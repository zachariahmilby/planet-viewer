from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import pytz
import spiceypy as spice
from astropy.coordinates import Latitude, Angle, Longitude, EarthLocation
from astropy.time import Time
from spiceypy.utils.callbacks import SpiceUDFUNS, SpiceUDFUNB
from timezonefinder import TimezoneFinder

from planetviewer.files import _project_directory
from planetviewer.spice_functions import get_azel, get_equatorial_radius, \
    get_flattening_coefficient
from planetviewer.spice_kernels import _earth_obs_kernel_from_coords, \
    _find_code
from planetviewer.objects import SolarSystemBody


class ObservingSite:
    """
    A general class for an observing site. Can be located on any SPICE
    ephemeris object (Earth-based sites have additional functionality).

    If no coordinates provided (latitude, longitude and altitude are all
    `None`), it will attempt to find the observatory in a list of pre-defined
    locations. If that fails (which would likely occur for minor observatories
    on Earth or any off-Earth site), you will need to provide coordinates.
    """
    def __init__(self,
                 name: str,
                 body: str = 'Earth',
                 latitude: u.Quantity | Angle | Latitude = None,
                 longitude: u.Quantity | Angle | Longitude = None,
                 altitude: u.Quantity = None,
                 overwrite_existing: bool = False):
        """
        Parameters
        ----------
        name : str
            The name of the observing site.
        body : str, optional
            The body on which the observing site is located. Defaults to
            'Earth'.
        latitude : u.Quantity or Angle or Latitude, optional
            The geodetic latitude of the observing site.
        longitude : u.Quantity or Angle or Latitude, optional
            The geodetic longitude of the observing site.
        altitude : u.Quantity, optional
            The altitude of the observing site.
        overwrite_existing : bool, optional
            Whether to overwrite observing site kernels (if they exist).
        """
        self._name = name
        self._body = body
        self._parse_site(latitude, longitude, altitude)
        self._set_timezone()
        self._furnish_kernel(overwrite=overwrite_existing)

    def _parse_site(self, latitude, longitude, altitude) -> None:
        """
        Set up an observing site. If no coordinates provided (latitude,
        longitude and altitude are all `None`), it will attempt to find the
        observatory in the list of pre-defined locations. If it can't find it
        there, it will attempt to resolve the location using Astropy's
        `EarthLocation.of_site()` method.

        Parameters
        ----------
        latitude : u.Quantity or Angle or Latitude, optional
            The geodetic latitude of the observing site.
        longitude : u.Quantity or Angle or Latitude, optional
            The geodetic longitude of the observing site.
        altitude : u.Quantity, optional
            The altitude of the observing site.

        Returns
        -------
        None
            None. Just sets up class attributes.
        """
        if (latitude is None) | (longitude is None) | (altitude is None):
            file = Path(_project_directory, 'anc', 'observatories.dat')
            sites = pd.read_csv(file)
            row = sites[sites['key'] == self._name]
            if len(row) == 0:
                try:
                    location = EarthLocation.of_site(self._name)
                    self._latitude = location.lat
                    self._longitude = location.lon
                    self._altitude = location.height.to(u.km)
                    self._earth_location = location
                except:
                    msg = ("Must supply either a pre-defined observatory site "
                           "or specify an observatory site's latitude, "
                           "longitude, altitude and the Solar System body on "
                           "which the site is found.")
                    raise ValueError(msg)
            self._latitude = Latitude(row['lat'].values[0] * u.degree)
            self._longitude = Longitude(row['lon'].values[0] * u.degree)
            self._altitude = row['alt'].values[0] * u.km
            self._location = self._make_astropy_object()

        else:
            self._latitude = Latitude(latitude)
            self._longitude = Longitude(longitude)
            self._altitude = altitude.to(u.km)
            self._location = self._make_astropy_object()
        self._code = _find_code(self._body)

    def _set_timezone(self) -> pytz.timezone:
        """
        If the observing site is located on Earth, determine its timezone.
        For non-Earth sites the timezone is set to UTC.

        Returns
        -------
        pytz.timezone
            The observing site's timezone.
        """
        if self._body == 'Earth':
            lon = self._longitude.degree
            if lon > 180:
                lon -= 360
            tf = TimezoneFinder().timezone_at(lat=self._latitude.degree,
                                              lng=lon)
            self._timezone = pytz.timezone(tf)
        else:
            self._timezone = pytz.utc

    def _make_astropy_object(self) -> EarthLocation | None:
        """
        If the observing site is located on Earth, set up an Astropy
        `EarthLocation` object to aid in calculating sidereal time.

        Returns
        -------
        EarthLocation or None
            The `EarthLocation` object. If the site is not on Earth, None is
            returned.
        """
        if self._body == 'Earth':
            return EarthLocation.from_geodetic(lat=self._latitude,
                                               lon=self._longitude,
                                               height=self._altitude)
        else:
            return None

    def _furnish_kernel(self,
                        overwrite: bool) -> None:
        """
        Make SPICE kernels necessary for this site to be used as an observer.

        Parameters
        ---------
        overwrite : bool
            Whether to overwrite observing site kernels (if they exist).

        Returns
        -------
        None
            None.
        """
        kernels = _earth_obs_kernel_from_coords(self._name,
                                                latitude=self._latitude,
                                                longitude=self._longitude,
                                                altitude=self._altitude,
                                                overwrite=overwrite)
        [kernel.furnish() for kernel in kernels]

    def _get_target_event(self,
                          target: str,
                          event: str,
                          time: Time,
                          which: str,
                          horizon: float = None) -> Time | None:
        """
        Use SPICE geometry finder to determine the precise time of a target's
        rise time, set time or meridian transits.

        Parameters
        ----------
        target : str
            The name of a Solar System ephemeris object. Must be available in
            SPICE.
        event : str
            The type of event to look for. Options are 'rise', 'set', 'transit'
            or 'antitransit'. Transit is the meridian crossing, antitransit is
            the antimeridian crossing.
        time : Time
            The reference time.
        which : str
            The timing of the event. Options are 'nearest' for the one closest
            to the reference time, 'previous' for the closest before the
            reference time and 'next' for the closest after the reference time.
        horizon : float, optional
            A specific horizon for calculating the time of the event. If None,
            then the refracted true horizon is used.

        Returns
        -------
        Time or None
            The UTC time of the target event. If no event is found (like no
            sunset near Earth's north pole during summer), None is returned.
        """
        et = spice.utc2et(time.isot)
        if which == 'nearest':
            left = et - 86400 / 2
            right = et + 86400 / 2
        elif which == 'previous':
            left = et - 86400
            right = et
        elif which == 'next':
            left = et
            right = et + 86400
        else:
            msg = ("Argument 'which' must be one of 'nearest' or 'previous' "
                   "or 'next'.")
            raise ValueError(msg)

        cnfine = spice.cell_double(86400)
        result = spice.cell_double(86400)
        spice.wninsd(left=left, right=right, window=cnfine)

        @SpiceUDFUNS
        def udfuns(et_in: float) -> float:
            return get_azel(target, et_in, self._name)[1] * spice.dpr()

        @SpiceUDFUNB
        def udqdec(udfunc, et_in: float) -> bool:
            return spice.uddc(udfunc, et_in, 1.0)

        # need to account for angular size of Sun, atmospheric refraction at
        # the horizon (for Earth) and altitude above the reference sphereoid

        # apparent angular size of Sun
        sun = SolarSystemBody(target)
        size = sun.get_angular_radius(time, self._name).to(u.degree).value

        # refraction angle used by NOAA
        # https://gml.noaa.gov/grad/solcalc/calcdetails.html
        if self._body == 'Earth':
            refraction = 0.833
        else:
            refraction = 0.0

        # additional angular requirement due to altitude above reference
        # sphereoid
        re = get_equatorial_radius(self._body, self._longitude.rad)
        f = get_flattening_coefficient(self._body, self._longitude.rad)
        rectan = spice.georec(lon=self._longitude.rad,
                              lat=self._latitude.rad,
                              alt=0,
                              re=re,
                              f=f)
        r0 = spice.vnorm(rectan)
        rectan = spice.georec(lon=self._longitude.rad,
                              lat=self._latitude.rad,
                              alt=self._altitude.value,
                              re=re,
                              f=f)
        r1 = spice.vnorm(rectan)
        altitude_offset = 90 - np.degrees(np.arcsin(r0 / r1))

        # use geometry finder to find when these events occur
        # determine geometry finder parameters for events
        if event in ['rise', 'set']:
            relate = '='
            if horizon is None:
                horizon = 0-size-refraction-altitude_offset
        elif event == 'antitransit':
            relate = 'ABSMIN'
            if horizon is None:
                horizon = 0
        elif event == 'transit':
            relate = 'ABSMAX'
            if horizon is None:
                horizon = 0
        else:
            msg = ("'event' must be one of 'rise', 'set', 'transit' or "
                   "'antitransit'.")
            raise ValueError(msg)

        result = spice.gfuds(udfuns=udfuns,
                             udqdec=udqdec,
                             relate=relate,
                             refval=horizon,
                             adjust=0,
                             step=60,
                             nintvls=2*86400,
                             cnfine=cnfine,
                             result=result)
        count = spice.wncard(result)
        if count == 0:
            return None
        for i in range(count):
            et = spice.wnfetd(result, i)[0]
            az1 = get_azel(target, et+0.5, self._name)[1] * spice.dpr()
            az0 = get_azel(target, et-0.5, self._name)[1] * spice.dpr()
            dalt = np.sign(az1 - az0)
            utc = spice.et2utc(et, format_str='ISOC', prec=3)
            out_time = Time(Time(utc).iso)
            if event == 'set':
                if dalt == -1.0:
                    return out_time
            elif event == 'rise':
                if dalt == 1.0:
                    return out_time
            elif event in ['transit', 'antitransit']:
                return out_time

    def convert_to_local_time(self,
                              utc: Time,
                              local_timezone: str = None) -> Time:
        """
        Convert a UTC time to a local time.

        Parameters
        ----------
        utc : Time
            UTC time to convert.
        local_timezone : str, optional
            The name of the local timezone. If None, defaults to the local
            timezone of the observing site.

        Returns
        -------

        """
        if local_timezone is None:
            timezone = self._timezone
        else:
            timezone = pytz.timezone(local_timezone)
        local_time = utc.to_datetime(timezone=pytz.utc).astimezone(timezone)
        fmt = '%Y-%m-%dT%H:%M:%S.%f'
        return Time(Time(local_time.strftime(fmt)).iso)

    def get_sunset_time(self,
                        reftime: Time,
                        which: str,
                        local_time: bool = False,
                        timezone: str = None) -> Time:
        """
        Get UTC sunset time relative to a UTC reference time. Accounts for
        the angular size of the Sun, atmospheric refraction at the horizon (for
        Earth locations only) and the relative altitude of an observer.

        Parameters
        ----------
        reftime : Time
            UTC reference time.
        which : str
            Which sunset you want to find. Options are 'previous' for the one
            preceeding the reference time, 'next' for the one following the
            reference time or 'nearest' for the one closest in time.
        local_time : bool, optional
            Whether or not to convert from UTC to a local time. Must specify a
            timezone. Default is False.
        timezone: str, optional
            If you want a local time instead of a UTC time, set the timezone
            here. If `local_time=True` and `timezone=None`, it will default to
            the local timezone at the observing size. If you want somewhere
            else, you can pass any `pytz` timezone.

        Returns
        -------
        Time
            The UTC sunset time as an Astropy `Time` object.
        """
        time = self._get_target_event(target='Sun',
                                      event='set',
                                      time=reftime,
                                      which=which)
        if local_time:
            time = self.convert_to_local_time(time, timezone)
        return time

    def get_solar_midnight_time(self,
                                reftime: Time,
                                which: str,
                                local_time: bool = False,
                                timezone: str = None) -> Time:
        """
        Get UTC solar midnight time relative to a UTC reference time. Accounts
        for the angular size of the Sun, atmospheric refraction at the horizon
        (for Earth locations only) and the relative altitude of an observer.

        Parameters
        ----------
        reftime : Time
            UTC reference time.
        which : str
            Which midnight you want to find. Options are 'previous' for the one
            preceeding the reference time, 'next' for the one following the
            reference time or 'nearest' for the one closest in time.
        local_time : bool, optional
            Whether or not to convert from UTC to a local time. Must specify a
            timezone. Default is False.
        timezone: str, optional
            If you want a local time instead of a UTC time, set the timezone
            here. If `local_time=True` and `timezone=None`, it will default to
            the local timezone at the observing size. If you want somewhere
            else, you can pass any `pytz` timezone.

        Returns
        -------
        Time
            The UTC midnight time as an Astropy `Time` object.
        """
        time = self._get_target_event(target='Sun',
                                      event='antitransit',
                                      time=reftime,
                                      which=which)
        if local_time:
            time = self.convert_to_local_time(time, timezone)
        return time

    def get_sunrise_time(self,
                         reftime: Time,
                         which: str,
                         local_time: bool = False,
                         timezone: str = None) -> Time:
        """
        Get UTC sunrise time relative to a UTC reference time. Accounts for
        the angular size of the Sun, atmospheric refraction at the horizon (for
        Earth locations only) and the relative altitude of an observer.

        Parameters
        ----------
        reftime : Time
            UTC reference time.
        which : str
            Which sunrise you want to find. Options are 'previous' for the one
            preceeding the reference time, 'next' for the one following the
            reference time or 'nearest' for the one closest in time.
        local_time : bool, optional
            Whether or not to convert from UTC to a local time. Must specify a
            timezone. Default is False.
        timezone: str, optional
            If you want a local time instead of a UTC time, set the timezone
            here. If `local_time=True` and `timezone=None`, it will default to
            the local timezone at the observing size. If you want somewhere
            else, you can pass any `pytz` timezone.

        Returns
        -------
        Time
            The UTC sunset time as an Astropy `Time` object.
        """
        time = self._get_target_event(target='Sun',
                                      event='rise',
                                      time=reftime,
                                      which=which)
        if local_time:
            time = self.convert_to_local_time(time, timezone)
        return time

    def get_solar_noon_time(self,
                            reftime: Time,
                            which: str,
                            local_time: bool = False,
                            timezone: str = None) -> Time:
        """
        Get UTC solar noon time relative to a UTC reference time. Accounts
        for the angular size of the Sun, atmospheric refraction at the horizon
        (for Earth locations only) and the relative altitude of an observer.

        Parameters
        ----------
        reftime : Time
            UTC reference time.
        which : str
            Which solar noon you want to find. Options are 'previous' for the
            one preceeding the reference time, 'next' for the one following the
            reference time or 'nearest' for the one closest in time.
        local_time : bool, optional
            Whether or not to convert from UTC to a local time. Must specify a
            timezone. Default is False.
        timezone: str, optional
            If you want a local time instead of a UTC time, set the timezone
            here. If `local_time=True` and `timezone=None`, it will default to
            the local timezone at the observing size. If you want somewhere
            else, you can pass any `pytz` timezone.

        Returns
        -------
        Time
            The UTC solar noon time as an Astropy `Time` object.
        """
        time = self._get_target_event(target='Sun',
                                      event='transit',
                                      time=reftime,
                                      which=which)
        if local_time:
            time = self.convert_to_local_time(time, timezone)
        return time

    def get_moonrise_time(self,
                          reftime: Time,
                          which: str,
                          local_time: bool = False,
                          timezone: str = None) -> Time:
        """
        Get UTC moonrise time relative to a UTC reference time. Accounts for
        the angular size of the Moon, atmospheric refraction at the horizon
        (for Earth locations only) and the relative altitude of an observer.

        Parameters
        ----------
        reftime : Time
            UTC reference time.
        which : str
            Which moon rise time you want to find. Options are 'previous' for
            the one preceeding the reference time, 'next' for the one following
            the reference time or 'nearest' for the one closest in time.
        local_time : bool, optional
            Whether or not to convert from UTC to a local time. Must specify a
            timezone. Default is False.
        timezone: str, optional
            If you want a local time instead of a UTC time, set the timezone
            here. If `local_time=True` and `timezone=None`, it will default to
            the local timezone at the observing size. If you want somewhere
            else, you can pass any `pytz` timezone.

        Returns
        -------
        Time
            The UTC moon rise time as an Astropy `Time` object.
        """
        time = self._get_target_event(target='Moon',
                                      event='rise',
                                      time=reftime,
                                      which=which)
        if local_time:
            time = self.convert_to_local_time(time, timezone)
        return time

    def get_moonset_time(self,
                         reftime: Time,
                         which: str,
                         local_time: bool = False,
                         timezone: str = None) -> Time:
        """
        Get UTC moonset time relative to a UTC reference time. Accounts for the
        angular size of the Moon, atmospheric refraction at the horizon (for
        Earth locations only) and the relative altitude of an observer.

        Parameters
        ----------
        reftime : Time
            UTC reference time.
        which : str
            Which moon set time you want to find. Options are 'previous' for
            the one preceeding the reference time, 'next' for the one following
            the reference time or 'nearest' for the one closest in time.
        local_time : bool, optional
            Whether or not to convert from UTC to a local time. Must specify a
            timezone. Default is False.
        timezone: str, optional
            If you want a local time instead of a UTC time, set the timezone
            here. If `local_time=True` and `timezone=None`, it will default to
            the local timezone at the observing size. If you want somewhere
            else, you can pass any `pytz` timezone.

        Returns
        -------
        Time
            The UTC moon set time as an Astropy `Time` object.
        """
        time = self._get_target_event(target='Moon',
                                      event='set',
                                      time=reftime,
                                      which=which)
        if local_time:
            time = self.convert_to_local_time(time, timezone)
        return time

    def get_moon_illumination(self,
                              time: Time) -> float:
        """
        Get fraction of Moon illuminated as observed from this site.

        Parameters
        ----------
        time : Time
            UTC reference time.

        Returns
        -------
        float
            The fraction of Moon illuminated.
        """
        moon = SolarSystemBody('Moon')
        phase_angle = moon.get_phase_angle(time, self._name)
        k = (1 + np.cos(phase_angle)) / 2.0
        return k.value

    def get_local_sidereal_time(self,
                                local_time: Time,
                                kind: str) -> Longitude:
        """
        Convert a local time to a local sidereal time.

        Parameters
        ----------
        local_time : Time
            The local UTC time.
        kind : str
            Type of calculation. If 'mean', then the calculation includes
            precession but not nutation. If 'apparent', it also includes
            nutation.

        Returns
        -------
        Longitude
            The local sidereal time as an Astropy `Longitude` object with units
            of hourangle.
        """
        return local_time.sidereal_time(kind=kind, longitude=self._longitude)

    def get_hour_angle(self,
                       target: str,
                       time: Time,
                       kind: str) -> Angle:
        """
        Get the hour angle for a given target/UTC time for this site.

        Parameters
        ----------
        target : str
            The name of a Solar System ephemeris object. Must be available in
            SPICE.
        time : Time
            The local UTC time.
        kind : str
            Type of calculation. If 'mean', then the calculation includes
            precession but not nutation. If 'apparent', it also includes
            nutation.

        Returns
        -------
        Angle
            The target's hour hangle as an Astropy `Angle` object.
        """
        lst = self.get_local_sidereal_time(local_time=time, kind=kind)
        coord = SolarSystemBody(target).get_skycoord(time, self._name,
                                                     kind=kind)
        hour_angle = Angle(lst - coord.ra)
        hour_angle.wrap_at('12h', inplace=True)
        return hour_angle

    def get_parallactic_angle(self,
                              target: str,
                              time: Time,
                              kind: str) -> Angle:
        """
        Get the parallactic angle for a given target/UTC time for this site.

        Parameters
        ----------
        target : str
            The name of a Solar System ephemeris object. Must be available in
            SPICE.
        time : Time
            The local UTC time.
        kind : str
            Type of calculation. If 'mean', then the calculation includes
            precession but not nutation. If 'apparent', it also includes
            nutation.

        Returns
        -------
        Angle
            The local sidereal time as an Astropy `Angle` object.
        """
        coord = SolarSystemBody(target).get_skycoord(time, self._name,
                                                     kind=kind)
        hour_angle = self.get_hour_angle(target=target, time=time, kind=kind)
        y = np.sin(hour_angle)
        x = (np.tan(self._latitude) * np.cos(coord.dec) -
             np.sin(coord.dec) * np.cos(hour_angle))
        q = (np.arctan2(y, x)).to(u.degree)
        return Angle(q)

    def get_twilight_time(self,
                          reftime: Time,
                          which: str,
                          kind: str,
                          local_time: bool = False,
                          timezone: str = None) -> Time:
        """
        Get the time for different kinds of twilight.

        Parameters
        ----------
        reftime : Time
            UTC reference time.
        which : str
            Which moon set time you want to find. Options are 'previous' for
            the one preceeding the reference time, 'next' for the one following
            the reference time or 'nearest' for the one closest in time.
        kind : str
            The type of twilight: options are 'civil dusk', 'nautical dusk',
            'astronomical dusk', 'astronomical dawn', 'nautical dawn' and
            'civil dawn'.
        local_time : bool, optional
            Whether or not to convert from UTC to a local time. Must specify a
            timezone. Default is False.
        timezone: str, optional
            If you want a local time instead of a UTC time, set the timezone
            here. If `local_time=True` and `timezone=None`, it will default to
            the local timezone at the observing size. If you want somewhere
            else, you can pass any `pytz` timezone.

        Returns
        -------
        Time
            The UTC twilight time as an Astropy `Time` object, optionally
            converted to local time if `local_time=True`.
        """
        if 'dusk' in kind:
            event = 'set'
        elif 'dawn' in kind:
            event = 'rise'
        else:
            msg = 'Must specify whether a twilight occurs at dusk or dawn.'
            raise ValueError(msg)

        if 'civil' in kind:
            horizon = -6
        elif 'nautical' in kind:
            horizon = -12
        elif 'astronomical' in kind:
            horizon = -18
        else:
            msg = ('Must specify a subcategory of twilight: civil, nautical or '
                   'astronomical.')
            raise ValueError(msg)

        time = self._get_target_event(target='Sun',
                                      event=event,
                                      time=reftime,
                                      which=which,
                                      horizon=horizon)
        if local_time:
            time = self.convert_to_local_time(time, timezone)
        return time

    @property
    def name(self) -> str:
        """The observing site's name."""
        return self._name

    @property
    def code(self) -> int:
        """The observing site's assigned SPICE integer ID code."""
        return self._code

    @property
    def body(self) -> str:
        """The body on which the observing site is located."""
        return self._body

    @property
    def latitude(self) -> Latitude:
        """The observing site's latitude."""
        return self._latitude

    @property
    def longitude(self) -> Longitude:
        """The observing site's east longitude."""
        return self._longitude

    @property
    def altitude(self) -> u.Quantity:
        """The observing site's altitude."""
        return self._altitude

    @property
    def earth_location(self) -> EarthLocation | None:
        """The Astropy `EarthLocation` object for Earth-based observatories."""
        return self._location
