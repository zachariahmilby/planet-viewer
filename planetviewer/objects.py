from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spiceypy
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, Latitude, Longitude
from astropy.time import Time
from astropy.visualization.wcsaxes.core import WCSAxes
from spiceypy.utils.exceptions import NotFoundError

from planetviewer.files import _project_directory
from planetviewer.plotting import (
    plot_limb, plot_disk, plot_nightside, plot_latlon, plot_ring, plot_arc,
    place_ring_pericenter_markers, plot_primary_shadow,
    _default_night_patch_kwargs)
from planetviewer.spice_functions import (
    determine_prograde_rotation, get_radii, get_sky_coordinates, angle_unit,
    get_sub_observer_latlon, get_sub_solar_latlon,
    get_sub_observer_phase_angle, get_state, length_unit, get_light_time,
    get_et, get_limb_radec, get_latlon_radec, get_terminator_radec,
    get_rotation_rate, abcorr, ref, angle_rate_unit,
    determine_eclipse_or_transit, get_shadow_intercept, get_dayside_radec)


def _r(semimajor_axis: u.Quantity,
       eccentricity: float,
       true_anomaly: u.Quantity) -> u.Quantity:
    """
    Calculate radial disance as a function of true anomaly and semimajor
    axis.

    Parameters
    ----------
    semimajor_axis : u.Quantity
        Ring semimajor axis. Set as an input variable rather than using the
        value in the class constructor so that you can calculate inner and
        outer edges of the ring rather than just the ring center.
    true_anomaly : u.Quantity
        Angle from ring pericenter.

    Returns
    -------
    The radial disance as a function of true anomaly and semimajor axis.
    """
    numerator = semimajor_axis * (1 - eccentricity ** 2)
    denominator = 1 + eccentricity * np.cos(true_anomaly)
    out = numerator / denominator
    return out.to(length_unit)


def _ring_xyz(time: Time,
              observer: str,
              planet: str,
              semimajor_axis: u.Quantity,
              width: u.Quantity,
              eccentricity: float,
              inclination: u.Quantity,
              longitude_of_periapsis: u.Quantity,
              longitude_of_ascending_node: u.Quantity,
              dlon_dt: u.Quantity,
              dnode_dt: u.Quantity,
              step: u.Quantity = 1 * u.degree,
              pericenter: bool = False,
              edge: str = None,
              minlon: u.Quantity = 0 * u.degree,
              dminlon_dt: u.Quantity = 0 * u.degree / u.day,
              delta_lon: u.Quantity = 0 * u.degree,
              boundaries: bool = False,
              ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate ring positions vectors in planetary reference frame.

    Parameters
    ----------
    time : Time
        UTC at the time of observation.
    observer : str
        The observer or observatory. Could be a Solar System body like
        "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
        like "Juno".
    planet : str
        Name of parent planet.
    semimajor_axis : u.Quantity
        Ring's semimajor axis.
    width : u.Quantity
        Ring's width.
    eccentricity : float
        Ring's eccentricity.
    inclination : u.Quantity
        Ring's inclination.
    longitude_of_periapsis : u.Quantity
        Ring's longitude of periapsis.
    longitude_of_ascending_node : u.Quantity
        Ring's longitude of ascending node.
    dlon_dt : u.Quantity
        Precession rate of ring's longitude of periapsis.
    dnode_dt : u.Quantity
        Precession rate of ring's longitude of ascending node.
    step : u.Quantity
        Angular precision of the calculation. Default is 1 degree.
    pericenter : bool
        If true, return just the position of the ring's pericenter. Applies
        only to some of the rings of Uranus with non-zero eccentricity.
    edge : str
        For rings with non-zero width, whether you want the 'inner' or 'outer'
        edge.
    minlon : u.Quantity
        If you want just a subsection of the ring, use this to set the starting
        longitude. This is useful for selecting ring arc segments.
    dminlon_dt : u.Quantity
        The precession rate of `minlon`. Again, useful for ring arc segments
        with a measured precession rate.
    delta_lon : u.Quantity
        The angular length of the segment starting from `minlon`.
    boundaries : bool
        If true, returns an array with size n+1 so that you get the left and
        right boundaries of each resolution element.

    Returns
    -------
    A tuple (x, y, z) containing the vector components.
    """

    if edge == 'inner':
        a = semimajor_axis - width / 2
    elif edge == 'outer':
        a = semimajor_axis + width / 2
    else:
        a = semimajor_axis

    et = get_et(time)
    light_travel_time = get_light_time(planet, et, observer)
    delta_t = (et - light_travel_time) * u.s
    prograde = determine_prograde_rotation(planet)

    varpi = deepcopy(longitude_of_periapsis)
    big_omega = deepcopy(longitude_of_ascending_node)
    minlon = deepcopy(minlon)

    # add angular precession
    varpi += dlon_dt * delta_t
    big_omega += dnode_dt * delta_t
    minlon += dminlon_dt * delta_t

    # flip for retrograde rotation (Uranus)
    if not prograde:
        varpi = 360 * u.degree - varpi
        big_omega = 360 * u.degree - big_omega
        minlon = 360 * u.degree - minlon

    # remove planetary rotation
    parent_rotation_rate = get_rotation_rate(
        planet) * angle_rate_unit
    varpi -= parent_rotation_rate * delta_t
    big_omega -= parent_rotation_rate * delta_t
    minlon -= parent_rotation_rate * delta_t

    varpi %= 360 * u.degree
    big_omega %= 360 * u.degree
    minlon %= 360 * u.degree

    # calculate argument of periapsis
    omega = varpi - big_omega

    # calculate true anomaly
    if not pericenter:
        if boundaries:
            theta = np.arange(
                minlon.to(u.degree).value-step.value/2,
                (minlon + delta_lon).to(u.degree).value+step.value/2,
                step.value) * u.degree
        else:
            theta = np.arange(minlon.to(u.degree).value,
                              (minlon + delta_lon).to(u.degree).value,
                              step.value) * u.degree
    else:
        theta = varpi
    f = theta - varpi

    r = _r(a, eccentricity, f)
    ang_arg = omega + f

    term1 = np.cos(big_omega) * np.cos(ang_arg)
    term2 = np.sin(big_omega) * np.sin(ang_arg) * np.cos(inclination)
    term3 = np.sin(big_omega) * np.cos(ang_arg)
    term4 = np.cos(big_omega) * np.sin(ang_arg) * np.cos(inclination)
    term5 = np.sin(ang_arg) * np.sin(inclination)

    x = r * (term1 - term2)
    y = r * (term3 + term4)
    z = r * term5

    x = x.to(u.km).value
    y = y.to(u.km).value
    z = z.to(u.km).value

    if pericenter and eccentricity == 0.0:
        return np.array([np.nan]), np.array([np.nan]), np.array([np.nan])
    elif pericenter and eccentricity != 0.0:
        return np.array([x]), np.array([y]), np.array([z])
    else:
        return x, y, z


class Arc:
    """
    An arc within a planetary ring. Currently applicable only to Neptune.
    Assumes provided parameters have been measured at the same epoch as the
    parent ring.
    """
    def __init__(self,
                 name: str,
                 min_longitude: u.Quantity,
                 max_longitude: u.Quantity,
                 dlon_dt: u.Quantity):
        """
        Parameters
        ----------
        name : str
            The name of the arc.
        min_longitude : u.Quantity
            The minimum longitude of the arc.
        max_longitude : u.Quantity
            The maximum longitude of the arc.
        dlon_dt : u.Quantity
            The precession rate of the arc in units of [angle/time].
        """
        self._name = name
        self._min_longitude = min_longitude
        if max_longitude < min_longitude:
            max_longitude += 360 * u.degree
        self._max_longitude = max_longitude
        self._dlon_dt = dlon_dt

    @property
    def name(self) -> str:
        return self._name

    @property
    def min_longitude(self) -> u.Quantity:
        return self._min_longitude

    @property
    def delta_longitude(self) -> u.Quantity:
        return self._max_longitude - self._min_longitude

    @property
    def dlon_dt(self) -> u.Quantity:
        return self._dlon_dt


class Ring:
    """
    A planetary ring. Assumes provided parameters have been measured or offset
    to the J2000 reference time of 2000-01-01 12:00 UT at the parent body (not
    the observer).
    """
    def __init__(self,
                 name: str,
                 planet: str,
                 planet_prograde: bool,
                 semimajor_axis: u.Quantity,
                 width: u.Quantity = 0 * u.km,
                 eccentricity: float = 0.0,
                 inclination: u.Quantity = 0 * u.degree,
                 longitude_of_periapsis: u.Quantity = 0.0 * u.degree,
                 longitude_of_ascending_node: u.Quantity = 0.0 * u.degree,
                 dlon_dt: u.Quantity = 0.0 * u.degree / u.day,
                 dnode_dt: u.Quantity = 0.0 * u.degree / u.day,
                 arcs: dict[str, Arc] = None):
        """
        Parameters
        ----------
        name : str
            The name of the ring.
        planet : str
            The name of the parent body.
        planet_prograde : bool
            Whether or not the parent body's rotation is prograde.
        semimajor_axis : u.Quantity
            Ring's semimajor axis.
        width : u.Quantity
            Ring's width.
        eccentricity : float
            Ring's eccentricity.
        inclination : u.Quantity
            Ring's inclination.
        longitude_of_periapsis : u.Quantity
            Ring's longitude of periapsis.
        longitude_of_ascending_node : u.Quantity
            Ring's longitude of ascending node.
        dlon_dt : u.Quantity
            Precession rate of ring's longitude of periapsis.
        dnode_dt : u.Quantity
            Precession rate of ring's longitude of ascending node.
        arcs : dict[str, Arc]
            A dictionary containing any ring arcs. The keys should be the name
            of the arcs.
        """
        self._name = name
        self._planet = planet
        self._fixref = f'IAU_{planet}'
        self._planet_prograde = planet_prograde
        self._semimajor_axis = semimajor_axis
        self._width = width
        self._eccentricity = eccentricity
        self._inclination = inclination
        self._longitude_of_periapsis = longitude_of_periapsis
        self._longitude_of_ascending_node = longitude_of_ascending_node
        self._dlon_dt = dlon_dt
        self._dnode_dt = dnode_dt
        if arcs is None:
            arcs = {}
        self._arcs = arcs

    def _get_xyz(self,
                 time: Time,
                 observer: str,
                 step: u.Quantity = 1 * u.degree,
                 edge: str = None,
                 pericenter: bool = False,
                 minlon: u.Quantity = 0 * u.degree,
                 dminlon_dt: u.Quantity = 0 * u.degree / u.day,
                 delta_lon: u.Quantity = 360 * u.degree,
                 boundaries: bool = False
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Get ring positions vectors in planetary reference frame.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        step : u.Quantity
            Angular precision of the calculation. Default is 1 degree.
        pericenter : bool
            If true, return just the position of the ring's pericenter. Applies
            only to some of the rings of Uranus with non-zero eccentricity.
        edge : str
            For rings with non-zero width, whether you want the 'inner' or
            'outer' edge.
        minlon : u.Quantity
            If you want just a subsection of the ring, use this to set the
            starting longitude. This is useful for selecting ring arc segments.
        dminlon_dt : u.Quantity
            The precession rate of `minlon`. Again, useful for ring arc
            segments with a measured precession rate.
        delta_lon : u.Quantity
            The angular length of the segment starting from `minlon`.
        boundaries : bool
            If true, returns an array with size n+1 so that you get the left
            and right boundaries of each resolution element.

        Returns
        -------
        A tuple (x, y, z) containing the vector components.
        """
        args = dict(
            time=time,
            observer=observer,
            planet=self._planet,
            semimajor_axis=self._semimajor_axis,
            width=self._width,
            eccentricity=self._eccentricity,
            inclination=self._inclination,
            longitude_of_periapsis=self._longitude_of_periapsis,
            longitude_of_ascending_node=self._longitude_of_ascending_node,
            dlon_dt=self._dlon_dt,
            dnode_dt=self._dnode_dt,
            step=step,
            pericenter=pericenter,
            edge=edge,
            minlon=minlon,
            dminlon_dt=dminlon_dt,
            delta_lon=delta_lon,
            boundaries=boundaries)
        return _ring_xyz(**args)

    # noinspection DuplicatedCode
    def _get_sky_coordinates(
            self,
            time: Time,
            observer: str,
            step: u.Quantity = 1 * u.degree,
            edge: str = None,
            pericenter: bool = False,
            minlon: u.Quantity = 0 * u.degree,
            dminlon_dt: u.Quantity = 0 * u.degree / u.day,
            delta_lon: u.Quantity = 360 * u.degree,
            boundaries: bool = False) -> list[SkyCoord] | None:
        """
        Get ring sky coordinates.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        step : u.Quantity
            Angular precision of the calculation. Default is 1 degree.
        pericenter : bool
            If true, return just the position of the ring's pericenter. Applies
            only to some of the rings of Uranus with non-zero eccentricity.
        edge : str
            For rings with non-zero width, whether you want the 'inner' or
            'outer' edge.
        minlon : u.Quantity
            If you want just a subsection of the ring, use this to set the
            starting longitude. This is useful for selecting ring arc segments.
        dminlon_dt : u.Quantity
            The precession rate of `minlon`. Again, useful for ring arc
            segments with a measured precession rate.
        delta_lon : u.Quantity
            The angular length of the segment starting from `minlon`.
        boundaries : bool
            If true, returns an array with size n+1 so that you get the left
            and right boundaries of each resolution element.

        Returns
        -------
        A list containing the ring sky coordinates as Astropy `SkyCoord`
        objects.
        """
        et = get_et(time)
        args = dict(time=time,
                    observer=observer,
                    step=step,
                    pericenter=pericenter,
                    edge=edge,
                    minlon=minlon,
                    dminlon_dt=dminlon_dt,
                    delta_lon=delta_lon,
                    boundaries=boundaries)
        x, y, z = self._get_xyz(**args)
        if (x.size == 1) & (np.isnan(x).any()):
            return None
        coords = []
        params = dict(method='ELLIPSOID',
                      target=self._planet,
                      ilusrc='Sun',
                      et=et,
                      fixref=self._fixref,
                      abcorr=abcorr,
                      obsrvr=observer)
        for i in range(x.size):
            spoint = spiceypy.vpack(float(x[i]), float(y[i]), float(z[i]))
            trgepc, vec, _, _, _, visible, lit = spiceypy.illumf(
                spoint=spoint, **params)
            rotation_matrix = spiceypy.pxfrm2(
                self._fixref, 'J2000', trgepc, et)
            point = spiceypy.mxv(rotation_matrix, vec)
            _, ra_i, dec_i = spiceypy.recrad(point)
            coord = SkyCoord(ra=Angle(ra_i, unit='rad'),
                             dec=Angle(dec_i, unit='rad'))
            coords.append(coord)
        return coords

    # noinspection DuplicatedCode
    def _get_distances(self,
                       time: Time,
                       observer: str,
                       step: u.Quantity = 1 * u.degree,
                       edge: str = None,
                       pericenter: bool = False,
                       minlon: u.Quantity = 0 * u.degree,
                       dminlon_dt: u.Quantity = 0 * u.degree / u.day,
                       delta_lon: u.Quantity = 360 * u.degree,
                       boundaries: bool = False) -> u.Quantity | None:
        """
        Get distances from observer to each ring segment.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        step : u.Quantity
            Angular precision of the calculation. Default is 1 degree.
        pericenter : bool
            If true, return just the position of the ring's pericenter. Applies
            only to some of the rings of Uranus with non-zero eccentricity.
        edge : str
            For rings with non-zero width, whether you want the 'inner' or
            'outer' edge.
        minlon : u.Quantity
            If you want just a subsection of the ring, use this to set the
            starting longitude. This is useful for selecting ring arc segments.
        dminlon_dt : u.Quantity
            The precession rate of `minlon`. Again, useful for ring arc
            segments with a measured precession rate.
        delta_lon : u.Quantity
            The angular length of the segment starting from `minlon`.
        boundaries : bool
            If true, returns an array with size n+1 so that you get the left
            and right boundaries of each resolution element.

        Returns
        -------
        An array of distances from the observer to each ring segment as an
        Astropy `Quantity` object.
        """
        et = spiceypy.str2et(time.to_string())
        args = dict(time=time,
                    observer=observer,
                    step=step,
                    pericenter=pericenter,
                    edge=edge,
                    minlon=minlon,
                    dminlon_dt=dminlon_dt,
                    delta_lon=delta_lon,
                    boundaries=boundaries)
        x, y, z = self._get_xyz(**args)
        if (x.size == 1) & (np.isnan(x).any()):
            return None
        distances = np.full(x.shape, np.nan)
        for i in range(x.size):

            # get distance from observer to parent body in J2000
            starg, lt = spiceypy.spkezr(self._planet, et, ref, abcorr,
                                        observer)
            vec2 = starg[:3]

            # convert vector to J2000
            rotate = spiceypy.pxform(self._fixref, ref, et-lt)
            vec1 = spiceypy.vpack(float(x[i]), float(y[i]), float(z[i]))
            vec1 = spiceypy.mxv(rotate, vec1)

            # add vectors to get distance from observer to ring point
            dist_vec = vec1 + vec2

            # calculate magnitude of distance vector
            distances[i] = spiceypy.vnorm(dist_vec)

        return distances * length_unit

    # noinspection DuplicatedCode
    def _determine_eclipsed(self,
                            time: Time,
                            step: u.Quantity = 1 * u.degree,
                            edge: str = None,
                            pericenter: bool = False,
                            minlon: u.Quantity = 0 * u.degree,
                            dminlon_dt: u.Quantity = 0 * u.degree / u.day,
                            delta_lon: u.Quantity = 360 * u.degree,
                            boundaries: bool = False):
        """
        Determine if a ring segment is eclipsed by the planet relative to the
        Sun.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        step : u.Quantity
            Angular precision of the calculation. Default is 1 degree.
        pericenter : bool
            If true, return just the position of the ring's pericenter. Applies
            only to some of the rings of Uranus with non-zero eccentricity.
        edge : str
            For rings with non-zero width, whether you want the 'inner' or
            'outer' edge.
        minlon : u.Quantity
            If you want just a subsection of the ring, use this to set the
            starting longitude. This is useful for selecting ring arc segments.
        dminlon_dt : u.Quantity
            The precession rate of `minlon`. Again, useful for ring arc
            segments with a measured precession rate.
        delta_lon : u.Quantity
            The angular length of the segment starting from `minlon`.
        boundaries : bool
            If true, returns an array with size n+1 so that you get the left
            and right boundaries of each resolution element.

        Returns
        -------
        A boolean array indicating whether a ring segment is eclipsed or not.
        `True` means it is eclipsed. `False` means it is lit.
        """
        observer = 'Sun'
        et = spiceypy.str2et(time.to_string())
        args = dict(time=time,
                    observer=observer,
                    step=step,
                    pericenter=pericenter,
                    edge=edge,
                    minlon=minlon,
                    dminlon_dt=dminlon_dt,
                    delta_lon=delta_lon,
                    boundaries=boundaries)
        starg, lt = spiceypy.spkezr(self._planet, et, ref, abcorr, observer)
        parent_distance = spiceypy.vnorm(starg[:3]) * length_unit
        distances = self._get_distances(**args)
        x, y, z = self._get_xyz(**args)
        shaded = np.full(distances.shape, False, dtype=bool)
        params = dict(method='ELLIPSOID',
                      target=self._planet,
                      ilusrc='Sun',
                      et=et,
                      fixref=self._fixref,
                      abcorr=abcorr,
                      obsrvr=observer)
        for i in range(distances.shape[0]):
            if distances[i] < parent_distance:
                continue
            else:
                spoint = spiceypy.vpack(float(x[i]), float(y[i]), float(z[i]))
                trgepc, vec, _, _, _, _, _ = spiceypy.illumf(
                    spoint=spoint, **params)
                rotation_matrix = spiceypy.pxfrm2(self._fixref, 'J2000',
                                                  trgepc, et)
                rotvec = spiceypy.mxv(rotation_matrix, vec)
                try:
                    spiceypy.sincpt('ELLIPSOID', self._planet, et,
                                    self._fixref, abcorr, observer, 'J2000',
                                    rotvec)
                except NotFoundError:
                    continue
                else:
                    shaded[i] = True
        return shaded

    def get_sky_coordinates(self,
                            time: Time,
                            observer: str,
                            step: u.Quantity = 1 * u.degree,
                            pericenter: bool = False,
                            edge: str = None,
                            boundaries: bool = False) -> list[SkyCoord]:
        """
        Get the ring's sky coordinates for a given time and observer.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        step : u.Quantity
            Angular precision of the calculation. Default is 1 degree.
        pericenter : bool
            If true, return just the position of the ring's pericenter. Applies
            only to some of the rings of Uranus with non-zero eccentricity.
        edge : str
            For rings with non-zero width, whether you want the 'inner' or
            'outer' edge.
        boundaries : bool
            If true, returns an array with size n+1 so that you get the left
            and right boundaries of each resolution element.

        Returns
        -------
        A list containing the ring sky coordinates as Astropy `SkyCoord`
        objects.
        """
        args = dict(time=time,
                    observer=observer,
                    step=step,
                    edge=edge,
                    pericenter=pericenter,
                    boundaries=boundaries)
        return self._get_sky_coordinates(**args)

    def get_distances(self,
                      time: Time,
                      observer: str,
                      step: u.Quantity = 1 * u.degree,
                      pericenter: bool = False,
                      edge: str = None) -> u.Quantity:
        """
        Get the distances from the observer to each ring segment for a given
        time and observer.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        step : u.Quantity
            Angular precision of the calculation. Default is 1 degree.
        pericenter : bool
            If true, return just the position of the ring's pericenter. Applies
            only to some of the rings of Uranus with non-zero eccentricity.
        edge : str
            For rings with non-zero width, whether you want the 'inner' or
            'outer' edge.

        Returns
        -------
        An array of distances from the observer to each ring segment as an
        Astropy `Quantity` object.
        """
        args = dict(time=time,
                    observer=observer,
                    step=step,
                    edge=edge,
                    pericenter=pericenter)
        return self._get_distances(**args)

    def get_eclipsed(self,
                     time: Time,
                     step: u.Quantity = 1 * u.degree,
                     pericenter: bool = False,
                     edge: str = None):
        """
        Determine if a ring segment is eclipsed by the planet relative to the
        Sun for a given time and observer.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        step : u.Quantity
            Angular precision of the calculation. Default is 1 degree.
        pericenter : bool
            If true, return just the position of the ring's pericenter. Applies
            only to some of the rings of Uranus with non-zero eccentricity.
        edge : str
            For rings with non-zero width, whether you want the 'inner' or
            'outer' edge.

        Returns
        -------
        A boolean array indicating whether a ring segment is eclipsed or not.
        `True` means it is eclipsed. `False` means it is lit.
        """
        args = dict(time=time,
                    step=step,
                    edge=edge,
                    pericenter=pericenter)
        return self._determine_eclipsed(**args)

    def get_arc_sky_coordinates(self,
                                name: str,
                                time: Time,
                                observer: str,
                                step: u.Quantity = 1 * u.degree,
                                edge: str = None,
                                boundaries: bool = False) -> list[SkyCoord]:
        """
        Get a ring arc's sky coordinates for a given time and observer.

        Parameters
        ----------
        name : str
            The name of the ring arc.
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        step : u.Quantity
            Angular precision of the calculation. Default is 1 degree.
        edge : str
            For rings with non-zero width, whether you want the 'inner' or
            'outer' edge.
        boundaries : bool
            If true, returns an array with size n+1 so that you get the left
            and right boundaries of each resolution element.

        Returns
        -------
        A list containing the ring sky coordinates as Astropy `SkyCoord`
        objects.
        """
        arc = self._arcs[name]
        args = dict(time=time,
                    observer=observer,
                    step=step,
                    edge=edge,
                    pericenter=False,
                    minlon=arc.min_longitude,
                    dminlon_dt=arc.dlon_dt,
                    delta_lon=arc.delta_longitude,
                    boundaries=boundaries)
        return self._get_sky_coordinates(**args)

    def get_arc_distances(self,
                          name: str,
                          time: Time,
                          observer: str,
                          step: u.Quantity = 1 * u.degree,
                          edge: str = None) -> u.Quantity:
        """
        Get the distances from the observer to each ring arc segment for a
        given time and observer.

        Parameters
        ----------
        name : str
            The name of the ring arc.
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        step : u.Quantity
            Angular precision of the calculation. Default is 1 degree.
        edge : str
            For rings with non-zero width, whether you want the 'inner' or
            'outer' edge.

        Returns
        -------
        An array of distances from the observer to each ring segment as an
        Astropy `Quantity` object.
        """
        arc = self._arcs[name]
        args = dict(time=time,
                    observer=observer,
                    step=step,
                    edge=edge,
                    pericenter=False,
                    minlon=arc.min_longitude,
                    dminlon_dt=arc.dlon_dt,
                    delta_lon=arc.delta_longitude)
        return self._get_distances(**args)

    @property
    def name(self) -> str:
        return self._name

    @property
    def semimajor_axis(self) -> u.Quantity:
        return self._semimajor_axis

    @property
    def width(self) -> u.Quantity:
        return self._width

    @property
    def eccentricity(self) -> float:
        return self._eccentricity

    @property
    def inclination(self) -> u.Quantity:
        return self._inclination

    @property
    def longitude_of_periapsis(self) -> u.Quantity:
        return self._longitude_of_periapsis

    @property
    def longitude_of_ascending_node(self) -> u.Quantity:
        return self._longitude_of_ascending_node

    @property
    def dlon_dt(self) -> u.Quantity:
        return self._dlon_dt

    @property
    def dnode_dt(self) -> u.Quantity:
        return self._dnode_dt

    @property
    def arcs(self) -> dict[str, dict]:
        return self._arcs


class Planet:

    def __init__(self, name: str):
        self._name = name
        self._code = spiceypy.bodn2c(name)
        self._prograde = determine_prograde_rotation(name)
        self._fixref = f'IAU_{name}'
        self._rings = self._get_rings()

    def get_ra(self,
               time: Time,
               observer: str,
               string: bool = False) -> str | Angle:
        """
        Get astrometric right ascension for a given UTC time and
        observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        string : bool
            If true, return the RA as an hms string with 4 decimal places of
            precision.

        Returns
        -------
        The RA as either a string or Astropy `Angle` object.
        """
        et = get_et(time)
        ra, _ = get_sky_coordinates(self._name, et, observer)
        ra = Angle(ra * angle_unit)
        if string:
            return ra.to_string(unit=u.hour, precision=4)
        else:
            return ra

    def get_dec(self,
                time: Time,
                observer: str,
                string: bool = False) -> str | Angle:
        """
        Get astrometric declination for a given UTC time and
        observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        string : bool
            If true, return the declination as an hms string with 3 decimal
            places of precision.

        Returns
        -------
        The declination as either a string or Astropy `Angle` object.
        """
        et = get_et(time)
        _, dec = get_sky_coordinates(self._name, et, observer)
        dec = Angle(dec * angle_unit)
        if string:
            return dec.to_string(unit=u.degree, precision=3)
        else:
            return dec

    def get_skycoord(self,
                     time: Time,
                     observer: str) -> SkyCoord:
        """
        Get astrometric right ascension and declination for a given UTC time
        and observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The sky coordinates as an Astropy `SkyCoord` object.
        """
        et = get_et(time)
        ra, dec = get_sky_coordinates(self._name, et, observer)
        return SkyCoord(ra=ra*angle_unit, dec=dec*angle_unit)

    def get_offset(self,
                   time: Time,
                   observer: str,
                   refcoord: SkyCoord = None) -> tuple[Angle, Angle]:
        """
        Get spherical offset between this body and a reference body for a given
        UTC time and observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        refcoord : SkyCoord
            The other body's coordinates.

        Returns
        -------
        The spherical offset between this body and a reference body as Astropy
        `Angle` objects in units of arcseconds. If `refcoord` is not provided,
        then the returned offset will be (0, 0).
        """
        if refcoord is None:
            return Angle('0d').to(u.arcsec), Angle('0d').to(u.arcsec)
        coord = self.get_skycoord(time, observer)
        offset = refcoord.spherical_offsets_to(coord)
        return offset[0].to(u.arcsec), offset[1].to(u.arcsec)

    def get_sub_observer_latitude(self,
                                  time: Time,
                                  observer: str) -> Latitude:
        """
        Get sub-observer geodetic latitude for a given UTC time and
        observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The sub-observer geodetic latitude as an Astropy `Latitude` object.
        """
        et = get_et(time)
        lat, _ = get_sub_observer_latlon(self._name, et, observer)
        return Latitude(lat * angle_unit).to(u.degree)

    def get_sub_observer_longitude(self,
                                   time: Time,
                                   observer: str) -> Longitude:
        """
        Get sub-observer geodetic longitude for a given UTC time and
        observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The sub-observer geodetic latitude as an Astropy `Longitude` object.
        """
        et = get_et(time)
        _, lon = get_sub_observer_latlon(self._name, et, observer)
        return Longitude(lon * angle_unit).to(u.degree)

    def get_subsolar_latitude(self,
                              time: Time,
                              observer: str) -> Latitude:
        """
        Get sub-solar geodetic latitude for a given UTC time and
        observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The sub-observer geodetic latitude as an Astropy `Latitude` object.
        """
        et = get_et(time)
        lat, _ = get_sub_solar_latlon(self._name, et, observer)
        return Latitude(lat * angle_unit).to(u.degree)

    def get_subsolar_longitude(self,
                               time: Time,
                               observer: str) -> Longitude:
        """
        Get sub-solar geodetic longitude for a given UTC time and
        observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The sub-observer geodetic latitude as an Astropy `Longitude` object.
        """
        et = get_et(time)
        _, lon = get_sub_solar_latlon(self._name, et, observer)
        return Longitude(lon * angle_unit).to(u.degree)

    def get_phase_angle(self,
                        time: Time,
                        observer: str) -> Angle:
        """
        Get the phase angle at the sub-observer point for a given UTC time and
        observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The phase angle at the sub-observer point as an Astropy `Angle` object.
        """
        et = get_et(time)
        phase = get_sub_observer_phase_angle(self._name, et, observer)
        return Angle(phase * angle_unit).to(u.degree)

    def get_distance(self,
                     time: Time,
                     observer: str) -> u.Quantity:
        """
        Get the distance to the object for a given UTC time and
        observer/observatory.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The distance to the object as an Astropy `Quantity` object.
        """
        et = get_et(time)
        vec = get_state(self._name, et, observer)
        return spiceypy.vnorm(vec) * length_unit

    def get_ring_subsolar_latitude(self,
                                   time: Time,
                                   observer: str) -> Latitude:
        """
        This is the same as the sub-solar latitude.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The sub-solar latitude as an Astropy `Latitude` object.
        """
        return self.get_subsolar_latitude(time, observer)

    def get_ring_subsolar_latitude_range(self,
                                         time: Time,
                                         observer: str
                                         ) -> tuple[Latitude, Latitude]:
        """
        Get the range in subsolar latitude due to the apparent angular size of
        the Sun.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------

        """
        _, radii = spiceypy.bodvcd(10, 'RADII', 3)
        et = get_et(time)
        vec = get_state(self._name, et, 'Sun')
        angular_radius = Angle(
            radii[0] / spiceypy.vnorm(vec), unit=u.rad).to(u.degree)
        lat0 = self.get_ring_subsolar_latitude(time, observer)
        return Latitude(lat0 - angular_radius), Latitude(lat0 + angular_radius)

    def get_ring_plane_opening_angle(self,
                                     time: Time,
                                     observer: str) -> Angle:
        """
        This is the same as the sub-observer latitude.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The ring plane opening angle as an Astropy `Angle` object.
        """
        return Angle(self.get_sub_observer_latitude(time, observer))

    def determine_if_rings_illuminated(self,
                                       time: Time,
                                       observer: str) -> bool:
        """
        Determine if the rings appear illuminated to the observer. Assumes the
        ring plane is at 0° latitude.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        A bool indicating if the rings are illuminated to the observer.
        """
        lat0, lat1 = self.get_ring_subsolar_latitude_range(time, observer)
        sun_on_one_side_of_plane = np.sign(lat0) == np.sign(lat1)

        # both sides of rings illuminated if Sun on both sides of ring plane
        if not sun_on_one_side_of_plane:
            return True

        # otherwise, the rings are illuminated only if the sub-solar latitude
        #  is the same sign for both the Sun and the observer
        lat_obs = self.get_sub_observer_longitude(time, observer)
        lat_sun = self.get_subsolar_latitude(time, observer)
        sun_on_same_side_as_observer = np.sign(lat_obs) == np.sign(lat_sun)
        return sun_on_same_side_as_observer

    def get_ring_center_phase_angle(self,
                                    time: Time,
                                    observer: str) -> Angle:
        """
        This is the same as the body phase angle.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The body phase angle as an Astropy `Angle` object.
        """
        return self.get_phase_angle(time, observer).to(u.degree)

    def get_ascending_node_longitude(self,
                                     time: Time) -> Longitude:
        """
        Get geodetic longitude of the ring plane ascending node (by definition,
        the ascending node of the primary body's equator on the J2000 equator).

        Parameters
        ----------
        time : Time
            UTC at the time of observation.

        Returns
        -------
        The geodetic longitude of the ring plane ascending node as an Astropy
        `Longitude` object.
        """
        et = get_et(time)
        ra, _, _, _ = spiceypy.bodeul(self._code, et)
        if not self._prograde:
            ra = (ra + spiceypy.pi()) % spiceypy.twopi()
        offset = Angle(90, unit='deg').rad
        node_vector = np.array([np.cos(ra + offset), np.sin(ra + offset), 0])
        radii = get_radii(self._name)
        re = float(radii[0])
        rp = float(radii[2])
        f = (re - rp)/re
        node_vector *= re
        rotation = spiceypy.pxform('J2000', self._fixref, et)
        node_vector = spiceypy.mxv(rotation, node_vector)
        reflon, _, _ = spiceypy.recgeo(node_vector, re, f)
        return Longitude(reflon * u.rad).to(u.degree)

    def _get_ring_sub_longitude(self,
                                time: Time,
                                observer: str,
                                desc: str) -> Longitude:
        """
        Calculate the sub-solar or sub-observer geodetic longitude relative to
        the ring plane ascending node.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        desc: str
            Whether you want 'observer' for sub-observer coordinates or 'solar'
            for sub-solar coordinates.

        Returns
        -------
        The sub-solar or sub-observer longitude geodetic longitude relative to
        the ring plane ascending node as an Astropy `Longitude` object.
        """
        ascnode = self.get_ascending_node_longitude(time)
        if desc == 'solar':
            sublon = self.get_subsolar_longitude(time, observer)
        else:
            sublon = self.get_sub_observer_longitude(time, observer)
        longitude = sublon - ascnode
        if longitude < 0 * u.degree:
            longitude += 360 * u.degree
        return Longitude(longitude).to(u.degree)

    def get_ring_subsolar_longitude(self,
                                    time: Time,
                                    observer: str) -> Longitude:
        """
        Get sub-solar longitude of the rings relative to the plane ascending
        node.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The sub-solar longitude of the rings relative to the plane ascending
        node as an Astropy `Longitude` object.
        """
        return self._get_ring_sub_longitude(time, observer, 'solar')

    def get_ring_sub_observer_longitude(self,
                                        time: Time,
                                        observer: str) -> Longitude:
        """
        Get sub-observer longitude of the rings relative to the plane ascending
        node.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The sub-observer longitude of the rings relative to the plane ascending
        node as an Astropy `Longitude` object.
        """
        return self._get_ring_sub_longitude(time, observer, 'observer')

    def get_sun_target_distance(self,
                                time: Time,
                                unit: u.Unit = u.au) -> u.Quantity:
        """
        Get the distance between the target and the Sun.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        unit : u.Unit
            Unit of output distance. Default is "au".

        Returns
        -------
        The distance between the Sun and the target as an Astropy `Quantity`
        object.
        """
        return self.get_distance(time, 'Sun').to(unit)

    def get_observer_target_distance(self,
                                     time: Time,
                                     observer: str) -> u.Quantity:
        """
        Get the distance between the target and the observer.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The distance between the observer and the target as an Astropy
        `Quantity` object.
        """
        return self.get_distance(time, observer)

    def get_light_travel_time(self,
                              time: Time,
                              observer: str) -> u.Quantity:
        """
        Get the light travel time between the target and the observer.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The light travel time between the observer and the target as an Astropy
        `Quantity` object.
        """
        et = get_et(time)
        return get_light_time(self._name, et, observer) * u.s

    def get_latlon_sky_coordinates(self,
                                   time: Time,
                                   observer: str,
                                   latitude: u.Quantity,
                                   longitude: u.Quantity,
                                   ) -> list[SkyCoord] or SkyCoord:
        """
        Get the right ascension and declination of of a latitude/longitude
        point on the body's surface for a given observer and observation time.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        latitude : u.Quantity
            A latitude.
        longitude : u.Quantity
            A longitude.

        Returns
        -------
        The right ascension and declination of the latitude/longitude point
        as an Astropy `SkyCoord` object.'
        """
        et = get_et(time)
        latitude = latitude.to(u.rad).value
        longitude = longitude.to(u.rad).value
        ra, dec = get_latlon_radec(
            self._name, et, observer, latitude, longitude)
        return SkyCoord(ra=Angle(ra*angle_unit), dec=Angle(dec*angle_unit))

    def get_longitude_line_coordinates(self,
                                       time: Time,
                                       observer: str,
                                       longitude: Angle,
                                       dlat: Angle = Angle(1, unit='deg'),
                                       ) -> list[SkyCoord]:
        """
        Get the right ascension and declination of a longitude line.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        longitude : Angle
            The longitude.
        dlat : Angle
            The latitudinal resolution of the longitude line. Default is 1
            degree.

        Returns
        -------
        The right ascension and declination of the longitude line as a list of
        Astropy `SkyCoord` objects.'
        """
        latitudes = np.linspace(-90, 90, int(180/dlat.value) + 1) * u.degree
        coords = []
        for latitude in latitudes:
            coords.append(self.get_latlon_sky_coordinates(time, observer,
                                                          latitude, longitude))
        return coords

    def get_latitude_line_coordinates(self,
                                      time: Time,
                                      observer: str,
                                      latitude: Angle,
                                      dlon: Angle = Angle(1, unit='deg'),
                                      ) -> list[SkyCoord]:
        """
        Get the right ascension and declination of a latitude line.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        latitude : Angle
            The latitude.
        dlon : Angle
            The longitudinal resolution of the latitude line. Default is 1
            degree.

        Returns
        -------
        The right ascension and declination of the latitude line as a list of
        Astropy `SkyCoord` objects.'
        """
        longitudes = np.linspace(-180, 180, int(360/dlon.value) + 1) * u.degree
        coords = []
        for longitude in longitudes:
            coords.append(self.get_latlon_sky_coordinates(time, observer,
                                                          latitude, longitude))
        return coords

    @staticmethod
    def _convert_to_sky_coords(time: Time,
                               radec_func,
                               **kwargs) -> list[SkyCoord]:
        """
        Convenience function to calculate RA/Dec with one of the SPICE
        functions and convert to SkyCoord objects.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        radec_func
            The function which returns RA and Dec in units of [rad].
        **kwargs
            The arguments to `radec_func` with the exception of `et` which is
            calculated internally.

        Returns
        -------
        The coordinates as a list of Astropy `SkyCoord` objects.
        """
        et = get_et(time)
        ras, decs = radec_func(et=et, **kwargs)
        coordinates = []
        for ra, dec in zip(ras, decs):
            coordinates.append(
                SkyCoord(ra=ra * angle_unit, dec=dec * angle_unit))
        return coordinates

    def get_limb_sky_coordinates(self,
                                 time: Time,
                                 observer: str,
                                 resolution: int = 360) -> list[SkyCoord]:
        """
        Get the right ascension and declination of the body's limb for a given
        observer and observation time.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        resolution : int
            The angular resolution of the limb points. The default is 360,
            which calculates the limb of a spherical body in 1-degree
            increments.

        Returns
        -------
        The right ascension and declination of the body's limb as a list of
        Astropy `SkyCoord` objects.'
        """
        kwargs = dict(time=time,
                      radec_func=get_limb_radec,
                      target=self._name,
                      obs=observer,
                      resolution=resolution)
        return self._convert_to_sky_coords(**kwargs)

    def get_terminator_sky_coordinates(
            self,
            time: Time,
            observer: str,
            resolution: int = 360,
            kind: str = 'umbral') -> list[SkyCoord]:
        """
        Get the right ascension and declination of the body's terminator for a
        given observer and observation time.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        resolution : int
            The angular resolution of the limb points. The default is 360,
            which calculates the limb of a spherical body in 1-degree
            increments.
        kind : str
            Type of terminator. Default is 'umbral' but can also be
            'penumbral'.

        Returns
        -------
        The right ascension and declination of the body's limb as a list of
        Astropy `SkyCoord` objects.'
        """
        kwargs = dict(time=time,
                      radec_func=get_terminator_radec,
                      target=self._name,
                      obs=observer,
                      resolution=resolution,
                      kind=kind)
        return self._convert_to_sky_coords(**kwargs)

    def get_dayside_sky_coordinates(self,
                                    time: Time,
                                    observer: str,
                                    resolution: int = 360) -> list[SkyCoord]:
        """
        Get the right ascension and declination for the dayside of the body's
        apparent disk.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        resolution : int
            The angular resolution of the limb points. The default is 360,
            which calculates the limb of a spherical body in 1-degree
            increments.

        Returns
        -------
        The right ascension and declination of the body's apparent dayside disk
        as a list of Astropy `SkyCoord` objects.'
        """
        kwargs = dict(time=time,
                      radec_func=get_dayside_radec,
                      target=self._name,
                      obs=observer,
                      resolution=resolution)
        return self._convert_to_sky_coords(**kwargs)

    def get_shadow_intersection_sky_coordinates(
            self,
            shadow_casting_body: str,
            time: Time,
            observer: str,
            resolution: int = 360) -> list[SkyCoord]:
        """
        Get the right ascension and declination for the intersection of the
        shadow cast by one body with the apparent disk of another.

        Parameters
        ----------
        shadow_casting_body : str
            The name of the body casting the shadow.
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        resolution : int
            The angular resolution of the limb points. The default is 360,
            which calculates the limb of a spherical body in 1-degree
            increments.

        Returns
        -------
        The right ascension and declination of the body's apparent dayside disk
        as a list of Astropy `SkyCoord` objects.'
        """
        kwargs = dict(time=time,
                      radec_func=get_shadow_intercept,
                      target1=self._name,
                      target2=shadow_casting_body,
                      obs=observer,
                      resolution=resolution)
        return self._convert_to_sky_coords(**kwargs)

    def get_angular_radius(self,
                           time: Time,
                           observer: str,) -> Angle:
        """
        Get the body's angular radius for a given observer and observation
        time. Calculated using the object's mean equatorial radius.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The angular radius of the body as an Astropy `Angle` object.
        """
        distance = self.get_distance(time, observer)
        re = np.mean(get_radii(self._name)[:2])
        return Angle(np.arctan(re * length_unit / distance)).to(u.arcsec)

    @staticmethod
    def _parse_unit(value, unit) -> u.Quantity:
        """
        Replace NaN values with zero.

        Parameters
        ----------
        value : int or float
            The quantity's value.
        unit : u.Unit
            The quantity's unit.

        Returns
        -------
        The parsed quantity as an Astropy `Quantity` object.
        """
        if np.isnan(value):
            return 0.0 * unit
        else:
            return value * unit

    def _add_units_to_ring_props(self, properties: dict) -> dict:
        """
        Add units to ring orbital elements.
        
        Parameters
        ----------
        properties : dict
            A dictionary of the ring orbital elements.

        Returns
        -------
        A dictionary like the input dictionary but now with units.
        """
        lengths = ['semimajor_axis', 'width']
        for i in lengths:
            properties[i] = self._parse_unit(properties[i], u.km)

        angles = ['inclination', 'longitude_of_periapsis',
                  'longitude_of_ascending_node', 'longitude_offset']
        for i in angles:
            properties[i] = self._parse_unit(properties[i], u.degree)

        angular_rates = ['dlon_dt', 'dnode_dt']
        for i in angular_rates:
            properties[i] = self._parse_unit(
                properties[i], u.degree/u.day).to(u.degree/u.s)

        properties['epoch'] = Time(properties['epoch'], scale='utc')

        if np.isnan(properties['eccentricity']):
            properties['eccentricity'] = 0.0

        for arc in properties['arcs']:
            arc['min_longitude'] = self._parse_unit(arc['min_longitude'],
                                                    u.degree)
            arc['max_longitude'] = self._parse_unit(arc['max_longitude'],
                                                    u.degree)
            arc['dlon_dt'] = self._parse_unit(arc['dlon_dt'], u.degree/u.day)

        return properties

    def _make_ring(self, properties: dict) -> Ring:
        """
        Convenience function to create `Ring` objects.
        
        Parameters
        ----------
        properties : dict
            A dictionary of the ring orbital elements.

        Returns
        -------
        The ring as a `Ring` object.
        """
        ascending_node = self.get_ascending_node_longitude(properties['epoch'])
        offset = properties['longitude_offset']
        dt = get_et(properties['epoch']) * u.s
        parent_rate = get_rotation_rate(self._name) * angle_rate_unit
        if self._prograde:
            scale = -1
        else:
            scale = 1
        parent_rotation_offset = parent_rate * dt * scale

        if properties['longitude_of_periapsis'] != 0.0 * u.degree:
            properties['longitude_of_periapsis'] -= ascending_node
            properties['longitude_of_periapsis'] += offset
            properties['longitude_of_periapsis'] -= properties['dlon_dt'] * dt
            properties['longitude_of_periapsis'] -= parent_rotation_offset
            properties['longitude_of_periapsis'] %= 360 * u.degree
        if properties['longitude_of_ascending_node'] != 0.0 * u.degree:
            properties['longitude_of_ascending_node'] -= ascending_node
            properties['longitude_of_ascending_node'] += offset
            properties['longitude_of_ascending_node'] -= (
                    properties['dnode_dt'] * dt)
            properties['longitude_of_ascending_node'] -= parent_rotation_offset
            properties['longitude_of_ascending_node'] %= 360 * u.degree
        arcs = {}
        for arc in properties['arcs']:
            arc['min_longitude'] -= ascending_node
            arc['max_longitude'] -= ascending_node
            arc['min_longitude'] += offset
            arc['max_longitude'] += offset
            arc['min_longitude'] -= arc['dlon_dt'] * dt
            arc['max_longitude'] -= arc['dlon_dt'] * dt
            arc['max_longitude'] -= parent_rotation_offset
            arc['min_longitude'] -= parent_rotation_offset
            arc['min_longitude'] %= 360 * u.degree
            arc['max_longitude'] %= 360 * u.degree
            arc = Arc(name=arc['name'],
                      min_longitude=arc['min_longitude'],
                      max_longitude=arc['max_longitude'],
                      dlon_dt=arc['dlon_dt'])
            arcs[arc.name] = arc
        properties['arcs'] = arcs

        # remove unneeded properties
        properties.pop('epoch')
        properties.pop('longitude_offset')

        return Ring(**properties)

    def _get_rings(self) -> dict[str, Ring]:
        """
        A convenience function to generate all rings for a given parent body.
        
        Returns
        -------
        A dictionary of the rings as `Ring` objects.
        """
        properties = pd.read_csv(Path(_project_directory, 'anc', 'rings.dat'))
        ring_arcs = pd.read_csv(Path(_project_directory, 'anc', 'arcs.dat'))
        rings = {}
        for row in properties.iloc:
            if row['planet'] == self._name:
                ring_properties = row.to_dict()
                if np.isnan(ring_properties['semimajor_axis']):
                    a = (ring_properties['outer_radius'] +
                         ring_properties['inner_radius']) / 2
                    width = (ring_properties['outer_radius'] -
                             ring_properties['inner_radius'])
                    ring_properties['semimajor_axis'] = a
                    ring_properties['width'] = width
                else:
                    ring_properties['width'] = 0.0

                # remove unneeded properties
                ring_properties.pop('inner_radius')
                ring_properties.pop('outer_radius')

                # add additional properties
                ring_properties['planet_prograde'] = self._prograde

                arcs = []
                if self._name in ring_arcs['planet'].tolist():
                    names = ring_arcs['ring_name'].tolist()
                    if ring_properties['name'] in names:
                        for arc_prop in ring_arcs.iloc:
                            props = arc_prop.to_dict()
                            props.pop('planet')
                            props.pop('ring_name')
                            arcs.append(props)
                ring_properties['arcs'] = arcs

                # add units
                ring_properties = self._add_units_to_ring_props(
                    ring_properties)
                rings[f"{ring_properties['name']}"] = (
                    self._make_ring(ring_properties))

        return rings

    def get_angular_offset_between_times(self,
                                         relative_time: Time,
                                         current_time: Time,
                                         observer: str) -> tuple[Angle, Angle]:
        """
        Calculate the RA/Dec offset necessary to shift this object's sky
        coordinates at  `current_time` to `relative_time`. This will let you
        plot an object relative to another object at a specific time centered
        in the field of view. Useful for time-series plots showing things like
        transits of one body across another.

        Parameters
        ----------
        relative_time : Time
            The time from which you want the offsets calculated.
        current_time : Time
            The current time.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        The offsets in RA and Dec as Astropy `Angle` objects.
        """
        relative_coord = self.get_skycoord(relative_time, observer)
        current_coord = self.get_skycoord(current_time, observer)
        dra, ddec = current_coord.spherical_offsets_to(relative_coord)
        return dra, ddec

    def parse_fov(self,
                  time: Time,
                  observer: str,
                  fov: float or int or Angle or u.Quantity) -> Angle:
        """
        Convert distances at the target body to angular sizes on the sky. This
        could be an `Angle` or angular quantity with units like [rad], [deg] or
        [arcsec], a length with units like [km] or an integer or float which is
        interpreted to be multiples of the planet's equatorial radius.

        This is useful for defining a WCS field of view. If you want an FOV
        with a width of 4 planetary radii, then calling
        `parse_fov(time, observer, 4)` will return the angle corresponding to
        that distance on the sky.

        Parameters
        ----------
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        fov : u.Quantity
            The width of the field-of-view.

        Returns
        -------
        The angular field-of-view as an Astropy `Angle` object.
        """
        distance = self.get_distance(time, observer)
        if isinstance(fov, int) or isinstance(fov, float):
            fov = self.get_angular_radius(time, observer) * fov
        else:
            try:
                fov = fov.to(length_unit)
                fov = Angle(np.arctan(fov/distance))
            except u.UnitConversionError:
                fov = Angle(fov)
        return fov.to(u.arcsec)

    @staticmethod
    def _parse_list_or_str(thing: list[str] or str) -> list[str]:
        """
        Determine if a thing is a list of strings or a string. Helps account
        for people who supply single ring names as a string rather than a list
        with one item.
        
        Parameters
        ----------
        thing : list[str] or str
            The thing in question.

        Returns
        -------
        A list of string(s).
        """
        if isinstance(thing, str):
            return [thing]
        else:
            return thing

    @staticmethod
    def _get_ring_pericenters(
            ring: Ring,
            time: Time,
            observer: str) -> tuple[SkyCoord, u.Quantity] or tuple[None, None]:
        """
        Convenience function to get the coordinate of and distance to the 
        ring's pericenter. Applies only to some of Uranus's rings.
         
        Parameters
        ----------
        ring : Ring
            A `Ring` object.
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".

        Returns
        -------
        A tuple (`SkyCoord`, `Quantity`) of the coordinate and the distance.
        """
        args = dict(time=time, observer=observer, pericenter=True)
        pericenter_coord = ring.get_sky_coordinates(**args)
        pericenter_distance = ring.get_distances(**args)
        if pericenter_coord is None:
            return None, None
        else:
            return pericenter_coord[0], pericenter_distance[0]

    def _parse_and_plot_ring(self,
                             axis: WCSAxes,
                             time: Time,
                             observer: str,
                             side: str,
                             dra: Angle,
                             ddec: Angle,
                             rings: list[str] or None,
                             arcs: list[str] or None) -> None:
        """
        Since the rings have to be plotted twice (front half and back half),
        this combines all of the plotting logic into a single function.

        Parameters
        ----------
        axis : WCSAxes
            The axis in which to plot the ring.
        time : Time
            UTC at the time of observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        side : str
            The side of the rings to plot, either 'front' or 'back'. You'll
            want to plot the back first, then plot the planet, then the front.
        dra : Angle
            Optional angular offset in right ascension.
        ddec : Angle
            Optional angular offset in declination.
        rings : list[str] or None
            The list of rings to plot.
        arcs : list[str] or None
            The list of arcs to plot.

        Returns
        -------
        None.
        """
        body_distance = self.get_distance(time, observer)
        edges = ['inner', 'outer']

        for name in rings:
            distance_range = None
            if name in self._rings.keys():
                ring = self._rings[name]
                args = dict(time=time, observer=observer)
                if ring.width != 0.0:
                    coords = [
                        ring.get_sky_coordinates(**args, edge=edge,
                                                 boundaries=True)
                        for edge in edges]
                    distances = [
                        ring.get_distances(**args, edge=edge)
                        for edge in edges]
                    eclipsed = [
                        ring.get_eclipsed(time, edge=edge)
                        for edge in edges]
                else:
                    coords = [ring.get_sky_coordinates(
                        **args, boundaries=True)]
                    distances = [ring.get_distances(time, observer)]
                    eclipsed = [ring.get_eclipsed(time)]
                if distance_range is None:
                    distance_range = (distances[0].max() -
                                      distances[0].min())
                plot_ring(axis, coords, distances, eclipsed, side, dra,
                          ddec)
                if arcs is None:
                    use_arcs = list(ring.arcs.keys())
                else:
                    use_arcs = self._parse_list_or_str(arcs)
                for arc in use_arcs:
                    arc_distances = ring.get_arc_distances(
                        arc, time, observer)
                    mean_distance = np.mean(arc_distances)
                    # check if near ansa, if so, skip for backside and only
                    # draw on the front side to avoid ring overlap issues
                    if side == 'back':
                        if ((np.abs(mean_distance - body_distance)
                             < 0.5 * distance_range) or
                                (mean_distance < body_distance)):
                            coords = ring.get_arc_sky_coordinates(
                                arc, time, observer, boundaries=True)
                            plot_arc(axis, coords, dra, ddec)
                    elif side == 'front':
                        if (np.abs(mean_distance - body_distance)
                                < 0.5 * distance_range):
                            pass
                        elif mean_distance > body_distance:
                            coords = ring.get_arc_sky_coordinates(
                                arc, time, observer, boundaries=True)
                            plot_arc(axis, coords, dra, ddec)

        # pericenter parkers
        for name in rings:
            if name in self._rings.keys():
                ring = self._rings[name]
                coord, distance = (
                    self._get_ring_pericenters(ring, time, observer))
                place_ring_pericenter_markers(axis, coord, distance,
                                              body_distance, side, dra,
                                              ddec)

    def draw(self,
             axis: WCSAxes,
             time: Time,
             observer: str,
             rings: list[str] or str = None,
             arcs: list[str] or str = None,
             dra: Angle = Angle(0, unit='deg'),
             ddec: Angle = Angle(0, unit='deg')) -> None:
        """
        Draw the planet (and optionally its rings and/or ring arcs) for a given
        time and observer.

        Parameters
        ----------
        axis : WCSAxes
            The axis in which to draw the planet. Must be an axis with an
            Astropy `WCS` projection.
        time : Time
            The UTC time of the observation.
        observer : str
            The observer or observatory. Could be a Solar System body like
            "Ganymede", an Earth-based observatory like "Keck" or a spacecraft
            like "Juno".
        rings : list[str] or None
            A list of rings to plot. If `None`, all rings will be shown. If an
            empty list, no rings will be shown.
        arcs : list[str] or None
            A list of arcs to plot. If `None`, all arcs will be shown. If an
            empty list, no arcs will be shown.
        dra : Angle
            Optional angular offset in right ascension.
        ddec : Angle
            Optional angular offset in declination.

        Returns
        -------
        None.
        """

        # if not the primary body, plot any shadow from the primary body
        occultation_code = 0
        parent_name = None
        et = get_et(time)
        if str(self._code)[1:3] != '99':
            parent_name = spiceypy.bodc2n(int(f'{str(self._code)[0]}99'))
            occultation_code = determine_eclipse_or_transit(
                self._name, parent_name, et, observer)

        # rings back side
        if rings is None:
            rings = list(self._rings.keys())
        else:
            rings = self._parse_list_or_str(rings)
        self._parse_and_plot_ring(axis, time, observer, side='back', dra=dra,
                                  ddec=ddec, rings=rings, arcs=arcs)

        # plot the disk
        coords = self.get_limb_sky_coordinates(time, observer)
        facecolor = 'white'
        if occultation_code == 3:
            facecolor = _default_night_patch_kwargs()['facecolor']
        plot_disk(axis, coords, dra, ddec, facecolor=facecolor)

        # plot night side
        limb_coords = self.get_limb_sky_coordinates(time, observer)
        terminator_coords = self.get_terminator_sky_coordinates(time, observer)
        plot_nightside(axis, limb_coords, terminator_coords, dra, ddec)

        if (occultation_code == 1) & (parent_name is not None):
            disk_coords = self.get_dayside_sky_coordinates(time, observer)
            shadow_coords = self.get_shadow_intersection_sky_coordinates(
                parent_name, time, observer)
            plot_primary_shadow(axis, disk_coords, shadow_coords, dra, ddec)

        # plot latitudes
        for latitude in np.arange(-60, 60+30, 30) * u.degree:
            coords = self.get_latitude_line_coordinates(
                time, observer, latitude)
            plot_latlon(axis, coords, dra, ddec)

        # plot longitudes
        linewidth = plt.rcParams['lines.linewidth'] / 2
        for longitude in np.arange(30, 360, 30) * u.degree:
            coords = self.get_longitude_line_coordinates(
                time, observer, longitude)
            plot_latlon(axis, coords, dra, ddec)
        longitude = 0 * u.degree
        coords = self.get_longitude_line_coordinates(time, observer, longitude)
        plot_latlon(axis, coords, dra, ddec, edgecolor='black',
                    linewidth=linewidth)

        # plot the limb
        coords = self.get_limb_sky_coordinates(time, observer)
        plot_limb(axis, coords, dra, ddec)

        # rings/arcs front side
        self._parse_and_plot_ring(axis, time, observer, side='front', dra=dra,
                                  ddec=ddec, rings=rings, arcs=arcs)

    @property
    def name(self) -> str:
        return self._name

    @property
    def code(self) -> int:
        return self._code

    @property
    def prograde(self) -> bool:
        return self._prograde

    @property
    def rings(self) -> dict[str, Ring]:
        return self._rings
