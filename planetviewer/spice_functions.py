import astropy.units as u
import numpy as np
import spiceypy as spice
from spiceypy.utils.exceptions import NotFoundError

# Set fundamental SPICE units for length, time and angle.
length_unit = u.km
time_unit = u.s
angle_unit = u.rad
angle_rate_unit = u.degree / u.day


# Set common SPICE parameters
abcorr = 'LT'
ref = 'J2000'
subpoint_method = 'ELLIPSOID'
surface_method = 'NEAR POINT/ELLIPSOID'
corloc = 'ELLIPSOID LIMB'
refvec = spice.vpack(0, 0, 1)

_planets = ['Mercury',
            'Venus',
            'Earth',
            'Mars',
            'Jupiter',
            'Saturn',
            'Uranus',
            'Neptune',
            'Pluto']


def determine_prograde_rotation(target: str) -> bool:
    """
    Determine if an object rotates prograde (most of the planets) or retrograde
    (Venus and Uranus).

    Parameters
    ----------
    target : str
        Name of target body.

    Returns
    -------
    bool
        A bool indicating whether or not the body's rotation is prograde or
        not.
    """
    code = spice.bodn2c(target)
    pm = spice.bodvcd(code, 'PM', 3)[1]
    if pm[1] < 0:
        return False
    else:
        return True


def get_rotation_rate(target: str) -> float:
    """
    Get parent body angular rotation rate.

    Returns
    -------
    float
        Parent body angular rotation rate in units of [deg/day].
    """
    code = spice.bodn2c(target)
    rate = float(spice.bodvcd(code, 'PM', 3)[1][1])
    return rate


def get_radii(target: str) -> np.ndarray:
    """
    Get triaxial ellipsoid radii for the target object.

    Parameters
    ----------
    target : str
        Name of target body.

    Returns
    -------
    np.ndarray
        An array of the triaxial ellipsoid axes in [km].
    """
    code = spice.bodn2c(target)
    _, radii = spice.bodvcd(code, 'RADII', 3)
    return radii


def get_equatorial_radius(target: str,
                          lon: float = 0.0) -> float:
    """
    Get equatorial ellipsoid radius for the target object. If the x and y-axis
    radii are not equal, then it has to be the ellipse radius, not just the
    x-axis radius! Learned that one the hard way...

    Parameters
    ----------
    target : str
        Name of target body.
    lon : float
        Planetographic longitude of target body.

    Returns
    -------
    float
        An array of the triaxial ellipsoid axes in [km].
    """
    radii = get_radii(target)
    re = np.sqrt(
        (radii[0] * np.cos(lon)) ** 2 + (radii[1] * np.sin(lon)) ** 2)
    return float(re)


def get_flattening_coefficient(target: str,
                               lon: float = 0.0) -> float:
    """
    Get ellipsoid flattening coefficient.

    Parameters
    ----------
    target : str
        Name of target body.
    lon : float
        Planetographic longitude of target body.

    Returns
    -------
    float
        The flattening coefficient (0 = spherical, 1 = pancake).
    """
    rp = get_radii(target)[2]
    re = get_equatorial_radius(target, lon)
    return float((re - rp) / re)


def get_sky_coordinates(target: str,
                        et: float,
                        obs: str) -> tuple[float, float]:
    """
    Get the RA and Dec of an object as viewed from a specified observatory at
    a given epoch.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    tuple[float, float]
        The object RA and Dec in [rad].
    """
    ptarg, _ = spice.spkpos(target, et, ref, abcorr, obs)
    _, ra, dec = spice.recrad(ptarg)
    return ra, dec


def get_light_time(target: str,
                   et: float,
                   obs: str,
                   direction: str = '<-') -> float:
    """
    Get the one-way light travel time between an object and an observer at a
    given epoch.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    direction : str
        Direction of light travel in units of [deg/day].

    Returns
    -------
    tuple[float, float]
        The apparent epoch at the target and the light travel time in [s].
    """
    code = spice.bods2c(target)
    _, lt = spice.ltime(et, code, direction, spice.bods2c(obs))
    return lt


def get_apparent_epoch(target: str,
                       et: float,
                       obs: str,
                       direction: str = '<-') -> float:
    """
    Get the apparent epoch at a target for an observer at a given epoch.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    direction : str
        Direction of light travel in units of [deg/day].

    Returns
    -------
    tuple[float, float]
        The apparent epoch at the target and the light travel time in [s].
    """
    code = spice.bods2c(target)
    et_target, _ = spice.ltime(et, code, direction, spice.bods2c(obs))
    return et_target


def get_target_frame(target: str) -> str:
    """
    Get IAU target frame name except for Earth and the Moon, which return the
    high-precision ITRF93 and MEAN_ME frames.
    
    Parameters
    ----------
    target : str
        Name of target body.

    Returns
    -------
    str
        The frame name as a string.
    """
    if target == 'Earth':
        return 'ITRF93'
    elif target == 'Moon':
        return 'MEAN_ME'
    else:
        return f'IAU_{target}'


def _get_subpnt(target: str,
                et: float,
                obs: str,
                desc: str) -> tuple[np.ndarray, float, np.ndarray]:
    """
    Wrapper function to get either sub-solar or sub-observer points.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    desc: str
        Whether you want 'observer' for sub-observer coordinates or 'solar' for
        sub-solar coordinates.

    Returns
    -------
    tuple[np.ndarray, float, np.ndarray]
        The sub-observer point on the target body, the sub-observer point
        epoch and the vector from the observer to sub-observer point.
    """
    kwargs = dict(method=surface_method,
                  target=target,
                  et=et,
                  fixref=get_target_frame(target),
                  abcorr=abcorr,
                  obsrvr=obs)
    if desc == 'observer':
        return spice.subpnt(**kwargs)
    elif desc == 'solar':
        return spice.subslr(**kwargs)
    else:
        raise SystemExit("You have to choose between 'observer' and 'solar'!")


def _get_sub_latlon(target: str,
                    et: float,
                    obs: str,
                    desc: str) -> tuple[float, float]:
    """
    Get sub-observer or sub-solar geodetic latitude and east longitude on the
    target body.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    desc: str
        Whether you want 'observer' for sub-observer coordinates or 'solar' for
        sub-solar coordinates.

    Returns
    -------
    tuple[float, float]
        The sub-observer or sub-solar latitude and longitude on the target body
        in [rad].
    """
    spoint, _, _ = _get_subpnt(target, et, obs, desc)
    lon = np.arctan2(spoint[1], spoint[0])
    re = get_equatorial_radius(target, lon)
    f = get_flattening_coefficient(target, lon)
    lon, lat, _ = spice.recpgr(target, spoint, re, f)
    return lat, lon


def _get_sub_distance(target: str,
                      et: float,
                      obs: str,
                      desc: str) -> float:
    """
    Get distance to the sub-observer point or sub-solar point on target body.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    desc: str
        Whether you want 'observer' for sub-observer coordinates or 'solar' for
        sub-solar coordinates.

    Returns
    -------
    float
        The distance to the sub-observer point or sub-solar point on target
        body in [km].
    """
    _, _, srfvec = _get_subpnt(target, et, obs, desc)
    return spice.vnorm(srfvec)


def get_sub_observer_latlon(target: str,
                            et: float,
                            obs: str) -> tuple[float, float]:
    """
    Get sub-solar geodetic latitude and east longitude on target body's
    surface.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    tuple[float, float]
        The sub-solar latitude and longitude on the target body in [rad].
    """
    return _get_sub_latlon(target, et, obs, 'observer')


def get_sub_observer_distance(target: str,
                              et: float,
                              obs: str) -> float:
    """
    Get distance from observer to sub-observer point on target body's surface.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    float
        The distance between the observer and the sub-observer point on the
        target body's surface in [km].
    """
    return _get_sub_distance(target, et, obs, 'observer')


def get_sub_observer_epoch(target: str,
                           et: float,
                           obs: str) -> float:
    """
    Get the apparent epoch at the sub-observer point on target body's surface.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    float
        The apparent time at the sub-observer point (effectively the epoch with
        the light travel time to the sub-observer point subtracted).
    """
    _, trgepc, _ = _get_subpnt(target, et, obs, 'observer')
    return trgepc


def get_sub_observer_phase_angle(target: str,
                                 et: float,
                                 obs: str) -> float:
    """
    Get the phase angle at the sub-observer point on the target body's surface.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    float
        The phase angle at the sub-observer point in [rad].
    """
    spoint, _, _ = _get_subpnt(target, et, obs, 'observer')
    _, _, phase, _, _ = spice.ilumin(
        subpoint_method, target, et, get_target_frame(target), abcorr, obs,
        spoint)
    return phase


def get_sub_solar_latlon(target: str,
                         et: float,
                         obs: str) -> tuple[float, float]:
    """
    Get sub-solar geodetic latitude and east longitude on the target body's
    surface.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.e
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    tuple[float, float]
        The sub-solar latitude and longitude on the target body in [rad].
    """
    return _get_sub_latlon(target, et, obs, 'solar')


def get_sub_solar_distance(target: str,
                           et: float,
                           obs: str) -> float:
    """
    Get distance from observer to sub-solar point on the target body's surface.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    float
        The distance between the observer and the sub-solar point on the target
        body's surface in [km].
    """
    return _get_sub_distance(target, et, obs, 'solar')


def get_anti_solar_radec(et: float,
                         obs: str) -> tuple[float, float]:
    """
    Get RA/Dec of the anti-solar point from a given observatory.

    Parameters
    ----------
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    tuple[float, float]
        The RA/Dec of the anti-solar point in [rad].
    """
    ra, dec = get_sky_coordinates('Sun', et, obs)
    ra = (ra + spice.pi()) % spice.twopi()
    return ra, -dec


def get_state(target: str,
              et: float,
              obs: str) -> np.ndarray:
    """
    Get the state (position and velocity) of the target relative to the
    observer. Also separate the components.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    np.ndarray
        The state of the body in the J2000 reference frame. The first three
        indices of the array are the (x, y, z) position vector components in
        [km], and the second three indices of the array are the (dx/dt, dy/dt,
        dz/dt) velocity vector components in [km/s].
    """
    starg, _ = spice.spkezr(target, et, ref, abcorr, obs)
    return starg


def get_limb_radec(target: str,
                   et: float,
                   obs: str,
                   resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Get the J2000 right ascension and declination of the target body's apparent
    limb for a given observer.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    resolution : int
        The number of limb points.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The right ascension and declination of the target body's apparent limb
        in units of [rad].
    """
    fixref = get_target_frame(target)
    params = dict(method='TANGENT/ELLIPSOID',
                  target=target,
                  et=et,
                  fixref=fixref,
                  abcorr=abcorr,
                  corloc=corloc,
                  obsrvr=obs,
                  refvec=refvec,
                  rolstp=spice.twopi()/resolution,
                  ncuts=resolution+1,
                  schstp=1e-4,
                  soltol=1e-7,
                  maxn=resolution+1)
    _, _, epochs, tangts = spice.limbpt(**params)
    ra = []
    dec = []
    for epoch, point in zip(epochs, tangts):
        rotation_matrix = spice.pxfrm2(fixref, 'J2000', epoch, et)
        point = spice.mxv(rotation_matrix, point)
        point = spice.recrad(point)
        ra.append(point[1])
        dec.append(point[2])
    return np.array(ra), np.array(dec)


def get_terminator_radec(target: str,
                         et: float,
                         obs: str,
                         resolution: int,
                         kind: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Get the J2000 right ascension and declination of the target body's
    terminator for a given observer.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    resolution : int
        The number of limb points.
    kind : str
        Type of terminator. Options are 'umbral' or 'penumbral'.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The right ascension and declination of the target body's terminator in
        units of [rad].
    """
    fixref = get_target_frame(target)
    termpt_params = dict(method=f'{kind.upper()}/TANGENT/ELLIPSOID',
                         ilusrc='Sun',
                         target=target,
                         et=et,
                         fixref=fixref,
                         abcorr=abcorr,
                         corloc='ELLIPSOID TERMINATOR',
                         obsrvr=obs,
                         refvec=refvec,
                         rolstp=spice.twopi()/resolution,
                         ncuts=resolution+1,
                         schstp=1e-4,
                         soltol=1e-7,
                         maxn=resolution+1)
    illumf_params = dict(method=subpoint_method,
                         target=target,
                         ilusrc='Sun',
                         et=et,
                         fixref=fixref,
                         abcorr=abcorr,
                         obsrvr=obs)
    _, spoints, epochs, _ = spice.termpt(**termpt_params)
    ra = []
    dec = []
    for epoch, spoint in zip(epochs, spoints):
        trgepc, vec, _, _, _, visible, _ = spice.illumf(spoint=spoint,
                                                           **illumf_params)
        if visible:
            rotation_matrix = spice.pxfrm2(fixref, 'J2000', epoch, et)
            point = spice.mxv(rotation_matrix, vec)
            _, ra_i, dec_i = spice.recrad(point)
            ra.append(ra_i)
            dec.append(dec_i)
        else:
            ra.append(np.nan)
            dec.append(np.nan)
    ra = np.array(ra)
    dec = np.array(dec)

    nans = np.where(np.isnan(ra))[0]
    if nans[0] != 0:
        ra = np.roll(ra, -nans[0])
        dec = np.roll(dec, -nans[0])
    ra = ra[~np.isnan(ra)]
    dec = dec[~np.isnan(dec)]
    return ra, dec


def get_latlon_radec(target: str,
                     et: float,
                     obs: str,
                     latitude: float,
                     longitude: float) -> tuple[float, float]:
    """
    Get the J2000 right ascension and declination of the target body's apparent
    limb for a given observer.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    latitude : float
        Planetodetic longitude on the target body in units of [rad].
    longitude : float
        Planetodetic latitude on the target body in units of [rad].

    Returns
    -------
    tuple[float, float]
        The right ascension and declination of the latitude/longitude pair.
        Returns NaN if the point is not visible to the observer at the
        specified time.
    """
    params = dict(method=subpoint_method,
                  target=target,
                  ilusrc='Sun',
                  et=et,
                  fixref=get_target_frame(target),
                  abcorr=abcorr,
                  obsrvr=obs)
    re = get_equatorial_radius(target, longitude)
    f = get_flattening_coefficient(target, longitude)
    spoint = spice.pgrrec(target, longitude, latitude, 0, re, f)
    trgepc, vec, _, _, _, visible, _ = spice.illumf(spoint=spoint, **params)
    if visible:
        rotation_matrix = spice.pxfrm2(get_target_frame(target), 'J2000',
                                       trgepc, et)
        point = spice.mxv(rotation_matrix, vec)
        _, ra, dec = spice.recrad(point)
    else:
        ra, dec = np.nan, np.nan
    return ra, dec


def determine_occultation(target1: str,
                          target2: str,
                          et: float,
                          obs: str,
                          aberration: str = abcorr) -> int:
    """
    Determine if one target occults another. The output code has the following
    meaning:

    +------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Code | Meaning                                                                                                                                                                            |
    +======+====================================================================================================================================================================================+
    | -3   | Total occultation of first target by second. This occurs when the second target totally blocks the first target.                                                                   |
    +------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | -2   | Annular occultation of first target by second. The second target does not block the limb of the first. For instance, this occurs when the first target is a small moon in transit. |
    +------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | -1   | Partial occultation of first target by second target. This is the transition phase when the second target begins to transit across the first target.                               |
    +------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    |  0   | No occultation or transit: both objects are completely visible to the observer.                                                                                                    |
    +------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    |  1   | Partial occultation of second target by first target. This is the transition phase when the first target begins to transit across the second target.                               |
    +------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    |  2   | Annular occultation of second target by first. For instance, this occurs when the first target is a small moon transiting the second target.                                       |
    +------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    |  3   | Total occultation of second target by first. This occurs when the first target totall blocks the second target.                                                                    |
    +------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

    Parameters
    ----------
    target1 : str
        The name of the first target body, assumed to be the primary target of
        interest.
    target2 : str
        The name of the second target body, probably the one casting the
        shadow.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    aberration: str
        Aberration correction. Default is 'LT'.

    Returns
    -------
    int
        The occultation identification code.
    """
    params = dict(target1=target1,
                  shape1='ELLIPSOID',
                  frame1=get_target_frame(target1),
                  target2=target2,
                  shape2='ELLIPSOID',
                  frame2=get_target_frame(target2),
                  abcorr=aberration,
                  observer=obs,
                  et=et)
    ocltid = spice.occult(**params)
    return ocltid


def determine_if_in_shadow(target1: str,
                           target2: str,
                           et: float,
                           obs: str) -> int:
    """
    Get the code indicating the occultation status of `target1` by `target2` as
    seen by the Sun. Interpret this as determining when `target1` is in eclipse
    behind `target2` (code > 0) or `target1` transiting `target1` (code < 0) as
    viewed from the Sun. If code = 0, then neither object occults the other.

    These codes seem to be opposite to those listed on SpiceyPy, but I don't
    know why. For instance, when `target1` passes behind `target2` as viewed by
    the Sun, the code should be < 0, but SpiceyPy returns > 0.

    Parameters
    ----------
    target1 : str
        The primary target of interest. Probably a satelite.
    target2 : str
        The body casting a shadow. Probably a primary Solar System planet.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.

    Returns
    -------
    int
        The occultation code.
    """
    # get apparent ephemeris times at target and Sun as viewed from target
    et_target = get_apparent_epoch(target1, et, obs)
    et_sun = get_apparent_epoch(target1, et_target, 'Sun')

    # get occultation code with Sun as the observer
    return determine_occultation(target1, target2, et_sun, 'Sun', 'XLT')


def get_dayside_radec(target: str,
                      et: float,
                      obs: str,
                      resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Get the J2000 right ascension and declination of the target body's apparent
    dayside disk for a given observer.

    Parameters
    ----------
    target : str
        Name of target body.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    resolution : int
        The number of disk points.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The right ascension and declination of the target body's apparent
        dayside disk in units of [rad].
    """
    fixref = get_target_frame(target)
    params = dict(method='TANGENT/ELLIPSOID',
                  target=target,
                  et=et,
                  fixref=fixref,
                  abcorr=abcorr,
                  corloc='ELLIPSOID LIMB',
                  obsrvr=obs,
                  refvec=refvec,
                  rolstp=spice.twopi()/resolution,
                  ncuts=resolution+1,
                  schstp=1e-4,
                  soltol=1e-7,
                  maxn=resolution+1)
    _, limb_points, _, _ = spice.limbpt(**params)

    params = dict(method='UMBRAL/TANGENT/ELLIPSOID',
                  target=target,
                  et=et,
                  fixref=fixref,
                  abcorr=abcorr,
                  corloc='ELLIPSOID TERMINATOR',
                  obsrvr=obs,
                  refvec=refvec,
                  ilusrc='Sun',
                  rolstp=spice.twopi()/resolution,
                  ncuts=resolution+1,
                  schstp=1e-4,
                  soltol=1e-7,
                  maxn=resolution+1)
    _, term_points, _, _ = spice.termpt(**params)

    ra = []
    dec = []
    _, radii = spice.bodvrd(target, 'RADII', 3)
    for lpoint, tpoint in zip(limb_points, term_points):
        trgepc, srfvec, _, _, _, visibl, _ = spice.illumf(
            'ELLIPSOID', target, 'Sun', et, fixref, abcorr, obs, tpoint)
        if visibl:
            rotation_matrix = spice.pxfrm2(fixref, 'J2000', trgepc, et)
            point = spice.mxv(rotation_matrix, srfvec)
            point = spice.recrad(point)
            ra.append(point[1])
            dec.append(point[2])
        else:
            trgepc, srfvec, _, _, _, visibl, _ = spice.illumf(
                'ELLIPSOID', target, 'Sun', et, fixref, abcorr, obs, lpoint)
            rotation_matrix = spice.pxfrm2(get_target_frame(target),
                                              'J2000', trgepc, et)
            point = spice.mxv(rotation_matrix, srfvec)
            point = spice.recrad(point)
            ra.append(point[1])
            dec.append(point[2])

    return np.array(ra), np.array(dec)


def get_shadow_intercept(target1: str,
                         target2: str,
                         et: float,
                         obs: str,
                         resolution: int,
                         scale: int | float = 1.0
                         ) -> tuple[np.ndarray, np.ndarray]:
    """
    Get the J2000 right ascension and declination of the intersection between
    `target2`'s shadow with `target1`'s disk.

    Parameters
    ----------
    target1 : str
        The primary target of interest. Probably a satelite.
    target2 : str
        The body casting a shadow. Probably a primary Solar System planet.
    et : float
        Observer epoch.
    obs : str
        Name of observing body, e.g., 'Earth' or 'Keck' or 'JWST'.
    resolution : int
        The number of disk points.
    scale : int | float
        Scaling factor to increase resolution for intersection.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, bool]
        The right ascension and declination of the `target2`'s shadow on the
        disk of `target1` in units of [rad].
    """

    # get fixed reference frames and target1 radii
    fixref1 = get_target_frame(target1)
    fixref2 = get_target_frame(target2)
    _, radii = spice.bodvrd(target1, 'RADII', 3)
    
    # get apparent ephemeris times at target and Sun as viewed from target
    et_target = get_apparent_epoch(target1, et, obs, '<-')
    et_sun = get_apparent_epoch('Sun', et_target, target1, '<-')

    # get vectors from Sun to parent's limb at much higher resolution to 
    # account for projected size differences
    hires = int(resolution * scale)
    params = dict(method='TANGENT/ELLIPSOID',
                  target=target2,
                  et=et_sun,
                  fixref=fixref2,
                  abcorr='XLT',
                  corloc=corloc,
                  obsrvr='Sun',
                  refvec=refvec,
                  rolstp=spice.twopi()/hires,
                  ncuts=hires+1,
                  schstp=1e-4,
                  soltol=1e-7,
                  maxn=hires+1)
    _, _, _, tangts = spice.limbpt(**params)

    # first get indices
    indices = []
    for i, tangt in enumerate(tangts):
        try:
            spice.sincpt(subpoint_method, target1, et_sun, fixref1, 'XLT',
                            'Sun', fixref2, tangt)
            indices.append(i)
        except NotFoundError:
            continue
    indices = np.concatenate(([indices[0] - 1], indices, [indices[-1] + 1]))
    indices %= tangts.shape[0]

    # calculate surface intercept on target of parent's limb vectors and
    # determine visibility from observer
    ra, dec = [], []
    for tangt in tangts[indices]:
        try:
            spoint, epoch, srfvec = spice.sincpt(
                subpoint_method, target1, et_sun, fixref1, 'XLT', 'Sun',
                fixref2, tangt)
            trgepc, srfvec, _, _, _, visibl, _ = spice.illumf(
                'ELLIPSOID', target1, 'Sun', et, fixref1, abcorr, obs, spoint)
            if visibl:
                rotation_matrix = spice.pxfrm2(fixref1, 'J2000', trgepc, et)
                tpoint = spice.mxv(rotation_matrix, srfvec)
                tpoint = spice.recrad(tpoint)
                ra.append(tpoint[1])
                dec.append(tpoint[2])
        except NotFoundError:
            _, alt, _, srfpt, trgepc, srfvec = spice.tangpt(
                'ELLIPSOID', target1, et_sun, fixref1, 'XLT',
                'SURFACE POINT', 'Sun', fixref2, tangt)
            trgepc, srfvec, _, _, _, visibl, _ = spice.illumf(
                'ELLIPSOID', target1, 'Sun', et, fixref1, abcorr, obs, srfpt)
            rotation_matrix = spice.pxfrm2(fixref1, 'J2000', trgepc, et)
            tpoint = spice.mxv(rotation_matrix, srfvec)
            tpoint = spice.recrad(tpoint)
            ra.append(tpoint[1])
            dec.append(tpoint[2])
            
    return np.flip(ra), np.flip(dec)
