import warnings

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, Angle
from astropy.visualization.wcsaxes.core import WCSAxes
from astropy.wcs import WCS
from matplotlib.patches import Polygon

from planetviewer.spice_functions import ref

wcs_shape = (3, 3)
"""Default WCS window shape."""

standard_colors = {
    'red': '#D62728',
    'orange': '#FF7F0E',
    'yellow': '#FDB813',
    'green': '#2CA02C',
    'blue': '#0079C1',
    'violet': '#9467BD',
    'cyan': '#17BECF',
    'magenta': '#D64ECF',
    'brown': '#8C564B',
    'lightestgrey': '#ECECEC',
    'lightergrey': '#E5E5E5',
    'lightgrey': '#C9C9C9',
    'grey': '#AEAEAE',
    'darkgrey': '#838383',
    'darkergrey': '#3A3A3A',
    'black': '#000000',
    'white': '#FFFFFF',
}
"""A set of standard colors."""


def _default_limb_kwargs() -> dict:
    return dict(closed=True,
                edgecolor=standard_colors['black'],
                facecolor='none',
                linewidth=plt.rcParams['lines.linewidth']/2)


def _default_disk_kwargs() -> dict:
    return dict(closed=True,
                edgecolor='none',
                facecolor=standard_colors['white'],
                linewidth=0.0)


def _default_night_patch_kwargs() -> dict:
    return dict(facecolor=standard_colors['grey'],
                edgecolor='none',
                linewidth=0,
                closed=True)


def _default_latlon_kwargs() -> dict:
    return dict(edgecolor=standard_colors['darkergrey'],
                facecolor='none',
                closed=False,
                linewidth=plt.rcParams['lines.linewidth']/4,
                capstyle='round')


def _default_lit_ring_edge_kwargs() -> dict:
    return dict(edgecolor=standard_colors['darkergrey'],
                facecolor='none',
                closed=False,
                linewidth=plt.rcParams['lines.linewidth']/4)


def _default_lit_ring_fill_kwargs() -> dict:
    return dict(edgecolor='none',
                facecolor=standard_colors['lightestgrey'],
                closed=True)


def _default_eclipsed_ring_edge_kwargs() -> dict:
    return dict(edgecolor=standard_colors['darkergrey'],
                facecolor='none',
                closed=False,
                linewidth=plt.rcParams['lines.linewidth']/4)


def _default_eclipsed_ring_fill_kwargs() -> dict:
    return dict(closed=True,
                edgecolor='none',
                facecolor=standard_colors['grey'])


def _default_pericenter_marker_kwargs() -> dict:
    return dict(color=standard_colors['black'],
                marker='.',
                s=(4*plt.rcParams['lines.linewidth'])**2,
                linewidth=plt.rcParams['lines.linewidth']/4,
                capstyle='round')


def _default_arc_kwargs() -> dict:
    return dict(edgecolor=standard_colors['black'],
                facecolor='none',
                closed=False,
                linewidth=plt.rcParams['lines.linewidth'],
                capstyle='round')


# noinspection PyUnresolvedReferences
def make_wcs(center: SkyCoord,
             fov: Angle) -> WCS:
    """
    Convenience function for generating a WCS projection for plotting and
    coordinate transforms.

    Parameters
    ----------
    center : SkyCoord
        The center of the field of view.
    fov : Angle
        The width/height of the field of view.

    Returns
    -------
    WCS
        The WCS for plotting and coordinate transforms.
    """
    fov *= 1.000001  # helps ensure the actual edge values are shown
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [wcs_shape[1] / 2 + 0.5, wcs_shape[0] / 2 + 0.5]
    wcs.wcs.crval = [center.ra.degree, center.dec.degree]
    wcs.wcs.cunit = ["deg", "deg"]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cdelt = [-fov.degree/wcs_shape[1], fov.degree/wcs_shape[0]]
    return wcs


def set_standard_axis_labels(axis: WCSAxes) -> None:
    """
    Set axis labels to right ascension and declination for the standard SPICE
    output reference vector (default is J2000).

    Parameters
    ----------
    axis : WCSAxes
        The axis in need of labelling.

    Returns
    -------
    None
        None.
    """
    axis.set_xlabel(f'Right Ascension ({ref})')
    axis.set_ylabel(f'Declination ({ref})')


def set_standard_axis_limits(axis: WCSAxes) -> None:
    """
    Reset the plot limits of a WCS axis. This is usually necessary after
    drawing objects in the field of view because Matplotlib will scale the axis
    to fit everything.

    Parameters
    ----------
    axis : WCSAxes
        The axis in need of limiting.

    Returns
    -------
    None
        None.
    """
    axis.set_xlim(-0.5, wcs_shape[0] - 0.5)
    axis.set_ylim(-0.5, wcs_shape[1] - 0.5)
    axis.set_aspect('equal')


def convert_to_relative_axis(center: SkyCoord,
                             axis: WCSAxes) -> None:
    """
    Convert an axis from absolute sky coordinates to relative offset from the
    center of the field of view.

    Parameters
    ----------
    center : SkyCoord
        The center of the field of view.
    axis : WCSAxes
        The WCS axis in need of converting.

    Returns
    -------
    None
        None.
    """
    # Remove the absolute coordinates
    ra = axis.coords['ra']
    dec = axis.coords['dec']
    ra.set_ticks_visible(False)
    ra.set_ticklabel_visible(False)
    dec.set_ticks_visible(False)
    dec.set_ticklabel_visible(False)
    ra.set_axislabel('')
    dec.set_axislabel('')

    # Create an overlay with relative coordinates
    aframe = center.skyoffset_frame()
    overlay = axis.get_coords_overlay(aframe)
    ra_offset = overlay['lon']
    dec_offset = overlay['lat']
    ra_offset.set_axislabel(f'Relative Right Ascension')
    dec_offset.set_axislabel(f'Relative Declination')
    ra_offset.set_ticks_position('bt')
    ra_offset.set_ticklabel_position('b')
    dec_offset.set_ticks_position('lr')
    dec_offset.set_ticklabel_position('l')
    ra_offset.set_axislabel_position('b')
    dec_offset.set_axislabel_position('l')
    ra_offset.coord_wrap = 180 * u.degree


def _parse_kwargs(kwargs: dict,
                  default_kwargs: dict) -> dict:
    """
    Combine any user-provided kwargs with existing default kwargs, overriding 
    the defaults when necessary.
    
    Parameters
    ----------
    kwargs : dict
        User provided kwargs.
    default_kwargs : dict
        Default kwargs.

    Returns
    -------
    dict
        The combined kwargs.
    """
    parsed_kwargs = dict()
    if kwargs is None:  # noqa
        kwargs = dict()
    for key in default_kwargs.keys():
        parsed_kwargs[key] = default_kwargs[key]
    for key in kwargs.keys():
        parsed_kwargs[key] = kwargs[key]
    return parsed_kwargs


def _offset_coord(coord: SkyCoord,
                  dra: Angle,
                  ddec: Angle) -> SkyCoord:
    """
    Calcualte the spherical offset of a coordinate given a change in RA and
    Dec.

    Parameters
    ----------
    coord : SkyCoord
        The original coordinate.
    dra : Angle
        The offset in right ascension.
    ddec : Angle
        The offset in declination.

    Returns
    -------
    SkyCoord
        The offset coordinate as an Astropy `SkyCoord` object.
    """
    return coord.spherical_offsets_by(dra, ddec)


def _retrieve_coords(coords: list[SkyCoord],
                     dra: Angle,
                     ddec: Angle) -> tuple[np.ndarray, np.ndarray]:
    """
    Get RA and Dec coordinates from a list of coordinates.

    Parameters
    ----------
    coords : list[SkyCoord]
        The list of coordinates.
    dra : Angle
        Optional angular offset in right ascension.
    ddec : Angle
        Optional angular offset in declination.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        A tuple of arrays containing the RA and Dec in degrees.
    """
    ra = []
    dec = []
    for coord in coords:
        new_coord = _offset_coord(coord, dra, ddec)
        ra.append(new_coord.ra.degree)
        dec.append(new_coord.dec.degree)
    return np.array(ra), np.array(dec)


# noinspection DuplicatedCode
def plot_limb(axis: WCSAxes,
              coords: list[SkyCoord],
              dra: Angle = Angle(0, unit='deg'),
              ddec: Angle = Angle(0, unit='deg'),
              **kwargs) -> None:
    """
    Draw the target body's limb.
    
    Parameters
    ----------
    axis : WCSAxes
        The axis in which to plot the limb.
    coords : list[SkyCoord]
        The coordinates of the target body's limb. Can be generated with the
        `Planet` method 'get_limb_sky_coordinates'.
    dra : Angle
        Optional angular offset in right ascension.
    ddec : Angle
        Optional angular offset in declination.
    **kwargs
        Additional Matplotlib `plot` keyword arguments like `linewidth`,
        `color`, etc.

    Returns
    -------
    None
        None.
    """
    transform = axis.get_transform('world')
    kwargs = _parse_kwargs(kwargs, _default_limb_kwargs())

    ra, dec = _retrieve_coords(coords, dra, ddec)
    xy = np.concatenate(([ra], [dec])).T

    limb = Polygon(xy, transform=transform, **kwargs)
    axis.add_patch(limb)


# noinspection DuplicatedCode
def plot_disk(axis: WCSAxes,
              coords: list[SkyCoord],
              dra: Angle = Angle(0, unit='deg'),
              ddec: Angle = Angle(0, unit='deg'),
              **kwargs) -> None:
    """
    Draw the target body's disk.

    Parameters
    ----------
    axis : WCSAxes
        The axis in which to plot the disk.
    coords : list[SkyCoord]
        The coordinates of the target body's limb. Can be generated with the
        `Planet` method 'get_limb_sky_coordinates'.
    dra : Angle, optional
        Angular offset in right ascension.
    ddec : Angle, optional
        Angular offset in declination.
    **kwargs
        Additional Matplotlib `Polygon` patch keyword arguments like
        `facecolor `, `linewidth`, etc.

    Returns
    -------
    None
        None.
    """
    transform = axis.get_transform('world')
    kwargs = _parse_kwargs(kwargs, _default_disk_kwargs())

    ra, dec = _retrieve_coords(coords, dra, ddec)
    xy = np.concatenate(([ra], [dec])).T

    disk = Polygon(xy, transform=transform, **kwargs)
    axis.add_patch(disk)


def plot_nightside(axis: WCSAxes,
                   limb_coords: list[SkyCoord],
                   terminator_coords: list[SkyCoord],
                   dra: Angle = Angle(0, unit='deg'),
                   ddec: Angle = Angle(0, unit='deg'),
                   **kwargs):
    """
    Add a shaded patch for the disk's nightside.

    Parameters
    ----------
    axis : WCSAxes
        The axis in which to plot the nightside shading.
    limb_coords : list[SkyCoord]
        The coordinates of the target body's limb. Can be generated with the
        `Planet` method 'get_limb_sky_coordinates'.
    terminator_coords : list[SkyCoord]
        The coordinates of the target body's terminator. Can be generated with
        the `Planet` method 'get_terminator_sky_coordinates'.
    dra : Angle, optional
        Angular offset in right ascension.
    ddec : Angle, optional
        Angular offset in declination.
    **kwargs
        Additional Matplotlib `Polygon` patch keyword arguments like
        `facecolor `, `linewidth`, etc.

    Returns
    -------
    None
        None.
    """
    transform = axis.get_transform('world')
    kwargs = _parse_kwargs(kwargs, _default_night_patch_kwargs())

    # first calculate the limb points
    limb_ra, limb_dec = _retrieve_coords(limb_coords, dra, ddec)

    # then calculate the terminator points
    terminator_ra, terminator_dec = _retrieve_coords(
        terminator_coords, dra, ddec)

    # find the intersection of the limb with the apparent terminator
    try:
        start = np.sqrt((limb_ra - terminator_ra[0])**2 +
                        (limb_dec - terminator_dec[0])**2).argmin()
        end = np.sqrt((limb_ra - terminator_ra[-1])**2 +
                      (limb_dec - terminator_dec[-1])**2).argmin()
        ind = np.arange(start, start+np.size(limb_ra), 1) % limb_ra.size
        ind0 = np.where(ind == start)[0][0]
        ind1 = np.where(ind == end)[0][0]
        ind = ind[ind0:ind1]
    # if the observer is the Sun, then there is no nightside!
    except IndexError:
        return

    # select the limb section
    limb_ra = np.flip(limb_ra[ind])
    limb_dec = np.flip(limb_dec[ind])

    # join the coordinates
    ra = np.concatenate((terminator_ra, limb_ra))
    dec = np.concatenate((terminator_dec, limb_dec))
    xy = np.concatenate(([ra], [dec])).T

    # plot the night side
    nightside = Polygon(xy, transform=transform, **kwargs)
    axis.add_patch(nightside)


# TODO: try to convert this to an intersection of two Polygons, that way it
#  won't matter what the shapes are or whether or not the shadow is smaller or
#  larger than the primary target.
def plot_primary_shadow(axis: WCSAxes,
                        disk_coords: list[SkyCoord],
                        shadow_coords: list[SkyCoord],
                        dra: Angle = Angle(0, unit='deg'),
                        ddec: Angle = Angle(0, unit='deg'),
                        **kwargs):
    """
    Add a shaded patch for the shadow cast by a primary body (like Jupiter's
    shadow cast on Ganymede).

    Parameters
    ----------
    axis : WCSAxes
        The axis in which to plot the nightside shading.
    disk_coords : list[SkyCoord]
        The coordinates of the target body's apparent dayside. Can be generated
        with the `Planet` method 'get_dayside_sky_coordinates'.
    shadow_coords : list[SkyCoord]
        The coordinates of the shadow casting body's shadow where it intersects
        the target body's disk. Can be generated with the
        `Planet` method 'get_shadow_intersection_sky_coordinates'.
    dra : Angle, optional
        Angular offset in right ascension.
    ddec : Angle, optional
        Angular offset in declination.
    **kwargs
        Additional Matplotlib `Polygon` patch keyword arguments like
        `facecolor `, `linewidth`, etc.

    Returns
    -------
    None
        None.
    """
    plot_nightside(axis, disk_coords, shadow_coords, dra, ddec,
                   **kwargs)


# TODO: add examples in each of these docstrings.
# noinspection DuplicatedCode
def plot_latlon(axis: WCSAxes,
                coords: list[SkyCoord],
                dra: Angle = Angle(0, unit='deg'),
                ddec: Angle = Angle(0, unit='deg'),
                **kwargs) -> None:
    """
    Plot latitude or longitude gridlines.

    Parameters
    ----------
    axis : WCSAxes
        The axis in which to plot the latitude or longitude gridlines.
    coords : list[SkyCoord]
        The coordinates of the target body's latitude or longitude line. Can be
        generated with the `Planet` method `get_longitude_line_coordinates` or
        `get_latitude_line_coordinates`.
    dra : Angle, optional
        Angular offset in right ascension.
    ddec : Angle, optional
        Angular offset in declination.
    **kwargs
        Additional Matplotlib `plot` keyword arguments like `linewidth`,
        `color`, etc.

    Returns
    -------
    None
        None.
    """
    transform = axis.get_transform('world')
    kwargs = _parse_kwargs(kwargs, _default_latlon_kwargs())

    ra, dec = _retrieve_coords(coords, dra, ddec)
    xy = np.concatenate(([ra], [dec])).T

    line = Polygon(xy, transform=transform, **kwargs)
    axis.add_patch(line)


def _expand_ind(ind: np.ndarray, size: int) -> np.ndarray:
    """
    Add an additional index at the end for slicing ring boundary arrays.
    """
    if len(ind) == 0:
        return ind
    else:
        end = ind[-1] + 1
        return np.concatenate((ind, [end])) % size


def _remove_jumps(ind: np.ndarray) -> np.ndarray:
    loc = np.where(np.diff(ind) > 1)[0]
    if len(loc) == 0:
        return ind
    else:
        loc = loc[0] + 1
        return np.roll(ind, -loc)


def _find_sets(ind: np.ndarray) -> list[np.ndarray]:
    locs = np.where(np.diff(ind) > 1)[0]
    if len(locs) == 0:
        return [ind]
    else:
        locs += 1
        return np.array_split(ind, locs)


def _get_side_indices(distances: np.ndarray,
                      side: str) -> np.ndarray:
    center = distances.argmin()
    count = distances.size
    if side == 'front':
        ind = np.arange(center - count / 4, center + count / 4) % count
        ind[np.where(ind < 0)] += count
    else:
        ind = np.arange(center + count / 4, center + 3 * count / 4) % count
    return ind.astype(int)


# noinspection DuplicatedCode
def plot_ring(axis: WCSAxes,
              coords: list[list[SkyCoord]],
              distances: list[np.ndarray],
              eclipsed: list[np.ndarray],
              side: str,
              dra: Angle = Angle(0, unit='deg'),
              ddec: Angle = Angle(0, unit='deg'),
              eclipse_edge_kwargs: dict = None,
              eclipse_fill_kwargs: dict = None,
              lit_edge_kwargs: dict = None,
              lit_fill_kwargs: dict = None) -> None:
    """
    Plot a planet's rings (if any). Rings with non-zero width are filled with
    light grey.

    Parameters
    ----------
    axis : WCSAxes
        The axis in which to plot the ring.
    coords : list[list[SkyCoord]]
        The coordinates of the ring. Can be generated with the `Ring` method
        `get_sky_coordinates` WITH the option `boundaries=True`. If the list
        contains a single list of coordinates, then the ring is assumed to have
        no width. If the list contains two lists of coordinates, they are
        assumed to correspond to the inner and outer radius.
    distances : list[np.ndarray]
        Distances from the observer to each ring point. Can be generated with
        the `Ring` method `get_distances`. Like the ring coordinates, if the
        list contains a single list of distances, then the ring is assumed to
        have no width. If the list contains two lists of distances, they are
        assumed to correspond to the inner and outer radius.
    eclipsed: list[np.ndarray]
        A boolean array indicating whether or not a ring point is eclipsed by
        the planet. Can be generated with the `Ring` method `get_eclipsed`.
        Like the ring coordinates, if the list contains a single list of
        boolean values, then the ring is assumed to have no width. If the list
        contains two arrays of boolean values, they are assumed to correspond
        to the inner and outer radius.
    side : str
        The side of the rings to plot, either 'front' or 'back'. You'll want to
        plot the back first, then plot the planet, then the front.
    dra : Angle, optional
        Angular offset in right ascension.
    ddec : Angle, optional
        Angular offset in declination.
    eclipse_edge_kwargs : dict, optional
        Additional Matplotlib `plot` keyword arguments like `linewidth`,
        `color`, etc. for the edges of the eclipsed portions of the rings.
    eclipse_fill_kwargs : dict, optional
        Additional Matplotlib `Polygon` patch keyword arguments like
        `facecolor `, `linewidth`, etc. for the fills of the eclipsed portions
        of the rings.
    lit_edge_kwargs : dict, optional
        Additional Matplotlib `plot` keyword arguments like `linewidth`,
        `color`, etc. for the edges of the eclipsed portions of the rings.
    lit_fill_kwargs : dict, optional
        Additional Matplotlib `Polygon` patch keyword arguments like
        `facecolor `, `linewidth`, etc. for the fills of the eclipsed portions
        of the rings.

    Returns
    -------
    None
        None.
    """
    transform = axis.get_transform('world')
    eclipse_edge_kwargs = _parse_kwargs(eclipse_edge_kwargs,
                                        _default_eclipsed_ring_edge_kwargs())
    eclipse_fill_kwargs = _parse_kwargs(eclipse_fill_kwargs,
                                        _default_eclipsed_ring_fill_kwargs())
    lit_edge_kwargs = _parse_kwargs(lit_edge_kwargs,
                                    _default_lit_ring_edge_kwargs())
    lit_fill_kwargs = _parse_kwargs(lit_fill_kwargs,
                                    _default_lit_ring_fill_kwargs())

    # get indices
    ind = [_get_side_indices(i, side) for i in distances]
    count = distances[0].size
    eclipsed = [np.where(i)[0] for i in eclipsed]

    # decide if eclisped portion should be drawn
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        eclipsed_distance = np.mean(distances[0][eclipsed[0]])
        avg_ring_distance = np.mean(distances[0])
    if (side == 'front') & (eclipsed_distance < avg_ring_distance):
        draw_eclipsed = True
    elif (side == 'back') & (eclipsed_distance > avg_ring_distance):
        draw_eclipsed = True
    else:
        draw_eclipsed = False

    # draw 2D ring
    if len(ind) == 2:
        # get coordinates
        ra_inner, dec_inner = _retrieve_coords(coords[0], dra, ddec)
        ra_outer, dec_outer = _retrieve_coords(coords[1], dra, ddec)

        # get eclipsed coordinates
        eclipsed_inner = _expand_ind(_remove_jumps(eclipsed[0]), count)
        eclipsed_outer = _expand_ind(_remove_jumps(eclipsed[1]), count)
        ra_eclipsed_inner = ra_inner[eclipsed_inner]
        dec_eclipsed_inner = dec_inner[eclipsed_inner]
        ra_eclipsed_outer = ra_outer[eclipsed_outer]
        dec_eclipsed_outer = dec_outer[eclipsed_outer]
        ra_eclipsed_combined = np.concatenate([ra_eclipsed_inner,
                                               np.flip(ra_eclipsed_outer)])
        dec_eclipsed_combined = np.concatenate([dec_eclipsed_inner,
                                                np.flip(dec_eclipsed_outer)])

        # draw eclipsed shapes if they exist
        if (ra_eclipsed_combined.size != 0) & draw_eclipsed:
            xy = np.concatenate(([ra_eclipsed_combined],
                                 [dec_eclipsed_combined])).T
            polygon = Polygon(xy=xy, transform=transform,
                              **eclipse_fill_kwargs)
            axis.add_patch(polygon)
        if (ra_eclipsed_inner.size != 0) & draw_eclipsed:
            xy = np.concatenate(([ra_eclipsed_inner], [dec_eclipsed_inner])).T
            line = Polygon(xy, transform=transform, **eclipse_edge_kwargs)
            axis.add_patch(line)
        if (ra_eclipsed_outer.size != 0) & draw_eclipsed:
            xy = np.concatenate(([ra_eclipsed_outer], [dec_eclipsed_outer])).T
            line = Polygon(xy, transform=transform, **eclipse_edge_kwargs)
            axis.add_patch(line)
        
        # get lit coordinates
        inner_sets = _find_sets(
            np.intersect1d(np.setxor1d(eclipsed[0], ind[0]), ind[0]))
        outer_sets = _find_sets(
            np.intersect1d(np.setxor1d(eclipsed[1], ind[1]), ind[1]))
        for i, j in zip(inner_sets, outer_sets):
            inner = _expand_ind(_remove_jumps(i), count)
            outer = _expand_ind(_remove_jumps(j), count)
            ra_lit_inner = ra_inner[inner]
            dec_lit_inner = dec_inner[inner]
            ra_lit_outer = ra_outer[outer]
            dec_lit_outer = dec_outer[outer]
            ra_lit_combined = np.concatenate([ra_lit_inner,
                                              np.flip(ra_lit_outer)])
            dec_lit_combined = np.concatenate([dec_lit_inner,
                                               np.flip(dec_lit_outer)])

            # draw lit shapes if they exist
            if ra_lit_combined.size != 0:
                xy = np.concatenate(([ra_lit_combined],
                                     [dec_lit_combined])).T
                polygon = Polygon(xy=xy, transform=transform,
                                  **lit_fill_kwargs)
                axis.add_patch(polygon)
            if ra_lit_inner.size != 0:
                xy = np.concatenate(([ra_lit_inner], [dec_lit_inner])).T
                line = Polygon(xy, transform=transform, **lit_edge_kwargs)
                axis.add_patch(line)
            if ra_lit_outer.size != 0:
                xy = np.concatenate(([ra_lit_outer], [dec_lit_outer])).T
                line = Polygon(xy, transform=transform, **lit_edge_kwargs)
                axis.add_patch(line)

    # draw 1D ring
    else:
        # get coordinates
        ra, dec = _retrieve_coords(coords[0], dra, ddec)
        ra_eclipsed, dec_eclipsed = ra[eclipsed[0]], dec[eclipsed[0]]

        # draw eclipsed section if it exists
        if (ra_eclipsed.size != 0) & draw_eclipsed:
            xy = np.concatenate(([ra_eclipsed], [dec_eclipsed])).T
            line = Polygon(xy, transform=transform, **eclipse_edge_kwargs)
            axis.add_patch(line)

        # lit portions
        sets = _find_sets(
            np.intersect1d(np.setxor1d(eclipsed[0], ind[0]), ind[0]))
        for i in sets:
            ind = _expand_ind(_remove_jumps(i), count)
            ra_lit, dec_lit = ra[ind], dec[ind]

            # add shapes if they exist
            if ra_lit.size != 0:
                xy = np.concatenate(([ra_lit], [dec_lit])).T
                line = Polygon(xy, transform=transform, **lit_edge_kwargs)
                axis.add_patch(line)


def place_ring_pericenter_markers(axis: WCSAxes,
                                  coord: SkyCoord,
                                  distance: u.Quantity,
                                  parent_body_distance: u.Quantity,
                                  side: str,
                                  dra: Angle = Angle(0, unit='deg'),
                                  ddec: Angle = Angle(0, unit='deg'),
                                  **kwargs) -> None:
    """
    Scatterplot markers at the pericenter or eccentric rings. Currently
    applicable only to some of the rings of Uranus.

    Parameters
    ----------
    axis : WCSAxes
        The axis in which to plot the ring.
    coord : SkyCoord
        The coordinates of the ring's pericenter. Can be found with the `Ring`
        method `get_sky_coordinates` with the option `pericenter=True`.
    distance : float
        The distance from the observer to the pericenter point. Can be
        calculated with the `Ring` method `get_distances` with the option
        `pericenter=True`.
    parent_body_distance : float
        Distance to the parent body.
    side : str
        The side of the rings to plot, either 'front' or 'back'. You'll want to
        plot the back first, then plot the planet, then the front.
    dra : Angle, optional
        Angular offset in right ascension.
    ddec : Angle, optional
        Angular offset in declination.
    **kwargs
        Additional Matplotlib `scatter` keyword arguments like `marker`,
        `size`, `color`, etc. for the markers of the ring pericenter. Currently
        applies only to some of Uranus's rings.

    Returns
    -------
    None
        None.
    """
    transform = axis.get_transform('world')
    kwargs = _parse_kwargs(kwargs, _default_pericenter_marker_kwargs())

    if distance is None:
        pass
    else:
        if (side == 'front') & (distance < parent_body_distance):
            place_marker = True
        elif (side == 'back') & (distance > parent_body_distance):
            place_marker = True
        else:
            place_marker = False
        if place_marker:
            new_coord = _offset_coord(coord, dra, ddec)
            ra = new_coord.ra.degree
            dec = new_coord.dec.degree
            axis.scatter(ra, dec, transform=transform, **kwargs)


def plot_arc(axis: WCSAxes,
             coords: list[SkyCoord],
             dra: Angle = Angle(0, unit='deg'),
             ddec: Angle = Angle(0, unit='deg'),
             **kwargs) -> None:
    """
    Plot a planet's ring arcs (if any). Currently applies only to Neptune's
    Adams ring.

    Parameters
    ----------
    axis : WCSAxes
        The axis in which to plot the ring.
    coords : list[SkyCoord]
        The coordinates of the ring arc. Can be generated with the `Ring`
        method `get_arc_sky_coordinates`.
    dra : Angle, optional
        Angular offset in right ascension.
    ddec : Angle, optional
        Angular offset in declination.
    **kwargs
        Additional Matplotlib `plot` keyword arguments like `linewidth`,
        `color`, etc.

    Returns
    -------
    None
        None.
    """
    transform = axis.get_transform('world')
    kwargs = _parse_kwargs(kwargs, _default_arc_kwargs())
    lw = kwargs['linewidth']
    fc = kwargs['edgecolor']

    ra, dec = _retrieve_coords(coords, dra, ddec)
    if len(ra) < 2:
        axis.scatter(ra, dec, marker='.', s=(2*lw)**2, color=fc, ec='none',
                     transform=transform)
    else:
        xy = np.concatenate(([ra], [dec])).T
        line = Polygon(xy, transform=transform, **kwargs)
        axis.add_patch(line)


def _parse_object_label_position(label_position: str) -> dict:
    """
    Take an object label position like "upper right" and determine various
    offsets and positions for Matplotlib `Annotation`.

    Parameters
    ----------
    label_position : str
        The position. Some combination of 'upper', 'center' or 'lower' for
        vertical position and 'left', 'center', or 'right' for horizontal
        position. For example, 'upper right'. If you want the body center, you
        can just supply a single 'center'.

    Returns
    -------
    dict
        A dictionary containing 'ha', 'va', 'dx', 'dy' and 'xytext'. 'ha' and
        'va' are the text alignments. 'dx' and 'dy' are the scale factors for
        the `xy` argument (for example, to place the label 1 planetary radius
        away from the coordinate `(x0, y0)`, you'd supply `xy=(x0+radius*dx,
        y0+radius*dy)`. Finally, 'xytext' is the text offset in points. Here
        I've chosen 1/4 of the fontsize. Seems to look good.
    """
    fs = plt.rcParams['font.size']

    va = {'upper': 'baseline', 'center': 'center', 'lower': 'top'}
    ha = {'left': 'right', 'center': 'center', 'right': 'left'}
    angle = {'center right': 0*u.deg, 'upper right': 45*u.deg,
             'upper center': 90*u.deg, 'upper left': 135*u.deg,
             'center left': 180*u.deg, 'lower left': 235*u.deg,
             'lower center': 270*u.deg, 'lower right': 315*u.deg,
             'center': 0*u.deg}
    offset = fs / 4
    xytext = {'center right': (offset, 0), 'upper right': (offset, offset),
              'upper center': (0, offset), 'upper left': (-offset, offset),
              'center left': (-offset, 0), 'lower left': (-offset, offset),
              'lower center': (0, -offset), 'lower right': (offset, -offset),
              'center': (0, 0)}
    dx = np.cos(angle[label_position]).value
    dy = np.sin(angle[label_position]).value

    if label_position == 'center':
        vertical, horizontal = label_position, label_position
        dx, dy = 0, 0
    else:
        vertical, horizontal = label_position.split(' ')

    out = {'ha': ha[horizontal],
           'va': va[vertical],
           'dx': dx,
           'dy': dy,
           'xytext': xytext[label_position]}

    return out


def place_label(axis: WCSAxes,
                label: str,
                position: str,
                coord: SkyCoord,
                body_radius: Angle,
                dra: Angle = Angle(0, unit='deg'),
                ddec: Angle = Angle(0, unit='deg'),
                **kwargs):
    """
    Convenience function to place a label relative to an object.

    Parameters
    ----------
    axis : WCSAxes
        The axis in which to plot the ring.
    label: str
        The label.
    position: str
        The label's position as a combination of the vertical and horizontal
        position. For example, 'top right' will work, but 'right top' is a bad
        ideal. Vertical position can be one of 'top', 'center' or 'bottom'.
        Horizontal position can be one of 'left', 'center' or 'right'. If you
        want it directly in the center, use a single 'center' instead of
        'center center'.
    coord : SkyCoord
        The coordinates of the body being labelled. Can be generated with the
        `Planet` method `get_skycoord`.
    body_radius : Angle
        The angular radius of the object.
    dra : Angle, optional
        Angular offset in right ascension.
    ddec : Angle, optional
        Angular offset in declination.
    **kwargs
        Additional Matplotlib `text` keyword arguments like `fontsize`,
        `color`, etc.

    Returns
    -------
    None
        None.
    """

    position = _parse_object_label_position(position)
    ra_offset = body_radius * position['dx']
    dec_offset = body_radius * position['dy']

    offset1 = coord.spherical_offsets_by(-ra_offset, dec_offset)
    offset2 = offset1.spherical_offsets_by(dra, ddec)

    xy = (offset2.ra.degree, offset2.dec.degree)
    axis.annotate(label, xy=xy, xycoords=axis.get_transform('world'),
                  xytext=position['xytext'], textcoords='offset points',
                  ha=position['ha'], va=position['va'],
                  clip_on=True, **kwargs)
