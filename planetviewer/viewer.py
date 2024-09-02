# import astropy.units as u
# import matplotlib.pyplot as plt
# import numpy as np
# from astropy.coordinates import Angle
# from astropy.time import Time
# from astropy.visualization.wcsaxes.core import WCSAxes
#
# from planetviewer.objects import Planet
# from planetviewer.plotting import _parse_fov, make_wcs, \
#     set_standard_axis_labels, set_standard_axis_limits, plot_limb, \
#     shade_nightside, plot_latitudes, plot_longitudes, plot_rings, \
#     place_label, add_pixel, place_time_offset_label, convert_to_relative_axis
#
#
# def _sort_bodies_by_distance(bodies: list[Planet],
#                              time: Time,
#                              observer: str) -> list[Planet]:
#     """
#     Sort a list of Planet objects by their distance from furthest to closest.
#
#     Parameters
#     ----------
#     bodies
#     time
#     observer
#
#     Returns
#     -------
#     The sorted list of bodies.
#     """
#     distances = np.zeros(len(bodies))
#     for i, body in enumerate(bodies):
#         distances[i] = body.get_distance(time, observer).value
#     sort = np.flip(distances.argsort())
#     return [bodies[i] for i in sort]
#
#
# class LocationNotFound(Exception):
#     """Raised if observer location not available."""
#     pass
#
#
# class OneSheet:
#
#     def __init__(self):
#         pass
#
#
# class Viewer:
#     """
#
#     """
#
#     def __init__(self, target: str, fov: Angle or u.Quantity, observer: str):
#         self._target = Planet(target)
#         self._fov = fov
#         self._observer = observer
#
#     def plot(self,
#              time: Time,
#              additional_targets: list[str] = None,
#              rings: list[str] = None,
#              figsize: tuple = None,
#              fontsize: int or float = None) -> tuple[plt.Figure, WCSAxes]:
#
#         center = self._target.get_skycoord(time, self._observer)
#
#         if additional_targets is None:
#             additional_targets = []
#         targets = [self._target]
#         targets.extend([Planet(i) for i in additional_targets])
#         sorted_targets = _sort_bodies_by_distance(targets, time,
#                                                   self._observer)
#
#         fov = _parse_fov(self._target, time, self._observer, self._fov)
#         wcs = make_wcs(center, fov)
#         fig, axis = plt.subplots(figsize=figsize,
#                                  subplot_kw={'projection': wcs})
#         axis.set_aspect('equal')
#
#         zoffset = 0
#         for target in sorted_targets:
#             args = dict(body=target, time=time,
#                         observer=self._observer, axis=axis,
#                         zoffset=zoffset, relative_body=self._target)
#             plot_limb(**args)
#             shade_nightside(**args)
#             plot_latitudes(**args)
#             plot_longitudes(**args)
#             plot_rings(rings=rings, **args)
#
#             zoffset += 100
#
#         axis.zorder = zoffset + 100
#         set_standard_axis_labels(axis)
#         set_standard_axis_limits(axis)
#
#         return fig, axis
#
#
# class TimeSeries:
#     """
#     Show a series of additional targets at discrete times relative to a central
#     target. For example, if you wanted to show Io transiting Jupiter at
#     5-minute intervals, this is the class for you.
#     """
#
#     def __init__(self,
#                  target: str,
#                  target_display_time: Time,
#                  fov: Angle or u.Quantity,
#                  observer: str,
#                  dra: Angle = None,
#                  ddec: Angle = None):
#         self._target = Planet(target)
#         self._target_display_time = target_display_time
#         self._fov = fov
#         self._observer = observer
#         self._dra = dra
#         self._ddec = ddec
#
#     def plot(self,
#              times: list[Time],
#              additional_targets: list[str] = None,
#              rings: list[str] = None,
#              figsize: tuple = None,
#              fontsize: int or float = None,
#              label_times: bool = True) -> tuple[plt.Figure, WCSAxes]:
#         center = self._target.get_skycoord(
#             self._target_display_time, self._observer)
#         center = center.spherical_offsets_by(self._dra, self._ddec)
#
#         if additional_targets is None:
#             additional_targets = []
#         targets = [self._target]
#         targets.extend([Planet(i) for i in additional_targets])
#         sorted_targets = _sort_bodies_by_distance(
#             targets, self._target_display_time, self._observer)
#
#         fov = _parse_fov(
#             self._target, self._target_display_time, self._observer, self._fov)
#         wcs = make_wcs(center, fov)
#         fig, axis = plt.subplots(figsize=figsize,
#                                  subplot_kw={'projection': wcs})
#         axis.set_aspect('equal')
#
#         zoffset = 0
#         for target in sorted_targets:
#             if target.name == self._target.name:
#                 args = dict(body=self._target, time=self._target_display_time,
#                             observer=self._observer, axis=axis,
#                             zoffset=zoffset)
#                 plot_limb(**args)
#                 shade_nightside(**args)
#                 plot_latitudes(**args)
#                 plot_longitudes(**args)
#                 plot_rings(rings=rings, **args)
#                 zoffset += 100
#             else:
#                 for time in times:
#                     args = dict(body=target, time=time,
#                                 observer=self._observer, axis=axis,
#                                 zoffset=zoffset, relative_body=self._target,
#                                 relative_time=self._target_display_time)
#                     plot_limb(**args)
#                     shade_nightside(**args)
#                     plot_latitudes(**args)
#                     plot_longitudes(**args)
#                     plot_rings(rings=rings, **args)
#                     if label_times:
#                         place_label(label_type='time',
#                                     label_position='upper center',
#                                     fontsize=fontsize, time_fmt='%H:%M',
#                                     **args)
#                     zoffset += 100
#         axis.zorder = zoffset + 100
#         set_standard_axis_labels(axis)
#
#         convert_to_relative_axis(center, axis)
#         set_standard_axis_limits(axis)
#
#         return fig, axis
