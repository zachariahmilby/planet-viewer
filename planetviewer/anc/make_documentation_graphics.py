import matplotlib.pyplot as plt
from astropy.coordinates import Angle
from astropy.time import Time

from planetviewer import (load_spice_kernels, Planet, make_wcs,
                          set_standard_axis_limits, set_standard_axis_labels,
                          get_available_observatories, get_available_spacecraft)

load_spice_kernels(download=False)

print(get_available_spacecraft())
print(get_available_observatories())
raise SystemExit(7)

time = Time('2021-06-08 14:00')
observer = 'Triton'
planet = Planet('Neptune')
center = planet.get_skycoord(time=time, observer=observer)
fov = planet.parse_fov(time=time, observer=observer, fov=Angle(30, unit='deg'))
wcs = make_wcs(center=center, fov=fov)
fig, axis = plt.subplots(subplot_kw=dict(projection=wcs))
set_standard_axis_limits(axis)
set_standard_axis_labels(axis)
plt.grid()

plt.savefig('docs/axis0.png')

planet.draw(axis=axis, time=time, observer=observer)
plt.savefig('docs/axis1.png')

satellites = ['Proteus', 'Larissa', 'Galatea', 'Despina', 'Thalassa']
for moon in satellites:
    satellite = Planet(moon)
    satellite.draw(axis=axis, time=time, observer=observer)
plt.savefig('docs/axis2.png')

transform = axis.get_transform('world')
for moon in satellites:
    satellite = Planet(moon)
    coord = satellite.get_skycoord(time=time, observer=observer)
    axis.annotate(moon, xy=(coord.ra.deg, coord.dec.deg), xycoords=transform,
                  annotation_clip=True)
plt.savefig('docs/axis3.png')
