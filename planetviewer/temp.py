import matplotlib.pyplot as plt
from astropy.time import Time

from funkyfresh import set_style

from planetviewer import make_wcs, Planet, load_spice_kernels, \
    set_standard_axis_limits, set_standard_axis_labels

load_spice_kernels(download=False)

set_style('AAS', silent=True)

observer = 'Keck'
time = Time('2021-06-08 12:44')

jupiter = Planet('Jupiter')
ganymede = Planet('Ganymede')

center = jupiter.get_skycoord(time, observer)
fov = jupiter.parse_fov(time, observer, 8)
wcs = make_wcs(center, fov)

fig, axis = plt.subplots(figsize=(5.2, 5), subplot_kw=dict(projection=wcs),
                         layout='constrained', clear=True)

jupiter.draw(axis, time, observer, rings='Main')
ganymede.draw(axis, time, observer)

set_standard_axis_limits(axis)
set_standard_axis_labels(axis)

plt.savefig(f'/Users/zachariahmilby/Desktop/test.pdf', dpi=200)
plt.close(fig)
