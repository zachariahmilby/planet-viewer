from planetviewer import set_kernel_path, furnish_spice_kernels, SolarSystemBody
from astropy.time import Time

set_kernel_path('/Users/zachariahmilby/SPICE')
furnish_spice_kernels()

satellite = SolarSystemBody('Europa')

time = Time('2021-06-12 19:30')
print(time.isot, satellite.get_occultation('Jupiter', time, 'Maunakea'))
