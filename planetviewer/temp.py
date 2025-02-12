from planetviewer import update_spice_kernels, SolarSystemBody, make_custom_observer
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u
import astropy.constants as c
import spiceypy as spice


satellite = SolarSystemBody('Pluto')

make_custom_observer('Home', 34.1547658820689*u.deg, -118.12965108301688*u.deg,
                     800*u.imperial.ft, 1701)

time = Time('2021-06-08 12:30')
print(satellite.get_apparent_epoch(time, 'Home'))
