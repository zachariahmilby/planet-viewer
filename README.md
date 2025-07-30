# `planetviewer`: Ephemeris and Visualization Utilities for Solar System Observers 

`planetviewer` is a Python-based planetary geometry viewer based on the
[JPL Horizons Ephemeris Service](https://ssd.jpl.nasa.gov/horizons/) and the 
NASA PDS Ring-Moon Systems Node 
[Planet Viewers](https://pds-rings.seti.org/tools/). It incorporates the 
functionality of both Matplotlib and Astropy as detailed below. In particular,
it requires the use of Astropy convenience classes like `Time` for times and 
dates, `Angle` (and the related classes `Longitude` and `Latitude`) for angular
quantities and `Quantity` for other general physical quantities with units like
distances, velocities, etc. These help to force users to be cognizant of the 
units of inputs they provide. Most of the functions and methods are decently 
described in their docstrings, but I'll provide some overviews and a few 
examples here.

>**CAUTION**<br>
> Because of all of the Matplotlib transforms and plotting, this can be 
> significantly slower than the PDS Planet Viewers.

## Installation
Here are some installation instructions for the average Anaconda user; if 
you're more advanced I'm sure you can figure it out from here. (Note: in the
instructions below I will assume that you are using a virtual environment named 
`myenv`.) I've tested this using Python 3.10.
1. Activate your virtual environment:<br>
    `% conda activate myenv`
2. Install the `planetviewer` package and its dependencies:<br>
    `% python -m pip install git+https://github.com/zachariahmilby/planet-viewer.git`

You're now ready to use the `planetviewer` package!

## SPICE Kernels

`planetviewer` use NASA's SPICE system for ephemeris calculations. The first 
time you create a `SolarSystemBody` (see "Ephemeris Calculations" section 
below) the package will download all data files (called "kernels") it needs to
make calculations. Some kernels are frequently updated, particularly for 
spacecraft with evolving orbits like JWST and HST. To update the kernels, 
import and call the function `update_spice_kernels()`.

You can also furnish additional kernels using the `furnsh()` function available
in `spiceypy`.

## Observers
Observers must be provided as a string. The choice of observer/observatory can 
be one of three types:
1. A planetary body defined in SPICE (major planets, minor planets, 
   satellites or asteroids),
2. A spacecraft, or
3. A ground-based observatory on Earth.

The lists below detail which of each type are currently available. I haven't 
figured out comets yet.

Available planets/satellites/asteroids:
- See the "NAMES" column in the "Barycenters", "Planets and Satellites", 
  "Asteroids" sections of 
  https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html.

Available spacecraft:
- `JWST`: James Webb Space Telescope
- `HST`: Hubble Space Telescope
- `Galileo`
- `Juno`
- `Cassini`
- `Voyager 1`
- `Voyager 2`

Available ground-based Earth observatories:
- `ALMA`: Atacama Large Millimeter/submillimeter Array
- `Apache Point`: Apache Point Observatory, New Mexico, United States
- `Cahill Roof`: Roof of Cahill Center for Astronomy and Astrophysics, Pasadena, California, United States
- `Keck`: W. M. Keck Observatory, Hawaii, United States
- `Kitt Peak`: Kitt Peak National Observatory, Arizona, United States
- `Lowell`: Lowell Observatory, Arizona, United States
- `Maunakea`: W. M. Keck Observatory, Hawaii, United States
- `McDonald`: McDonald Observatory, Texas, United States
- `Palomar`: Palomar Observatory, California, United States
- `Sommers Bausch`: Sommers-Bausch Observatory, Colorado, United States
- `VLT`: Very Large Telescope, Chile

You can also make a custom observer for anywhere on Earth using the function
`make_custom_observer()` which will then be available as an observer using the
name provided.

## Times
Observation epochs (the time of an observation at the observer) must be 
provided as Astropy `Time` objects.

## Ephemeris Calculations
The fundamental object in `planetviewer` is `SolarSystemBody`. You can choose 
from any Solar System planet, dwarf planet or satellite known to the SPICE 
system. To create a `SolarSystemBody` object, simply provide its name:

```
from planetviewer import SolarSystemBody
planet = SolarSystemBody('Jupiter')
```

This object has a variety of built-in methods which allow you to calculate 
observer-dependent ephemeris information for the planet. A lot of this can of 
course be done using JPL Horizons (and/or its implementation in Astroquery). 
However, a part of my porting of the PDS Ring-Moon Systems Node Planet Viewers 
was recreating all of the various ephemeris information provided in the output.
These methods are also often used internally to produce the visualizations 
detailed below.

> **NOTE**<br>
> I don't believe in west longitudes. Neither does SPICE. All longitudes 
> returned by methods or otherwise reported from this package are east 
> longitudes. I wish I could change the definition of the north pole on planets
> like Venus and Uranus to the angular momentum pole, but I'll save that for 
> another time. As a Boulder mom's bumper sticker might say, "be the change you 
> want to see in the world."

### Small Bodies
By default, the SPICE system does not contain most of the minor bodies in the
Solar System, like comets and asteroids. However, you can get the necessary 
SPICE data for them from JPL Horizons using the function `get_small_body`. For
example, to get the ephemeris for the third interstellar object I3/Atlas
(called C/2025 N1 officially) discovered during summer 2025, use the following: 
```
from planetviewer import get_small_body
get_small_body('C/2025 N1')
```

You can now create a `SolarSystemObject` for C/2025 N1. You only need to call 
this function once; once the ephemeris has been downloaded it will be available
for all future uses of `planetviewer` within this virtual environment. However, 
if the ephemeris has been updated since you last downloaded it, you can pass 
the option `update=True`. You can also specify time ranges for which to 
download the ephemeris data.

> **CAUTION**<br>
> Don't include parentheses in the names of small bodies. For instance, if you
> look up "C/2025 N1" on the 
> [Small-Body Database Lookup web app](https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html),
> the name will return as "C/2025 N1 (ATLAS)", however if you try to search 
> using that name, it will tell you your search string contains one or more bad
> characters. So just use "C/2025 N1" instead.

Some small bodies will have multiple ephemerides, so if you specify the object
by name or by NAIF ID number (the `spk_id` kwarg), the program may exit and 
tell you that you need to specify a record number instead. Choose from those
listed in the terminal output, pass that to `record_no` and try again.

### Ephemeris Calculation Methods
The following table lists all of the ephemeris methods available in the 
`SolarSystemBody` object. These should calculate all of the information (and 
more) returned by the various PDS planet viewers. See each method's docstring 
for specifics on arguments and outputs.

Almost every method requires the fundamental arguments `time` (as an Astropy 
`Time` object) and `observer` (as a string). Other methods have additional 
arguments detailed in their docstrings.

| Method                                     | Description                                                                                                 |
|:-------------------------------------------|:------------------------------------------------------------------------------------------------------------|
| `get_skycoord`                             | Get J2000 sky coordinate (RA/Dec).                                                                          |
| `get_offset`                               | Get offset between this object and another `SkyCoord` object.                                               |
| `get_radec_rates`                          | Get rate at which target's Ra/Dec is changing.                                                              |
| `get_altaz`                                | Get alt/az coordinates for a given observer.                                                                |
| `get_altaz_rates`                          | Get rate at which target's alt/az is changing.                                                              |
| `get_sub_observer_latitude`                | Get apparent sub-observer latitude on the planet.                                                           |
| `get_sub_observer_longitude`               | Get apparent sub-observer longitude on the planet.                                                          |
| `get_subsolar_latitude`                    | Get apparent sub-solar latitude on the planet.                                                              |
| `get_subsolar_longitude`                   | Get apparent sub-solar longitude on the planet.                                                             |
| `get_phase_angle`                          | Get the phase angle at the sub-observer point.                                                              |
| `get_distance`                             | Get the distance from a SPICE ephemeris object to the object's barycenter.                                  |
| `get_relative_velocity`                    | Get the relative velocity of the object.                                                                    |
| `get_ring_subsolar_latitude`               | Get the sub-solar latitude on the rings (the same as the planet's sub-solar latitude).                      |
| `get_ring_subsolar_latitude_range`         | Get the range in sub-solar latitude on the rings due to the apparent angular size of the Sun.               |
| `get_ring_plane_opening_angle`             | Get the ring plane opening angle (the same as the planet's sub-observer latitude).                          |
| `determine_if_rings_illuminated`           | Determine if the rings appear illuminated to the observer.                                                  |
| `get_ring_center_phase_angle`              | Get the ring phase angle (the same as the phase angle at the sub-observer point).                           |
| `get_ascending_node_longitude`             | Get the planetographic longitude of the ascending node of the body's equator on the J2000 frame equator.    |
| `get_ring_subsolar_longitude`              | Get the sub-solar longitude measured relative to the ring plane ascending node.                             |
| `get_ring_sub_observer_longitude`          | Get the sub-observer longitude measured relative to the ring plane ascending node.                          |
| `get_sun_target_distance`                  | Get the distance from the Sun to the object's barycenter.                                                   |
| `get_observer_target_distance`             | Get the distance from the observer to the object's barycenter.                                              |
| `get_light_travel_time`                    | Get the light travel time between a SPICE ephemeris object and the object's barycenter.                     |
| `get_apparent_epoch`                       | Get the apparent epoch at the object as seen by the observer.                                               |
| `get_necessary_epoch`                      | Get the observer epoch necessary to observe the target at a given epoch.                                    |
| `get_occultation`                          | Determine if this object is occulted by another object or vice-versa.                                       |
| `get_eclipsed`                             | Determine if this body casts a shadow on another body or vice-versa.                                        |
| `find_eclipses`                            | Find when object is eclipse (in the shadow) of another body for a given time range.                         |
| `get_latlonalt_sky_coordinates`            | Get the apparent RA/Dec of a latitude/longitude point at a given altitude above the object's surface.       |
| `get_longitude_line_coordinates`           | Get the apparent RA/Dec of an entire meridian.                                                              |
| `get_latitude_line_coordinates`            | Get the apparent RA/Dec of an entire parallel.                                                              |
| `get_limb_sky_coordinates`                 | Get the apparent RA/Dec of an object's limb.                                                                |
| `get_terminator_sky_coordinates`           | Get the apparent RA/Dec of an object's terminator.                                                          |
| `get_dayside_sky_coordinates`              | Get the apaprent RA/Dec of an object's dayside disk (the region between the dayside limb and terminator).   |
| `get_shadow_intersection_sky_coordinates`  | Get the apparent RA/Dec of the intersection of another object's shadow on the apparent disk of this object. |
| `get_angular_radius`                       | Convert the object's mean equatorial radius to angular radius on the sky.                                   |
| `get_angular_offset_between_times`         | Calculate the RA/Dec offset necessary to shift this object's sky coordinates between two times.             |
| `parse_fov`                                | Convert an angle, length or scalar to an angular field-of-view.                                             |
| `check_for_shadows`                        | Determine if a body casts a shadow on another body or vice-versa.                                           |
| `draw`                                     | A convenience function to draw a planet in an Astropy `WCSAxes` axis.                                       |

## Visualization with Matplotlib
The primary purpose of this package is to provide a convenient way to calculate 
ephemeris information without resorting to JPL Horizons (which can be slow and 
inaccurate) and also visualize the appearance of a Solar System object(s) to an 
observer at a particular time. I imagine most people will use this to generate 
figures for papers or presentations. Astropy provides a variety of convenient 
extensions to Matplotlib with produce excellent astronomical graphics, and the 
methods I've created here are designed to interface with both Astropy's 
extensions.

`SolarSystemBody` objects can be drawn in parts (using the various functions 
described below) or using the convenient `draw` method. Rings can be drawn for 
Outer Solar System planets, and Neptune's Adams ring arcs can also be drawn. 
Available rings can be found in the `rings` property of the `SolarSystemBody` 
class.

The World Coordinate System (WCS) projection and corresponding transforms are 
the core of the visualization system. Several functions conveniently calculate 
these for you (see below).

Ths visualization methods will also incorporate any changes you make to 
Matplotlib's runtime configuration (rcParams), so you can style the graphics 
however you'd like.
> **CAUTION**<br>
> There is a [known bug](https://github.com/astropy/astropy/issues/15344) 
> between Astropy's extensions to Matpltolib and pgf rendering. This will cause 
> an exception if you try to use pgf. Hopefully this will be fixed soon.

### Plotting Functions
The following table lists each of the plotting functions available. Their 
individual docstrings provide more detailed instructions for their use.

| Function                        | Description                                                                                      |
|:--------------------------------|:-------------------------------------------------------------------------------------------------|
| `plot_limb`                     | Draw the planet's apparent limb.                                                                 |
| `plot_disk`                     | Draw the planet's apparent disk as a filled polygon.                                             |
| `plot_nightside`                | Draw the planet's apparent night side as a filled polyon.                                        |
| `plot_primary_shadow`           | If the planet is a moon, draw the primary planet's shadow if it intersects with the moon's disk. |
| `plot_latlon`                   | Draw a latitude or longitude line.                                                               |
| `plot_ring`                     | Draw a particular ring.                                                                          |
| `place_ring_pericenter_markers` | Place markers at eccentric ring pericenters (applies only to some rings of Uranus).              |
| `plot_arc`                      | Draw ring arc (applies only to Neptune's Adams ring).                                            |
| `place_label`                   | Place a label at a given coordinate.                                                             |
| `set_standard_axis_labels`      | Set standard J2000 axis labels.                                                                  |
| `set_standard_axis_limits`      | Ensures axis limits match the chosen center and FOV.                                             |
| `convert_to_relative_axis`      | Convert axis from absolute RA/Dec to relative angle from the center.                             |
| `sort_by_distance`              | Sort a list of ephemeris objects by their distance from furthest to closest.                     |
