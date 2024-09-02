import numpy as np
import spiceypy as spice
from astropy.coordinates import Angle, SkyCoord

from planetviewer.spice_kernels import load_spice_kernels

load_spice_kernels()

# the time doesn't really matter here since these frames aren't time-dependent
et = 0
b1950_to_j2000 = spice.pxform('FK4', 'J2000', et)
j2000_to_b1950 = spice.pxform('J2000', 'FK4', et)
axis, angle = spice.raxisa(j2000_to_b1950)

offsets = []

for body in ['Jupiter', 'Saturn', 'Uranus', 'Neptune (invariable pole)',
             'Neptune (planet pole)']:
    print(body)

    # get RA/Dec from SPICE PCK
    code = spice.bodn2c(body.split(' ')[0])
    ra_j2000, dec_j2000, _, _ = spice.bodeul(code, et)
    ra_j2000 = Angle(ra_j2000, unit='rad')
    if body == 'Uranus':  # need angular momentum pole!
        ra_j2000 -= Angle(180, unit='deg')
        dec_j2000 = -dec_j2000
    dec_j2000 = Angle(dec_j2000, unit='rad')

    # from Jacobson (2009)
    if body == 'Neptune (invariable pole)':
        ra_j2000 = Angle(299.46, unit='deg')
        dec_j2000 = Angle(43.40, unit='deg')

    # calculate angular momentum pole unit vector
    j2000_pole_vector = np.array([np.cos(ra_j2000) * np.cos(dec_j2000),
                                  np.sin(ra_j2000) * np.cos(dec_j2000),
                                  np.sin(dec_j2000)])

    # rotate J2000 pole unit vector to B1950
    b1950_pole_vector = spice.mxv(j2000_to_b1950, j2000_pole_vector)

    # calculate pole RA/Dec in B1950 frame
    dec_b1950 = Angle(np.arcsin(b1950_pole_vector[2]), unit='rad')
    ra_b1950 = Angle(np.arctan2(b1950_pole_vector[1], b1950_pole_vector[0]),
                     unit='rad')

    # determine ascending node in B1950 frame (node is 90° ahead of RA)
    offset = Angle(90, unit='deg')
    b1950_ascending_node_vector = np.array([np.cos(ra_b1950 + offset),
                                            np.sin(ra_b1950 + offset),
                                            0])

    # rotate B1950 ascending node into J2000 frame
    j2000_ascending_node_vector = spice.mxv(
        b1950_to_j2000, b1950_ascending_node_vector)

    # calculate node J2000 RA/Dec
    dec_j2000_node = Angle(np.arcsin(j2000_ascending_node_vector[2]),
                           unit='rad')
    ra_j2000_node = Angle(np.arctan2(j2000_ascending_node_vector[1],
                                     j2000_ascending_node_vector[0]),
                          unit='rad')

    # rotate B1950 node in J2000 frame to planetary equatorial coordinates
    # where x-axis points toward the J2000 ascending node
    vec = spice.rotvec(j2000_ascending_node_vector, (ra_j2000 + offset).rad, 3)
    vec = spice.rotvec(vec, (offset - dec_j2000).rad, 1)

    # calculate the angular offset between x-axis (J2000 ascending node) and
    # the B1950 ascending node vector
    dlon = Angle(np.arctan2(vec[1], vec[0]), unit='rad')

    # save to list
    offsets.append(f'{body},{dlon.degree}')

    # Print out to compare with Phil Nicholson's results; they differ slightly
    # but that is because of my choice of RA/Dec. If I use his they match
    # perfectly, but I'm not sure where his came from since they are more
    # precise than any published values I've found.
    print(f'   B1950 pole vector:  {b1950_pole_vector[0]:.8f}  {b1950_pole_vector[1]:.8f}  {b1950_pole_vector[2]:.8f}')
    print(f'   RA, DEC (B1950):  {ra_b1950.degree:.6f}  {dec_b1950.degree:.6f}')
    print(f'   B1950 asc node vector:  {b1950_ascending_node_vector[0]:.8f}  {b1950_ascending_node_vector[1]:.8f}  {b1950_ascending_node_vector[2]:.8f}')
    print('')
    print(f'   J2000 pole vector:  {j2000_pole_vector[0]:.8f}  {j2000_pole_vector[1]:.8f}  {j2000_pole_vector[2]:.8f}')
    print(f'   RA, DEC (J2000):  {ra_j2000.degree:.6f}  {dec_j2000.degree:.6f}')
    print('')
    print(f'   precessed B1950 node:  {j2000_ascending_node_vector[0]:.8f}  {j2000_ascending_node_vector[1]:.8f}  {j2000_ascending_node_vector[2]:.8f}')
    print(f'   RA, DEC (J2000):  {ra_j2000_node.degree:.6f}  {dec_j2000_node.degree:.6f}')
    print('')
    print(f'   J2000_long of B1950 node:  {dlon.degree:.8f}')
    print('\n')

# save offsets to file
out = 'planet,offset_deg' + '\n'
out += '\n'.join(offsets)
with open('b1950-j2000_ring_longitude_offsets.dat', 'w') as file:
    file.write(out)
