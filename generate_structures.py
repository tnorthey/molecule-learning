import numpy as np
# import my own modules
import molecule

# create class object
m = molecule.Molecule()

# starting coordinates
xyzheader, comment, atomlist, chargelist, xyz = m.read_xyz("xyz/equilibrium.xyz")

# read normal modes
nmfile = "nm/test_normalmodes.txt"
natom = 3
nmodes = 3
a = 0.2
displacements = m.read_nm_displacements(nmfile, natom)

# generate random structures
nstructures = 100
for i in range(nstructures):
    #mass = mi[ int( (i - 1) / 3 ) ]  # int rounds down by default
    #displacement_constant = (mass**.5 * 0.172*freqcm1**.5)**-1
    factors = np.random.rand(nmodes)*2*a - a 
    print(factors)

    displaced_xyz = m.nm_displacer(xyz, displacements, factors)

    fname = "xyz/generated/%s.xyz" % str(i).zfill(3)
    comment = "generated: %s" % str(i).zfill(3)
    m.write_xyz(fname, comment, atomlist, displaced_xyz)

