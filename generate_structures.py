import random

# import my own modules
import molecule

# create class object
m = molecule.Molecule()

# starting coordinates
xyzheader, comment, atomlist, chargelist, xyz = m.read_xyz("xyz/equilibrium.xyz")

# read normal modes
nmfile = "nm/test_normalmodes.txt"
natom = 3
displacements = m.read_nm_displacements(nmfile, natom)

# generate random structures
nstructures = 100
for i in range(nstructures):
    factors = [random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)]
    print(factors)

    displaced_xyz = m.nm_displacer(xyz, displacements, factors)

    fname = "xyz/generated/%s.xyz" % str(i).zfill(3)
    comment = "generated: %s" % str(i).zfill(3)
    m.write_xyz(fname, comment, atomlist, displaced_xyz)


# make Coulomb matrix
# tcm, fcm = m.triangle_cm(chargelist, xyz, dim)
