"""
molecule.
### Read & write xyz files, convert to Z-matrix, Coulomb matrix,
    transform / sample by normal modes, Z-matrix manipulation
### Goals (done: -, not done: x)
- read xyz
- write xyz
x convert to Z-matrix
- convert to Coulomb matrix (CM)
x ability to sort atomlist and xyz by charge (consistently ordered CM)
- reduced CM
x normal mode displacement
x normal mode sampling
x Z-matrix displacement
x Z-matrix sampling
"""
######
import numpy as np

######
class Molecule:
    """methods to manipulate molecules"""

    def __init__(self):
        pass

    def periodicfunc(self, element):
        """Outputs atomic number for each element in the periodic table"""
        with open("pt.txt") as pt_file:
            for line in pt_file:
                if line.split()[0] == element:
                    atomnum = line.split()[1]
                    break
        return int(atomnum)

    # read/write xyz files

    def read_xyz(self, fname):
        """Read a .xyz file"""
        with open(fname, "r") as xyzfile:
            xyzheader = int(xyzfile.readline())
            comment = xyzfile.readline()
            # chargearray = zeros((xyzheader, 1))
            xyzmatrix = np.loadtxt(fname, skiprows=2, usecols=[1, 2, 3])
            atominfoarray = np.loadtxt(fname, skiprows=2, dtype=str, usecols=[0])
            chargearray = [self.periodicfunc(symbol) for symbol in atominfoarray]
        return xyzheader, comment, atominfoarray, chargearray, xyzmatrix

    def write_xyz(self, fname, comment, atoms, xyz):
        """Write .xyz file"""
        natom = len(atoms)
        xyz = np.transpose(xyz)
        atoms_xyz = np.transpose(np.append([atoms], xyz, axis=0))
        with open(fname, "w") as xyzfile:
            np.savetxt(
                fname,
                atoms_xyz,
                fmt='%s',
                delimiter="  ",
                header=str(natom) + "\n" + comment,
                footer="",
                comments="",
            )
        return

    def sort_array(self, tosort, sortbyarray):
        """sort tosort by sortbyarray (have to be same size)"""
        indices = np.argsort(sortbyarray)
        indices = indices[::-1]
        sorted_array = tosort[indices]
        return sorted_array

    # Coulomb matrix

    def triangle_cm(self, charges, xyz, dim):
        """Computes the triangle Coulomb matrix from charges and xyz arrays"""

        tcm = np.zeros((dim, dim))  # the CM of size dim**2
        fcm = np.zeros((dim, dim))  # must make sure to np.zeros; fcm=tcm doesn't work.
        natom = len(charges)  # number of atoms

        for i in range(natom):
            diag_element = 0.5 * charges[i] ** 2.4  # diagonal elements
            tcm[i, i] = diag_element
            fcm[i, i] = diag_element
            for j in range(i + 1, natom):
                dist = np.linalg.norm(xyz[i, :] - xyz[j, :])
                reps = charges[i] * charges[j] / dist  # Pair-wise repulsion
                tcm[i, j] = reps
                fcm[i, j] = reps
                fcm[j, i] = reps  # opposite elements are equal
        return tcm, fcm

    def reduced_cm(self, cm, size):
        """change CM to reduced CM"""
        # only 1st row of CM
        rcm = cm[0:size, 0]
        return rcm
