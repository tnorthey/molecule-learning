"""
molecule.
### Read & write xyz files, convert to Z-matrix, Coulomb matrix,
    transform / sample by normal modes, Z-matrix manipulation
### Goals:
✔ read xyz
- write xyz
- convert to Z-matrix
✔ convert to Coulomb matrix (CM)
- ordered CM
- reduced CM
- normal mode displacement
- normal mode sampling
- Z-matrix displacement
- Z-matrix sampling
"""

import numpy as np


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

    def coulombmat(self, fname, dim):
        """Reads xyz file, computes Coulomb matrix"""

        xyzheader, _, _, chargearray, xyzmatrix = self.read_xyz(fname)
        cij = np.zeros((dim, dim))

        for i in range(xyzheader):
            for j in range(xyzheader):
                if i == j:
                    cij[i, j] = (
                        0.5 * chargearray[i] ** 2.4
                    )  # Diagonal term described by Potential energy of isolated atom
                else:
                    dist = np.linalg.norm(xyzmatrix[i, :] - xyzmatrix[j, :])
                    cij[i, j] = (
                        chargearray[i] * chargearray[j] / dist
                    )  # Pair-wise repulsion
        return cij
