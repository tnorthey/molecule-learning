"""
molecule, T. Northey, 2022
### Read & write xyz files, convert to Z-matrix, Coulomb matrix,
    transform / sample by normal modes, Z-matrix manipulation
### Goals (done: -, not done: x)
- read xyz
- write xyz
x convert to Z-matrix
- convert to Coulomb matrix (CM)
- ability to sort atomlist and xyz by charge (consistently ordered CM)
- reduced CM
- normal mode displacement
x normal mode sampling
x Z-matrix displacement
x Z-matrix sampling
"""
######
import numpy as np
import pandas as pd

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
                    return int(line.split()[1])

    def atomic_mass(self, element):
        """Outputs atomic mass for each element in the periodic table"""
        with open("atomic_masses.txt") as am_file:
            for line in am_file:
                if line.split()[0] == element:
                    return line.split()[1]

    # read/write xyz files

    def read_xyz(self, fname):
        """Read a .xyz file"""
        with open(fname, "r") as xyzfile:
            xyzheader = int(xyzfile.readline())
            comment = xyzfile.readline()
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
                fmt="%s",
                delimiter=" ",
                header=str(natom) + "\n" + comment,
                footer="",
                comments="",
            )
        return

    # Bagel stuff
    def write_bagel_dyson(self, xyzfile, outfile="bagel_inp.json"):
        """writes a bagel dyson norm input file based on 
        bagel_dyson.template with atoms and geometry from xyzfile"""
        _, _, atoms, _, xyzmatrix = self.read_xyz(xyzfile)
        bagel_df = pd.read_json(
            "templates/bagel_dyson.template"
        )  # read bagel input template as pandas dataframe
        for k in range(len(atoms)):
            bagel_df["bagel"][0]["geometry"][k]["atom"] = atoms[k]
            bagel_df["bagel"][0]["geometry"][k]["xyz"] = xyzmatrix[k, :]
        bagel_df.to_json(outfile, indent=4)  # this runs in bagel!
        return

    def read_bagel_dyson(self, bagel_dyson_output, max_rows):
        """read dyson norms and ionisation energies from a bagel dyson output file"""
        str_find = "Norms^2 of Dyson orbitals approximately indicate the strength of an inization transitions."
        energy, norm = [], []  # define here to avoid return error if str_find isn't found
        with open(bagel_dyson_output, "r") as f:
            for line in f:
                if str_find in line:    # go to line containing str
                    out_array = np.loadtxt(     # numpy loadtxt into an array
                        f,
                        dtype={
                            "names": ("from", "-", "to", "energy", "norm"),
                            "formats": ("i4", "a2", "i4", "f4", "f4"),
                        },
                        skiprows = 4,
                        max_rows = max_rows
                    )
                    energy = out_array["energy"]
                    norm = out_array["norm"]
        return energy, norm
    ### End Bagel stuff 

    def sort_array(self, tosort, sortbyarray):
        """sort tosort by sortbyarray (have to be same size)"""
        indices = np.argsort(sortbyarray)
        indices = indices[::-1]
        sorted_array = tosort[indices]
        return sorted_array

    ### distances array

    def distances_array(self, xyz):
        """Computes matrix of distances from xyz"""
        natom = xyz.shape[0]  # number of atoms
        dist_array = np.zeros((natom, natom))  # the array of distances
        for i in range(natom):
            dist_array[i, i] = 0
            for j in range(i + 1, natom):
                dist = np.linalg.norm(xyz[i, :] - xyz[j, :])
                dist_array[i, j] = dist
                dist_array[j, i] = dist  # opposite elements are equal
        return dist_array


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

    ### Normal modes section
    def read_nm_displacements(self, fname, natoms):
        """ read_nm_displacements: Reads displacement vector from file=fname e.g. 'normalmodes.txt'
        Inputs: 	natoms (int), total number of atoms
        Outputs:	displacements, array of displacements, size: (nmodes, natoms, 3) """
        if natoms == 2:
            nmodes = 1
        elif natoms > 2:
            nmodes = 3 * natoms - 6
        else:
            print("ERROR: natoms. Are there < 2 atoms?")
            return False
        with open(fname, "r") as xyzfile:
            tmp = np.loadtxt(fname)
        displacements = np.zeros((nmodes, natoms, 3))
        for i in range(3 * natoms):
            for j in range(nmodes):
                if i % 3 == 0:  # Indices 0,3,6,...
                    dindex = int(i / 3)
                    displacements[j, dindex, 0] = tmp[i, j]  # x coordinates
                elif (i - 1) % 3 == 0:  # Indices 1,4,7,...
                    displacements[j, dindex, 1] = tmp[i, j]  # x coordinates
                elif (i - 2) % 3 == 0:  # Indices 2,5,8,...
                    displacements[j, dindex, 2] = tmp[i, j]  # x coordinates
        return displacements

    def displace_xyz(self, xyz, displacement, factor):
        """displace xyz by displacement * factor
        xyz and displacement should be same size"""
        return xyz + displacement * factor

    def nm_displacer(self, xyz, displacements, modes, factors):
        """displace xyz along all displacements by factors array"""
        natoms = xyz.shape[0]
        summed_displacement = np.zeros(displacements[0, :, :].shape)
        c = 0
        for mode in modes:
            summed_displacement += displacements[mode, :, :] * [factors][c]
            c += 1
        displaced_xyz = self.displace_xyz(xyz, summed_displacement, 1)
        return displaced_xyz

    def animate_mode(self, mode, xyz_start_file, nmfile, natoms):
        """make xyz file animation along normal mode"""
        displacements = self.read_nm_displacements(nmfile, natoms)
        a = 0.4
        factor = np.linspace(-a, a, 20, endpoint=True)
        factor = np.append(factor, np.linspace(a, -a, 20, endpoint=True))
        _, _, atoms, _, xyz_start = self.read_xyz(xyz_start_file)
        for k in range(len(factor)):
            xyz = self.nm_displacer(xyz_start, displacements, [mode], factor[k])
            xyzfile_out = "animate/mode%i_%s.xyz" % (mode, str(k).zfill(2))
            self.write_xyz(xyzfile_out, str(factor[k]), atoms, xyz)

    ### End normal modes section
