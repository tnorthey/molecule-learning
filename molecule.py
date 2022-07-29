import numpy as np

class molecule:
    def __init__(self):
        pass

    def periodicfunc(self, element):
        """ Outputs atomic number for each element in the periodic table """
        with open("pt.txt") as f:
            for line in f:
                if line.split()[0] == element:
                    atomnum = line.split()[1]
                    break
        return int(atomnum)


    def read_xyz(self, fname):
        """Read a .xyz file """
        with open(fname, 'r') as xyzfile:
            xyzheader = int(xyzfile.readline())
            comment = xyzfile.readline()
            #chargearray = zeros((xyzheader, 1))
            xyzmatrix = np.loadtxt(fname, skiprows=2, usecols=[1, 2, 3])
            atominfoarray = np.loadtxt(fname, skiprows=2, dtype=str, usecols=[0])
            chargearray = [self.periodicfunc(symbol) for symbol in atominfoarray]
        return xyzheader, comment, atominfoarray, chargearray, xyzmatrix


    def write_xyz(self, AtomList, Coords, fname, comment):
        """
        ##################################################
        # write_xyz: Write xyz file 'fname' using atom list 'AtomList'
        #            and Cartesian coordinates 'Coords'
        # Inputs:	AtomList (string list), list of atomic labels
        # 		Coords (float list), coordinates as column vector 
        #           with the format X1, Y1, Z1, X2, Y2, Z2,...
        # 		fname (str), output xyz file name
        # 		comment (str), comment line in the xyz file
        # Outputs:	None, but it creates the fname.xyz file
        ##################################################
        """
        Nat = len(AtomList)  # Number of atoms
        with open(fname, "w") as f:  # Open file for writing
            f.write(
                str(Nat) + "\n"
            )  # The first line of an xyz file contains only the number of atoms
            f.write(
                comment + "\n"
            )  # The second line is blank or contains a title string (convention)
            for i in range(3 * Nat):  # Loop over vectors of lentgh 3*Nat
                if i % 3 == 0:  # Indices 0,3,6,...
                    Atom = AtomList[
                        int(i / 3)
                    ]  # Read atom labels (at indices 0,1,2,... every i=0,3,6,...)
                    x = Coords[i]  # x coordinate
                elif (i - 1) % 3 == 0:  # Indices 1,4,7,...
                    y = Coords[i]  # y coordinate
                elif (i - 2) % 3 == 0:  # Indices 2,5,8,...
                    z = Coords[i]  # z coordinate
                    # f.write( Atom + '  ' + str(x) + '  ' + str(y) + '  ' + str(z) + '\n')	# Write out Lines containing atom labels and x,y,z coordinates
                    f.write(
                        "%2s %12.8f %12.8f %12.8f \n" % (Atom, x, y, z)
                    )  # Write out Lines containing atom labels and x,y,z coordinates
        ##################################################
        return

    def cmcalc(AtomList, Coords):
        """
        CALCULATE REDUCED CM OR MODIFIED FULL CM
        integer :: i,j,m,natm,mx,c
        character(len=2) :: atmtyp(natm)
        real*8 :: xc(natm),yc(natm),zc(natm),z(natm),dist
        real*8, dimension (natm,natm) :: cm
        real*8, dimension (m) :: rcm
        logical :: rcm_on
        """
        natom = len(AtomList)

        xc = Coords
        # convert atom type list to atomic number Z
        z = atomtype2z(AtomList)

        # sort by atomic number (z) then loop over sorted indices to create CM
        # Implement this here <----------

        # diagonal elements
        for i in range(natom):
            cm[i, i] = 0.5 * z[i] ** 2.4

        # OFF DIAGONALS
        for i in range(natom):
            for j in range(i, natom):
                dist = (
                    (xc[i] - xc[j]) ** 2 + (yc[i] - yc[j]) ** 2 + (zc[i] - zc[j]) ** 2
                ) ** 0.5
                if dist > 1.0e-8:
                    cm[i, j] = z[i] * z[j] / dist
                else:
                    cm[i, j] = 0.0
                # other side of diagonal is equivalent
                cm[j, i] = cm[i, j]

        """ 
        106 format(f12.8)
        !OPEN(7,FILE='rcm.dat',FORM="FORMATTED",STATUS="replace")
        !write(7,*) 'CM: ',idx
        if (rcm_on) then
          rcm = cm(1:m,1)
        else
          c=0
          do i=1,mx-1
            do j=i+1,mx
              c=c+1
              rcm(c) = cm(i,j)
            enddo
          enddo
        endif
        rcm(1) = 0.0d0  ! first element is absorbing atom Z^2, set to 0 (for reduced or full CM)
      
        !do i=1,m
          !write(7,106) rcm[i]
        !enddo
        !close(7)
        """

        return cm


    def coulombmat(self, fname, dim):
        """
        This function takes in an xyz input file for a molecule, number of atoms in the biggest molecule  to computes the corresponding coulomb Matrix
        """
        xyzheader, comment, atominfoarray, chargearray, xyzmatrix = self.read_xyz(fname)
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
