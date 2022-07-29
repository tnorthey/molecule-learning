class molecule:
    def __init__():
        pass
    
    def read_xyz(fname):
        """
        ##################################################
        # read_xyz:	Reads xyz file 'fname'
        # Inputs:	fname, the file name of the xyz file to read
        # Outputs: 	AtomList (string list), list of atomic labels;
        # 		Coords (float list), coordinates as column vector with the format X1,Y1,Z1,X2,Y2,Z2,...
        ##################################################
        """
        with open(fname, "r") as f:  # open file
            AtomList = []
            Coords = []
            c = 0
            for line in f:
                c += 1
                if c == 1:
                    pass
                elif c == 2:
                    comment = line
                else:
                    AtomList.append(line.split()[0])  # List of atomic labels
                    for i in range(1, 4):
                        Coords.append(float(line.split()[i]))
        ##################################################dd
        return AtomList, Coords, comment
    
    
    def write_xyz(AtomList, Coords, fname, comment):
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


    def cmcalc(natm, atmtyp, xc, yc, zc, cm, rcm, mx, m, rcm_on):
        """
        integer :: i,j,m,natm,mx,c
        character(len=2) :: atmtyp(natm)
        real*8 :: xc(natm),yc(natm),zc(natm),z(natm),dist
        real*8, dimension (natm,natm) :: cm
        real*8, dimension (m) :: rcm
        logical :: rcm_on
        """
    
        # convert atom type list to atomic number Z
        z = atomtype2z(atomtype)
    
        # sort by atomic number (z) then loop over sorted indices to create CM
        # Implement this here <----------
    
        # diagonal elements
        for i in range(natm):
            cm[i, i] = 0.5 * z[i] ** 2.4
    
        # OFF DIAGONALS
        for i in range(natm):
            for j in range(i, natm):
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
        #! CALCULATE REDUCED CM OR MODIFIED FULL CM
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
