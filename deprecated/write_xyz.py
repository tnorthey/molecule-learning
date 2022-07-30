
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

