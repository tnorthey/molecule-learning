import molecule

# create class object
m = molecule.molecule()

def test_read_xyz():
    xyzheader, comment, atomlist, chargelist, xyzmatrix = m.read_xyz("xyz/test.xyz")
    assert xyzheader == 3, "xyzheader should be 3"
    assert comment.__contains__("test"), "comment should be 'test'"
    assert atomlist[0] == "O", "1st atom should be O"
    assert chargelist[0] == 8, "1st atomic charge should be 8"
    assert xyzmatrix[0, 0] == 0.0, "Upper left coordinate should be 0.0"

def test_periodicfunc():
    h = m.periodicfunc("H")
    he = m.periodicfunc("He")
    c = m.periodicfunc("C")
    assert h == 1, "H should have atom number 1"
    assert he == 2, "He should have atom number 2"
    assert c == 6, "C should have atom number 2"

def test_coulombmat():
    dim = 3
    cm = m.coulombmat('xyz/test.xyz', dim)
    # at the moment coulombmat reads an xyz file instead of e.g. an array, I could update
    # it also has 2nd argument "dim" to pad with 0s, this is ok
    assert round(cm[0, 0]) == 74, "upper left element doesn't round to 74"
    assert cm[dim - 1, dim - 1] == 0.5, "bottom right element isn't 0.5"

