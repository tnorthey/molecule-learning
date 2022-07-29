import molecule

# create class object
m = molecule.molecule()

def test_read_xyz():
    AtomList, Coords, comment = m.read_xyz("xyz/equilibrium.xyz")
    assert AtomList[0] == "8", "1st atom should be 8"
    assert Coords[0] == 0.0, "1st Coord should be 0.0"

def test_periodicfunc():
    h = m.periodicfunc("H")
    he = m.periodicfunc("He")
    c = m.periodicfunc("C")
    assert h == 1, "H should have atom number 1"
    assert he == 2, "He should have atom number 2"
    assert c == 6, "C should have atom number 2"

"""
def test_cmcalc():
    AtomList, Coords, comment = m.read_xyz('linear.xyz')
    print(AtomList)
    print(Coords)
    cm = m.cmcalc(AtomList, Coords)

test_cmcalc()

"""
