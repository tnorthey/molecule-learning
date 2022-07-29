from rw_xyz import read_xyz

def test_read_xyz():
    AtomList, Coords, comment = read_xyz('equilibrium.xyz')
    assert AtomList[0] == "8", "1st atom should be 8"
    assert Coords[0] == 0.0, "1st Coord should be 0.0"

