import sys
sys.path.append('../')
import molecule as m

def test_read_xyz():
    AtomList, Coords, comment = m.read_xyz('equilibrium.xyz')
    assert AtomList[0] == "8", "1st atom should be 8"
    assert Coords[0] == 0.0, "1st Coord should be 0.0"

