import chemcoord as cc
#from chemcoord.xyz_functions import get_rotation_matrix

water = cc.Cartesian.read_xyz('test.xyz', start_index=1)
water = cc.Cartesian.read_xyz('test.xyz')

print(water)
coords = water[['x', 'y', 'z']]
print('coords')
print(coords)

# or explicit label based indexing
print(water.loc[:, 'x'])
# or explicit integer based indexing
print(water.iloc[:, 1])

zmat = water.get_zmat()
print(zmat)
