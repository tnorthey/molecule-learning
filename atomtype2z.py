import numpy as np

def atomtype2z(atomtyp):
  """ 
  convert atom type (atomtype) i.e. H, He, Li, etc. to atomic number Z
  inputs: list of atom types of length natm
  output: list of atomic numbers of length natm
  """

  natom = len(atomtyp)
  z = np.zeros(natom)
  
  for i in range(0, natom):
      if (atomtype[i] == ' H' or atomtype[i] == 'H '):
         z[i] = 1
      elif (atomtype[i] == ' B' or atomtype[i] == 'B '):
         z[i] = 5
      elif (atomtype[i] == ' C' or atomtype[i] == 'C '):
         z[i] = 6
      elif (atomtype[i] == ' N' or atomtype[i] == 'N '):
         z[i] = 7
      elif (atomtype[i] == ' O' or atomtype[i] == 'O '):
         z[i] = 8
      elif (atomtype[i] == ' P' or atomtype[i] == 'P '):
         z[i] = 15
      elif (atomtype[i] == ' F' or atomtype[i] == 'F '):
         z[i] = 9
      elif (atomtype[i] == ' I' or atomtype[i] == 'I '):
         z[i] = 53
      elif (atomtype[i] == ' K' or atomtype[i] == 'K '):
         z[i] = 19
      elif (atomtype[i] == ' S' or atomtype[i] == 'S '):
         z[i] = 16
      elif (atomtype[i] == ' Y' or atomtype[i] == 'Y '):
         z[i] = 39
      elif (atomtype[i] == ' W' or atomtype[i] == 'W '):
         z[i] = 74
      elif (atomtype[i] == ' U' or atomtype[i] == 'U '):
         z[i] = 92
      elif (atomtype[i] == ' V' or atomtype[i] == 'V '):
         z[i] = 23
      elif (atomtype[i] == 'Si'):
         z[i] = 14
      elif (atomtype[i] == 'Zn'):
         z[i] =30
      elif (atomtype[i] == 'Nb'):
         z[i] =41
      elif (atomtype[i] == 'Mo'):
         z[i] =42
      elif (atomtype[i] == 'Ru'):
         z[i] =44
      elif (atomtype[i] == 'Rh'):
         z[i] =45
      elif (atomtype[i] == 'Pd'):
         z[i] =46
      elif (atomtype[i] == 'Ag'):
         z[i] =47
      elif (atomtype[i] == 'Cd'):
         z[i] =48
      elif (atomtype[i] == 'In'):
         z[i] =49
      elif (atomtype[i] == 'Hf'):
         z[i] =72
      elif (atomtype[i] == 'Ta'):
         z[i] =73
      elif (atomtype[i] == 'Re'):
         z[i] =75
      elif (atomtype[i] == 'Os'):
         z[i] =76
      elif (atomtype[i] == 'Ir'):
         z[i] =77
      elif (atomtype[i] == 'Pt'):
         z[i] =78
      elif (atomtype[i] == 'Au'):
         z[i] =79
      elif (atomtype[i] == 'Pb'):
         z[i] =82
      elif (atomtype[i] == 'Bi'):
         z[i] =83
      elif (atomtype[i] == 'Po'):
         z[i] =84
      elif (atomtype[i] == 'At'):
         z[i] =85
      elif (atomtype[i] == 'Rb'):
         z[i] =37
      elif (atomtype[i] == 'Cs'):
         z[i] =55
      elif (atomtype[i] == 'Fr'):
         z[i] =87
      elif (atomtype[i] == 'Ce'):
         z[i] =58
      elif (atomtype[i] == 'Pr'):
         z[i] =59
      elif (atomtype[i] == 'Nd'):
         z[i] =60
      elif (atomtype[i] == 'Sm'):
         z[i] =62
      elif (atomtype[i] == 'Eu'):
         z[i] =63
      elif (atomtype[i] == 'Gd'):
         z[i] =64
      elif (atomtype[i] == 'Tb'):
         z[i] =65
      elif (atomtype[i] == 'Dy'):
         z[i] =66
      elif (atomtype[i] == 'Ho'):
         z[i] =67
      elif (atomtype[i] == 'Er'):
         z[i] =68
      elif (atomtype[i] == 'Tm'):
         z[i] =69
      elif (atomtype[i] == 'Yb'):
         z[i] =70
      elif (atomtype[i] == 'Lu'):
         z[i] =71
      elif (atomtype[i] == 'Th'):
         z[i] =90
      elif (atomtype[i] == 'Pa'):
         z[i] =91
      elif (atomtype[i] == 'Np'):
         z[i] =93
      elif (atomtype[i] == 'Pu'):
         z[i] =94
      elif (atomtype[i] == 'Am'):
         z[i] =95
      elif (atomtype[i] == 'Cm'):
         z[i] =96
      elif (atomtype[i] == 'Zr'):
         z[i] =40
      elif (atomtype[i] == 'Cr'):
         z[i] = 24
      elif (atomtype[i] == 'Mn'):
         z[i] = 25
      elif (atomtype[i] == 'Ni'):
         z[i] = 28 
      elif (atomtype[i] == 'Cu'):
         z[i] = 29
      elif (atomtype[i] == 'As'):
         z[i] = 33
      elif (atomtype[i] == 'Na'):
         z[i] = 11
      elif (atomtype[i] == 'Be'):
         z[i] = 4
      elif (atomtype[i] == 'Li'):
         z[i] = 3
      elif (atomtype[i] == 'Mg'):
         z[i] = 21
      elif (atomtype[i] == 'Ca'):
         z[i] = 20
      elif (atomtype[i] == 'Co'):
         z[i] = 27
      elif (atomtype[i] == 'Ba'):
         z[i] = 56
      elif (atomtype[i] == 'Sr'):
         z[i] = 38
      elif (atomtype[i] == 'Sc'):
         z[i] = 21
      elif (atomtype[i] == 'Ra'):
         z[i] = 88
      elif (atomtype[i] == 'Ge'):
         z[i] = 32
      elif (atomtype[i] == 'Sn'):
         z[i] = 50
      elif (atomtype[i] == 'Ti'):
         z[i] = 81
      elif (atomtype[i] == 'Ga'):
         z[i] = 31
      elif (atomtype[i] == 'Cl'):
         z[i] = 17
      elif (atomtype[i] == 'Br'):
         z[i] = 35
      elif (atomtype[i] == 'Se'):
         z[i] = 34
      elif (atomtype[i] == 'Al'):
          z[i] = 13
      elif (atomtype[i] == 'Sb'):        
          z[i] = 51
      elif (atomtype[i] == 'Te'):        
          z[i] = 52
      elif (atomtype[i] == 'Fe'):
          z[i] = 26
      else 
          print('Element undefined. Manually add to cmcalc.')
          print(atomtype[i])
          print(i)
    return z
