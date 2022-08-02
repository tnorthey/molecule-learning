import numpy as np


class Xray:
    def __init__(self):
        pass

    def atomic_factor(self, atom_number, qvector):
        """ returns atomic x-ray scattering factor for atom_number, and qvector """
        # np.array([[1, 2], [3, 4]])
        aa = np.array(
            [
                [0.489918, 0.262003, 0.196767, 0.049879],  # hydrogen
                [0.8734, 0.6309, 0.3112, 0.1780],  # helium
                [1.1282, 0.7508, 0.6175, 0.4653],  # lithium
                [1.5919, 1.1278, 0.5391, 0.7029],  # berylium
                [2.0545, 1.3326, 1.0979, 0.7068],  # boron
                [2.3100, 1.0200, 1.5886, 0.8650],  # carbon
                [12.2126, 3.1322, 2.0125, 1.1663], # nitrogen
                [3.0485, 2.2868, 1.5463, 0.8670],  # oxygen
                [3.5392, 2.6412, 1.5170, 1.0243],  # fluorine
                [3.9553, 3.1125, 1.4546, 1.1251],  # neon
                [4.7626, 3.1736, 1.2674, 1.1128],  # sodium
                [5.4204, 2.1735, 1.2269, 2.3073],  # magnesium
                [6.4202, 1.9002, 1.5936, 1.9646],  # aluminium
                [6.2915, 3.0353, 1.9891, 1.5410],  # Siv
                [6.4345, 4.1791, 1.7800, 1.4908],  # phosphorus
                [6.9053, 5.2034, 1.4379, 1.5863],  # sulphur
            ]
        )

        bb = np.array(
            [
                [20.6593, 7.74039, 49.5519, 2.20159],  # hydrogen
                [9.1037, 3.3568, 22.9276, 0.9821],  # helium
                [3.9546, 1.0524, 85.3905, 168.261],  # lithium
                [43.6427, 1.8623, 103.483, 0.5420],  # berylium
                [23.2185, 1.0210, 60.3498, 0.1403],  # boron
                [20.8439, 10.2075, 0.5687, 51.6512],  # carbon
                [0.00570, 9.8933, 28.9975, 0.5826],  # nitrogen
                [13.2771, 5.7011, 0.3239, 32.9089],  # oxygen
                [10.2825, 4.2944, 0.2615, 26.1476],  # fluorine
                [8.4042, 3.4262, 0.2306, 21.7184],  # Ne
                [3.2850, 8.8422, 0.3136, 129.424],  # Na
                [2.8275, 79.2611, 0.3808, 7.1937],  # Mg
                [3.0387, 0.7426, 31.5472, 85.0886],  # Al
                [2.4386, 32.3337, 0.6785, 81.6937],  # Siv
                [1.9067, 27.1570, 0.5260, 68.1645],  # P
                [1.4679, 22.2151, 0.2536, 56.1720],
            ]
        )  # S

        cc = np.array(
            [
                0.001305,  # hydrogen
                0.0064,  # helium
                0.0377,  # lithium
                0.0385,  # berylium
                -0.1932,  # boron
                0.2156,  # carbon
                -11.529,  # nitrogen
                0.2508,  # oxygen
                0.2776,  # fluorine
                0.3515,
                0.6760,
                0.8584,
                1.1151,
                1.1407,
                1.1149,
                0.8669,
            ]
        )

        qlen = len(qvector)
        atomfactor = np.zeros(qlen)
        for j in range(qlen):
            for i in range(4):
                atomfactor[j] += aa[atom_number, i] * np.exp(
                    -bb[atom_number, i] * (0.25 * qvector[j] / np.pi) ** 2
                )
        atomfactor += cc[atom_number]
        return atomfactor

    def iam_calc(self, atoms, xyz, qlength):
        """ calculate IAM scattering curve for molecule atoms, xyz """
        sfactors = np.load("xray/Scattering_Factors.npy", allow_pickle=True)
        print(sfactors)
        # to do
        qvector = 0
        iam = 0
        return qvector, iam

    """
    Nq=length(q)],        # length of q
Nat=size(Atoms,1)],   # number of atoms
ZZ=Atoms(:,1)],       # atomic number
X=Atoms(:,2)],
Y=Atoms(:,3)],
Z=Atoms(:,4)],

Imol=zeros(1,Nq)],
Iat=Imol],
for j1=1:Nat
    for j2=1:Nat
        rij=sqrt((X(j1)-X(j2))^2+(Y(j1)-Y(j2))^2+(Z(j1)-Z(j2))^2)],
        for i=1:Nq
            brac=q(i)*rij],
            if brac>1e-9
               Imol(i) = Imol(i) + formfact(ZZ(j1),q(i)).*formfact(ZZ(j2),q(i)).*sin(brac)./brac],
            else
               Iat(i) = Iat(i) + formfact(ZZ(j1),q(i)).*formfact(ZZ(j2),q(i))],
            end
        end
    end
end
Ideb=Iat+Imol],



function[fq]=formfact(an,q)
# Obtain form-factor from:
# http://lamp.tu-graz.ac.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

# input atom number (an) and vector q in Angstroms
# outputs form factor for atom number and q-range

au2ang = 0.52917721092d0], # convert Ang coordinates to au at the end
iam_constants],  # H to F constants, if needed obtain more from website.
q=q/au2ang],

Nq=length(q)],
fq=zeros(Nq,1)],
for j=1:Nq
    for i=1:4
        fq(j)=fq(j)+aa(an,i)*exp(-bb(an,i)*(.25*q(j)/pi)^2)],
    end
end
fq=fq+cc(an)],
"""
