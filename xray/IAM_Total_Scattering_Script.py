import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from scipy.interpolate import InterpolatedUnivariateSpline

#Load in all the important premade files, and Adjust any tunables
geometry = np.loadtxt('SF6_geometry.txt') #Should be a 3 columned XYZ Array.  Column 0 = X, Column 1 = Y, Column 2 = Z
Scattering_Factor_Coefficients = np.load('Scattering_Factors.npy', allow_pickle=True) #For calculating atomic scattering factors for elastic scattering
Compton_Factors = np.load('Compton_Factors.npy') #For calculating the inelastic scattering
atomic_numbers = np.array([16, 9, 9, 9, 9, 9, 9]) #These should be the atomic numbers corresponding to the rows in the geometry file.  This indicates row 0 is the XYZ coordinates for S, rows 1-6 are the XYZ coordinates for the F atoms
natoms = atomic_numbers.size 
q_bins = 1000
phi_bins = 720

#Define Normalization and other variables
q = np.linspace(0, 10, q_bins, endpoint=True)
q = np.reshape(np.repeat(q, phi_bins), (q_bins, phi_bins))
phi = np.linspace(0, 2*np.pi, phi_bins, endpoint=False)
phi = np.transpose(np.reshape(np.repeat(phi, q_bins), (phi_bins, q_bins)))

    
#This loop calculates the atomic scattering factors from the coefficients, and calculates the inelastic scattering contribution.  The atomic scattering factors will be used when the elastic scattering is calculated later.
q_inelastic = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.5, 2., 8., 15.]) #These go with the Compton factors array- don't change them unless Compton Array changes.
q_inelastic *= 4*np.pi
inelastic_scattering = np.zeros_like(q_inelastic)
scattering_factors = np.zeros((natoms, q_bins, phi_bins))
for i in range(0, natoms):
    atom = atomic_numbers[i] - 1
    factor_coefficient = Scattering_Factor_Coefficients[atom,:]
    inelastic_scattering += Compton_Factors[atom,:]
    atomic_factor = factor_coefficient[0]*np.exp(-factor_coefficient[4]*(q/(4*np.pi))**2) + factor_coefficient[1]*np.exp(-factor_coefficient[5]*(q/(4*np.pi))**2) + factor_coefficient[2]*np.exp(-factor_coefficient[6]*(q/(4*np.pi))**2) + factor_coefficient[3]*np.exp(-factor_coefficient[7]*(q/(4*np.pi))**2) + factor_coefficient[8]
    scattering_factors[i,:,:] = atomic_factor
    
scat_spline = InterpolatedUnivariateSpline(q_inelastic, inelastic_scattering)
Inelastic_Jungfrau_Shaped = scat_spline(q) 

#This part of the script sets up the plotting
#The saveplot function shows what a detector image may look like, without polarization effects accounted for! In real detector images, the X-ray polarization causes anisotropy in the image.
cmap = 'jet'
def saveplot(data,mmin,mmax,name): 
        ax1 = pl.subplots()
        ax1 = pl.subplot(111, polar=True)
        ax1.set_yticklabels([])
        image = ax1.pcolormesh(phi, q, data, cmap=cmap)
        image.set_clim(mmin,mmax)
        pl.colorbar(image)
        pl.savefig('{}.png'.format(name))
        pl.clf()
def saveplot_1D(data1,name):
        pl.figure()
        pl.plot(data1[:,0], data1[:,1])
        pl.savefig('{}.png'.format(name))

#Defining every r vector between 2 atoms- should be N**2 total
r = np.zeros((natoms, 3, natoms))
#This loop calculates the elastic scattering
molecular_contribution = np.zeros((q_bins, phi_bins))
atomic_contribution = np.sum(scattering_factors**2, axis = 0)
debye_image = atomic_contribution
for l in range(0,natoms):
    r[:,:,l] = geometry - geometry[l,:]
for i in range(0,natoms):
    for j in range(i+1, natoms):
        a = scattering_factors[i,:,:]
        b = scattering_factors[j,:,:]
        molecular_contribution = a*b*2*np.sinc(q*np.sqrt(r[j,0,i]**2 + r[j,1,i]**2 + r[j,2,i]**2)/np.pi)
        debye_image += molecular_contribution

debye_image += Inelastic_Jungfrau_Shaped


#The Debye result has been produced- the next session plots and saves the result in various ways.
saveplot(debye_image, 0,4900, 'SF6_test_w_inelastic')  #Save the Debye result as an image
#The next lines save the result as an array
Debye_Result = np.zeros((q[:,0].size,2))
Debye_Result[:,0] = q[:,0]
Debye_Result[:,1] = debye_image[:,0]
saveplot_1D(Debye_Result, 'Debye_example_w_inelastic') #Saves an image
np.save('SF6_Debye_Total.npy', Debye_Result) #Saves as a numpy array to be loaded in and used in other scripts
