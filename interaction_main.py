import numpy as np
import interaction_calculation as int_calc
import interaction_dependencies as int_dep
import interaction_plotting as int_plt
import time

#This version currently supports only rectangular birth surfaces.
a = 5 #rect. width
b = 30 #rect. length
n = 100 #amount of photons to be modelled.
#This version currently supports only cylindrical volumes for interaction.
rad = 20 #radius of the cylinder
heigth = 10 #heigth of the cylinder

#This version currently support only a single detector.
det_x = 0
det_y = 0
det_z = 11

E_min = 0.1 #Minimum energy allowed for a photon.
E_max = 3 #Initial energy for the photons.
step = 0.005 #Step for dividing energy into energy zones.
filename = 'energy2.txt' #File from which to read energies and interaction cross-sections.

#Execute energy zone creation.
print('Creating energy zones...')
z_energy, z_sigma_t, z_sigma_k = int_dep.create_zones(filename, E_min, E_max, step)

#Execute model computation. all_phot contains history of all interactions, flux_dens
#contains information of contributions to the detector.
start_time = time.time()
print('Executing model computation...')
all_phot, flux_dens = int_calc.model_computation(a, b, n, rad, heigth, det_x, det_y, det_z, z_energy, z_sigma_t, z_sigma_k)
print("\n--- Done in %s seconds ---\n" % (time.time() - start_time))
#Plotting the results.
print('Plotting flux density...')
int_plt.plot_contribution(flux_dens, z_energy)

usr_inp = input('Plot the hedgehog? (y/[n]): ')
if usr_inp == 'y':
    int_plt.plot_hedgehog(r, all_phot)

print('Done.')
