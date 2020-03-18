import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from scipy.interpolate import interp1d
import time

class photon_prop:
    #class variable here
    m0c2 = 0.511

    def __init__(self, x0, y0, z0, E):
        self.x = x0
        self.y = y0
        self.z = z0
        self.energy = E
        self.inter_amount = 0
        self.hist = []
        self.contr_hist = []
        self.weight = 1
        self.sigma_norm = 1

    def calc_length(self):
        for i in range(len(z_energy)-1):
            if self.energy > z_energy[i] and self.energy <= z_energy[i+1]:
                self.zone_num = i+1

        sigma_t = (z_sigma_t[self.zone_num-1] + z_sigma_t[self.zone_num])/2
        sigma_k = (z_sigma_k[self.zone_num-1] + z_sigma_k[self.zone_num])/2

        gamma = np.random.uniform(0,1,1)
        self.length = -np.log(gamma)/sigma_t
        self.sigma_t = sigma_t
        self.sigma_k = sigma_k
        self.sigma_norm = sigma_k/sigma_t
        self.weight = self.sigma_norm * self.weight

    def calc_coor(self):
        self.psi = np.random.uniform(0,2*np.pi,1)
        if self.inter_amount == 0:
            self.theta = np.random.uniform(0,np.pi/2,1)
            self.om1 = np.sin(self.theta) * np.cos(self.psi)
            self.om2 = np.sin(self.theta) * np.sin(self.psi)
            self.om3 = np.cos(self.theta)
        else:
            self.om1 = np.sin(np.arccos(self.mu)) * np.cos(self.psi)
            self.om2 = np.sin(np.arccos(self.mu)) * np.sin(self.psi)
            self.om3 = np.cos(np.arccos(self.mu))

        self.x = self.x + self.om1 * self.length
        self.y = self.y + self.om2 * self.length
        self.z = self.z + self.om3 * self.length

    def write_hist(self):
        self.hist.append([self.x, self.y, self.z])

    def prop_energy(self):
        self.inter_amount += 1

        alpha0 = self.energy / (photon_prop.m0c2)
        i=0
        while i < 1:
            gamma1 = np.random.uniform(0,1,1)
            gamma2 = np.random.uniform(0,1,1)

            x_e = alpha0 * (1+2*alpha0*gamma1) / (1+2*alpha0)
            p = x_e/alpha0 + alpha0/x_e + (1/alpha0 - 1/x_e)*(2+1/alpha0 - 1/x_e)

            if gamma2*(1+2*alpha0+1/(1+2*alpha0)) < p:
                alpha = x_e
                i += 1
        mu = 1 - 1/alpha + 1/alpha0
        self.mu = mu
        self.energy = alpha*photon_prop.m0c2
        self.diff_s = 1/2 * (1+alpha0*(1-mu))**(-2) * (1+mu**2+alpha0**2*(1-mu)**2/(1+alpha0*(1-mu)))
        self.int_s = 2*np.pi*((1+alpha0)/alpha0**2 *(2*(1+alpha0)/(1+2*alpha0)-np.log(1+2*alpha0)/alpha0)+np.log(1+2*alpha0)/(2*alpha0) - (1+3*alpha0)/(1+2*alpha0)**2)

    def calc_contrib(self):
        self.dist = np.sqrt((det_x - self.x)**2 + (det_y - self.y)**2 + (det_z - self.z)**2)
        self.eta = self.weight * np.exp(-self.sigma_norm*self.dist) * self.diff_s / (self.dist**2 * self.int_s)
        self.contr_hist.append([self.eta[0], self.zone_num])

energ, S_k, S_t = np.array_split(np.loadtxt("energy2.txt", dtype=float), 3, axis=1)
ln = len(energ)
energ = energ.reshape(1,ln)
energ = energ[0]
S_k = S_k.reshape(1,ln)
S_k = S_k[0]
S_t = S_t.reshape(1,ln)
S_t = S_t[0]

f_t = interp1d(energ, S_t)
f_k = interp1d(energ, S_k)

new_energy = np.linspace(0.1,3,290)
new_t = f_t(new_energy)
new_k = f_k(new_energy)

z_sigma_t = [new_t[0]]
z_sigma_k = [new_k[0]]
z_energy = [new_energy[0]]

i = 1
while i < len(new_t):
    if (z_sigma_t[-1] - new_t[i]) >= 0.005:
        z_sigma_t.append(new_t[i])
        z_sigma_k.append(new_k[i])
        z_energy.append(new_energy[i])
    i+=1
z_energy.append(3.0)
z_sigma_t.append(new_t[-1])
z_sigma_k.append(new_k[-1])

det_x = 0
det_y = 0
det_z = 11

a = 5
b = 30
n = 10000

x0 = np.random.uniform(-b/2,b/2,n)
x0 = x0.reshape(n,1)
y0 = np.random.uniform(-a/2,a/2,n)
y0 = y0.reshape(n,1)
z0 = np.zeros(n)
z0 = z0.reshape(n,1)

print('Calculating interactions for %i photons.\n'%n)
start_time = time.time()
all_phot = []
for i in range(n):
    obj = photon_prop(x0[i], y0[i], z0[i], 3)
    all_phot.append(obj)
    gamma = np.random.uniform(0,1,1)
    while np.sqrt(obj.x[0]**2 + obj.y[0]**2) < 20 and obj.z[0] < 10 and obj.z[0] >= 0 and obj.energy > 0.1 and gamma <= obj.sigma_norm and obj.weight > 10e-11:
        obj.write_hist()
        obj.calc_length()
        obj.calc_coor()
        if np.sqrt(obj.x[0]**2 + obj.y[0]**2) < 20 and obj.z[0] < 10 and obj.z[0] >= 0 and obj.energy > 0.1 and gamma <= obj.sigma_norm and obj.weight > 10e-11:
            obj.prop_energy()
            obj.calc_contrib()
print("--- Done in %s seconds ---" % (time.time() - start_time))

sum_contributions = np.zeros(len(z_energy)-1)
interactions = []
for i in range(n):
    interactions.append(len(all_phot[i].contr_hist))
    if len(all_phot[i].contr_hist) > 0:
        for j in range(len(all_phot[i].contr_hist)):
            sum_contributions[all_phot[i].contr_hist[j][1]-1] += all_phot[i].contr_hist[j][0]
flux_dens = sum_contributions/n

print('Most interactions =', max(interactions), '\n')

inter_stage = n
for i in range(max(interactions)):
    inter_stage = inter_stage - interactions.count(i)
    print('Phot. survived after inter. #%i ='%(i+1), inter_stage)

#--------------------Plotting------------------------------------------------------------------
#Figure 1.
switch = input("\nPlot the hedgehog? : ")
if switch == 'y':
    print('Plotting.')
    circ_x = np.linspace(-20,20,100)
    circ_y1 = np.sqrt(400 - circ_x**2)
    circ_y2 = -np.sqrt(400 - circ_x**2)

    fig = plt.figure(1)

    ax = fig.gca(projection='3d')

    colors = ['green','mediumseagreen', 'orange', 'red', 'darkred']
    for i in range(n):
        grn = 0.5
        for j in range(len(all_phot[i].hist) - 1):
            xs = np.array([all_phot[i].hist[j][0][0], all_phot[i].hist[j+1][0][0]])
            ys = np.array([all_phot[i].hist[j][1][0], all_phot[i].hist[j+1][1][0]])
            zs = np.array([all_phot[i].hist[j][2][0], all_phot[i].hist[j+1][2][0]])

            r = 0.2 + j*0.1
            if r > 1:
                r = 1
            g = 0.5 - j*0.05
            if g > 1:
                g = 1
            ax.plot(xs,ys,zs,color =  (r, g, 0.2))
            ax.plot(all_phot[i].hist[j+1][0],all_phot[i].hist[j+1][1],all_phot[i].hist[j+1][2], '.', color = 'C' + '%i'%(j+1))

    ax.plot(circ_x, circ_y1, 'b', alpha = 0.4)
    ax.plot(circ_x, circ_y2, 'b', alpha = 0.4)
    ax.plot(circ_x, circ_y1, 10, 'b', alpha = 0.4)
    ax.plot(circ_x, circ_y2, 10, 'b', alpha = 0.4)
    ax.plot(np.array([det_x]), np.array([det_y]), np.array([det_z]), '*')
    ax.plot(x0,y0, '.')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
else:
    print('Plotting next Figure.')

#Figure 2.------------------------------------------------------------------------
plt.figure(2)
plt.clf()
plt.grid(linestyle='--')
widths = []
centers = []
for i in range(len(z_energy)-1):
    widths.append(z_energy[i+1] - z_energy[i])
    centers.append((z_energy[i+1] + z_energy[i])/2)

plt.bar(z_energy[0:-1], list(flux_dens) ,align = 'edge', width = widths, edgecolor = 'r', alpha = 0.7)
plt.plot(centers, flux_dens, color = 'y', linestyle = '-.')
plt.plot(centers, flux_dens, '.', color = 'g')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlim(0,3)
plt.xlabel('Energy, MeV')
plt.ylabel('Flux density')
plt.title('Contribution vs Energy.')
plt.show()
#Figure 3.------------------------------------------------------------------------
plt.figure(3)
plt.clf()
plt.grid(linestyle = '--')
plt.plot(z_energy,z_sigma_t, label = 'Sigma Total.')
plt.plot(z_energy,z_sigma_k, '--', label = 'Sigma Kompt.')
for i in range(len(z_energy)):
    plt.axvline(z_energy[i], color = 'r', linewidth = '0.5', alpha = 0.7)
plt.xlim(0.05,3.05)
plt.legend()
plt.xlabel('Energy, MeV')
plt.ylabel('Sigma')
plt.title('Energy zones.')
plt.show()
