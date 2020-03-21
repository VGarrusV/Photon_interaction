import numpy as np
import random
from math import sqrt
from scipy.interpolate import interp1d

def rect_surface(a,b,n):
    x0 = np.random.uniform(-b/2,b/2,n)
    y0 = np.random.uniform(-a/2,a/2,n)
    z0 = np.zeros(n)
    return x0,y0,z0

def circ_surface(r,n):
    x0 = np.zeros(n)
    y0 = np.zeros(n)
    for i in range(n):
        while True:
            beta1 = 2*r*random.random() - r
            beta2 = 2*r*random.random() - r
            if (beta1*beta1 + beta2*beta2) < r*r:
                break
        x0[i] = beta1
        y0[i] = beta2
    return x0, y0, np.zeros(n)

def ellipt_surface(a,b,n):
    x0 = np.zeros(n)
    y0 = np.zeros(n)
    for i in range(n):
        while True:
            #beta1 = 2*a*random.random() - a
            #beta2 = 2*b*random.random() - b
            beta1 = 2*a*random.random() - a
            beta2 = b*random.random()
            if beta1*beta1 / (a*a) + beta2*beta2 / (b*b) < 1:
                break
        x0[i] = beta1
        y0[i] = beta2
    return x0,y0, np.zeros(n)

def in_or_out(x,y,z,r,h):
    if sqrt(x*x+y*y) < r and z >= 0 and z < h:
        return True

def create_zones(filename, E_min, E_max, step):
    energ, S_k, S_t = np.array_split(np.loadtxt(filename, dtype=float), 3, axis=1)
    ln = len(energ)
    energ = energ.reshape(1,ln)
    energ = energ[0]
    S_k = S_k.reshape(1,ln)
    S_k = S_k[0]
    S_t = S_t.reshape(1,ln)
    S_t = S_t[0]

    f_t = interp1d(energ, S_t)
    f_k = interp1d(energ, S_k)

    new_energy = np.linspace(E_min,E_max,int((E_max-E_min)*100))
    new_t = f_t(new_energy)
    new_k = f_k(new_energy)

    z_sigma_t = [new_t[0]]
    z_sigma_k = [new_k[0]]
    z_energy = [new_energy[0]]

    i = 1
    while i < len(new_t):
        if (z_sigma_t[-1] - new_t[i]) >= step:
            z_sigma_t.append(new_t[i])
            z_sigma_k.append(new_k[i])
            z_energy.append(new_energy[i])
        i+=1
    z_energy.append(E_max)
    z_sigma_t.append(new_t[-1])
    z_sigma_k.append(new_k[-1])
    return z_energy, z_sigma_t, z_sigma_k

class photon_prop:
    m0c2 = 0.511

    def __init__(self, x, y, z, energy):
        self.x = x
        self.y = y
        self.z = z
        self.energy = energy
        self.inter_amount = 0
        self.weight = 1
        self.sigma_norm = 1

    def calc_length(self, z_energy, mean_sigma_k, mean_sigma_t):
        # for i in range(len(z_energy)-1):
        #     if self.energy > z_energy[i] and self.energy <= z_energy[i+1]:
        #         self.zone_num = i+1
        #         break
        for i in range(len(z_energy)-1):
            if self.energy <= z_energy[-(i+1)] and self.energy >= z_energy[-(i+2)]:
                self.zone_num = len(z_energy) - i - 1
                break

        # sigma_t = (z_sigma_t[self.zone_num-1] + z_sigma_t[self.zone_num])/2
        # sigma_k = (z_sigma_k[self.zone_num-1] + z_sigma_k[self.zone_num])/2
        sigma_t = mean_sigma_t[self.zone_num-1]
        sigma_k = mean_sigma_k[self.zone_num-1]

        gamma = random.random()
        self.length = -np.log(gamma)/sigma_t
        self.sigma_t = sigma_t
        self.sigma_k = sigma_k
        self.sigma_norm = sigma_k/sigma_t
        self.weight = self.sigma_norm * self.weight

    def calc_coor(self):
        while True:
            a = 1 - 2*random.random()
            b = 1 - 2*random.random()
            d = a*a + b*b
            if d <= 1:
                break
        cos_psi = a/sqrt(d)
        sin_psi = b/sqrt(d)
        if self.inter_amount == 0:
            cos_theta = 2*random.random() - 1
            sin_theta = sqrt(1-cos_theta*cos_theta)
            #self.theta = np.random.uniform(0,np.pi/2,1)
            # self.om1 = np.sin(self.theta) * cos_psi
            # self.om2 = np.sin(self.theta) * sin_psi
            # self.om3 = np.cos(self.theta)
            self.om1 = sin_theta * cos_psi
            self.om2 = sin_theta * sin_psi
            self.om3 = cos_theta
        else:
            # self.om1 = np.sin(np.arccos(self.mu)) * cos_psi
            # self.om2 = np.sin(np.arccos(self.mu)) * sin_psi
            # self.om3 = np.cos(np.arccos(self.mu))
            om3 = self.om3*self.mu + sqrt((1-self.mu*self.mu)*(1-self.om3*self.om3))*cos_psi
            om2 = (self.om2*(self.mu - self.om3*om3) + self.om1 * sin_psi *\
            sqrt((1-self.mu*self.mu)*(1-self.om3*self.om3)))/(1-self.om3*self.om3)
            om1 = (self.om1*(self.mu - self.om3*om3) - self.om2 * sin_psi *\
            sqrt((1-self.mu*self.mu)*(1-self.om3*self.om3)))/(1-self.om3*self.om3)
            self.om3 = om3
            self.om2 = om2
            self.om1 = om1

        self.x = self.x + self.om1 * self.length
        self.y = self.y + self.om2 * self.length
        self.z = self.z + self.om3 * self.length
        return self.x, self.y, self.z

    def prop_energy(self):
        self.inter_amount += 1

        alpha0 = self.energy / (photon_prop.m0c2)

        while True:
            gamma1 = random.random()
            gamma2 = random.random()

            x_e = alpha0 * (1+2*alpha0*gamma1) / (1+2*alpha0)
            p = x_e/alpha0 + alpha0/x_e + (1/alpha0 - 1/x_e)*(2+1/alpha0 - 1/x_e)

            if gamma2*(1+2*alpha0+1/(1+2*alpha0)) < p:
                break
        alpha = x_e

        mu = 1 - 1/alpha + 1/alpha0
        self.mu = mu
        self.energy = alpha*photon_prop.m0c2
        self.diff_s = 1/2 * (1+alpha0*(1-mu))**(-2) * (1+mu**2+alpha0**2*\
        (1-mu)**2/(1+alpha0*(1-mu)))
        self.int_s = 2*np.pi*((1+alpha0)/alpha0**2 *(2*(1+alpha0)/(1+2*alpha0)-\
        np.log(1+2*alpha0)/alpha0)+np.log(1+2*alpha0)/(2*alpha0) - (1+3*alpha0)/(1+2*alpha0)**2)

class detector:
    def __init__(self, det_x, det_y, det_z):
        self.det_x = det_x
        self.det_y = det_y
        self.det_z = det_z

    def calc_contrib(self, x, y, z, weight, sigma_t, diff_s, int_s):
        self.dist = sqrt((self.det_x - x)**2 + (self.det_y - y)**2 + (self.det_z - z)**2)
        self.eta = weight * np.exp(-sigma_t*self.dist) * diff_s / (self.dist**2 * int_s)
        return self.eta
