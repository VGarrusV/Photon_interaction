import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from scipy.interpolate import interp1d

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

#Рассчет длины свободного пробега.
    def calc_length(self):
        sigma_t = f1(self.energy)
        sigma_k = f2(self.energy)
        gamma = np.random.uniform(0,1,1)
        self.length = -np.log(gamma)/sigma_t
        self.sigma_t = sigma_t
        self.sigma_k = sigma_k
        self.sigma_norm = sigma_k/sigma_t
        self.weight = self.sigma_norm * self.weight
#Подсчет новых координат.
    def calc_coor(self):
        self.psi = np.random.uniform(0,2*np.pi,1)
        if self.inter_amount == 0:
            self.theta = np.random.uniform(0,np.pi/2,1)
            om1 = np.sin(self.theta) * np.cos(self.psi)
            om2 = np.sin(self.theta) * np.sin(self.psi)
            om3 = np.cos(self.theta)
        else:
            om1 = np.sin(np.arccos(self.mu)) * np.cos(self.psi)
            om2 = np.sin(np.arccos(self.mu)) * np.sin(self.psi)
            om3 = np.cos(np.arccos(self.mu))

        self.x = self.x + om1 * self.length
        self.y = self.y + om2 * self.length
        self.z = self.z + om3 * self.length
#Запись координат в историю если необходимо.
    def write_hist(self):
        self.hist.append([self.x, self.y, self.z])
#Подсчет новой энергии.
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
#Подсчет вклада в детектор.
    def calc_contrib(self):
        self.dist = np.sqrt((det_x - self.x)**2 + (det_y - self.y)**2 + (det_z - self.z)**2)
        self.eta = self.weight * np.exp(-self.sigma_norm*self.dist) * self.diff_s / (self.dist**2 * self.int_s)
        self.contr_hist.append(self.eta[0])

#Координаты детектора
det_x = 0
det_y = 0
det_z = 10

#параметры пов. испускания
a = 5
b = 30
n = 100

#Начвальные координаты испускания фотонов
x0 = np.random.uniform(-b/2,b/2,n)
x0 = x0.reshape(n,1)
y0 = np.random.uniform(-a/2,a/2,n)
y0 = y0.reshape(n,1)
z0 = np.zeros(n)
z0 = z0.reshape(n,1)

#Считывание энергий и сечений из файла
energ, S_k, S_t = np.array_split(np.loadtxt("energy2.txt", dtype=float), 3, axis=1)
energ = energ.reshape(1,82)
energ = energ[0]
S_k = S_k.reshape(1,82)
S_k = S_k[0]
S_t = S_t.reshape(1,82)
S_t = S_t[0]

#Интерполяция сечений
f1 = interp1d(energ, S_t)
f2 = interp1d(energ, S_k)

#Создание фотонов и прохождение через взаимодействия. Жуткие ifы.
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

#Подсчет вклада по первому взаимодействию
sum_contribution = 0
layer = 1
for i in range(n):
    if len(all_phot[i].contr_hist) > 0:
        sum_contribution += all_phot[i].contr_hist[layer-1]
print(sum_contribution)

#Вывод на экран
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
        #ax.plot(xs,ys,zs, color = colors[j])
        #ax.plot(xs,ys,zs, 'g')

        ax.plot(xs,ys,zs,color =  (0.2 + j*0.1, 0.5 - j*0.05, 0.2))
        ax.plot(all_phot[i].hist[j+1][0],all_phot[i].hist[j+1][1],all_phot[i].hist[j+1][2], '.', color = 'C' + '%i'%(j+1))

ax.plot(circ_x, circ_y1, 'b', alpha = 0.4)
ax.plot(circ_x, circ_y2, 'b', alpha = 0.4)
ax.plot(circ_x, circ_y1, 10, 'b', alpha = 0.4)
ax.plot(circ_x, circ_y2, 10, 'b', alpha = 0.4)
ax.plot(x0,y0, '.')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
