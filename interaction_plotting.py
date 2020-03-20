import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

def plot_contribution(flux_dens, z_energy):
    plt.figure(1)
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

def plot_hedgehog(r, all_phot):
    circ_x = np.linspace(-r,r,100)
    circ_y1 = np.sqrt(r*r - circ_x**2)
    circ_y2 = -np.sqrt(r*r - circ_x**2)

    fig = plt.figure(2)

    ax = fig.gca(projection='3d')

    for i in range(len(all_phot)):
        curr_phot = all_phot[i]
        for j in range(len(curr_phot) - 1):
            xs = np.array([curr_phot[j][0], curr_phot[j+1][0]])
            ys = np.array([curr_phot[j][1], curr_phot[j+1][1]])
            zs = np.array([curr_phot[j][2], curr_phot[j+1][2]])

            r = 0.2 + j*0.1
            if r > 1:
                r = 1
            g = 0.5 - j*0.05
            if g > 1:
                g = 1
            ax.plot(xs,ys,zs,color =  (r, g, 0.2))
            ax.plot(np.array([curr_phot[j+1][0]]), np.array([curr_phot[j+1][1]]), np.array([curr_phot[j+1][2]]), '.', color = 'C' + '%i'%(j+1))
        ax.plot(np.array([curr_phot[0][0]]), np.array([curr_phot[0][1]]), np.array([curr_phot[0][2]]), '.', color = 'blue')

    ax.plot(circ_x, circ_y1, 'b', alpha = 0.4)
    ax.plot(circ_x, circ_y2, 'b', alpha = 0.4)
    ax.plot(circ_x, circ_y1, 10, 'b', alpha = 0.4)
    ax.plot(circ_x, circ_y2, 10, 'b', alpha = 0.4)

    #ax.plot(, np.array([det_y]), np.array([det_z]), '*')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    plt.show()
