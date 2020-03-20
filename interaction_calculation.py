import numpy as np
import interaction_dependencies as int_dep
import random

def model_computation(a,b,n,r,h,det_x,det_y,det_z, z_energy, mean_sigma_k, mean_sigma_t):
    x0, y0, z0 = int_dep.rect_surface(a,b,n)

    det = int_dep.detector(det_x,det_y,det_z)

    all_phot = []

    for i in range(n):
        history = [np.array([x0[i],y0[i],z0[i]])]
        obj = int_dep.photon_prop(x0[i],y0[i],z0[i],max(z_energy))

        while True:
            curr_int_num = obj.inter_amount
            obj.calc_length(z_energy, mean_sigma_k, mean_sigma_t)
            x,y,z = obj.calc_coor()
            if int_dep.in_or_out(x,y,z,r,h):
                if obj.energy > 0.1:
                    gamma = random.random()
                    if gamma < obj.sigma_norm:
                        if obj.weight > 10e-11:
                            obj.prop_energy()
                            history.append(np.array([x,y,z, obj.zone_num, det.calc_contrib(x,y,z, obj.weight, obj.sigma_t, obj.diff_s, obj.int_s)]))
            if curr_int_num == obj.inter_amount:
                del obj
                break
        all_phot.append(history)

    sum_contributions = np.zeros(len(z_energy)-1)
    interactions = []
    for i in range(n):
        curr_phot = all_phot[i]
        interactions.append(len(curr_phot) - 1)
        for j in range(len(curr_phot)-1):
            sum_contributions[int(curr_phot[j+1][3]) - 1] += curr_phot[j+1][4]
    flux_dens = sum_contributions/n

    print('Most interactions =', max(interactions), '\n')

    inter_stage = n
    for i in range(max(interactions)):
        inter_stage = inter_stage - interactions.count(i)
        print('Phot. survived after inter. #%i ='%(i+1), inter_stage)

    return all_phot, flux_dens
