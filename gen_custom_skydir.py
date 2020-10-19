import os

import numpy as np
import matplotlib.pyplot as plt

def find_nearest(array, value):

    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

def P_ani(phi, theta, alpha):

    norm_const = 1/(4*np.pi*(1+alpha**2))
    b = 1+(alpha*2*np.sqrt(3)*np.cos(theta))+(3*alpha**2*(np.cos(theta)**2))
    return norm_const*b

def gen_custom_sample_sphere(n_samples, alpha, savetofile = False, fname = 'none'):


    phi = np.linspace(-np.pi, np.pi, 100)
    theta = np.linspace(-np.pi/2, np.pi/2, 100)
    PHI, THETA = np.meshgrid(phi, theta)

    theta_out = []
    phi_out = []

    #P_ANI = P_ani(PHI,THETA,alpha)

    while len(theta_out)<n_samples:
        u = np.random.uniform(-1,1)
        v = np.random.uniform(0,1)
        
        p_test = np.pi*u
        t_test = np.arccos(2*v-1)
        t_test = t_test#-np.pi/2

        
        weight_test= P_ani(p_test,t_test,alpha)#P_ANI[th_idx][ph_idx]

        q = np.random.uniform(0,1)
        if weight_test>q:
            #conversion for mollweide in theta
            theta_out.append(t_test-np.pi/2)
            phi_out.append(p_test)
        else:
            pass

    angles_out = np.column_stack([phi_out, theta_out])
    return angles_out

#example of usage:
if __name__ == '__main__':
	#look at the underlying distribution
	# fig, ax = plt.subplots(subplot_kw=dict(projection='mollweide'), figsize=(10,8))
	# im = ax.pcolormesh(X, Y, P_ANI, vmin = 0, vmax = 0.4)
	# ax.grid()
	# fig.colorbar(im, orientation='horizontal');
	# plt.savefig('dipole.png')

    #generate some sky directions to use
    angles= gen_custom_sample_sphere(1000, 0.5, savetofile = False, fname='anglesout_alpha05.txt')

    #plot the sky positions you have generated 
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.scatter(angles[:,0], angles[:,1], marker='+', s=50)
    ax.grid(True)
    plt.show()
    #plt.savefig('uniform.png')

