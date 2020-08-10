import os

import numpy as np
import matplotlib.pyplot as plt

def find_nearest(array, value):
    """
        Finds the nearest value to one given in an array
    """
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

def P_ani(phi, theta, alpha):
"""
	Custom sky distribution from which to draw samples

    Parameters
    ----------
    phi : np.meshgrid
            grid of phi values over which to calculate distribution
    theta: np.meshgrid
           grid of theta values over which to calculate distribution
    alpha: float
           dipole-y-ness
            
    Returns
    -------
    numpy ndarray
        meshgrid of probabilities that a source occurs in given direction
"""
    norm_const = 1/(4*np.pi*(1+alpha**2))
    b = 1+(alpha*2*np.sqrt(3))+(3*alpha**2*(np.cos(theta)**2))
    return norm_const*b

def gen_custom_sample_sphere(n_samples, alpha, savetofile = False, fname = 'none'):
    """
    generate a user-defined random number of points on the celestial sphere

    Parameters
    ----------
    n_samples : int
            the number of samples to generate
    alpha: float
    		alpha for the distribution
    savetofile: Bool
            save to file in current directory, true or false
    fname: string
            filename to save as, 'none' if savetofile = False
            
    Returns
    -------
    numpy ndarray
        table of RA and DEC in radians
    """

	x = np.linspace(-np.pi, np.pi, 100)
	y = np.linspace(-np.pi/2, np.pi/2, 100)
	X, Y = np.meshgrid(x, y)

	# Spherical coordinate arrays derived from x, y
	# Necessary conversions to get Mollweide right
	phi = x.copy()    # physical copy
	phi[x < 0] = 2 * np.pi + x[x<0]
	theta = np.pi/2 - y
	PHI, THETA = np.meshgrid(phi, theta)

	theta_out = []
	phi_out = []

	P_ANI = P_ani(PHI,THETA,alpha)

	while len(theta_out)<n_samples:
	    u = np.random.uniform(-1,1)
	    v = np.random.uniform(0,1)
	    
	    p_test = np.pi*u
	    t_test = np.arccos(2*v-1)
	    t_test = t_test-np.pi/2
	    t_in.append(t_test)
	    p_in.append(p_test)
	    
	    coord = [p_test, t_test]
	    ph_idx = find_nearest(X[0][:], coord[0])[0]
	    th_idx = find_nearest(Y[:,0], coord[1])[0]
	    ph_val = find_nearest(X[0][:], coord[0])[1]
	    th_val = find_nearest(Y[:,0], coord[1])[1]

	    
	    ### DO NOT TOUCH THIS LINE SPECIFICALLY I DON'T CARE IT LOOKS WRONG IT IS CORRECT###
	    weight_test= P_ANI[th_idx][ph_idx]
	    #ssh no tears only dreams now

	    q = np.random.uniform(0,1)
	    if weight_test>=q:
	        theta_out.append(t_test)
	        phi_out.append(p_test)
	    else:
	        pass

	angles_out = np.column_stack([phi_out, theta_out])

	if savetofile==True:
        if fname == 'none':
            raise ValueError('invalid file name')
        else:
            np.savetxt(fname, angles_out)
    else:
        print('generated {} sky direction samples'.format(n_samples))
        
    return angles_out

#example of usage:
if __name__ == '__main__':
	#look at the underlying distribution
	fig, ax = plt.subplots(subplot_kw=dict(projection='mollweide'), figsize=(10,8))
	im = ax.pcolormesh(X, Y, P_ANI, vmin = 0, vmax = 0.4)
	ax.grid()
	fig.colorbar(im, orientation='horizontal');
	plt.savefig('dipole.png')

    #generate some sky directions to use
    angles= gen_custom_sample_sphere(100, 0.5, savetofile = True, fname='anglesout_alpha05.txt')

    #plot the sky positions you have generated 
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.scatter(angles[:,0], angles[:,1], marker='+', s=50)
    ax.grid(True)
    plt.savefig('skydirections.png')

