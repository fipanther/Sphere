### generate random samples of the celestial sphere. (C) F H Panther 2020

from __future__ import print_function, division

import os

import numpy as np

from astropy import units
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt

CWD = os.getcwd

def gen_uniform_sphere_sample(n_samples, savetofile = False, fname = 'none'):
    """
    generate a user-defined random number of points on the celestial sphere

    Parameters
    ----------
    n_samples : int
            the number of samples to generate
    savetofile: Bool
            save to file in current directory, true or false
    fname: string
            filename to save as, 'none' if savetofile = False
            
    Returns
    -------
    numpy ndarray
        table of RA and DEC in radians
    """
    
    # generate two random numbers that are uniformly distributed between 0 and 1
    u = np.random.uniform(0,1, size=n_samples)
    v = np.random.uniform(0,1, size=n_samples)
    
    #generate the angle around the equator
    phi = 2*np.pi*u
    #generate the polar angle - must be weighted at the equator or you will draw
    #more points at the poles
    theta = np.arccos(2*v-1)
    theta = [i-np.pi/2 for i in theta]

    c = SkyCoord(ra=phi*units.radian, dec=theta*units.radian, frame='icrs')
    ra = c.ra
    #the sky is weird so wrap the angles at 180 degrees
    ra = ra.wrap_at(180*units.degree)
    dec = c.dec
    out = np.column_stack([ra.radian, dec.radian])
    
    if savetofile==True:
        if fname == 'none':
            raise ValueError('invalid file name')
        else:
            np.savetxt(fname, out)
    else:
        print('generated {} sky direction samples'.format(n_samples))
        
    return out

#example of usage:
# if __name__ == '__main__':
# 	angles= gen_uniform_sphere_sample(100, savetofile = True, fname='anglesout.txt')

# 	#plot the sky positions you have generated 
# 	fig = plt.figure(figsize=(16, 12))
# 	ax = fig.add_subplot(111, projection="mollweide")
# 	ax.scatter(angles[:,0], angles[:,1], marker='+', s=50)
# 	ax.grid(True)
# 	plt.savefig('skydirections.png')