### generate random samples of the celestial sphere. (C) F H Panther 2020

from __future__ import print_function, division

import os

import numpy as np
import healpy as hp

from astropy import units
from astropy.coordinates import SkyCoord
from scipy.special import sph_harm


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

def make_map_vec(theta, phi, data, NSIDE=4):
    """
    Turn a vector map into a healpix map

    Parameters
    ----------
    theta : float, list
            values of theta over which data is defined
    phi: float, list
            values of phi over which data is defined
    data: float, array
            weights of the data at each theta, phi 
    NSIDE: int, must be a power of 2 no larger than 8
            for healpix        
            
    Returns
    -------
    list
        weights of each healpix pixel
    """
    assert len(theta) == len(phi) == len(data)
    e1map = np.full(hp.nside2npix(NSIDE), hp.UNSEEN, dtype=np.float)
    index = hp.ang2pix(NSIDE, theta, phi)
    values = np.fromiter((np.sum(data[index==i]) for i in np.unique(index)), float, count=len(np.unique(index)))
    values = [i/np.max(values) for i in values]
    e1map[np.unique(index)] = values
    return e1map

def gen_uniform_sample_sph_harm(n_samples, alpha, l, m, NSIDE = 4, savetofile = False, fname ='none'):
    """
    generate a user-defined random number of points on the celestial sphere distributed 
    according to spherical harmonics parameterized by l, m and alpha

    Parameters
    ----------
    n_samples : int
            the number of samples to generate
    alpha: float, array
            the weight to apply to spherical harmonic i
    l: int, array
            the l value for spherical harmonic i
    m: int, array
            the m value for spherical harmonic i, limited by l   
    NSIDE: int, must be a power of 2 and no greater than 8
            for healpix     
    savetofile: Bool
            save to file in current directory, true or false
    fname: string
            filename to save as, 'none' if savetofile = False
            
    Returns
    -------
    numpy ndarray
        table of RA and DEC in radians sampled from the desired linear combination of spherical harmonics
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

    SH_SP = np.zeros((100,100))

    assert len(alpha) == len(l) == len(m)

    for i in range(len(alpha)):
        sp_temp = alpha[i]*sph_harm(m[i],l[i], PHI, THETA).real
        SH_SP = SH_SP+sp_temp

    #turn the spherical harmonic map into a healpix map
    hpmap = make_map_vec(THETA, PHI, abs(SH_SP))

    theta_out = []
    phi_out = []

    while len(theta_out)<n_samples:

    #generate an ra and dec - should be an underlying uniform distribution
        u = np.random.uniform(0,1)
        v = np.random.uniform(0,1)
        ph = 2*np.pi*u
        th = np.arccos(2*v-1)

        #calculate the hp index of this theta and phi, and get it's value
        idx = hp.ang2pix(NSIDE, th, ph)
        pix_val = hpmap[idx]

        q = np.random.uniform(min(hpmap), max(hpmap))
        print(q)

        if pix_val>q:
            theta_out.append(th-np.pi/2)
            phi_out.append(ph)
        else:
            pass

    c = SkyCoord(ra=phi_out*units.radian, dec=theta_out*units.radian, frame='icrs')
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
if __name__ == '__main__':
    #generate some sky directions to use
    angles= gen_uniform_sphere_sample(100, savetofile = True, fname='anglesout.txt')

    #plot the sky positions you have generated 
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.scatter(angles[:,0], angles[:,1], marker='+', s=50)
    ax.grid(True)
    plt.savefig('skydirections.png')


    angles_dipole = gen_uniform_sample_sph_harm(100, [1], [1], [0], NSIDE = 4, savetofile = True, fname ='dipole_out.txt')

    #plot the sky directions from the dipole
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.scatter(angles_dipole[:,0], angles_dipole[:,1], marker='+', s=50)
    ax.grid(True)
    plt.savefig('skydirectionsdipole.png')

