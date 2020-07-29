
from __future__ import print_function, division

import os

import numpy as np
import healpy as hp

def make_map_vec(theta, phi, data, NSIDE=4, normalize = False):
    """
    Turn a vector map of data for each theta and phi into a healpix map

    from https://stackoverflow.com/questions/59492213/turn-a-fits-file-into-healpix-map
    The user uncovers an interesting bug where you get lots of unseen pixels if you crank NSIDE higher than 8, I have not yet
    figured out the origin of this bug, so leave NSIDE to be low for now

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
    normalize: Bool, to normalize the map or not
            
    Returns
    -------
    list
        weights of each healpix pixel
    """
    assert len(theta) == len(phi) == len(data)
    #make a healpix map of the sky
    e1map = np.full(hp.nside2npix(NSIDE), hp.UNSEEN, dtype=np.float)
    #convert the theta and phi coordinates into healpix pixel indicies
    index = hp.ang2pix(NSIDE, theta, phi)
    # bin the data into the healpix pixels: count up the data for all thetas and phis in bin i
    values = np.fromiter((np.sum(data[index==i]) for i in np.unique(index)), float, count=len(np.unique(index)))
    if normalize == True:
    	#if you want to normalize the map, do it here, not usually necessary
    	values = [i/np.max(values) for i in values]
    	e1map[np.unique(index)] = values
    else:
    	e1map[np.unique(index)] = values

    #you can not plot this map using hp.mollview(e1map)

    return e1map
