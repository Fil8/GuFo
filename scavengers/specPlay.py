#!/usr/bin/env python3.6
'''

Set of tools for modifying spectra


'''

import os, sys
import yaml


import numpy as np


class specplay:
    '''Modules to modify spectra
    
    - rebinSpec
        rebins the spectrum given a binsize     
    '''

    def rebinSpec(self,spec,binSize,method='sum'):
        
        '''Rebins spectrum
        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        spec: np.array()
            input array (2,N) np.array

        binSize: float
            size of the new spectrum

        function: str
            how to rebin the data: sum, median, mean

        Returns
        -------

            newSpec: np.array()
                binned spectrum
        '''
        


        minSpec = spec[0,0]
        maxSpec = spec[0,-1]
        print(minSpec,maxSpec)
        binSpec=np.arange(minSpec,maxSpec+binSize,binSize)

        vels = np.hstack([1.5*spec[0,0]-0.5*spec[0,1],
                      0.5*(spec[0,0:-1]+spec[0,1:]),
                      1.5*spec[0,-1]-0.5*spec[0,-2]])

        newSpec = np.zeros([len(vels),len(vels)])
        newSpec[0,:] = vels

        for i in range(0, newSpec.shape[1]-1):
            #look for the right velocity bin
            index = (newSpec[0,i] <= spec[0,:]) & (spec[0,:] < newSpec[0,i+1])

            #update the flux bin
            if method == 'sum':
                newSpec[1,i] = np.nansum(spec[1,index])
            elif method == 'mean':
                newSpec[1,i] = np.nanmean(spec[1,index])
            elif method == 'median':
                newSpec[1,i] = np.nanmedian(spec[1,index])

        return newSpec
















