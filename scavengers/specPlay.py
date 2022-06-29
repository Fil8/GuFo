#!/usr/bin/env python3.6
'''

Set of tools for modifying spectra


'''

import os, sys
import yaml


from astropy.io import fits, ascii

import numpy as np

from scavengers import absInt
from astropy.table import Table

aBs=absInt.absint()

class specplay:
    '''Modules to modify spectra
    
    - rebinSpec
        rebins the spectrum given a binsize     
    '''

    def rebinSpec(self,spec,binSize,fluxCont=None,method='sum'):
        
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
        

        if isinstance(spec,str):  
            specName=spec          
            spec_vec = ascii.read(spec)
            vel = np.array(spec_vec[spec_vec.colnames[0]], dtype=float)
            flux = np.array(spec_vec[spec_vec.colnames[1]], dtype=float)
            noise = np.array(spec_vec[spec_vec.colnames[2]], dtype=float)
                        
            spec=np.empty([3,len(vel)])
            spec[0,:] = vel[::-1]/1e3
            spec[1,:] = flux[::-1]
            spec[2,:] = noise[::-1]
        



        minSpec = spec[0,0]
        maxSpec = spec[0,-1]
        binSpec=np.arange(minSpec,maxSpec,binSize)
        vels = np.hstack([1.5*spec[0,0]-0.5*spec[0,1],
                      0.5*(spec[0,0:-1]+spec[0,1:]),
                      1.5*spec[0,-1]-0.5*spec[0,-2]])

        newSpec = np.zeros((3,len(binSpec)))
        newSpec[0,:] = binSpec
        for i in range(0, newSpec.shape[1]-1):
            #look for the right velocity bin
            index = np.where(np.logical_and(binSpec[i] <= spec[0,:], spec[0,:] < binSpec[i+1]))

            #update the flux bin
            if method == 'sum':
                newSpec[1,i] = np.nansum(spec[1,index])
                newSpec[2,i] = np.nansum(spec[2,index])

            elif method == 'mean':
                newSpec[1,i] = np.nanmean(spec[1,index])
                newSpec[2,i] = np.nanmean(spec[2,index])
            elif method == 'median':
                newSpec[1,i] = np.nanmedian(spec[1,index])
                newSpec[2,i] = np.nanmedian(spec[2,index])

        if specName is not None:    

            # measure noise in the spectrum outside of the line
            end_spec = float(len(vels))
            end_spec_th = int(end_spec/3.)
            end_spec = int(end_spec)
            mean_rms = (np.std(newSpec[1,0:end_spec_th]) +
                                np.std(newSpec[1,end_spec-end_spec_th:end_spec])) / 2.
            mean_rms_arr = np.zeros(newSpec.shape[1])+mean_rms
            if fluxCont is not None:
                tau = aBs.optical_depth(newSpec[1,::-1], fluxCont)
                tau_noise=aBs.optical_depth(newSpec[2,::-1],fluxCont)
                nhi=aBs.nhiAbs(-tau,binSize)
                nhi_noise=aBs.nhiAbs(-tau_noise,binSize)
            else:
                tau = np.zeros(len(vels))
                tau_noise =  np.zeros(len(vels))
                nhi = np.zeros(len(vels))
                nhi_noise =  np.zeros(len(vels))
            xcol = 'Velocity [km/s]'

            t = Table([newSpec[0,::-1]*1e3, newSpec[1,::-1], newSpec[2,::-1], -tau, -tau_noise, nhi, nhi_noise], names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'NHI', 'noise_NHI'),meta={'name': 'Spectrum'})
            #t = Table([tau, tau_noise], names=('Optical depth','Noise optical depth'),meta={'name': 'Spectrum'})
           

            out_spec=specName.replace('.txt','_rb.txt')
            ascii.write(t,out_spec,overwrite=True)
        
            return out_spec
        else:
            return newSpec
















