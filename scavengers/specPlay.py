#!/usr/bin/env python3.6
'''

Set of tools for modifying spectra


'''

import os, sys
import yaml


from astropy.io import fits, ascii
from astropy import units as u 
from astropy.table import Table

import numpy as np

from scavengers import absInt, cubePlay, cvPlay, headPlay, hiPlay

aBs=absInt.absint()
cP = cubePlay.cubeplay()
cvP = cvPlay.convert()
hP = headPlay.headplay()
hiP = hiPlay.hiplay()



class specplay:
    '''Modules to modify spectra
    
    - rebinSpec
        rebins the spectrum given a binsize     
    '''

    def rebinSpec(self,spec,binSize,fluxCont=None,method='sum'):
        
        '''Rebins spectrum
        Parameters
        ----------

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
        else:
            spec_vec=spec
            spec=np.empty([3,len(spec_vec[0,:])])

            spec[0,:]=spec_vec[0,::-1]/1e3
            spec[1,:]=spec_vec[1,::-1]
            spec[2,:]=spec_vec[1,::-1]*np.nan
            specName=None



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

    def writeSpec(self,cfg_par,source,cubeName,channels,flux,n_pix=None,nameOpt='SoFiA'):

        '''Writes extracted spectrum to table

        Parameters
        ----------
        source : astropy table (single_row)
            properties of the source extracted by SoFiA

        chan: np.array
            array with channel numbers of spectrum

       cubeName: str
            path-to-datacube. Needed to generate array of frequencies and velocities

        freq: np.array
            array with velocities of each channel

        flux: np.array
            array with fluxes of each channel

        nPix: np.array (optional)
            _default=None_, array with fluxes of each channel

        Returns
        -------

            outSpec: str
                full path to output spectrum


        '''

        if 'freq' in source.colnames:
            velocities = cvP.chan2freq(channels, cubeName) 
            if n_pix is not None:
                t = Table([channels,velocities, flux, n_pix], 
                    names=('chan','freq','f_sum', 'n_pix'),
                    meta={'name': 'Spectrum'})
            else:
                t = Table([channels,velocities, flux], 
                    names=('chan','freq','f_sum'),
                    meta={'name': 'Spectrum'})                
            dV = np.abs(np.nanmean(np.diff(cvP.chan2vel(channels,cubeName))))

        elif 'FELO' in source.colnames:
            velocities = cvP.felo2vel(channels, cubeName) 
            if n_pix is not None:

                t = Table([channels,velocities, flux, n_pix], 
                    names=('chan','vel','f_sum', 'n_pix'),
                    meta={'name': 'Spectrum'})
            else:
                t = Table([channels,velocities, flux], 
                    names=('chan','freq','f_sum'),
                    meta={'name': 'Spectrum'})                
            dV = np.abs(np.nanmean(np.diff(velocities)))

        else:
            velocities = cvP.chan2vel(channels,cubeName)
            if n_pix is not None:
                t = Table([channels,velocities, flux, n_pix], 
                    names=('chan','vel','f_sum', 'n_pix'),
                    meta={'name': 'Spectrum'})

            else:
                t = Table([channels,velocities, flux], 
                    names=('chan','freq','f_sum'),
                    meta={'name': 'Spectrum'})                

            dV = np.abs(np.nanmean(np.diff(velocities)))


        outSpec = cfg_par['HIabs']['absDir']+cfg_par['HIabs']['specName']+nameOpt+'.txt'

        t.write(outSpec,overwrite=True,format='ascii')

        return outSpec, dV/1e3



    def intSpecExt(self,cfg_par,source,cubeName,outSpec=None,zRange=None,hiMass=False):
        
        '''Extracts integrated spectrum from source in SoFiA catalogue

        Parameters
        ----------
        source : astropy table (single_row)
            properties of the source extracted by SoFiA

        cubeName: str
            path-to-datacube. The mask produced by Sofia (datacube_mask.fits must be located in the same directory)

        outSpec : str (optional)
            _default=None_, path-to-output spectrum. Default will save the spectrum in the same folder of the datacube (with postfix _specInt.txt)

        zRange: list (optional)
            _default=None_, [chan_min,chan_max] range where to extract the spectrum (default is all channels)
        
        hiMass : bool (optional)
            _default=False_, if set to true HI mass from integrated spectrum and HI mass in science format is returned.

        Returns
        -------

            intSpec: astropy Table
                integrated spectrum : chanNumb, freq (or velocity), integrated flux, number of pixel per channel
        '''
        
        baseCubeName=str.split(cubeName,'.fits')[0]
        cube,xCtr,yCtr = cP.makeSubCube(source,cubeName)
        mask,mxCtr,myCtr = cP.makeSubCube(source, baseCubeName + '_mask.fits')
        pxBeam = hP.pxBeam(cubeName)
        channels = np.asarray(range(cube.shape[0]))
        intFlux = np.zeros(cube.shape[0])
        n_pix = np.zeros(cube.shape[0])

        
        for i in range(0,mask.shape[0]):
            chanFov=cube[i,:,:]
            chanMask=mask[i,:,:]
            idxMask = chanMask>0
            intFlux[i] = np.nansum(chanFov[idxMask])/pxBeam
            n_pix[i] = np.nansum(idxMask)

        outSpecSof,dV = self.writeSpec(cfg_par,source,cubeName,channels,intFlux,n_pix)
        totFlux=np.nansum(intFlux)*dV*u.Jy
        if hiMass==True:

            hiMassVl, hiMassScience = hiP.hiMass(totFlux,cfg_par['galaxy']['dL'],cfg_par['galaxy']['z'])
        

        areaMaskMed =np.median(n_pix[np.nonzero(n_pix)])
        hlfSideMaskMed =int(np.rint(np.sqrt(areaMaskMed)/2.))
        xCtr = int(xCtr)
        yCtr = int(yCtr)
        sqrMaskIdx_X = [xCtr-hlfSideMaskMed,xCtr+hlfSideMaskMed] 
        sqrMaskIdx_Y = [yCtr-hlfSideMaskMed,yCtr+hlfSideMaskMed] 


        for i in range(0,cube.shape[0]):

            if intFlux[i] == 0.:
                chanFov=cube[i,:,:]
                intFlux[i] = np.nansum(chanFov[sqrMaskIdx_Y[0]:sqrMaskIdx_Y[1],sqrMaskIdx_X[0]:sqrMaskIdx_X[1]])/pxBeam
                n_pix[i] = int(areaMaskMed)

        outSpecFull = self.writeSpec(cfg_par,source,cubeName,channels,intFlux,n_pix,nameOpt='Full')

        if hiMass == True:
            return totFlux, hiMassVl, hiMassScience
        else:
            return totFlux













