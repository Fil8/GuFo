#!/usr/bin/env python3.6

import os,sys
import numpy as np
from astropy.io import ascii
from scipy import interpolate
from scipy.signal import find_peaks

from astropy.coordinates import Angle
from astropy.io import fits
from astropy.nddata import StdDevUncertainty

from astropy import units as u 
from astropy.io.votable import parse_single_table

from matplotlib import pyplot as plt

from specutils import Spectrum1D
from specutils.analysis import fwzi
from scavengers import headPlay, absInt, specPlay

hP = headPlay.headplay()
aBs = absInt.absint()
sP = specPlay.specplay()

class hiplay(object):

    def __init__(self):

        self.C = 2.99792458e8
        self.Ccm = 2.99792458e10
        self.mH = 1.6735575e-24
        self.mSun = 1.989e33
        self.pc = 3.08577758149137e18
        self.nuHI = 1.42040575177e+09
        self.asecSr = 4.25e10
        self.kb = 1.380649e-16 #erg K-1
        self.Jy = 1e-23
        self.HImassEmission=2.35E5


    def nhi(self,value,bMaj,bMin,dV,z=0.,vunit='m/s'):
        '''

        Module to convert flux Jy beam-1 to column density [cm^-2]

        Parameters
        ----------
        
        value: float
            full path to input mom0 map. 
            Units are in Jy beam-1.

        bMaj: astropy.Coords.Angle()
            size of beam major axis

        bMin: astropy.Coords.Angle()
            size of beam minor axis

        dV: float
            channel width in km/s

        Returns
        -------
            
            nhi: float
                column density value
        
        Notes
        -----

        Conversion formula:

            N_{HI} = 1.10e24*np.power((1+z),2) {SdV}/{Theta^2}
            
            - SdV : integrated flux
            - $Theta^2$ : beam size in **arcseconds**
            
        '''        

        if vunit == 'm/s':
            conversionFactor = 1.10e21*np.power((1+z),2)/(bMaj.arcsecond*bMin.arcsecond)
        elif vunit == 'km/s':
            conversionFactor = 1.10e24*np.power((1+z),2)/(bMaj.arcsecond*bMin.arcsecond)
        nhi = np.multiply(value*dV,conversionFactor)
        nhiScience="{:.2e}".format(nhi)

        return nhi, nhiScience


    def surfBrightHI(self,value):
        '''

        Module to convert column density into surface brightness

        Parameters
        ----------
        
        value: float

            column density value


        Returns
        -------
            
            sB: float
                surface brightness value
        
        Notes
        -----

        Conversion formula:

            N_{HI} = 3.1x10^{17} {SdV}/{Theta^2}
            
            - SdV : integrated flux
            - $Theta^2$ : beam size in **arcminutes**
            
        '''        

        conversionFactor = self.mH*np.power(self.pc,2)/(self.mSun)


        sB = np.multiply(value,conversionFactor)

        return sB

    def HISurfBright(self,value):
        '''

        Module to convert column density into surface brightness

        Parameters
        ----------
        
        value: float

            HI surface brightness value Msun / pc^2


        Returns
        -------
            
            nHI: float
                N(HI)
        
        Notes
        -----

        Conversion formula:

            N_{HI} = 3.1x10^{17} {SdV}/{Theta^2}
            
            - SdV : integrated flux
            - $Theta^2$ : beam size in **arcminutes**
            
        '''        

        conversionFactor = self.mH*np.power(self.pc,2)/(self.mSun)


        nHI = np.divide(value,conversionFactor)

        return nHI  

    def hiMass(self,value,bMaj,bMin,dV,pxSize,DL,z=0.):
        '''

        Estimate HI mass from a given flux density, velocity width, beam size and DL

        Parameters
        ----------
        
        value: float
            flux density value 
            Units are in Jy beam-1 channel

        bMaj: astropy.Coords.Angle()
            size of beam major axis

        bMin: astropy.Coords.Angle()
            size of beam minor axis

        dV: float
            channel width in km/s

        Returns
        -------
            
            nhi: float
                column density value
        
        Notes
        -----

        Conversion formula:

            N_{HI} = 1.10e24*np.power((1+z),2) {SdV}/{Theta^2}
            
            - SdV : integrated flux
            - $Theta^2$ : beam size in **arcseconds**
            
        '''
        
        beamcorr=2.*np.pi*(bMaj.deg*bMin.deg)/(2.35482**2)/(np.power(pxSize.deg,2))

        factor=value*dV/beamcorr
        
        mhi=self.HImassEmission/np.power(1+z,2)*(DL**2)*factor
        mhiScience="{:.2e}".format(mhi)

        return mhi,mhiScience



    def hiMassFromSofiaTable(self,cfg_par):
        '''

        Estimate total HI mass from moment 0 map

        Parameters
        ----------
            
        inMom0: str
            full path to moment zero

        cutoff: float
            threshold within which compute the mass in _Jy*km/s_
 
        z: float
        redshift the source

        DL: float
            luminosity distance of the source in Mpc

        Returns
        -------
            
            mhi: float
                HI mass in Msun

        Notes
        -----

            Conversion formula:
            M(HI) = 2.35x10^5/(1+z)^2*DL^2[Mpc]*Sdv[Jy*km/s]

        '''

        tab = parse_single_table(cfg_par['compute']['inTable'])
        dataTab = tab.array

        if cfg_par['compute']['hiMassFromTable']['vunit'] == 'm/s':
            flux = dataTab['f_sum']/1e3
        else:
            flux = dataTab['f_sum']
        #totFlux = np.sum(flux)

        # bx = Angle(cfg_par['compute']['hiMassFromTable']['beamX'], u.arcsec)
        # by = Angle(cfg_par['compute']['hiMassFromTable']['beamY'], u.arcsec)

        DL=cfg_par['compute']['hiMassFromTable']['dL']
        z=cfg_par['compute']['hiMassFromTable']['z']
        # pixSize= Angle(cfg_par['compute']['hiMassFromTable']['pixSize'], u.arcsec)

        # beamcorr=2.*np.pi*(bx.deg*by.deg/(2.35482**2))/(np.power(pixSize.deg,2))

        mhi=self.HImassEmission/np.power(1+z,2)*(DL**2)*flux    
        
        for i in range(0,len(mhi)):
                
                print('M(HI) = {} x10^9 Msun,'.format(np.round(mhi[i],3)/1e9))

        return mhi

    def hiMassFromMom0(self,inMom0,cutoff,z,DL,zunit):
        '''

        Estimate total HI mass from moment 0 map

        Parameters
        ----------
            
        inMom0: str
            full path to moment zero

        cutoff: float
            threshold within which compute the mass in _Jy*km/s_
 
        z: float
        redshift the source

        DL: float
            luminosity distance of the source in Mpc

        Returns
        -------
            
            mhi: float
                HI mass in Msun

        Notes
        -----

            Conversion formula:
            M(HI) = 2.35x10^5/(1+z)^2*DL^2[Mpc]*Sdv[Jy*km/s]

        '''

        totFlux = self.totalFluxHI(inMom0,cutoff,z,DL,zunit)
        mhi=self.HImassEmission/np.power(1+z,2)*(DL**2)*totFlux

        #print('Sdv [Jy*km/s] = {},'.format(np.round(totFlux,3)))

        #print('M(HI) = {} x10^9 Msun,'.format(np.round(mhi,3)/1e9))
        
        return mhi

    def totalFluxHI(self,inMom0,cutoff,z,DL,zunit='km/s',verbose=True):
        '''

        Estimate total flux of a source given a cutoff in flux density

        Parameters
        ----------
            
        inMom0: str
            full path to moment zero [Jy beam-1*km/s]

        cutoff: float
            threshold within which compute the mass in _Jy beam-1*km/s_


        Returns
        -------
            
            mhi: float
                HI mass in Msun

        Notes
        -----

            Conversion formula:
            S(HI) = Sum(flux_within_cutoff) * 1/beamArea * pixelArea

        '''
        data = fits.getdata(inMom0,0)
        head = fits.getheader(inMom0, 0)
        if zunit=='m/s':
            data/=1e3

        idx = np.where(data<cutoff)
        data[idx] = np.nan
        totFlux = np.nansum(data)
        #totFlux=164452.9274609777e-3
        bx = Angle(head['BMAJ'],u.deg)
        by = Angle(head['BMIN'],u.deg)
        pixSize = Angle(head['CDELT2'], u.deg)

        beamcorr=2.*np.pi*(bx.deg*by.deg)/(2.35482**2)/(np.power(pixSize.deg,2))
        factor=totFlux/beamcorr
        #print(factor,beamcorr)
        #factor=totFlux
        if verbose==True:
            print('Sdv [Jy*km/s] = {},'.format(np.round(factor,3)))
        return factor


    def surfBrightCO(self,value,bmaj,bmin,nu_obs,dV,facMass=4.3):
        '''

        Module converting CO flux density [Jy*km/s] to surface brightness units.
        
        Parameters
        ----------
        
        value: str
            full path to input flux density map.

        bmaj: float
            size of beam major axis in arcseconds

        bmin: float
            size of beam minor axis in arcseconds

        nu_obs: float
            observed frequency of the CO line
 
        Returns
        -------
            
            sB: float
                surface brightness value


        Notes
        -----

            Conversion formula:
            The **brightness temperature** corresponding to an observed flux density S_nu, 
            measured with a telescope of main beam solid angle Omega_b, is the **blackbody 
            temperature an extended object would need to have to produce the observed flux 
            in the Rayleigh-Jeans limit** (h*nu << k_B*T). From the Rayleigh-Jeans law, the
            luminosity emitted by a blackbody per unit area into a unit solid angle, B_nu is 
            given by:


            B_nu = I_nu = S_nu/Omega_beam = 2k_B*nu_line^2*T/c^2
            Rightarrow
            T_B = c^2 S_nu/2 k_B*nu^2*Omega_beam[sr]


            taking into account the constants, that 1 sr = 1 rad^2 = 4.25 x 10^10, and that the flux 
            is typically given in Jy and not in erg.


            Follows that for the CO (1-0) line the surface brightness is 
            
            Sigma_H_2 = X_2 T_{B,CO} [K*km/s]


            where X_2 is the CO J = 1 → 0-to-H_2 conversion factor normalized to 2 x 10^{20} cm^{−2} (K*km/s)
            X_2=4.3 (default) in the Galaxy (Bolatto et al. 2012)


        '''



        fac = 2.*np.log(2.)/(np.pi)*(self.Ccm*self.Ccm)/(self.kb)*self.asecSr*self.Jy


        beamAreaFreq = bmaj*bmin*np.power(nu_obs,2)

        fac = fac*facMass/beamAreaFreq

        sB=np.multiply(value*dV,fac)

        return sB

    def nhiMap(self,inMap,outMap=None,z=0,vunit='m/s',corrFlux=None):
        '''

        Module converting moment 0 map in column density units [cm^-2]

        Parameters
        ----------
        
        inMap: str
            full path to input mom0 map. 
            Units are in Jy beam-1*m/s.

        outMap: str, optional
            full path to output map.

        Returns
        -------
            
            outMap: str
                full path to output mom0 map.
        
        Notes
        -----

        Conversion formula:

            N_{HI} = 3.1x10^{17} {SdV}/{Theta^2}
            
            - SdV : integrated flux in **Jy beam$^{-1}$ m s$^{-1}$**
            - $Theta^2$ : beam size in **arcminutes**
            
        '''
        hh,dd = hP.cleanHead(inMap,writeFile=False)

        bMaj = Angle(hh['BMAJ'], u.deg)
        bMin = Angle(hh['BMIN'], u.deg)
        if vunit == 'm/s':
            conversionFactor = 1.10e21*np.power((1+z),2)/(bMaj.arcsecond*bMin.arcsecond)
        elif vunit == 'km/s':
            conversionFactor = 1.10e24*np.power((1+z),2)/(bMaj.arcsecond*bMin.arcsecond)

        #multiply by conversion factor
        hh['BUNIT']='atoms cm-2'
        
        dd = np.multiply(dd,conversionFactor)
        if corrFlux is not None:
            dd/=corrFlux
        if outMap == None:
            outMap=str.split(inMap,'.fits')[0]
            outMap=outMap+'_nhi.fits'


        fits.writeto(outMap,dd,hh,overwrite=True)

        return outMap

    def nhiCube(self,inCube,z=0.,vunit='m/s',outCube=None):
        '''

        Module converting datacube in column density units [cm^-2]

        Parameters
        ----------
        
        inCube: str
            full path to input datacube. 
            Units are in Jy beam-1.

        outCube: str, optional
            full path to output map.

        Returns
        -------
            
            outCube: str
                full path to output mom0 map.
        
        Notes
        -----

        Conversion formula:

            N_{HI} = 1.10e24*np.power((1+z),2) S/{Theta^2}
            
            - SdV : flux in **Jy beam$^{-1}$            
        '''

        base = fits.open(inCube)
        hh = base[0].header
        dd = base[0].data
        
        bMaj = Angle(hh['BMAJ'], u.deg)
        bMin = Angle(hh['BMIN'], u.deg)
        if vunit == 'm/s':
            dv = hh['CDELT3']/1e3
        elif vunit =='km/s':
            dv = hh['CDELT3']

        conversionFactor = 1.10e24*np.power((1+z),2)/(bMaj.arcsecond*bMin.arcsecond)*dv

        #multiply by conversion factor
        hh['BUNIT']='atoms cm-2'
        
        dd = np.multiply(dd,conversionFactor)

        if outCube == None:
            outCube=str.split(inCube,'.fits')[0]
            outCube=outCube+'_nhi.fits'


        fits.writeto(outCube,dd,hh,overwrite=True)

        return outCube


    def surfBrightHIMap(self,inMap,outMap=None):
        '''

        Module converting column density map to surface brightness units [Msun/pc^2].

        Parameters
        ----------
        
        inMap: str
            full path to input column density map. 
            Units are in cm^-2.

        outMap: str, optional
            full path to output map.

        Returns
        -------
            
            outMap: str
                full path to output surface brightness map.
        
        Notes
        -----

        Conversion formula:

            Sigma_HI [M_sun/pc^2] = N_HI[atom/cm^2] x 1.67*10^-24[g]/atom x(3.08*10^18)^2[cm^2]/pc^2xM_sun/1.989*10^33[g]= 
            8.119*10^-21 x N_HI x M_sun/pc^2
            

            since the column density is N_HI = 3.1*10^17x SdV/Theta^2 then:

            
            Sigma_HI [M_sun/pc^2] = 2.483*10^-3 x SdV/Theta^2 [mJy beam^-1 km s^-1/arcmin^2] 
            
            - m_H = 1.6735575x10^-24 g
            - M_sun = 1.989x10^33 g
            - pc = 3.08577758149137 x 10^18 cm

        '''
        

        hh,dd = hP.cleanHead(inMap,writeFile=False)

        conversionFactor = self.mH*np.power(self.pc,2)/(self.mSun)

        dd = np.multiply(dd,conversionFactor)

        hh['BUNIT'] = 'Msun pc-2'


        if outMap == None:

            outMap=str.split(inMap,'.fits')[0]
            outMap=outMap+'_surfBright.fits'
        fits.writeto(outMap,dd,hh,overwrite=True)

        return outMap

    def fromHztokms(self,inMap,kind):
        '''

        Module converting flux density map from Jy*Hz to Jy*km/s.
        
        Parameters
        ----------
        
        inMap: str
            full path to input flux density map.

        delta_freq: float
            channel-widht in m/s.
 
        Returns
        -------
            
            outMap: str
                full path to output flux density map.        
        '''

        hh,dd = hP.cleanHead(inMap,writeFile=False)
        
        kms=self.C/float(hh['RESTFRQ'])/1e3

        if kind == 'mom0':


            dd=np.multiply(dd,kms)
            hh['BUNIT']='Jy/beam*km/s'

            outMap=str.split(inMap,'.fits')[0]
            outMap=outMap+'_Jybeamkms.fits'
        elif kind=='mom1':
            dd-=float(hh['RESTFRQ'])
            dd = np.multiply(dd,-kms)

            outMap=str.split(inMap,'.fits')[0]
            outMap=outMap+'_kms.fits'
        elif kind =='mom2':
            dd = np.multiply(dd,kms)

            outMap=str.split(inMap,'.fits')[0]
            outMap=outMap+'_kms.fits'

        fits.writeto(outMap,dd,hh,overwrite=True)

        return outMap


    def surfBrightCOMap(self,inMap,facMass=4.3):
        '''

        Module converting CO (1-0) flux density map [Jy*km/s] to surface brightness units [Msun/pc^2].
        
        Parameters
        ----------
        
        inMap: str
            full path to input flux density map.

 
        Returns
        -------
            
            outMap: str
                full path to output surface brightness map.


        Notes
        -----

            Conversion formula:
            The **brightness temperature** corresponding to an observed flux density S_nu, 
            measured with a telescope of main beam solid angle Omega_b, is the **blackbody 
            temperature an extended object would need to have to produce the observed flux 
            in the Rayleigh-Jeans limit** (h*nu << k_B*T). From the Rayleigh-Jeans law, the
            luminosity emitted by a blackbody per unit area into a unit solid angle, B_nu is 
            given by:


            B_nu = I_nu = S_nu/Omega_beam = 2k_B*nu_line^2*T/c^2
            Rightarrow
            T_B = c^2 S_nu/2 k_B*nu^2*Omega_beam[sr]


            taking into account the constants, that 1 sr = 1 rad^2 = 4.25 x 10^10, and that the flux 
            is typically given in Jy and not in erg.


            Follows that for the CO (1-0) line the surface brightness is 
            
            Sigma_H_2 = X_2 T_{B,CO} [K*km/s]


            where X_2 is the CO J = 1 → 0-to-H_2 conversion factor normalized to 2 x 10^{20} cm^{−2} (K*km/s)
            X_2=4.3 (default) in the Galaxy (Bolatto et al. 2012)


        '''
        

        hh,dd = hP.cleanHead(inMap,writeFile=False)

        fac = 2.*np.log(2.)/(np.pi)*(self.Ccm*self.Ccm)/(self.kb)*self.asecSr*self.Jy

        bmaj = hh['BMAJ']
        bmin = hh['BMIN']

        beamAreaFreq = bmaj*3600.*bmin*3600.*np.power(hh['RESTFRQ'],2)

        fac = fac*facMass/beamAreaFreq


        outMap=str.split(inMap,'.fits')[0]
        outMap=outMap+'_surfBright.fits'

        dd = np.multiply(dd,fac)

        hh['BUNIT'] = 'Msun pc-2'

        fits.writeto(outMap,dd,hh,overwrite=True)

        return outMap


    def vRadCube(self,inCube,line='HI'):

        '''
        Module converting cube in frequency to radial velocity units.

        Parameters
        ----------
        
        inCube: str
            full path to input datacube. 
            Units are in Hz

        Returns
        -------
            
            outCube: str
                full path to output surface brightness map.
        
        Notes
        -----

        Conversion formula:

            V(radio) = (nu_0-nu/nu_0) * c = (lambda-lambda_0/lambda) * c 

        To complete for HI and CO

        '''
        f=fits.open(inCube)
        hh=f[0].header
        dd=f[0].data

        if 'RESTFREQ' in hh:
            vRadBlue = (hh['RESTFREQ'] - hh['CRVAL3']) /hh['RESTFREQ'] *self.C
            vRadRed = (hh['RESTFREQ'] - (hh['CRVAL3']+hh['CDELT3']*hh['NAXIS3'])) /hh['RESTFREQ'] *self.C        
        else:
            vRadBlue = (self.nuHI - hh['CRVAL3']) /self.nuHI *self.C
            vRadRed = (self.nuHI - (hh['CRVAL3']+hh['CDELT3']*hh['NAXIS3'])) /self.nuHI *self.C
                

        if line == 'CO':

            hh['CRVAL3']=vRadBlue
            hh['CDELT3']=-(vRadBlue-vRadRed)/hh['NAXIS3']
        #hh['CUNIT3']='km/s'
            hh['CTYPE3']='VRAD'
        
        if line=='HI':

            hh['CRVAL3']=vRadRed
            hh['CDELT3']=(vRadBlue-vRadRed)/hh['NAXIS3']
        #hh['CUNIT3']='km/s'
            hh['CTYPE3']='VRAD'
            dd=np.flip(dd,0)
        
        outCube=str.split(inCube,'.fits')[0]
        outCube=outCube+'_vrad.fits'
        
        fits.writeto(outCube,dd,hh,overwrite=True)

        return outCube



    def computeStats(self,cfg_par,contFlux,specName,dV,bMaj,bMin):

        spec_vec = ascii.read(specName)

        vel = np.array(spec_vec[spec_vec.colnames[0]], dtype=float)
        flux = np.array(spec_vec[spec_vec.colnames[1]], dtype=float)
        tau = np.array(spec_vec[spec_vec.colnames[3]], dtype=float)
        noise = np.array(spec_vec[spec_vec.colnames[2]], dtype=float)
        tauNoise = np.array(spec_vec[spec_vec.colnames[4]], dtype=float)
        if cfg_par['HIabs']['vunit']=='m/s':
            vel/=1e3
            dV/=1e3
        
        peakFlux, idxPeak = aBs.findPeak(cfg_par,vel,flux)

        width,zeroLeft,zeroRight= aBs.fwzi(cfg_par,vel,flux)
       

        fluxArr=np.nansum(flux[idxPeak-zeroLeft:idxPeak+zeroRight])       
        tauArr=np.nansum(tau[idxPeak-zeroLeft:idxPeak+zeroRight])       

        niHIarr = aBs.nhiAbs(tau,np.nanmean(np.diff(vel)))

        
        #niHIRange = aBs.nhiAbs(np.nansum(tau[idxPeak-zeroLeft:idxPeak+zeroRight]),width)

        intNHI =-np.nansum(niHIarr[idxPeak-zeroLeft:idxPeak+zeroRight])
        
        absMass =aBs.mhi_abs(intNHI,bMaj,bMin,cfg_par['galaxy']['dL'],cfg_par['galaxy']['z'])

        return fluxArr,tauArr,width,intNHI,absMass


    def hiatz(self,z):
        hobs = self.nuHI/(1+z)/1e06   #MHz
    
        return hobs

