#!/usr/bin/env python3.6

import os
import numpy as np
from astropy.io import ascii
from scipy import interpolate

from astropy.coordinates import Angle
from astropy.io import fits
from astropy import units as u 
from astropy.io.votable import parse_single_table

from matplotlib import pyplot as plt


import headPlay

hP = headPlay.headplay()



class convert(object):
    '''This class converts units in useful various ways.
    
    '''
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



    def pccm(self,dist,option='cmTopc'):
        '''Converts parsecs to cm and viceversa
        
        Parameters
        ----------

        pc: np.array()
            array with distances in parsecs
        
        option: str
            cmTopc for centimeter to parsec conversion
            pcTocm for parsec to centimeter conversion

        Returns
        -------

        distConv: np.array()
            array with distances in cm

        Notes
        ------

        1 pc = 3.08577758149137e18

        '''        

        if option=='cmTopc':

            distConv = np.divide(dist,self.pc)

        elif option=='pcTocm':

            distConv = np.multiply(dist,self.pc)

        return distConv

    def lambdaVRad(self,lambdaWave,lambdaRest):
        '''Converts wavelength to radial velocity
        
        Parameters
        ----------

        lambdaWave: np.array()
            wavelengths in Angstrom to be converted
        
        lambdaRest: float
            rest wavelenght of reference

        Returns
        -------

        vRad: np.array()
            array with radial velocities in km/s

        ''' 

        vRad = ((lambdaWave-lambdaRest)/lambdaRest)*self.C/1e3

        return vRad

    def vRadLambda(self,vRad,lambdaRest):
        '''Converts radial velocity to wavelength
        
        Parameters
        ----------

        vRad: np.array()
            array with radial velocities in km/s
        
        lambdaRest: float
            rest wavelenght of reference

        Returns
        -------

        lambdaWave: np.array()
            wavelengths in Angstrom

        '''         
        # velocity in kms
        lambdaWave = (vRad*1e3*lambdaRest*1e-10)/(self.C)/1e-10 +lambdaRest

        # wavelenght in Angstrom
        return lambdaWave

    def specRes(self,cfg_par):
        '''Computes the spectral resolution given the MUSE specs
        
        Parameters
        ----------

        cfg_par: configuration file
        Returns
        -------

        dVlambda: np.array()
            spectral resolution in the MUSE wavelength range

        ''' 

        MUSEfile = cfg_par['general']['gfitPath'] + 'MUSEspecs.csv'

        MUSETable = ascii.read(MUSEfile)

        x = MUSETable['Lambda']
        y = MUSETable['R']

        tck = interpolate.splrep(x, y, s=0)
        xnew = np.linspace(np.min(x),np.max(x),1e4)
        ynew = interpolate.splev(xnew, tck, der=0)
        dVlambda = 1./ynew *xnew

        return dVlambda/(2.*np.sqrt(2.*np.log(2.)))
        #return dVlambda


    def hms2deg(self,ra_hms):
        '''
        Converts Right Ascension from HH:MM:SS to degrees

        Parameters
        ----------
        
        ra_hms: str
            right ascension 'hh:mm:ss'

        Returns
        -------
            
            raDeg: float
                right ascension in degrees
            
        '''        
        ra = str.split(ra_hms, ':')

        hh = float(ra[0])*15
        mm = (float(ra[1])/60)*15
        ss = (float(ra[2])/3600)*15
        
        raDeg = hh+mm+ss

        return raDeg


    def dms2deg(self,dec_dms):
        '''
        Converts Declination from DD:MM:SS to degrees

        Parameters
        ----------
        
        dec_hms: str
            right ascension 'dd:mm:ss'

        Returns
        -------
            
            decDeg: float
                right ascension in degrees
            
        '''        
        dec = str.split(dec_dms, ':')
        
        dd = abs(float(dec[0]))
        mm = float(dec[1])/60
        ss = float(dec[2])/3600
        
        if float(dec[0])>= 0:
            decDeg = dd+mm+ss
        else:
            decDeg = -(dd+mm+ss)
        
        return decDeg

    def deg2dms(dec_deg):
        '''
        Converts Declination from degrees to DD:MM:SS

        Parameters
        ----------
        
        decDeg: float
            right ascension in degrees

        Returns
        -------
            
        dec_hms: str
            right ascension 'dd:mm:ss'
            
        '''          
        negative = dec_deg < 0
        dec_deg = abs(dec_deg)
        minutes,seconds = divmod(dec_deg*3600,60)
        degrees,minutes = divmod(minutes,60)
        if negative:
            if degrees > 0:
                degrees = -degrees
            elif minutes > 0:
                minutes = -minutes
            else:
                seconds = -seconds
        return (degrees,minutes,seconds)


    def electronDensity(self,R):
        '''Computes the electron density from the [SII]6716/[SII]6731 ratio
        Parameters
        ----------

        R: np.array()
            array witth the [SII]6716/[SII]6731 ratio

        Returns
        -------

        nE: np.array()
            array with the electron density in log10 values

        Notes
        ------

        Empirical formula for electron density taken from : 10.1051/0004-6361/201322581
        previously published in Osterbork & Farland 2003, Astrophysics of Gaseous Nebulae and Active Galactic Nuclei

        '''

        R2=np.power(R,2)
        R3=np.power(R,3)

        C1=0.0543*np.tan(-3.0553*R+2.8506)+6.98
        C2=-10.6905*R
        C3=+9.9186*R2
        C4=-3.5442*R3

        nE=C1+C2+C3+C4

        indexMin = np.power(10,nE)<40.
        nE[indexMin] = np.log10(40.)

        indexMax = np.power(10,nE)>1e4

        nE[indexMax] = np.log10(1e4)

        return nE

    def dustExtinction(self,R):
        '''Computes the electron density from the [SII]6716/[SII]6731 ratio
        Parameters
        ----------

        R: np.array()
            array witth the [SII]6716/[SII]6731 ratio

        Returns
        -------

        nE: np.array()
            array with the electron density in log10 values

        Notes
        ------

        Empirical formula for electron density taken from : 10.1051/0004-6361/201322581
        previously published in Osterbork & Farland 2003, Astrophysics of Gaseous Nebulae and Active Galactic Nuclei

        '''
        EBV = 1.97*np.log10(R/2.86)
        AHalpha = 3.33*EBV
        Av = 4.05*EBV
        
        return AHalpha, Av

    def electronDensityX(self,Lx,Lambda,V):
        '''Computes the electron density from X ray luminosity and cooling function (from Schure et al. 2009)
        Parameters
        ----------

        Lx: np.array()
            array with the X-ray thermal luminosity in erg s^-1

        Lambda: np.array()
            corresponding cooling function values

        V: np.array()
            corresponding volumes

        Returns
        -------

        nE: electron density in 1/cm^-3 ?

        Notes
        ------

        n_e ~ (frac{Delta L_X}{Delta V}frac{mu_i/mu_e}{Lambda(T_X,Z)})^{1/2}
        
        mu_i = 1.30, mu_e = 1.18
        Delta V = 4/3 pi (R_out^3-R_in^3)
        '''

        muRatio = 1.30/1.18

        constants = np.divide(muRatio,Lambda)
        LfracV = np.divide(Lx,V)

        nE = np.sqrt(np.multiply(LfracV,constants))

        return nE

    def temperature(self,TkeV):
        '''Converts keV to Kelvin dividing by converting to Joules and dividing by the Boltzmann's constant in J/K
        Parameters
        ----------

        Tkev: float
            temperature to convert

        Returns
        -------

        TK: np.array()
            temperature in Kelvin

        Notes
        ------
        1/kB = 1.602e-19 J/eV / 1.380e-23 J/K = 11604.51812 K/eV
            
        '''

        TK = np.multiply(TkeV, 11604.51812)

        return TK

    def tcool(self,ne,T,Lambda):
        '''Computes the cooling time from the electron density, temperature and cooling function of the gas
        Parameters
        ----------

        ne: np.array()
            electron density in cm^-3

        T: np.array()
            temperatures in K

        Lambda: np.array()
            cooling function in erg s-1 cm3

        Returns
        -------

        TK: np.array()
            temperature in Kelvin

        Notes
        ------
        t_cool = 3/2(n_e+n_i)k_B T/(n_e*n_i Lambda(T))
            
        '''        

        ni= np.multiply(0.92,ne)
        
        density=3./2.*(ni+ne)

        densSq = np.multiply(ne,ni)

        kT = np.multiply(self.kb, T)

        
        denominator = np.multiply(densSq,Lambda)

        tcool = np.divide(np.multiply(kT,density),denominator)

        yrinsec = 3600.*24.*365.2425
        
        tcool = np.divide(tcool,yrinsec)

        return tcool

    def nhi(self,value,bMaj,bMin,dV,z=0.):
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
            channel width in m/s

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

#        conversionFactor = 1.10e21*np.power((1+z),2)/(bMaj.arcsecond*bMin.arcsecond)
        conversionFactor = 1.10e21*np.power((1+z),2)/(bMaj*bMin)

        nhi = np.multiply(value*dV,conversionFactor)

        return nhi


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

    def hiMassFromMom0(self,inMom0,cutoff,z,DL):
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


        data = fits.getdata(inMom0,0)
        head = fits.getheader(inMom0, 0)
        idx = np.where(data<cutoff)
        data[idx] = np.nan
        totFlux = np.nansum(data/1e3)
        bx = Angle(head['BMAJ'],u.deg)
        by = Angle(head['BMIN'],u.deg)
        pixSize = Angle(head['CDELT2'], u.deg)
        #dl=c.lum_dist(cfg_par['compute']['hiMassFromTable']['dL'])/3.085678e24
        beamcorr=2.*np.pi*(bx.deg*by.deg)/(2.35482**2)/(np.power(pixSize.deg,2))
        factor=totFlux/beamcorr
        mhi=self.HImassEmission/np.power(1+z,2)*(DL**2)*factor

        print('Sdv [Jy*km/s] = {},'.format(np.round(factor,3)))

        print('M(HI) = {} x10^9 Msun,'.format(np.round(mhi,3)/1e9))
        
        return mhi

    def totalFlux(self,inMom0,cutoff,z,DL):
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
            S(HI) = Sum(flux_within_cutoff) / beamArea * pixelArea

        '''



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

    def nhiMap(self,inMap,outMap=None,z=0,vunit='m/s'):
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
        print(inMap)
        hh,dd = hP.cleanHead(inMap,writeFile=False)

        bMaj = Angle(hh['BMAJ'], u.deg)
        bMin = Angle(hh['BMIN'], u.deg)
        if vunit == 'm/s':
            conversionFactor = 1.10e21*np.power((1+z),2)/(bMaj.arcsecond*bMin.arcsecond)
        elif vunit == 'km/s':
            conversionFactor = 1.10e24*np.power((1+z),2)/(bMaj.arcsecond*bMin.arcsecond)

        print(conversionFactor)
        #multiply by conversion factor
        hh['BUNIT']='atoms cm-2'
        
        dd = np.multiply(dd,conversionFactor)

        if outMap == None:
            outMap=str.split(inMap,'.fits')[0]
            outMap=outMap+'_nhi.fits'


        fits.writeto(outMap,dd,hh,overwrite=True)

        return outMap

    def nhiCube(self,inCube,dv,z=0.,outCube=None):
        '''

        Module converting datacube in column density units [cm^-2]

        Parameters
        ----------
        
        inCube: str
            full path to input datacube. 
            Units are in Jy beam-1.

        outCube: str, optional
            full path to output map.

        dV : float
            channel width in km/s

        Returns
        -------
            
            outCube: str
                full path to output mom0 map.
        
        Notes
        -----

        Conversion formula:

            N_{HI} = 1.10e24*np.power((1+z),2) S/{Theta^2}
            
            - SdV : flux in **Jy beam$^{-1} m/s$            
        '''

        base = fits.open(inCube)
        hh = base[0].header
        dd = base[0].data
        
        bMaj = Angle(hh['BMAJ'], u.deg)
        bMin = Angle(hh['BMIN'], u.deg)

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
        
        print(hh)
        outCube=str.split(inCube,'.fits')[0]
        outCube=outCube+'_vrad.fits'
        
        fits.writeto(outCube,dd,hh,overwrite=True)

        return outCube


    def computeNoise(self,imName):

        niter = 4 # nr of noise measurements
        clip = 4 
        medNoise = []

        unit_dict = {
            'vrad'   :   [1e+3, 'Velocity [km/s]'],
            'freq'   :   [1e+6, 'Frequency [MHz]'],
        }
        
        f=fits.open(imName)
        head=f[0].header
        cube=f[0].data

        cube=cube[:,cube.shape[1]//4:3*cube.shape[1]//4,cube.shape[2]//4:3*cube.shape[2]//4]

        iter=1
        while iter<niter:
          std=np.nanmedian(np.nanstd(cube,axis=(1,2)))
          print('# iter {1:d}, median std = {0:.2e} Jy/beam'.format(std,iter))
          cube[np.abs(cube)>clip*std]=np.nan
          iter+=1

        noise=np.nanstd(cube,axis=(1,2))
        print('# iter {1:d}, median std = {0:.2e} Jy/beam'.format(np.nanmedian(noise),iter))

        freqs=(np.arange(head['naxis3'])-(head['crpix3']-1))*head['cdelt3']+head['crval3']

        BaseName = os.path.basename(imName)

        plt.plot(freqs/unit_dict[head['CTYPE3'].lower()][0],noise*1e+3,'k-')
        plt.axhline(y=np.median(noise)*1e+3,linestyle=':',color='k')
        plt.title(BaseName)
        plt.xlabel(unit_dict[head['CTYPE3'].lower()][1])
        plt.ylabel('Noise [mJy/beam]')

        plt.savefig(imName.replace('.fits','rms.png'))
        plt.clf()





