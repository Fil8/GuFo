#!/usr/bin/env python3.6

import numpy as np
from astropy.io import ascii
from scipy import interpolate

class convert(object):

    def __init__(self):

        self.C = 2.99792458e8


    def lambdaVRad(self,lambdaWave,lambdaRest):

        # wavelenght in Angstrom
        vRad = ((lambdaWave-lambdaRest)/lambdaRest)*self.C/1e3

        #velocity in km/s
        return vRad

    def vRadLambda(self,vRad,lambdaRest):
        
        # velocity in kms
        lambdaWave = (vRad*1e3*lambdaRest*1e-10)/(self.C)/1e-10 +lambdaRest

        # wavelenght in Angstrom
        return lambdaWave

    def specRes(self,cfg_par):


        MUSEfile = cfg_par['general']['gfitPath'] + 'MUSEspecs.csv'

        MUSETable = ascii.read(MUSEfile)

        x = MUSETable['Lambda']
        y = MUSETable['R']

        tck = interpolate.splrep(x, y, s=0)
        xnew = np.linspace(np.min(x),np.max(x),1e4)
        ynew = interpolate.splev(xnew, tck, der=0)
        dVlambda = 1./ynew *xnew
        return dVlambda


    def hms2deg(self,ra_hms):
        
        ra = str.split(ra_hms, ':')

        hh = float(ra[0])*15
        mm = (float(ra[1])/60)*15
        ss = (float(ra[2])/3600)*15
        
        raDeg = hh+mm+ss

        return raDeg


    # DMS -> degrees
    def dms2deg(self,dec_dms):
        dec = str.split(dec_dms, ':')
        
        dd = abs(float(dec[0]))
        mm = float(dec[1])/60
        ss = float(dec[2])/3600
        
        if float(dec[0])>= 0:
            decDeg = dd+mm+ss
        else:
            decDeg = -(dd+mm+ss)
        
        return decDeg