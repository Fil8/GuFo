#!/usr/bin/env python

class convert:

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