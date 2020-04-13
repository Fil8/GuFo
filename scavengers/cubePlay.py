#!/usr/bin/env python3.6

import os, sys
import yaml

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult
from lmfit.model import load_modelresult


from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np
#import numpy.ma as ma

import tPlay,cvPlay,bptPlot,momPlot

tP = tPlay.tplay()
cvP = cvPlay.convert()
mPl = momPlot.MOMplot()

class cubeplay:
    '''Modules to create cubelets of real fitted lines and residuals
    - makeHeader
        make header of line in vrad velocity frame 
    - makeLineCube
        make cubelets for each line marked in lineList.txt
    '''
    def makeHeader(self,cfg_par,lineWave,header,waveAng):
        '''
        Defines header of output line cubelets
        Puts wcs coordinates in datacube from information provided in the configuration file
        Sets velocity axis to vrad

        Parameters:
            - cfg_par: configuration file
            - header of datacube
        
        Returns:
            - header of datacubes
        '''       
            
        vel=cvP.lambdaVRad(waveAng,lineWave)+float(cfg_par['general']['velsys'])
        cdelt3=np.mean(np.ediff1d(vel))

        if 'CRDER3' in header:
            del header['CRDER3']

        header['CRPIX3'] = len(vel)
        header['CRVAL3'] = vel[0]*1e3
        header['CDELT3'] = -cdelt3*1e3
        header['CTYPE3'] = "VRAD"
        header['SPECSYS'] = "barycent"
        header['CELLSCAL'] = 'constant'
        return header
    

    def makeLineCubes(self,cfg_par):

        cubeletsDir = cfg_par['general']['cubeletsDir']
        cubeDir = cfg_par['general']['cubeDir']
        modName = cfg_par['gFit']['modName']
        momModDir = cfg_par['gFit']['modName']

        if cfg_par['cubelets']['cube'] == 'vorLine': 
            f = fits.open(cfg_par['general']['dataCubeName'])
            dd = f[0].data
            hh = f[0].header

        elif cfg_par['cubelets']['cube'] == 'fitLine': 
            f = fits.open(cubeDir+'fitCube_'+modName+'.fits')
            mom = fits.open(momModDir+'mom0_OIII5006.fits')
            mm=mom[0].data
            indexFltr = np.where(mm==np.nan)
            dd=dd[:,[indexFltr]]

        else:
            f = fits.open(cubeDir+cfg_par['cubelets']['cube'])
            dd = f[0].data
            hh = f[0].header


        lineInfo = tP.openLineList(cfg_par)
        index = np.where(lineInfo['Cubelets'] == 0)
        fltr =  np.array(index)[0]
        lineInfo.remove_rows(list(fltr))
        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec = tP.openVorLineOutput(cfg_par,cfg_par['general']['outVorLineTableName'],
            cfg_par['general']['outVorSpectra'])


        lambdaMin = np.log(cfg_par['gFit']['lambdaMin'])
        lambdaMax = np.log(cfg_par['gFit']['lambdaMax'])
        idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
        idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0])

        wave=wave[idxMin:idxMax]
        dd=dd[idxMin:idxMax,:,:]
 
        velRange = float(cfg_par['cubelets']['velRange'])
       
        for ii in range(0,len(lineInfo['ID'])):

            lineNameStr = str(lineInfo['Name'][ii])

            if '[' in lineNameStr:
                lineName = lineNameStr.replace("[", "")
                lineName = lineName.replace("]", "")
                lineName = lineName+str(int(lineInfo['Wave'][ii]))
            else:
                lineName = lineNameStr+str(int(lineInfo['Wave'][ii]))

            lineNameStr=lineNameStr+str(int(lineInfo['Wave'][ii]))


            print('\n\t         +++\t\t    '+lineName+'\t\t +++')
            
            waveMin =  np.log(lineInfo['Wave'][ii] - lineInfo['lineRangeAng'][ii])
            waveMax =  np.log(lineInfo['Wave'][ii] + lineInfo['lineRangeAng'][ii])
            idxMin = int(np.where(abs(wave-waveMin)==abs(wave-waveMin).min())[0]) 
            idxMax = int(np.where(abs(wave-waveMax)==abs(wave-waveMax).min())[0] )
            
            dCbl = dd[idxMin:idxMax,:,:]

            waveAng=np.exp(wave[idxMin:idxMax])

            header = self.makeHeader(cfg_par,lineInfo['Wave'][ii],hh,waveAng)
            
            if cfg_par['cubelets']['cube'] == 'vorLine': 
                outCubelet = cubeletsDir+str(lineNameStr)+'_measVor.fits'
            elif cfg_par['cubelets']['cube'] == 'fitLine': 
                outCubelet = cubeletsDir+str(lineNameStr)+'_fit.fits'          
            else:
                outCubelet = cubeletsDir+str(lineNameStr)+'_meas.fits'

            fits.writeto(outCubelet,np.flip(dCbl,axis=0),header,overwrite=True)

        return