#!/usr/bin/env python
import os, sys
import yaml

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult

from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np
import numpy.ma as ma

import shutil
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator
from matplotlib import transforms as mtransforms
from matplotlib.ticker import LogFormatter 
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import gufo as gf
import cvPlay
import specPlot
import tPlay

#gf = gufo()
cvP = cvPlay.convert()
sP = specPlot.specplot()
tP = tPlay.tplay()

class gplay:
    

    def modDef(self,modName):

        
        Gmod = GaussianModel()

        gauss1 = GaussianModel(prefix='g1_ln'+str(i)+'_')

        pars = gauss1.make_params()

        pars['g1_ln'+0+'_'+'center'].set(value=np.argmax[y[indexMin:indexMax]],
            min=waveMin,max=waveMax)
        pars['g1_ln'+str(i)+'_'+'height'].set(value=ampIn1)
        pars['g1_ln'+str(i)+'_'+'sigma'].set(expr='g1_sigma')


        gauss2 = GaussianModel(prefix='g2_')
        pars.update(gauss2.make_params())

        pars['g2_center'].set(value=-10, min=wMin, max=wMax)
        pars['g2_sigma'].set(value=100, min=0)
        pars['g2_amplitude'].set(value=2000, min=1)

        cenIn3 = cenIn1 + dltV13
        sigmaIn3 = pars['g1_ln'+str(i)+'_'+'sigma'] + dltSgm12
        ampIn3 = ampIn1 + dltAmp13

        gauss3 = GaussianModel(prefix='g3_')
        pars.update(gauss3.make_params())

        pars['g3_center'].set(value=-20, min=wMin, max=wMax)
        pars['g3_sigma'].set(value=200, min=0)
        pars['g3_amplitude'].set(value=2000, min=1)

        gName1 = ['g1']
        gName2 = ['g1', 'g2']   
        gName3 = ['g1', 'g2', 'g3']
        
        mod1 = gauss1 
        mod2 = gauss1+gauss2
        mod3 = gauss1+gauss2+gauss3

        if modName == 'g1':
            return mod1,gName1
        if modName == 'g2':
            return mod2,gName2
        if modName == 'g3':
            return mod3,gName3

    def lineModDef(self,cfg_par,wave,y,lineInfo):

        gName = cfg_par['gFit']['modName']
        for i in xrange(0,len(lineInfo['ID'])):
            
            #waveMin =  np.log(lineInfo['Wave'][i] - lineInfo['lineRangeAng'][i])
            #waveMax =  np.log(lineInfo['Wave'][i] + lineInfo['lineRangeAng'][i])
           
            waveAmpIn1Min = np.log(lineInfo['Wave'][i]-lineInfo['cenRangeAng'][i])
            indexMin = int(np.where(abs(wave-waveAmpIn1Min)==abs(wave-waveAmpIn1Min).min())[0]) 
    
            waveAmpIn1Max = np.log(lineInfo['Wave'][i]+lineInfo['cenRangeAng'][i])
            indexMax = int(np.where(abs(wave-waveAmpIn1Max)==abs(wave-waveAmpIn1Max).min())[0])
                       
            Gmod = GaussianModel()

            gauss1 = GaussianModel(prefix='g1ln'+str(i)+'_')

            sigmaIn1 = lineInfo['deltaSigmaAng_In1'][i]
            ampIn1 = np.max(y[indexMin:indexMax])*max(2.220446049250313e-16, sigmaIn1)/0.3989423
            smallWave = wave[indexMin:indexMax]
            cenIn1 = smallWave[np.argmax(y[indexMin:indexMax])]

            if i == 0:
                pars = gauss1.make_params()
                pars['g1ln'+str(i)+'_'+'sigma'].set(value=sigmaIn1,
                    min=sigmaIn1/5.,max=sigmaIn1*5.)
                pars['g1ln'+str(i)+'_'+'center'].set(value=cenIn1,
                min=waveAmpIn1Min,max=waveAmpIn1Max)
                pars['g1ln'+str(i)+'_'+'amplitude'].set(value=ampIn1,min=0,max=None)
                mod = gauss1
          
            else:
            
                pars.update(gauss1.make_params())    
                pars['g1ln'+str(i)+'_'+'center'].set(value=cenIn1,
                        min=waveAmpIn1Min,max=waveAmpIn1Max)
                pars['g1ln'+str(i)+'_'+'amplitude'].set(value=ampIn1,min=0,max=None)
                if cfg_par['gFit']['fixSigma'] == True:
                    pars['g1ln'+str(i)+'_'+'sigma'].set(expr='g1ln'+str(0)+'_'+'sigma')
                else:
                    pars['g1ln'+str(i)+'_'+'sigma'].set(value=sigmaIn1,
                    min=sigmaIn1/5.,max=sigmaIn1*5.)    

                mod += gauss1            
            
            if gName != 'g1':

                Gmod = GaussianModel()
                gauss2 = GaussianModel(prefix='g2ln'+str(i)+'_')
                pars.update(gauss2.make_params())

                if i == 0:
                    sigmaIn2 = pars['g1ln'+str(i)+'_'+'sigma'] +lineInfo['deltaSigmaAng_12'][i]

                    pars['g2ln'+str(i)+'_'+'sigma'].set(value=sigmaIn2,min=sigmaIn2/5.,max=sigmaIn2*5.)
                elif cfg_par['gFit']['fixSigma'] == True:
                    pars['g2ln'+str(i)+'_'+'sigma'].set(expr='g2ln'+str(0)+'_'+'sigma')
                else:
                    pars['g2ln'+str(i)+'_'+'sigma'].set(value=sigmaIn2,min=sigmaIn2/5.,max=sigmaIn2*5.)
                
                ampIn2 = ampIn1*cfg_par['gFit']['dltAmp12']               
                cenIn2Pos = cenIn1 + lineInfo['deltaVAng_12'][i]
                cenIn2Neg = cenIn1 - lineInfo['deltaVAng_12'][i]
              
                if cfg_par['gFit']['fixCentre'] == True and i >0:

                    pars['g2ln'+str(i)+'_'+'center'].set(value=cenIn2Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i],vary=True)
                    
                    pars.add(name='g2ln'+str(i)+'Split_'+'center', expr='g2ln'+str(0)+'_'+'center - g1ln'+str(0)+'_'+'center')
                    pars.add(name='g2ln'+str(i)+'Pos_'+'center', value=cenIn2Pos,
                        max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i],min=cenIn1, vary=True)
                    pars.add(name='g2ln'+str(i)+'Neg_'+'center', value=cenIn2Neg,max=cenIn1,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i], vary=True)

                    pars['g2ln'+str(i)+'_'+'center'].set(expr='g2ln'+str(i)+'Pos_center if g2ln'+str(i)+'Split_center >= 0 else g2ln'+str(i)+'Neg_center' )                    

                else:
                    pars['g2ln'+str(i)+'_'+'center'].set(value=cenIn2Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i])                   

                pars['g2ln'+str(i)+'_'+'amplitude'].set(value=ampIn2,min=0,max=None)
                

                mod += gauss2


                if gName == 'g3':

                    Gmod = GaussianModel()

                    gauss3 = GaussianModel(prefix='g3ln'+str(i)+'_')

                    pars.update(gauss3.make_params())

                    
                    if i == 0:
                        sigmaIn3 = pars['g1ln'+str(i)+'_'+'sigma'] + lineInfo['deltaSigmaAng_13'][i]
                        pars['g3ln'+str(i)+'_'+'sigma'].set(value=sigmaIn3,min=sigmaIn3/5.,max=sigmaIn3*5.)
                    elif cfg_par['gFit']['fixSigma'] == True:
                        pars['g3ln'+str(i)+'_'+'sigma'].set(expr='g3ln'+str(0)+'_'+'sigma')
                    else:
                        pars['g3ln'+str(i)+'_'+'sigma'].set(value=sigmaIn3,min=sigmaIn3/5.,max=sigmaIn3*5.)

                    ampIn3 = ampIn1*cfg_par['gFit']['dltAmp13']
                    
                    cenIn3Pos = cenIn1 + lineInfo['deltaVAng_13'][i]
                    cenIn3Neg = cenIn1 - lineInfo['deltaVAng_13'][i]
                    
                    pars['g3ln'+str(i)+'_'+'center'].set(value=cenIn3Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i])
                
                    if cfg_par['gFit']['fixCentre'] == True and i >0:

                        pars['g3ln'+str(i)+'_'+'center'].set(value=cenIn3Pos,
                            min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i],vary=True)
                        
                        pars.add(name='g3ln'+str(i)+'Split_'+'center', expr='g3ln'+str(0)+'_'+'center - g1ln'+str(0)+'_'+'center')
                        pars.add(name='g3ln'+str(i)+'Pos_'+'center', value=cenIn3Pos,
                            max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i],min=cenIn1, vary=True)
                        pars.add(name='g3ln'+str(i)+'Neg_'+'center', value=cenIn3Neg,max=cenIn1,
                            min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i], vary=True)

                        pars['g3ln'+str(i)+'_'+'center'].set(expr='g3ln'+str(i)+'Pos_center if g3ln'+str(i)+'Split_center >= 0 else g3ln'+str(i)+'Neg_center' )                    

                    else:
                        pars['g3ln'+str(i)+'_'+'center'].set(value=cenIn3Pos,
                            min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i])                   

                    pars['g3ln'+str(i)+'_'+'amplitude'].set(value=ampIn3,min=0,max=None)

                    mod += gauss3

        #pars.pretty_print()
        return mod,pars

    def gFit(self,cfg_par):
        
        key = 'general'

        workDir = cfg_par[key]['workdir']
        
        #open line lineList
        lineInfo = tP.openLineList(cfg_par)
        diffusion = 1e-5
        
        #open table for bins
        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo = tP.openTablesPPXF(cfg_par,workDir+cfg_par[key]['tableBinName'],
            workDir+cfg_par[key]['tableSpecName'])
        #open datacube
        f = fits.open(workDir+cfg_par[key]['dataCubeName'])
        hh = f[0].header
        dd = f[0].data
        



        #define x-axis array
        lambdaMin = np.log(cfg_par['gFit']['lambdaMin'])
        lambdaMax = np.log(cfg_par['gFit']['lambdaMax'])


        idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
        idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0] )


        Ydim = dd.shape[1]
        Xdim = dd.shape[2]
        
        binID, binArr, fitResArr, lineArr = tP.makeInputArrays(cfg_par,lineInfo, Xdim, Ydim)
       
        counter = 0
        #for j in xrange(205,208):
        #    for i in xrange(250,252):
        for j in xrange(0,dd.shape[1]):
            for i in xrange(0,dd.shape[2]):
                y = dd[idxMin:idxMax,j,i]

                waveCut = wave[idxMin:idxMax]
                #check if spectrum is not empty                   
                if np.sum(y)>0:

                    gMod,gPars = self.lineModDef(cfg_par,waveCut,y,lineInfo)

                    # identify voronoi bin
                    xVal = xAxis[i]
                    yVal = yAxis[j]
                    
                    index = np.where((vorBinInfo['X'] < (xVal+pxSize/2.+diffusion)) & 
                    ((xVal-pxSize/2.-diffusion) < vorBinInfo['X']) & (vorBinInfo['Y'] < (yVal+pxSize/2.+diffusion)) & 
                    ((yVal-pxSize/2.-diffusion) < vorBinInfo['Y']))
                    
                    if np.sum(index)>0: 
                        binArr = tP.updateBinArray(cfg_par,binArr,vorBinInfo,index,i,j,counter)
                        binIDName = binArr['BIN_ID'][counter]     
                    else:
                        #fitResArr = np.delete(fitResArr,counter,0)
                        #lineArr = np.delete(lineArr,counter,0)  
                        counter+=1
                        continue
                    
                    #check if it is first time in bin
                    if binIDName not in binID[:,:] and np.sum(index)>0:
 
                        binID[j,i] = binIDName
                        noiseVec = noiseBin[binIDName][:]

                        # FIT
                        result = gMod.fit(y, gPars, x=waveCut)
                        save_modelresult(result, self.cfg_par['general']['modNameDir']+str(binIDName)+'_'+cfg_par['gFit']['modName']+'.sav')

                        fitResArr = tP.updateFitArray(cfg_par,fitResArr,result,binIDName,counter)
                        lineArr = tP.updateLineArray(cfg_par,lineArr,result,lineInfo,binIDName,counter)

                        #plot Fit
                        if cfg_par['gPlot']['enable'] == True:
                        #self.plotSpecFit(waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
                            sP.plotLineZoom(cfg_par,waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
                          
                counter+=1
        
        match_indices = np.where(binArr['BIN_ID'] == 0.0)[0]
        binArr = np.delete(binArr,match_indices,0)                                
        match_indices = np.where(fitResArr['BIN_ID'] == 0.0)[0]
        fitResArr = np.delete(fitResArr,match_indices,0)                                
        match_indices = np.where(lineArr['BIN_ID'] == 0)[0]
        lineArr = np.delete(lineArr,match_indices,0)                                

                #print 'end_for'
        tP.saveOutputTable(cfg_par, binArr, fitResArr, lineArr)
    
        print('''\t+---------+\n\t gFit done\n\t+---------+''')
    
        return 0

    def gPlot(self,cfg_par):
        
        key = 'general'

        workDir = cfg_par[key]['workdir']
        
        #open line lineList
        lineInfo = tP.openLineList(cfg_par)
        diffusion = 1e-5
        
        #open table for bins
        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo = tP.openTablesPPXF(cfg_par,workDir+cfg_par[key]['tableBinName'],
            workDir+cfg_par[key]['tableSpecName'])
        #open datacube
        f = fits.open(workDir+cfg_par[key]['dataCubeName'])
        hh = f[0].header
        dd = f[0].data
        

        modNameDir = cfg_par[key]['runNameDir']+'/myMods/'
        cfg_par[key]['modNameDir'] = modNameDir
        if not os.path.exists(modNameDir):
            os.mkdir(modNameDir)

        #define x-axis array
        lambdaMin = np.log(cfg_par['gFit']['lambdaMin'])
        lambdaMax = np.log(cfg_par['gFit']['lambdaMax'])


        idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
        idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0] )


        Ydim = dd.shape[1]
        Xdim = dd.shape[2]
        
        binID, binArr, fitResArr, lineArr = tP.makeInputArrays(cfg_par,lineInfo, Xdim, Ydim)
       
        counter = 0
        print 'culo'
        #for j in xrange(205,208):
        #    for i in xrange(250,252):
        for j in xrange(0,dd.shape[1]):
            for i in xrange(0,dd.shape[2]):
                y = dd[idxMin:idxMax,j,i]

                waveCut = wave[idxMin:idxMax]
                #check if spectrum is not empty                   
                if np.sum(y)>0:

                    # identify voronoi bin
                    xVal = xAxis[i]
                    yVal = yAxis[j]
                    
                    index = np.where((vorBinInfo['X'] < (xVal+pxSize/2.+diffusion)) & 
                    ((xVal-pxSize/2.-diffusion) < vorBinInfo['X']) & (vorBinInfo['Y'] < (yVal+pxSize/2.+diffusion)) & 
                    ((yVal-pxSize/2.-diffusion) < vorBinInfo['Y']))
                    
                    if np.sum(index)>0: 
                        binIDName = binArr['BIN_ID'][counter]     
                    else:
                        #fitResArr = np.delete(fitResArr,counter,0)
                        #lineArr = np.delete(lineArr,counter,0)  
                        counter+=1
                        continue
                    print binIDName
                    print binID
                    #check if it is first time in bin
                    if binIDName not in binID[:,:] and np.sum(index)>0:
 
                        binID[j,i] = binIDName
                        noiseVec = noiseBin[binIDName][:]

                        # FIT
                        result = load_modelresult(modNameDir+str(binIDName)+'_'+cfg_par['gFit']['modName']+'.sav')
                        print result
                        #plot Fit
                        if cfg_par['gPlot']['enable'] == True:
                            print 'i'
                        #self.plotSpecFit(waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
                            sP.plotLineZoom(cfg_par,waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
                          
                counter+=1

        print('''\t+---------+\n\t gPlot done\n\t+---------+''')

        return 0
