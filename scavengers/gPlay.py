#!/usr/bin/env python
import os, sys
import yaml

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult
from lmfit.model import load_modelresult


from numba import njit, prange

from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np
import numpy.ma as ma

import shutil

import gufo as gf
import cvPlay
import specPlot
import tPlay

#gf = gufo()
cvP = cvPlay.convert()
sP = specPlot.specplot()
tP = tPlay.tplay()

class gplay(object):
    

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

        dLambda = cvP.specRes(cfg_par)

        gName = cfg_par['gFit']['modName']
        for i in xrange(0,len(lineInfo['ID'])):
            
            if lineInfo['Wave'][i] == 6547.96:
                kk = i
            if lineInfo['Wave'][i] == 6716.31:
                zz = i
            else:
                kk = 0
                zz = 0

            waveInRed = cfg_par['general']['redshift']*lineInfo['Wave'][i]+lineInfo['Wave'][i]

            indexWaveInRed = int(np.where(abs(np.exp(wave)-waveInRed)==abs(np.exp(wave)-waveInRed).min())[0])
            dLIn = dLambda[indexWaveInRed]
            dLIn = np.log(waveInRed+dLIn/2.)-np.log(waveInRed-dLIn/2.)

            indexWaveIn = int(np.where(abs(np.exp(wave)-lineInfo['Wave'][i]).min())[0])
            waveIn = np.exp(wave[indexWaveIn])
            
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
                pars.add(name = 'Wintln'+str(i), value=dLIn,vary=False)
                pars.add(name = 'g1intln'+str(i), value=sigmaIn1,
                    min=sigmaIn1/5.,max=sigmaIn1*5.,vary=True)
                
                pars['g1ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g1intln'+str(i)+',2))')
                pars['g1ln'+str(i)+'_'+'center'].set(value=cenIn1,
                min=waveAmpIn1Min,max=waveAmpIn1Max,vary=True)

                pars.add(name='lineWave'+str(i),value=lineInfo['Wave'][i],vary=False)
                pars.add(name='cenDist',expr ='((exp(g1ln'+str(i)+'_'+'center)-lineWave'+str(i)+')/lineWave'+str(0)+')*2.99792458e8/1e3',vary=False)

                pars['g1ln'+str(i)+'_'+'amplitude'].set(value=ampIn1,min=0,max=None)
                mod = gauss1
          
            else:
                pars.update(gauss1.make_params())    
                
                if cfg_par['gFit']['fixCentre'] == True:
                    pars.add(name='lineWave'+str(i),value=lineInfo['Wave'][i],vary=False)

                    pars.add(name='cenDistAng'+str(i),expr ='log((cenDist*1e3*lineWave'+str(i)+'*1e-10)/2.99792458e8/1e-10+lineWave'+str(i)+')')

                    cenDist = cvP.lambdaVRad(np.exp(pars['g1ln'+str(0)+'_'+'center']),lineInfo['Wave'][0])
                    cenDistAng = np.log(cvP.vRadLambda(cenDist,lineInfo['Wave'][i]))
                    #pars.add(name='cenDistln'+str(i), value=cenDistAng,vary=False)
                    pars['g1ln'+str(i)+'_'+'center'].set(expr='cenDistAng'+str(i))
                else:
                    pars['g1ln'+str(i)+'_'+'center'].set(value=cenIn1,
                    min=waveAmpIn1Min,max=waveAmpIn1Max,vary=True) 

                pars['g1ln'+str(i)+'_'+'amplitude'].set(value=ampIn1,min=0,max=None)
                
                if lineInfo['Wave'][i] == 6583.34:
                    ampMin = pars['g1ln'+str(kk)+'_'+'height'] * 1./cfg_par['gFit']['ampRatioNII']
                    ampIn1 = np.max(y[indexMin:indexMax])
                    pars['g1ln'+str(i)+'_'+'height'].set(value=ampIn1,min=ampMin,max=None,vary=True)

                if lineInfo['Wave'][i] == 6730.68:
                    ampMin = pars['g1ln'+str(zz)+'_'+'height'] * cfg_par['gFit']['ampRatioSII']
                    ampIn1 = np.max(y[indexMin:indexMax])
                    pars['g1ln'+str(i)+'_'+'height'].set(value=ampIn1,min=ampMin,max=None,vary=True)
               

                pars.add(name = 'Wintln'+str(i), value=dLIn,vary=False)
                if cfg_par['gFit']['fixSigma'] == True:
                    pars.add(name = 'g1intln'+str(i), expr='g1intln'+str(0))  
                    pars['g1ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g1intln'+str(i)+',2))')
                else:
                    pars.add(name = 'g1intln'+str(i), value=sigmaIn1,
                    min=sigmaIn1/5.,max=sigmaIn1*5.,vary=True)
                    pars['g1ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g1intln'+str(i)+',2))')
 
#                    pars['g1ln'+str(i)+'_'+'sigma'].set(value=sigmaIn1,
#                    min=sigmaIn1/5.,max=sigmaIn1*5.)    

                mod += gauss1            
            
            if gName != 'g1':

                Gmod = GaussianModel()
                gauss2 = GaussianModel(prefix='g2ln'+str(i)+'_')
                pars.update(gauss2.make_params())

                if i == 0:
                    sigmaIn2 = pars['g1ln'+str(i)+'_'+'sigma'] +lineInfo['deltaSigmaAng_12'][i]
                #    pars['g2ln'+str(i)+'_'+'sigma'].set(value=sigmaIn2,min=sigmaIn2/5.,max=sigmaIn2*5.)
                    pars.add('g2intln'+str(i), value=sigmaIn2,
                    min=sigmaIn2/5.,max=sigmaIn2*5.,vary=True)
                    pars['g2ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g2intln'+str(i)+',2))')
                elif cfg_par['gFit']['fixSigma'] == True:
                    #pars['g2ln'+str(i)+'_'+'sigma'].set(expr='g2ln'+str(0)+'_'+'sigma')
                    pars.add(name = 'g2intln'+str(i), expr='g2intln'+str(0))  
                    pars['g2ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g2intln'+str(i)+',2))')
                else:
                    pars.add(name= 'g2intln'+str(i), value=sigmaIn2,
                    min=sigmaIn2/5.,max=sigmaIn2*5.,vary=True)
                    pars['g2ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g2intln'+str(i)+',2))')

#                else:
#                    pars['g2ln'+str(i)+'_'+'sigma'].set(value=sigmaIn2,min=sigmaIn2/5.,max=sigmaIn2*5.)
                
                ampIn2 = ampIn1*cfg_par['gFit']['dltAmp12']               
                cenIn2Pos = cenIn1 + lineInfo['deltaVAng_12'][i]
                cenIn2Neg = cenIn1 - lineInfo['deltaVAng_12'][i]

                if i == 0:
                    pars['g2ln'+str(i)+'_'+'center'].set(value=cenIn2Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i])

                elif cfg_par['gFit']['posCentre'] == True and i >0:
                    pars['g2ln'+str(i)+'_'+'center'].set(value=cenIn2Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i],vary=True)
                    pars.add(name='g2ln'+str(i)+'Split_'+'center', expr='g2ln'+str(0)+'_'+'center - g1ln'+str(0)+'_'+'center')
                    pars.add(name='g2ln'+str(i)+'Pos_'+'center', value=cenIn2Pos,
                        max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i],min=cenIn1, vary=True)
                    pars.add(name='g2ln'+str(i)+'Neg_'+'center', value=cenIn2Neg,max=cenIn1,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i], vary=True)
                    pars['g2ln'+str(i)+'_'+'center'].set(expr='g2ln'+str(i)+'Pos_center if g2ln'+str(i)+'Split_center >= 0 else g2ln'+str(i)+'Neg_center' )                    
                
                elif cfg_par['gFit']['fixCentre'] == True and i >0:
                    
                    cenDist = cvP.lambdaVRad(np.exp(pars['g2ln'+str(0)+'_'+'center']),lineInfo['Wave'][0])
                    cenDistAng = np.log(cvP.vRadLambda(cenDist,lineInfo['Wave'][i]))
                    pars.add(name='cenDistln'+str(i), value=cenDistAng,vary=False)
                    pars['g2ln'+str(i)+'_'+'center'].set(expr='cenDistln'+str(i))

                else:
                    pars['g2ln'+str(i)+'_'+'center'].set(expr='g2ln'+str(0)+'_'+'center')                   

                pars['g2ln'+str(i)+'_'+'amplitude'].set(value=ampIn2,min=0,max=None)
                

                mod += gauss2


                if gName == 'g3':

                    Gmod = GaussianModel()

                    gauss3 = GaussianModel(prefix='g3ln'+str(i)+'_')

                    pars.update(gauss3.make_params())

                    
                    if i == 0:
                        sigmaIn3 = pars['g1ln'+str(i)+'_'+'sigma'] + lineInfo['deltaSigmaAng_13'][i]
                        pars.add(name = 'g3intln'+str(i)+'_'+'sigma', value=sigmaIn3,
                            min=sigmaIn3/5.,max=sigmaIn3*5.,vary=True)               
                        pars['g3ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g3intln'+str(i)+'_'+'sigma,2))')
                        #pars['g3ln'+str(i)+'_'+'sigma'].set(value=sigmaIn3,min=sigmaIn3/5.,max=sigmaIn3*5.)
                    
                    elif cfg_par['gFit']['fixSigma'] == True:
                    #    pars['g3ln'+str(i)+'_'+'sigma'].set(expr='g3ln'+str(0)+'_'+'sigma')
                        pars.add(name = 'g3intln'+str(i)+'_'+'sigma', expr='g3intln'+str(0)+'_'+'sigma')  
                        pars['g3ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g3intln'+str(i)+'_'+'sigma,2))')
                    else:
                        pars.add(name = 'g3intln'+str(i)+'_'+'sigma',value=sigmaIn3,
                        min=sigmaIn3/5.,max=sigmaIn3*5.,vary=True)
                        pars['g3ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g3intln'+str(i)+'_'+'sigma,2))')

                    
                #    else:
                #        pars['g3ln'+str(i)+'_'+'sigma'].set(value=sigmaIn3,min=sigmaIn3/5.,max=sigmaIn3*5.)

                    ampIn3 = ampIn1*cfg_par['gFit']['dltAmp13']
                    
                    cenIn3Pos = cenIn1 + lineInfo['deltaVAng_13'][i]
                    cenIn3Neg = cenIn1 - lineInfo['deltaVAng_13'][i]
                    
                    pars['g3ln'+str(i)+'_'+'center'].set(value=cenIn3Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i])
                
                    if cfg_par['gFit']['posCentre'] == True and i >0:
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
        f = fits.open(cfg_par[key]['cubeDir']+cfg_par[key]['dataCubeName'])
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
        for j in prange(0,dd.shape[1]):
            for i in prange(0,dd.shape[2]):
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
                        save_modelresult(result, cfg_par['general']['modNameDir']+str(binIDName)+'_'+cfg_par['gFit']['modName']+'.sav')
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
        f = fits.open(cfg_par[key]['cubeDir']+cfg_par[key]['dataCubeName'])
        hh = f[0].header
        dd = f[0].data
        

        modNameDir = cfg_par[key]['modNameDir']

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
                        fitResArr = np.delete(fitResArr,counter,0)
                        lineArr = np.delete(lineArr,counter,0)  
                        counter+=1
                        continue
                    #check if it is first time in bin
                    if binIDName not in binID[:,:] and np.sum(index)>0:
 
                        binID[j,i] = binIDName
                        noiseVec = noiseBin[binIDName][:]

                        # FIT
                        result = load_modelresult(cfg_par[key]['modNameDir']+'_'+cfg_par['gFit']['modName']+'.sav')
                        cfg_par['gPlot']['loadModel'] = True

                        #plot Fit
                        if cfg_par['gPlot']['enable'] == True:
                        #self.plotSpecFit(waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
                            sP.plotLineZoom(cfg_par,waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
                          
                counter+=1

        f.close()
        print('''\t+---------+\n\t gPlot done\n\t+---------+''')

        return 0

    def plotSingleBin(self,cfg_par,binID):

        key = 'general'
        workDir = cfg_par[key]['workdir']
        cubeDir = cfg_par[key]['cubeDir']


        #open datacube
        f = fits.open(cubeDir+cfg_par[key]['dataCubeName'])
        dd = f[0].data
        f.close()
        
        lineInfo = tP.openLineList(cfg_par)


        #open table for bins
        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo = tP.openTablesPPXF(cfg_par,workDir+cfg_par[key]['tableBinName'],
            workDir+cfg_par[key]['tableSpecName'])

        hdul = fits.open(cfg_par['general']['outTableName'])
        tabGen = hdul['BinInfo'].data


        lambdaMin = np.log(cfg_par['gFit']['lambdaMin'])
        lambdaMax = np.log(cfg_par['gFit']['lambdaMax'])
        idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
        idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0])

        print idxMin, idxMax
        idxTable = int(np.where(tabGen['BIN_ID'] == int(binID))[0][0])
        
        print idxTable
        print int(tabGen['PixY'][idxTable]),int(tabGen['PixX'][idxTable])

        y = dd[idxMin:idxMax,int(tabGen['PixY'][idxTable]),int(tabGen['PixX'][idxTable])]

        waveCut = wave[idxMin:idxMax]

        result = load_modelresult(cfg_par[key]['modNameDir']+str(binID)+'_'+cfg_par['gFit']['modName']+'.sav')
        
        cfg_par['gPlot']['loadModel'] = True

        noiseVec = noiseBin[int(binID)][:]

        #plot Fit
        #self.plotSpecFit(waveCut, y,result,noiseVec[idxMin:idxMax],i,j,tab,vorBinInfo[index])
        sP.plotLineZoom(cfg_par,waveCut, y,result,noiseVec[idxMin:idxMax],int(tabGen['PixX'][idxTable]),int(tabGen['PixY'][idxTable]),lineInfo,tabGen[:][idxTable])

        print('''\t+---------+\n\t bin Plotted\n\t+---------+''')


        return 0
