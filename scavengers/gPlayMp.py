#!/usr/bin/env python3.6

import os, sys
import yaml

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult
from lmfit.model import load_modelresult

#import pickle4MPplay
import multiprocessing as mp
from multiprocessing import Queue, Manager, Process

#needed to change the pickle
#ctx = mp.get_context()
#ctx.reducer = pickle4MPplay.Pickle4Reducer()
#mp.context._default_context.reducer = pickle4MPplay.Pickle2Reducer()
#import multiprocessing.connection
#import ctypes
#from tqdm import tqdm


from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np
import numpy.ma as ma

import cvPlay, specPlot, tPlay
import specPlot
import tPlay

cvP = cvPlay.convert()
sP = specPlot.specplot()
tP = tPlay.tplay()


def workerGFitMp(cfg_par,dd,rank,nprocs):
    '''
    Worker of gaussian loop for multiprocess. 
    Fits in the same time N spectra of the Voronoi binned datacube.
    N is the number of cores used.

    Parameters:
        cfg_par - configuration file
        dd - datacube
        rank - starting point of for loop
        nprocs - number of cores where to send spectraters of the multi-gaussian model
    
    Returns:
        binArr
        fitResArr
        lineArr
    '''   

    key = 'general'

    workDir = cfg_par[key]['workdir']
    
    lineInfo = tP.openLineList(cfg_par)
    diffusion = 1e-5
    
    wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec = tP.openVorLineOutput(cfg_par,cfg_par['general']['outVorLineTableName'],
        cfg_par['general']['outVorSpectra'])
        
    ubins = np.unique(vorBinInfo['BIN_ID'])

    counter = 0

    binArr, fitResArr, lineArr = tP.makeInputArraysMP(cfg_par,lineInfo, vorBinInfo)

    for ii in range(rank,len(ubins), nprocs):
        
        counter,binArr,fitResArr,lineArr = gFitMp(cfg_par,lineInfo,vorBinInfo,wave,dd,noiseBin,counter,ii,ubins,binArr,
            fitResArr,lineArr)

    match_indices = np.where(binArr['BIN_ID'] == 0.0)[0]
    binArr = np.delete(binArr,match_indices,0)                                
    match_indices = np.where(fitResArr['BIN_ID'] == 0.0)[0]
    fitResArr = np.delete(fitResArr,match_indices,0)                                
    match_indices = np.where(lineArr['BIN_ID'] == 0)[0]
    lineArr = np.delete(lineArr,match_indices,0)    

    return binArr, lineArr, fitResArr  

# def show_prog(q, total_bytes):
#     prog = tqdm(total=total_bytes, desc="Total", unit='B', unit_scale=True)
#     while 1:
#         try:
#             to_add = q.get(timeout=1)
#             prog.n += to_add
#             prog.update(0)
#             if prog.n >= total_bytes:
#                 break
#         except:
#             continue


def gFitMp(cfg_par,lineInfo,vorBinInfo,wave,dd,noiseBin,counter,ii,ubins,binArr,fitResArr,lineArr):
    '''
    Fits the spectrum of a single Voronoi bin a model of 1 or multiple gaussians

    Parameters:
        cfg_par - configuration file
        lineInfo - array with columns of lineList.txt file
        vorBinInfo - array of table of voronoi binned datacube
        wave - x-axis of spectrum in log(lambda) units
        dd - datacube
        noiseBin - voronoi binned noise cube
        counter - step of multiprocess loop
        ubins - array of BINID
        binArr - array where to store information about bins
        fitResArr - array where to store information about fit results
        lineArr - array where to store best fit parameters of the multi-gaussian model
    
    Returns:
        binArr
        fitResArr
        lineArr
    '''

    match_bin = np.where(ubins[ii]==vorBinInfo['BIN_ID'])[0]

    index = match_bin[0]

    j = int(vorBinInfo['PixY'][index])
    i = int(vorBinInfo['PixX'][index])

    lambdaMin = np.log(cfg_par['gFit']['lambdaMin'])
    lambdaMax = np.log(cfg_par['gFit']['lambdaMax'])


    idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
    idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0])

    Ydim = int(dd.shape[1])
    Xdim = int(dd.shape[2])

    y = dd[idxMin:idxMax,j,i]
    waveCut = wave[idxMin:idxMax]

    if not np.isnan(y).all():
        gMod,gPars = lineModDefMp(cfg_par,waveCut,y,lineInfo)

        if np.sum(index)>0: 
            binArr = tP.updateBinArray(cfg_par,binArr,vorBinInfo,index,i,j,counter)
            binIDName = vorBinInfo['BIN_ID'][index]     
        else:
            fitResArr = np.delete(fitResArr,counter,0)
            lineArr = np.delete(lineArr,counter,0)  
            counter+=1
            return counter,binArr,fitResArr,lineArr

        noiseVec = noiseBin[binIDName][:]

        # FIT
        result = gMod.fit(y, gPars, x=waveCut)
        
        save_modelresult(result, cfg_par['general']['modNameDir']+str(binIDName)+'_'+cfg_par['gFit']['modName']+'.sav')
        fitResArr = tP.updateFitArray(cfg_par,fitResArr,result,binIDName,counter)
        lineArr = tP.updateLineArray(cfg_par,waveCut,lineArr,result,noiseVec[idxMin],lineInfo,binIDName,counter)
        
        #plot Fit
        if cfg_par['gPlot']['enable'] == True:
            sP.plotLineZoom(cfg_par,waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])

    counter+=1
    
    return counter,binArr,fitResArr,lineArr

# this is likely not used
# def update(i, ans):
#     '''
#     Defines composite gaussian model for lmfit fitting
    
#     Parameters:
#         wave: x-axis of spectrum in log(lambda) to be fitted
#         y: spectrum
#         lineInfo: lines to be fitted, wavelenght and other parameters defined in lineList.txt
#     Returns:
#         model of 1 or multiple gaussians to use for fit
#     '''  
#     print(i,ans)
#     res[i] = ans  # put answer into correct index of result list
#     pbar.update()



def lineModDefMp(cfg_par,wave,y,lineInfo):
    '''
    Defines composite gaussian model for lmfit fitting
    
    Parameters:
        wave: x-axis of spectrum in log(lambda) to be fitted
        y: spectrum
        lineInfo: lines to be fitted, wavelenght and other parameters defined in lineList.txt
    Returns:
        model of 1 or multiple gaussians to use for fit in gFitMp module
    '''
    dLambda = cvP.specRes(cfg_par)

    gName = cfg_par['gFit']['modName']
    for i in range(0,len(lineInfo['ID'])):
        
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
        
        waveAmpIn1Min = np.log(lineInfo['Wave'][i]-lineInfo['cenRangeAng'][i])
        indexMin = int(np.where(abs(wave-waveAmpIn1Min)==abs(wave-waveAmpIn1Min).min())[0]) 

        waveAmpIn1Max = np.log(lineInfo['Wave'][i]+lineInfo['cenRangeAng'][i])
        indexMax = int(np.where(abs(wave-waveAmpIn1Max)==abs(wave-waveAmpIn1Max).min())[0])
                   
        Gmod = GaussianModel()

        gauss1 = GaussianModel(prefix='g1ln'+str(i)+'_')

        sigmaMin = lineInfo['deltaSigmaAng_Min'][i]
        sigmaMaxG1 = lineInfo['deltaSigmaAng_MaxG1'][i]
        sigmaMaxG2 = lineInfo['deltaSigmaAng_MaxG2'][i]
        sigmaMaxG3 = lineInfo['deltaSigmaAng_MaxG3'][i]
        ampIn1 = np.max(y[indexMin:indexMax])*max(2.220446049250313e-16, sigmaMin)/0.3989423
        smallWave = wave[indexMin:indexMax]
        cenIn1 = smallWave[np.argmax(y[indexMin:indexMax])]


        if i == 0:

            pars = gauss1.make_params()
            pars.add(name = 'Wintln'+str(i), value=dLIn,vary=False)
            
            #if gName=='g1':
            pars.add(name = 'g1intln'+str(i), value=sigmaMin*5.,vary=True,min=sigmaMin,max=sigmaMaxG1)
            #else:
            #    pars.add(name = 'g1intln'+str(i), value=sigmaMin*5.,vary=True,min=sigmaMin,max=sigmaMaxG1)

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

                pars['g1ln'+str(i)+'_'+'center'].set(expr='cenDistAng'+str(i))
            else:
                pars['g1ln'+str(i)+'_'+'center'].set(value=cenIn1,
                min=waveAmpIn1Min,max=waveAmpIn1Max,vary=True) 

            pars['g1ln'+str(i)+'_'+'amplitude'].set(value=ampIn1,min=0,max=None)
           
            pars.add(name = 'Wintln'+str(i), value=dLIn,vary=False)
            if cfg_par['gFit']['fixSigma'] == True:
                pars.add(name = 'g1intln'+str(i), expr='g1intln'+str(0))  
                pars['g1ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g1intln'+str(i)+',2))')
            else:
                pars['g1ln'+str(i)+'_'+'sigma'].set(value=sigmaMin,
                    min=sigmaMin,max=sigmaMaxG1,vary=True)    

            mod += gauss1            
        
        if gName != 'g1':

            Gmod = GaussianModel()
            gauss2 = GaussianModel(prefix='g2ln'+str(i)+'_')
            pars.update(gauss2.make_params())
            cenIn2Pos = cenIn1

            ampIn2 = ampIn1*cfg_par['gFit']['dltAmp12']      
            pars['g2ln'+str(i)+'_'+'amplitude'].set(value=ampIn2,min=0,max=None)

            if i == 0:
                sigmaIn2 = pars['g1ln'+str(i)+'_'+'sigma'] +lineInfo['deltaSigmaAng_12'][i]
            #    pars['g2ln'+str(i)+'_'+'sigma'].set(value=sigmaIn2,min=sigmaIn2/5.,max=sigmaIn2*5.)
                pars.add('g2intln'+str(i), value=sigmaMin*5,
                    min=pars['g1intln'+str(i)].value,vary=True,max=sigmaMaxG2)
                pars['g2ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g2intln'+str(i)+',2))')

                pars['g2ln'+str(i)+'_'+'center'].set(value=cenIn2Pos,
                    min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i],vary=True)

                pars.add(name='lineWave'+str(i),value=lineInfo['Wave'][i],vary=False)
            
                pars.add(name='cenDistg2',expr ='((exp(g2ln'+str(i)+'_'+'center)-lineWave'+str(i)+')/lineWave'+str(0)+')*2.99792458e8/1e3',vary=False)


            else:
                if cfg_par['gFit']['fixSigma'] == True:
                    #pars['g2ln'+str(i)+'_'+'sigma'].set(expr='g2ln'+str(0)+'_'+'sigma')
                    pars.add(name = 'g2intln'+str(i), expr='g2intln'+str(0))  
                    pars['g2ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g2intln'+str(i)+',2))')
                
                if cfg_par['gFit']['fixCentre'] == True:
                    pars.add(name='lineWave'+str(i),value=lineInfo['Wave'][i],vary=False)

                    pars.add(name='cenDistAngg2'+str(i),expr ='log((cenDistg2*1e3*lineWave'+str(i)+'*1e-10)/2.99792458e8/1e-10+lineWave'+str(i)+')')
                   
                    cenDist = cvP.lambdaVRad(np.exp(pars['g2ln'+str(0)+'_'+'center']),lineInfo['Wave'][0])
                    cenDistAng = np.log(cvP.vRadLambda(cenDist,lineInfo['Wave'][i]))
                    pars.add(name='cenDistln'+str(i), value=cenDistAng,vary=False)
                    pars['g2ln'+str(i)+'_'+'center'].set(expr='cenDistAngg2'+str(i))

            mod += gauss2


            if gName == 'g3':

                Gmod = GaussianModel()

                gauss3 = GaussianModel(prefix='g3ln'+str(i)+'_')

                pars.update(gauss3.make_params())
                
                if i == 0:
                    sigmaIn3 = pars['g1ln'+str(i)+'_'+'sigma'] + lineInfo['deltaSigmaAng_13'][i]
                    pars.add(name = 'g3intln'+str(i)+'_'+'sigma', value=sigmaMin*6.,min=pars['g1intln'+str(i)].value,vary=True,max=sigmaMax*3)               
                    pars['g3ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g3intln'+str(i)+'_'+'sigma,2))')
                elif cfg_par['gFit']['fixSigma'] == True:
                    pars.add(name = 'g3intln'+str(i)+'_'+'sigma', expr='g3intln'+str(0)+'_'+'sigma')  
                    pars['g3ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g3intln'+str(i)+'_'+'sigma,2))')
                else:
                    pars.add(name = 'g3intln'+str(i)+'_'+'sigma',value=sigmaIn3,
                    min=sigmaIn3/5.,max=sigmaIn3*5.,vary=True)
                    pars['g3ln'+str(i)+'_'+'sigma'].set(expr='sqrt(pow(Wintln'+str(i)+',2)+pow(g3intln'+str(i)+'_'+'sigma,2))')

                ampIn3 = ampIn1*cfg_par['gFit']['dltAmp13']
                
                cenIn3Pos = cenIn1 + lineInfo['deltaVAng_13'][i]
                cenIn3Neg = cenIn1 - lineInfo['deltaVAng_13'][i]
                
                if i == 0:
                    pars['g3ln'+str(i)+'_'+'center'].set(value=cenIn3Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i])
            

                if cfg_par['gFit']['fixCentre'] == True and i>0:
                    cenDist = cvP.lambdaVRad(np.exp(pars['g3ln'+str(0)+'_'+'center']),lineInfo['Wave'][0])
                    cenDistAng = np.log(cvP.vRadLambda(cenDist,lineInfo['Wave'][i]))
                    pars.add(name='cenDistln'+str(i), value=cenDistAng,vary=False)
                    pars['g3ln'+str(i)+'_'+'center'].set(expr='cenDistln'+str(i))

                else:
                    pars['g3ln'+str(i)+'_'+'center'].set(value=cenIn3,
                    min=waveAmpIn1Min,max=waveAmpIn1Max,vary=True)                   
                    mod += gauss3

    return mod,pars

def log_result(result):
    '''
    This is called whenever foo_pool(i) returns a result.
    result_list is modified only by the main process, not the pool workers.
    '''
    result = [p.get() for p in result]
    
           
def main(cfg_par):
    result_list = []

    #freeze_support()
    key = 'general'

    workDir = cfg_par[key]['workdir']
    
    #open line lineList
    #lineInfo = tP.openLineList(cfg_par)
    #diffusion = 1e-5
    
    #open table for bins
    wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec = tP.openVorLineOutput(cfg_par,cfg_par['general']['outVorLineTableName'],
        cfg_par['general']['outVorSpectra'])

    #open datacube
    f = fits.open(cfg_par['general']['outVorLines'])
    hh = f[0].header
    dd = f[0].data


    if mp.current_process().name == "MainProcess":

        if cfg_par[key]['nCores']:
            nprocs = cfg_par[key]['nCores']
        else:
            nprocs = mp.cpu_count()
 
        inputs = [(cfg_par,dd,rank, nprocs) for rank in range(nprocs)]
        print('\n\t         +++\t\t  multiProcess\t\t +++')


        pool = mp.Pool(processes=nprocs)

        multi_result = [pool.apply_async(workerGFitMp, args=(inp)) for inp in inputs]
        
        result = [p.get() for p in multi_result]

        binArr = np.array(result[0][0])
        lineArr = np.array(result[0][1])
        fitResArr = np.array(result[0][2])

        for i in range(1,nprocs):
            binArr = np.hstack([binArr,np.array(result[i][0])])
            lineArr = np.hstack([lineArr,np.array(result[i][1])])
            fitResArr = np.hstack([fitResArr,np.array(result[i][2])])

        tP.saveOutputTable(cfg_par, binArr, fitResArr, lineArr)