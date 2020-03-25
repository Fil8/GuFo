#!/usr/bin/env python3.6

import os, sys
import yaml

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult
from lmfit.model import load_modelresult

import multiprocessing as mp
from multiprocessing import Queue, Manager, Process

from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np
import numpy.ma as ma


import tPlay,cvPlay

tP = tPlay.tplay()
cvP = cvPlay.convert()

def workerAncels(cfg_par,lines,wave,lineInfo,dLambda,sigmaCen,rank,nprocs):

 
    counter = 0

    for ii in range(rank,len(lines['BIN_ID']), nprocs):
        
        counter,sigmamaCen = w50Centroid(cfg_par,lines,wave,lineInfo,dLambda,sigmaCen,counter)
    

    print(sigmaCen['BIN_ID'])
    match_indices = np.where(sigmaCen['BIN_ID'] == 0.0)[0]
    sigmaCen = np.delete(sigmaCen,match_indices,0)                                
    print(sigmaCen)
    return sigmaCen  
    

def w50Centroid(cfg_par,lines,wave,lineInfo,dLambda,sigmaCen,counter):

    #tt=Table([lines['BIN_ID']])
    binID = lines['BIN_ID'][counter]
    modName = cfg_par['gFit']['modName']
    result = load_modelresult(cfg_par['general']['modNameDir']+str(binID)+'_'+cfg_par['gFit']['modName']+'.sav')
    comps = result.eval_components()
    
    sigmaCen['BIN_ID'][counter]=binID
    for ii in range(0,len(lineInfo['ID'])):
        lambdaRest = lineInfo['Wave'][ii]

        waveInRed = cfg_par['general']['redshift']*lambdaRest+lambdaRest
        indexWaveInRed = int(np.where(abs(np.exp(wave)-waveInRed)==abs(np.exp(wave)-waveInRed).min())[0])
        dLIn1 = dLambda[indexWaveInRed]  
        
        lineName = str(lineInfo['Name'][ii])+str(int(lineInfo['Wave'][ii]))
        if '[' in lineName:
            lineName = lineName.replace("[", "")
            lineName = lineName.replace("]", "")
        
        lineNameID = lineName
        lineThresh = float(lineInfo['SNThresh'][ii])


        lineFit= comps['g1ln'+str(ii)+'_']+comps['g2ln'+str(ii)+'_']
        centreFit=result.params['g1ln'+str(ii)+'_center']
        ampSpax = lines[modName+'-AmpSpax_'+lineNameID][counter]

        if ampSpax >= lineThresh:
            height =np.max(lineFit)

            #print(dLIn,dLIn1,centreFit)

            centroidG1 = lines['g1_Centre_'+str(lineNameID)][counter]
            #print(centroidG1)
            centroidG2 = lines['g2_Centre_'+str(lineNameID)][counter]
            ampG1 = lines['g1_Amp_'+str(lineNameID)][counter]
            ampG2 = lines['g2_Amp_'+str(lineNameID)][counter]
            centroidG1Lambda = cvP.vRadLambda(centroidG1,lambdaRest)
            centroidG2Lambda = cvP.vRadLambda(centroidG2,lambdaRest)

            if np.sum(lineFit)==0:
                w50 = np.nan
            else:    
                indexHeight =np.abs(lineFit-height).argmin()

                height /=2.

                lineFitLeft = lineFit[:indexHeight]
                lineFitRight = lineFit[indexHeight:]
                indexWaveLeft = (np.abs(lineFitLeft-height)).argmin()
                indexWaveRight = (np.abs(lineFitRight-height)).argmin()+indexHeight

                waveDist = np.exp(wave[indexWaveRight])-np.exp(wave[indexWaveLeft])
                #print(indexWaveLeft,indexWaveRight,np.exp(wave[indexWaveLeft]),
                #    np.exp(wave[indexWaveRight]),waveDist)

                sigmaLambda = waveDist/(2*np.sqrt(2*np.log(2)))
                sigmaInt = np.sqrt(np.power(sigmaLambda,2)-np.power(dLIn1,2))    
                sigmaCen['sigma_'+lineName][counter] =  cvP.lambdaVRad(lambdaRest+sigmaInt,lambdaRest)
                
            centroidToT = np.divide(np.sum([np.multiply(centroidG1Lambda,ampG1)+np.multiply(centroidG2Lambda,ampG2)]),
                np.sum([ampG1,ampG2]))          

            sigmaCen['centroid_'+lineName][counter] = cvP.lambdaVRad(centroidToT,lambdaRest)
            #print(centroidToTVel[j],sigmaInt,sigma[j])
        else:

            sigmaCen['sigma_'+lineName][counter]=np.nan
            sigmaCen['centroid_'+lineName][counter] = np.nan

    counter +=1
    
    return counter, sigmaCen

def main(cfg_par):
    
    result_list = []

    #freeze_support()
    key = 'general'

    workDir = cfg_par[key]['workdir']

    lineNameID=[]
    lineThresh=[]


    lineInfo = tP.openLineList(cfg_par)



    hdul = fits.open(cfg_par['general']['outTableName'])
    lines = hdul['LineRes_'+cfg_par['gFit']['modName']].data 

    dLambda = cvP.specRes(cfg_par)


    wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec = tP.openVorLineOutput(cfg_par,cfg_par['general']['outVorLineTableName'],
        cfg_par['general']['outVorSpectra'])

    sigmaNameList = []
    sigmaList = []
    sigmaNameList.append('BIN_ID')
    sigmaList.append('i4')
    for ii in range(0,len(lineInfo['ID'])):
        lineRest = lineInfo['Wave'][ii]
        lineName = str(lineInfo['Name'][ii])+str(int(lineInfo['Wave'][ii]))
        if '[' in lineName:
            lineName = lineName.replace("[", "")
            lineName = lineName.replace("]", "")

        sigmaNameList.append('sigma_'+lineName)
        sigmaList.append('f8')
        sigmaNameList.append('centroid_'+lineName)
        sigmaList.append('f8')
    
    sigmaCenArr = np.zeros([len(lines['BIN_ID'])],dtype={'names':(sigmaNameList), 'formats':(sigmaList)})


    print(sigmaCenArr.dtype.names)

    if mp.current_process().name == "MainProcess":

        if cfg_par[key]['nCores']:
            nprocs = cfg_par[key]['nCores']
        else:
            nprocs = mp.cpu_count()

        inputs = [(cfg_par,lines,wave,lineInfo,dLambda,sigmaCenArr,rank, nprocs) for rank in range(nprocs)]
        print('''\t+---------+\n\t sigma & centroid\n\t+---------+''')

        pool = mp.Pool(processes=nprocs)
        multi_result = [pool.apply_async(workerAncels, args=(inp)) for inp in inputs]        
        result = [p.get() for p in multi_result]
        print(result[0][0])
        print(result)
        sigmaCen = np.asarray(result)
    
        #for i in range(1,nprocs):
        #    sigmaCen = np.hstack([sigmaCen,np.array(result[i][0])])

        print(sigmaCen)
        print(sigmaCen.shape)
        tP.saveAncelsTable(cfg_par, sigmaCen)

        print('''\t+---------+\n\t sigma & centroid done\n\t+---------+''')



