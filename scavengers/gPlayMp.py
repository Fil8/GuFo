#!/usr/bin/env python3.6

import os, sys
import yaml

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult
from lmfit.model import load_modelresult

import pickle4MPplay
import multiprocessing as mp
from multiprocessing import Queue, Manager, Process

ctx = mp.get_context()
ctx.reducer = pickle4MPplay.Pickle4Reducer()
mp.context._default_context.reducer = pickle4MPplay.Pickle4Reducer()

import ctypes
from tqdm import tqdm


from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np
import numpy.ma as ma

import shutil

#import gufo as gf
import cvPlay
import specPlot
import tPlay

import util as pretty

#gf = gufo()
cvP = cvPlay.convert()
sP = specPlot.specplot()
tP = tPlay.tplay()


#lock = mp.Lock()
def workerGFitMp(cfg_par,dd,rank,nprocs):
    

    #existing_shm = shared_memory.SharedMemory(name=binIDShareName)
    key = 'general'

    workDir = cfg_par[key]['workdir']
    
    #open line lineList
    lineInfo = tP.openLineList(cfg_par)
    diffusion = 1e-5
    
    wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec = tP.openVorLineOutput(cfg_par,cfg_par['general']['outVorLineTableName'],
        cfg_par['general']['outVorSpectra'])
        
    ubins = np.unique(vorBinInfo['BIN_ID'])
    #for i in range(0,len(ubins)):
    #    print(ubins[i])
    #
    counter = 0

    binArr, fitResArr, lineArr = tP.makeInputArraysMP(cfg_par,lineInfo, vorBinInfo)
    #for j in range(205,208):
    #    for i in range(250,252):
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

def show_prog(q, total_bytes):
    prog = tqdm(total=total_bytes, desc="Total", unit='B', unit_scale=True)
    while 1:
        try:
            to_add = q.get(timeout=1)
            prog.n += to_add
            prog.update(0)
            if prog.n >= total_bytes:
                break
        except:
            continue


def gFitMp(cfg_par,lineInfo,vorBinInfo,wave,dd,noiseBin,counter,ii,ubins,binArr,fitResArr,lineArr):

            #pretty.printProgress( ii, len(ubins), barLength=50 )

        #for i in range(rank, nsteps, nprocs):
            #0,dd.shape[2]):
            match_bin = np.where(ubins[ii]==vorBinInfo['BIN_ID'])[0]

            index = match_bin[0]
        
            j = int(vorBinInfo['PixY'][index])
            i = int(vorBinInfo['PixX'][index])

            lambdaMin = np.log(cfg_par['gFit']['lambdaMin'])
            lambdaMax = np.log(cfg_par['gFit']['lambdaMax'])


            idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
            idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0])

            #dd.shape[2]=12
            Ydim = int(dd.shape[1])
            Xdim = int(dd.shape[2])

            y = dd[idxMin:idxMax,j,i]
            waveCut = wave[idxMin:idxMax]



                #check if spectrum is not empty                   
                #if np.nansum(y)>0:
            if not np.isnan(y).all():
                gMod,gPars = lineModDefMp(cfg_par,waveCut,y,lineInfo)

                # identify voronoi bin
                #xVal = xAxis[i]
                #yVal = yAxis[j]
                
                #index = np.where((vorBinInfo['X'] < (xVal+pxSize/2.+diffusion)) & 
                #((xVal-pxSize/2.-diffusion) < vorBinInfo['X']) & (vorBinInfo['Y'] < (yVal+pxSize/2.+diffusion)) & 
                #((yVal-pxSize/2.-diffusion) < vorBinInfo['Y']))
                
                if np.sum(index)>0: 
                    binArr = tP.updateBinArray(cfg_par,binArr,vorBinInfo,index,i,j,counter)
                    binIDName = vorBinInfo['BIN_ID'][index]     
                else:
                    fitResArr = np.delete(fitResArr,counter,0)
                    lineArr = np.delete(lineArr,counter,0)  
                    counter+=1
                    return counter,binArr,fitResArr,lineArr

                
                #print(binIDName,rank)

                #binIDShare = np.ndarray((Ydim, Xdim), dtype=np.int16, buffer=existing_shm.buf)
                #print(np.sum(binIDShare))
                #lock.acquire()    
                #print('inBin')
                #check if it is first time in bin
                
                #if int(binIDName) not in binIDShare[:,:] and np.sum(index)>0:
                #binIDShare[j,i] = binIDName
                #print('io')
                #print(binIDShare[j,i], str(rank))

                # Always release the lock!
                #lock.release()

                noiseVec = noiseBin[binIDName][:]

                # FIT
                result = gMod.fit(y, gPars, x=waveCut)
                
                save_modelresult(result, cfg_par['general']['modNameDir']+str(binIDName)+'_'+cfg_par['gFit']['modName']+'.sav')
                fitResArr = tP.updateFitArray(cfg_par,fitResArr,result,binIDName,counter)
                lineArr = tP.updateLineArray(cfg_par,lineArr,result,noiseVec[idxMin],lineInfo,binIDName,counter)
                
                #plot Fit
                if cfg_par['gPlot']['enable'] == True:
                #self.plotSpecFit(waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
                    sP.plotLineZoom(cfg_par,waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
            #else:
        
                    #lock.release()
        
            #existing_shm.close()

            counter+=1
            
            return counter,binArr,fitResArr,lineArr


def update(i, ans):
    # note: input comes from async `wrapMyFunc`
    print(i,ans)
    res[i] = ans  # put answer into correct index of result list
    pbar.update()


def lineModDefMp(cfg_par,wave,y,lineInfo):

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

            elif cfg_par['gFit']['fixCentre'] == True and i >0:
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

# def create_shared_block(dd):

#     smallD = np.zeros([dd.shape[1],dd.shape[2]])

#     a = np.ones(shape=smallD.shape, dtype=np.int8)

#     shm = shared_memory.SharedMemory(create=True, size=smallD.nbytes)
    
#     np_array = np.ndarray(smallD.shape, dtype=np.int16, buffer=shm.buf)
#     np_array[:] = a[:]  
    
#     return shm, np_array
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result = [p.get() for p in result]
    #print(result)

    #result_list.append(result)

           
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

    #binIDShare_base = mp.Array(ctypes.c_int, dd.shape[1]*dd.shape[2])
    #binIDShare = np.ctypeslib.as_array(binIDShare_base.get_obj())
    #binIDShare = binIDShare.reshape(*(dd.shape[1],dd.shape[2]))
    #print binIDShare.shape
    #lock = mp.Lock()
    #define x-axis array

        #binIDShare, array = create_shared_block(dd)
        if cfg_par[key]['nCores']:
            nprocs = cfg_par[key]['nCores']
        else:
            nprocs = mp.cpu_count()
        #nprocs -= 2

        inputs = [(cfg_par,dd,rank, nprocs) for rank in range(nprocs)]
        #print inputs 
        print('''\t+---------+\n\t going to process\n\t+---------+''')
        
        #test = mpFit()

        #processes = []
        #for i in range(nprocs):
            #print(inputs[i])
        #    _process = mp.Process(target=gFitMp, args=(inputs[i]))
        #    processes.append(_process)
        #    _process.start()

        #for _process in processes:
        #    _process.join()

        pool = mp.Pool(processes=nprocs)
        #pathos = pp.ProcessPool(nprocs)
        
        #print(_process[0])
        #   freeze_support()
        
        #q = mp.Queue()
        #workers = [mp.Process(target=gFitMp, args = (binIDShare.name,cfg_par,lineInfo,dd,rank, nprocs, 12,))
        #for rank in range(nprocs)]

        #for p in workers:
        #    p.start()
        #for p in workers:
        #    p.join()
        #while not q.empty():
        #    print(q.get())
        #manager = Manager()
        #queue = manager.Queue()
        #:
        #for inp in range(10):

    # tqdm.write(str(a))
        #for i in range(pbar.total):
        #pbar = tqdm(total=len(inputs))
        #:
        multi_result = [pool.apply_async(workerGFitMp, args=(inp)) for inp in inputs]
        
        result = [p.get() for p in multi_result]

        #progress = Process(target=show_prog, args=(queue, final_bytes))

        #progress.start()

        #multi_result = [pool.apply_async(gFitMp, inputs[0]) ]
        #tqdm.write('scheduled')
        #print(result)
        #sys.exit(0)
        
        #pool.close()
        #pool.start()
        #pool.join()
        #progress.join()
        
        #]                
        #print(result.shape())
        #result = [p.get() for p in _process]                
        #binIDShare.close()
        #binIDShare.unlink()

        binArr = np.array(result[0][0])
        lineArr = np.array(result[0][1])
        fitResArr = np.array(result[0][2])

        #print binArr.shape
        for i in range(1,nprocs):
            #print result[i][0].shape
            binArr = np.hstack([binArr,np.array(result[i][0])])
            lineArr = np.hstack([lineArr,np.array(result[i][1])])
            fitResArr = np.hstack([fitResArr,np.array(result[i][2])])

            #fitResArr = np.vstack([fitResArr,np.array(result[i][1])])
            #lineArr = np.vstack([lineArr,np.array(result[i][2])])

        #for i in range(0,nprocs):
        #    binArr = np.concatenate(result[0])
        #print fitResArr.shape
        #print  result[0][0]
        tP.saveOutputTable(cfg_par, binArr, fitResArr, lineArr)

        print('''\t+---------+\n\t gFit done\n\t+---------+''')

#if __name__ == '__main__':  
#    main()