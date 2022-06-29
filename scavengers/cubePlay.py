#!/usr/bin/env python3.6

import os, sys, shutil
import yaml
from math import floor,ceil
from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult
from lmfit.model import load_modelresult

from scipy.ndimage import map_coordinates

from reproject import reproject_interp as rp

from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from reproject import reproject_interp, reproject_exact


from matplotlib import pyplot as plt

#from MontagePy.main    import *
#from MontagePy.archive import *

import numpy as np
#import numpy.ma as ma

import tPlay,cvPlay,bptPlot,momPlot, headPlay

tP = tPlay.tplay()
cvP = cvPlay.convert()
mPl = momPlot.MOMplot()
hP = headPlay.headplay()

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
        header['CUNIT3'] = 'm/s'
        header['BUNIT'] = 'Jy/beam'
        return header

    def delete_header(self, header, keyword):
        '''
        from SoFiA
        '''

        if keyword in header:
            del header[keyword]
            return True
        return False

    def delete_3rd_axis(self,header):
        '''
        from SoFiA
        '''
        self.delete_header(header, "NAXIS3")
        self.delete_header(header, "CTYPE3")
        self.delete_header(header, "CRPIX3")
        self.delete_header(header, "CRVAL3")
        self.delete_header(header, "CDELT3")
        self.delete_header(header, "CUNIT3")
        if 'WCSAXES' in header:
            header['WCSAXES'] = 2

        if 'NAXIS' in header:
            header['NAXIS'] = 2

        return

    def measRMS(self,cubePath):
        niter = 4 # nr of noise measurements
        clip = 4  # clipping at each iteration (in units of std)

        medNoise = []


        unit_dict = {
            'vrad'   :   [1e+3, 'Velocity [km/s]'],
            'freq'   :   [1e+6, 'Frequency [MHz]'],
            'velo-hel': [1e+3, 'Velocity [km/s]'],
            'vopt-f2w': [1e+3, 'Velocity [km/s]'],
        }

        f=fits.open(cubePath)
        head=f[0].header
        cube=f[0].data
        cube = np.squeeze(cube)
        cube=cube[:,cube.shape[1]//4:3*cube.shape[1]//4,cube.shape[2]//4:3*cube.shape[2]//4]
        iter=1
        while iter<niter:
          std=np.nanmedian(np.nanstd(cube,axis=(1,2)))
          #print('# iter {1:d}, median std = {0:.2e} Jy/beam'.format(std,iter))
          cube[np.abs(cube)>clip*std]=np.nan
          iter+=1

        noise=np.nanstd(cube,axis=(1,2))
        #print('# iter {1:d}, median std = {0:.2e} Jy/beam'.format(np.nanmedian(noise),iter))

        freqs=(np.arange(head['naxis3'])-(head['crpix3']-1))*head['cdelt3']+head['crval3']


        noisePath = os.path.dirname(cubePath)
        baseName = os.path.basename(cubePath)

        if not os.path.exists(noisePath+'/noisePlots/'):
            os.mkdir(noisePath+'/noisePlots/')
        noisePath=noisePath+'/noisePlots/'

        plt.plot(freqs/unit_dict[head['CTYPE3'].lower()][0],noise*1e+3,'k-')
        plt.axhline(y=np.median(noise)*1e+3,linestyle=':',color='k')
        plt.xlabel(unit_dict[head['CTYPE3'].lower()][1])
        plt.ylabel('Noise [mJy/beam]')
        outPath=noisePath+baseName.replace('.fits','_rms.png')
        plt.savefig(outPath)
        plt.clf()

        return ((np.nanmedian(noise))),outPath

    def makeBFLineCube(self,cfg_par):

        cubeletsDir = cfg_par['general']['cubeletsDir']
        cubeDir = cfg_par['general']['cubeDir']
        momDir = cfg_par['general']['momDir']

        if cfg_par['bestFitSel']['BFcube']['rotationID'] == True:
            modelCube = fits.open(cfg_par['bestFitSel']['BFcube']['modelCube'])
            mdC = modelCube[0].data
            rotMoM = np.zeros([mdC.shape[1],mdC.shape[2]])*np.nan



        hdul = fits.open(cfg_par['general']['outTableName'])
        tabGen = hdul['BININFO'].data

        residuals = hdul['Residuals_'+cfg_par['gFit']['modName']].data
        ancels = hdul['Ancels'+cfg_par['gFit']['modName']].data
        rotArr = np.empty(len(ancels['BIN_ID']))

        bF = np.array(residuals['bestFit'],dtype=int)

        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec = tP.openVorLineOutput(cfg_par,cfg_par['general']['outVorLineTableName'],
            cfg_par['general']['outVorSpectra'])

        lineInfo = tP.openLineList(cfg_par)
        lineInfoAll = lineInfo.copy()
        index = np.where(lineInfo['BFCUbe'] == 0)
        fltr =  np.array(index)[0]
        lineInfo.remove_rows(list(fltr))

        lineNameStr = str(lineInfo['Name'][0])

        if '[' in lineNameStr:
            lineName = lineNameStr.replace("[", "")
            lineName = lineName.replace("]", "")
            lineName = lineName+str(int(lineInfo['Wave'][0]))
        else:
            lineName = lineNameStr+str(int(lineInfo['Wave'][0]))
        
        lineNameName = lineNameStr+str(int(lineInfo['Wave'][0]))

        lineNameStr=np.array([lineNameStr+str(int(lineInfo['Wave'][0]))],dtype='a16')
        
        lineNamesStrAll = np.empty(len(lineInfoAll),dtype='a16')
        for ii in range(0,len(lineInfoAll['ID'])):
        #for ii in range(0,1):

            lineNameStrAll = str(lineInfoAll['Name'][ii])

            if '[' in lineNameStrAll:
                lineNameAll = lineNameStrAll.replace("[", "")
                lineNameAll = lineNameAll.replace("]", "")
                lineNameAll= lineNameAll+str(int(lineInfoAll['Wave'][ii]))
            else:
                lineNameAll = lineNameStrAll+str(int(lineInfoAll['Wave'][ii]))

            lineNamesStrAll[ii] = str(lineNameStrAll+str(int(lineInfoAll['Wave'][ii])))

        indexLine =  np.where(lineNamesStrAll == lineNameStr)[0]
        if cfg_par['bestFitSel']['BFcube']['rotationID'] == True:
            modName='BF'
            f = fits.open(cubeDir+'fitCube_'+modName+'.fits')
            dd = f[0].data
            hh = f[0].header
        else:
            modName='g2'
            f = fits.open(cubeDir+'fitCube_'+modName+'.fits')
            dd = f[0].data
            hh = f[0].header            
        mom = fits.open(momDir+'g2/mom0_tot-'+lineName+'.fits')
        mm=mom[0].data
        indexFltr = np.broadcast_to(np.isnan(mm), dd.shape)
        dd[indexFltr] = np.nan
        #lambdaMin = np.log(cfg_par['gFit']['lambdaMin'])
        #lambdaMax = np.log(cfg_par['gFit']['lambdaMax'])
        #idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
        #idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0])

        lambdaMin = np.log(cfg_par['gFit']['lambdaMin'])
        lambdaMax = np.log(cfg_par['gFit']['lambdaMax'])
        idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
        idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0])

        wave=wave[idxMin:idxMax]

        dd=dd[idxMin:idxMax,:,:]

        velRangeMin = cvP.vRadLambda(-cfg_par['bestFitSel']['BFcube']['velRange'][0],
            lineInfo['Wave'][0] )-lineInfo['Wave'][0] 
        velRangeMax= cvP.vRadLambda(cfg_par['bestFitSel']['BFcube']['velRange'][1],
            lineInfo['Wave'][0] )-lineInfo['Wave'][0] 
                    
        waveMin =  np.log(lineInfo['Wave'][0] - velRangeMin)
        waveMax =  np.log(lineInfo['Wave'][0] + velRangeMax)
        

        idxMin1 = int(np.where(abs(wave-waveMin)==abs(wave-waveMin).min())[0]) 
        idxMax1 = int(np.where(abs(wave-waveMax)==abs(wave-waveMax).min())[0])

        wave=wave[idxMin1:idxMax1]
        
        waveAng=np.exp(wave)
        
        vel=cvP.lambdaVRad(waveAng,lineInfo['Wave'][0])+float(cfg_par['general']['velsys'])
        
        dd=dd[idxMin1:idxMax1,:,:]
        fitCube = np.empty([dd.shape[0],dd.shape[1],dd.shape[2]])
        fitCubeMask = np.zeros([dd.shape[0],dd.shape[1],dd.shape[2]])
        fitCubeMD = np.zeros([dd.shape[0],dd.shape[1],dd.shape[2]])

        fitCubeMaskInter = np.zeros([dd.shape[0],dd.shape[1],dd.shape[2]])
        vecSumMap = np.zeros([dd.shape[1],dd.shape[2]])*np.nan
        lenghtLineMap = np.zeros([dd.shape[1],dd.shape[2]])*np.nan

        for i in range(0,len(ancels['BIN_ID'])):
            
            match_bin = np.where(tabGen['BIN_ID']==ancels['BIN_ID'][i])[0]
            binID = ancels['BIN_ID'][i]
            indexBinIDres = np.where(residuals['BIN_ID']==binID)[0]
            if bF[indexBinIDres] == 0:
                modName = 'g1'
   
        
                for index in match_bin:
                    #print('figa')
                    if np.sum(~np.isnan(dd[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])])) != 0: 
           
                        result = load_modelresult(cfg_par['general']['runNameDir']+'models/'+modName+'/'+str(int(binID))+'_'+modName+'.sav')
                        comps = result.eval_components()                    
                        #if modName=='g1':
                        #elif modName =='g2':
                        fit = comps['g1ln'+str(indexLine[0])+'_']
                        #print('cazzo')
                        fitCube[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = fit[idxMin1:idxMax1]

                        if cfg_par['bestFitSel']['BFcube']['rotationID'] == True:
                            mdSpec = mdC[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]
                            
                            mdSpec[mdSpec!=0]=1.
                            mdSpec = np.flipud(mdSpec)
                            fitCubeMD[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = mdSpec
                            #centroid = ancels['centroid_'+lineName][i]
                            #width = ancels['w80_'+lineName][i]
                            fitSmall = fit[idxMin1:idxMax1]
                            idxPeak = np.nanargmax(fitSmall)
                            #print(idxPeak)
                            idxLeft = int(np.where(abs(fitSmall[:idxPeak]-10.)==abs(fitSmall[:idxPeak]-10.).min())[0]) 
                            idxRight = int(np.where(abs(fitSmall[idxPeak:]-10.)==abs(fitSmall[idxPeak:]-10.).min())[0]) +idxPeak

                            #print(idxLeft,idxRight)
                            #velMin = centroid-(width/2.)
                            #velMax = centroid+(width/2.)
                            #indexVelMin = int(np.where(abs(vel-velMin)==abs(vel-velMin).min())[0]) 
                            #indexVelMax = int(np.where(abs(vel-velMax)==abs(vel-velMax).min())[0]) 
                            
                            fitMask = np.zeros(len(fit[idxMin1:idxMax1]))
                            fitMaskIntercect = np.zeros(len(fit[idxMin1:idxMax1]))
                            #fitMask[indexVelMin:indexVelMax] = 1.
                            
                            fitMask[idxLeft:idxRight] = 1.
                            lenghtLine = np.count_nonzero(fitMask == 1.)
                            
                            fitCubeMask[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = fitMask
                            #print(fitMask,mdSpec)
                            vecCount = np.where((fitMask==1.)& (mdSpec==1.))
                            fitMaskIntercect[vecCount] = 1.
                            
                            vecSum = np.count_nonzero(fitMaskIntercect == 1.)
                            fitCubeMaskInter[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = fitMaskIntercect
                            
                            vecSumMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = vecSum
                            
                            lenghtLineMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lenghtLine/100.*cfg_par['bestFitSel']['BFcube']['rotationPercent']
                            lengthLineSpec= np.count_nonzero(mdSpec == 1.)
                            
                            if vecSum>(lengthLineSpec/100.*cfg_par['bestFitSel']['BFcube']['rotationPercent']):
                                rotArr[i]=1.
                                rotMoM[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]=1.

                            else:
                                rotArr[i]=0.
                                rotMoM[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]=0.
                    else:
                        fitCube[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        fitCubeMask[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        fitCubeMD[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan

                        fitCubeMaskInter[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        vecSumMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        lenghtLineMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        rotArr[i]=np.nan
                        
                        if cfg_par['bestFitSel']['BFcube']['rotationID'] == True:
                            rotMoM[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]=np.nan

                
            elif bF[indexBinIDres] == 1:
                modName = 'g2' #we consider only the first component as rotation
                #print(bF[i],residuals['BIN_ID'][i],residuals['bestFit'][i])
                #sys.exit(0)        
            #print('culo')
            #print(rotArr[i],print(ancels['BIN_ID'][i]))
                for index in match_bin:
                    if np.sum(~np.isnan(dd[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])])) != 0: 
                        result = load_modelresult(cfg_par['general']['runNameDir']+'models/'+modName+'/'+str(int(binID))+'_'+modName+'.sav')
                        comps = result.eval_components()
                        fit = comps['g1ln'+str(indexLine[0])+'_']+comps['g2ln'+str(indexLine[0])+'_']
                #fit = comps['g1ln'+str(indexLine[0])+'_'] 
                        fitCube[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = fit[idxMin1:idxMax1]           
                    else:
                        fitCube[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        fitCubeMask[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        fitCubeMD[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        fitCubeMaskInter[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        vecSumMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        lenghtLineMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                        rotArr[i]=np.nan
                        
                        if cfg_par['bestFitSel']['BFcube']['rotationID'] == True:
                            rotMoM[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]=np.nan


                if cfg_par['bestFitSel']['BFcube']['rotationID'] == True:

                    rotArr[i]=0.
                    rotMoM[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]=0.


            else:
                fitCube[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                fitCubeMask[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                fitCubeMD[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan

                fitCubeMaskInter[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                vecSumMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                lenghtLineMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                rotArr[i]=np.nan
                
                if cfg_par['bestFitSel']['BFcube']['rotationID'] == True:
                    rotMoM[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]=np.nan

        waveAng=np.exp(wave)

        header = self.makeHeader(cfg_par,lineInfo['Wave'][0],hh,waveAng)
        
        outCubelet = cubeletsDir+str(lineNameName)+'_BF.fits'            
        outCubeletMask = cubeletsDir+str(lineNameName)+'_BFMask.fits'        
        outCubeletMaskfl = cubeletsDir+str(lineNameName)+'_BFMaskInter.fits'        
        outCubeletMaskMD = cubeletsDir+str(lineNameName)+'_BFMD.fits'        

        fits.writeto(outCubelet,np.flip(fitCube,axis=0),header,overwrite=True)
        fits.writeto(outCubeletMask,np.flip(fitCubeMask,axis=0),header,overwrite=True)
        fits.writeto(outCubeletMaskfl,np.flip(fitCubeMaskInter,axis=0),header,overwrite=True)
        fits.writeto(outCubeletMaskMD,np.flip(fitCubeMD,axis=0),header,overwrite=True)

        if cfg_par['bestFitSel']['BFcube']['rotationID'] == True:

            outMomRot =  momDir+str(lineNameName)+'_RotMom.fits'        
            outMomSum =  momDir+str(lineNameName)+'_SumInter.fits'        
            outMomLength =  momDir+str(lineNameName)+'_LengthLine.fits'        
            outMomDiff =  momDir+str(lineNameName)+'_diffInter.fits'        

            if 'CUNIT3' in header:
                del header['CUNIT3']
            if 'CTYPE3' in header:
                del header['CTYPE3']
            if 'CDELT3' in header:
                del header['CDELT3']
            if 'CRVAL3' in header:  
                del header['CRVAL3']
            if 'CRPIX3' in header:
                del header['CRPIX3'] 
            if 'NAXIS3' in header:
                del header['NAXIS3']            
            header['WCSAXES'] = 2
            header['SPECSYS'] = 'topocent'
            header['BUNIT'] = 'Jy/beam'

            fits.writeto(outMomRot,rotMoM,header,overwrite=True)
            fits.writeto(outMomSum,vecSumMap,header,overwrite=True)
            fits.writeto(outMomLength,lenghtLineMap,header,overwrite=True)
            fits.writeto(outMomDiff,vecSumMap-lenghtLineMap,header,overwrite=True)

            t=Table(ancels)

            if 'RotMod' not in ancels.dtype.names: 
                t.add_column(Column(rotArr,name='RotMod'))
            else:
                t.replace_column('RotMod',Column(rotArr,name='RotMod'))        

            try:
                tt = Table(hdul['AncelsBF'].data)
                hdul['AncelsBF'] = fits.BinTableHDU(t.as_array(),name='AncelsBF')
            except KeyError as e:
                tt=fits.BinTableHDU.from_columns(t.as_array(),name='AncelsBF')   
                hdul.append(tt)    

        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)
            

        return

    def selectRotation(self,cfg_par):

        cubeletsDir = cfg_par['general']['cubeletsDir']
        momDir = cfg_par['general']['momDir']


        hdul = fits.open(cfg_par['general']['outTableName'])
        tabGen = hdul['BININFO'].data
        ancels = hdul['Ancels'+cfg_par['gFit']['modName']].data

        modelCube = fits.open(cfg_par['otherGasKinAnalysis']['rotID']['modelCube'])
        mdC = modelCube[0].data
        #indexFltrMod = mdC ==0.0
        #print(mdC[65,64,43])
        mdC[mdC>=float(cfg_par['otherGasKinAnalysis']['rotID']['maskValue'])] = 1.
        #print(mdC[65,64,43])

        mdC[mdC==0] = np.nan        

        maskCube = fits.open(cfg_par['otherGasKinAnalysis']['rotID']['gasMask'])
        maskD = np.array(maskCube[0].data,dtype=float)
       

        maskD[maskD!=0] = 1.
        maskD[maskD==0] = np.nan        

        rotMoM = np.zeros([mdC.shape[1],mdC.shape[2]])*np.nan

        rotArr = np.empty(len(ancels['BIN_ID']))*np.nan

        f = fits.open(cfg_par['otherGasKinAnalysis']['rotID']['gasCube'])
        dd = f[0].data
        header = f[0].header
        mom = fits.open(cfg_par['otherGasKinAnalysis']['rotID']['gasMoment'])
        mm=mom[0].data
        indexFltr = np.broadcast_to(mm==0.0, dd.shape)
        dd[indexFltr] = np.nan

        lineNameName = cfg_par['otherGasKinAnalysis']['Name']


        vel = ((np.linspace(1, dd.shape[0], dd.shape[0]) - header['CRPIX3']) 
            * header['CDELT3'] + header['CRVAL3'])/1e3
        velMin = cfg_par['otherGasKinAnalysis']['rotID']['velRange'][1]+float(cfg_par['general']['velsys'])
        velMax = cfg_par['otherGasKinAnalysis']['rotID']['velRange'][0]+float(cfg_par['general']['velsys'])
        idxMin1 = int(np.where(abs(vel-velMin)==abs(vel-velMin).min())[0]) 
        idxMax1 = int(np.where(abs(vel-velMax)==abs(vel-velMax).min())[0])        

        dd=dd[idxMin1:idxMax1,:,:]
        mdC=mdC[idxMin1:idxMax1,:,:]
        maskD=maskD[idxMin1:idxMax1,:,:]


        fitCube = np.empty([dd.shape[0],dd.shape[1],dd.shape[2]])
        fitCubeMask = np.zeros([dd.shape[0],dd.shape[1],dd.shape[2]])
        fitCubeMD = np.zeros([dd.shape[0],dd.shape[1],dd.shape[2]])

        fitCubeMaskInter = np.zeros([dd.shape[0],dd.shape[1],dd.shape[2]])
        vecSumMap = np.zeros([dd.shape[1],dd.shape[2]])*np.nan
        lenghtLineMap = np.zeros([dd.shape[1],dd.shape[2]])*np.nan


        for i in range(0,len(ancels['BIN_ID'])):
            
            match_bin = np.where(tabGen['BIN_ID']==ancels['BIN_ID'][i])[0]
                
            for index in match_bin:
                
                if np.sum(~np.isnan(dd[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])])) != 0 and ancels['sigma_'+lineNameName][i] !=np.nan: 

                    fitCube[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = dd[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]
                    fit = dd[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]
                    

                    mdSpec = mdC[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]
                    #if cfg_par['otherGasKinAnalysis']['rotID']['enable']==True:
                        #idxModFalse = mdSpec<float(cfg_par['otherGasKinAnalysis']['rotID']['maskValue'])
                        
                        #idxMod = mdSpec>=float(cfg_par['otherGasKinAnalysis']['rotID']['maskValue'])
                        #mdSpec[idxMod] =1.
                        #mdSpec[idxMod] =0.
                    #indexMod = np.logical_and(mdSpec!=0,mdSpec!=np.nan)

                    #mdSpec[indexMod ]=1.
#                    mdSpec = np.flipud(mdSpec)
 
                    fitCubeMD[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = mdSpec
                    #centroid = ancels['centroid_'+lineName][i]
                    #width = ancels['w80_'+lineName][i]
                    fitSmall = fit[:]
                    specMask = maskD[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]
                    velMinPeak = -100.+float(cfg_par['general']['velsys'])
                    velMaxPeak = +100.+float(cfg_par['general']['velsys'])
                    idxMinPeak = int(np.where(abs(vel-velMin)==abs(vel-velMin).min())[0]) 
                    idxMaxPeak = int(np.where(abs(vel-velMax)==abs(vel-velMax).min())[0])    
                    idxPeak = np.nanargmax(fitSmall[idxMinPeak:idxMaxPeak])
                    idxPeak1 = np.nanargmax(fitSmall)

                    #print(fitSmall[idxPeak],fitSmall[idxPeak1])
  

                    peakValue = fitSmall[idxPeak]
                    noiseValue= cfg_par['otherGasKinAnalysis']['rotID']['sigmaNoise']*float(cfg_par['otherGasKinAnalysis']['rotID']['noiseValue'])
                   
                    idxMask = np.logical_or(specMask==0.,specMask==np.nan)
                    idxMaskTrue = specMask==1.
    
                    #print(np.sum(idxMask))
                    fitSmall[idxMask] = np.nan

                    #idxLeft = int(np.where(abs(fitSmall[:idxPeak]-noiseValue)==abs(fitSmall[:idxPeak]-noiseValue).min())[0]) 
                    #idxRight = int(np.where(abs(fitSmall[idxPeak:]-noiseValue)==abs(fitSmall[idxPeak:]-noiseValue).min())[0]) +idxPeak

                    #print(idxLeft,idxRight)
                    #velMin = centroid-(width/2.)
                    #velMax = centroid+(width/2.)
                    #indexVelMin = int(np.where(abs(vel-velMin)==abs(vel-velMin).min())[0]) 
                    #indexVelMax = int(np.where(abs(vel-velMax)==abs(vel-velMax).min())[0]) 
                    
                    #fitMask = np.zeros(len(fit[:]))
                    #fitMask[idxMaskTrue] = 1.
                    fitMaskIntercect = np.zeros(len(fit[:]))
                    #fitMask[indexVelMin:indexVelMax] = 1.
                    
                    #fitMask[idxLeft:idxRight] = 1.
                    
                    vecMask = np.logical_and(fitSmall > noiseValue, specMask==1.)
                    
                    fitMaskOne = np.zeros(len(fit[:]))
                    fitMaskOne[vecMask]=1.

                    lenghtLine = np.nansum(fitMaskOne == 1.)
                    
                    #fitCubeMask[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = fitMask
                    fitCubeMask[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = fitMaskOne
                    #print(fitMask,mdSpec)
                    #vecCount = np.where((fitMask==1.)& (mdSpec!=0.0))
                    vecCount = np.where((fitMaskOne==1.)& (mdSpec==1.))
                    
                    fitMaskIntercect[vecCount] = 1.
                    
                    vecSum = np.count_nonzero(fitMaskIntercect == 1.)
                    #print(vecSum,lenghtLine,lenghtLine/100.*cfg_par['otherGasKinAnalysis']['rotID']['rotationPercent'])
                    fitCubeMaskInter[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = fitMaskIntercect
                    
                    vecSumMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = vecSum
                    
                    lenghtLineMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lenghtLine/100.*cfg_par['otherGasKinAnalysis']['rotID']['rotationPercent']
                    #lengthLineSpec= np.count_nonzero(mdSpec == 1.)

                    if vecSum>(lenghtLine/100.*cfg_par['otherGasKinAnalysis']['rotID']['rotationPercent']):
                        rotArr[i]=1.
                        rotMoM[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]=1.

                    else:
                        rotArr[i]=0.
                        rotMoM[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]=0.


                else:
                    fitCube[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                    fitCubeMask[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                    fitCubeMD[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan

                    fitCubeMaskInter[:,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                    vecSumMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan
                    lenghtLineMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = np.nan

        #waveAng=np.exp(wave)

        #header = self.makeHeader(cfg_par,lineInfo['Wave'][0],hh,waveAng)
        
        header['CRVAL3'] = (float(cfg_par['otherGasKinAnalysis']['rotID']['velRange'][0])+float(cfg_par['general']['velsys']))*1e3

        outCubelet = cubeletsDir+str(lineNameName)+'_BF.fits'            
        outCubeletMask = cubeletsDir+str(lineNameName)+'_BFMask.fits'        
        outCubeletMaskfl = cubeletsDir+str(lineNameName)+'_BFMaskInter.fits'        
        outCubeletMaskMD = cubeletsDir+str(lineNameName)+'_BFMD.fits'        

        fits.writeto(outCubelet,fitCube,header,overwrite=True)
        fits.writeto(outCubeletMask,fitCubeMask,header,overwrite=True)
        fits.writeto(outCubeletMaskfl,fitCubeMaskInter,header,overwrite=True)
        fits.writeto(outCubeletMaskMD,fitCubeMD,header,overwrite=True)

        outMomRot =  momDir+str(lineNameName)+'_RotMom.fits'        
        outMomSum =  momDir+str(lineNameName)+'_SumInter.fits'        
        outMomLength =  momDir+str(lineNameName)+'_LengthLine.fits'        
        outMomDiff =  momDir+str(lineNameName)+'_diffInter.fits'        

        if 'CUNIT3' in header:
            del header['CUNIT3']
        if 'CTYPE3' in header:
            del header['CTYPE3']
        if 'CDELT3' in header:
            del header['CDELT3']
        if 'CRVAL3' in header:  
            del header['CRVAL3']
        if 'CRPIX3' in header:
            del header['CRPIX3'] 
        if 'NAXIS3' in header:
            del header['NAXIS3']            
        header['WCSAXES'] = 2
        header['SPECSYS'] = 'topocent'
        header['BUNIT'] = 'Jy/beam'

        fits.writeto(outMomRot,rotMoM,header,overwrite=True)
        fits.writeto(outMomSum,vecSumMap,header,overwrite=True)
        fits.writeto(outMomLength,lenghtLineMap,header,overwrite=True)
        fits.writeto(outMomDiff,vecSumMap-lenghtLineMap,header,overwrite=True)

        t=Table(ancels)

        if 'RotMod' not in ancels.dtype.names: 
            t.add_column(Column(rotArr,name='RotMod'))
        else:
            t.replace_column('RotMod',Column(rotArr,name='RotMod'))        

        try:
            tt = Table(hdul['Ancels'+cfg_par['gFit']['modName']].data)
            hdul['AncelsBF'] = fits.BinTableHDU(t.as_array(),name='AncelsBF')
        except KeyError as e:
            tt=fits.BinTableHDU.from_columns(t.as_array(),name='AncelsBF')   
            hdul.append(tt)    

        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)
            

        return

    def makeLineCubes(self,cfg_par):

        cubeletsDir = cfg_par['general']['cubeletsDir']
        cubeDir = cfg_par['general']['cubeDir']
        

        modName = cfg_par['gFit']['modName']
        momDir = cfg_par['general']['momDir']


        if cfg_par['cubePlay']['cubelets']['cube'] == 'vorLine': 
            f = fits.open(cfg_par['general']['dataCubeName'])
            dd = f[0].data
            hh = f[0].header

        elif cfg_par['cubePlay']['cubelets']['cube'] == 'fitLine': 
            f = fits.open(cubeDir+'fitCube_'+modName+'.fits')
            dd = f[0].data
            hh = f[0].header
            mom = fits.open(momDir+'g1/mom0_g1-OIII5006.fits')
            mm=mom[0].data
            print(mm)
            indexFltr = np.broadcast_to(np.isnan(mm), dd.shape)
            print(indexFltr.shape)
            dd[indexFltr] = np.nan
            print(dd.shape)

        elif cfg_par['cubePlay']['cubelets']['cube'] == 'residuals':
            f = fits.open(cubeDir+'resCube_'+modName+'.fits')
            dd = f[0].data
            hh = f[0].header
            mom = fits.open(momDir+'g1/mom0_g1-OIII5006.fits')
            mm=mom[0].data
            indexFltr = np.broadcast_to(np.isnan(mm), dd.shape)
            dd[indexFltr] = np.nan
        else:
            f = fits.open(cfg_par['general']['outLines'])
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
            
            velRangeMin = cvP.vRadLambda(-cfg_par['cubePlay']['cubelets']['velRange'][0],
                lineInfo['Wave'][ii] )-lineInfo['Wave'][ii] 
            velRangeMax= cvP.vRadLambda(cfg_par['cubePlay']['cubelets']['velRange'][1],
                lineInfo['Wave'][ii] )-lineInfo['Wave'][ii] 
            #velRangeLeft = cvP.vRadLambda(cfg_par['cubePlay']['velRange'],
            #    lineInfo['Wave'][ii] )+lineInfo['Wave'][ii] 

            #waveMin =  np.log(lineInfo['Wave'][ii] - lineInfo['lineRangeAng'][ii])
            #waveMax =  np.log(lineInfo['Wave'][ii] + lineInfo['lineRangeAng'][ii])
            
            waveMin =  np.log(lineInfo['Wave'][ii] - velRangeMin)
            waveMax =  np.log(lineInfo['Wave'][ii] + velRangeMax)
            

            idxMin = int(np.where(abs(wave-waveMin)==abs(wave-waveMin).min())[0]) 
            idxMax = int(np.where(abs(wave-waveMax)==abs(wave-waveMax).min())[0] )
            
            dCbl = dd[idxMin:idxMax,:,:]

            waveAng=np.exp(wave[idxMin:idxMax])

            header = self.makeHeader(cfg_par,lineInfo['Wave'][ii],hh,waveAng)
            
            if cfg_par['cubePlay']['cubelets']['cube'] == 'vorLine': 
                outCubelet = cubeletsDir+str(lineNameStr)+'_measVor.fits'
            elif cfg_par['cubePlay']['cubelets']['cube'] == 'fitLine': 
                outCubelet = cubeletsDir+str(lineNameStr)+'_fit_'+modName+'.fits'   
            elif cfg_par['cubePlay']['cubelets']['cube'] == 'residuals':        
                outCubelet = cubeletsDir+str(lineNameStr)+'_res_'+modName+'.fits'        
            else:
                outCubelet = cubeletsDir+str(lineNameStr)+'_meas.fits'

            fits.writeto(outCubelet,np.flip(dCbl,axis=0),header,overwrite=True)

        return


    # def makeBFLineCube(self,cfg_par):



    def circPoint(self,alpha,beta,m,q,pa,ww,keyAlpha):
        a = -2.*alpha
        b = -2.*beta 
        if ww != 0.:
            r = (ww*1.5)/np.tan(np.radians(float(pa)))
        else:
            return alpha,beta
        c = np.power(alpha,2)+np.power(beta,2)-np.power(r,2)

        bb = 2.*m*q+a+b*m
        delta = np.sqrt(np.power(bb,2)-4.*(1.+np.power(m,2))*(np.power(q,2)+b*q+c))
        dueA = 2.*(1.+np.power(m,2))
        x = np.divide(-bb+delta,dueA)
        if keyAlpha == 'min':
            if x<alpha:
                x = np.divide(-bb-delta,dueA)
        elif keyAlpha == 'max':
            if x>alpha:
                x = np.divide(-bb-delta,dueA)
        y = m*x+q
        return x,y
    
    def findExtremes(self,ww,pa,Xc,Yc,shape):

        m = np.tan(np.radians(float(pa)))
        q = Yc-m*Xc
        x0 = -q/m
        
        y0 = m*shape[2]+q
        if m >= 0:
            if q < 0:
                if y0 <= shape[1]:
                    alpha= x0
                    beta= 0.
                    patmp = pa
                    x1,y1 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='min')
                    alpha= shape[2]
                    beta= y0 
                    patmp = 90.-pa
                    x2,y2 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='max')
                    #extrs = [[x0,0],[data.shape[2],y0]]                         
                else:
                    alpha= x0
                    beta= 0.
                    patmp = pa
                    x1,y1 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='min')
                    xtmp = (shape[1]-q)/m                            
                    alpha= xtmp
                    beta= shape[1]
                    patmp = pa
                    x2,y2 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='max')
                    #extrs = [[x0,0],[xtmp,data.shape[1]]]                         
            elif q >= 0 and q < shape[1]:
                if y0 <= shape[1]:
                    alpha= 0.
                    beta= q
                    patmp = pa
                    x1,y1 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='min')
                    alpha= shape[2]
                    beta= y0
                    patmp = 90.-pa
                    x2,y2 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='max')
                    #extrs = [[0,q],[data.shape[2],y0]]
                else:
                    alpha= 0.
                    beta= q
                    patmp = pa
                    x1,y1 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='min')
                    xtmp = (shape[1]-q)/m                            
                    alpha= xtmp
                    beta= shape[1]
                    patmp = pa
                    x2,y2 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='max')
                    #extrs = [[0,q],[x1,data.shape[1]]]
            elif q >=data.shape[1]:
                    print('ERROR: cutting outside of datacube, adjust PA or centre')
                    sys.exit(0)
        elif m < 0:
            if q > shape[1]:
                if y0 >= 0:
                    print('ciao')
                    alpha = (shape[1]-q)/m
                    beta= shape[1]
                    patmp=180.-pa
                    x1,y1 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='min')
                    print(x1,y1)
                    alpha= shape[2]
                    beta= y0
                    patmp =(180.-pa)
                    x2,y2 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='max')                              
                    print(x2,y2)

                    #extrs = [[x1,data.shape[1]],[data.shape[2],y0]]
                else:
                    x1 = (shape[1]-q)/m
                    alpha= x1
                    beta= shape[1]
                    patmp = 180.-pa
                    x1,y1 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='min')
                    alpha= x0
                    beta= 0
                    patmp = 180.-pa
                    x2,y2 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='max')
                    #extrs = [[x1,data.shape[1]],[x0,0]]
            elif q >=0 and q < shape[1]:
                if y0 >= 0:
                    alpha= 0.
                    beta= q
                    patmp = 90.-(180.-pa)
                    x1,y1 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='min')
                    alpha= shape[2]
                    beta= y0
                    patmp=90. - (180.-pa)                
                    x2,y2 = self.circPoint(alpha,beta,m,q,patmp,ww,keyAlpha='max')  
                    #extrs = [[0,q],[data.shape[2],y0]]
                else:
                    alpha= 0.
                    beta= q
                    patmp = 90.-(180.-pa)
                    x1,y1 = self.circPoint(alpha,beta,m,q,pa,ww,keyAlpha='min')
                    alpha = x0
                    beta = 0.
                    patmp = 180.-pa
                    x2,y2 = self.circPoint(alpha,beta,m,q,pa,ww,keyAlpha='max')
                    #extrs = [[0,q],[x0,0]]
            elif q <=0:
                    print('ERROR: cutting outside of datacube, adjust PA or centre')
                    sys.exit(0)
        
        extrs = [[np.rint(x1),np.rint(y1)],[np.rint(x2),np.rint(y2)]]
        return extrs,m

    def pvDiagramLines(self,cfg_par,pa):
        # -------------------------
        # Position-velocity diagram
        # -------------------------

        lineInfo = tP.openLineList(cfg_par)
 
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


            if cfg_par['cubePlay']['cube'] == 'vorLine': 
                inCubelet = cubeletsDir+str(lineNameStr)+'_measVor.fits'
            elif cfg_par['cubePlay']['cube'] == 'fitLine': 
                inCubelet = cubeletsDir+str(lineNameStr)+'_fit_'+modName+'.fits'   
            elif cfg_par['cubePlay']['cube'] == 'residuals':        
                inCubelet = cubeletsDir+str(lineNameStr)+'_res_'+modName+'.fits' 
            else:
                inCubelet = cubeletsDir+str(lineNameStr)+'_meas.fits'

            self.pvDiagram(cfg_par,inCubelet,pa)


            return

    def mosaicCubes(self, cfg_par):
        '''
        calls Montage() to mosaic a set of datacubes

        Parameters:
            - cfg_par: configuration file
        Uses:
            - cfg_par['cubePlay']['mosaic']['inDirectory']: directory of input cubes
        Creates: 
            - cfg_par['cubePlay']['mosaic']['inDirectory']+ mosaicProjections: directory with
                datacubes projected to common frame of reference
            - cfg_par['cubePlay']['mosaic']['inDirectory']+ mosaicProjections: directory with
                - projected datacube: mosaicCube.fits
                - tableInputs.tbl: table of input dataCubes
                - tableMosaic.tbl: table of projected dataCubes
                - hdrMosaic.hdr : header of common reference frame for mosaic
        Options:
            - cfg_par['cubePlay']['mosaic']['cleanup'] = True
                intermediate files are deleted and final mosaic is moved in base directory of cfg_par['cubePlay']['mosaic']['inDirectory']   
        '''

        pathCubes=cfg_par['cubePlay']['mosaic']['inDirectory']

        workDir = os.path.normpath(os.path.join(pathCubes, os.pardir))

        pathProj = workDir+'/mosaicProjections/'
        if not os.path.exists(pathProj):
            os.mkdir(pathProj)

        pathMosaic =  workDir+'/mosaicCube/'
        if not os.path.exists(pathMosaic):
            os.mkdir(pathMosaic)

        mImgtbl(pathCubes,pathMosaic+'tableMontage.tbl')
        mMakeHdr(pathMosaic+'tableMontage.tbl',pathMosaic+'hdrMosaic.hdr')
        
        with os.scandir(pathCubes) as it:
            for entry in it:
                if entry.name.endswith(".fits") and entry.is_file(): 
                    outProjection = str.split(entry.name,'.fits')[0]+'proj.fits'
                    mProjectCube(pathCubes+entry.name,pathProj+outProjection,pathMosaic+'hdrMosaic.hdr')

        mImgtbl(pathProj,pathMosaic+'tableMosaic.tbl')
        mAddCube(pathProj,pathMosaic+'tableMosaic.tbl',pathMosaic+'hdrMosaic.hdr',pathMosaic+'mosaicCube.fits')

        if cfg_par['cubePlay']['mosaic']['cleanup'] == True:
            ## Try to delete the file ##
            shutil.rmtree(pathProj)
            shutil.move(pathMosaic+'mosaicCube.fits',workDir+'/mosaicCube.fits')
            shutil.rmtree(pathMosaic)

        return

    def regridCube(self,cfg_par,header,inCube):

        
        outCube = str.split((inCube),'.')[0]
        outCube=outCube+'-reg.fits'
 
        
        #print('culo')
        #mGetHdr(tCube,'template.hrd')
        mProjectCube(inCube,outCube,header)

        return 0


    def regridCubeExact(self,inCube,pixSize):


        d = fits.getdata(inCube)
        h = fits.getheader(inCube)

        h2D , d2D = hP.cleanHead(inCube)

        h2D_or = h2D.copy()
   
        wcsh2D_or = WCS(h2D_or)

        h2D['NAXIS1'] = int(np.floor(-h2D['NAXIS1']*h2D['CDELT1']/(float(pixSize)/3600.)))
        h2D['NAXIS2'] = int(np.floor(h2D['NAXIS2']*h2D['CDELT2']/(float(pixSize)/3600.)))


        h2D['CRPIX1'] = int(np.floor(-h2D['CRPIX1']*h2D['CDELT1']/(float(pixSize)/3600.)))
        h2D['CRPIX2'] = int(np.floor(h2D['CRPIX2']*h2D['CDELT2']/(float(pixSize)/3600.)))

        h2D['CDELT1'] = float(-pixSize)/3600.
        h2D['CDELT2'] = float(pixSize)/3600.


        wcs2D = WCS(h2D)

        newCube = np.empty([d.shape[0],h2D['NAXIS2'],h2D['NAXIS1']])
        newIm = np.empty([h2D['NAXIS2'],h2D['NAXIS1']])
        print(d.shape)
        print(newIm.shape)
        for i in range(d.shape[0]):

            array, footprint = reproject_exact((d[i,:,:], h2D_or) ,
                                            wcs2D, shape_out=newIm.shape)
            newCube[i,:,:] = array



        newCubeName=inCube.replace('fits','_reg_pix_'+str(int(pixSize))+'asec.fits') 


        h3D = h2D.copy()

        h3D['NAXIS']= 3
        
        h3D['CRPIX3'] = h['CRPIX3'] 
        h3D['CDELT3'] = h['CDELT3'] 
        h3D['CRVAL3'] = h['CRVAL3'] 
        
        if 'CUNIT3' in h:
            h3D['CUNIT3'] = h['CUNIT3']
        if 'CTYPE3' in h:
            h3D['CTYPE3'] = h['CTYPE3']

        fits.writeto(newCubeName,newCube,h3D,overwrite=True)

        return newCubeName

    def rebinCube(self,templateFile,inputFile):

        tFile = fits.open(templateFile)
        iFile = fits.open(inputFile)

        tHead = tFile[0].header
        tData = tFile[0].data

        tVel = ((np.linspace(1,tData.shape[0],tData.shape[0])-tHead['CRPIX3'])*tHead['CDELT3']+tHead['CRVAL3'])/1e3
        tVel = tVel[::-1]
        iHead = iFile[0].header
        iData = iFile[0].data

        iVel = ((np.linspace(1,iData.shape[0],iData.shape[0])-iHead['CRPIX3'])*iHead['CDELT3']+iHead['CRVAL3'])/1e3
        iVel = iVel[::-1]
        data = np.zeros([tData.shape[0],tData.shape[1],tData.shape[2]])

        if iData.shape[0] != tData.shape[0]:
            for i in range(0,tData.shape[0]-1):
                index = (tVel[i] <= iVel) & (iVel < tVel[i+1])
                print(index)
                print(tVel[i])
                data[i,:,:] = np.sum(iData[index,:,:])
        else:
            data = np.copy(iData)
        iHead['CRVAL3'] = tHead['CRVAL3']
        iHead['CDELT3'] = tHead['CDELT3']
        iHead['CRPIX3'] = tHead['CRPIX3']
        iHead['CTYPE3'] = 'm/s'


        rebinFileName = str.split((inputFile),'.')[0]
        rebinFileName=rebinFileName+'-rebin.fits'
        print(rebinFileName)
        fits.writeto(rebinFileName,data,iHead,overwrite=True)



        return 0

    def medianSubtract(self,cfg_par):


        inCubelet=cfg_par['cubePlay']['medSub']['inCube']
        f = fits.open(inCubelet)
        data = f[0].data
        headerCubelets = f[0].header

        for i in range(0,data.shape[2]):
            for j in range(0,data.shape[1]):
                median1 = np.nanmedian(data[cfg_par['cubePlay']['medSub']['interval'][0][0]:cfg_par['cubePlay']['medSub']['interval'][0][1],j,i])
                median2 = np.nanmedian(data[cfg_par['cubePlay']['medSub']['interval'][1][0]:cfg_par['cubePlay']['medSub']['interval'][1][1],j,i])

                median = (median1+median2)/2.
                data[:,j,i] =np.subtract(data[:,j,i],median)
        
        medSubFile = str.split((inCubelet),'.')[0]
        medSubFile=medSubFile+'-medSub.fits'
        fits.writeto(medSubFile,data,headerCubelets,overwrite=True)

        return 0



