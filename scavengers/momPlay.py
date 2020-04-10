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
bpt = bptPlot.BPTplot()
mPl = momPlot.MOMplot()

class momplay:
    '''Modules to create moment maps, residual maps, line ratios maps
    - makeMoments
        load line list, datacube and create loop for moments module
    - makeSigmaCentroidMap
        load line list, datacube and create loop for momSigmaCentroid module
    - makeMomPlots
        load line list, call momPlot for each line for mom0, mom1 maps from gaussian components
    - momSigmaCentroid
        mom0, mom1 from centroid and sigma of fitted line
    - moments
        mom0, mom1, mom2 of the fitted gaussian components of the line
    - resCube
        cube of residuals of fit
    - resLines
        for each line residuals are computed as the standard deviation of line-fit 
        within vrange and as the sum of the absolute value of line-fit 
    - makeLineRatioMaps
        load line list, and create loop for momLineRatio
    - momLineRatio
        maps of the line ratios (OIII/Hbeta, NII/Halpha, SII/Halpha)
    - momCDist
        map of the eta-parameter (distance from Kauffmann and Kewley SF curves)        
    '''

    def makeMoments(self,cfg_par):

        workDir = cfg_par['general']['cubeDir']

        f = fits.open(cfg_par['general']['dataCubeName'])
        dd = f[0].header

        lineInfo = tP.openLineList(cfg_par)
        for ii in range(0,len(lineInfo['ID'])):
        #for ii in range(0,1):

            lineNameStr = str(lineInfo['Name'][ii])

            if '[' in lineNameStr:
                lineName = lineNameStr.replace("[", "")
                lineName = lineName.replace("]", "")
                lineName = lineName+str(int(lineInfo['Wave'][ii]))
            else:
                lineName = lineNameStr+str(int(lineInfo['Wave'][ii]))

            lineNameStr=lineNameStr+str(int(lineInfo['Wave'][ii]))
            lineThresh = float(lineInfo['SNThresh'][ii])
            cenRange = float(lineInfo['cenRange'][ii])

            print('\n\t         +++\t\t    '+lineName+'\t\t +++')

            
            if ii==0:
                doBinMap=True
            else:
                doBinMap=True
            
            self.moments(cfg_par,lineName,lineNameStr,dd,cfg_par['general']['outTableName'],lineThresh,doBinMap,cenRange)

        return

    def makeSigmaCentroidMaps(self,cfg_par):

        workDir = cfg_par['general']['cubeDir']

        f = fits.open(cfg_par['general']['dataCubeName'])
        dd = f[0].header

        lineInfo = tP.openLineList(cfg_par)
        for ii in range(0,len(lineInfo['ID'])):
        #for ii in range(0,1):

            lineNameStr = str(lineInfo['Name'][ii])

            if '[' in lineNameStr:
                lineName = lineNameStr.replace("[", "")
                lineName = lineName.replace("]", "")
                lineName = lineName+str(int(lineInfo['Wave'][ii]))
            else:
                lineName = lineNameStr+str(int(lineInfo['Wave'][ii]))

            lineNameStr=lineNameStr+str(int(lineInfo['Wave'][ii]))
            lineThresh = float(lineInfo['SNThresh'][ii])
            cenRange = float(lineInfo['cenRange'][ii])

            print('\t         +++\t\t'+lineName+'\t\t+++')
                        
            self.momSigmaCentroid(cfg_par,lineName,lineNameStr,dd,lineThresh,cenRange)

        return

    def makeMomPlots(self,cfg_par):

        workDir = cfg_par['general']['cubeDir']
        modName = cfg_par['gFit']['modName']
        momModDir = cfg_par['general']['momDir']+modName+'/'

        lineInfo = tP.openLineList(cfg_par)
        for ii in range(0,len(lineInfo['ID'])):
        #for ii in range(0,1):

            lineNameStr = str(lineInfo['Name'][ii])

            if '[' in lineName:
                lineName = lineNameStr.replace("[", "")
                lineName = lineName.replace("]", "")

            lineName = lineName+str(int(lineInfo['Wave'][ii]))
            lineThresh = float(lineInfo['SNThresh'][ii])
            cenRange = float(lineInfo['cenRange'][ii])

            print('\n\t *********** --- Plot Moms: '+lineName+' --- ***********\n')
            
            mom0Name = momModDir+'mom0_g1-'+lineName+'.fits'          
            mom1Name = momModDir+'mom1_g1-'+lineName+'.fits'
            
            mPl.mom0Plot(cfg_par, mom0Name,lineName,lineNameStr,lineThresh)

            mPl.mom1Plot(cfg_par, mom1Name,lineName,lineThresh,lineNameStr, 'moments',vRange=[-cenRange,cenRange])

        return

    def momSigmaCentroid(self,cfg_par,lineName,lineNameStr,header,lineThresh,cenRange):

        modName = cfg_par['gFit']['modName']
        momModDir = cfg_par['general']['momDir']+modName+'/'

        if not os.path.exists(momModDir):
            os.mkdir(momModDir)

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

        momSigmaHead = header.copy()
        momCentroidHead = header.copy()
        momW80Head = header.copy()

        hdul = fits.open(cfg_par['general']['outTableName'])
        lines = hdul['Ancels'+cfg_par['gFit']['modName']].data
        
        hduGen = fits.open(cfg_par['general']['outVorLineTableName'])
        tabGen = hduGen[1].data

        momSigma = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        momCentroid = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        momW80 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        
        for i in range(0,len(lines['BIN_ID'])):

            match_bin = np.where(tabGen['BIN_ID']==lines['BIN_ID'][i])[0]

            for index in match_bin:
                
                momW80[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['w80_'+lineName][i]
                momSigma[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['sigma_'+lineName][i]
                momCentroid[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['centroid_'+lineName][i]

        del momSigmaHead['CRDER3']
        del momCentroidHead['CRDER3']
        del momW80Head['CRDER3']

        momSigmaHead['WCSAXES'] = 2
        momSigmaHead['SPECSYS'] = 'topocent'
        momSigmaHead['BUNIT'] = 'km/s'

        fits.writeto(momModDir+'momSigma-'+lineName+'.fits',momSigma,momSigmaHead,overwrite=True)

        mPl.mom2Plot(cfg_par, momModDir+'momSigma-'+lineName+'.fits',lineName,lineThresh,lineNameStr,'ancillary')

        momCentroidHead['WCSAXES'] = 2
        momCentroidHead['SPECSYS'] = 'topocent'
        momCentroidHead['BUNIT'] = 'km/s'
        fits.writeto(momModDir+'momCentroid-'+lineName+'.fits',momCentroid,momCentroidHead,overwrite=True)
        mPl.mom1Plot(cfg_par, momModDir+'momCentroid-'+lineName+'.fits',lineName,lineThresh,
            lineNameStr,'ancillary',vRange=[-cenRange,cenRange],modName=cfg_par['gFit']['modName'])

        momW80Head['WCSAXES'] = 2
        momW80Head['SPECSYS'] = 'topocent'
        momW80Head['BUNIT'] = 'km/s'
        fits.writeto(momModDir+'momW80-'+lineName+'.fits',momW80,momW80Head,overwrite=True)
        mPl.mom2Plot(cfg_par, momModDir+'momW80-'+lineName+'.fits',lineName,lineThresh,lineNameStr,'ancillary')

        return



    def moments(self,cfg_par,lineName,lineNameStr,header,outTableName,lineThresh,doBinMap,cenRange):

        modName = cfg_par['gFit']['modName']
        momModDir = cfg_par['general']['momDir']+modName+'/'

        if not os.path.exists(momModDir):
            os.mkdir(momModDir)

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

        mom0Head = header.copy()
        mom1Head = header.copy()
        mom2Head = header.copy()
        binHead  = header.copy()


        hdul = fits.open(cfg_par['general']['outTableName'])

        lines = hdul['LineRes_'+cfg_par['gFit']['modName']].data
        

        hduGen = fits.open(cfg_par['general']['outVorLineTableName'])
        tabGen = hduGen[1].data

        ampSpax = np.empty(len(tabGen['BIN_ID']))


        mom0G1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        mom1G1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        mom2G1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        heightG1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan

        if doBinMap==True:
            binMap = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        
        if modName != 'g1':
            mom0Tot = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            mom0G2 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            mom1G2 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            mom2G2 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            if modName == 'g3':
                mom0G3 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
                mom1G3 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
                mom2G3 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        
        for i in range(0,len(lines['BIN_ID'])):
            #if lines['BIN_ID'][i]< 0:
            #    continue
            #else:

            match_bin = np.where(tabGen['BIN_ID']==lines['BIN_ID'][i])[0]
            for index in match_bin:
                
                if modName=='g1':
                    thresHold = lines['g1_Amp_Hb4861'][i]/tabGen['NSPAX'][index]
                    ampSpax[index] = lines['g1_Amp_'+lineName][i]/tabGen['NSPAX'][index]                   
                elif modName=='g2':
                    thresHold = (lines['g1_Amp_Hb4861'][i]+lines['g2_Amp_Hb4861'][i])/tabGen['NSPAX'][index]
                    ampSpax[index] = (lines['g1_Amp_'+lineName][i]+lines['g2_Amp_Hb4861'][i])/tabGen['NSPAX'][index]                   
                elif modName=='g3':
                    thresHold = (lines['g1_Amp_Hb4861'][i]+lines['g2_Amp_Hb4861'][i]+lines['g3_Amp_Hb4861'][i])/tabGen['NSPAX'][index]
                    ampSpax[index] = (lines['g1_Amp_'+lineName][i]+lines['g2_Amp_Hb4861'][i]+lines['g3_Amp_Hb4861'][i])/tabGen['NSPAX'][index]  

                #thresHold = lines['g1_Height_'+lineName][i]/0.3989423*lines['g1_Sigma_'+lineName][i]/noise[0,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]
                #print(lines['g1_Height_'+lineName][i]/0.3989423*lines['g1_Sigma_'+lineName][i],lines['g1_Sigma_'+lineName][i],lines['g1_Height_'+lineName][i])
                #print(thresHold,lineThresh)
                if thresHold >= lineThresh:
                    mom0G1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = ampSpax[index]
#                        mom0G1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g1_Height_'+lineName][i]/tabGen['NSPAX'][index]
                    mom1G1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g1_Centre_'+lineName][i]
                    mom2G1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g1_SigIntr_'+lineName][i]
                    heightG1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g1_Height_'+lineName][i]                

                    if doBinMap==True:
                        binMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['BIN_ID'][i]
                    
                    if modName != 'g1':
                        mom0G2[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = ampSpax[index]
                        mom1G2[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g2_Centre_'+lineName][i]
                        mom2G2[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g2_SigIntr_'+lineName][i]
                    
                        if modName == 'g3':
                            mom0G3[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = ampSpax[index]
                            mom1G3[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g3_Centre_'+lineName][i]
                            mom2G3[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g3_SigIntr_'+lineName][i]
                        
                        mom0Tot[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = ampSpax[i]


                    #else#:
                    #    print(mom1G1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])])
                    #    print(int(tabGen['PixY'][index]),int(tabGen['PixX'][index]))
        
        if doBinMap==True:

            binHead['SPECSYS'] = 'topocent'
            binHead['BUNIT'] = 'Flux'
            fits.writeto(momModDir+'binMapMom0_'+lineName+'.fits',binMap, binHead,overwrite=True)

        del mom0Head['CRDER3']
        del mom1Head['CRDER3']
        del mom2Head['CRDER3']

        mom0Head['WCSAXES'] = 2

        mom0Head['SPECSYS'] = 'topocent'
        mom0Head['BUNIT'] = 'Jy/beam.km/s'
        fits.writeto(momModDir+'mom0_g1-'+lineName+'.fits',mom0G1,mom0Head,overwrite=True)
        fits.writeto(momModDir+'height_g1-'+lineName+'.fits',heightG1,mom0Head,overwrite=True)
        
        mPl.mom0Plot(cfg_par, momModDir+'mom0_g1-'+lineName+'.fits',lineName,lineNameStr,lineThresh)

        mom1Head['WCSAXES'] = 2
        mom1Head['SPECSYS'] = 'topocent'
        mom1Head['BUNIT'] = 'km/s'
        fits.writeto(momModDir+'mom1_g1-'+lineName+'.fits',mom1G1,mom1Head,overwrite=True)
        mPl.mom1Plot(cfg_par, momModDir+'mom1_g1-'+lineName+'.fits',lineName,
            lineThresh, lineNameStr,'moments', vRange=[-cenRange,cenRange],modName='g1')

        mom2Head['WCSAXES'] = 2
        mom2Head['SPECSYS'] = 'topocent'
        mom2Head['BUNIT'] = 'km/s'
        fits.writeto(momModDir+'mom2_g1-'+lineName+'.fits',mom2G1,mom2Head,overwrite=True)
        
        if modName != 'g1':
            fits.writeto(momModDir+'mom0_g2-'+lineName+'.fits',mom0G2,mom0Head,overwrite=True)
            fits.writeto(momModDir+'mom1_g2-'+lineName+'.fits',mom1G2,mom1Head,overwrite=True)
            mPl.mom0Plot(cfg_par, momModDir+'mom0_g2-'+lineName+'.fits',lineName,lineNameStr,lineThresh)
            mPl.mom1Plot(cfg_par, momModDir+'mom1_g2-'+lineName+'.fits',lineName,
                lineNameStr,lineThresh,'moments',vRange=[-cenRange,cenRange],
                modName='g2')

            if modName == 'g2':
                fits.writeto(momModDir+'mom0_tot-'+lineName+'.fits',mom0G1+mom0G2,mom0Head,overwrite=True)
                
                mPl.mom0Plot(cfg_par, momModDir+'mom0_tot-'+lineName+'.fits',lineName,lineNameStr,lineThresh)


            if modName == 'g3':
                fits.writeto(momModDir+'mom0_g3-'+lineName+'.fits',mom0G3,mom0Head,overwrite=True)
                fits.writeto(momModDir+'mom1_g3-'+lineName+'.fits',mom1G3,mom1Head,overwrite=True)
                fits.writeto(momModDir+'mom2_g3-'+lineName+'.fits',mom2G3,mom2Head,overwrite=True)
                fits.writeto(momModDir+'mom0_tot-'+lineName+'.fits',mom0G1+mom0G2+mom0G3,mom0Head,overwrite=True)


        t=Table(tabGen)
        if modName+'-AmpSpax_'+lineName not in tabGen.dtype.names: 
            t.add_column(Column(ampSpax,name=modName+'-AmpSpax_'+lineName))
        else:
            t.replace_column(modName+'-AmpSpax_'+lineName,Column(ampSpax,name=modName+'-AmpSpax_'+lineName))        
        
        #try:
        #    tt = Table(hduGen['VORBININFO'].data)
        hduGen['VORBININFO'] = fits.BinTableHDU(t.as_array(),name='VORBININFO')

        #except KeyError as e:
        #    tt=fits.BinTableHDU(t.as_array(),name='VORBININFO')   
        #    hdul.append(tt)          
        
        hduGen.writeto(cfg_par['general']['outVorLineTableName'],overwrite=True)


        return


    def resCube(self,cfg_par):

        key = 'general'
        cubeDir = cfg_par['general']['cubeDir']
        workDir = cfg_par['general']['workdir']


        modName = cfg_par['gFit']['modName']

        resModDir = cfg_par['general']['resDir']+modName+'/'


        # if not os.path.exists(resModDir):
        #     os.mkdir(momModDir)

        f = fits.open(cfg_par['general']['dataCubeName'])
        dd = f[0].data
        resHead = f[0].header
        
        hdul = fits.open(cfg_par['general']['outTableName'])
        lines = hdul['LineRes_'+cfg_par['gFit']['modName']].data

        hduGen = fits.open(cfg_par['general']['outVorLineTableName'])
        tabGen = hduGen[1].data

        resG1 = np.zeros([resHead['NAXIS3'],resHead['NAXIS2'],resHead['NAXIS1']])*np.nan
        fitCube = np.zeros([resHead['NAXIS3'],resHead['NAXIS2'],resHead['NAXIS1']])*np.nan
    

        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec = tP.openVorLineOutput(cfg_par,cfg_par['general']['outVorLineTableName'],
            cfg_par['general']['outVorSpectra'])

        #hdul = fits.open(cfg_par['general']['outTableName'])
        #tabGen = hdul['BinInfo'].data

        lambdaMin = np.log(cfg_par['gFit']['lambdaMin'])
        lambdaMax = np.log(cfg_par['gFit']['lambdaMax'])
        idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
        idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0])
        
        if modName != 'g1':
            resTot = np.zeros([resHead['NAXIS3'],resHead['NAXIS2'],resHead['NAXIS1']])*np.nan
            resG2 = np.zeros([resHead['NAXIS3'],resHead['NAXIS2'],resHead['NAXIS1']])*np.nan
            if modName == 'g3':
                res0G3 = np.zeros([resHead['NAXIS3'],resHead['NAXIS2'],resHead['NAXIS1']])*np.nan

        for i in range(0,len(lines['BIN_ID'])):

            match_bin = np.where(tabGen['BIN_ID']==lines['BIN_ID'][i])[0]

            result = load_modelresult(cfg_par[key]['modNameDir']+str(lines['BIN_ID'][i])+'_'+cfg_par['gFit']['modName']+'.sav')

            for index in match_bin:
                yy = dd[idxMin:idxMax,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]
                fit = result.best_fit
                residuals = result.best_fit-yy
                resG1[idxMin:idxMax,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = residuals
                fitCube[idxMin:idxMax,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = fit

        resHead['SPECSYS'] = 'topocent'
        resHead['BUNIT'] = 'Flux'
        fits.writeto(resModDir+'resAllLines_'+modName+'.fits',resG1,resHead,overwrite=True)

        return


    def resLines(self,cfg_par): 
        '''
        Computes for each the residuals of the fit. Within which velocity range? 
            At the moment is within 6*sigmaG1 and 3*sigmaG2

        Parameters:
            cfg_par: parameter file
                gFit_modName: specifies # of gaussian components used for the fit

        Uses:
            - voroni binned line subtracted datacube
            - table of voronoi binned datacube and spectra
            -
        Returns (located in /moments/modName/):
            - resAbs_linename:  residuals computed as sum of the absolute value of line-fit
                                within a velocity range given by 6*sigmag1 weighted on the fitted amplitude of the line
            - resSTD_linename:  residuals computed as the standard deviation of line-fit
                                within a velocity range given by 6*sigmag1 weighted on the fitted amplitude of the line
        '''


        key = 'general'
        workDir = cfg_par[key]['workdir'] 
        cubeDir = cfg_par[key]['cubeDir'] 
        modName = cfg_par['gFit']['modName']
        resModDir = cfg_par['general']['resDir']+modName+'/'
        noiseDir = cfg_par['general']['noiseDir']
        resName = resModDir+'resAllLines_'+modName+'.fits'
        fitCubeName = cubeDir+'fitCube_'+modName+'.fits'

        if not os.path.exists(resName):
            self.resCube(cfg_par)
        else:
            pass
        
        f = fits.open(resName)
        resCube = f[0].data
        resHead = f[0].header
        if 'CUNIT3' in resHead:
            del resHead['CUNIT3']
        if 'CTYPE3' in resHead:
            del resHead['CTYPE3']
        if 'CDELT3' in resHead:
            del resHead['CDELT3']
        if 'CRVAL3' in resHead:  
            del resHead['CRVAL3']
        if 'CRPIX3' in resHead:
            del resHead['CRPIX3'] 
        if 'NAXIS3' in resHead:
            del resHead['NAXIS3']
        if 'CRDER3'in resHead:
            del resHead['CRDER3']
        

        f = fits.open(cfg_par['general']['outVorLines'])
        dd = f[0].data
        
        #to load Voronoi Bin noise : noiseBin
        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec = tP.openVorLineOutput(cfg_par,cfg_par['general']['outVorLineTableName'],
            cfg_par['general']['outVorSpectra'])

        f = fits.open(cfg_par['general']['dataCubeName'])
        dd = f[0].data
        header = f[0].header

        hdul = fits.open(cfg_par['general']['outTableName'])
        lines = hdul['LineRes_'+cfg_par['gFit']['modName']].data

        resNameList=['BIN_ID']
        frmList=['i4']
        tot = lines['BIN_ID']

        hduGen = fits.open(cfg_par['general']['outVorLineTableName'])
        tabGen = hduGen[1].data
        
        lineInfo = tP.openLineList(cfg_par)

        tableSpec = workDir+cfg_par[key]['tableSpecName']
        tab = fits.open(tableSpec)
        dataSpec = tab[1].data
        specExp = tab[2].data
        wave = [item for t in specExp for item in t] 
        
        noiseMapName =noiseDir+'noiseMap.fits'
        noiseMap = np.empty([resHead['NAXIS2'],resHead['NAXIS1']])*np.nan

        for ii in range(0,len(lineInfo['ID'])):
            
            stdArr = np.empty(len(lines['BIN_ID']))
            noiseArr = np.empty(len(lines['BIN_ID']))
            SNValues = np.empty(len(lines['BIN_ID']))
            SNStdValues = np.empty(len(lines['BIN_ID']))

            lineName = str(lineInfo['Name'][ii])
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")

            lineName = lineName+str(int(lineInfo['Wave'][ii]))            
            lineThresh = float(lineInfo['SNThresh'][ii])
            print('\n\t         +++\t\t    '+lineName+'\t\t +++')
          
            resG1Abs = np.empty([resHead['NAXIS2'],resHead['NAXIS1']])*np.nan
            resG1Std = np.empty([resHead['NAXIS2'],resHead['NAXIS1']])*np.nan
            noiseLine = np.empty([resHead['NAXIS2'],resHead['NAXIS1']])*np.nan
            SNLineMap = np.empty([resHead['NAXIS2'],resHead['NAXIS1']])*np.nan
            SNStdLineMap = np.empty([resHead['NAXIS2'],resHead['NAXIS1']])*np.nan


            resNameOutAbs =resModDir+'resAbs_'+lineName+'.fits'
            resNameOutStd =resModDir+'resStd_'+lineName+'.fits'
            noiseNameLine =noiseDir+'noise_'+lineName+'.fits'
            SNMapName =noiseDir+'SN_'+lineName+'.fits'
            SNStdMapName =resModDir+'SNRes_'+lineName+'.fits'
   
            for i in range(0,len(lines['BIN_ID'])):

            
                #lineHeigth = np.max(y[indexMin:indexMax])    
        
                amp = lines['g1_Amp_'+lineName][i]
                
                cenKmsG1 = lines['g1_Centre_'+lineName][i]
                sigKmsG1 = lines['g1_SigIntr_'+lineName][i]
                
                if sigKmsG1 >=2.e3:
                    sigKmsG1=2.e3

                cenG1 = np.log(cvP.vRadLambda(cenKmsG1,lineInfo['Wave'][ii]))
                leftG1 = np.log(cvP.vRadLambda(cenKmsG1-6.*sigKmsG1,lineInfo['Wave'][ii]))
                rightG1 = np.log(cvP.vRadLambda(cenKmsG1+6.*sigKmsG1,lineInfo['Wave'][ii]))

                idxLeft = int(np.where(abs(wave-leftG1)==abs(wave-leftG1).min())[0])
                idxRight = int(np.where(abs(wave-rightG1)==abs(wave-rightG1).min())[0])

                #define interval where to measure maximum of real line from centroid of 1G-fit
                peakLeft = np.log(cvP.vRadLambda(cenKmsG1-140.,lineInfo['Wave'][ii]))
                peakRight = np.log(cvP.vRadLambda(cenKmsG1+140.,lineInfo['Wave'][ii]))

                idxPeakLeft = int(np.where(abs(wave-peakLeft)==abs(wave-peakLeft).min())[0])
                idxPeakRight = int(np.where(abs(wave-peakRight)==abs(wave-peakRight).min())[0])
 
                if cfg_par['residuals']['computeNoise']==True:
            
                    leftNoise = np.log(lineInfo['Wave'][ii]-60.)
                    leftleftNoise = np.log(lineInfo['Wave'][ii]-80.)
                    rightNoise = np.log(lineInfo['Wave'][ii]+60.)
                    rightrightNoise = np.log(lineInfo['Wave'][ii]+80.)

                    idxLeftLeftNoise = int(np.where(abs(wave-leftleftNoise)==abs(wave-leftleftNoise).min())[0])
                    idxLeftNoise = int(np.where(abs(wave-leftNoise)==abs(wave-leftNoise).min())[0])

                    idxRightRightNoise = int(np.where(abs(wave-rightrightNoise)==abs(wave-rightrightNoise).min())[0])
                    idxRightNoise = int(np.where(abs(wave-rightNoise)==abs(wave-rightNoise).min())[0])
                    idxTable = int(np.where(tabGen['BIN_ID'] == int(lines['BIN_ID'][i]))[0][0])
                    y = dd[:,int(tabGen['PixY'][idxTable]),int(tabGen['PixX'][idxTable])]                

                if modName == 'g2':
                    amp = lines['g1_Amp_'+lineName][i]+lines['g2_Amp_'+lineName][i]
                    cenKmsG2 = lines['g2_Centre_'+lineName][i]
                    sigKmsG2 = lines['g2_SigMeas_'+lineName][i]
            
                    cenG2 = np.log(cvP.vRadLambda(cenKmsG2,lineInfo['Wave'][ii]))
                    leftG2 = np.log(cvP.vRadLambda(cenKmsG2-3.*sigKmsG2,lineInfo['Wave'][ii]))
                    rightG2 = np.log(cvP.vRadLambda(cenKmsG2+3.*sigKmsG2,lineInfo['Wave'][ii]))
                    
                    idxLeftG2 = int(np.where(abs(wave-leftG2)==abs(wave-leftG2).min())[0])
                    idxRightG2 = int(np.where(abs(wave-rightG2)==abs(wave-rightG2).min())[0])
                    
                    idxLeft = np.min([idxLeft,idxLeftG2])
                    idxRight = np.max([idxRight,idxRightG2])

                    if modName =='g3':

                        cenKmsG3 = lines['g3_Centre_'+lineName][i]
                        sigKmsG3 = lines['g3_SigMeas_'+lineName][i]
                
                        cenG2 = np.log(cvP.vRadLambda(cenKmsG1,lineInfo['Wave'][ii]))
                        leftG2 = np.log(cvP.vRadLambda(cenKmsG1-3.*sigKmsG3,lineInfo['Wave'][ii]))
                        rightG2 = np.log(cvP.vRadLambda(cenKmsG1+3.*sigKmsG3,lineInfo['Wave'][ii]))
                        
                        idxLeftG3 = int(np.where(abs(wave-leftG3)==abs(wave-leftG3).min())[0])
                        idxRightG3 = int(np.where(abs(wave-rightG3)==abs(wave-rightG3).min())[0])
                        
                        idxLeft = np.min([idxLeft,idxLeftG3])
                        idxRight = np.max([idxRight,idxRightG3])

                if ii==0 and cfg_par['residuals']['computeNoise']==True:
                    noiseValue = noiseBin[lines['BIN_ID'][i]][idxLeft]*amp

                match_bin = np.where(tabGen['BIN_ID']==lines['BIN_ID'][i])[0]


                #result = load_modelresult(cfg_par[key]['modNameDir']+str(lines['BIN_ID'][i])+'_'+cfg_par['gFit']['modName']+'.sav')
                if idxRight-idxLeft <2.:
                    idxLeft-=4
                    idxRight+=4
                for index in match_bin:

                    # if modName=='g1':
                    #     thresHold = lines['g1_Amp_Hb4861'][i]/tabGen['NSPAX'][index]
                    # elif modName=='g2':
                    #     thresHold = (lines['g1_Amp_Hb4861'][i]+lines['g2_Amp_Hb4861'][i])/tabGen['NSPAX'][index]
                    # elif modName=='g3':
                    #     thresHold = (lines['g1_Amp_Hb4861'][i]+lines['g2_Amp_Hb4861'][i]+lines['g3_Amp_Hb4861'][i])/tabGen['NSPAX'][index]
                    
                    # if thresHold >= lineThresh:
                    absValue = np.multiply(np.nansum(np.abs(resCube[idxLeft:idxRight,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]),axis=0),amp)
                    stdValue = np.multiply(np.nanstd(resCube[idxLeft:idxRight,int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]),amp)

                    resG1Abs[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = absValue
                    resG1Std[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = stdValue
                    
                    if cfg_par['residuals']['computeNoise']==True:
                        noise = np.nanstd(np.concatenate([y[idxLeftLeftNoise:idxLeftNoise],y[idxRightNoise:idxRightRightNoise]]))
                        linePeak = np.max(y[idxPeakLeft:idxPeakRight])
                        sn = np.divide(linePeak,noise)
                        snStd = np.divide(linePeak,stdValue)

                        noiseLine[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = noise

                        SNLineMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = sn
                        SNStdLineMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = snStd


                        if ii==0: 
                            noiseMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = noiseValue

                stdArr[i] = stdValue
                noiseArr[i] = noise
                SNValues[i] = sn
                SNStdValues[i] = snStd
 
            tot = np.column_stack((tot,stdArr))
            resNameList.append('res_'+lineName)
            frmList.append('f8')


            resHead['WCSAXES'] = 2
                      
            fits.writeto(resNameOutAbs,resG1Abs,resHead,overwrite=True)
            fits.writeto(resNameOutStd,resG1Std,resHead,overwrite=True)

            if cfg_par['residuals']['computeNoise']==True:
                fits.writeto(noiseNameLine,noiseLine,resHead,overwrite=True)
                fits.writeto(SNMapName,SNLineMap,resHead,overwrite=True)
                fits.writeto(SNStdMapName,SNStdLineMap,resHead,overwrite=True)

                tot = np.column_stack((tot,noiseArr))
                resNameList.append('noise_'+lineName)
                frmList.append('f8')

                tot = np.column_stack((tot,SNValues))
                resNameList.append('SN_'+lineName)
                frmList.append('f8')
                
                tot = np.column_stack((tot,SNStdValues))
                resNameList.append('snRes_'+lineName)
                frmList.append('f8')

                if ii==0:
                    fits.writeto(noiseMapName,noiseMap,resHead,overwrite=True)

        t = Table(tot, names=(resNameList))
        hdul.append(fits.BinTableHDU(t.as_array(), name='Residuals_'+modName))

        try:
            tt = Table(hdul['Residuals_'+modName].data)
            hdul['Residuals_'+modName] = fits.BinTableHDU(tt.as_array(),name='Residuals_'+modName)

        except KeyError as e:
            tt=fits.BinTableHDU.from_columns(t.as_array(),name='Residuals_'+modName)   
            hdul.append(tt)          
        
        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)



        return 0

    def makeLineRatioMaps(self,cfg_par):

        workDir = cfg_par['general']['cubeDir']

        f = fits.open(cfg_par['general']['dataCubeName'])
        dd = f[0].header

        #lineInfo = self.openLineList()

        #for ii in range(0,len(lineInfo['ID'])):
        #    lineName = str(lineInfo['Name'][ii])
        #    if '[' in lineName:
        #        lineName = lineName.replace("[", "")
        #        lineName = lineName.replace("]", "")

            #lineName = lineName+str(int(lineInfo['Wave'][ii]))

        self.momLineRatio(cfg_par,dd,cfg_par['general']['outTableName'])

        return


    def momLineRatio(self,cfg_par,header,outTableName):

        modName = cfg_par['gFit']['modName']
        bptDir = cfg_par['general']['bptDir']+'/'
        momModDir = cfg_par['general']['momDir']+modName+'/'

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
        if 'WCSAXES' in header:
            del header['WCSAXES']
        if 'CRDER3' in header:
            del header['CRDER3']

        lineMapHead = header.copy()

        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data

        hduGen = fits.open(cfg_par['general']['outVorLineTableName'])
        tabGen = hduGen[1].data

        hbetaMap = fits.open(momModDir+'mom0_'+modName+'-Hb4861.fits')
        hbetaData = hbetaMap[0].data

        numCols = len(lineBPT.dtype.names)

        if modName == 'g2':
            numCols = int((numCols-1)/3)
            numCols +=1 
        if modName == 'g3':
            numCols = int((numCols-1)/4)
            numCols +=1 

        for i in range(1,numCols):
            lineMapG1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        
            if modName != 'g1':
                lineMapToT = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
                lineMapG2 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            
            if modName == 'g3':
            
                lineMapG3 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            
            for j in range(0,len(lineBPT['BIN_ID'])):
                match_bin = np.where(tabGen['BIN_ID']==lineBPT['BIN_ID'][j])[0]

                for index in match_bin:
                    if ~np.isnan(hbetaData[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]):
                        lineMapG1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j][i]

                    if modName != 'g1': 

                        lineMapToT[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j][i+numCols*2-2]
                        lineMapG2[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j][i+numCols-1]

                        if modName == 'g3':
                            lineMapG3[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j][i+numCols+2] #TOREVIEW!!!!

            lineMapHead['BUNIT'] = 'Flux'
            outBPT = bptDir+'BPT-'+str(lineBPT.dtype.names[i])+'.fits'
            fits.writeto(bptDir+'BPT-'+str(lineBPT.dtype.names[i])+'.fits',lineMapG1,lineMapHead,overwrite=True)
            
            if modName != 'g1':
                outBPTg2 = bptDir+'BPT-'+str(lineBPT.dtype.names[i+numCols-1])+'.fits'
                outBPTtot = bptDir+'BPT-'+str(lineBPT.dtype.names[i+numCols*2-2])+'.fits'

                fits.writeto(outBPTg2,lineMapG2,lineMapHead,overwrite=True)
                fits.writeto(outBPTtot,lineMapToT,lineMapHead,overwrite=True)

                if modName == 'g3':
                    outBPTg3 = bptDir+'BPT-'+str(lineBPT.dtype.names[i+numCols+2])+'.fits'
                    fits.writeto(bptDir+'BPT-'+str(lineBPT.dtype.names[i+numCols+2])+'.fits',lineMapG3,lineMapHead,overwrite=True)

            if cfg_par['lineRatios']['bptMap'] == True:
                bpt.bptIM(cfg_par,outBPT)
                if modName != 'g1':
                    bpt.bptIM(cfg_par,outBPTg2)
                    bpt.bptIM(cfg_par,outBPTtot)
                elif modName=='g3':
                    bpt.bptIM(cfg_par,outBPTg3)

        return

    def momCDist(self,cfg_par):

        f = fits.open(cfg_par['general']['dataCubeName'])
        header = f[0].header
        f.close()
        modName = cfg_par['gFit']['modName']
        bptDir = cfg_par['general']['bptDir']+'/'
        momModDir = cfg_par['general']['momDir']+modName+'/'

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
        if 'WCSAXES' in header:
            del header['WCSAXES']
        if 'CRDER3' in header:
            del header['CRDER3']

        lineMapHead = header.copy()

        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data
        hbetaMap = fits.open(momModDir+'mom0_'+modName+'-Hb4861.fits')
        hbetaData = hbetaMap[0].data

        hduGen = fits.open(cfg_par['general']['outVorLineTableName'])
        tabGen = hduGen[1].data

        numCols = len(lineBPT.dtype.names)

        if modName == 'g2':
            numCols = int((numCols-1)/3)
            numCols +=1 
        if modName == 'g3':
            numCols = int((numCols-1)/4)
            numCols +=1 

        lineMapG1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
    
        if modName != 'g1':
            lineMapToT = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            lineMapG2 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        
        if modName == 'g3':
        
            lineMapG3 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        
        for j in range(0,len(lineBPT['BIN_ID'])):
            match_bin = np.where(tabGen['BIN_ID']==lineBPT['BIN_ID'][j])[0]

            for index in match_bin:
                if ~np.isnan(hbetaData[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])]):
                    
                    lineMapG1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j]['cDist-OIIIG1']

                    if modName != 'g1':

                        lineMapG2[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j]['cDist-OIIIG2']

                    if modName == 'g3':
                        lineMapG3[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j]['cDist-OIIIG3'] #TOREVIEW!!!!

                    lineMapToT[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j]['cDist-OIIIToT']
                    

        lineMapHead['BUNIT'] = 'cDistance'
        outBPT = bptDir+'cDist-OIIIG1.fits'
        fits.writeto(outBPT,lineMapG1,lineMapHead,overwrite=True)
    
        if modName != 'g1':
            
            #print('\n\t************* --- GuFo : ERROR --- **************\n')
            #outBPTg2 = bptDir+'BPT-'+str(lineBPT.dtype.names[i+numCols-1])+'.fits'
            outBPTG2 = bptDir+'BPT-cDist-OIIIG2.fits'
            outBPTToT = bptDir+'BPT-cDist-OIIIToT.fits'

            fits.writeto(outBPTG2,lineMapG2,lineMapHead,overwrite=True)
            fits.writeto(outBPTToT,lineMapToT,lineMapHead,overwrite=True)

            if modName == 'g3':
                outBPTG3 = bptDir+'BPT-cDist-OIIIG3.fits'
                fits.writeto(outBPTG3,lineMapG3,lineMapHead,overwrite=True)

        if cfg_par['lineRatios']['cDistPlot'] == True:
            
            bpt.cDistIM(cfg_par,outBPT)

            if modName != 'g1':
                bpt.cDistIM(cfg_par,outBPTG2)
                bpt.cDistIM(cfg_par,outBPTToT)
            elif modName=='g3':

                bpt.cDistIM(cfg_par,outBPTG3)

        return