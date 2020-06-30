#!/usr/bin/env python3.6
import os, sys, math
import yaml

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult

from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np
import numpy.ma as ma

import shutil



import cvPlay
import gufo as gf

#gf = gufo.gufo()
cvP = cvPlay.convert()


class tplay(object):

    def openLineList(self,cfg_par):
        
        workDir = cfg_par['general']['workdir']
       
        lineList = workDir+cfg_par['general']['lineListName']
        lineInfo = ascii.read(lineList) 

        #mask line list 
        index = np.where(lineInfo['Fit'] == 0)

        indexLines = np.where(lineInfo['Fit'] == 1)
        #idxMin = 

        fltr =  np.array(index)[0]
        lineInfo.remove_rows(list(fltr))

        lenTable = len(lineInfo['ID'])
        dltSigmaMinAng = np.zeros([lenTable])
        dltSigmaMaxAngG1 = np.zeros([lenTable])
        dltSigmaMaxAngG2 = np.zeros([lenTable])
        dltSigmaMaxAngG3 = np.zeros([lenTable])

        dltV12Ang = np.zeros([lenTable])
        dltSigma12Ang = np.zeros([lenTable])

        dltV13Ang = np.zeros([lenTable])
        dltSigma13Ang = np.zeros([lenTable])
        lineRange = np.zeros([lenTable])
        cenRange = np.zeros([lenTable])

        #ampThresh = np.zeros([lenTable])


        for i in range(0,lenTable):

            lambdaRest = lineInfo['Wave'][i]
            
            lineRange[i] = cvP.vRadLambda(lineInfo['lineRange'][i],
                lambdaRest)-lambdaRest    
            cenRange[i] = cvP.vRadLambda(lineInfo['cenRange'][i],
                lambdaRest)-lambdaRest    
            deltaV12 = np.log(cvP.vRadLambda(cfg_par['gFit']['dltV12'],
                lambdaRest))
            deltaV12 -= np.log(lambdaRest)       
            deltaV13 =np.log(cvP.vRadLambda(cfg_par['gFit']['dltV13'],
                lambdaRest))
            deltaV13 -= np.log(lambdaRest)


            deltaSigmaMin = np.log(cvP.vRadLambda(cfg_par['gFit']['sigmaMin'],
                lambdaRest))
            deltaSigmaMin -= np.log(lambdaRest)   
            deltaSigmaMaxG1 = np.log(cvP.vRadLambda(cfg_par['gFit']['sigmaMaxG1'],
                lambdaRest))
            deltaSigmaMaxG1 -= np.log(lambdaRest)   
            deltaSigmaMaxG2 = np.log(cvP.vRadLambda(cfg_par['gFit']['sigmaMaxG2'],
                lambdaRest))
            deltaSigmaMaxG2 -= np.log(lambdaRest)   
            deltaSigmaMaxG3 = np.log(cvP.vRadLambda(cfg_par['gFit']['sigmaMaxG3'],
                lambdaRest))
            deltaSigmaMaxG3 -= np.log(lambdaRest)            
            
            deltaSigma12 = np.log(cvP.vRadLambda(cfg_par['gFit']['dltSigma12'],
                lambdaRest))
            deltaSigma12 -=  np.log(lambdaRest)
            deltaSigma13 = np.log(cvP.vRadLambda(cfg_par['gFit']['dltSigma13'],
                lambdaRest))
            deltaSigma13 -= np.log(lambdaRest)


            dltSigmaMinAng[i] = deltaSigmaMin
            dltSigmaMaxAngG1[i] = deltaSigmaMaxG1
            dltSigmaMaxAngG2[i] = deltaSigmaMaxG2
            dltSigmaMaxAngG3[i] = deltaSigmaMaxG3

            dltV12Ang[i] = deltaV12            
            dltSigma12Ang[i] = deltaSigma12
            dltV13Ang[i] = deltaV13
            dltSigma13Ang[i] = deltaSigma13

            #ampThresh[i] = lineInfo['ampThresh'][i]

        dltSigmaMinCol = Column(name='deltaSigmaAng_Min', data=dltSigmaMinAng)        
        dltSigmaMaxColG1 = Column(name='deltaSigmaAng_MaxG1', data=dltSigmaMaxAngG1)        
        dltSigmaMaxColG2 = Column(name='deltaSigmaAng_MaxG2', data=dltSigmaMaxAngG2)        
        dltSigmaMaxColG3 = Column(name='deltaSigmaAng_MaxG3', data=dltSigmaMaxAngG3)        

        dltV12Col = Column(name='deltaVAng_12', data=dltV12Ang)
        dltSigma12Col = Column(name='deltaSigmaAng_12', data=dltSigma12Ang)
        dltV13Col = Column(name='deltaVAng_13', data=dltV13Ang)
        dltSigma13Col = Column(name='deltaSigmaAng_13', data=dltSigma13Ang)
        lineRangeCol = Column(name='lineRangeAng', data=lineRange)
        cenRangeCol = Column(name='cenRangeAng', data=cenRange)
        #ampThreshCol = Column(name='ampThresh', data=ampThresh)


        lineInfo.add_column(dltSigmaMinCol)
        lineInfo.add_column(dltSigmaMaxColG1)
        lineInfo.add_column(dltSigmaMaxColG2)
        lineInfo.add_column(dltSigmaMaxColG3)

        lineInfo.add_column(dltV12Col)
        lineInfo.add_column(dltSigma12Col)
        lineInfo.add_column(dltV13Col)
        lineInfo.add_column(dltSigma13Col)
        lineInfo.add_column(lineRangeCol)
        lineInfo.add_column(cenRangeCol)
        #lineInfo.add_column(ampThreshCol)


        return lineInfo

    def openTablesPPXF(self,cfg_par,tableBin,tableSpec):
        
        crPix1=cfg_par['starSub']['pixX']
        crPix2=cfg_par['starSub']['pixY']
      
        tab = fits.open(tableBin)
        head = tab[0].header
        headTab = tab[1].header
        dataTab = tab[1].data    
        head['CRPIX1'] = crPix1
        head['CRPIX2'] = crPix2 
        
        xMin = np.min(dataTab['X'])
        xMax = np.max(dataTab['X'])

        shapeX = (xMax-xMin)/head['PIXSIZE']

        yMin = np.min(dataTab['Y'])
        yMax = np.max(dataTab['Y'])

        shapeY = (yMax-yMin)/head['PIXSIZE']

        xAxis = np.arange(xMin, xMax,head['PIXSIZE'])
        yAxis = np.arange(yMin, yMax,head['PIXSIZE'])
        
        tab = fits.open(tableSpec)
        tab.info()
        dataSpec = tab[1].data
        specExp = tab[2].data
        wave = [item for t in specExp for item in t] 
        noiseBin = dataSpec['ESPEC']
        pxSize = head['PIXSIZE']/3600.

        return wave,xAxis,yAxis,pxSize,noiseBin,dataTab


    def openPPXFforSubtraction(self,cfg_par,tableBin,tableSpec,tableStar):
        
        crPix1=cfg_par['starSub']['pixX']
        crPix2=cfg_par['starSub']['pixY']
      
        tab = fits.open(tableBin)
        head = tab[0].header
        headTab = tab[1].header
        dataTab = tab[1].data    
        head['CRPIX1'] = crPix1
        head['CRPIX2'] = crPix2 
        
        xMin = np.min(dataTab['X'])
        xMax = np.max(dataTab['X'])

        shapeX = (xMax-xMin)/head['PIXSIZE']

        yMin = np.min(dataTab['Y'])
        yMax = np.max(dataTab['Y'])

        shapeY = (yMax-yMin)/head['PIXSIZE']

        xAxis = np.arange(xMin, xMax,head['PIXSIZE'])
        yAxis = np.arange(yMin, yMax,head['PIXSIZE'])
        tab = fits.open(tableSpec)
        dataSpec = tab[1].data
        specExp = tab[2].data
        wave = [item for t in specExp for item in t] 

        pxSize = head['PIXSIZE']/3600.

        noiseBin = dataSpec['ESPEC']

        tabStar = fits.open(tableStar)
        dataStar = tabStar[1].data


        return wave,xAxis,yAxis,pxSize,noiseBin,dataTab,dataSpec,dataStar

    def openVorLineOutput(self,cfg_par,tableBin,tableSpec):
        
        crPix1=cfg_par['starSub']['pixX']
        crPix2=cfg_par['starSub']['pixY']
      
        tab = fits.open(tableBin)
        head = tab[0].header
        headTab = tab[1].header
        dataTab = tab[1].data    
        head['CRPIX1'] = crPix1
        head['CRPIX2'] = crPix2 
        head['PIXSIZE'] = head['PIXSIZE']
        xMin = np.min(dataTab['X'])
        xMax = np.max(dataTab['X'])

        shapeX = (xMax-xMin)/head['PIXSIZE']

        yMin = np.min(dataTab['Y'])
        yMax = np.max(dataTab['Y'])

        shapeY = (yMax-yMin)/head['PIXSIZE']

        xAxis = np.arange(xMin, xMax,head['PIXSIZE'])
        yAxis = np.arange(yMin, yMax+head['PIXSIZE'],head['PIXSIZE'])
        tab = fits.open(tableSpec)
        dataSpec = tab[1].data
        specExp = tab[2].data
        wave = [item for t in specExp for item in t] 

        pxSize = head['PIXSIZE']/3600.

        noiseBin = dataSpec['ESPEC']

        return wave,xAxis,yAxis,pxSize,noiseBin,dataTab,dataSpec

    def makePixelTable(self,cfg_par):

        tab = fits.open(cfg_par['general']['outVorTableName'])
        
        headTab = tab[0].header
        dataTab = tab[0].data  

        NSPAX = np.zeros(0,len(dataTab['NSPAX']))+1.
        BIN_ID = dataTab['ID'].copy()

        nam = tuple (['ID', 'BIN_ID', 'X', 'Y', 'PixX', 'PixY', 'NSPAX'])
        tableArr = np.array([dataTab['ID'],BIN_ID,dataTab['X'],dataTab['Y'],dataTab['PixX'],dataTab['PixY'],NSPAX], dtype={'names':nam,
                          'formats':( 'i4', 'i4', 'f8', 'f8', 'i4', 'i4', 'i4', 'i4')})

        cols = []

        cols.append(fits.Column(name='ID',        format='J',   array=dataTab['ID']     ))
        cols.append(fits.Column(name='BIN_ID',    format='J',   array=BIN_ID            ))
        cols.append(fits.Column(name='X',         format='D',   array=dataTab['X']      ))
        cols.append(fits.Column(name='Y',         format='D',   array=dataTab['Y']      ))
        cols.append(fits.Column(name='PixX',      format='D',   array=dataTab['PixX']   ))
        cols.append(fits.Column(name='PixY',      format='D',   array=dataTab['PixY']   ))
        cols.append(fits.Column(name='XBIN',      format='D',   array=dataTab['PixX']   ))
        cols.append(fits.Column(name='YBIN',      format='D',   array=dataTab['PixY']   ))        
        cols.append(fits.Column(name='NSPAX',     format='J',   array=NSPAX             ))

        tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        tbhdu.writeto(cfg_par['general']['outVorLineTableName'], overwrite=True)

        return 0

    def makeInputArrays(self,cfg_par,lineInfo, Xdim,Ydim):

        binID = np.zeros([Ydim,Xdim],dtype=int)

        nam = tuple (['ID', 'BIN_ID', 'X', 'Y', 'PixX', 'PixY'])
        binArr = np.zeros([Ydim*Xdim], dtype={'names':nam,
                          'formats':('i4', 'i4', 'i4', 'f8', 'f8', 'i4', 'i4')})
        nam = tuple(['BIN_ID', 'fitSuccess', 'redChi', 'aic', 'bic', 'nData', 'nVariables', 'nFev'])
        fitResArr = np.zeros([Ydim*Xdim], dtype={'names':nam,
                          'formats':( 'i4', '?', 'f8', 'f8', 'f8', 'i4', 'i4', 'i4')})

        lineNameList = []
        frmList = []
        lineNameList.append('BIN_ID')
        frmList.append('i4')
        for i in range (0,len(lineInfo['ID'])):
            lineName = str(lineInfo['Name'][i])+str(int(lineInfo['Wave'][i]))
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")
            

            lineNameList.append(lineName)
            frmList.append('i4')
            lineNameList.append('g1_Amp_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_Height_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_Centre_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_Sigma_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_FWHM_'+lineName)
            frmList.append('f8')

            lineNameList.append('g1_SN_'+lineName)
            frmList.append('f8')
            
            if cfg_par['gFit']['modName'] == 'g2':
                
                lineNameList.append('g2_Amp_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_Height_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_Centre_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_Sigma_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_FWHM_'+lineName)
                frmList.append('f8')

                lineNameList.append('g2_SN_'+lineName)
                frmList.append('f8')

            
                if cfg_par['gFit']['modName'] == 'g3':

                    lineNameList.append('g3_Amp_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_Height_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_Centre_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_Sigma_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_FWHM_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_SN_'+lineName)
                    frmList.append('f8')
        
        if cfg_par['gFit']['modName'] == 'g1':
            lineArr = np.zeros([Ydim*Xdim], dtype={'names':(lineNameList), 'formats':(frmList)})
        elif cfg_par['gFit']['modName'] == 'g2':
            lineArr = np.zeros([Ydim*Xdim], dtype={'names':(lineNameList), 'formats':(frmList)})
        elif cfg_par['gFit']['modName'] == 'g3':
            lineArr = np.zeros([Ydim*Xdim], dtype={'names':(lineNameList), 'formats':(frmList)})

        return binID, binArr, fitResArr, lineArr

    def makeInputArraysMP(self,cfg_par,lineInfo,vorBinInfo):

        nam = tuple (['ID', 'BIN_ID', 'X', 'Y', 'PixX', 'PixY'])
        binArr = np.zeros([len(vorBinInfo['ID'])], dtype={'names':nam,
                          'formats':('i4', 'i4', 'i4', 'f8', 'f8', 'i4', 'i4')})
        nam = tuple(['BIN_ID', 'fitSuccess', 'redChi', 'aic', 'bic', 'nData', 'nVariables', 'nFev'])
        fitResArr = np.zeros([len(vorBinInfo['ID'])], dtype={'names':nam,
                          'formats':( 'i4', '?', 'f8', 'f8', 'f8', 'i4', 'i4', 'i4')})

        lineNameList = []
        frmList = []
        lineNameList.append('BIN_ID')
        frmList.append('i4')
        lineNameList.append('noiseBin')
        frmList.append('f8')
        for i in range (0,len(lineInfo['ID'])):
            lineName = str(lineInfo['Name'][i])+str(int(lineInfo['Wave'][i]))
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")
            

            lineNameList.append(lineName)
            frmList.append('i4')
            lineNameList.append('g1_Amp_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_Height_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_Centre_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_SigMeas_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_SigIntr_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_dLambda_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_FWHM_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_SN_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_centre_'+lineName)         
            frmList.append('f8')
            lineNameList.append('g1_sigLambda_'+lineName)         
            frmList.append('f8')
            
            if cfg_par['gFit']['modName'] == 'g2':
                
                lineNameList.append('g2_Amp_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_Height_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_Centre_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_SigMeas_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_SigIntr_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_FWHM_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_SN_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_centre_'+lineName)         
                frmList.append('f8')
                lineNameList.append('g2_sigLambda_'+lineName)         
                frmList.append('f8')
            
                if cfg_par['gFit']['modName'] == 'g3':

                    lineNameList.append('g3_Amp_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_Height_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_Centre_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_SigMeas_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_SigIntr_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_FWHM_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_SN_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_centre_'+lineName)         
                    frmList.append('f8')
                    lineNameList.append('g3_sigLambda_'+lineName)         
                    frmList.append('f8')

        if cfg_par['gFit']['modName'] == 'g1':
            lineArr = np.zeros([len(vorBinInfo['ID'])], dtype={'names':(lineNameList), 'formats':(frmList)})
        elif cfg_par['gFit']['modName'] == 'g2':
            lineArr = np.zeros([len(vorBinInfo['ID'])], dtype={'names':(lineNameList), 'formats':(frmList)})
        elif cfg_par['gFit']['modName'] == 'g3':
            lineArr = np.zeros([len(vorBinInfo['ID'])], dtype={'names':(lineNameList), 'formats':(frmList)})

        return binArr, fitResArr, lineArr


    def updateBinArray(self,cfg_par,binArr,vorBinInfo,index,i,j,counter):
  
        binArr['BIN_ID'][counter] = vorBinInfo['BIN_ID'][index]
        binArr['ID'][counter] = vorBinInfo['ID'][index]
        binArr['X'][counter] = vorBinInfo['X'][index]
        binArr['Y'][counter] = vorBinInfo['Y'][index]
        binArr['PixX'][counter] = int(i)
        binArr['PixY'][counter] = int(j)

        return binArr

    def updateFitArray(self,cfg_par,fitResArr,result,binIDName,counter):

        aic = result.aic
        bic = result.bic
        redchi = result.redchi
        success = result.success
        ndata = result.ndata
        nvarys = result.nvarys
        nfev = result.nfev
        success = result.success
        fitResArr['BIN_ID'][counter] = binIDName
        fitResArr['fitSuccess'][counter] = success
        fitResArr['redChi'][counter] = redchi
        fitResArr['aic'][counter] = aic
        fitResArr['nData'][counter] = bic
        fitResArr['nVariables'][counter] = nvarys
        fitResArr['nFev'][counter] = nfev
        fitResArr['nData'][counter] = ndata
        

        return fitResArr


    def updateLineArray(self,cfg_par,wave,lineArr,result,noiseValue,lineInfo,binIDName,counter):
        
        fitRes = result.params.valuesdict()
        dLambda = cvP.specRes(cfg_par)

        modName = cfg_par['gFit']['modName']
        lineArr['BIN_ID'][counter] = binIDName
        lineArr['noiseBin'][counter] = noiseValue

        for ii in range(0,len(lineInfo['ID'])):

            lineName = str(lineInfo['Name'][ii])
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")
            
            lineName = lineName+str(int(lineInfo['Wave'][ii]))


            intR = fitRes['Wintln'+str(ii)]

            waveInRed = cfg_par['general']['redshift']*lineInfo['Wave'][ii]+lineInfo['Wave'][ii]

            indexWaveInRed = int(np.where(abs(np.exp(wave)-waveInRed)==abs(np.exp(wave)-waveInRed).min())[0])
            

            dLIn = dLambda[indexWaveInRed]
            dLIn = np.log(waveInRed+dLIn/2.)-np.log(waveInRed-dLIn/2.)

            if modName == 'g1':

                amp = fitRes['g1ln'+str(ii)+'_amplitude']
                ctr = fitRes['g1ln'+str(ii)+'_center']
                sig = fitRes['g1ln'+str(ii)+'_sigma']

                fwhm = fitRes['g1ln'+str(ii)+'_fwhm']
                height = fitRes['g1ln'+str(ii)+'_height']

                sigmaInt = np.sqrt(np.power(sig,2)-np.power(dLIn,2))

                g1Ctr = cvP.lambdaVRad(np.exp(ctr),lineInfo['Wave'][ii])

                g1SigmaInt = cvP.lambdaVRad(np.exp(ctr+sigmaInt),lineInfo['Wave'][ii])-g1Ctr
                g1Sigma = cvP.lambdaVRad(np.exp(ctr+sig),lineInfo['Wave'][ii])-g1Ctr            
                g1FWHM = cvP.lambdaVRad(np.exp(ctr+fwhm),lineInfo['Wave'][ii])-g1Ctr
                g1dL = cvP.lambdaVRad(np.exp(ctr+dLambda[indexWaveInRed]),lineInfo['Wave'][ii])-g1Ctr

                #amp_err = result.params[modName+'ln'+str(i)+'_amplitude'].stderr
                #sig_err = result.params[modName+'ln'+str(i)+'_sigma'].stderr
                #g1SigmaErr = self.lambdaVRad(np.exp(sig_err),lineInfo['Wave'][i])
                #cen_err = result.params[modName+'ln'+str(i)+'_center'].stderr  
                #g1CtrErr = self.lambdaVRad(np.exp(cen_err),lineInfo['Wave'][i])
                lineArr[lineName][counter] = int(lineInfo['Wave'][ii])     
                
                lineArr['g1_Amp_'+lineName][counter] = amp
                lineArr['g1_Height_'+lineName][counter] = height
                lineArr['g1_SN_'+lineName][counter]=height/noiseValue         
                lineArr['g1_Centre_'+lineName][counter] = g1Ctr
                lineArr['g1_SigMeas_'+lineName][counter] = g1Sigma
                lineArr['g1_SigIntr_'+lineName][counter] = g1SigmaInt
                lineArr['g1_FWHM_'+lineName][counter] = g1FWHM

                lineArr['g1_dLambda_'+lineName][counter] = g1dL
                lineArr['g1_centre_'+lineName][counter]=ctr         
                lineArr['g1_sigLambda_'+lineName][counter]=sig       

            elif modName == 'g2':

                sigTmp1 = fitRes['g1ln'+str(ii)+'_sigma']
                sigTmp2 = fitRes['g2ln'+str(ii)+'_sigma']

                if sigTmp1 <= sigTmp2:

                    amp1 = fitRes['g1ln'+str(ii)+'_amplitude']
                    ctr1 = fitRes['g1ln'+str(ii)+'_center']
                    sig1 = fitRes['g1ln'+str(ii)+'_sigma']

                    fwhm1 = fitRes['g1ln'+str(ii)+'_fwhm']
                    height1 = fitRes['g1ln'+str(ii)+'_height']

                    amp2 = fitRes['g2ln'+str(ii)+'_amplitude']
                    ctr2 = fitRes['g2ln'+str(ii)+'_center']
                    sig2 = fitRes['g2ln'+str(ii)+'_sigma']
                    fwhm2 = fitRes['g2ln'+str(ii)+'_fwhm']
                    height2 = fitRes['g2ln'+str(ii)+'_height']
                
                else:
                    
                    amp2 = fitRes['g1ln'+str(ii)+'_amplitude']
                    ctr2 = fitRes['g1ln'+str(ii)+'_center']
                    sig2 = fitRes['g1ln'+str(ii)+'_sigma']
                    fwhm2 = fitRes['g1ln'+str(ii)+'_fwhm']
                    height2 = fitRes['g1ln'+str(ii)+'_height']

                    amp1 = fitRes['g2ln'+str(ii)+'_amplitude']
                    ctr1 = fitRes['g2ln'+str(ii)+'_center']
                    sig1 = fitRes['g2ln'+str(ii)+'_sigma']
                    fwhm1 = fitRes['g2ln'+str(ii)+'_fwhm']
                    height1 = fitRes['g2ln'+str(ii)+'_height']

                sigmaInt1 = np.sqrt(np.power(sig1,2)-np.power(dLIn,2))

                g1Ctr = cvP.lambdaVRad(np.exp(ctr1),lineInfo['Wave'][ii])

                g1SigmaInt = cvP.lambdaVRad(np.exp(ctr1+sigmaInt1),lineInfo['Wave'][ii])-g1Ctr
                g1Sigma = cvP.lambdaVRad(np.exp(ctr1+sig1),lineInfo['Wave'][ii])-g1Ctr            
                g1FWHM = cvP.lambdaVRad(np.exp(ctr1+fwhm1),lineInfo['Wave'][ii])-g1Ctr
                g1dL = cvP.lambdaVRad(np.exp(ctr1+dLambda[indexWaveInRed]),lineInfo['Wave'][ii])-g1Ctr                    


                sigmaInt2 = np.sqrt(np.power(sig2,2)-np.power(dLIn,2))

                g2Ctr = cvP.lambdaVRad(np.exp(ctr2),lineInfo['Wave'][ii])
                g2SigmaInt = cvP.lambdaVRad(np.exp(ctr2+sigmaInt2),lineInfo['Wave'][ii])-g2Ctr
                g2Sigma = cvP.lambdaVRad(np.exp(ctr2+sig2),lineInfo['Wave'][ii])-g2Ctr
                g2FWHM = cvP.lambdaVRad(np.exp(ctr2+fwhm2),lineInfo['Wave'][ii])-g2Ctr

                #amp_err = result.params[modName+'ln'+str(i)+'_amplitude'].stderr
                #sig_err = result.params[modName+'ln'+str(i)+'_sigma'].stderr
                #g1SigmaErr = self.lambdaVRad(np.exp(sig_err),lineInfo['Wave'][i])
                #cen_err = result.params[modName+'ln'+str(i)+'_center'].stderr  
                #g1CtrErr = self.lambdaVRad(np.exp(cen_err),lineInfo['Wave'][i])

                lineArr[lineName][counter] = int(lineInfo['Wave'][ii])     
                
                lineArr['g1_Amp_'+lineName][counter] = amp1
                lineArr['g1_Height_'+lineName][counter] = height1
                lineArr['g1_SN_'+lineName][counter]=height1/noiseValue         
                lineArr['g1_Centre_'+lineName][counter] = g1Ctr
                lineArr['g1_SigMeas_'+lineName][counter] = g1Sigma
                lineArr['g1_SigIntr_'+lineName][counter] = g1SigmaInt
                lineArr['g1_FWHM_'+lineName][counter] = g1FWHM

                lineArr['g1_dLambda_'+lineName][counter] = g1dL
                lineArr['g1_centre_'+lineName][counter]=ctr1         
                lineArr['g1_sigLambda_'+lineName][counter]=sig1  
          
                lineArr['g2_Amp_'+lineName][counter] = amp2
                lineArr['g2_Height_'+lineName][counter] = height2
                lineArr['g2_Centre_'+lineName][counter] = g2Ctr
                lineArr['g2_SigMeas_'+lineName][counter] = g2Sigma
                lineArr['g2_SigIntr_'+lineName][counter] = g2SigmaInt
                lineArr['g2_FWHM_'+lineName][counter] = g2FWHM
                lineArr['g2_SN_'+lineName][counter]=height2/noiseValue         
                lineArr['g2_sigLambda_'+lineName][counter]=sig2  

            elif modName == 'g3':

                    amp = fitRes['g3ln'+str(ii)+'_amplitude']
                    ctr = fitRes['g3ln'+str(ii)+'_center']
                    sig = fitRes['g3intln'+str(ii)+'_sigma']
                    fwhm = fitRes['g3ln'+str(ii)+'_fwhm']
                    height = fitRes['g3ln'+str(ii)+'_height']

                    sigmaInt = np.sqrt(np.power(sig,2)-np.power(dLIn,2))


                    g3Ctr = cvP.lambdaVRad(np.exp(ctr),lineInfo['Wave'][ii])
                    g3Sigma = cvP.lambdaVRad(np.exp(ctr+sig),lineInfo['Wave'][ii])-g3Ctr
                    g3SigmaInt = cvP.lambdaVRad(np.exp(ctr+sigmaInt),lineInfo['Wave'][ii])-g3Ctr
                    g3Sigma = cvP.lambdaVRad(np.exp(ctr+sig),lineInfo['Wave'][ii])-g3Ctr
                    g3FWHM = cvP.lambdaVRad(np.exp(ctr+fwhm),lineInfo['Wave'][ii])-g3Ctr

                    #amp_err = result.params[modName+'ln'+str(i)+'_amplitude'].stderr
                    #sig_err = result.params[modName+'ln'+str(i)+'_sigma'].stderr
                    #g1SigmaErr = self.lambdaVRad(np.exp(sig_err),lineInfo['Wave'][i])
                    #cen_err = result.params[modName+'ln'+str(i)+'_center'].stderr  
                    #g1CtrErr = self.lambdaVRad(np.exp(cen_err),lineInfo['Wave'][i])
              
                    lineArr['g3_Amp_'+lineName][counter] = amp
                    lineArr['g3_Height_'+lineName][counter] = height
                    lineArr['g3_Centre_'+lineName][counter] = g3Ctr
                    lineArr['g3_SigMeas_'+lineName][counter] = g3Sigma
                    lineArr['g3_SigIntr_'+lineName][counter] = g3SigmaInt
                    lineArr['g3_dLambda_'+lineName][counter] = g3dL
                    lineArr['g3_FWHM_'+lineName][counter] = g3FWHM
                    lineArr['g3_SN_'+lineName][counter]=height/noiseValue         
    

        return lineArr

    def saveOutputTable(self,cfg_par, binArr, fitResArr, lineArr):
        
        #outTableName = cfg_par['general']['runNameDir']+'/gPlayOut1.fits'
        modNameList = cfg_par['gFit']['modName']

        if os.path.exists(cfg_par['general']['outTableName']):
            hdul = fits.open(cfg_par['general']['outTableName'])
            t2 = fits.BinTableHDU.from_columns(fitResArr,name='FitRes_'+modNameList)
            hdul.append(t2)  
            t3 = fits.BinTableHDU.from_columns(lineArr,name='LineRes_'+modNameList)
            hdul.append(t3)  
        else:    
            hdr = fits.Header()
            hdr['COMMENT'] = "Here are the outputs of gPlay"
            hdr['COMMENT'] = "Ext 1 = binInfo Ext 2 = fit result Ext 3 = line parameters"
            
            empty_primary = fits.PrimaryHDU(header=hdr)
           
            t1 = fits.BinTableHDU.from_columns(binArr,name='BinInfo')  
            hdul = fits.HDUList([empty_primary,t1])        

            t2 = fits.BinTableHDU.from_columns(fitResArr,name='FitRes_'+modNameList)
            hdul.append(t2)  

            t3 = fits.BinTableHDU.from_columns(lineArr,name='LineRes_'+modNameList)
            hdul.append(t3)  

        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)

        return

    def cleanTable(self,cfg_par):
        
        hdul = fits.open(cfg_par['general']['outTableName'])

        if cfg_par['gFit']['modName'] == 'g1':
            hdl = fits.HDUList([hdul[0],hdul['BININFO'],hdul['FitRes_g1'],hdul['LineRes_g1']])
        elif cfg_par['gFit']['modName'] == 'g2':
            hdl = fits.HDUList([hdul[0],hdul['BININFO'],hdul['FitRes_g1'],hdul['LineRes_g1'],hdul['FitRes_g2'],hdul['LineRes_g2']])

        hdl.writeto(cfg_par['general']['outTableName'],overwrite=True)

        return


    def reorderTable(self,cfg_par):

        lineInfo = self.openLineList(cfg_par)

        hdul = fits.open(cfg_par['general']['outTableName'])

        lines = hdul['LineRes_g2'].data

        for i in range(0,len(lines['BIN_ID'])):

            sigTmp1 = lines['g1_SigIntr_OIII5006'][i]
            sigTmp2 = lines['g2_SigIntr_OIII5006'][i]

            if sigTmp1 <= sigTmp2:
                pass
            else:
                
                for ii in range(0,len(lineInfo['ID'])):

                    lineName = str(lineInfo['Name'][ii])
                    if '[' in lineName:
                        lineName = lineName.replace("[", "")
                        lineName = lineName.replace("]", "")
                    
                    lineName = lineName+str(int(lineInfo['Wave'][ii]))

                    amp2 = lines['g1_Amp_'+lineName][i]
                    height2 = lines['g1_Height_'+lineName][i]
                    g2Ctr = lines['g1_Centre_'+lineName][i]
                    g2Sigma = lines['g1_SigMeas_'+lineName][i]
                    g2SigmaInt = lines['g1_SigIntr_'+lineName][i]
                    g2FWHM = lines['g1_FWHM_'+lineName][i]
                    g2SN = lines['g1_SN_'+lineName][i]
                    g2Centre = lines['g1_centre_'+lineName][i]
                    g2SigL = lines['g1_sigLambda_'+lineName][i]
                    
                    lines['g1_Amp_'+lineName][i] = lines['g2_Amp_'+lineName][i]
                    lines['g1_Height_'+lineName][i] = lines['g2_Height_'+lineName][i]
                    lines['g1_SN_'+lineName][i]=lines['g2_SN_'+lineName][i]        
                    lines['g1_Centre_'+lineName][i] = lines['g2_Centre_'+lineName][i]
                    lines['g1_SigMeas_'+lineName][i] = lines['g2_SigMeas_'+lineName][i]
                    lines['g1_SigIntr_'+lineName][i] = lines['g2_SigIntr_'+lineName][i]
                    lines['g1_FWHM_'+lineName][i] = lines['g2_FWHM_'+lineName][i]

                    lines['g1_centre_'+lineName][i]=lines['g2_centre_'+lineName][i]       
                    lines['g1_sigLambda_'+lineName][i]=lines['g2_sigLambda_'+lineName][i]
          
                    lines['g2_Amp_'+lineName][i] = amp2
                    lines['g2_Height_'+lineName][i] = height2
                    lines['g2_Centre_'+lineName][i] = g2Ctr
                    lines['g2_SigMeas_'+lineName][i] = g2Sigma
                    lines['g2_SigIntr_'+lineName][i] = g2SigmaInt
                    lines['g2_FWHM_'+lineName][i] = g2FWHM
                    lines['g2_SN_'+lineName][i]=g2SN         
                    lines['g2_sigLambda_'+lineName][i]= g2SigL 
                    lines['g2_centre_'+lineName][i]=g2Centre    

        hdl = fits.HDUList([hdul[0],hdul['BININFO'],hdul['FitRes_g1'],hdul['LineRes_g1'],hdul['Residuals_g1'],
            hdul['FitRes_g2'],hdul['LineRes_g2'],hdul['Residuals_g2']])
        
        hdl.writeto(cfg_par['general']['runNameDir']+'gPlayOutReord.fits',overwrite=True)


    def selectBestFit(self,cfg_par):

        print(cfg_par['bestFitSel']['tableNames'])
        
        tableNames = np.array(cfg_par['bestFitSel']['tableNames'])
        
        hdul = fits.open(cfg_par['general']['runNameDir']+tableNames[0])
        bins = hdul['lineRes_g2'].data['BIN_ID']
        nrows = hdul['lineRes_g2'].data.shape[0]

        linesG1 = hdul['lineRes_g1'].data
        resG1 = hdul['residuals_g1'].data
        fitG1 = hdul['fitRes_g1'].data
        ancG1 = hdul['ancelsg1'].data

        hdul = fits.open(cfg_par['general']['runNameDir']+tableNames[1])

        linesG2R1 = hdul['lineRes_g2'].data
        resG2R1 = hdul['residuals_g2'].data
        fitG2R1 = hdul['fitRes_g2'].data
        ancG2R1 = hdul['ancelsg2'].data

        hdul = fits.open(cfg_par['general']['runNameDir']+tableNames[2])       
        linesG2R2 = hdul['lineRes_g2'].data
        resG2R2 = hdul['residuals_g2'].data
        fitG2R2 = hdul['fitRes_g2'].data
        ancG2R2 = hdul['ancelsg2'].data


        hdul = fits.open(cfg_par['general']['runNameDir']+tableNames[3])
        linesG2R3 = hdul['lineRes_g2'].data
        resG2R3 = hdul['residuals_g2'].data
        fitG2R3 = hdul['fitRes_g2'].data
        ancG2R3 = hdul['ancelsg2'].data
        
        res=np.zeros([len(tableNames),nrows])


        #for i in range(len(tableNames)):

        res[0,:] = np.array(resG1['res_NII6583'])
        res[1,:] = np.array(resG2R1['res_NII6583'])
        res[2,:] = np.array(resG2R2['res_NII6583'])
        res[3,:] = np.array(resG2R3['res_NII6583'])

        bestFitTable = fits.BinTableHDU.from_columns(hdul['lineRes_g2'].columns, nrows=nrows,name='lineRes_g2')
        resTable = fits.BinTableHDU.from_columns(hdul['residuals_g2'].columns, nrows=nrows,name='residuals_g2')
        fitResTable = fits.BinTableHDU.from_columns(hdul['fitres_g2'].columns, nrows=nrows,name='fitRes_g2')
        ancTable = fits.BinTableHDU.from_columns(hdul['ancelsg2'].columns, nrows=nrows,name='ancelsg2')
        #print(fitResTable.columns.names,linesG1.columns.names)
        bestres = []
        
        for i in range(nrows):
            bestres.append(np.argmin(res[:,i]))
            print(res[:,i])
            if bestres[i] == 0:
                for colname in linesG1.columns.names:
                    bestFitTable.data[colname][i] = linesG1[colname][i]
                fitResTable.data[:][i] = fitG1[:][i]
                resTable.data[:][i] = resG1[:][i]
                ancTable.data[:][i] = ancG1[:][i]
            elif bestres[i] ==1:
                bestFitTable.data[:][i] = linesG2R1[:][i]
                fitResTable.data[:][i] = fitG2R1[:][i]
                resTable.data[:][i] = resG2R1[:][i]
                ancTable.data[:][i] = ancG2R2[:][i]
            elif bestres[i] ==2:
                bestFitTable.data[:][i] = linesG2R2[:][i]
                fitResTable.data[:][i] = fitG2R2[:][i]
                resTable.data[:][i] = resG2R2[:][i]
                ancTable.data[:][i] = ancG2R2[:][i]
            elif bestres[i] ==3:
                bestFitTable.data[:][i] = linesG2R3[:][i]
                fitResTable.data[:][i] = fitG2R3[:][i]
                resTable.data[:][i] = resG2R3[:][i]
                ancTable.data[:][i] = ancG2R3[:][i]

        #tot = np.column_stack(( resTable.data.columns,bestres))

        new_col = fits.ColDefs([fits.Column(name='bestFit', format='D', array=bestres)])
        
        orig_cols = resTable.data.columns
        
        hduBF = fits.BinTableHDU.from_columns(orig_cols + new_col)
        
        hdl = fits.HDUList([hdul[0],hdul['BININFO'],fitResTable,hduBF,resTable,ancTable])

        hdl.writeto(cfg_par['general']['runNameDir']+'gPlayOutBF.fits',overwrite=True)


    def binLineRatio(self,cfg_par,lineInfo):

        lineNameID=[]
        modName = cfg_par['gFit']['modName']
                #open line lineList

        for ii in range(0,len(lineInfo['ID'])):

            lineName = str(lineInfo['Name'][ii])
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")
            
            lineNameID.append(lineName+str(int(lineInfo['Wave'][ii])))
        
        lineThresh = float(lineInfo['SNThresh'][2])

        hdul = fits.open(cfg_par['general']['outTableName'])
        lines = hdul['LineRes_'+cfg_par['gFit']['modName']].data 
        linesG1 = hdul['LineRes_G1'].data

        residuals = hdul['Residuals_'+cfg_par['gFit']['modName']].data
        sigmaTable = hdul['LineRes_'+cfg_par['gFit']['modName']].data

        thresHold = residuals['SN_OIII5006']
        sigmaThresh = linesG1['g1_SigIntr_OIII5006']


        lineNameList=['BIN_ID']
        frmList=['i4']

        tot = lines['BIN_ID']

        lineNameID = np.array(lineNameID)
       

        index = np.where(thresHold<=lineThresh)[0]
        indexSigma = np.where(sigmaThresh>cfg_par['moments']['sigmaThresh'])

        if 'OIII5006' in lineNameID and 'Hb4861' in lineNameID:
            
            oIII = np.copy(lines['g1_Amp_'+'OIII5006'])
            oIII[index] = np.nan
            oIII[indexSigma] = np.nan

            hBeta = np.copy(lines['g1_Amp_'+'Hb4861'])
            hBeta[index] = np.nan
            hBeta[indexSigma] = np.nan
            
            lrOHbG1 = np.divide(oIII,hBeta)
            logOHbG1 = np.log10(lrOHbG1)

            tot = np.column_stack((tot,lrOHbG1,logOHbG1))

            lineNameList.append('G1-OIII5006/Hb4861')
            lineNameList.append('log_G1-OIII5006/Hb4861')

            frmList.append('f8')
            frmList.append('f8')


        if 'NII6583' in lineNameID and 'Ha6562' in lineNameID:
            
            NII = np.copy(lines['g1_Amp_'+'NII6583'])
            NII[index]=np.nan
            NII[indexSigma] = np.nan

            Halpha = np.copy(lines['g1_Amp_'+'Ha6562'])
            Halpha[index]=np.nan
            Halpha[indexSigma]=np.nan

            lrNIIHaG1 = np.divide(NII,Halpha)
            logNIIHaG1 = np.log10(lrNIIHaG1)
            tot = np.column_stack((tot,lrNIIHaG1,logNIIHaG1))
            lineNameList.append('G1-NII6583/Ha6562')
            lineNameList.append('log_G1-NII6583/Ha6562')

            frmList.append('f8')
            frmList.append('f8')

        if 'OI6300' in lineNameID and 'Ha6562' in lineNameID:
            
            OI = np.copy(lines['g1_Amp_'+'OI6300'])
            OI[index] = np.nan
            OI[indexSigma] = np.nan

            Halpha = np.copy(lines['g1_Amp_'+'Ha6562'])
            Halpha[index]=np.nan
            Halpha[indexSigma]=np.nan

            lrOIHaG1 = np.divide(OI,Halpha)
            logOIHaG1 = np.log10(lrOIHaG1)            
            tot = np.column_stack((tot,lrOIHaG1,logOIHaG1))
            lineNameList.append('G1-OI6300/Ha6562')
            lineNameList.append('log_G1-OI6300/Ha6562')

            frmList.append('f8')
            frmList.append('f8')

        if 'SII6716' in lineNameID and 'Ha6562' in lineNameID:
            
            SII1 = np.copy(lines['g1_Amp_'+'SII6716'])
            SII1[index] = np.nan
            SII1[indexSigma] = np.nan

            SII2 = np.copy(lines['g1_Amp_'+'SII6730'])
            SII2[index] = np.nan
            SII2[indexSigma] = np.nan

            Halpha = np.copy(lines['g1_Amp_'+'Ha6562'])
            Halpha[index] = np.nan
            Halpha[indexSigma]=np.nan


            lrSIIHaG1 = np.divide((SII1+SII2),Halpha)
            logSIIHaG1 = np.log10(lrSIIHaG1)            
            
            tot = np.column_stack((tot,lrSIIHaG1,logSIIHaG1))
            
            lineNameList.append('G1-SII6716/Ha6562')
            lineNameList.append('log_G1-SII6716/Ha6562')

            frmList.append('f8')
            frmList.append('f8')


        if modName != 'g1':


            if 'OIII5006' in lineNameID and 'Hb4861' in lineNameID:

                oIIIG2 = np.copy(lines['g2_Amp_'+'OIII5006'])
                oIIIG2[index] = np.nan
                oIIIG2[indexSigma] = np.nan

                hBetaG2 = np.copy(lines['g2_Amp_'+'Hb4861'])
                hBetaG2[index] = np.nan
                hBetaG2[indexSigma] = np.nan

                lrOHbG2 = np.divide(oIIIG2,hBetaG2)
                logOHbG2 = np.log10(lrOHbG2)            

                lrOHb = np.divide((oIII+oIIIG2),(hBeta+hBetaG2))
                logOHb = np.log10(lrOHb)            
                
                tot = np.column_stack((tot,lrOHbG2,logOHbG2,lrOHb,logOHb))
                
                lineNameList.append('G2-OIII5006/Hb4861')
                lineNameList.append('log_G2-OIII5006/Hb4861')

                lineNameList.append('ToT-OIII5006/Hb4861')
                lineNameList.append('log_ToT-OIII5006/Hb4861')

                frmList.append('f8')
                frmList.append('f8')
                frmList.append('f8')
                frmList.append('f8')
            
            if 'NII6583' in lineNameID and 'Ha6562' in lineNameID:

                NIIG2 = np.copy(lines['g2_Amp_'+'NII6583'])
                NIIG2[index] = np.nan
                NIIG2[indexSigma] = np.nan

                HalphaG2 = np.copy(lines['g2_Amp_'+'Ha6562'])
                HalphaG2[index] = np.nan
                HalphaG2[indexSigma] = np.nan

                lrNIIHaG2 = np.divide(NIIG2,HalphaG2)
                logNIIHaG2 = np.log10(lrNIIHaG2)
                lrNIIHa = np.divide((NII+NIIG2),(Halpha+HalphaG2))
                logNIIHa = np.log10(lrNIIHa)
                
                tot = np.column_stack((tot,lrNIIHaG2,logNIIHaG2,lrNIIHa,logNIIHa))
                lineNameList.append('G2-NII6583/Ha6562')
                lineNameList.append('log_G2-NII6583/Ha6562')

                lineNameList.append('ToT-NII6583/Ha6562')
                lineNameList.append('log_ToT-NII6583/Ha6562')

                frmList.append('f8')
                frmList.append('f8')
                frmList.append('f8')
                frmList.append('f8')


            if 'OI6300' in lineNameID and 'Ha6562' in lineNameID:

                OIG2 = np.copy(lines['g2_Amp_'+'OI6300'])
                OIG2[index] = np.nan
                OIG2[indexSigma] = np.nan

                HalphaG2 = np.copy(lines['g2_Amp_'+'Ha6562'])
                HalphaG2[index] = np.nan
                HalphaG2[indexSigma] = np.nan


                lrOIHaG2 = np.divide(OIG2,HalphaG2)
                logOIHaG2 = np.log10(lrOIHaG2)            

                lrOIHa = np.divide((OI+OIG2),(Halpha+HalphaG2))
                logOIHa = np.log10(lrOIHa)            
                

                tot = np.column_stack((tot,lrOIHaG2,logOIHaG2,lrOIHa,logOIHa))
                lineNameList.append('G2-OI6300/Ha6562')
                lineNameList.append('log_G2-OI6300/Ha6562')

                lineNameList.append('ToT-OI6300/Ha6562')
                lineNameList.append('log_ToT-OI6300/Ha6562')

                frmList.append('f8')
                frmList.append('f8')
                frmList.append('f8')
                frmList.append('f8')

            if 'SII6716' in lineNameID and 'Ha6562' in lineNameID:

                SII1G2 = np.copy(lines['g2_Amp_'+'SII6716'])
                SII1G2[index] = np.nan
                SII1G2[indexSigma] = np.nan
                SII2G2 = np.copy(lines['g2_Amp_'+'SII6730'])
                SII2G2[index] = np.nan
                SII2G2[indexSigma] = np.nan

                HalphaG2 = np.copy(lines['g2_Amp_'+'Ha6562'])
                HalphaG2[index] = np.nan
                HalphaG2[indexSigma] = np.nan

                lrSIIHaG2 = np.divide((SII1G2+SII2G2),HalphaG2)
                logSIIHaG2 = np.log10(lrSIIHaG2)            
                
                lrSIIHa = np.divide((SII1+SII2+SII1G2+SII1G2),(Halpha+HalphaG2))
                logSIIHa = np.log10(lrSIIHa)            
                
                tot = np.column_stack((tot,lrSIIHaG2,logSIIHaG2,lrSIIHaG2,logSIIHa))
                lineNameList.append('G2-SII6716/Ha6562')
                lineNameList.append('log_G2-SII6716/Ha6562')
    
                lineNameList.append('ToT-SII6716/Ha6562')
                lineNameList.append('log_ToT-SII6716/Ha6562')

                frmList.append('f8')
                frmList.append('f8')
                frmList.append('f8')
                frmList.append('f8')

            if modName == 'g3':
                
                if 'OIII5006' in lineNameID and 'Hb4861' in lineNameID:

                    oIIIG3 = np.copy(lines['g3_Amp_'+'OIII5006'])
                    oIIIG3[index] = np.nan
                    oIIIG3[indexSigma] = np.nan

                    hBetaG3 = lines['g3_Amp_'+'Hb4861']
                    hBetaG3[index] = np.nan
                    hBetaG3[indexSigma] = np.nan

                    lrOHbG3 = np.divide(oIIIG3,HBetaG3)
                    logOHbG3 = np.log10(lrOHbG3)            

                    lrOHb = np.divide((lines['g1_Amp_'+'OIII5006']+lines['g2_Amp_'+'OIII5006']+lines['g3_Amp_'+'OIII5006']),(lines['g1_Amp_'+'Hb4861']+
                        lines['g2_Amp_'+'Hb4861']+lines['g3_Amp_'+'Hb4861']))
                    logOHb = np.log10(lrOHb)            
                    
                    tot = np.column_stack((tot,logOHbG3,lrOHbG3,logOHb))
                    
                    lineNameList.append('G3-OIII5006/Hb4861')
                    lineNameList.append('log_G3-OIII5006/Hb4861')

                    lineNameList.append('ToT-OIII5006/Hb4861')
                    lineNameList.append('log_ToT-OIII5006/Hb4861')

                    frmList.append('f8')
                    frmList.append('f8')
                    frmList.append('f8')
                    frmList.append('f8')
                
                if 'NII6583' in lineNameID and 'Ha6562' in lineNameID:

                    NIIG3 = np.copy(lines['g3_Amp_'+'NII6583'])
                    NIIG3[index] = np.nan
                    NIIG3[indexSigma] = np.nan
                    
                    HalphaG3 = np.copy(lines['g3_Amp_'+'Ha6562'])
                    HalphaG3[index] = np.nan
                    HalphaG3[indexSigma] = np.nan

                    
                    lrNIIHaG3 = np.divide(NIIG3,HalphaG3)
                    logNIIHaG3 = np.log10(lrNIIHaG3)            

                    lrNIIHa = np.divide((NII+NIIG2+NIIG3),(Halpha+HalphaG2+HalphaG3))
                    logNIIHa = np.log10(lrNIIHa)            

                    tot = np.column_stack((tot,lrNIIHaG3,logNIIHaG3,lrNIIHa,logNIIHa))
                    lineNameList.append('G3-NII6583/Hb4861')
                    lineNameList.append('log_G3-NII6583/Hb4861')

                    lineNameList.append('ToT-NII6583/Hb4861')
                    lineNameList.append('log_ToT-NII6583/Hb4861')

                    frmList.append('f8')
                    frmList.append('f8')
                    frmList.append('f8')
                    frmList.append('f8')

                if 'OI6300' in lineNameID and 'Ha6562' in lineNameID:

                    OIG3 = np.copy(lines['g3_Amp_'+'OI6300'])
                    OIG3[index] = np.nan
                    OIG3[indexSigma] = np.nan

                    HalphaG3= np.copy(lines['g3_Amp_'+'Ha6562'])
                    HalphaG3[index] = np.nan
                    HalphaG3[indexSigma] = np.nan

                    lrOIHaG3 = np.divide(OIG3,HalphaG3)
                    logOIHaG3 = np.log10(lrOIHaG3)            

                    lrOIHa = np.divide((OI+OIG2+OIG3),(Halpha+HalphaG2+HalphaG3))
                    logOIHa = np.log10(lrOIHa)            

                    tot = np.column_stack((tot,lrOIHaG3,logOIHaG3,lrOIHa,logOIHa))
                    lineNameList.append('G3-OI6300/Ha6562')
                    lineNameList.append('log_G3-OI6300/Ha6562')

                    lineNameList.append('ToT-OI6300/Ha6562')
                    lineNameList.append('log_ToT-OI6300/Ha6562')

                    frmList.append('f8')
                    frmList.append('f8')
                    frmList.append('f8')
                    frmList.append('f8')


                if 'SII6716' in lineNameID and 'Ha6562' in lineNameID:

                    SII1G3 = np.copy(lines['g3_Amp_'+'SII6716'])
                    SII1G3[index] = np.nan
                    SII1G3[indexSigma] = np.nan
                    
                    SII2G3 = np.copy(lines['g3_Amp_'+'SII6730'])
                    SII2G3[index] = np.nan
                    SII2G3[indexSigma] = np.nan

                    HalphaG3 = lines['g3_Amp_'+'Ha6562']
                    HalphaG3[index] = np.nan
                    HalphaG3[indexSigma] = np.nan

                    lrSIIHaG3 = np.divide((SII1G3+SII2G3),HalphaG3)
                    logSIIHaG3 = np.log10(lrSIIHaG3)            

                    lrSIIHa = np.divide((SII1+SII2+SII1G2+SII2G2+SII1G2+SII2G2),
                        (Halpha+HalphaG2+HalphaG3))
                    logSIIHa = np.log10(lrSIIHa)            

                    tot = np.column_stack((tot,lrSIIHaG3,logSIIHaG3,lrSIIHa,logSIIHa))

                    lineNameList.append('G3-SII6716/Ha6562')
                    lineNameList.append('log_G3-SII6716/Ha6562')

                    lineNameList.append('ToT-SII6716/Ha6562')
                    lineNameList.append('log_ToT-SII6716/Ha6562')

                    frmList.append('f8')
                    frmList.append('f8')
                    frmList.append('f8')
                    frmList.append('f8')


        t = Table(tot, names=(lineNameList))

        try:
            tb = Table(hdul['LineRatios_'+modName].data)
            hdul['LineRatios_'+modName] = fits.BinTableHDU(t.as_array(),name='LineRatios_'+modName)
        except KeyError as e:
            tt=fits.BinTableHDU.from_columns(t.as_array(),name='LineRatios_'+modName)   
            hdul.append(tt) 



        #hdul.append(fits.BinTableHDU(t.as_array(), name='LineRatios_'+modName))

        indexSFK = np.where(np.logical_and(np.logical_and(t['log_G1-OIII5006/Hb4861'] < 0.61 / (t['log_G1-NII6583/Ha6562'] - 0.05) + 1.3,
            t['log_G1-OIII5006/Hb4861']<3),
            t['log_G1-NII6583/Ha6562']<-0.05)
            )


        indexSF = np.where(np.logical_and(np.logical_and(np.logical_and(t['log_G1-OIII5006/Hb4861'] < 0.61 / (t['log_G1-NII6583/Ha6562'] - 0.47) + 1.19, 
            t['log_G1-OIII5006/Hb4861'] >= 0.61 / (t['log_G1-NII6583/Ha6562'] - 0.05) + 1.3),
            t['log_G1-OIII5006/Hb4861']<3.),
            t['log_G1-NII6583/Ha6562']<0.47))
        
        indexAGN  = np.where(np.logical_and(np.logical_and(t['log_G1-OIII5006/Hb4861'] >= 0.61 / (t['log_G1-NII6583/Ha6562'] - 0.47) + 1.19,
            t['log_G1-OIII5006/Hb4861']<3.),
            t['log_G1-NII6583/Ha6562']<=0.5))
        

        indexBadFit = np.where(np.logical_or.reduce((t['log_G1-OIII5006/Hb4861']>3,t['log_G1-OIII5006/Hb4861']<-3,t['log_G1-NII6583/Ha6562']>0.8,
            t['log_G1-NII6583/Ha6562']<-2.)))
        #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-OIII5006/Hb4861']),np.isnan(t['G1-NII6583/Ha6562'])))
        

        LrOIII  = np.zeros(len(lines['BIN_ID']))*np.nan
        LrOIII[indexSFK] = 0.
        LrOIII[indexSF] = 1.
        LrOIII[indexAGN] = 2.

        LrOIII[indexBadFit] = -1.
        #LrOIII[indexBadFit2] = -2.


        indexSF = np.where(np.logical_and.reduce((t['log_G1-OIII5006/Hb4861'] < 0.72 / (t['log_G1-SII6716/Ha6562'] - 0.32) + 1.30,
                            t['log_G1-OIII5006/Hb4861']<3.,
                            t['log_G1-SII6716/Ha6562']<0.32,
                            t['log_G1-OIII5006/Hb4861']>=-2.,
                            t['log_G1-SII6716/Ha6562']>=-2.)))
        indexSey = np.where(np.logical_and.reduce((t['log_G1-OIII5006/Hb4861'] >= 0.72 / (t['log_G1-SII6716/Ha6562'] - 0.32) + 1.30, 
            t['log_G1-OIII5006/Hb4861'] > 1.89* t['log_G1-SII6716/Ha6562'] + 0.76,
                            t['log_G1-OIII5006/Hb4861']<3.,
                            t['log_G1-SII6716/Ha6562']<1.,
                            t['log_G1-OIII5006/Hb4861']>=-2.,
                            t['log_G1-SII6716/Ha6562']>=-2.)))
        
        indexLIN = np.where(np.logical_and.reduce((t['log_G1-OIII5006/Hb4861'] >= 0.72 / (t['log_G1-SII6716/Ha6562'] - 0.32) + 1.30,
            t['log_G1-OIII5006/Hb4861'] < 1.89*t['log_G1-SII6716/Ha6562'] + 0.76,
                            t['log_G1-OIII5006/Hb4861']<3.,
                            t['log_G1-SII6716/Ha6562']<1.,
                            t['log_G1-OIII5006/Hb4861']>=-2.,
                            t['log_G1-SII6716/Ha6562']>=-2.)))

        indexBadFit = np.where(np.logical_or.reduce((t['log_G1-OIII5006/Hb4861']<=-2.,t['log_G1-OIII5006/Hb4861']>3.,
            t['log_G1-SII6716/Ha6562']<=-2.,t['log_G1-SII6716/Ha6562']>0.5)))
        #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-SII6716/Ha6562']),np.isnan(t['G1-NII6583/Ha6562'])))

        LrSII  = np.zeros(len(lines['BIN_ID']))*np.nan
        LrSII[indexSF] = 0.
        LrSII[indexSey] = 1.
        LrSII[indexLIN] = 2.
        LrSII[indexBadFit] = -1.
        #LrSII[indexBadFit2] = -2.

        indexSF = np.where(np.logical_and.reduce((t['log_G1-OIII5006/Hb4861'] < 0.73 / (t['log_G1-OI6300/Ha6562'] + 0.59) + 1.33,
            t['log_G1-OIII5006/Hb4861']<3.,
            t['log_G1-OI6300/Ha6562']<-0.59,
            t['log_G1-OIII5006/Hb4861']>=-2.,            
            t['log_G1-OI6300/Ha6562']>=-3.)))
        
        indexSey = np.where(np.logical_and.reduce((t['log_G1-OIII5006/Hb4861'] >= 0.73 / (t['log_G1-OI6300/Ha6562'] + 0.59) +1.33,
            t['log_G1-OIII5006/Hb4861'] >= 1.18* t['log_G1-OI6300/Ha6562'] + 1.30,
            t['log_G1-OIII5006/Hb4861']<3.,
            t['log_G1-OI6300/Ha6562']<2.,
            t['log_G1-OIII5006/Hb4861']>=-2.,            
            t['log_G1-OI6300/Ha6562']>=-3.)))
        
        indexLIN = np.where(np.logical_and.reduce((t['log_G1-OIII5006/Hb4861'] >= 0.73 / (t['log_G1-OI6300/Ha6562'] + 0.59)+1.33, 
            t['log_G1-OIII5006/Hb4861'] < 1.18* t['log_G1-OI6300/Ha6562'] + 1.30,
            t['log_G1-OIII5006/Hb4861']<3.,
            t['log_G1-OI6300/Ha6562']<2.,
            t['log_G1-OIII5006/Hb4861']>=-2.,            
            t['log_G1-OI6300/Ha6562']>=-3.)))

        indexBadFit = np.where(np.logical_or.reduce((t['log_G1-OIII5006/Hb4861']<-2.,t['log_G1-OIII5006/Hb4861']>3.,
            t['log_G1-OI6300/Ha6562']<-3.,t['log_G1-OI6300/Ha6562']>0.)))
        
        #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-OI6300/Ha6562']),np.isnan(t['G1-OIII5006/Hb4861'])))

        LrOI  = np.zeros(len(lines['BIN_ID']))*np.nan
        LrOI[indexSF] = 0.
        LrOI[indexSey] = 1.
        LrOI[indexLIN] = 2.
        LrOI[indexBadFit] = -1.
       # LrOI[indexBadFit2] = -2.

        tt=Table([lines['BIN_ID'],LrOIII,LrSII,LrOI],names=('BIN_ID','G1-BPT_OIII','G1-BPT_SII','G1-BPT_OI'))

        if modName != 'g1':

            indexSFK = np.where(np.logical_and(np.logical_and(t['log_G2-OIII5006/Hb4861'] < 0.61 / (t['log_G2-NII6583/Ha6562'] - 0.05) + 1.3,
                t['log_G2-OIII5006/Hb4861']<3),
                t['log_G2-NII6583/Ha6562']<-0.05)
                )

            indexSF = np.where(np.logical_and(np.logical_and(np.logical_and(t['log_G2-OIII5006/Hb4861'] < 0.61 / (t['log_G2-NII6583/Ha6562'] - 0.47) + 1.19, 
                t['log_G2-OIII5006/Hb4861'] >= 0.61 / (t['log_G2-NII6583/Ha6562'] - 0.05) + 1.3),
                t['log_G2-OIII5006/Hb4861']<3.),
                t['log_G2-NII6583/Ha6562']<0.47))
            
            indexAGN  = np.where(np.logical_and(np.logical_and(t['log_G2-OIII5006/Hb4861'] >= 0.61 / (t['log_G2-NII6583/Ha6562'] - 0.47) + 1.19,
                t['log_G2-OIII5006/Hb4861']<3.),
                t['log_G2-NII6583/Ha6562']<=0.5))
            

            indexBadFit = np.where(np.logical_or.reduce((t['log_G2-OIII5006/Hb4861']>3,t['log_G2-OIII5006/Hb4861']<-3,t['log_G2-NII6583/Ha6562']>0.8,
                t['log_G2-NII6583/Ha6562']<-2.)))
            #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-OIII5006/Hb4861']),np.isnan(t['G1-NII6583/Ha6562'])))
            

            LrOIII  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrOIII[indexSFK] = 0.
            LrOIII[indexSF] = 1.
            LrOIII[indexAGN] = 2.

            LrOIII[indexBadFit] = -1.
            #LrOIII[indexBadFit2] = -2.

            indexSF = np.where(np.logical_and.reduce((t['log_G2-OIII5006/Hb4861'] < 0.72 / (t['log_G2-SII6716/Ha6562'] - 0.32) + 1.30,
                                t['log_G2-OIII5006/Hb4861']<3.,
                                t['log_G2-SII6716/Ha6562']<0.32,
                                t['log_G2-OIII5006/Hb4861']>=-2.,
                                t['log_G2-SII6716/Ha6562']>=-2.)))
            indexSey = np.where(np.logical_and.reduce((t['log_G2-OIII5006/Hb4861'] >= 0.72 / (t['log_G2-SII6716/Ha6562'] - 0.32) + 1.30, 
                t['log_G2-OIII5006/Hb4861'] > 1.89* t['log_G2-SII6716/Ha6562'] + 0.76,
                                t['log_G2-OIII5006/Hb4861']<3.,
                                t['log_G2-SII6716/Ha6562']<1.,
                                t['log_G2-OIII5006/Hb4861']>=-2.,
                                t['log_G2-SII6716/Ha6562']>=-2.)))
            
            indexLIN = np.where(np.logical_and.reduce((t['log_G2-OIII5006/Hb4861'] >= 0.72 / (t['log_G2-SII6716/Ha6562'] - 0.32) + 1.30,
                t['log_G2-OIII5006/Hb4861'] < 1.89*t['log_G2-SII6716/Ha6562'] + 0.76,
                                t['log_G2-OIII5006/Hb4861']<3.,
                                t['log_G2-SII6716/Ha6562']<1.,
                                t['log_G2-OIII5006/Hb4861']>=-2.,
                                t['log_G2-SII6716/Ha6562']>=-2.)))

            #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-SII6716/Ha6562']),np.isnan(t['G1-NII6583/Ha6562'])))
            indexBadFit = np.where(np.logical_or.reduce((t['log_G2-OIII5006/Hb4861']<=-2.,t['log_G2-OIII5006/Hb4861']>3.,
                t['log_G2-SII6716/Ha6562']<=-2.,t['log_G2-SII6716/Ha6562']>0.5)))



            LrSII  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrSII[indexSF] = 0.
            LrSII[indexSey] = 1.
            LrSII[indexLIN] = 2.
            LrSII[indexBadFit] = -1.
            #LrSII[indexBadFit2] = -2.

            indexSF = np.where(np.logical_and.reduce((t['log_G2-OIII5006/Hb4861'] < 0.73 / (t['log_G2-OI6300/Ha6562'] + 0.59) + 1.33,
                t['log_G2-OIII5006/Hb4861']<3.,
                t['log_G2-OI6300/Ha6562']<-0.59,
                t['log_G2-OIII5006/Hb4861']>=-2.,            
                t['log_G2-OI6300/Ha6562']>=-3.)))
            
            indexSey = np.where(np.logical_and.reduce((t['log_G2-OIII5006/Hb4861'] >= 0.73 / (t['log_G2-OI6300/Ha6562'] + 0.59) +1.33,
                t['log_G2-OIII5006/Hb4861'] >= 1.18* t['log_G2-OI6300/Ha6562'] + 1.30,
                t['log_G2-OIII5006/Hb4861']<3.,
                t['log_G2-OI6300/Ha6562']<2.,
                t['log_G2-OIII5006/Hb4861']>=-2.,            
                t['log_G2-OI6300/Ha6562']>=-3.)))
            
            indexLIN = np.where(np.logical_and.reduce((t['log_G2-OIII5006/Hb4861'] >= 0.73 / (t['log_G2-OI6300/Ha6562'] + 0.59)+1.33, 
                t['log_G2-OIII5006/Hb4861'] < 1.18* t['log_G2-OI6300/Ha6562'] + 1.30,
                t['log_G2-OIII5006/Hb4861']<3.,
                t['log_G2-OI6300/Ha6562']<2.,
                t['log_G2-OIII5006/Hb4861']>=-2.,            
                t['log_G2-OI6300/Ha6562']>=-3.)))

            indexBadFit = np.where(np.logical_or.reduce((t['log_G2-OIII5006/Hb4861']<-2.,t['log_G2-OIII5006/Hb4861']>3.,
                t['log_G2-OI6300/Ha6562']<-3.,t['log_G2-OI6300/Ha6562']>0.)))
            
            #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-OI6300/Ha6562']),np.isnan(t['G1-OIII5006/Hb4861'])))

            LrOI  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrOI[indexSF] = 0.
            LrOI[indexSey] = 1.
            LrOI[indexLIN] = 2.
            LrOI[indexBadFit] = -1.

            tt.add_column(Column(LrOIII,name='G2-BPT_OIII'))
            tt.add_column(Column(LrSII,name='G2-BPT_SII'))
            tt.add_column(Column(LrOI,name='G2-BPT_OI'))

            indexSFK = np.where(np.logical_and(np.logical_and(t['log_ToT-OIII5006/Hb4861'] < 0.61 / (t['log_ToT-NII6583/Ha6562'] - 0.05) + 1.3,
                t['log_ToT-OIII5006/Hb4861']<3),
                t['log_ToT-NII6583/Ha6562']<-0.05)
                )

            indexSF = np.where(np.logical_and(np.logical_and(np.logical_and(t['log_ToT-OIII5006/Hb4861'] < 0.61 / (t['log_ToT-NII6583/Ha6562'] - 0.47) + 1.19, 
                t['log_ToT-OIII5006/Hb4861'] >= 0.61 / (t['log_ToT-NII6583/Ha6562'] - 0.05) + 1.3),
                t['log_ToT-OIII5006/Hb4861']<3.),
                t['log_ToT-NII6583/Ha6562']<0.47))
            
            indexAGN  = np.where(np.logical_and(np.logical_and(t['log_ToT-OIII5006/Hb4861'] >= 0.61 / (t['log_ToT-NII6583/Ha6562'] - 0.47) + 1.19,
                t['log_ToT-OIII5006/Hb4861']<3.),
                t['log_ToT-NII6583/Ha6562']<=0.5))
            

            indexBadFit = np.where(np.logical_or.reduce((t['log_ToT-OIII5006/Hb4861']>3,t['log_ToT-OIII5006/Hb4861']<-3,t['log_ToT-NII6583/Ha6562']>0.8,
                t['log_ToT-NII6583/Ha6562']<-2.)))
            #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-OIII5006/Hb4861']),np.isnan(t['G1-NII6583/Ha6562'])))
            

            LrOIII  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrOIII[indexSFK] = 0.
            LrOIII[indexSF] = 1.
            LrOIII[indexAGN] = 2.

            LrOIII[indexBadFit] = -1.
            #LrOIII[indexBadFit2] = -2.


            indexSF = np.where(np.logical_and.reduce((t['log_ToT-OIII5006/Hb4861'] < 0.72 / (t['log_ToT-SII6716/Ha6562'] - 0.32) + 1.30,
                                t['log_ToT-OIII5006/Hb4861']<3.,
                                t['log_ToT-SII6716/Ha6562']<0.32,
                                t['log_ToT-OIII5006/Hb4861']>=-2.,
                                t['log_ToT-SII6716/Ha6562']>=-2.)))
            indexSey = np.where(np.logical_and.reduce((t['log_ToT-OIII5006/Hb4861'] >= 0.72 / (t['log_ToT-SII6716/Ha6562'] - 0.32) + 1.30, 
                t['log_ToT-OIII5006/Hb4861'] > 1.89* t['log_ToT-SII6716/Ha6562'] + 0.76,
                                t['log_ToT-OIII5006/Hb4861']<3.,
                                t['log_ToT-SII6716/Ha6562']<1.,
                                t['log_ToT-OIII5006/Hb4861']>=-2.,
                                t['log_ToT-SII6716/Ha6562']>=-2.)))
            
            indexLIN = np.where(np.logical_and.reduce((t['log_ToT-OIII5006/Hb4861'] >= 0.72 / (t['log_ToT-SII6716/Ha6562'] - 0.32) + 1.30,
                t['log_ToT-OIII5006/Hb4861'] < 1.89*t['log_ToT-SII6716/Ha6562'] + 0.76,
                                t['log_ToT-OIII5006/Hb4861']<3.,
                                t['log_ToT-SII6716/Ha6562']<1.,
                                t['log_ToT-OIII5006/Hb4861']>=-2.,
                                t['log_ToT-SII6716/Ha6562']>=-2.)))

            #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-SII6716/Ha6562']),np.isnan(t['G1-NII6583/Ha6562'])))
            indexBadFit = np.where(np.logical_or.reduce((t['log_ToT-OIII5006/Hb4861']<=-2.,t['log_ToT-OIII5006/Hb4861']>3,
                t['log_ToT-NII6583/Ha6562']<=2.,t['log_ToT-NII6583/Ha6562']>0.5)))



            LrSII  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrSII[indexSF] = 0.
            LrSII[indexSey] = 1.
            LrSII[indexLIN] = 2.
            LrSII[indexBadFit] = -1.
            #LrSII[indexBadFit2] = -2.

            indexSF = np.where(np.logical_and.reduce((t['log_ToT-OIII5006/Hb4861'] < 0.73 / (t['log_ToT-OI6300/Ha6562'] + 0.59) + 1.33,
                t['log_ToT-OIII5006/Hb4861']<3.,
                t['log_ToT-OI6300/Ha6562']<-0.59,
                t['log_ToT-OIII5006/Hb4861']>=-2.,            
                t['log_ToT-OI6300/Ha6562']>=-3.)))
            
            indexSey = np.where(np.logical_and.reduce((t['log_ToT-OIII5006/Hb4861'] >= 0.73 / (t['log_ToT-OI6300/Ha6562'] + 0.59) +1.33,
                t['log_ToT-OIII5006/Hb4861'] >= 1.18* t['log_ToT-OI6300/Ha6562'] + 1.30,
                t['log_ToT-OIII5006/Hb4861']<3.,
                t['log_ToT-OI6300/Ha6562']<2.,
                t['log_ToT-OIII5006/Hb4861']>=-2.,            
                t['log_ToT-OI6300/Ha6562']>=-3.)))
            
            indexLIN = np.where(np.logical_and.reduce((t['log_ToT-OIII5006/Hb4861'] >= 0.73 / (t['log_ToT-OI6300/Ha6562'] + 0.59)+1.33, 
                t['log_ToT-OIII5006/Hb4861'] < 1.18* t['log_ToT-OI6300/Ha6562'] + 1.30,
                t['log_ToT-OIII5006/Hb4861']<3.,
                t['log_ToT-OI6300/Ha6562']<2.,
                t['log_ToT-OIII5006/Hb4861']>=-2.,            
                t['log_ToT-OI6300/Ha6562']>=-3.)))

            indexBadFit = np.where(np.logical_or.reduce((t['log_ToT-OIII5006/Hb4861']<-2.,t['log_ToT-OIII5006/Hb4861']>3.,
                t['log_ToT-OI6300/Ha6562']<-3.,t['log_ToT-OI6300/Ha6562']>0.)))
            
            #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-OI6300/Ha6562']),np.isnan(t['G1-OIII5006/Hb4861'])))

            LrOI  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrOI[indexSF] = 0.
            LrOI[indexSey] = 1.
            LrOI[indexLIN] = 2.
            LrOI[indexBadFit] = -1.

            tt.add_column(Column(LrOIII,name='ToT-BPT_OIII'))
            tt.add_column(Column(LrSII,name='ToT-BPT_SII'))
            tt.add_column(Column(LrOI,name='ToT-BPT_OI'))


            if modName == 'g3':

                indexSFK = np.where(np.logical_and(np.logical_and(t['log_G3-OIII5006/Hb4861'] < 0.61 / (t['log_G3-NII6583/Ha6562'] - 0.05) + 1.3,
                    t['log_G3-OIII5006/Hb4861']<3),
                    t['log_G3-NII6583/Ha6562']<-0.05)
                    )

                indexSF = np.where(np.logical_and(np.logical_and(np.logical_and(t['log_G3-OIII5006/Hb4861'] < 0.61 / (t['log_G3-NII6583/Ha6562'] - 0.47) + 1.19, 
                    t['log_G3-OIII5006/Hb4861'] >= 0.61 / (t['log_G3-NII6583/Ha6562'] - 0.05) + 1.3),
                    t['log_G3-OIII5006/Hb4861']<3.),
                    t['log_G3-NII6583/Ha6562']<0.47))
                
                indexAGN  = np.where(np.logical_and(np.logical_and(t['log_G3-OIII5006/Hb4861'] >= 0.61 / (t['log_G3-NII6583/Ha6562'] - 0.47) + 1.19,
                    t['log_G3-OIII5006/Hb4861']<3.),
                    t['log_G3-NII6583/Ha6562']<=0.5))
                

                indexBadFit = np.where(np.logical_or.reduce((t['log_G3-OIII5006/Hb4861']>3,t['log_G3-OIII5006/Hb4861']<-3,t['log_G3-NII6583/Ha6562']>0.8,
                    t['log_G3-NII6583/Ha6562']<-2.)))
                #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-OIII5006/Hb4861']),np.isnan(t['G1-NII6583/Ha6562'])))
                

                LrOIII  = np.zeros(len(lines['BIN_ID']))*np.nan
                LrOIII[indexSFK] = 0.
                LrOIII[indexSF] = 1.
                LrOIII[indexAGN] = 2.

                LrOIII[indexBadFit] = -1.
                #LrOIII[indexBadFit2] = -2.

                indexSF = np.where(np.logical_and.reduce((t['log_G3-OIII5006/Hb4861'] < 0.72 / (t['log_G3-SII6716/Ha6562'] - 0.32) + 1.30,
                                    t['log_G3-OIII5006/Hb4861']<3.,
                                    t['log_G3-SII6716/Ha6562']<0.32,
                                    t['log_G3-OIII5006/Hb4861']>=-2.,
                                    t['log_G3-SII6716/Ha6562']>=-2.)))
                indexSey = np.where(np.logical_and.reduce((t['log_G3-OIII5006/Hb4861'] >= 0.72 / (t['log_G3-SII6716/Ha6562'] - 0.32) + 1.30, 
                    t['log_G3-OIII5006/Hb4861'] > 1.89* t['log_G3-SII6716/Ha6562'] + 0.76,
                                    t['log_G3-OIII5006/Hb4861']<3.,
                                    t['log_G3-SII6716/Ha6562']<1.,
                                    t['log_G3-OIII5006/Hb4861']>=-2.,
                                    t['log_G3-SII6716/Ha6562']>=-2.)))
                
                indexLIN = np.where(np.logical_and.reduce((t['log_G3-OIII5006/Hb4861'] >= 0.72 / (t['log_G3-SII6716/Ha6562'] - 0.32) + 1.30,
                    t['log_G3-OIII5006/Hb4861'] < 1.89*t['log_G3-SII6716/Ha6562'] + 0.76,
                                    t['log_G3-OIII5006/Hb4861']<3.,
                                    t['log_G3-SII6716/Ha6562']<1.,
                                    t['log_G3-OIII5006/Hb4861']>=-2.,
                                    t['log_G3-SII6716/Ha6562']>=-2.)))

                indexBadFit = np.where(np.logical_or.reduce(t['log_G3-OIII5006/Hb4861']<=-2.,t['log_G3-OIII5006/Hb4861']>3.,
                    t['log_G3-SII6716/Ha6562']<=-2.,t['log_G3-SII6716/Ha6562']>0.5))
                #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-SII6716/Ha6562']),np.isnan(t['G1-NII6583/Ha6562'])))

                LrSII  = np.zeros(len(lines['BIN_ID']))*np.nan
                LrSII[indexSF] = 0.
                LrSII[indexSey] = 1.
                LrSII[indexLIN] = 2.
                LrSII[indexBadFit] = -1.
                #LrSII[indexBadFit2] = -2.

                indexSF = np.where(np.logical_and.reduce((t['log_G3-OIII5006/Hb4861'] < 0.73 / (t['log_G3-OI6300/Ha6562'] + 0.59) + 1.33,
                    t['log_G3-OIII5006/Hb4861']<3.,
                    t['log_G3-OI6300/Ha6562']<-0.59,
                    t['log_G3-OIII5006/Hb4861']>=-2.,            
                    t['log_G3-OI6300/Ha6562']>=-3.)))
                
                indexSey = np.where(np.logical_and.reduce((t['log_G3-OIII5006/Hb4861'] >= 0.73 / (t['log_G3-OI6300/Ha6562'] + 0.59) +1.33,
                    t['log_G3-OIII5006/Hb4861'] >= 1.18* t['log_G3-OI6300/Ha6562'] + 1.30,
                    t['log_G3-OIII5006/Hb4861']<3.,
                    t['log_G3-OI6300/Ha6562']<2.,
                    t['log_G3-OIII5006/Hb4861']>=-2.,            
                    t['log_G3-OI6300/Ha6562']>=-3.)))
                
                indexLIN = np.where(np.logical_and.reduce((t['log_G3-OIII5006/Hb4861'] >= 0.73 / (t['log_G3-OI6300/Ha6562'] + 0.59)+1.33, 
                    t['log_G3-OIII5006/Hb4861'] < 1.18* t['log_G3-OI6300/Ha6562'] + 1.30,
                    t['log_G3-OIII5006/Hb4861']<3.,
                    t['log_G3-OI6300/Ha6562']<2.,
                    t['log_G3-OIII5006/Hb4861']>=-2.,            
                    t['log_G3-OI6300/Ha6562']>=-3.)))

                indexBadFit = np.where(np.logical_or.reduce((t['log_G3-OIII5006/Hb4861']<-2.,t['log_G3-OIII5006/Hb4861']>3.,
                    t['log_G3-OI6300/Ha6562']<-3.,t['log_G3-OI6300/Ha6562']>0.)))
                
                #indexBadFit2 = np.where(np.logical_or(np.isnan(t['G1-OI6300/Ha6562']),np.isnan(t['G1-OIII5006/Hb4861'])))

                LrOI  = np.zeros(len(lines['BIN_ID']))*np.nan
                LrOI[indexSF] = 0.
                LrOI[indexSey] = 1.
                LrOI[indexLIN] = 2.
                LrOI[indexBadFit] = -1.

                tt.add_column(Column(LrOIII,name='G3-BPT_OIII'))
                tt.add_column(Column(LrSII,name='G3-BPT_SII'))
                tt.add_column(Column(LrOI,name='G3-BPT_OI'))


        try:
            tt = Table(hdul['BPT_'+modName].data)
            hdul['BPT_'+modName] = fits.BinTableHDU(tt.as_array(),name='BPT_'+modName)

        except KeyError as e:
            tt=fits.BinTableHDU.from_columns(tt.as_array(),name='BPT_'+modName)   
            hdul.append(tt)          
        
        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)

        #hdul.append(fits.BinTableHDU(tt.as_array(), name='BPT_'+modName))
        #hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)

        return



    def saveAncelsTable(self,cfg_par, sigmaCenArr):
        modName = cfg_par['gFit']['modName']
 
        hdul = fits.open(cfg_par['general']['outTableName'])

        try:
            tt = Table(hdul['Ancels'+modName].data)
            hdul['Ancels'+modName] = fits.BinTableHDU.from_columns(sigmaCenArr,name='Ancels'+modName)

        except KeyError as e:
            tt=fits.BinTableHDU.from_columns(sigmaCenArr,name='Ancels'+modName)   
            hdul.append(tt)  
        
        
        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)


        return

    def cardan(self,a,b,c,d):
        J=np.exp(2j*np.pi/3)
        Jc=1/J
        u=np.empty(2,np.complex128)
        z0=b/3/a
        a2,b2 = a*a,b*b    
        p=-b2/3/a2 +c/a
        q=(b/27*(2*b2/a2-9*c/a)+d)/a
        D=-4*p*p*p-27*q*q
        r=np.sqrt(-D/27+0j)        
        u=((-q-r)/2)**0.33333333333333333333333
        v=((-q+r)/2)**0.33333333333333333333333
        w=u*v
        w0=np.abs(w+p/3)
        w1=np.abs(w*J+p/3)
        w2=np.abs(w*Jc+p/3)
        if w0<w1: 
            if w2<w0 : v*=Jc
        elif w2<w1 : v*=Jc
        else: v*=J        
        return u+v-z0, u*J+v*Jc-z0,u*Jc+v*J-z0

    def findferrari(self,a,b,c,d,e):
        "resolution of P=ax^4+bx^3+cx^2+dx+e=0"
        "CN all coeffs real."
        "First shift : x= z-b/4/a  =>  P=z^4+pz^2+qz+r"
        z0=b/4./a
        a2,b2,c2,d2 = a*a,b*b,c*c,d*d 
        p = -3.*b2/(8.*a2)+c/a
        q = b*b2/8/a/a2 - 1./2.*b*c/a2 + d/a
        r = -3./256.*b2*b2/a2/a2 +c*b2/a2/a/16.-b*d/a2/4.+e/a
        "Second find X so P2=AX^3+BX^2+C^X+D=0"
        A=8.
        B=-4.*p
        C=-8.*r
        D=4.*r*p-q*q
        y0,y1,y2=self.cardan(A,B,C,D)
        if np.abs(y1.imag)<np.abs(y0.imag): y0=y1 
        if np.abs(y2.imag)<np.abs(y0.imag): y0=y2 
        a0=(-p+2*y0.real)**.5
        if a0==0 : b0=y0**2.-r
        else : b0=-q/2./a0
        r0,r1=self.roots2(1.,a0,y0+b0)
        r2,r3=self.roots2(1.,-a0,y0-b0)
        return r0-z0,r1-z0,r2-z0,r3-z0

    def roots2(self,a,b,c):
        bp=b/2    
        delta=bp*bp-a*c
        u1=(-bp-delta**.5)/a
        u2=-u1-b/a
        return u1,u2  

    def carolloDistOIII(self,cfg_par):
        
        modName = cfg_par['gFit']['modName']

        hdul = fits.open(cfg_par['general']['outTableName'])

        if cfg_par['gFit']['modName'] == 'g1':
            modString = ['G1']
            gNameTable = ['G1']
            gName = ['G1']
        elif cfg_par['gFit']['modName'] == 'g2':
            modString = ['G1','G2','ToT']
            gNameTable = ['G2','G2','G2']
            gName = ['G1','G2','ToT']
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']
            gNameTable = ['G3','G3','G3','G3']
            gName = ['G1','G2','G3','ToT']

        for j in range(0,len(modString)):

            bptInfo = hdul['BPT_'+gNameTable[j]].data
            lineRatio = hdul['LINERATIOS_'+gNameTable[j]].data 
            LrOIII  = bptInfo[gName[j]+'-BPT_OIII']
            xArr = np.log10(lineRatio[gName[j]+'-NII6583/Ha6562'])
            yArr = np.log10(lineRatio[gName[j]+'-OIII5006/Hb4861'])

            nu = np.zeros(len(bptInfo['BIN_ID']))*np.nan

            for i in range(0,len(xArr)):
                x1 = xArr[i]
                y1 = yArr[i]
                if np.isnan(LrOIII[i]):
                    continue

                if LrOIII[i] == 0:
                    a = 40000
                    b = -40000*x1-6000
                    c = 6000*x1+300
                    d = -300.*x1+24400*y1-31725
                    e = -1220.*y1+5.*x1-13298.
                    
                    sols = self.findferrari(a,b,c,d,e)

                    if len(sols) == 0 : 
                        nu[i] = np.nan
                        continue

                    d1=[]
                    for xS in sols:
                        dist = np.sqrt(np.power(x1-xS.real,2)+np.power(y1-(np.divide(0.61,xS.real-0.05)+1.3),2))
                        d1.append(dist)

                    if len(d1) !=0: 
                        d1=np.min(d1)
                    else:
                        nu[i] = np.nan
                        continue

                    nu[i] = -0.5 - d1


                elif LrOIII[i] == 2:

                    a = 1e6
                    b = -1e6*x1-1.41e6
                    c = 1.41e6*x1 + 6.627e5
                    d = -6.627e5*x1+6.1e5*y1-8.29723e5
                    e = -2.867e5*y1+1.03823e5*x1-3.0927e4            
                    
                    sols = self.findferrari(a,b,c,d,e)
                    if len(sols) == 0 : 
                        nu[i] = np.nan
                        continue

                    d1=[]
                    for xS in sols:
                        dist = np.sqrt(np.power(x1-xS.real,2)+np.power(y1-(np.divide(0.61,xS.real-0.47)+1.19),2))
                        d1.append(dist)
     
                    if len(d1) !=0: 
                        d1=np.min(d1)
                    else:
                        nu[i] = np.nan
                        continue           

                    nu[i] = 0.5 + d1

                elif LrOIII[i] == 1:

                    aSF = 40000.
                    bSF = -40000.*x1-6000.
                    cSF = 6000.*x1+300.
                    dSF = -300.*x1+24400.*y1-31725
                    eSF = -1220.*y1+5.*x1-13298.

                    solsSF = self.findferrari(aSF,bSF,cSF,dSF,eSF)

                    if len(solsSF) == 0 : 
                        nu[i] = np.nan
                        continue

                    d1=[]
                    for xS in solsSF:
                        dist = np.sqrt(np.power(x1-xS.real,2)+np.power(y1-(np.divide(0.61,xS.real-0.05)+1.3),2))
                        d1.append(dist)

                    if len(d1) !=0: 
                         d1=np.min(d1)
                    else:
                         nu[i] = np.nan
                         continue

                    aAG = 1e6
                    bAG = -1e6*x1-1.41e6
                    cAG = 1.41e6*x1 + 6.627e5
                    dAG = -6.627e5*x1+6.1e5*y1-8.29723e5
                    eAG = -2.867e5*y1+1.03823e5*x1-3.0927e4
                    
                    solsAG = self.findferrari(aAG,bAG,cAG,dAG,eAG)
                    
                    if len(solsAG) == 0 : 
                        nu[i] = np.nan
                        continue
                    
                    d2=[]
                    for xS in solsSF:
                        dist = np.sqrt(np.power(x1-xS.real,2)+np.power(y1-(np.divide(0.61,xS.real-0.05)+1.3),2))
                        d2.append(dist)

                    if len(d2) !=0: 
                         d2=np.min(d2)
                    else:
                         nu[i] = np.nan
                         continue

                    nu[i] = -0.5+d1

            t2 = Table(bptInfo)
     
            if 'cDist-OIII'+modString[j] not in bptInfo.dtype.names: 
                t2.add_column(Column(nu,name='cDist-OIII'+modString[j]))
            else:
                t2.replace_column('cDist-OIII'+modString[j],Column(nu,name='cDist-OIII'+modString[j]))
    
            try:
                tt = Table(hdul['BPT_'+modName].data)
                hdul['BPT_'+modName] = fits.BinTableHDU(t2.as_array(),name='BPT_'+modName)
            except KeyError as e:
                tt=fits.BinTableHDU(t2.as_array(),name='LineRes_'+modName)   
                hdul.append(tt)          
        
        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)
        return



