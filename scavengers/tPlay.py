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
        dltSigmaIn1Ang = np.zeros([lenTable])

        dltV12Ang = np.zeros([lenTable])
        dltSigma12Ang = np.zeros([lenTable])

        dltV13Ang = np.zeros([lenTable])
        dltSigma13Ang = np.zeros([lenTable])
        lineRange = np.zeros([lenTable])
        cenRange = np.zeros([lenTable])

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


            deltaSigmaIn1 = np.log(cvP.vRadLambda(cfg_par['gFit']['sigmaIn1'],
                lambdaRest))
            deltaSigmaIn1 -= np.log(lambdaRest)            
            deltaSigma12 = np.log(cvP.vRadLambda(cfg_par['gFit']['dltSigma12'],
                lambdaRest))
            deltaSigma12 -=  np.log(lambdaRest)
            deltaSigma13 = np.log(cvP.vRadLambda(cfg_par['gFit']['dltSigma13'],
                lambdaRest))
            deltaSigma13 -= np.log(lambdaRest)


            dltSigmaIn1Ang[i] = deltaSigmaIn1
            dltV12Ang[i] = deltaV12            
            dltSigma12Ang[i] = deltaSigma12
            dltV13Ang[i] = deltaV13
            dltSigma13Ang[i] = deltaSigma13

        dltSigmaIn1Col = Column(name='deltaSigmaAng_In1', data=dltSigmaIn1Ang)        
        dltV12Col = Column(name='deltaVAng_12', data=dltV12Ang)
        dltSigma12Col = Column(name='deltaSigmaAng_12', data=dltSigma12Ang)
        dltV13Col = Column(name='deltaVAng_13', data=dltV13Ang)
        dltSigma13Col = Column(name='deltaSigmaAng_13', data=dltSigma13Ang)
        lineRangeCol = Column(name='lineRangeAng', data=lineRange)
        cenRangeCol = Column(name='cenRangeAng', data=cenRange)


        lineInfo.add_column(dltSigmaIn1Col)
        lineInfo.add_column(dltV12Col)
        lineInfo.add_column(dltSigma12Col)
        lineInfo.add_column(dltV13Col)
        lineInfo.add_column(dltSigma13Col)
        lineInfo.add_column(lineRangeCol)
        lineInfo.add_column(cenRangeCol)


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

        xAxis = (np.linspace(1, shapeX+1, shapeX+1)-head['CRPIX1']+1) *head['PIXSIZE']
        yAxis = (np.linspace(1, shapeY+1, shapeY+1)-head['CRPIX2']+1) *head['PIXSIZE']
        tab = fits.open(tableSpec)
        dataSpec = tab[1].data
        specExp = tab[2].data
        wave = [item for t in specExp for item in t] 

        noiseBin = dataSpec['ESPEC']

        pxSize = head['PIXSIZE']

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

        xAxis = (np.linspace(1, shapeX+1, shapeX+1)-head['CRPIX1']+1) *head['PIXSIZE']
        yAxis = (np.linspace(1, shapeY+1, shapeY+1)-head['CRPIX2']+1) *head['PIXSIZE']
        
        tab = fits.open(tableSpec)
        dataSpec = tab[1].data
        specExp = tab[2].data
        wave = [item for t in specExp for item in t] 

        pxSize = head['PIXSIZE']

        noiseBin = dataSpec['ESPEC']

        tabStar = fits.open(tableStar)
        dataStar = tabStar[1].data


        return wave,xAxis,yAxis,pxSize,noiseBin,dataTab,dataSpec,dataStar

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


    def updateLineArray(self,cfg_par,lineArr,result,lineInfo,binIDName,counter):
        
        fitRes = result.params.valuesdict()


        modName = cfg_par['gFit']['modName']
        lineArr['BIN_ID'][counter] = binIDName

        for ii in range(0,len(lineInfo['ID'])):

            lineName = str(lineInfo['Name'][ii])
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")
            
            lineName = lineName+str(int(lineInfo['Wave'][ii]))


            intR = fitRes['Wintln'+str(ii)]

            amp = fitRes['g1ln'+str(ii)+'_amplitude']
            ctr = fitRes['g1ln'+str(ii)+'_center']
            sig = fitRes['g1intln'+str(ii)]

            fwhm = fitRes['g1ln'+str(ii)+'_fwhm']
            height = fitRes['g1ln'+str(ii)+'_height']
            

            g1Ctr = cvP.lambdaVRad(np.exp(ctr),lineInfo['Wave'][ii])
            g1Sigma = cvP.lambdaVRad(np.exp(ctr+sig),lineInfo['Wave'][ii])-g1Ctr
            g1FWHM = cvP.lambdaVRad(np.exp(ctr+fwhm),lineInfo['Wave'][ii])-g1Ctr
            #amp_err = result.params[modName+'ln'+str(i)+'_amplitude'].stderr
            #sig_err = result.params[modName+'ln'+str(i)+'_sigma'].stderr
            #g1SigmaErr = self.lambdaVRad(np.exp(sig_err),lineInfo['Wave'][i])
            #cen_err = result.params[modName+'ln'+str(i)+'_center'].stderr  
            #g1CtrErr = self.lambdaVRad(np.exp(cen_err),lineInfo['Wave'][i])
            lineArr[lineName][counter] = int(lineInfo['Wave'][ii])     
            lineArr['g1_Amp_'+lineName][counter] = amp
            lineArr['g1_Height_'+lineName][counter] = height
            lineArr['g1_Centre_'+lineName][counter] = g1Ctr
            lineArr['g1_Sigma_'+lineName][counter] = g1Sigma
            lineArr['g1_FWHM_'+lineName][counter] = g1FWHM

            if modName != 'g1':

                amp = fitRes['g2ln'+str(ii)+'_amplitude']
                ctr = fitRes['g2ln'+str(ii)+'_center']
                sig = fitRes['g2intln'+str(ii)]
                fwhm = fitRes['g2ln'+str(ii)+'_fwhm']
                height = fitRes['g2ln'+str(ii)+'_height']

                g2Ctr = cvP.lambdaVRad(np.exp(ctr),lineInfo['Wave'][ii])
                g2Sigma = cvP.lambdaVRad(np.exp(ctr+sig),lineInfo['Wave'][ii])-g2Ctr
                g2FWHM = cvP.lambdaVRad(np.exp(ctr+fwhm),lineInfo['Wave'][ii])-g3Ctr

                #amp_err = result.params[modName+'ln'+str(i)+'_amplitude'].stderr
                #sig_err = result.params[modName+'ln'+str(i)+'_sigma'].stderr
                #g1SigmaErr = self.lambdaVRad(np.exp(sig_err),lineInfo['Wave'][i])
                #cen_err = result.params[modName+'ln'+str(i)+'_center'].stderr  
                #g1CtrErr = self.lambdaVRad(np.exp(cen_err),lineInfo['Wave'][i])
          
                lineArr['g2_Amp_'+lineName][counter] = amp
                lineArr['g2_Height_'+lineName][counter] = height
                lineArr['g2_Centre_'+lineName][counter] = g2Ctr
                lineArr['g2_Sigma_'+lineName][counter] = g2Sigma
                lineArr['g2_FWHM_'+lineName][counter] = g2FWHM

                if modName == 'g3':

                    amp = fitRes['g3ln'+str(ii)+'_amplitude']
                    ctr = fitRes['g3ln'+str(ii)+'_center']
                    sig = fitRes['g3intln'+str(ii)]
                    fwhm = fitRes['g3ln'+str(ii)+'_fwhm']
                    height = fitRes['g3ln'+str(ii)+'_height']

                    g3Ctr = cvP.lambdaVRad(np.exp(ctr),lineInfo['Wave'][ii])
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
                    lineArr['g3_Sigma_'+lineName][counter] = g3Sigma
                    lineArr['g3_FWHM_'+lineName][counter] = g3FWHM
            

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

    def binLineRatio(self,cfg_par,lineInfo):

        lineNameID=[]
        modName = cfg_par['gFit']['modName']

        for ii in range(0,len(lineInfo['ID'])):

            lineName = str(lineInfo['Name'][ii])
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")
            
            lineNameID.append(lineName+str(int(lineInfo['Wave'][ii])))
        
        hdul = fits.open(cfg_par['general']['outTableName'])
        lines = hdul['LineRes_'+cfg_par['gFit']['modName']].data 


        lineNameList=['BIN_ID']
        frmList=['i4']

        tot = lines['BIN_ID']

        #tot = np.column_stack((binColumn))

        if 'OIII5006' in lineNameID and 'Hb4861' in lineNameID:
            lrOHbG1 = np.divide(lines['g1_Amp_'+'OIII5006'],lines['g1_Amp_'+'Hb4861'])
            tot = np.column_stack((tot,lrOHbG1))
            lineNameList.append('G1-OIII5006/Hb4861')
            frmList.append('f8')

        if 'NII6583' in lineNameID and 'Ha6562' in lineNameID:
            lrNIIHaG1 = np.divide(lines['g1_Amp_'+'NII6583'],lines['g1_Amp_'+'Ha6562'])
            tot = np.column_stack((tot,lrNIIHaG1))
            lineNameList.append('G1-NII6583/Ha6562')
            frmList.append('f8')

        if 'OI6300' in lineNameID and 'Ha6562' in lineNameID:
            lrOIHaG1 = np.divide(lines['g1_Amp_'+'OI6300'],lines['g1_Amp_'+'Ha6562'])
            tot = np.column_stack((tot,lrOIHaG1))
            lineNameList.append('G1-OI6300/Ha6562')
            frmList.append('f8')

        if 'SII6716' in lineNameID and 'Ha6562' in lineNameID:
            lrSIIHaG1 = np.divide((lines['g1_Amp_'+'SII6716']+lines['g1_Amp_'+'SII6730']),lines['g1_Amp_'+'Ha6562'])
            tot = np.column_stack((tot,lrSIIHaG1))
            lineNameList.append('G1-SII6716/Ha6562')
            frmList.append('f8')

        if modName != 'g1':


            if 'OIII5006' in lineNameID and 'Hb4861' in lineNameID:
                lrOHbG2 = np.divide(lines['g2_Amp_'+'OIII5006'],lines['g2_Amp_'+'Hb4861'])
                lrOHb = np.divide((lines['g1_Amp_'+'OIII5006']+lines['g2_Amp_'+'OIII5006']),(lines['g1_Amp_'+'Hb4861']+lines['g2_Amp_'+'Hb4861']))
                tot = np.column_stack((tot,lrOHbG2,lrOHb))
                lineNameList.append('G2-OIII5006/Hb4861')
                lineNameList.append('ToT-OIII5006/Hb4861')
                frmList.append('f8')
                frmList.append('f8')

            
            if 'NII6583' in lineNameID and 'Ha6562' in lineNameID:
                lrNIIHaG2 = np.divide(lines['g2_Amp_'+'NII6583'],lines['g2_Amp_'+'Ha6562'])
                lrNIIHa = np.divide((lines['g1_Amp_'+'NII6583']+lines['g2_Amp_'+'NII6583']),(lines['g1_Amp_'+'Ha6562']+lines['g2_Amp_'+'Ha6562']))
                tot = np.column_stack((tot,lrNIIHaG2,lrNIIHa))
                lineNameList.append('G2-NII6583/Ha6562')
                lineNameList.append('ToT-NII6583/Ha6562')
                frmList.append('f8')
                frmList.append('f8')


            if 'OI6300' in lineNameID and 'Ha6562' in lineNameID:
                lrOIHaG2 = np.divide(lines['g2_Amp_'+'OI6300'],lines['g2_Amp_'+'Ha6562'])
                lrOIHa = np.divide((lines['g1_Amp_'+'OI6300']+lines['g2_Amp_'+'OI6300']),(lines['g1_Amp_'+'Ha6562']+lines['g2_Amp_'+'Ha6562']))
                tot = np.column_stack((tot,lrOIHaG2,lrOIHa))
                lineNameList.append('G2-OI6300/Ha6562')
                lineNameList.append('ToT-OI6300/Ha6562')
                frmList.append('f8')
                frmList.append('f8')


            if 'SII6716' in lineNameID and 'Ha6562' in lineNameID:
                lrSIIHaG2 = np.divide((lines['g2_Amp_'+'SII6716']+lines['g2_Amp_'+'SII6730']),lines['g2_Amp_'+'Ha6562'])
                lrSIIHa = np.divide((lines['g1_Amp_'+'SII6716']+lines['g1_Amp_'+'SII6730']+lines['g2_Amp_'+'SII6716']+lines['g2_Amp_'+'SII6730']),(lines['g1_Amp_'+'Ha6562']+lines['g2_Amp_'+'Ha6562']))
                tot = np.column_stack((tot,lrSIIHaG2,lrSIIHa))
                lineNameList.append('G2-SII6716/Ha6562')
                lineNameList.append('ToT-SII6716/Ha6562')
                frmList.append('f8')
                frmList.append('f8')

            if modName == 'g3':

                
                if 'OIII5006' in lineNameID and 'Hb4861' in lineNameID:
                    lrOHbG3 = np.divide(lines['g3_Amp_'+'OIII5006'],lines['g3_Amp_'+'Hb4861'])
                    lrOHb = np.divide((lines['g1_Amp_'+'OIII5006']+lines['g2_Amp_'+'OIII5006']+lines['g3_Amp_'+'OIII5006']),(lines['g1_Amp_'+'Hb4861']+
                        lines['g2_Amp_'+'Hb4861']+lines['g3_Amp_'+'Hb4861']))
                    tot = np.column_stack((tot,lrOHbG3))
                    lineNameList.append('G3-OIII5006/Hb4861')
                    lineNameList.append('ToT-OIII5006/Hb4861')
                    frmList.append('f8')
                    frmList.append('f8')

                
                if 'NII6583' in lineNameID and 'Ha6562' in lineNameID:
                    lrNIIHaG3 = np.divide(lines['g3_Amp_'+'NII6583'],lines['g3_Amp_'+'Ha6562'])
                    lrNIIHa = np.divide((lines['g1_Amp_'+'NII6583']+lines['g2_Amp_'+'NII6583']+lines['g3_Amp_'+'NII6583']),(lines['g1_Amp_'+'Ha6562']+
                        lines['g2_Amp_'+'Ha6562']+lines['g3_Amp_'+'Ha6562']))
                    tot = np.column_stack((tot,lrNIIHaG3,lrNIIHa))
                    lineNameList.append('G3-NII6583/Hb4861')
                    lineNameList.append('ToT-NII6583/Hb4861')
                    frmList.append('f8')
                    frmList.append('f8')

                if 'OI6300' in lineNameID and 'Ha6562' in lineNameID:
                    lrOIHaG3 = np.divide(lines['g3_Amp_'+'OI6300'],lines['g3_Amp_'+'Ha6562'])
                    lrOIHa = np.divide((lines['g1_Amp_'+'OI6300']+lines['g2_Amp_'+'OI6300']+lines['g3_Amp_'+'OI6300']),(lines['g1_Amp_'+'Ha6562']+
                        lines['g2_Amp_'+'Ha6562']+lines['g3_Amp_'+'Ha6562']))
                    tot = np.column_stack((tot,lrOIHaG3,lrOIHa))
                    lineNameList.append('G3-OI6300/Ha6562')
                    lineNameList.append('ToT-OI6300/Ha6562')
                    frmList.append('f8')
                    frmList.append('f8')


                if 'SII6716' in lineNameID and 'Ha6562' in lineNameID:
                    lrSIIHaG3 = np.divide((lines['g2_Amp_'+'SII6716']+lines['g2_Amp_'+'SII6730']),lines['g2_Amp_'+'Ha6562'])
                    lrSIIHa = np.divide((lines['g1_Amp_'+'SII6716']+lines['g1_Amp_'+'SII6730']+lines['g2_Amp_'+'SII6716']+lines['g2_Amp_'+'SII6730']+lines['g3_Amp_'+'SII6716']+lines['g3_Amp_'+'SII6730']),
                        (lines['g1_Amp_'+'Ha6562']+lines['g2_Amp_'+'Ha6562']+lines['g3_Amp_'+'Ha6562']))
                    tot = np.column_stack((tot,lrSIIHaG3,lrSIIHa))
                    lineNameList.append('G3-SII6716/Ha6562')
                    lineNameList.append('ToT-SII6716/Ha6562')
                    frmList.append('f8')
                    frmList.append('f8')


        t = Table(tot, names=(lineNameList))

        hdul.append(fits.BinTableHDU(t.as_array(), name='LineRatios_'+modName))


        indexSFK = np.where(np.logical_and(np.logical_and(np.log10(t['G1-OIII5006/Hb4861']) < 0.61 / (np.log10(t['G1-NII6583/Ha6562']) - 0.05) + 1.3,
            np.log10(t['G1-OIII5006/Hb4861'])<3),
            np.log10(t['G1-NII6583/Ha6562'])<-0.05)
            )
        indexSF = np.where(np.logical_and(np.logical_and(np.logical_and(np.log10(t['G1-OIII5006/Hb4861']) < 0.61 / (np.log10(t['G1-NII6583/Ha6562']) - 0.47) + 1.19, 
            np.log10(t['G1-OIII5006/Hb4861']) >= 0.61 / (np.log10(t['G1-NII6583/Ha6562']) - 0.05) + 1.3),
            np.log10(t['G1-OIII5006/Hb4861'])<3.),
            np.log10(t['G1-NII6583/Ha6562'])<0.47))
        indexAGN  = np.where(np.logical_and(np.logical_and(np.log10(t['G1-OIII5006/Hb4861']) >= 0.61 / (np.log10(t['G1-NII6583/Ha6562']) - 0.47) + 1.19,
            np.log10(t['G1-OIII5006/Hb4861'])<3.),
            np.log10(t['G1-NII6583/Ha6562'])<2.))
        
        LrOIII  = np.zeros(len(lines['BIN_ID']))*np.nan
        LrOIII[indexSFK] = 0.
        LrOIII[indexSF] = 1.
        LrOIII[indexAGN] = 2.


        indexSF = np.where(np.logical_and.reduce((np.log10(t['G1-OIII5006/Hb4861']) < 0.72 / (np.log10(t['G1-SII6716/Ha6562']) - 0.32) + 1.30,
                            np.log10(t['G1-OIII5006/Hb4861'])<3.,
                            np.log10(t['G1-SII6716/Ha6562'])<0.32,
                            np.log10(t['G1-OIII5006/Hb4861'])>=-2.,
                            np.log10(t['G1-SII6716/Ha6562'])>=-2.)))
        indexSey = np.where(np.logical_and.reduce((np.log10(t['G1-OIII5006/Hb4861']) >= 0.72 / (np.log10(t['G1-SII6716/Ha6562']) - 0.32) + 1.30, 
            np.log10(t['G1-OIII5006/Hb4861']) > 1.89* np.log10(t['G1-SII6716/Ha6562']) + 0.76,
                            np.log10(t['G1-OIII5006/Hb4861'])<3.,
                            np.log10(t['G1-SII6716/Ha6562'])<1.,
                            np.log10(t['G1-OIII5006/Hb4861'])>=-2.,
                            np.log10(t['G1-SII6716/Ha6562'])>=-2.)))
        indexLIN = np.where(np.logical_and.reduce((np.log10(t['G1-OIII5006/Hb4861']) >= 0.72 / (np.log10(t['G1-SII6716/Ha6562']) - 0.32) + 1.30,
            np.log10(t['G1-OIII5006/Hb4861']) < 1.89*np.log10(t['G1-SII6716/Ha6562']) + 0.76,
                            np.log10(t['G1-OIII5006/Hb4861'])<3.,
                            np.log10(t['G1-SII6716/Ha6562'])<1.,
                            np.log10(t['G1-OIII5006/Hb4861'])>=-2.,
                            np.log10(t['G1-SII6716/Ha6562'])>=-2.)))

        LrSII  = np.zeros(len(lines['BIN_ID']))*np.nan
        LrSII[indexSF] = 0.
        LrSII[indexSey] = 1.
        LrSII[indexLIN] = 2.

        indexSF = np.where(np.logical_and.reduce((np.log10(t['G1-OIII5006/Hb4861']) < 0.73 / (np.log10(t['G1-OI6300/Ha6562']) + 0.59) + 1.33,
            np.log10(t['G1-OIII5006/Hb4861'])<3.,
            np.log10(t['G1-OI6300/Ha6562'])<-0.59,
            np.log10(t['G1-OIII5006/Hb4861'])>=-2.,            
            np.log10(t['G1-OI6300/Ha6562'])>=-3.)))
        
        indexSey = np.where(np.logical_and.reduce((np.log10(t['G1-OIII5006/Hb4861']) >= 0.73 / (np.log10(t['G1-OI6300/Ha6562']) + 0.59) +1.33,
            np.log10(t['G1-OIII5006/Hb4861']) >= 1.18* np.log10(t['G1-OI6300/Ha6562']) + 1.30,
            np.log10(t['G1-OIII5006/Hb4861'])<3.,
            np.log10(t['G1-OI6300/Ha6562'])<2.,
            np.log10(t['G1-OIII5006/Hb4861'])>=-2.,            
            np.log10(t['G1-OI6300/Ha6562'])>=-3.)))
        
        indexLIN = np.where(np.logical_and.reduce((np.log10(t['G1-OIII5006/Hb4861']) >= 0.73 / (np.log10(t['G1-OI6300/Ha6562']) + 0.59)+1.33, 
            np.log10(t['G1-OIII5006/Hb4861']) < 1.18* np.log10(t['G1-OI6300/Ha6562']) + 1.30,
            np.log10(t['G1-OIII5006/Hb4861'])<3.,
            np.log10(t['G1-OI6300/Ha6562'])<2.,
            np.log10(t['G1-OIII5006/Hb4861'])>=-2.,            
            np.log10(t['G1-OI6300/Ha6562'])>=-3.)))

        LrOI  = np.zeros(len(lines['BIN_ID']))*np.nan
        LrOI[indexSF] = 0.
        LrOI[indexSey] = 1.
        LrOI[indexLIN] = 2.

        tt=Table([lines['BIN_ID'],LrOIII,LrSII,LrOI],names=('BIN_ID','G1-BPT_OIII','G1-BPT_SII','G1-BPT_OI'))

        if modName != 'g1':

            indexSFK = np.where(np.log10(t['G2-OIII5006/Hb4861']) < 0.61 / (np.log10(t['G2-NII6583/Ha6562']) - 0.05) + 1.3)
            indexSF = np.where(np.logical_and(np.log10(t['G2-OIII5006/Hb4861']) < 0.61 / (np.log10(t['G2-NII6583/Ha6562']) - 0.47) + 1.19, 
                np.log10(t['G2-OIII5006/Hb4861']) >= 0.61 / (np.log10(t['G2-NII6583/Ha6562']) - 0.05) + 1.3))
            indexAGN  = np.where((np.log10(t['G2-OIII5006/Hb4861']) >= 0.61 / (np.log10(t['G2-NII6583/Ha6562']) - 0.47) + 1.19))    
            
            LrOIII  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrOIII[indexSFK] = 0.
            LrOIII[indexSF] = 1.
            LrOIII[indexAGN] = 2.


            indexSF = np.where(np.log10(t['G2-OIII5006/Hb4861']) < 0.72 / (np.log10(t['G2-SII6716/Ha6562']) - 0.32) + 1.30)
            indexSey = np.where(np.logical_and(np.log10(t['G2-OIII5006/Hb4861']) >= 0.72 / (np.log10(t['G2-SII6716/Ha6562']) - 0.32) + 1.30, 
                np.log10(t['G2-OIII5006/Hb4861']) > 1.89* np.log10(t['G1-SII6716/Ha6562']) + 0.76))
            indexLIN = np.where(np.logical_and(np.log10(t['G2-OIII5006/Hb4861']) >= 0.72 / (np.log10(t['G2-SII6716/Ha6562']) - 0.32) + 1.30,
                np.log10(t['G2-OIII5006/Hb4861']) < 1.89*np.log10(t['G2-SII6716/Ha6562']) + 0.76))

            LrSII  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrSII[indexSF] = 0.
            LrSII[indexSey] = 1.
            LrSII[indexLIN] = 2.

            indexSF = np.where(np.log10(t['G2-OIII5006/Hb4861']) < 0.73 / (np.log10(t['G2-OI6300/Ha6562']) + 0.59) + 1.33)
            indexSey = np.where(np.logical_and( np.log10(t['G2-OIII5006/Hb4861']) >= 0.73 / (np.log10(t['G2-OI6300/Ha6562']) + 0.59) +1.33,
                np.log10(t['G2-OIII5006/Hb4861']) >= 1.18* np.log10(t['G2-OI6300/Ha6562']) + 1.30))
            indexLIN = np.where(np.logical_and( np.log10(t['G2-OIII5006/Hb4861']) >= 0.73 / (np.log10(t['G2-OI6300/Ha6562']) + 0.59)+1.33, 
                np.log10(t['G2-OIII5006/Hb4861']) < 1.18* np.log10(t['G2-OI6300/Ha6562']) + 1.30))

            LrOI  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrOI[indexSF] = 0.
            LrOI[indexSey] = 1.
            LrOI[indexLIN] = 2.

            tt.add_column(Column(LrOIII,name='G2-BPT_OIII'))
            tt.add_column(Column(LrSII,name='G2-BPT_SII'))
            tt.add_column(Column(LrOI,name='G2-BPT_OI'))

            indexSFK = np.where(np.log10(t['ToT-OIII5006/Hb4861']) < 0.61 / (np.log10(t['ToT-NII6583/Ha6562']) - 0.05) + 1.3)
            indexSF = np.where(np.logical_and(np.log10(t['ToT-OIII5006/Hb4861']) < 0.61 / (np.log10(t['ToT-NII6583/Ha6562']) - 0.47) + 1.19, 
                np.log10(t['ToT-OIII5006/Hb4861']) >= 0.61 / (np.log10(t['ToT-NII6583/Ha6562']) - 0.05) + 1.3))
            indexAGN  = np.where((np.log10(t['ToT-OIII5006/Hb4861']) >= 0.61 / (np.log10(t['ToT-NII6583/Ha6562']) - 0.47) + 1.19))    
            
            LrOIII  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrOIII[indexSFK] = 0.
            LrOIII[indexSF] = 1.
            LrOIII[indexAGN] = 2.


            indexSF = np.where(np.log10(t['ToT-OIII5006/Hb4861']) < 0.72 / (np.log10(t['ToT-SII6716/Ha6562']) - 0.32) + 1.30)
            indexSey = np.where(np.logical_and(np.log10(t['ToT-OIII5006/Hb4861']) >= 0.72 / (np.log10(t['ToT-SII6716/Ha6562']) - 0.32) + 1.30, 
                np.log10(t['ToT-OIII5006/Hb4861']) > 1.89* np.log10(t['G1-SII6716/Ha6562']) + 0.76))
            indexLIN = np.where(np.logical_and(np.log10(t['ToT-OIII5006/Hb4861']) >= 0.72 / (np.log10(t['ToT-SII6716/Ha6562']) - 0.32) + 1.30,
                np.log10(t['ToT-OIII5006/Hb4861']) < 1.89*np.log10(t['ToT-SII6716/Ha6562']) + 0.76))

            LrSII  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrSII[indexSF] = 0.
            LrSII[indexSey] = 1.
            LrSII[indexLIN] = 2.

            indexSF = np.where(np.log10(t['ToT-OIII5006/Hb4861']) < 0.73 / (np.log10(t['ToT-OI6300/Ha6562']) + 0.59) + 1.33)
            indexSey = np.where(np.logical_and( np.log10(t['ToT-OIII5006/Hb4861']) >= 0.73 / (np.log10(t['ToT-OI6300/Ha6562']) + 0.59) +1.33,
                np.log10(t['ToT-OIII5006/Hb4861']) >= 1.18* np.log10(t['ToT-OI6300/Ha6562']) + 1.30))
            indexLIN = np.where(np.logical_and( np.log10(t['ToT-OIII5006/Hb4861']) >= 0.73 / (np.log10(t['ToT-OI6300/Ha6562']) + 0.59)+1.33, 
                np.log10(t['ToT-OIII5006/Hb4861']) < 1.18* np.log10(t['ToT-OI6300/Ha6562']) + 1.30))

            LrOI  = np.zeros(len(lines['BIN_ID']))*np.nan
            LrOI[indexSF] = 0.
            LrOI[indexSey] = 1.
            LrOI[indexLIN] = 2.

            tt.add_column(Column(LrOIII,name='ToT-BPT_OIII'))
            tt.add_column(Column(LrSII,name='ToT-BPT_SII'))
            tt.add_column(Column(LrOI,name='ToT-BPT_OI'))

            if modName == 'g3':

                indexSFK = np.where(np.log10(t['G3-OIII5006/Hb4861']) < 0.61 / (np.log10(t['G3-NII6583/Ha6562']) - 0.05) + 1.3)
                indexSF = np.where(np.logical_and(np.log10(t['G3-OIII5006/Hb4861']) < 0.61 / (np.log10(t['G3-NII6583/Ha6562']) - 0.47) + 1.19, 
                    np.log10(t['G3-OIII5006/Hb4861']) >= 0.61 / (np.log10(t['G3-NII6583/Ha6562']) - 0.05) + 1.3))
                indexAGN  = np.where((np.log10(t['G3-OIII5006/Hb4861']) >= 0.61 / (np.log10(t['G3-NII6583/Ha6562']) - 0.47) + 1.19))    
                
                LrOIII  = np.zeros(len(lines['BIN_ID']))*np.nan
                LrOIII[indexSFK] = 0.
                LrOIII[indexSF] = 1.
                LrOIII[indexAGN] = 2.


                indexSF = np.where(np.log10(t['G3-OIII5006/Hb4861']) < 0.72 / (np.log10(t['G3-SII6716/Ha6562']) - 0.32) + 1.30)
                indexSey = np.where(np.logical_and(np.log10(t['G3-OIII5006/Hb4861']) >= 0.72 / (np.log10(t['G3-SII6716/Ha6562']) - 0.32) + 1.30, 
                    np.log10(t['G3-OIII5006/Hb4861']) > 1.89* np.log10(t['G1-SII6716/Ha6562']) + 0.76))
                indexLIN = np.where(np.logical_and(np.log10(t['G3-OIII5006/Hb4861']) >= 0.72 / (np.log10(t['G3-SII6716/Ha6562']) - 0.32) + 1.30,
                    np.log10(t['G3-OIII5006/Hb4861']) < 1.89*np.log10(t['G3-SII6716/Ha6562']) + 0.76))

                LrSII  = np.zeros(len(lines['BIN_ID']))*np.nan
                LrSII[indexSF] = 0.
                LrSII[indexSey] = 1.
                LrSII[indexLIN] = 2.

                indexSF = np.where(np.log10(t['G3-OIII5006/Hb4861']) < 0.73 / (np.log10(t['G3-OI6300/Ha6562']) + 0.59) + 1.33)
                indexSey = np.where(np.logical_and( np.log10(t['G3-OIII5006/Hb4861']) >= 0.73 / (np.log10(t['G3-OI6300/Ha6562']) + 0.59) +1.33,
                    np.log10(t['G3-OIII5006/Hb4861']) >= 1.18* np.log10(t['G3-OI6300/Ha6562']) + 1.30))
                indexLIN = np.where(np.logical_and( np.log10(t['G3-OIII5006/Hb4861']) >= 0.73 / (np.log10(t['G3-OI6300/Ha6562']) + 0.59)+1.33, 
                    np.log10(t['G3-OIII5006/Hb4861']) < 1.18* np.log10(t['G3-OI6300/Ha6562']) + 1.30))

                LrOI  = np.zeros(len(lines['BIN_ID']))*np.nan
                LrOI[indexSF] = 0.
                LrOI[indexSey] = 1.
                LrOI[indexLIN] = 2.

                tt.add_column(Column(LrOIII,name='G3-BPT_OIII'))
                tt.add_column(Column(LrSII,name='G3-BPT_SII'))
                tt.add_column(Column(LrOI,name='G3-BPT_OI'))


        hdul.append(fits.BinTableHDU(tt.as_array(), name='BPT_'+modName))
        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)

        
        return

