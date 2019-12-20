#!/usr/bin/env python3.8
 
from astropy.io import fits
from astropy import wcs
import os
import numpy as np
from astropy.table import Table, Column

from scavengers import cvPlay, tPlay
cvP = cvPlay.convert()
tP = tPlay.tplay()

class starsub(object):

    def makeCubes(self,cfg_par):

        workDir = cfg_par['general']['workdir']

        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec,dataStar = tP.openPPXFforSubtraction(cfg_par,workDir+cfg_par['general']['tableBinName'],
            workDir+cfg_par['general']['tableSpecName'],workDir+cfg_par['general']['tableStarName'])



        data=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        Stars=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        Lines=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])

        header = self.makeHeader(cfg_par, wave, pxSize)

        hdu = fits.PrimaryHDU(data=data,header=header)

        diffusion = 1e-5
        del header['LATPOLE']
        del header['LONPOLE']
        xxVec = []
        yyVec = []
        
        #create AllSpectra datacube
        for i in range(0,vorBinInfo['ID'].shape[0]):
            #print xAxis
            #print yAxis
            indexX = ((xAxis < (np.round(vorBinInfo['X'][i],1)+diffusion)) & 
                      ((np.round(vorBinInfo['X'][i],1)-diffusion) < xAxis))
            indexY = ((yAxis < (np.round(vorBinInfo['Y'][i],1)+diffusion)) & 
                      ((np.round(vorBinInfo['Y'][i],1)-diffusion) < yAxis))       

            xx = np.where(indexX)[0]
            yy = np.where(indexY)[0]

            #print yy

            indexBin =  vorBinInfo['BIN_ID'][i]
            
            if indexBin>0 and xx and yy: 
                   
                xxVec.append(xx[0])
                yyVec.append(yy[0])

                tmp = np.array(dataSpec[indexBin][0][:])
                tmp = tmp.tolist()
                data[:,yy[0],xx[0]] = tmp

                tmp = np.array(dataStar[indexBin][1][:])
                tmp = tmp.tolist()
                Stars[:,yy[0],xx[0]] = tmp

            else:
                pass

        xxVecArr= Column(np.array(xxVec), name='PixX')
        yyVecArr= Column(np.array(yyVec), name='PixY')

        #yyVecArr=np.array(yyVec,dtype={'names':('PixY')})
        #print(vorBinInfo.shape)
        #print(xxVecArr.shape)

        t = Table(vorBinInfo)
        t.add_column(xxVecArr,index=-1)
        t.add_column(yyVecArr,index=-1) 
        #vorBinInfo = np.column_stack((vorBinInfo,xxVec,yyVec))
        #vorBinInfo = np.vstack([vorBinInfo,yyVecArr])
        print(vorBinInfo.dtype.names)
        tab = fits.open(workDir+cfg_par['general']['tableBinName'])
        head = tab[0].header
        #data = tab[0].data
        #tab[1] = vorBinInfo    

        empty_primary = fits.PrimaryHDU(header=head)           

        #t2 = fits.BinTableHDU.from_columns(t,name='vorBinInfo')
        hdul = fits.HDUList([empty_primary])      

        hdul.append(fits.BinTableHDU(t.as_array(), name='vorBinInfo'))


        #hdul.append(t2)  
        hdul.writeto(workDir+cfg_par['general']['outVorTableName'],overwrite=True)

        Lines = np.subtract(data,Stars)

        outputs = cfg_par['starSub']['outputs']
        
        if outputs == 'lines':
            fits.writeto(cfg_par['general']['outLines'],Lines,header,overwrite=True)
            print('''\t+---------+\n\t Line Cube saved\n\t+---------+''')     
            return 
        else: 
            if outputs == 'stars':
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                print('''\t+---------+\n\t Star Cube saved\n\t+---------+''')    
                return 
            elif outputs == 'data':
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('''\t+---------+\n\t Data Cube saved\n\t+---------+''')                 
            elif outputs=='all':
                fits.writeto(cfg_par['general']['outLines'],Lines,header,overwrite=True)
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('''\t+---------+\n\t All Cubes saved\n\t+---------+''')    
                return


    def makeHeader(self,cfg_par,wave,pxSize):
        
        # these are arbitrary but must correspond, I choose the centre of the galaxy(usually offset with fov)
        crPix3 =cfg_par['starSub']['pixZ']
        crPix1=cfg_par['starSub']['pixX']
        crPix2=cfg_par['starSub']['pixY']

        crVal3 = np.min(wave)
        deltaLambda = (np.max(wave)-np.min(wave))/len(wave)

        #only for Fornax E
        #crVal3 = np.exp(dataSpecStraight[0])
        #crVal3 =4750.2734375
        crpix3 = 1
        crVal1=cfg_par['starSub']['ra']
        crVal2=cfg_par['starSub']['dec']

        ra = cvP.hms2deg(crVal1)
        dec = cvP.dms2deg(crVal2)

        w = wcs.WCS(naxis=3)


        # Set up an "Airy's zenithal" projection
        # Vector properties may be set with Python lists, or Numpy arrays
        w.wcs.crpix = [crPix1, crPix2, crPix3]
        w.wcs.cdelt = np.array([-pxSize,pxSize,deltaLambda])
        w.wcs.crval = [ra, dec, crVal3]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN","AWAV"]

        header = w.to_header()

        #w.wcs.set_pv([(2, 1, 45.0)])
        header['EQUINOX'] = 2000
        header['CRDER3'] = 0.026
        header['CTYPE3'] = "AWAV"
        header['CUNIT3'] = "Angstrom"
        #header['CRPIX1'] = crPix1
        #head['CRPIX2'] = crPix2

        return header
