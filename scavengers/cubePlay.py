#!/usr/bin/env python3.6

import os, sys
import yaml

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult
from lmfit.model import load_modelresult

from scipy.ndimage import map_coordinates


from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy import wcs

import numpy as np
#import numpy.ma as ma

import tPlay,cvPlay,bptPlot,momPlot

tP = tPlay.tplay()
cvP = cvPlay.convert()
mPl = momPlot.MOMplot()

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
        print(vel)
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


    def makeLineCubes(self,cfg_par):

        cubeletsDir = cfg_par['general']['cubeletsDir']
        cubeDir = cfg_par['general']['cubeDir']
        modName = cfg_par['gFit']['modName']
        momDir = cfg_par['general']['momDir']

        if cfg_par['cubelets']['cube'] == 'vorLine': 
            f = fits.open(cfg_par['general']['dataCubeName'])
            dd = f[0].data
            hh = f[0].header

        elif cfg_par['cubelets']['cube'] == 'fitLine': 
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

        elif cfg_par['cubelets']['cube'] == 'residuals':
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
 
        velRange = float(cfg_par['cubelets']['velRange'])
       
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
            
            velRange = cvP.vRadLambda(cfg_par['cubelets']['velRange'],
                lineInfo['Wave'][ii] )-lineInfo['Wave'][ii] 

            #velRangeLeft = cvP.vRadLambda(cfg_par['cubelets']['velRange'],
            #    lineInfo['Wave'][ii] )+lineInfo['Wave'][ii] 

            #waveMin =  np.log(lineInfo['Wave'][ii] - lineInfo['lineRangeAng'][ii])
            #waveMax =  np.log(lineInfo['Wave'][ii] + lineInfo['lineRangeAng'][ii])
            
            waveMin =  np.log(lineInfo['Wave'][ii] - velRange)
            waveMax =  np.log(lineInfo['Wave'][ii] + velRange)
            

            idxMin = int(np.where(abs(wave-waveMin)==abs(wave-waveMin).min())[0]) 
            idxMax = int(np.where(abs(wave-waveMax)==abs(wave-waveMax).min())[0] )
            
            dCbl = dd[idxMin:idxMax,:,:]

            waveAng=np.exp(wave[idxMin:idxMax])

            header = self.makeHeader(cfg_par,lineInfo['Wave'][ii],hh,waveAng)
            
            if cfg_par['cubelets']['cube'] == 'vorLine': 
                outCubelet = cubeletsDir+str(lineNameStr)+'_measVor.fits'
            elif cfg_par['cubelets']['cube'] == 'fitLine': 
                outCubelet = cubeletsDir+str(lineNameStr)+'_fit_'+modName+'.fits'   
            elif cfg_par['cubelets']['cube'] == 'residuals':        
                outCubelet = cubeletsDir+str(lineNameStr)+'_res_'+modName+'.fits'        
            else:
                outCubelet = cubeletsDir+str(lineNameStr)+'_meas.fits'

            fits.writeto(outCubelet,np.flip(dCbl,axis=0),header,overwrite=True)

        return

    def pvDiagram(self,cfg_par,inCubelet):
        '''
        from SoFiA
        '''


        f = fits.open(inCubelet)
        subcube = f[0].data
        headerCubelets = f[0].header

        # Centres and bounding boxes
        Xc = cfg_par['cubelets']['raCentre']
        Yc = cfg_par['cubelets']['decCentre']
        Xc = cvP.hms2deg(Xc)
        Yc = cvP.dms2deg(Yc)

        header2d = headerCubelets.copy()
        self.delete_3rd_axis(header2d)

        w = wcs.WCS(header2d)
        Xc,Yc=w.wcs_world2pix(Xc,Yc,0)
        if cfg_par['cubelets']['raMin'] != False:
            Xmin = cfg_par['cubelets']['raMin']
            Ymin = cfg_par['cubelets']['decMin']
            Xmax = cfg_par['cubelets']['raMax']
            Ymax = cfg_par['cubelets']['decMax']
        else:
            Xmin = 0
            Ymin = 0
            Xmax = subcube.shape[2]
            Ymax = subcube.shape[1]

        kin_pa = np.radians(float(cfg_par['cubelets']['pa']))

        pv_sampling = 10
        pv_r = np.arange(-max(subcube.shape[1:]), max(subcube.shape[1:]) - 1 + 1.0 / pv_sampling, 1.0 / pv_sampling)
        #pv_y = Yc - float(YminNew) + pv_r * np.cos(kin_pa)
        #pv_x = Xc - float(XminNew) - pv_r * np.sin(kin_pa)

        pv_y = Yc  + pv_r * np.cos(kin_pa)
        pv_x = Xc  - pv_r * np.sin(kin_pa)

        pv_x, pv_y = pv_x[(pv_x >= 0) * (pv_x <= subcube.shape[2] - 1)], pv_y[(pv_x >= 0) * (pv_x <= subcube.shape[2] - 1)]
        pv_x, pv_y = pv_x[(pv_y >= 0) * (pv_y <= subcube.shape[1] - 1)], pv_y[(pv_y >= 0) * (pv_y <= subcube.shape[1] - 1)]
        pv_x.resize((1, pv_x.shape[0]))
        pv_y.resize((pv_x.shape))
        pv_coords = np.concatenate((pv_y, pv_x), axis=0)
        pv_array=[]
        for jj in range(subcube.shape[0]):
            plane = map_coordinates(subcube[jj], pv_coords)
            plane = [plane[ii::pv_sampling] for ii in range(pv_sampling)]
            plane = np.array([ii[:plane[-1].shape[0]] for ii in plane])
            pv_array.append(plane.mean(axis=0))
        pv_array = np.array(pv_array)
        
        hdu = fits.PrimaryHDU(data=pv_array, header=headerCubelets)
        hdulist = fits.HDUList([hdu])
        hdulist[0].header["CTYPE1"] = "PV--DIST"
        hdulist[0].header["CDELT1"] = hdulist[0].header["CDELT2"]
        hdulist[0].header["CRVAL1"] = 0
        hdulist[0].header["CRPIX1"] = pv_array.shape[1] / 2
        hdulist[0].header["CTYPE2"] = hdulist[0].header["CTYPE3"]
        hdulist[0].header["CDELT2"] = hdulist[0].header["CDELT3"]
        hdulist[0].header["CRVAL2"] = hdulist[0].header["CRVAL3"]
        hdulist[0].header["CRPIX2"] = hdulist[0].header["CRPIX3"]
        hdulist[0].header["ORIGIN"] = 'GaNGiaLF'
        
        self.delete_3rd_axis(hdulist[0].header)
        
        outCubelet=str.split(os.path.basename(inCubelet),'.')[0]
        outCubelet=outCubelet+'_'+str(cfg_par['cubelets']['pa'])+'-pv.fits'
        name = cfg_par['general']['pvDir'] + outCubelet
        print(name)
        hdulist.writeto(name,overwrite=True)

        return

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


            if cfg_par['cubelets']['cube'] == 'vorLine': 
                inCubelet = cubeletsDir+str(lineNameStr)+'_measVor.fits'
            elif cfg_par['cubelets']['cube'] == 'fitLine': 
                inCubelet = cubeletsDir+str(lineNameStr)+'_fit_'+modName+'.fits'   
            elif cfg_par['cubelets']['cube'] == 'residuals':        
                inCubelet = cubeletsDir+str(lineNameStr)+'_res_'+modName+'.fits' 
            else:
                inCubelet = cubeletsDir+str(lineNameStr)+'_meas.fits'

            self.pvDiagram(cfg_par,inCubelet,pa)



            return
