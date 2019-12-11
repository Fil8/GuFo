from astropy.io import fits
from astropy import wcs
import os
import numpy as np

from scavengers import cvPlay, tPlay
cvP = cvPlay.convert()
tP = tPlay.tplay()

class starsub:

    def makeCubes(self,cfg_par):

        workDir = cfg_par['general']['workdir']

        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo = tP.openTablesPPXF(cfg_par,workDir+cfg_par['general']['tableBinName'],
            workDir+cfg_par['general']['tableSpecName'])

        data=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        Stars=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        Lines=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])

        self.makeHeader(cfg_par, wave)

        hdu = fits.PrimaryHDU(data=data,header=header)

        del header['LATPOLE']
        del header['LONPOLE']

        #create AllSpectra datacube
        for i in xrange(0,dataTab['ID'].shape[0]):
            #print xAxis
            #print yAxis
            indexX = ((xAxis < (np.round(dataTab['X'][i],1)+diffusion)) & 
                      ((np.round(dataTab['X'][i],1)-diffusion) < xAxis))
            indexY = ((yAxis < (np.round(dataTab['Y'][i],1)+diffusion)) & 
                      ((np.round(dataTab['Y'][i],1)-diffusion) < yAxis))       

            xx = np.where(indexX)[0]
            yy = np.where(indexY)[0]

            indexBin =  dataTab['BIN_ID'][i]
            
            if indexBin>0 and xx and yy: 
                       
                tmp = np.array(dataSpec[indexBin][0][:])
                tmp = tmp.tolist()
                data[:,yy[0],xx[0]] = tmp

                tmp = np.array(dataStar[indexBin][1][:])
                tmp = tmp.tolist()
                Stars[:,yy[0],xx[0]] = tmp
            
            else:
                pass

        Lines = np.subtract(data,Stars)

        outputs = cfg_par['starSub']['outputs']
        
        if outputs == 'lines':
            fits.writeto(self.cfg_par['general']['outLines'],Lines,header,overwrite=True)
            print('''\t+---------+\n\t Line Cube saved\n\t+---------+''')     
            return 
        else: 
            if outputs == 'stars':
                fits.writeto(self.cfg_par['general']['outStars'],Stars,header,overwrite=True)
                print('''\t+---------+\n\t Star Cube saved\n\t+---------+''')    
                return 
            elif outputs == 'data':
                fits.writeto(self.cfg_par['general']['outCube'],data,header,overwrite=True)
                print('''\t+---------+\n\t Data Cube saved\n\t+---------+''')                 
            elif outputs=='all':
                fits.writeto(self.cfg_par['general']['outLines'],Lines,header,overwrite=True)
                fits.writeto(self.cfg_par['general']['outStars'],Stars,header,overwrite=True)
                fits.writeto(self.cfg_par['general']['outCube'],data,header,overwrite=True)
                print('''\t+---------+\n\t All Cubes saved\n\t+---------+''')    
                return


    def makeHeader(self,cfg_par,wave):
        
        #workDir = cfg_par['general']['workdir']
        #tab = fits.open(workDir+cfg_par['general']['tableBinName'])
#
        #head = tab[0].header
        #fits.close()

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

        pixelSize = head['PIXSIZE']/3600.

        # Set up an "Airy's zenithal" projection
        # Vector properties may be set with Python lists, or Numpy arrays
        w.wcs.crpix = [crPix1, crPix2, crPix3]
        w.wcs.cdelt = np.array([-pixelSize,pixelSize,deltaLambda])
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
