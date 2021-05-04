#!/usr/bin/env python

'''

Set of generic tools to play with fits files.

'''


import sys, os, string
import numpy as np

from scavengers import headPlay, cvPlay
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from reproject import reproject_interp, reproject_exact
from astropy.io import fits
from astropy.coordinates import Angle


hP=headPlay.headplay()
cP=cvPlay.convert()
class fitsplay():
    '''This class plays with .fits files in useful ways.
    
    '''


    def divFits(self,cfg_par,inFile,divFile,output=False):
        '''Divides one fit file for another

            Parameters
            ----------
            inFile: str
                full path of input .fits file

            divFile: str
                full path of file to divide for input

            outMap: str, optional
                full path to output map.

            Returns
            -------
                
                outMap: str
                    full path to output surface brightness map.
            

        '''

        
        
        hh,dd = hP.cleanHead(inFile,writeFile=False)

        subh,sub = hP.cleanHead(divFile,writeFile=False)

        centre = SkyCoord(ra=hh['CRVAL1']*u.degree, dec=hh['CRVAL2']*u.degree, frame='fk5')
        size = u.Quantity((hh['CRVAL1'],cfg_par['moments']['sizePlots']), u.arcmin)

        dd = np.divide(dd,sub)
        dd[np.isnan(dd)] = 0.0

        dd[np.isinf(dd)] = 0.0

        if output==False:
            aaa = str.split(fileName, '.fits')
            output=aaa[0]+'_div.fits'
        
        fits.writeto(output,dd,hh,overwrite=True)

        return 0


    def coordCut(self,filename,rap1,decp1,rap2,decp2):

   
        rap1 = cP.hms2deg(rap1)
        rap2 = cP.hms2deg(rap2)
        decp1 = cP.dms2deg(decp1)
        decp2 = cP.dms2deg(decp2)
        hh,dd = hP.cleanHead(filename,writeFile=False)

        w = WCS(hh)    

        xmin,ymin=w.wcs_world2pix(rap1,decp1,0)
        xmin=int(np.round(xmin,0))
        ymin=int(np.round(ymin,0))
        
        xmax,ymax=w.wcs_world2pix(rap2,decp2,0)
        xmax=int(np.round(xmax,0))
        ymax=int(np.round(ymax,0))
        naxis1=xmax-xmin
        naxis2=ymax-ymin
        
        raCtr,decCtr=w.wcs_pix2world(xmin+(xmax-xmin)/2,ymin+(ymax-ymin)/2,0)


        hh['NAXIS1']=naxis1
        hh['NAXIS2']=naxis2
        hh['CRPIX1']=naxis1/2
        hh['CRPIX2']=naxis2/2
        hh['CRVAL1']=float(raCtr)
        hh['CRVAL2']=float(decCtr)  
             
        aaa = str.split(filename, '.fits')

        output=aaa[0]+'_coordCut.fits'
        print(output)
        fits.writeto(output,dd[ymin:ymax,xmin:xmax],hh,overwrite=True)

        return output

    def coordCutCube(self,filename,rap1,decp1,rap2,decp2):

   
        rap1 = cP.hms2deg(rap1)
        rap2 = cP.hms2deg(rap2)
        decp1 = cP.dms2deg(decp1)
        decp2 = cP.dms2deg(decp2)
        hh,dd = hP.cleanHead(filename,writeFile=False)

        w = WCS(hh)    

        xmin,ymin=w.wcs_world2pix(rap1,decp1,0)
        xmin=int(np.round(xmin,0))
        ymin=int(np.round(ymin,0))
        
        xmax,ymax=w.wcs_world2pix(rap2,decp2,0)
        xmax=int(np.round(xmax,0))
        ymax=int(np.round(ymax,0))
        naxis1=xmax-xmin
        naxis2=ymax-ymin
        
        raCtr,decCtr=w.wcs_pix2world(xmin+(xmax-xmin)/2,ymin+(ymax-ymin)/2,0)


        hh['NAXIS1']=naxis1
        hh['NAXIS2']=naxis2
        hh['CRPIX1']=naxis1/2
        hh['CRPIX2']=naxis2/2
        hh['CRVAL1']=float(raCtr)
        hh['CRVAL2']=float(decCtr)  
             
        aaa = str.split(filename, '.fits')

        output=aaa[0]+'_coordCut.fits'
        print(output)
        fits.writeto(output,dd[ymin:ymax,xmin:xmax],hh,overwrite=True)

        return output


    def coordCentreCut(self,filename,centreCoords,size,ychans=None):

   
        raC = cP.hms2deg(centreCoords[0])
        decC = cP.dms2deg(centreCoords[1])
        hh,dd = hP.cleanHead(filename,writeFile=False)

        w = WCS(hh)    

        xC,yC=w.wcs_world2pix(raC,decC,0)
        xC=int(np.round(xC,0))
        yC=int(np.round(yC,0))
        
        pixSize = Angle(hh['CDELT2'],u.deg)
        HalfSize = int(round(size/pixSize.arcminute)/2)
        #HalfSize = int(round(size/pixSize.arcminute)/2)
        
        xmax = xC+HalfSize
        xmin = xC-HalfSize
        ymax = yC+HalfSize
        ymin = yC-HalfSize
        
        naxis1=xmax-xmin
        naxis2=ymax-ymin
        
        raCtr,decCtr=w.wcs_pix2world(xmin+(xmax-xmin)/2,ymin+(ymax-ymin)/2,0)

        hh=fits.getheader(filename)
        hh['NAXIS1']=naxis1
        hh['NAXIS2']=naxis2
        hh['CRPIX1']=naxis1/2
        hh['CRPIX2']=naxis2/2
        hh['CRVAL1']=float(raC)
        hh['CRVAL2']=float(decC)  
             
        aaa = str.split(filename, '.fits')

        output=aaa[0]+'_coordCut.fits'
        print(output)
        if fits.getdata(filename).shape !=2:
            fits.writeto(output,fits.getdata(filename)[:,ymin:ymax,xmin:xmax],hh,overwrite=True)
        if ychans is not None:
            crval3=hh['CRVAL3']+(hh['CDELT3']*ychans[0])
            hh['CRVAL3']=crval3
            print(crval3)
            fits.writeto(output,fits.getdata(filename)[ychans[0]:ychans[1],ymin:ymax,xmin:xmax],hh,overwrite=True)
        if fits.getdata(filename).shape ==2:
            fits.writeto(output,fits.getdata(filename)[ymin:ymax,xmin:xmax],hh,overwrite=True)
        return output





