#!/usr/bin/env python

'''

Set of generic tools to play with fits files.

'''


import sys, os, string
import numpy as np

from scavengers import headPlay, cvPlay
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.nddata import Cutout2D
from reproject import reproject_interp, reproject_exact
from astropy.io import fits
from astropy.coordinates import Angle
from astropy.wcs import WCS

hP=headPlay.headplay()
cP=cvPlay.convert()

class fitsplay():
    '''This class plays with .fits files in useful ways.
    
    '''


    def divFits2D(self,cfg_par,inFile,divFile,output=False):
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

    def divFits3D(self,inFile,divFile,cfg_par=None,output=False):
        '''Divides one fitsFile file for another. Good for pbcorr datacubes.

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

        dd = fits.getdata(inFile)
        hh = fits.getheader(inFile)
        sub = fits.getdata(divFile)
        print(dd.shape[0],sub.shape[0])
        #for i in (0,np.min(dd.shape[0],sub.shape[0])):
        divide= np.divide(dd,sub)
        #dd[np.isnan(dd)] = 0.0

        #dd[np.isinf(dd)] = 0.0

        if output==False:
            aaa = str.split(inFile, '.fits')
            output=aaa[0]+'_pbcorr.fits'
        
        if dd.shape[0]<sub.shape[0]:
            hh['CDELT3'] = dd.shape[0]

        fits.writeto(output,divide,hh,overwrite=True)

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
        fits.writeto(output,dd[ymin:ymax,xmin:xmax],hh,overwrite=True)

        return output

    def coordCutCube(self,filename,rap1,decp1,rap2,decp2):

   
        rap1 = cP.hms2deg(rap1)
        rap2 = cP.hms2deg(rap2)
        decp1 = cP.dms2deg(decp1)
        decp2 = cP.dms2deg(decp2)
        hh,dd = hP.cleanHead(filename,writeFile=False)
        dC=fits.getdata(filename)
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
        hc=fits.getheader(filename)
        
        hc['NAXIS1']=naxis1
        hc['NAXIS1']=naxis1
        hc['NAXIS2']=naxis2
        hc['CRPIX1']=naxis1/2
        hc['CRPIX2']=naxis2/2
        hc['CRVAL1']=float(raCtr)
        hc['CRVAL2']=float(decCtr)  
             
        aaa = str.split(filename, '.fits')

        output=aaa[0]+'_coordCut.fits'
        print(dC.shape,ymin,ymax,xmin,xmax)
        print(hc)
        fits.writeto(output,dC[:,ymin:ymax,xmin:xmax],hc,overwrite=True)

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
        print(fits.getdata(filename).shape)

        if len(fits.getdata(filename).shape) !=2:
            fits.writeto(output,fits.getdata(filename)[:,ymin:ymax,xmin:xmax],hh,overwrite=True)
        if ychans is not None:
            crval3=hh['CRVAL3']+(hh['CDELT3']*ychans[0])
            hh['CRVAL3']=crval3
            print(crval3)
            fits.writeto(output,fits.getdata(filename)[ychans[0]:ychans[1],ymin:ymax,xmin:xmax],hh,overwrite=True)
        if len(fits.getdata(filename).shape) ==2:
            fits.writeto(output,fits.getdata(filename)[ymin:ymax,xmin:xmax],hh,overwrite=True)
        return output


    def coordToPix(self,imagename,ra,dec,verbose=False):
        '''
        
        Module called by abs_ex
        Converts ra,dec of continuum sources
        into pixel coordinates of the datacube
        
        '''

        #I load the WCS coordinate system:
        #open file

        hdulist = fits.open(imagename)  # read input

        # read data and header
        #what follows works for wcs, but can be written better
        # RS added some additional if clauses
        prihdr = hdulist[0].header
        prihdr, dats = hP.cleanHead(imagename)
        #if prihdr['NAXIS'] == 4:
        # if 'CTYPE4' in prihdr:
        #     del prihdr['CTYPE4']
        # if 'CDELT4' in prihdr:
        #     del prihdr['CDELT4']
        # if 'CRVAL4' in prihdr:
        #     del prihdr['CRVAL4']
        # if 'CRPIX4' in prihdr:
        #     del prihdr['CRPIX4']
        # if 'CUNIT4' in prihdr:
        #     del prihdr['CUNIT4']   
        # if 'NAXIS4' in prihdr:
        #     del prihdr['NAXIS4']
        # if 'CTYPE3' in prihdr:
        #     del prihdr['CTYPE3']
        # if 'CDELT3' in prihdr:
        #     del prihdr['CDELT3']
        # if 'CRVAL3' in prihdr:
        #     del prihdr['CRVAL3']
        # if 'CRPIX3' in prihdr:
        #     del prihdr['CRPIX3'] 
        # if 'NAXIS3' in prihdr:
        #     del prihdr['NAXIS3']
        # if 'CUNIT3' in prihdr:
        #     del prihdr['CUNIT3']

        # del prihdr['NAXIS']
        # prihdr['NAXIS']=2
        w=wcs.WCS(prihdr)    

        pixels=np.zeros([2])
        count_out = 0
        count_flag = 0 
        #for i in range(0,len(ra)):

        ra_deg = cP.hms2deg(ra)
        dec_deg = cP.dms2deg(dec)
        print(ra_deg,dec_deg)
        #px,py=w.wcs_world2pix(ra_deg,dec_deg,0)
        
        px,py=w.wcs_world2pix(ra_deg,dec_deg,0)
        print(px,py)
        if (0 < np.round(px) < prihdr['NAXIS1'] and
                0 < np.round(py) < prihdr['NAXIS2']): 
            pixels[0]= np.round(px)
            pixels[1]= np.round(py)
        else :
            pixels[0]= np.nan
            pixels[1]= np.nan
            count_out +=1
            if verbose == True:
                print('# Source # '+str([i])+ ' lies outside the fov of the data cube #')

        #print('# Total number of sources: \t'+str(len(ra)))
        #print('# Sources outside f.o.v.:\t'+str(count_out))
        #print('# Sources to analyze: \t\t'+str(len(ra)-count_flag-count_out))

        return pixels

    def takeExt(self,fitsName,ext=1):

        f=fits.open('fitsName')
        d=f[ext].data
        h=f[ext].to_header

        newName=fitsName.replace('.fits','_ext'+str(ext)+'.fits')

        fits.writeto(newName,d,h,overwrite=True)

        return newName

    def to32Bits(self,fileName):

        ff=fits.open(fileName)
        dd=ff[0].data
        hh=ff[0].header
        dd=dd.astype('float32')

        outfile=fileName.replace('.fits','_bt32.fits')
        fits.writeto(outfile,dd,hh,overwrite=True)

        return 0

    def regridmaskfits(self, galaxy, gal1, gal2, rob):


        beamsize = {0.0: [2, 1800], 0.5: [2.5, 1440], 1.5: [3, 1200]}
        maskhdu = astfit.open("{datapath}{galaxy}/cleanmask/{gal1}_{gal2}_cleanmask_r15_mask.fits".format(datapath=datapath, galaxy=galaxy, gal1=gal1, gal2=gal2))

        nhdr = (maskhdu[0].header).copy()
        nhdr['NAXIS1'] = beamsize[rob][1]
        nhdr['NAXIS2'] = beamsize[rob][1]
        nhdr['CRPIX1'] = int(beamsize[rob][1]/2)
        nhdr['CRPIX2'] = int(beamsize[rob][1]/2)
        nhdr['CDELT1'] = -beamsize[rob][0]/3600.
        nhdr['CDELT2'] = beamsize[rob][0]/3600.
        nwcs = astwcs.WCS(nhdr).sub(2)
        wcs = astwcs.WCS(maskhdu[0].header).sub(2)

        mask = maskhdu[0].data
        mask[ mask> 0] = 1
        nmask = np.zeros((nhdr['NAXIS3'], beamsize[rob][1],beamsize[rob][1] ))

        for i in range(nhdr['NAXIS3']):
            nmask[i], footprint = reproject.reproject_interp((mask[i], wcs), output_projection=nwcs, shape_out=(beamsize[rob][1],beamsize[rob][1]))

        nmask[nmask >= 0.5] = 1
        nmask[nmask < 0.5] = 0

        nhdu = astfit.HDUList([astfit.PrimaryHDU(nmask, header=nhdr)])
        nhdu.writeto('{datapath}{galaxy}/galx/output/masking/{gal1}_{gal2}_cleanmask_r15_mask_{npix}.fits'.format(datapath=datapath, galaxy=galaxy, gal1=gal1, gal2=gal2, npix=beamsize[rob][1]), overwrite=True)
        return




