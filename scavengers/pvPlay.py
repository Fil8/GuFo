#!/usr/bin/env python3.6
'''

Set of tools to draw pv-diagrams and plot them.

Requirements
------------
Datacube must exists.
'''
import os, sys
from astropy.io import fits
import numpy as np

from math import floor,ceil


from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import gridspec
from matplotlib import patches as mpatches
from matplotlib import colorbar
from matplotlib.patches import Rectangle, Ellipse
from matplotlib import colors
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator, FormatStrFormatter, ScalarFormatter
from matplotlib import transforms as mtransforms
from matplotlib.ticker import LogFormatter 
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import matplotlib.cm as cm

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import astropy.visualization as astviz
from reproject import reproject_interp, reproject_exact

from pvextractor import PathFromCenter
from spectral_cube import SpectralCube
from pvextractor import extract_pv_slice
from scavengers import cvPlay,tPlay, cubePlay
from scavengers import util as ut

cvP = cvPlay.convert()
tP = tPlay.tplay()
cP = cubePlay.cubeplay()

class pvplay(object):
    '''This class makes pv-diagrams and plots them.

    '''

    def pvCut(self,cfg_par,pa,width):
        
        objCoordsRA = cfg_par['pvDiagram']['centreRA']
        objCoordsDec = cfg_par['pvDiagram']['centreDec']
        objCoordsRA=cvP.hms2deg(objCoordsRA)
        objCoordsDec=cvP.dms2deg(objCoordsDec)
        centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree,frame='fk5')

        pvPath = PathFromCenter(center=centre,length=cfg_par['pvDiagram']['length'] * u.arcmin,angle=pa * u.deg,width=width * u.arcsec)
        #print(width,cfg_par['pvDiagram']['length'])
        #print(pvPath.values)
        cubeName=cfg_par['pvDiagram']['cubeName']
        cube = SpectralCube.read(cfg_par['pvDiagram']['cubeDir']+cubeName)
        print(cubeName)

        cube=fits.getdata(cfg_par['pvDiagram']['cubeDir']+cubeName)        
        header=fits.getheader(cfg_par['pvDiagram']['cubeDir']+cubeName)
        print(header['CDELT1'],header["CDELT2"])
        # del header['PC2_1']
        # del header['PC2_2']
        # del header['PC1_2']
        # del header['PC1_1']
        header['PC3_1']=0.0
        header['PC3_2']=0.0
        header['PC3_3']=1.0
        cubeWCS=WCS(header)

        print(cubeWCS)
        slice1 = extract_pv_slice(cube, pvPath,wcs=cubeWCS)  
        outName=cubeName.replace('.fits','_pv'+str(pa)+'_'+str(int(width))+'_asec.fits')
        print(outName)
        outPv = cfg_par['pvDiagram']['pvDir']+outName
        slice1.writeto(outPv,overwrite=True)  


        header = fits.getheader(outPv)
        header["CRVAL1"] = 0
        header["CRPIX1"] = header['NAXIS1']/2
        header["CUNIT2"] = 'm/s'

        fits.writeto(cfg_par['pvDiagram']['pvDir']+outName,fits.getdata(outPv),header,overwrite=True)
        return cfg_par['pvDiagram']['pvDir']+outName
    
    def pvCutCentreCoords(self,cfg_par,centre):
        centre = SkyCoord(ra=centre[0]*u.degree, dec=centre[1]*u.degree, frame='fk5')
        centreX = centre.ra
        centreY = centre.dec

        pvPath = PathFromCenter(center=centre,length=cfg_par['pvDiagram']['length'] * u.arcmin,angle=cfg_par['pvDiagram']['pa'] * u.deg,width=cfg_par['pvDiagram']['width'] * u.arcsec,sample=10)
        cubeName=cfg_par['pvDiagram']['cubeName']
        cube = SpectralCube.read(cfg_par['pvDiagram']['cubeDir']+cubeName)
        slice1 = extract_pv_slice(cube, pvPath)  
        
        header = slice1.header
        header["CRVAL1"] = 0
        header["CRPIX1"] = header['NAXIS1']/2
        header["CUNIT2"] = 'm/s'
        outName=cfg_par['general']['outPrefix'][0]+'.fits'
        fits.writeto(cfg_par['pvDiagram']['pvDir']+outName,slice1.data,header,overwrite=True)
        
        return str(cfg_par['pvDiagram']['pvDir']+outName)



    def circPoint(self,alpha,beta,m,q,pa,ww,keyAlpha):
        '''
        module called by findExtremes

        '''

        a = -2.*alpha
        b = -2.*beta 
        return alpha,beta
        # if ww != 0.:
        #     r = (ww*1.5)/np.tan(np.radians(float(pa)))
        # else:
        #     return alpha,beta
        # c = np.power(alpha,2)+np.power(beta,2)-np.power(r,2)

        # bb = 2.*m*q+a+b*m
        # delta = np.sqrt(np.power(bb,2)-4.*(1.+np.power(m,2))*(np.power(q,2)+b*q+c))
        # dueA = 2.*(1.+np.power(m,2))
        # x = np.divide(-bb+delta,dueA)
        # if keyAlpha == 'min':
        #     if x<alpha:
        #         x = np.divide(-bb-delta,dueA)
        # elif keyAlpha == 'max':
        #     if x>alpha:
        #         x = np.divide(-bb-delta,dueA)
        # y = m*x+q
        # return x,y

    def findExtremes(self,ww,pa,Xc,Yc,shape):
        '''
        module called by pvDiagram

        '''
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

    def pvDiagram(self,cfg_par):
        '''
        adapted from from SoFiA (https://github.com/SoFiA-Admin/SoFiA).

        Draws a pv-diagram and saves it to a .fits file. Centre (HH:MM:SS, DD:MM:SS), PA (N is up, E is left) and width (in arcseconds) 
        of the cut must be provided in the configuration file.
        
        Parameters
        ----------

        cfg_par['pvDiagram']: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.
        
        Returns
        -------

        name: str
            full path to pv-diagram

        Requirements:
        -------------

        cube must exist, thickness of cut must be lowered if error is given (problem with extremes).

        '''
        inCubelet=cfg_par['pvDiagram']['inCube']
        f = fits.open(inCubelet)
        data = f[0].data
        headerCubelets = f[0].header

        # Centres and bounding boxes
        Xc = cfg_par['pvDiagram']['raCentre']
        Yc = cfg_par['pvDiagram']['decCentre']
        Xc = cvP.hms2deg(Xc)
        Yc = cvP.dms2deg(Yc)
        header2d = headerCubelets.copy()
        cP.delete_3rd_axis(header2d)

        pa = cfg_par['pvDiagram']['pa']+90.
        kin_pa = np.radians(float(pa))

        w = WCS(header2d)

        Xc,Yc=w.wcs_world2pix(Xc,Yc,0)

        print(Xc,Yc)

        ww = np.rint(cfg_par['pvDiagram']['width']/3600./headerCubelets['CDELT2'])
        print(ww)
        if ww < 1 :     
            ww = 0.0
        
        if float(pa) == 90.:
            extrs = [[Xc,0],[Xc,data.shape[1]]]            
        else:

            extrs,m = self.findExtremes(ww,pa,Xc,Yc,data.shape)
        # if cfg_par['cubePlay']['cubelets']['raMin'] != False:
        #     Xmin = cfg_par['cubePlay']['cubelets']['raMin']
        #     Ymin = cfg_par['cubePlay']['cubelets']['decMin']
        #     Xmax = cfg_par['cubePlay']['cubelets']['raMax']
        #     Ymax = cfg_par['cubePlay']['cubelets']['decMax']
        # else:
        #     Xmin = 0
        #     Ymin = 0
        #     Xmax = subcube.shape[2]
        #     Ymax = subcube.shape[1]
        #sys.exit(0)
        print(Xc,Yc)
        N=5.
        px,py,dt,tmp,ss,tmpp=extrs[0][0],extrs[0][1],1./N,np.arange(0),np.arange(0),np.arange(0)
        Dx,Dy=extrs[1][0]-extrs[0][0],extrs[1][1]-extrs[0][1]
        dist=np.sqrt((Dx)**2+(Dy)**2)
        print(np.rint((px)) , np.rint(Xc) , np.rint((py)) , np.rint(Yc),dist)
        print(dt*(extrs[1][0]-extrs[0][0])/dist)

        # find perpendicular
        if ww and ww > 1:
            if Dx:
                Dyp=ww/np.sqrt(1+Dy**2/Dx**2)
                Dxp=-Dyp*Dy/Dx
            else:
                Dxp=ww/np.sqrt(1+Dx**2/Dy**2)
                Dyp=-Dxp*Dx/Dy
            distp=np.sqrt((Dxp)**2+(Dyp)**2)
        else: Dxp,Dyp,distp=0,0,1
        # move from one extreme to the other of the line
        # append 1 average spectrum 'tmp' per pixel to 'ss' sampling the pixel with N sub-pixels
        if m >= 0:
            while np.rint((px))<extrs[1][0] and np.rint((py))<extrs[1][1]:
                dw=-ww

                while dw<=ww+1e-9:
                    pxp,pyp=(px+dw*Dxp/distp),(py+dw*Dyp/distp)
                    #sys.exit(0)
                    if 0. <= pxp <= data.shape[2] and 0.<= pyp <=data.shape[1]:
                        if not tmpp.sum():
                            tmpp=data[:,int(floor(pyp)),int(floor(pxp))]
                            tmpp=np.atleast_2d(tmpp)
                        else: 
                            tmpp=np.vstack((tmpp,data[:,int(floor(pyp)),int(floor(pxp))]))
                    dw+=1./N
                if not tmp.sum():
                    tmp=tmpp.mean(axis=0)
                    tmp=np.atleast_2d(tmp)
                else: tmp=np.vstack((tmp,tmpp.mean(axis=0)))
                tmpp=np.arange(0)
                if tmp.shape[0]==N:
                    if not ss.sum():
                        ss=tmp.mean(axis=0)
                        ss=np.atleast_2d(ss)
                    else: ss=np.vstack((ss,tmp.mean(axis=0)))
                    tmp=np.arange(0)
                #print(ss.shape)
                if np.rint((px)) == np.rint(Xc) and np.rint((py)) == np.rint(Yc) :
                    centreX = ss.shape[0]
                px+=dt*(extrs[1][0]-extrs[0][0])/dist
                py+=dt*(extrs[1][1]-extrs[0][1])/dist
        elif m<0:
            while np.rint(floor(px))<extrs[1][0] and np.rint(floor(py))>extrs[1][1]:
                dw=-ww

                while dw<=ww+1e-9:
                    pxp,pyp=(px+dw*Dxp/distp),(py+dw*Dyp/distp)
                    
                    #sys.exit(0)
                    if not tmpp.sum():
                        tmpp=data[:,int(floor(pyp)),int(floor(pxp))]
                        tmpp=np.atleast_2d(tmpp)
                    else: tmpp=np.vstack((tmpp,data[:,int(floor(pyp)),int(floor(pxp))]))
                    dw+=1./N
                if not tmp.sum():
                    tmp=tmpp.mean(axis=0)
                    tmp=np.atleast_2d(tmp)
                else: tmp=np.vstack((tmp,tmpp.mean(axis=0)))
                tmpp=np.arange(0)
                if tmp.shape[0]==N:
                    if not ss.sum():
                        ss=tmp.mean(axis=0)
                        ss=np.atleast_2d(ss)
                    else: ss=np.vstack((ss,tmp.mean(axis=0)))
                    tmp=np.arange(0)
                if np.rint((px)) == np.rint(Xc) and np.rint((py)) == np.rint(Yc) :
                    centreX = ss.shape[0]
                px+=dt*(extrs[1][0]-extrs[0][0])/dist
                py+=dt*(extrs[1][1]-extrs[0][1])/dist            
            #print(pv_slice)

        # vel = ((np.linspace(1, ss.shape[0], ss.shape[0]) - headerCubelets['CRPIX3']) 
        #     * headerCubelets['CDELT3'] + headerCubelets['CRVAL3'])

        # print(vel[-1],vel[0],headerCubelets['CRPIX3'])
        # # print(vel[0],vel[-1])
        # if np.float(vel[0]) > np.float(vel[-1]):
        #     headerCubelets['CRPIX3'] = ss.shape[0]
        #     headerCubelets['CRVAL3'] = vel[-1]/1e3
        #     headerCubelets['CDELT3'] = -headerCubelets['CDELT3']
        #     ss = np.fliplr(ss).T
        #     print('ciao')
        # else:
        #     ss = ss.T


        hdu = fits.PrimaryHDU(data=ss.T, header=headerCubelets)
        hdulist = fits.HDUList([hdu])
        # print(hdulist[0].header)
        #if hdulist[0].header["CRPIX3"]
        # print(hdulist[0].header['CRVAL3'])

        hdulist[0].header["CTYPE1"] = "PV--DIST"
        hdulist[0].header["CDELT1"] = hdulist[0].header["CDELT2"]
        hdulist[0].header["CRVAL1"] = 0
        hdulist[0].header["CRPIX1"] = centreX
        hdulist[0].header["CTYPE2"] = hdulist[0].header["CTYPE3"]
        hdulist[0].header["CDELT2"] = hdulist[0].header["CDELT3"]
        hdulist[0].header["CRVAL2"] = hdulist[0].header["CRVAL3"]
        hdulist[0].header["CRPIX2"] = hdulist[0].header["CRPIX3"]
        if "CUNIT3" in hdulist[0].header:
            hdulist[0].header["CUNIT2"] = hdulist[0].header["CUNIT3"]
        else:
            hdulist[0].header["CUNIT2"] = 'm/s'
        hdulist[0].header["ORIGIN"] = 'GaNGiaLF'

        cP.delete_3rd_axis(hdulist[0].header)
        
        outCubelet=str.split(os.path.basename(inCubelet),'.')[0]
        outCubelet=outCubelet+'_'+str(cfg_par['pvDiagram']['pa'])+'-pv.fits'
        name = cfg_par['pvDiagram']['pvDir'] + outCubelet
        hdulist.writeto(name,overwrite=True)
        
        # hdu1 = fits.PrimaryHDU(data=pv_slice1[2,:,:], header=headerCubelets)
        # hdul = fits.HDUList([hdu1])
        # hdul[0].header["CTYPE1"] = "PV--DIST"
        # hdul[0].header["CDELT1"] = hdul[0].header["CDELT2"]
        # hdul[0].header["CRVAL1"] = 0
        # hdul[0].header["CRPIX1"] = pv_slice1.shape[1] / 2
        # hdul[0].header["CTYPE2"] = hdul[0].header["CTYPE3"]
        # hdul[0].header["CDELT2"] = hdul[0].header["CDELT3"]
        # hdul[0].header["CRVAL2"] = hdul[0].header["CRVAL3"]
        # hdul[0].header["CRPIX2"] = hdul[0].header["CRPIX3"]
        # hdul[0].header["ORIGIN"] = 'GaNGiaLF'

        # self.delete_3rd_axis(hdul[0].header)
        
        # outCubelet=str.split(os.path.basename(inCubelet),'.')[0]
        # outCubelet=outCubelet+'_'+str(cfg_par['cubePlay']['pvDiagram']['pa'])+'-pvSingle.fits'
        # name = cfg_par['general']['pvDir'] + outCubelet
        # print(name)
        # hdul.writeto(name,overwrite=True)

        return name

    def pvPlot(self,cfg_par,imName,cMap,
        imLevels=None,imColors=None,vsys=None,velRange=None,cRange=None,nanToZero=None,
        zeroToNan=None,interpMethod=None,cScale='linear',pvUnit='mJy',linthresh=None,base=10.,imContours=True):
        '''Draws a pv-diagram. Multiple contours are overlaid if more than one image is given.

        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        imName: list
            full paths to pvdiagrams

        cRange: np.array, optional
            [2,2] array with the min max range of the two colorscales

        imLevels: np.array, optional
            [2,N] array with N contours to plot (referring to im1, im2 or imMom0, if given)

        imColors: list, optional
            _default=['black','black']_ color of contours

        vsys: float, optional,
            if given, an horizontal dashed line is drawn at the systemic velocity of the source

        velRange: np.array, optional
            [2, N] array with min and max ranges of x and y axes.

        nanToZero: bool, optional
             _default=False_: converts nans values to zeros 

        zeroToNan: bool, optional
             _default=False_: converts zero values to np.nan 

        interpMethod: str, optional
            kind of interpolation for `matplotlib.imshow`

        title: list, optional
            title of plot

        cScale: str, optional
            _default=linear_: plot in 'linear', 'log' or 'sqrt' scale. If 'log' scale is used base and linthresh must be given.

        linthresh: float, optional
            _default=0.1_: modulates the logarithmic colorscale

        base: float, optional
            _default=10._: base of the logarithmic scale

        imContours: bool
            _default=True_: draw contours given in imLevels

        Returns
        ----------
        outFig: str
            full path to output figure

        Notes
        ----------
        Useful for plotting moment maps of different phases of the gas.

        '''

        params = ut.loadRcParams()
        plt.rcParams.update(params)
        hduIm = fits.open(imName[0])[0]
        wcsIm = WCS(hduIm.header)

        if nanToZero is not False:
            index=np.isnan(hduIm.data)
            hduIm.data[index] = 0.0
            mapName=nanToZero
        elif zeroToNan is not False:
            dd = np.array(hduIm.data,dtype=float)
            index=np.where(dd==0.)
            dd[index] = np.nan
            mapName=zeroToNan
            hduIm.data= dd
        

        hduIm.data*=1e3
        
        vel = ((np.linspace(1, hduIm.data.shape[0], hduIm.data.shape[0]) - hduIm.header['CRPIX2']) 
            * hduIm.header['CDELT2'] + hduIm.header['CRVAL2'])/1e3
        #print(vel)
        #print(vel)
        #print(hduIm.header)
        xS = ((np.linspace(1, hduIm.data.shape[1], hduIm.data.shape[1]) - hduIm.header['CRPIX1']) 
            * hduIm.header['CDELT1'] + hduIm.header['CRVAL1'])
        if velRange[1] is not None:
            ext_ymin=np.where(abs(np.asarray(vel,dtype=int)-velRange[1,0])==abs(np.asarray(vel,dtype=int)-velRange[1,0]).min())[0][0]
            ext_ymax =np.where(abs(np.asarray(vel,dtype=int)-velRange[1,1])==abs(np.asarray(vel,dtype=int)-velRange[1,1]).min())[0][0]
        else:
            ext_ymin = 0
            ext_ymax = hduIm.data.shape[0]
        if ext_ymin>ext_ymax:
            ext_ymin_tmp=ext_ymax
            ext_ymax=ext_ymin
            ext_ymin=ext_ymin_tmp
        if velRange[0] is not None:
            ext_xmin=np.where(abs(xS-velRange[0,0])==abs(xS-velRange[0,0]).min())[0][0]
            ext_xmax = np.where(abs(xS-velRange[0,1])==abs(xS-velRange[0,1]).min())[0][0]
        else:
            ext_xmin = 0
            ext_xmax = hduIm.data.shape[1]
        fig = plt.figure(figsize=(7.24409,4.074800625),constrained_layout=False)
        fig.set_tight_layout(False)

        ax1 = plt.subplot(projection=wcsIm)    

        divider = make_axes_locatable(ax1)

        current_cmap = cm.get_cmap(cMap)
        
        if nanToZero is not None or zeroToNan is not None:
            current_cmap.set_bad(color=mapName)
        

        if cScale == 'linear':
            if cRange is None:
                cRangeMin=np.nanmin(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax])
                cRangeMax=np.nanmax(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax])
            else:
                cRangeMin = cRange[0]
                cRangeMax = cRange[1]
            #cRangeMin=np.nanmin(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax])
            print('crange')
            print(cRangeMin,cRangeMax)
            norm = astviz.ImageNormalize(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], vmin=np.nanpercentile(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], 0.0), vmax=np.nanpercentile(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], 0.1))
            
            # img = ax1.imshow(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax],norm=norm, cmap=current_cmap,origin='lower',
            #     interpolation=interpMethod,vmin=cRangeMin,vmax=cRangeMax,aspect='auto',
            #     extent=[ext_xmin-np.abs(hduIm.header['CDELT1']/2),ext_xmax+np.abs(hduIm.header['CDELT1']/2),
            #     ext_ymin-np.abs(hduIm.header['CDELT1']/2),ext_ymax+np.abs(hduIm.header['CDELT1']/2)]) 


            img = ax1.imshow(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax],norm=norm, cmap=current_cmap,origin='lower',
                interpolation=interpMethod,aspect='auto',
                extent=[ext_xmin-np.abs(hduIm.header['CDELT1']/2),ext_xmax+np.abs(hduIm.header['CDELT1']/2),
                ext_ymin-np.abs(hduIm.header['CDELT1']/2),ext_ymax+np.abs(hduIm.header['CDELT1']/2)]) 




        elif cScale == 'sqrt':
            img = ax1.imshow(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], cmap=current_cmap,vmin=cRange[0],vmax=cRange[1])
        elif cScale == 'log':
            img = ax1.imshow(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], cmap=current_cmap,norm=SymLogNorm(linthresh=linthresh,linscale=1.,
                vmin=cRange[0],vmax=cRange[1],base=base),interpolation=interpMethod)

        if imLevels is not None and imContours==True:


            csPos = ax1.contour(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax],levels=imLevels[0,~np.isnan(imLevels[0,:,0]),0], 
                colors=imColors[0,0],linestyles = '-',linewidths=2.,
                origin='lower',aspect='auto',extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax],interpMethod='nearest')
            if np.nansum(imLevels[0,~np.isnan(imLevels[0,:,1]),1])!=0.0:
                print(imLevels[0,~np.isnan(imLevels[0,:,1]),1])
                csNeg = ax1.contour(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax],levels=imLevels[0,~np.isnan(imLevels[0,:,1]),1],
                    colors=imColors[0,1],linestyles = 'dashed',linewidths=1.5,
                    origin='lower',aspect='auto',extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])
            # if contValues[0]==1:
            #     ax1.clabel(csPos, inline=1, fontsize=14)
            #     ax1.clabel(csNeg, inline=1, fontsize=14)

        # if np.sum(cTicks) == 0.0:
        cRange=[np.nanpercentile(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], 5.),np.nanpercentile(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], 95.)]


        cTicks = np.linspace(float(cRange[0]),float(cRange[1]),5)    
        ax1.coords[1].get_format_unit()
        ax1.coords[1].set_format_unit(u.km/u.s)

        ax1.coords[1].set_axislabel(r'Velocity [km s$^{-1}$]')
        ax1.coords[0].set_axislabel(r'Radial offset [deg]')
        

        if len(imName) ==1:       
            cax = divider.append_axes("right", size='2%', pad=0.1)
            
            cbar = plt.colorbar(img, cax=cax,ticks =cTicks,
                        orientation='vertical', format='%.3e')   
            cax.coords[0].grid(False)
            cax.coords[1].grid(False)
            cax.tick_params(direction='in')
            cax.coords[0].tick_params(top=False, bottom=False,
                           labeltop=False, labelbottom=False,style='sci')
            cax.coords[1].set_ticklabel_position('r')
            if pvUnit == 'mJy':
                cBarLabel= r'mJy beam$^{-1}$'
            elif pvUnit == 'nhi':
                cBarLabel= r'$\times 10^{18}$ cm$^{-2}$'
            cax.coords[1].set_axislabel(cBarLabel)
            cax.coords[1].set_axislabel_position('r')        

        # if titleName is not None:
        #     ax1.set_title(titleName)
        
        if vsys is not None:
            # if corrAxes==None:
            #     corrAxes=[0,0]
            nearest_idx = np.where(abs(vel-cfg_par['pvDiagram']['pvPlots']['vsys'])==abs(vel-cfg_par['pvDiagram']['pvPlots']['vsys']).min())[0][0]
            ax1.axhline(y=nearest_idx,color='black',lw=1,ls='-.')
            
        nearest_idx = np.where(abs(xS-0.0)==abs(xS-0.0).min())[0][0]

        ax1.axvline(x=nearest_idx,color='black',lw=1,ls='-.')        

#        ax1.set_autoscale_on(False)    
 
        #SaveOutput
 
        if len(imName)>1:
            for i in range(1,len(imName)):
                hduCont = fits.open(imName[i])[0]
                wcsCont = WCS(hduCont.header)
                # hduCont.data*=1e3

                vell = ((np.linspace(1, hduCont.data.shape[0], hduCont.data.shape[0]) - hduCont.header['CRPIX2'])* hduCont.header['CDELT2'] + hduCont.header['CRVAL2'])/1e3
        #print(vel)
        #print(vel)
            #if i==3:
                #vell=np.flipud(vell)
                #print(vel,hduCont.header['CDELT2'],hduCont.header['CRVAL2'])
                xSS = ((np.linspace(1, hduCont.data.shape[1], hduCont.data.shape[1]) - hduCont.header['CRPIX1']) 
                    * hduCont.header['CDELT1'] + hduCont.header['CRVAL1'])
                print(xSS )


                if velRange[1] is not None:
                    ext_yminCont=np.where(abs(np.asarray(vell,dtype=int)-velRange[1,0])==abs(np.asarray(vell,dtype=int)-velRange[1,0]).min())[0][0]
                    ext_ymaxCont =np.where(abs(np.asarray(vell,dtype=int)-velRange[1,1])==abs(np.asarray(vell,dtype=int)-velRange[1,1]).min())[0][0]
                else:
                    ext_yminCont = 0
                    ext_ymaxCont = hduIm.data.shape[0]
                
                if ext_yminCont>ext_ymaxCont:
                    ext_ymin_tmp=ext_ymaxCont
                    ext_ymaxCont=ext_yminCont
                    ext_yminCont=ext_ymin_tmp
                
                if velRange[0] is not None:
                    ext_xminCont=np.where(abs(xSS-velRange[0,0])==abs(xSS-velRange[0,0]).min())[0][0]
                    ext_xmaxCont = np.where(abs(xSS-velRange[0,1])==abs(xSS-velRange[0,1]).min())[0][0]
                    print(velRange[0,0],ext_xminCont,ext_xmaxCont)
                    print('here')
                else:
                    ext_xminCont = 0
                    ext_xmaxCont = hduIm.data.shape[1]       

                # if velRange[0] is not None:
                #     ext_xminCont=np.where(abs(xSS-velRange[0,0])==abs(xSS-velRange[0,0]).min())[0][0]
                #     ext_xmaxCont = np.where(abs(xSS-velRange[0,1])==abs(xSS-velRange[0,1]).min())[0][0]
                #     ext_xminIm=np.where(abs(xS-velRange[0,0])==abs(xS-velRange[0,0]).min())[0][0]
                #     ext_xmaxIm = np.where(abs(xS-velRange[0,1])==abs(xS-velRange[0,1]).min())[0][0]
                #     if  i==2:
                #         ext_xminIm=np.where(abs(xS-xSS[0])==abs(xS-xSS[0]).min())[0][0]-30
                #         ext_xmaxIm = np.where(abs(xS-xSS[-1])==abs(xS-xSS[-1]).min())[0][0]-30
                #     else:
                #         ext_xminIm=ext_xmin
                #         ext_xmaxIm =ext_xmax                     
                # else:
                # ext_xminCont = 0
                # ext_xmaxCont = hduCont.data.shape[1]
                # ext_xminIm=ext_xmin
                # ext_xmaxIm =ext_xmax
                
                #print(ext_xminCont,ext_xmaxCont)
                #print(ext_yminCont,ext_ymaxCont)

                # size= (ext_ymaxCont-ext_yminCont,ext_xmaxCont-ext_xminCont)
                # centre = (ext_yminCont+size[1]/2.,ext_xminCont+size[0]/2.)
                # print(size,centre)
                # hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont, mode='partial')   
                # hduContCut.data *=1e3 

                # size1= (ext_ymax-ext_ymin,ext_xmax-ext_xmin)
                # centre1= (ext_ymin+size1[1]/2.,ext_xmin+size1[0]/2.)
                # print(hduIm.shape)
                # print(centre1,size1,hduIm.data.shape)
                # hduImCut = Cutout2D(hduIm.data, centre1, size1, wcs=wcsIm, mode='partial') 
                

                # array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
                #                             hduImCut.wcs,shape_out=hduImCut.shape)
                # print(array.data.shape)
                #print(hduCont.data[ext_yminCont:ext_ymaxCont,ext_xminCont:ext_xmaxCont].shape)
                csPos = ax1.contour(hduCont.data[ext_yminCont:ext_ymaxCont,ext_xminCont:ext_xmaxCont]*1e3,levels=imLevels[i,~np.isnan(imLevels[i,:,0]),0],
                    colors=imColors[i,0],linestyles = '-',linewidths=1,
                    origin='lower',aspect='auto',extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])
                print(csPos)
                
                if np.nansum(imLevels[0,~np.isnan(imLevels[0,:,1]),1]) != np.nan:
                    print('NEG')
                    csNeg = ax1.contour(hduCont.data[ext_yminCont:ext_ymaxCont,ext_xminCont:ext_xmaxCont]*1e3,levels=imLevels[i,~np.isnan(imLevels[i,:,1]),1], 
                        colors=imColors[i,1],linestyles = 'dashed',linewidths=1,origin='lower',aspect='auto',extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])   
                
                # if contValues[i]==1:
                #     ax1.clabel(csPos, inline=1, fontsize=14)
                #     ax1.clabel(csNeg, inline=1, fontsize=14)
        ax1.set_autoscale_on(False)    
        plotFormat= cfg_par['pvDiagram']['pvPlots']['plotFormat']
        if len(imName)==1:
            baseName = os.path.basename(imName[0])
        else: 
            baseName=cfg_par['general']['outPrefix'][0]+'pvPlot.fits'
        outFigName= cfg_par['pvDiagram']['pvPlots']['plotDir']+'/'+baseName.replace('.fits','.'+plotFormat)
        print(outFigName)
        fig.savefig(outFigName,format=plotFormat, bbox_inches = "tight",overwrite=True,dpi=300,transparent=False)#,
                    #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)

        return outFigName



    def pvMultiPlot(self,cfg_par,imName,cMap,
        imLevels=None,imColors=None,vsys=None,velRange=None,cRange=None,nanToZero=None,
        zeroToNan=None,interpMethod=None,cScale='linear',linthresh=None,base=10.,imContours=True):

        '''Draws a pv-diagram. Multiple contours are overlaid if more than one image is given.

        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        imName: list
            full paths to pvdiagrams

        cRange: np.array, optional
            [2,2] array with the min max range of the two colorscales

        imLevels: np.array, optional
            [2,N] array with N contours to plot (referring to im1, im2 or imMom0, if given)

        imColors: list, optional
            _default=['black','black']_ color of contours

        vsys: float, optional,
            if given, an horizontal dashed line is drawn at the systemic velocity of the source

        velRange: np.array, optional
            [2, N] array with min and max ranges of x and y axes.

        nanToZero: bool, optional
             _default=False_: converts nans values to zeros 

        zeroToNan: bool, optional
             _default=False_: converts zero values to np.nan 

        interpMethod: str, optional
            kind of interpolation for `matplotlib.imshow`

        title: list, optional
            title of plot

        cScale: str, optional
            _default=linear_: plot in 'linear', 'log' or 'sqrt' scale. If 'log' scale is used base and linthresh must be given.

        linthresh: float, optional
            _default=0.1_: modulates the logarithmic colorscale

        base: float, optional
            _default=10._: base of the logarithmic scale

        imContours: bool
            _default=True_: draw contours given in imLevels

        Returns
        ----------
        outFig: str
            full path to output figure

        Notes
        ----------
        Useful for plotting moment maps of different phases of the gas.

        '''


        numPlots=len(imName)

        n_rows = int(np.ceil(numPlots/2.))
        fig = plt.figure(figsize=(8.25, 11.67), constrained_layout=False)
        fig.tight_layout()

        #fig.set_tight_layout(False)
        fig.subplots_adjust(hspace=0.)

        #gs_top = plt.GridSpec(nrows=n_rows+1, ncols=2,  figure=fig, top=0.95)
        gs_base = plt.GridSpec(nrows=n_rows, ncols=2,  figure=fig, hspace=0.,top=0.99)

        # gs = fig.add_gridspec(nrows=n_rows+1, ncols=3, left=0.05, figsize=(8.25, 11.67))
        #gs = gridspec.GridSpec(nrows=n_rows+1, ncols=3,  figure=fig)
        
        #wave_ax = self.addFullSubplot(cfg_par,fig,gs_top,vel,y,result,noise,i,j,lineInfo,singleVorBinInfo)
        params = ut.loadRcParams()
        plt.rcParams.update(params)
        
        #for plot_count in range(n_plots):
        k=0
        for i in range(0,numPlots):

            if i == 0:
                j = 0
                k = 0
            elif i % 2 == 0:
                j +=1 
                k = 0
            else:
                k = 1
            # if i ==1:
            #     k=2
            # elif i==2:
            #     k=1

            hduIm = fits.open(imName[i])[0]
            wcsIm = WCS(hduIm.header)       
 

            if nanToZero is not False:
                index=np.isnan(hduIm.data)
                hduIm.data[index] = 0.0
                mapName=nanToZero
            elif zeroToNan is not False:
                dd = np.array(hduIm.data,dtype=float)
                index=np.where(dd==0.)
                dd[index] = np.nan
                mapName=zeroToNan
                hduIm.data= dd
            # hduIm.data=1e3
            vel = ((np.linspace(1, hduIm.data.shape[0], hduIm.data.shape[0]) - hduIm.header['CRPIX2']) 
                * hduIm.header['CDELT2'] + hduIm.header['CRVAL2'])/1e3
            #print(vel)
            #print(vel)
            #print(hduIm.header)
            xS = ((np.linspace(1, hduIm.data.shape[1], hduIm.data.shape[1]) - hduIm.header['CRPIX1']) 
                * hduIm.header['CDELT1'] + hduIm.header['CRVAL1'])
            if velRange[1] is not None:
                ext_ymin=np.where(abs(np.asarray(vel,dtype=int)-velRange[1,0])==abs(np.asarray(vel,dtype=int)-velRange[1,0]).min())[0][0]
                ext_ymax =np.where(abs(np.asarray(vel,dtype=int)-velRange[1,1])==abs(np.asarray(vel,dtype=int)-velRange[1,1]).min())[0][0]
            else:
                ext_ymin = 0
                ext_ymax = hduIm.data.shape[0]
            if ext_ymin>ext_ymax:
                ext_ymin_tmp=ext_ymax
                ext_ymax=ext_ymin
                ext_ymin=ext_ymin_tmp
            if velRange[0] is not None:
                ext_xmin=np.where(abs(xS-velRange[0,0])==abs(xS-velRange[0,0]).min())[0][0]
                ext_xmax = np.where(abs(xS-velRange[0,1])==abs(xS-velRange[0,1]).min())[0][0]
            else:
                ext_xmin = 0
                ext_xmax = hduIm.data.shape[1]
            print(j,k)
            ax1 = fig.add_subplot(gs_base[j,k],projection=wcsIm)


            divider = make_axes_locatable(ax1)

            current_cmap = cm.get_cmap(cMap)
            
            if nanToZero is not None or zeroToNan is not None:
                current_cmap.set_bad(color=mapName)
            

            if cScale == 'linear':
                if cRange is None:
                    cRangeMin=np.nanmin(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax])
                    cRangeMax=np.nanmax(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax])
                else:
                    cRangeMin = cRange[0]
                    cRangeMax = cRange[1]
                # norm = astviz.ImageNormalize(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], vmin=np.nanpercentile(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], 2.5.), vmax=np.nanpercentile(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], 98.))
                img = ax1.imshow(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], cmap=current_cmap,origin='lower',interpolation=interpMethod,vmin=cRangeMin,vmax=cRangeMax,aspect='auto',extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax]) 
            elif cScale == 'sqrt':
                img = ax1.imshow(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], cmap=current_cmap,vmin=cRange[0],vmax=cRange[1],norm=PowerNorm(gamma=0.5))
            elif cScale == 'log':
                img = ax1.imshow(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], cmap=current_cmap,norm=SymLogNorm(linthresh=linthresh,linscale=1.,
                    vmin=cRange[0],vmax=cRange[1],base=base),interpolation=interpMethod)

            if imLevels is not None and imContours==True:


                csPos = ax1.contour(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax],levels=imLevels[0,~np.isnan(imLevels[0,:,0]),0], 
                    colors=imColors[0,0],linestyles = '-',linewidths=1.,
                    origin='lower',aspect='auto',extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax],interpMethod='nearest')
                if np.nansum(imLevels[0,~np.isnan(imLevels[0,:,1]),1])!=0.0:
                    print(imLevels[0,~np.isnan(imLevels[0,:,1]),1])
                    csNeg = ax1.contour(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax],levels=imLevels[0,~np.isnan(imLevels[0,:,1]),1],
                        colors=imColors[0,1],linestyles = 'dashed',linewidths=1.5,
                        origin='lower',aspect='auto',extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])
                # if contValues[0]==1:
                #     ax1.clabel(csPos, inline=1, fontsize=14)
                #     ax1.clabel(csNeg, inline=1, fontsize=14)

            # if np.sum(cTicks) == 0.0:
            # cRange=[np.nanpercentile(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], 5.),np.nanpercentile(hduIm.data[ext_ymin:ext_ymax,ext_xmin:ext_xmax], 95.)]
            

            cTicks = np.linspace(float(cRange[0]),float(cRange[1]),5)    
            ax1.coords[1].get_format_unit()
            ax1.coords[1].set_format_unit(u.km/u.s)

            ax1.coords[1].set_axislabel(r'Velocity [km s$^{-1}$]')
            ax1.coords[0].set_axislabel(r'Radial offset [deg]')
            

            if len(imName[i]) ==1:       
                cax = divider.append_axes("right", size='2%', pad=0.1)
                
                cbar = plt.colorbar(img, cax=cax,ticks =cTicks,
                            orientation='vertical', format='%.3e')   
                cax.coords[0].grid(False)
                cax.coords[1].grid(False)
                cax.tick_params(direction='in')
                cax.coords[0].tick_params(top=False, bottom=False,
                               labeltop=False, labelbottom=False,style='sci')
                cax.coords[1].set_ticklabel_position('r')
                cBarLabel= r'$\times 10^{-3}$ mJy beam$^{-1}$'
                cax.coords[1].set_axislabel(cBarLabel)
                cax.coords[1].set_axislabel_position('r')        

            # if titleName is not None:
            #     ax1.set_title(titleName)
            
            if vsys is not None:
                # if corrAxes==None:
                #     corrAxes=[0,0]
                nearest_idx = np.where(abs(vel-cfg_par['pvDiagram']['pvPlots']['vsys'])==abs(vel-cfg_par['pvDiagram']['pvPlots']['vsys']).min())[0][0]
                ax1.axhline(y=nearest_idx,color='black',lw=1,ls='-.')
                nearest_idx = np.where(abs(xS-0.0)==abs(xS-0.0).min())[0][0]

                ax1.axvline(x=nearest_idx,color='black',lw=1,ls='-.')        

    #        ax1.set_autoscale_on(False)    
     
            #SaveOutput
     
            # if len(imName[i])>1:
            #     for ii in range(1,len(imName[i])):
            #         hduCont = fits.open(imName[i][ii])[0]
            #         wcsCont = WCS(hduCont.header)
            #         # hduCont.data*=1e3

            #         vell = ((np.linspace(1, hduCont.data.shape[0], hduCont.data.shape[0]) - hduCont.header['CRPIX2'])* hduCont.header['CDELT2'] + hduCont.header['CRVAL2'])/1e3
            # #print(vel)
            # #print(vel)
            #     #if i==3:
            #         #vell=np.flipud(vell)
            #         #print(vel,hduCont.header['CDELT2'],hduCont.header['CRVAL2'])
            #         xSS = ((np.linspace(1, hduCont.data.shape[1], hduCont.data.shape[1]) - hduCont.header['CRPIX1']) 
            #             * hduCont.header['CDELT1'] + hduCont.header['CRVAL1'])

            #         if velRange[1] is not None:
            #             ext_yminCont=np.where(abs(np.asarray(vell,dtype=int)-velRange[1,0])==abs(np.asarray(vell,dtype=int)-velRange[1,0]).min())[0][0]
            #             ext_ymaxCont =np.where(abs(np.asarray(vell,dtype=int)-velRange[1,1])==abs(np.asarray(vell,dtype=int)-velRange[1,1]).min())[0][0]
            #         else:
            #             ext_yminCont = 0
            #             ext_ymaxCont = hduIm.data.shape[0]
                    
            #         if ext_yminCont>ext_ymaxCont:
            #             ext_ymin_tmp=ext_ymaxCont
            #             ext_ymaxCont=ext_yminCont
            #             ext_yminCont=ext_ymin_tmp
                    
            #         if velRange[0] is not None:
            #             print('figa')
            #             ext_xminCont=np.where(abs(xSS-velRange[0,0])==abs(xSS-velRange[0,0]).min())[0][0]
            #             ext_xmaxCont = np.where(abs(xSS-velRange[0,1])==abs(xSS-velRange[0,1]).min())[0][0]
            #         else:
            #             ext_xminCont = 0
            #             ext_xmaxCont = hduIm.data.shape[1]                
            #         # if velRange[0] is not None:
            #         #     ext_xminCont=np.where(abs(xSS-velRange[0,0])==abs(xSS-velRange[0,0]).min())[0][0]
            #         #     ext_xmaxCont = np.where(abs(xSS-velRange[0,1])==abs(xSS-velRange[0,1]).min())[0][0]
            #         #     ext_xminIm=np.where(abs(xS-velRange[0,0])==abs(xS-velRange[0,0]).min())[0][0]
            #         #     ext_xmaxIm = np.where(abs(xS-velRange[0,1])==abs(xS-velRange[0,1]).min())[0][0]
            #         #     if  i==2:
            #         #         ext_xminIm=np.where(abs(xS-xSS[0])==abs(xS-xSS[0]).min())[0][0]-30
            #         #         ext_xmaxIm = np.where(abs(xS-xSS[-1])==abs(xS-xSS[-1]).min())[0][0]-30
            #         #     else:
            #         #         ext_xminIm=ext_xmin
            #         #         ext_xmaxIm =ext_xmax                     
            #         # else:
            #         # ext_xminCont = 0
            #         # ext_xmaxCont = hduCont.data.shape[1]
            #         # ext_xminIm=ext_xmin
            #         # ext_xmaxIm =ext_xmax
                    
            #         #print(ext_xminCont,ext_xmaxCont)
            #         #print(ext_yminCont,ext_ymaxCont)

            #         # size= (ext_ymaxCont-ext_yminCont,ext_xmaxCont-ext_xminCont)
            #         # centre = (ext_yminCont+size[1]/2.,ext_xminCont+size[0]/2.)
            #         # print(size,centre)
            #         # hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont, mode='partial')   
            #         # hduContCut.data *=1e3 

            #         # size1= (ext_ymax-ext_ymin,ext_xmax-ext_xmin)
            #         # centre1= (ext_ymin+size1[1]/2.,ext_xmin+size1[0]/2.)
            #         # print(hduIm.shape)
            #         # print(centre1,size1,hduIm.data.shape)
            #         # hduImCut = Cutout2D(hduIm.data, centre1, size1, wcs=wcsIm, mode='partial') 
                    

            #         # array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
            #         #                             hduImCut.wcs,shape_out=hduImCut.shape)
            #         # print(array.data.shape)
            #         #print(hduCont.data[ext_yminCont:ext_ymaxCont,ext_xminCont:ext_xmaxCont].shape)
            #         csPos = ax1.contour(hduCont.data[ext_yminCont:ext_ymaxCont,ext_xminCont:ext_xmaxCont],levels=imLevels[i,~np.isnan(imLevels[i,:,0]),0],
            #             colors=imColors[i,0],linestyles = '-',linewidths=1,
            #             origin='lower',aspect='auto',extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])
            #         if np.nansum(imLevels[0,~np.isnan(imLevels[0,:,1]),1]) != np.nan:
            #             csNeg = ax1.contour(hduCont.data[ext_yminCont:ext_ymaxCont,ext_xminCont:ext_xmaxCont],levels=imLevels[i,~np.isnan(imLevels[i,:,1]),1], 
            #                 colors=imColors[i,1],linestyles = 'dashed',linewidths=1,origin='lower',aspect='auto',extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])   
                    
            #         # if contValues[i]==1:
            #         #     ax1.clabel(csPos, inline=1, fontsize=14)
            #         #     ax1.clabel(csNeg, inline=1, fontsize=14)
            

            ax1.set_autoscale_on(False)    

            # fig.savefig(outFigName,format=plotFormat, bbox_inches = "tight",overwrite=True,dpi=300,transparent=False)#,
            #             #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)


           #  ax.axvline(color='k', linestyle=':', zorder=0)                           
           #  legend = ax.legend(loc='best',handlelength=0.0, handletextpad=0.0,frameon=False)
           #  legend.get_frame().set_facecolor('none')

 
           #  divider = make_axes_locatable(ax)
           #  ax2 = divider.append_axes("bottom", size='20%',pad=0)
           #  ax2.minorticks_on()

           #  ax.figure.add_axes(ax2)
           #  # Calculate axis limits
           #  x_min = np.nanmin(x_data_plot)
           #  x_max = np.nanmax(x_data_plot)
           #  if cfg_par['gPlot']['Res-fixed_scale']:
           #      y1_min = np.nanmin([-200.,np.nanmin(-y_sigma_plot)*1.5,np.nanmin(-y_Res_plot)*1.5])
           #      y1_max = np.nanmax([+200.,np.nanmax(+y_sigma_plot)*1.5,np.nanmax(+y_Res_plot)*1.5])
           #  else:
           #      y1_min = np.nanmin(y_Res_plot)*1.1
           #      y1_max = np.nanmax(y_Res_plot)*1.1

           #  ax2.set_xticks(xTicks)

           #  #ax2.set_yticks([-150,0,150])
           #  #ax2.set_yticklabels([])

           #  # Set axis limits
           #  ax2.set_xlim(x_min, x_max)
           #  ax2.set_ylim(y1_min, y1_max) 



           #  #ax2.plot(vel, amp(x, p(x,m,n)))
           #  ax2.step(x_data_plot, y_Res_plot, 'g-', label='residuals')
           #  ax2.axhline(color='k', linestyle=':', zorder=0)                           
           #  ax2.axvline(color='k', linestyle=':', zorder=0)                           
           #  ax2.fill_between(x_data_plot, -y_sigma_plot, y_sigma_plot,
           #                   facecolor='grey', alpha=0.5,step='mid')

           # # for the last plot add the x-axis label
           #  if i >= len(lineInfo['Wave'])-3:                
           #      ax2.set_xticks(xTicks)
           #      ax2.set_xlabel(
           #              r'v\,[$\mathrm{km}\,\mathrm{s}^{-1}$]', labelpad=2)
           #  else:
           #      ax2.set_xticklabels([])


            k+=1
        #delete unused subplots
        #i+=1
        #while not i % 3 ==0:   
        #    fig.delaxes(ax)

        #    gs[j-1][k].axes.get_xaxis().set_ticklabels(xTicksStr)
        #    gs[j-1][k].axes.get_xaxis().set_ticklabels(xTicksStr)
        #    k +=1
        #    i +=1
        plotFormat= cfg_par['pvDiagram']['pvPlots']['plotFormat']
    
        baseName=cfg_par['general']['outPrefix'][0]+'pvPlotMultu.fits'
        outFigName= cfg_par['pvDiagram']['pvPlots']['plotDir']+'/'+baseName.replace('.fits','.'+plotFormat)
        print(outFigName)

        plt.savefig(outFigName,dpi=300,bbox_inches='tight',
                    format='png',overwrite=True) # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        #plt.show()


