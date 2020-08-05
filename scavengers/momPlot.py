#!/usr/bin/env python3.6

import os, sys
from astropy.io import fits
import numpy as np

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult

from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator
from matplotlib import transforms as mtransforms
from matplotlib.ticker import LogFormatter 
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from reproject import reproject_interp, reproject_exact


#import aplpy
import cvPlay
cvP = cvPlay.convert()

class MOMplot(object):

#----------------------#
# rc param initialize
#----------------------#
  def loadRcParams(self):
  
      params = {'figure.figsize'      : '10,10',
        'figure.autolayout' : True,
        'font.family'         :'serif',
        'font.serif'          :'times',
        'font.style'          : 'normal',
        'font.weight'         : 'book',
        'font.size'           : 24,
        'axes.linewidth'      : 2.2,
        'lines.linewidth'     : 2,
        'xtick.labelsize'     : 22,
        'ytick.labelsize'     : 22, 
        'xtick.direction'     :'in',
        'ytick.direction'     :'in',
        'xtick.major.size'    : 6,
        'xtick.major.width'   : 2,
        'xtick.minor.size'    : 3,
        'xtick.minor.width'   : 1,
        'ytick.major.size'    : 6,
        'ytick.major.width'   : 2,
        'ytick.minor.size'    : 3,
        'ytick.minor.width'   : 1, 
        'text.usetex'         : True,
        #'text.latex.unicode'  : True
         }
      
      return params

  #def mom0Plot(self,cfg_par,imageName):
  def mom0Plot(self,cfg_par,imageName,lineName,lineNameStr,lineThresh,
    vRange=None,contourColors='black',nameFigLabel=None,overlayContours=False,
    contName=None,contLevels=None,contColors=None):
      
    objCoordsRA = cfg_par['moments']['centreRA']
    objCoordsDec = cfg_par['moments']['centreDec']
    
    centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree, frame='fk5')
    size = u.Quantity((cfg_par['moments']['sizePlots'],cfg_par['moments']['sizePlots']), u.arcmin)
  
    params = self.loadRcParams()
    plt.rcParams.update(params)


    hduIm = fits.open(imageName)[0]
    wcsIm = WCS(hduIm.header)

    hduImCut = Cutout2D(hduIm.data, centre, size, wcs=wcsIm)
    
    fig = plt.figure()
    
    ax1 = plt.subplot(projection=wcsIm)    

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size='2%', pad=0.1)
    if vRange == None:
      vRange=np.array([1,2],dtype=float)
      vRange[0] = lineThresh
      vRange[1] = np.nanmax(hduImCut.data)/5.

    img = ax1.imshow(hduImCut.data, cmap=cfg_par['moments']['colorMap'][0],vmin=vRange[0],vmax=vRange[1])

    colorTickLabels = np.linspace(vRange[0],vRange[1],9.)    

    ax1.coords[1].set_axislabel(r'Dec (J2000)')
    ax1.coords[0].set_axislabel(r'RA (J2000)')
    
    cax.coords[0].grid(False)
    cax.coords[1].grid(False)
    cax.tick_params(direction='in')
    cax.coords[0].tick_params(top=False, bottom=False,
                   labeltop=False, labelbottom=False)
    cax.coords[1].set_ticklabel_position('r')
    cax.coords[1].set_axislabel(cfg_par['moments']['cBarLabel'][0])
    cax.coords[1].set_axislabel_position('r')
    cbar = plt.colorbar(img, cax=cax,ticks =colorTickLabels,
                    orientation='vertical', format='%d')   

    if lineNameStr=='Hb4861':
      lineNameStr=r'H$_\beta$4861'
    elif lineNameStr=='Ha6562':
      lineNameStr=r'H$_\alpha$6562'
    
    ax1.set_title(lineNameStr)

    ax1.set_autoscale_on(False)    
    #SaveOutput
    outMom = os.path.basename(imageName)
    outMom= str.split(outMom, '.fits')[0]  
    modName = cfg_par['gFit']['modName'] 

    if nameFigLabel==None:
        nameFigLabel='' 
    
    if overlayContours:
        imLevels =lineThresh*1.2*(np.arange(1,10,2))
        #contLevels = np.linspace(lineThresh*1.2,np.nanmax(hduImCut.data)*0.95,step)
        ax1.contour(hduImCut.data,levels=imLevels, colors=contourColors)
        nameFigLabel = nameFigLabel+'_cs'
    if contName:
      if nameFigLabel=='':
        nameFigLabel='over_'
      for i in range(0,len(contName)):

        hduCont = fits.open(contName[i])[0]
        wcsCont = WCS(hduCont.header)
        hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
                                            hduImCut.wcs, shape_out=hduImCut.shape)

        ax1.contour(array.data,levels=contLevels[i], colors=contColors[i])
      

    outFig = cfg_par['general']['plotMomModDir']+outMom+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)
    return 0

  def momAncPlot(self,cfg_par,imageName,lineName,lineNameStr,lineThresh,
  contourColors='black',nameFigLabel=None,overlayContours=False,
  contName=None,contLevels=None,contColors=None):


    if cfg_par['ancillary']['plotRotation'] == True:
      CustomCmap = ListedColormap(['blue','darkseagreen'])
      cBarTickLabels= ['CCA','NII6583']

    else:
      CustomCmap = ListedColormap(['blue','darkseagreen','crimson'])
      cBarTickLabels= ['NII6583','CCA','Rotation']


    hduIm = fits.open(imageName)[0]
    wcsIm = WCS(hduIm.header)

    #sn = resTable['SN_OIII5006']
    #idx  = np.where(sn<=lineThresh)

    #hduIm.data[idx] = np.nan 

    #sigmaThresh = sigmaTable['g1_SigIntr_OIII5006']
    #idx  = np.where(sigmaThresh>=cfg_par['moments']['sigmaThresh'])

    #hduIm.data[idx] = np.nan 
    
    objCoordsRA = cfg_par['moments']['centreRA']
    objCoordsDec = cfg_par['moments']['centreDec']
    
    centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree, frame='fk5')
    size = u.Quantity((cfg_par['moments']['sizePlots'],cfg_par['moments']['sizePlots']), u.arcmin)
  
    params = self.loadRcParams()
    plt.rcParams.update(params)

    hduImCut = Cutout2D(hduIm.data, centre, size, wcs=wcsIm)
    
    #idxLin = np.where(hduImCut==2.)
    #idxSey = np.where(hduImCut==1.)
    #idxKew = np.where(hduImCut==0.)
    #idxBad = np.where(hduImCut==-1.)        
    
    fig = plt.figure()
    
    ax1 = plt.subplot(projection=wcsIm)    

    divider = make_axes_locatable(ax1)
    #cax = divider.append_axes("right", size='2%', pad=0.1)
    
    #if vRange == None:
    #  vRange=np.array([1,2])
    #  vRange[0] = lineThresh
    #  vRange[1] = np.nanmax(hduImCut.data)

    img = ax1.imshow(hduImCut.data, cmap=CustomCmap)

    #cBarTicks = [-1,0,1,2]

    ax1.coords[1].set_axislabel(r'Dec (J2000)')
    ax1.coords[0].set_axislabel(r'RA (J2000)')
    
    #cax.coords[0].grid(False)
    #cax.coords[1].grid(False)
    #cax.tick_params(direction='in')
    #cax.coords[0].tick_params(top=False, bottom=False,
    #              labeltop=False, labelbottom=False)
    #cax.coords[1].set_ticklabel_position('r')
    #cax.coords[1].tick_params(top=False, bottom=False,
    #               labeltop=False, labelbottom=False)
    #cax.coords[1].set_axislabel('$\eta$-parameter')
    #cax.coords[1].set_axislabel_position('r')
    #cbar = plt.colorbar(img, cax=cax,ticks =cBarTicks,
    #               orientation='vertical', format='%d')   
    
    ax1.set_autoscale_on(False)    
    #SaveOutput
    outBPT = os.path.basename(imageName)
    outBPT= str.split(outBPT, '.fits')[0]  
    modName = cfg_par['gFit']['modName'] 

    if nameFigLabel==None:
        nameFigLabel='' 
    
    #if overlayContours:
    #    imLevels =[-0.5,0.5]
        #contLevels = np.linspace(lineThresh*1.2,np.nanmax(hduImCut.data)*0.95,step)
    #    ax1.contour(hduImCut.data,levels=imLevels, colors=contourColors)
    #    nameFigLabel = nameFigLabel+'_cs'
    if contName:
      if nameFigLabel=='':
        nameFigLabel='over_'
      for i in range(0,len(contName)):

        hduCont = fits.open(contName[i])[0]
        wcsCont = WCS(hduCont.header)
        hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
                                            hduImCut.wcs, shape_out=hduImCut.shape)

        ax1.contour(array.data,levels=contLevels[i], colors=contColors[0])

    outFig = cfg_par['general']['plotMomModDir']+outBPT+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)
    return 0

  def mom1Plot(self,cfg_par,imageName,lineName,lineThresh,lineNameStr,keyword,
    vRange=None,modName='g1',contourColors='black',nameFigLabel=None,overlayContours=False,
    contName=None,contLevels=None,contValues=None,contColors=None):

    objCoordsRA = cfg_par['moments']['centreRA']
    objCoordsDec = cfg_par['moments']['centreDec']
    
    centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree, frame='fk5')
    size = u.Quantity((cfg_par['moments']['sizePlots'],cfg_par['moments']['sizePlots']), u.arcmin)
  
    params = self.loadRcParams()
    plt.rcParams.update(params)

    #mom0Map = fits.open(cfg_par['general']['momModDir']+'mom0_'+modName+'-'+lineName+'.fits')
    #hBetaData = mom0Map[0].data


    hduIm = fits.open(imageName)[0]
    wcsIm = WCS(hduIm.header)

    #idx = np.where(np.isnan(hBetaData))
    #hduIm.data[idx] = np.nan


    hduImCut = Cutout2D(hduIm.data, centre, size, wcs=wcsIm)

    
    fig = plt.figure()
    
    ax1 = plt.subplot(projection=wcsIm)    

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size='2%', pad=0.1)
    
    if vRange == None:
      vRange=np.array([1,2])
      vRange[0] = np.nanmin(hduImCut.data)
      vRange[1] = np.nanmax(hduImCut.data)

    if len(cfg_par[keyword]['cBarLabel'])>1:
      mom1BarLabel=r'velocity [km s$^{-1}$]'
      cMap = 'jet'
    else:
      mom1BarLabel = r+str(cfg_par[keyword]['cBarLabel'][1])
      cMap = cfg_par[keyword]['colorMap'][1]

    img = ax1.imshow(hduImCut.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5)

    colorTickLabels = np.linspace(vRange[0],vRange[1],9)    

    ax1.coords[1].set_axislabel(r'Dec (J2000)')
    ax1.coords[0].set_axislabel(r'RA (J2000)')
    
    cax.coords[0].grid(False)
    cax.coords[1].grid(False)
    cax.tick_params(direction='in')
    cax.coords[0].tick_params(top=False, bottom=False,
                   labeltop=False, labelbottom=False)
    cax.coords[1].set_ticklabel_position('r')

    cax.coords[1].set_axislabel(mom1BarLabel)
    cax.coords[1].set_axislabel_position('r')
    cbar = plt.colorbar(img, cax=cax,ticks =colorTickLabels,
                    orientation='vertical', format='%d')   
    
    if lineNameStr=='Hb4861':
      lineNameStr=r'H$_\beta$4861'
    elif lineNameStr=='Ha6562':
      lineNameStr=r'H$_\alpha$6562'

    ax1.set_title(lineNameStr)

    ax1.set_autoscale_on(False)    
    #SaveOutput
    outMom = os.path.basename(imageName)
    outMom= str.split(outMom, '.fits')[0]  
    modName = cfg_par['gFit']['modName'] 



    if nameFigLabel==None:
        nameFigLabel='' 
    if overlayContours:
        imLevels =lineThresh*1.2*(np.arange(1,10,2))
        #contLevels = np.linspace(lineThresh*1.2,np.nanmax(hduImCut.data)*0.95,step)
        cs = ax1.contour(hduImCut.data,levels=imLevels, colors=contourColors)
        nameFigLabel = nameFigLabel+'_cs'
        if contValues[0]==1:
            ax1.clabel(cs, inline=1, fontsize=14)

    if contName:
      if nameFigLabel=='':
        nameFigLabel='over_'
      for i in range(0,len(contName)):

        hduCont = fits.open(contName[i])[0]
        wcsCont = WCS(hduCont.header)
        hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
                                            hduImCut.wcs, shape_out=hduImCut.shape)

        cs = ax1.contour(array.data,levels=contLevels[i], colors=contColors[i])
        if contValues[i]==1:
            ax1.clabel(cs, inline=1, fontsize=14)

    outFig = cfg_par['general']['plotMomModDir']+outMom+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)


    return 0

  def mom2Plot(self,cfg_par,imageName,lineName,lineThresh,lineNameStr,keyword,
    modName='g1',contourColors='black',nameFigLabel=None,overlayContours=False,
    contName=None,contLevels=None,contValues=None,contColors=None):

    objCoordsRA = cfg_par['moments']['centreRA']
    objCoordsDec = cfg_par['moments']['centreDec']
    
    centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree, frame='fk5')
    size = u.Quantity((cfg_par['moments']['sizePlots'],cfg_par['moments']['sizePlots']), u.arcmin)
  
    params = self.loadRcParams()
    plt.rcParams.update(params)

    #mom0Map = fits.open(cfg_par['general']['momModDir']+'mom0_'+modName+'-'+lineName+'.fits')
    #hBetaData = mom0Map[0].data

    hduIm = fits.open(imageName)[0]
    wcsIm = WCS(hduIm.header)

    #idx = np.where(np.isnan(hBetaData))
    #hduIm.data[idx] = np.nan


    hduImCut = Cutout2D(hduIm.data, centre, size, wcs=wcsIm)
    
    fig = plt.figure()
    
    ax1 = plt.subplot(projection=wcsIm)    

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size='2%', pad=0.1)
  
    vRange=np.array([1,2])
    vRange[0] = np.nanmin(hduImCut.data)
    vRange[1] = np.nanmax(hduImCut.data)

    if len(cfg_par[keyword]['cBarLabel'])>1:
      mom1BarLabel=r' [km s$^{-1}$]'
      cMap = 'jet'
    else:
      mom1BarLabel = r+str(cfg_par[keyword]['cBarLabel'][2])
      cMap = cfg_par[keyword]['colorMap'][2]

    img = ax1.imshow(hduImCut.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5)

    colorTickLabels = np.linspace(vRange[0],vRange[1],9)    

    ax1.coords[1].set_axislabel(r'Dec (J2000)')
    ax1.coords[0].set_axislabel(r'RA (J2000)')
    
    cax.coords[0].grid(False)
    cax.coords[1].grid(False)
    cax.tick_params(direction='in')
    cax.coords[0].tick_params(top=False, bottom=False,
                   labeltop=False, labelbottom=False)
    cax.coords[1].set_ticklabel_position('r')

    cax.coords[1].set_axislabel(mom1BarLabel)
    cax.coords[1].set_axislabel_position('r')
    cbar = plt.colorbar(img, cax=cax,ticks =colorTickLabels,
                    orientation='vertical', format='%d')   
    
    if lineNameStr=='Hb4861':
      lineNameStr=r'H$_\beta$4861'
    elif lineNameStr=='Ha6562':
      lineNameStr=r'H$_\alpha$6562'

    ax1.set_title(lineNameStr)

    ax1.set_autoscale_on(False)    
    #SaveOutput
    outMom = os.path.basename(imageName)
    outMom= str.split(outMom, '.fits')[0]  
    modName = cfg_par['gFit']['modName'] 



    if nameFigLabel==None:
        nameFigLabel='' 
    if overlayContours:
        imLevels =lineThresh*1.2*(np.arange(1,10,2))
        #contLevels = np.linspace(lineThresh*1.2,np.nanmax(hduImCut.data)*0.95,step)
        cs = ax1.contour(hduImCut.data,levels=imLevels, colors=contourColors)
        nameFigLabel = nameFigLabel+'_cs'
        if contValues[0]==1:
            ax1.clabel(cs, inline=1, fontsize=14)

    if contName:
      if nameFigLabel=='':
        nameFigLabel='over_'
      for i in range(0,len(contName)):

        hduCont = fits.open(contName[i])[0]
        wcsCont = WCS(hduCont.header)
        hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
                                            hduImCut.wcs, shape_out=hduImCut.shape)

        cs = ax1.contour(array.data,levels=contLevels[i], colors=contColors[i])
        if contValues[i]==1:
            ax1.clabel(cs, inline=1, fontsize=14)

    outFig = cfg_par['general']['plotMomModDir']+outMom+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)


    return 0


      # params = self.loadRcParams()
      # plt.rcParams.update(params)
      # fig = plt.figure()

      # f1 = aplpy.FITSFigure(imageName, figure=fig)
      # f1.set_theme('publication')

      # f1.frame.set_linewidth(2)

      # f1.show_colorscale(aspect='equal', cmap='nipy_spectral',stretch = 'linear',vmin= -cenRange, vmax= cenRange)

      # f=fits.open(imageName)
      # dd=f[0].data
      
      # #print(outBPT)
      # #print(idxLin,idxSey,idxKew,idxBad)

      # #f1.show_contour(imageName,levels=[1, 5, 8, 11, 15], colors='black')
      # #f1.show_contour(imageName,levels=[1e-3, 1e-2, 1e-1,5e-1], colors='black')

      # f1.axis_labels.set_font( weight='book',size='medium', 
      #                          stretch='normal', family='serif', 
      #                          style='normal', variant='normal')
      # f1.axis_labels.set_xtext('RA (J2000)')
      # f1.axis_labels.set_ytext('Dec (J2000)')
      # f1.tick_labels.set_xformat('hh:mm:ss')
      # f1.tick_labels.set_yformat('dd:mm:ss')
      # f1.tick_labels.set_font( weight='book', size='small',
      #                          stretch='normal', family='serif', 
      #                          style='normal', variant='normal') 
      # f1.ticks.set_color('k')
      # f1.ticks.set_length(6)  # points
      # f1.ticks.set_linewidth(2)  # points
      # f1.ticks.set_minor_frequency(3)
      # f1.ticks.show()

      # outMom = os.path.basename(imageName)
      # outMom= str.split(outMom, '.fits')[0]  
      # modName = cfg_par['gFit']['modName']
      
      # outMom = cfg_par['general']['momDir']+modName+'/plots/'+outMom+'.'+cfg_par['moments']['plotFormat']
      
      # if os.path.exists(cfg_par['general']['momDir']+modName+'/plots/') == False:
      #     os.mkdir(cfg_par['general']['momDir']+modName+'/plots/')

      # fig.savefig(outMom,format=cfg_par['moments']['plotFormat'])
      #         #if pdf, dpi=300,bbox_inches='tight',transparent=False,overwrite=True)

  def resPlot(self,cfg_par,imageName,lineName,lineThresh,
    vRange=None,contourColors='black',nameFigLabel=None,overlayContours=False,
    contName=None,contLevels=None,contColors=None):

    objCoordsRA = cfg_par['moments']['centreRA']
    objCoordsDec = cfg_par['moments']['centreDec']
    
    centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree, frame='fk5')
    size = u.Quantity((cfg_par['moments']['sizePlots'],cfg_par['moments']['sizePlots']), u.arcmin)
  
    params = self.loadRcParams()
    plt.rcParams.update(params)


    hduIm = fits.open(imageName)[0]
    wcsIm = WCS(hduIm.header)

    hduImCut = Cutout2D(hduIm.data, centre, size, wcs=wcsIm)

    if contName: 
        hduCont = fits.open(contName)[0]
        idx = np.isnan(hduCont.data)
        hduImCut.data[idx] = np.nan
    
    fig = plt.figure()
    
    ax1 = plt.subplot(projection=wcsIm)    

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size='2%', pad=0.1)
    
    if vRange == None:
      vRange=np.array([1,2])
      vRange[0] = lineThresh
      vRange[1] = np.nanmax(hduImCut.data)
    
    #normalizeData
    normalized = (hduImCut.data-np.nanmin(hduImCut.data))/(np.nanmax(hduImCut.data)-np.nanmin(hduImCut.data))
    img = ax1.imshow(normalized, cmap=cfg_par['moments']['colorMap'][0])#,,vmin=-0.001,vmax=1.001

    colorTickLabels = np.linspace(vRange[0],vRange[1],9.)    

    ax1.coords[1].set_axislabel(r'Dec (J2000)')
    ax1.coords[0].set_axislabel(r'RA (J2000)')
    
    cax.coords[0].grid(False)
    cax.coords[1].grid(False)
    cax.tick_params(direction='in')
    cax.coords[0].tick_params(top=False, bottom=False,
                   labeltop=False, labelbottom=False)
    cax.coords[1].set_ticklabel_position('r')
    cax.coords[1].set_axislabel(cfg_par['moments']['cBarLabel'][0])
    cax.coords[1].set_axislabel_position('r')
    cbar = plt.colorbar(img, cax=cax,ticks =colorTickLabels,
                    orientation='vertical', format='%d')   
    
    ax1.set_title(lineName)

    ax1.set_autoscale_on(False)    
    #SaveOutput
    outMom = os.path.basename(imageName)
    outMom= str.split(outMom, '.fits')[0]  
    modName = cfg_par['gFit']['modName'] 

    if nameFigLabel==None:
        nameFigLabel='' 
    
    if overlayContours:
        #contLevels = np.linspace(lineThresh*1.2,np.nanmax(hduImCut.data)*0.95,step)
        nameFigLabel = nameFigLabel+'_cs'
        hduCont = fits.open(contName)[0]
        wcsCont = WCS(hduCont.header)
        hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
                                            hduImCut.wcs, shape_out=hduImCut.shape)
        if contLevels==None:
          contLevels = lineThresh*1.2*(np.arange(1,10,2))

        ax1.contour(array.data,levels=contLevels, colors=contColors)
      

    outFig = cfg_par['general']['momDir']+outMom+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)



    return 0


