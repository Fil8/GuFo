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
from matplotlib import patches as mpatches
from matplotlib import colors
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator
from matplotlib import transforms as mtransforms
from matplotlib.ticker import LogFormatter 
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import matplotlib.cm as cm

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
        'pdf.fonttype'        : 3,
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
        'text.latex.preamble' : r'\usepackage{amsmath}'
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

    colorTickLabels = np.linspace(vRange[0],vRange[1],9)    

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
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=300)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)
    return 0

  def momAncPlot(self,cfg_par,imageName,Labels,CustomCmap,colorList,
  contourColors='black',nameFigLabel=None,overlayContours=False,
  contName=None,contLevels=None,contColors=None):

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

    img = ax1.imshow(hduImCut.data, cmap=CustomCmap,interpolation='nearest')

    # get the colors of the values, according to the 
    # colormap used by imshow
    colours = [colors.to_rgba(colorList[i]) for i in range(0,len(Labels))]

    #colors = [img.cmap(img(value)) for value in values]
    #print(colors)
    # create a patch (proxy artist) for every color 
    patches = [ mpatches.Patch(color=colours[i],  label=Labels[i] ) for i in range(len(Labels)) ]
    # put those patched as legend-handles into the legend
    ax1.legend(handles=patches, loc=3, borderaxespad=0.)

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
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=300)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)
    return 0

  def momAncPlotOver(self,cfg_par,imageName,secondImageName,thirdImageName,Labels,CustomCmapOne,CustomCmapThree,colorList,
  contourColors='black',nameFigLabel=None,overlayContours=False,
  contName=None,contLevels=None,contColors=None):

    hduIm = fits.open(imageName)[0]
    hduImTwo = fits.open(secondImageName)[0]
    hduImThree = fits.open(thirdImageName)[0]
    hduImData= np.array(hduIm.data,dtype='float64')
    hduImTwoData= np.array(hduIm.data,dtype='float64')
    hduImTwoData= np.array(hduIm.data,dtype='float64')

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
    hduImTwoCut = Cutout2D(hduImTwo.data, centre, size, wcs=wcsIm)
    hduImThreeCut = Cutout2D(hduImThree.data, centre, size, wcs=wcsIm)

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

    #CustomCmapOne.set_under(color='white')
    #ccaMap = cm.Greens_r
    #ccaMap.set_over(color='white')
    #CustomCmapThree.set_under(color='white')

    img = ax1.imshow(hduImCut.data, cmap=CustomCmapOne,alpha=1,interpolation='nearest')
    
    img2 = ax1.imshow(hduImTwoCut.data, cmap='YlGn_r',alpha=1,interpolation='nearest',vmin=-1e-1,vmax=1.1)
    img3 = ax1.imshow(hduImThreeCut.data, cmap=CustomCmapThree,alpha=1,interpolation='nearest')

    # get the colors of the values, according to the 
    # colormap used by imshow
    colours = [colors.to_rgba(colorList[i]) for i in range(0,len(Labels))]

    #colors = [img.cmap(img(value)) for value in values]
    #print(colors)
    # create a patch (proxy artist) for every color 
    patches = [ mpatches.Patch(color=colours[i],  label=Labels[i] ) for i in range(len(Labels)) ]
    # put those patched as legend-handles into the legend
    ax1.legend(handles=patches, loc=3, borderaxespad=0.)

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

    outFig = cfg_par['general']['plotMomModDir']+outBPT+nameFigLabel+'-shade.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], #bbox_inches = "tight",overwrite=True,dpi=100)#,
            dpi=300,bbox_inches='tight',transparent=False,overwrite=True)
    return 0


  def momAncPlotMultiple(self,cfg_par):

    tableNames = cfg_par['multipleRegions']['tableNames']
    regionNames = cfg_par['multipleRegions']['regions']

    colorNames = cfg_par['multipleRegions']['colors']
    
    momModDir = cfg_par['general']['momDir']+'multipleRegions'+'/'

    colorMaps = ListedColormap(colorNames[0])
    print(colorMaps)

    imageName=[]
    for i in range(0,len(tableNames)):

        imageName.append(momModDir+'mom-'+regionNames[i]+'.fits')
    

    # hduImTwoCut = Cutout2D(hduImTwo.data, centre, size, wcs=wcsIm)
    # hduImThreeCut = Cutout2D(hduImThree.data, centre, size, wcs=wcsIm)

    hduIm = fits.open(imageName[0])[0]
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
    hduImCut = Cutout2D(hduIm.data, centre, size, wcs=wcsIm)

    params = self.loadRcParams()
    plt.rcParams.update(params)


    #idxLin = np.where(hduImCut==2.)
    #idxSey = np.where(hduImCut==1.)
    #idxKew = np.where(hduImCut==0.)
    #idxBad = np.where(hduImCut==-1.)        
    
    fig = plt.figure()
    
    ax1 = plt.subplot(projection=hduImCut.wcs)    

    divider = make_axes_locatable(ax1)


    #cax = divider.append_axes("right", size='2%', pad=0.1)
    
    #if vRange == None:
    #  vRange=np.array([1,2])
    #  vRange[0] = lineThresh
    #  vRange[1] = np.nanmax(hduImCut.data)

    #CustomCmapOne.set_under(color='white')
    #ccaMap = cm.Greens_r
    #ccaMap.set_over(color='white')
    #CustomCmapThree.set_under(color='white'

    
    for i in range(0,len(imageName)):
        colorMaps = ListedColormap(colorNames[i])

        hduIm = fits.open(imageName[i])[0]
    
        hduImCut = Cutout2D(hduIm.data, centre, size, wcs=wcsIm)

        img = ax1.imshow(hduImCut.data, cmap=colorMaps,alpha=1,interpolation='nearest',vmin=-1e-1,vmax=1.1)
    
    # get the colors of the values, according to the 
    # colormap used by imshow

    #colors = [img.cmap(img(value)) for value in values]
    #print(colors)
    # create a patch (proxy artist) for every color 
    patches = [ mpatches.Patch(color=colorNames[i],  label=regionNames[i] ) for i in range(len(regionNames)) ]
    # put those patched as legend-handles into the legend
    ax1.legend(handles=patches, loc=3, borderaxespad=0.,prop={'size': 18})

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

    
    #if overlayContours:
    #    imLevels =[-0.5,0.5]
        #contLevels = np.linspace(lineThresh*1.2,np.nanmax(hduImCut.data)*0.95,step)
    #    ax1.contour(hduImCut.data,levels=imLevels, colors=contourColors)
    #    nameFigLabel = nameFigLabel+'_cs'
    if cfg_par['multipleRegions']['contOver'] == True:
        contName=cfg_par['multipleRegions']['contName']
        contLevels=cfg_par['multipleRegions']['contLevels']
        hduCont = fits.open(contName)[0]
        wcsCont = WCS(hduCont.header)
        hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
                                            hduImCut.wcs, shape_out=hduImCut.shape)

        ax1.contour(array.data,levels=contLevels, colors='black')

    outFigDir = momModDir+'plots/'
    if not os.path.exists(outFigDir):
            os.mkdir(outFigDir)
    outFig = outFigDir +'KMultipleMom'+cfg_par['multipleRegions']['Name']+'.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], #bbox_inches = "tight",overwrite=True,dpi=100)#,
            dpi=300,bbox_inches='tight',transparent=False,overwrite=True)
    return 0




  def mom1Plot(self,cfg_par,imageName,lineName,lineThresh,lineNameStr,keyword,
    vRange=None,modName='g1',contourColors='black',nameFigLabel=None,overlayContours=False,
    contName=None,imLevels=None,contValues=None,contColors=None,interpolation=None):
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
      mom1BarLabel=r'[km s$^{-1}$]'
      cMap = 'jet'
    else:
      mom1BarLabel = r+str(cfg_par[keyword]['cBarLabel'][1])
      cMap = cfg_par[keyword]['colorMap'][1]

    img = ax1.imshow(hduImCut.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5,interpolation=interpolation)

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

    #ax1.set_title(lineNameStr)

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

        cs = ax1.contour(array.data,levels=contValues[i], colors=contColors[i])

        if contValues[i]==1:
            ax1.clabel(cs, inline=1, fontsize=14)

    outFig = cfg_par['general']['plotMomModDir']+outMom+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=300)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)


    return 0

  def mom2Plot(self,cfg_par,imageName,lineName,lineThresh,lineNameStr,keyword,vRange=None,
    modName='g1',contourColors='black',nameFigLabel=None,overlayContours=False,
    contName=None,contValues=None,contColors=None,interpolation=None):

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
  
    #vRange=np.array([1,2])
    #vRange[0] = np.nanmin(hduImCut.data)
    #vRange[1] = np.nanmax(hduImCut.data)

    if len(cfg_par[keyword]['cBarLabel'])>1:
      mom1BarLabel=r'[km s$^{-1}$]'
      cMap = plt.cm.jet
    else:
      mom1BarLabel = r+str(cfg_par[keyword]['cBarLabel'][2])
      cMap = plt.cm.cfg_par[keyword]['colorMap'][2]

    if cfg_par[keyword]['cBarLabel'][0].split(' ',maxsplit=2)[0]=='w80':
        mom1BarLabel=r'$w80$ [km s$^{-1}$]'
        cMap = plt.cm.Reds
        cMap.set_under(color='grey')

#    img = ax1.imshow(hduImCut.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5)
    img = ax1.imshow(hduImCut.data, cmap=cMap,vmin=vRange[0],vmax=vRange[1],
        interpolation=interpolation)

    colorTickLabels = np.arange(vRange[0],vRange[1]+100,100)    
    print(colorTickLabels)

    ax1.coords[1].set_axislabel(r'Dec (J2000)')
    ax1.coords[0].set_axislabel(r'RA (J2000)')
    ax1.coords[1].set_major_formatter('dd:mm:ss')
    ax1.coords[0].set_major_formatter('hh:mm:ss')

    
    cax.coords[0].grid(False)
    cax.coords[1].grid(False)
    cax.tick_params(direction='in')
    cax.coords[0].tick_params(top=False, bottom=False,
                   labeltop=False, labelbottom=False)
    cax.coords[1].set_ticklabel_position('r')

    cax.coords[1].set_axislabel(mom1BarLabel)
    cax.coords[1].set_axislabel_position('r')
    cax.coords[1].set_ticklabel_visible(True)     

    cbar = plt.colorbar(img, cax=cax,ticks =colorTickLabels,
                    orientation='vertical', format='%d')   
    
    if lineNameStr=='Hb4861':
      lineNameStr=r'H$_\beta$4861'
    elif lineNameStr=='Ha6562':
      lineNameStr=r'H$_\alpha$6562'

    #ax1.set_title(lineNameStr)

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

        cs = ax1.contour(array.data,levels=contValues[i], colors=contColors[i])
        if contValues[i]==1:
            ax1.clabel(cs, inline=1, fontsize=14)

    outFig = cfg_par['general']['plotMomModDir']+outMom+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=300,transparent=False)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)

    plt.close()
    return 0


  def mom12Triplet(self,cfg_par,im1,im2,im3,lineName,lineThresh,lineNameStr,
    vRange=None,modName='g1',contourColors='black',nameFigLabel=None,overlayContours=False,
    contName=None,imLevels=None,contValues=None,contColors=None,interpolation=None,title=False,kind=None):


    objCoordsRA = cfg_par['moments']['centreRA']
    objCoordsDec = cfg_par['moments']['centreDec']
    
    centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree, frame='fk5')
    size = u.Quantity((cfg_par['moments']['sizePlots'],cfg_par['moments']['sizePlots']), u.arcmin)
  
    params = self.loadRcParams()
    plt.rcParams.update(params)

    #mom0Map = fits.open(cfg_par['general']['momModDir']+'mom0_'+modName+'-'+lineName+'.fits')
    #hBetaData = mom0Map[0].data

    hduIm1 = fits.open(im1)[0]
    wcsIm1 = WCS(hduIm1.header)
    hduImCut1 = Cutout2D(hduIm1.data, centre, size, wcs=wcsIm1)

    hduIm2 = fits.open(im2)[0]
    wcsIm2 = WCS(hduIm2.header)
    hduImCut2 = Cutout2D(hduIm2.data, centre, size, wcs=wcsIm2)

    hduIm3 = fits.open(im3)[0]
    wcsIm3 = WCS(hduIm3.header)
    hduImCut3 = Cutout2D(hduIm3.data, centre, size, wcs=wcsIm3)
    #idx = np.where(np.isnan(hBetaData))
    #hduIm.data[idx] = np.nan

    fig = plt.figure(figsize =(12,10), constrained_layout=False)
    
    gs = plt.GridSpec(nrows=1, ncols=3,  figure=fig,wspace=0.0,hspace=0.0)

    ax1 = fig.add_subplot(gs[0,0],projection=wcsIm1)
    
    #divider = make_axes_locatable(ax1)
    #cax = divider.append_axes("right", size='2%', pad=0.1)
    if vRange==None:
        vRange=np.array([1,2])
        vRange[0] = np.nanmin(hduImCut.data)
        vRange[1] = np.nanmax(hduImCut.data)


    cMap = 'jet'


    c = SkyCoord('00:00:02.0','00:00:20.0',unit=(u.hourangle,u.deg))
    ax1.coords[0].set_ticks(spacing=c.ra.degree*u.degree)
    ax1.coords[1].set_ticks(spacing=c.dec.degree*u.degree)

    ax1.coords[1].set_axislabel(r'Dec (J2000)')

    ax1.coords[0].set_axislabel(r'RA (J2000)')

#    img = ax1.imshow(hduImCut.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5)
    img = ax1.imshow(hduImCut1.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5,
        interpolation=interpolation)

    if kind == 'mom2' or kind == None:
        colorTickLabels = np.arange(vRange[0],vRange[1]+100.,100.)    
    elif kind == 'mom1':
        colorTickLabels = np.linspace(vRange[0],vRange[1],5)  

    
    # cax.coords[0].grid(False)
    # cax.coords[1].grid(False)
    # cax.tick_params(direction='in')
    # cax.coords[0].tick_params(top=False, bottom=False,
    #                labeltop=False, labelbottom=False)
    # cax.coords[1].set_ticklabel_position('r')

    # cax.coords[1].set_axislabel(mom1BarLabel)
    # cax.coords[1].set_axislabel_position('r')
    # cbar = plt.colorbar(img, cax=cax,ticks =colorTickLabels,
    #                 orientation='vertical', format='%d')   
    
    if lineNameStr=='Hb4861':
      lineNameStr=r'H$_\beta$4861'
    elif lineNameStr=='Ha6562':
      lineNameStr=r'H$_\alpha$6562'

    if title==True:
        ax1.set_title('g1')

    ax1.set_autoscale_on(False)    
    #SaveOutput

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
                                            hduImCut1.wcs, shape_out=hduImCut1.shape)

        cs = ax1.contour(array.data,levels=contValues[i], colors=contColors[i])
        if contValues[i]==1:
            ax1.clabel(cs, inline=1, fontsize=14)

    ax2 = fig.add_subplot(gs[0,1],projection=wcsIm2)
    
    #divider = make_axes_locatable(ax1)
    #cax = divider.append_axes("right", size='2%', pad=0.1)


    cMap = 'jet'

#    img = ax1.imshow(hduImCut.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5)
    img = ax2.imshow(hduImCut2.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5,
        interpolation=interpolation)


    ax2.coords[0].set_ticks(spacing=c.ra.degree*u.degree)
    ax2.coords[1].set_ticks(spacing=c.dec.degree*u.degree)
    ax2.coords[1].set_axislabel(r'Dec (J2000)')
    ax2.coords[1].set_ticklabel_visible(False)
    ax2.coords[0].set_axislabel(r'RA (J2000)')
    
    # cax.coords[0].grid(False)
    # cax.coords[1].grid(False)
    # cax.tick_params(direction='in')
    # cax.coords[0].tick_params(top=False, bottom=False,
    #                labeltop=False, labelbottom=False)
    # cax.coords[1].set_ticklabel_position('r')

    # cax.coords[1].set_axislabel(mom1BarLabel)
    # cax.coords[1].set_axislabel_position('r')
    # cbar = plt.colorbar(img, cax=cax,ticks =colorTickLabels,
    #                 orientation='vertical', format='%d')   
    
    if lineNameStr=='Hb4861':
      lineNameStr=r'H$_\beta$4861'
    elif lineNameStr=='Ha6562':
      lineNameStr=r'H$_\alpha$6562'


    if title==True:
        ax2.set_title('g2')

    ax2.set_autoscale_on(False)    
    #SaveOutput
 
    modName = cfg_par['gFit']['modName'] 

    if nameFigLabel==None:
        nameFigLabel='' 
    if overlayContours:
        imLevels =lineThresh*1.2*(np.arange(1,10,2))
        #contLevels = np.linspace(lineThresh*1.2,np.nanmax(hduImCut.data)*0.95,step)
        cs = ax2.contour(hduImCut.data,levels=imLevels, colors=contourColors)
        nameFigLabel = nameFigLabel+'_cs'
        if contValues[0]==1:
            ax2.clabel(cs, inline=1, fontsize=14)

    if contName:
      if nameFigLabel=='':
        nameFigLabel='over_'
      for i in range(0,len(contName)):

        hduCont = fits.open(contName[i])[0]
        wcsCont = WCS(hduCont.header)
        hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
                                            hduImCut2.wcs, shape_out=hduImCut2.shape)

        cs = ax2.contour(array.data,levels=contValues[i], colors=contColors[i])
        if contValues[i]==1:
            ax2.clabel(cs, inline=1, fontsize=14)

    ax3 = fig.add_subplot(gs[0,2],projection=wcsIm3)
    
    #divider = make_axes_locatable(ax3)
    #cax = divider.append_axes("right", size='2%', pad=0.1)
  
    if kind=='mom2':
        mom1BarLabel=r'$\sigma_{\rm los}$ [km s$^{-1}$]'
        momName='mom2'

    elif kind=='mom1':
        mom1BarLabel=r'$v_{\rm los}-v_{\rm sys}$ [km s$^{-1}$]'
        momName='mom1'
    else:
        mom1BarLabel=r'[km s$^{-1}$]'
        momName='mom'

    cMap = 'jet'

#    img = ax1.imshow(hduImCut.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5)
    img = ax3.imshow(hduImCut3.data, cmap=cMap,vmin=vRange[0]-0.5,vmax=vRange[1]+0.5,
        interpolation=interpolation)

    ax3.coords[0].set_ticks(spacing=c.ra.degree*u.degree)
    ax3.coords[1].set_ticks(spacing=c.dec.degree*u.degree)
    ax3.coords[1].set_axislabel(r'Dec (J2000)')
    ax3.coords[1].set_ticklabel_visible(False)
    ax3.coords[0].set_axislabel(r'RA (J2000)')


    axins = inset_axes(ax3,
                   width="5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.05, 0., 1, 1),
                   bbox_transform=ax3.transAxes,
                   borderpad=0,
                   )

    #cax.coords[0].grid(False)
    #cax.coords[1].grid(False)
    #cax.tick_params(direction='in')
    #cax.coords[0].tick_params(top=False, bottom=False,
    #               labeltop=False, labelbottom=False)
    #cax.coords[1].set_ticklabel_position('r')

    #cax.coords[1].set_axislabel(mom1BarLabel)
    #cax.coords[1].set_axislabel_position('r')
    cbar = fig.colorbar(img, cax=axins,ticks =colorTickLabels,
                  orientation='vertical', format='%d')   
    cbar.set_label(mom1BarLabel, rotation=-90, va="bottom")
    if lineNameStr=='Hb4861':
      lineNameStr=r'H$_\beta$4861'
    elif lineNameStr=='Ha6562':
      lineNameStr=r'H$_\alpha$6562'

    if title==True:
        ax3.set_title('total')

    ax3.set_autoscale_on(False)    


    #outMom= str.split(outMom, '.fits')[0]+  
    modName = cfg_par['gFit']['modName'] 

    if nameFigLabel==None:
        nameFigLabel='' 
    if overlayContours:
        imLevels =lineThresh*1.2*(np.arange(1,10,2))
        #contLevels = np.linspace(lineThresh*1.2,np.nanmax(hduImCut.data)*0.95,step)
        cs = ax3.contour(hduImCut.data,levels=imLevels, colors=contourColors)
        nameFigLabel = nameFigLabel+'_cs'
        if contValues[0]==1:
            ax3.clabel(cs, inline=1, fontsize=14)

    if contName:
      if nameFigLabel=='':
        nameFigLabel='over_'
      for i in range(0,len(contName)):

        hduCont = fits.open(contName[i])[0]
        wcsCont = WCS(hduCont.header)
        hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
                                            hduImCut3.wcs, shape_out=hduImCut3.shape)

        cs = ax3.contour(array.data,levels=contValues[i], colors=contColors[i])
        if contValues[i]==1:
            ax3.clabel(cs, inline=1, fontsize=14)


    #SaveOutput
    if title==True:
        outMom = momName+'g1g2TotTitle'
    else:
        outMom = momName+'g1g2Tot'


    fig.subplots_adjust(hspace=0.,wspace=0.)


    outFig = cfg_par['general']['plotMomModDir']+outMom+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
    fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], dpi=300,bbox_inches = "tight",overwrite=True)#,
            #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)

    plt.close()
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


