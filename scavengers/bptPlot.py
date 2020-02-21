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
import matplotlib.axes as maxes

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from reproject import reproject_interp, reproject_exact

import cvPlay, tPlay

cvP = cvPlay.convert()
tP = tPlay.tplay()
class BPTplot(object):

#----------------------#
# rc param initialize
#----------------------#
    def loadRcParams(self):
    
        params = {'figure.figsize'      : '10,10',
          'font.family'         :' serif',
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
          'text.latex.unicode'  : True
           }
        
        return params

    def bptOIII(self,cfg_par,outPlotDir=None,xyRange=None):

        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data
        lineRatios = hdul['LineRatios_'+cfg_par['gFit']['modName']].data

        ampTable = hdul['LineRes_'+cfg_par['gFit']['modName']].data
        amps = ampTable[cfg_par['gFit']['modName']+'-AmpSpax'+'_Hb4861']

        lineInfo = tP.openLineList(cfg_par)
        lineThresh = float(lineInfo['SNThresh'][0])

        idx  = np.where(amps>=lineThresh)


        if cfg_par['gFit']['modName'] == 'g1':
            modString = ['G1']
        elif cfg_par['gFit']['modName'] == 'g2':
            modString = ['G1','G2','ToT']
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']


        for i in range (0, len(modString)):
            
            # initialize figure
            params = self.loadRcParams()
            plt.rcParams.update(params)
            fig = plt.figure(figsize =(10,8))
            fig.subplots_adjust(hspace=0.0)
            gs = gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0])

            y = lineRatios['log_'+modString[i]+'-OIII5006/Hb4861'][idx]
            x = lineRatios['log_'+modString[i]+'-NII6583/Ha6562'][idx]
            k = lineBPT[modString[i]+'-BPT_OIII'][idx]
            
            #ax.set_xticks([])
            
            ax1.set_xlabel(r'log([NII] 6583/H$_\alpha$ 6562)')
            ax1.set_ylabel(r'log([OIII]5006/H$_\beta$4861)')


            # Calculate axis limits and aspect ratio
            if xyRange:
                xMin = xyRange[0]
                xMax = xyRange[1]
                yMin = xyRange[2]
                yMax = xyRange[3]
            else:
                # Calculate axis limits and aspect ratio
                xMin = -2.
                xMax = 0.6
                yMin = -2.
                yMax = 1.5

            # Set axis limits
            ax1.set_xlim(xMin, xMax)
            ax1.set_ylim(yMin, yMax)

            ax1.tick_params(axis='both', which='major', pad=5)
            #ax1.xaxis.set_minor_locator()
            #ax1.yaxis.set_minor_locator()      
            
            idxAGN = np.where(k==2.)
            idxKauf = np.where(k==0.)
            idxKew = np.where(k==1.)
            idxBad = np.where(k==-1.)

            #print((idxAGN),(idxKew),(idxKauf),(idxBad))

            ax1.scatter(x[idxAGN], y[idxAGN], c='red', marker='+', s=80, linewidths=4, label='AGN')
            ax1.scatter(x[idxKew], y[idxKew], c='cyan', marker='+', s=80, linewidths=4, label='SF-Kewley')
            ax1.scatter(x[idxKauf], y[idxKauf], c='blue', marker='+', s=80, linewidths=4, label ='SF-Kauffmann')
            ax1.scatter(x[idxBad], y[idxBad], c='grey', marker='+', s=80, linewidths=4)

            kaX = np.log10(np.linspace(np.power(10,-1.),np.power(10,0.),1e4))
            kaY = 0.61 / (kaX - 0.05) + 1.3
            ax1.plot(kaX, kaY, ls='--',c='black', label='Kewley et al. 2001')

            keX = np.log10(np.linspace(np.power(10,-2.),np.power(10,0.5),1e4))

            keY= 0.61 / (keX - 0.47) + 1.19
            ax1.plot(keX, keY, ls=':',c='black', label='Kauffmann et al. 2003')


            ax1.legend = plt.legend(loc=1, prop={'size': 12})
            ax1.legend.get_frame().set_edgecolor('black')
            ax1.legend.get_frame().set_facecolor('white')

            if outPlotDir==None:
                outPlotDir = cfg_par['general']['bptDir']+cfg_par['gFit']['modName']+'/plots/'

            if not os.path.exists(outPlotDir):
                os.mkdir(outPlotDir)

            outPlot = outPlotDir+'BPT-'+modString[i]+'-OIII.png'
           
            plt.savefig(outPlot,format=cfg_par['lineRatios']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
                    # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
            plt.show()
            plt.close()
               
        return 0

    def bptSII(self,cfg_par,outPlotDir=None,xyRange=None):


        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data
        lineRatios = hdul['LineRatios_'+cfg_par['gFit']['modName']].data


        ampTable = hdul['LineRes_'+cfg_par['gFit']['modName']].data
        amps = ampTable[cfg_par['gFit']['modName']+'-AmpSpax'+'_Hb4861']

        lineInfo = tP.openLineList(cfg_par)
        lineThresh = float(lineInfo['SNThresh'][0])

        idx  = np.where(amps>=lineThresh)

        if cfg_par['gFit']['modName'] == 'g1':
            modString = ['G1']
        elif cfg_par['gFit']['modName'] == 'g2':
            modString = ['G1','G2','ToT']
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']


        for i in range (0, len(modString)):

            # initialize figure
            params = self.loadRcParams()
            plt.rcParams.update(params)
            fig = plt.figure(figsize =(10,8))
            fig.subplots_adjust(hspace=0.0)
            gs = gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0])

            y = lineRatios['log_'+modString[i]+'-OIII5006/Hb4861'][idx]
            x = lineRatios['log_'+modString[i]+'-SII6716/Ha6562'][idx]
            k = lineBPT[modString[i]+'-BPT_SII'][idx]
            #ax.set_xticks([])

            ax1.set_xlabel(r'log([SII] 6716,6730/H$_\alpha$ 6562)')
            ax1.set_ylabel(r'log([OIII]5006/H$_\beta$4861)')


            # Calculate axis limits and aspect ratio
            if xyRange:
                xMin = xyRange[0]
                xMax = xyRange[1]
                yMin = xyRange[2]
                yMax = xyRange[3]
            else:
                # Calculate axis limits and aspect ratio
                xMin = -2.
                xMax = 0.6
                yMin = -2.
                yMax = 1.5

            # Set axis limits
            ax1.set_xlim(xMin, xMax)
            ax1.set_ylim(yMin, yMax)

            ax1.tick_params(axis='both', which='major', pad=5)
            #ax1.xaxis.set_minor_locator()
            #ax1.yaxis.set_minor_locator()      
            
            idxLin = np.where(k==2.)
            idxSey = np.where(k==1.)
            idxKew = np.where(k==0.)
            idxBad = np.where(k==-1.)

            #print((idxLin),(idxSey),(idxKew),(idxBad))

            ax1.scatter(x[idxLin], y[idxLin], c='green', marker='+', s=80, linewidths=4, label='LINER')
            ax1.scatter(x[idxKew], y[idxKew], c='blue', marker='+', s=80, linewidths=4, label='SF')
            ax1.scatter(x[idxSey], y[idxSey], c='red', marker='+', s=80, linewidths=4, label='Seyfert')
            ax1.scatter(x[idxBad], y[idxBad], c='grey', marker='+', s=80, linewidths=4)

            keX = np.log10(np.linspace(np.power(10,-2.),np.power(10,0.5),1e4))
            keY= 0.72 / (keX - 0.32) + 1.30
            ax1.plot(keX, keY, ls=':',c='black', label='Kauffmann et al. 2003')
            
            seyX = np.log10(np.linspace(np.power(10,-0.4),np.power(10,0.5),1e4))
            seyLine = 1.89*seyX + 0.76
            ax1.plot(seyX, seyLine, ls=':',c='black', label='Seyfert-LINER')
        
            ax1.legend = plt.legend(loc=1, prop={'size': 12})
            ax1.legend.get_frame().set_edgecolor('black')
            ax1.legend.get_frame().set_facecolor('white')

            if outPlotDir==None:
                outPlotDir = cfg_par['general']['bptDir']+cfg_par['gFit']['modName']+'/plots/'

            if not os.path.exists(outPlotDir):
                os.mkdir(outPlotDir)

            outPlot = outPlotDir+'BPT-'+modString[i]+'-SII.png'
           
            plt.savefig(outPlot,format=cfg_par['lineRatios']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
                    # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
            plt.show()
            plt.close()
               
        return 0

    def bptOI(self,cfg_par,outPlotDir=None,xyRange=None):


        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data
        lineRatios = hdul['LineRatios_'+cfg_par['gFit']['modName']].data

        ampTable = hdul['LineRes_'+cfg_par['gFit']['modName']].data
        amps = ampTable[cfg_par['gFit']['modName']+'-AmpSpax'+'_Hb4861']

        lineInfo = tP.openLineList(cfg_par)
        lineThresh = float(lineInfo['SNThresh'][0])

        idx  = np.where(amps>=lineThresh)


        if cfg_par['gFit']['modName'] == 'g1':
            modString = ['G1']
        elif cfg_par['gFit']['modName'] == 'g2':
            modString = ['G1','G2','ToT']
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']


        for i in range (0, len(modString)):

            # initialize figure
            params = self.loadRcParams()
            plt.rcParams.update(params)
            fig = plt.figure(figsize =(10,8))
            fig.subplots_adjust(hspace=0.0)
            gs = gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0])            

            y = lineRatios['log_'+modString[i]+'-OIII5006/Hb4861'][idx]
            x = lineRatios['log_'+modString[i]+'-OI6300/Ha6562'][idx]
            k = lineBPT[modString[i]+'-BPT_OI'][idx]
            #ax.set_xticks([])

            ax1.set_xlabel(r'log([OI] 6300/H$_\alpha$ 6562)')
            ax1.set_ylabel(r'log([OIII]5006/H$_\beta$4861)')


            # Calculate axis limits and aspect ratio
            if xyRange:
                xMin = xyRange[0]
                xMax = xyRange[1]
                yMin = xyRange[2]
                yMax = xyRange[3]
            else:
                # Calculate axis limits and aspect ratio
                xMin = -2.
                xMax = 0.6
                yMin = -2.
                yMax = 1.5

            # Set axis limits
            ax1.set_xlim(xMin, xMax)
            ax1.set_ylim(yMin, yMax)

            ax1.tick_params(axis='both', which='major', pad=5)
            #ax1.xaxis.set_minor_locator()
            #ax1.yaxis.set_minor_locator()      
            
            idxLin = np.where(k==2.)
            idxSey = np.where(k==1.)
            idxKew = np.where(k==0.)
            idxBad = np.where(k==-1.)

            #print((idxLin),(idxSey),(idxKew),(idxBad))

            ax1.scatter(x[idxLin], y[idxLin], c='green', marker='+', s=80, linewidths=4, label='LINER')
            ax1.scatter(x[idxKew], y[idxKew], c='blue', marker='+', s=80, linewidths=4, label='SF')
            ax1.scatter(x[idxSey], y[idxSey], c='red', marker='+', s=80, linewidths=4, label= 'Seyfert')
            ax1.scatter(x[idxBad], y[idxBad], c='grey', marker='+', s=80, linewidths=4)


            kaX = np.log10(np.linspace(np.power(10,-3.),np.power(10,-0.595),1e4))
            kaY = 0.73 / (kaX + 0.59) + 1.33
            ax1.plot(kaX, kaY, ls='--',c='black', label='Kewley et al. 2001')

            seyX = np.log10(np.linspace(np.power(10,-1.),np.power(10,0.0),1e4))
            seyLine = 1.18*seyX + 1.30
            ax1.plot(seyX, seyLine, ls=':',c='black', label='Seyfert-LINER')
     
            ax1.legend = plt.legend(loc=1, prop={'size': 12})
            ax1.legend.get_frame().set_edgecolor('black')
            ax1.legend.get_frame().set_facecolor('white')
            
            if outPlotDir==None:
                outPlotDir = cfg_par['general']['bptDir']+cfg_par['gFit']['modName']+'/plots/'

            if not os.path.exists(outPlotDir):
                os.mkdir(outPlotDir)

            outPlot = outPlotDir+'BPT-'+modString[i]+'-OI.png'
            fig.savefig(outPlot,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
                #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)
            #plt.show()

        return 0

    def bptIM(self,cfg_par,outBPT,
    contourColors='black',nameFigLabel=None,overlayContours=False,
    contName=None,contLevels=None,contColors=None):

        if (os.path.basename(outBPT) == 'BPT-G1-BPT_OIII.fits' or os.path.basename(outBPT) == 'BPT-G2-BPT_OIII.fits' or os.path.basename(outBPT) == 'BPT-G3-BPT_OIII.fits' 
            or os.path.basename(outBPT) == 'BPT-ToT-BPT_OIII.fits' ):
            CustomCmap = ListedColormap(['grey','blue','cyan','red'])
            cBarTickLabels= ['bad fit','SF-Kewley', 'SF-Kauffmann', 'AGN']
        else:
            CustomCmap = ListedColormap(['grey','blue','red','green'])
            cBarTickLabels = ['bad fit','SF','AGN','LINER']


        objCoordsRA = cfg_par['moments']['centreRA']
        objCoordsDec = cfg_par['moments']['centreDec']
        
        centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree, frame='fk5')
        size = u.Quantity((cfg_par['moments']['sizePlots'],cfg_par['moments']['sizePlots']), u.arcmin)
      
        params = self.loadRcParams()
        plt.rcParams.update(params)


        hduIm = fits.open(outBPT)[0]
        wcsIm = WCS(hduIm.header)

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
        outBPT = os.path.basename(outBPT)
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

        outFig = cfg_par['general']['bptDir']+outBPT+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
        fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
                #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)
        return 0


        # params = self.loadRcParams()
        # plt.rcParams.update(params)
        # fig = plt.figure()




        # f1 = aplpy.FITSFigure(outBPT, figure=fig)
        # f1.set_theme('publication')

        # f1.frame.set_linewidth(2)

        # f1.show_colorscale(aspect='equal', cmap=CustomCmap,stretch = 'linear')

        # f=fits.open(outBPT)
        # dd=f[0].data
        # idxLin = np.where(dd==2.)
        # idxSey = np.where(dd==1.)
        # idxKew = np.where(dd==0.)
        # idxBad = np.where(dd==-1.)
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

        # outBPT = os.path.basename(outBPT)
        # outBPT= str.split(outBPT, '.fits')[0]
        # modName = cfg_par['gFit']['modName']
        # outBPT = cfg_par['general']['bptDir']+modName+'/plots/'+outBPT+'.'+cfg_par['lineRatios']['plotFormat']
        # if os.path.exists(cfg_par['general']['bptDir']+modName+'/plots/') == False:
        #     os.mkdir(cfg_par['general']['bptDir']+modName+'/plots/')

        # fig.savefig(outBPT,format=cfg_par['lineRatios']['plotFormat'])
        #         #if pdf, dpi=300,bbox_inches='tight',transparent=False,overwrite=True)

        # return 0

    def cDistIM(self,cfg_par,outBPT,
    vRange=None,contourColors='black',nameFigLabel=None,overlayContours=False,
    contName=None,contLevels=None,contColors=None):
      
        objCoordsRA = cfg_par['moments']['centreRA']
        objCoordsDec = cfg_par['moments']['centreDec']
        
        centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree, frame='fk5')
        size = u.Quantity((cfg_par['moments']['sizePlots'],cfg_par['moments']['sizePlots']), u.arcmin)
      
        params = self.loadRcParams()
        plt.rcParams.update(params)


        hduIm = fits.open(outBPT)[0]
        wcsIm = WCS(hduIm.header)

        hduImCut = Cutout2D(hduIm.data, centre, size, wcs=wcsIm)
        
        fig = plt.figure()
        
        ax1 = plt.subplot(projection=wcsIm)    

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size='2%', pad=0.1)
        
        if vRange == None:
          vRange=np.array([1,2])
          vRange[0] = lineThresh
          vRange[1] = np.nanmax(hduImCut.data)

        img = ax1.imshow(hduImCut.data, cmap=cfg_par['lineRatios']['cDistColorMap'],vmin=vRange[0],vmax=vRange[1])

        colorTickLabels = np.linspace(vRange[0],vRange[1],9.)    

        ax1.coords[1].set_axislabel(r'Dec (J2000)')
        ax1.coords[0].set_axislabel(r'RA (J2000)')
        
        cax.coords[0].grid(False)
        cax.coords[1].grid(False)
        cax.tick_params(direction='in')
        cax.coords[0].tick_params(top=False, bottom=False,
                       labeltop=False, labelbottom=False)
        cax.coords[1].set_ticklabel_position('r')
        cax.coords[1].set_axislabel('$\eta$-parameter')
        cax.coords[1].set_axislabel_position('r')
        cbar = plt.colorbar(img, cax=cax,ticks =colorTickLabels,
                        orientation='vertical', format='%d')   
        
        #ax1.set_title('')

        ax1.set_autoscale_on(False)    
        #SaveOutput
        outBPT = os.path.basename(outBPT)
        outBPT= str.split(outBPT, '.fits')[0]  
        modName = cfg_par['gFit']['modName'] 

        if nameFigLabel==None:
            nameFigLabel='' 
        
        if overlayContours:
            imLevels =[-0.5,0.5]
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

            ax1.contour(array.data,levels=contLevels[i], colors=contColors[0])

        outFig = cfg_par['general']['bptDir']+outBPT+nameFigLabel+'.'+cfg_par['moments']['plotFormat']
        fig.savefig(outFig,format=cfg_par['moments']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
                #dpi=300,bbox_inches='tight',transparent=False,overwrite=True)
        return 0  

    def bptCDist(self,cfg_par,outPlotDir=None,vRange=None):


        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data
        lineRatios = hdul['LineRatios_'+cfg_par['gFit']['modName']].data

        ampTable = hdul['LineRes_'+cfg_par['gFit']['modName']].data
        amps = ampTable[cfg_par['gFit']['modName']+'-AmpSpax'+'_Hb4861']

        lineInfo = tP.openLineList(cfg_par)
        lineThresh = float(lineInfo['SNThresh'][0])

        if cfg_par['gFit']['modName'] == 'g1':
            modString = ['G1']
        elif cfg_par['gFit']['modName'] == 'g2':
            modString = ['G1','G2','ToT']
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']

        idx  = np.where(amps>=lineThresh)

        for i in range (0, len(modString)):
            
            # initialize figure
            params = self.loadRcParams()
            plt.rcParams.update(params)
            fig = plt.figure(figsize =(10,8))
            fig.subplots_adjust(hspace=0.0)
            gs = gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0])

            y = lineRatios['log_'+modString[i]+'-OIII5006/Hb4861']
            x = lineRatios['log_'+modString[i]+'-NII6583/Ha6562']
            k = lineBPT['cDist-OIII']
            
            #ax.set_xticks([])
            
            ax1.set_xlabel(r'log([NII] 6583/H$_\alpha$ 6562)')
            ax1.set_ylabel(r'log([OIII]5006/H$_\beta$4861)')


            # Calculate axis limits and aspect ratio
            xMin = -2.
            xMax = 0.8
            yMin = -2.
            yMax = 1.5

            # Set axis limits
            ax1.set_xlim(xMin, xMax)
            ax1.set_ylim(yMin, yMax)

            ax1.tick_params(axis='both', which='major', pad=5)
            #ax1.xaxis.set_minor_locator()
            #ax1.yaxis.set_minor_locator()      
            
            #print((idxAGN),(idxKew),(idxKauf),(idxBad))

            if vRange==None:
                vRangeMin=-1.
                vRangeMax=1.
            else:
                vRangeMin=vRange[0]
                vRangeMax=vRange[1]

            ax1.scatter(x[idx], y[idx], c=k[idx], cmap='nipy_spectral', marker='+', s=80, linewidths=4, 
                label='Fornax A',vmin=vRangeMin,vmax=vRangeMax)

            kaX = np.log10(np.linspace(np.power(10,-1.),np.power(10,0.),1e4))
            kaY = 0.61 / (kaX - 0.05) + 1.3
            ax1.plot(kaX, kaY, ls='--',c='black', label='Kewley et al. 2001')

            keX = np.log10(np.linspace(np.power(10,-2.),np.power(10,0.5),1e4))

            keY= 0.61 / (keX - 0.47) + 1.19
            ax1.plot(keX, keY, ls=':',c='black', label='Kauffmann et al. 2003')


            ax1.legend = plt.legend(loc=1, prop={'size': 12})
            ax1.legend.get_frame().set_edgecolor('black')
            ax1.legend.get_frame().set_facecolor('white')
            if outPlotDir==None:
                outPlotDir = cfg_par['general']['bptDir']+cfg_par['gFit']['modName']+'/plots/'

            if not os.path.exists(outPlotDir):
                os.mkdir(outPlotDir)

            outPlot = outPlotDir+'BPT-'+cfg_par['gFit']['modName']+'etadist.png'
            plt.savefig(outPlot,
                        format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
            plt.show()
            plt.close()
               
        return 0      