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
import aplpy
import cvPlay
cvP = cvPlay.convert()

class MOMplot(object):

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

    def mom0Plot(self, cfg_par,imageName):
        
        params = self.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure()

        f1 = aplpy.FITSFigure(imageName, figure=fig)
        f1.set_theme('publication')

        f1.frame.set_linewidth(2)

        f1.show_colorscale(aspect='equal', cmap='afmhot',stretch = 'linear')

        f=fits.open(imageName)
        dd=f[0].data
        
        #print(outBPT)
        #print(idxLin,idxSey,idxKew,idxBad)

        #f1.show_contour(imageName,levels=[1, 5, 8, 11, 15], colors='black')
        #f1.show_contour(imageName,levels=[1e-3, 1e-2, 1e-1,5e-1], colors='black')

        f1.axis_labels.set_font( weight='book',size='medium', 
                                 stretch='normal', family='serif', 
                                 style='normal', variant='normal')
        f1.axis_labels.set_xtext('RA (J2000)')
        f1.axis_labels.set_ytext('Dec (J2000)')
        f1.tick_labels.set_xformat('hh:mm:ss')
        f1.tick_labels.set_yformat('dd:mm:ss')
        f1.tick_labels.set_font( weight='book', size='small',
                                 stretch='normal', family='serif', 
                                 style='normal', variant='normal') 
        f1.ticks.set_color('k')
        f1.ticks.set_length(6)  # points
        f1.ticks.set_linewidth(2)  # points
        f1.ticks.set_minor_frequency(3)
        f1.ticks.show()

        outMom = os.path.basename(imageName)
        print(imageName)

        outMom= str.split(outMom, '.fits')[0]  
        modName = cfg_par['gFit']['modName']
        
        outMom = cfg_par['general']['momDir']+modName+'/plots/'+outMom+'.'+cfg_par['moments']['plotFormat']
        
        if os.path.exists(cfg_par['general']['momDir']+modName+'/plots/') == False:
            os.mkdir(cfg_par['general']['momDir']+modName+'/plots/')

        fig.savefig(outMom,format=cfg_par['moments']['plotFormat'])
                #if pdf, dpi=300,bbox_inches='tight',transparent=False,overwrite=True)

        return 0

    def mom1Plot(self,cfg_par, imageName, cenRange):


        params = self.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure()

        f1 = aplpy.FITSFigure(imageName, figure=fig)
        f1.set_theme('publication')

        f1.frame.set_linewidth(2)

        f1.show_colorscale(aspect='equal', cmap='nipy_spectral',stretch = 'linear',vmin= -cenRange, vmax= cenRange)

        f=fits.open(imageName)
        dd=f[0].data
        
        #print(outBPT)
        #print(idxLin,idxSey,idxKew,idxBad)

        #f1.show_contour(imageName,levels=[1, 5, 8, 11, 15], colors='black')
        #f1.show_contour(imageName,levels=[1e-3, 1e-2, 1e-1,5e-1], colors='black')

        f1.axis_labels.set_font( weight='book',size='medium', 
                                 stretch='normal', family='serif', 
                                 style='normal', variant='normal')
        f1.axis_labels.set_xtext('RA (J2000)')
        f1.axis_labels.set_ytext('Dec (J2000)')
        f1.tick_labels.set_xformat('hh:mm:ss')
        f1.tick_labels.set_yformat('dd:mm:ss')
        f1.tick_labels.set_font( weight='book', size='small',
                                 stretch='normal', family='serif', 
                                 style='normal', variant='normal') 
        f1.ticks.set_color('k')
        f1.ticks.set_length(6)  # points
        f1.ticks.set_linewidth(2)  # points
        f1.ticks.set_minor_frequency(3)
        f1.ticks.show()

        outMom = os.path.basename(imageName)
        outMom= str.split(outMom, '.fits')[0]  
        modName = cfg_par['gFit']['modName']
        
        outMom = cfg_par['general']['momDir']+modName+'/plots/'+outMom+'.'+cfg_par['moments']['plotFormat']
        
        if os.path.exists(cfg_par['general']['momDir']+modName+'/plots/') == False:
            os.mkdir(cfg_par['general']['momDir']+modName+'/plots/')

        fig.savefig(outMom,format=cfg_par['moments']['plotFormat'])
                #if pdf, dpi=300,bbox_inches='tight',transparent=False,overwrite=True)

        return 0