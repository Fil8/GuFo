#!/usr/bin/env python
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

import cvPlay
cvP = cvPlay.convert()

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

    def bptOIII(self,cfg_par):


        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data
        lineRatios = hdul['LineRatios_'+cfg_par['gFit']['modName']].data

        # initialize figure
        params = self.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure(figsize =(10,8))
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])

        if cfg_par['gFit']['modName'] == 'g1':
            modString = ['G1']
        elif cfg_par['gFit']['modName'] == 'g2':
            modString = ['G1','G2','ToT']
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']


        for i in xrange (0, len(modString)):
            y = np.log10(lineRatios[modString[i]+'-OIII5006/Hb4861'])
            x = np.log10(lineRatios[modString[i]+'-NII6583/Ha6562'])
            k = lineBPT[modString[i]+'-BPT_OIII']
        #ax.set_xticks([])

        ax1.set_xlabel(r'log([NII] 6583/H$_\alpha$ 6562)')
        ax1.set_ylabel(r'log([OIII]5006/H$_\beta$4861)')


        # Calculate axis limits and aspect ratio
        xMin = -2.
        xMax = 2.
        yMin = -2.
        yMax = 4

        # Set axis limits
        ax1.set_xlim(xMin, xMax)
        ax1.set_ylim(yMin, yMax)

        ax1.tick_params(axis='both', which='major', pad=5)
        #ax1.xaxis.set_minor_locator()
        #ax1.yaxis.set_minor_locator()      
        
        idxAGN = np.where(k==2.)
        idxKauf = np.where(k==0.)
        idxKew = np.where(k==1.)


        ax1.scatter(x[idxAGN], y[idxAGN], c='red', marker='+', s=80, linewidths=4)
        ax1.scatter(x[idxKew], y[idxKew], c='cyan', marker='+', s=80, linewidths=4)
        ax1.scatter(x[idxKauf], y[idxKauf], c='blue', marker='+', s=80, linewidths=4)

        kaX = np.log10(np.linspace(np.power(10,-1.),np.power(10,0.),1e4))
        kaY = 0.61 / (kaX - 0.05) + 1.3
        ax1.plot(kaX, kaY, ls='--',c='black', label='Kewley et al. 2001')

        keX = np.log10(np.linspace(np.power(10,-2.),np.power(10,0.5),1e4))

        keY= 0.61 / (keX - 0.47) + 1.19
        ax1.plot(keX, keY, ls=':',c='black', label='Kauffmann et al. 2003')

    
        #ax1.text(xText, y1_max*0.90, r'BIN ID:\t'+str(singleVorBinInfo['BIN_ID'][0]), {'color': 'k', 'fontsize': 20})
        #ax1.text(xText, y1_max*0.94, r'X,Y:\t'+str(xx)+','+str(yy), {'color': 'k', 'fontsize': 20})

        #ax1.text(xText, x_max*0.85, r'Success:\t'+successStr, {'color': 'b'})
        #ax1.text(xText, y1_max*0.88, r'$\tilde{\chi}^2$:\t'+redchiStr, {'color': 'k', 'fontsize': 20})
        #ax1.text(xText, y1_max*0.82, r'aic:\t'+aicStr, {'color': 'k', 'fontsize': 20})
        #ax1.text(xText, y1_max*0.76, r'bic:\t'+bicStr, {'color': 'k', 'fontsize': 20})

        #ax1.fill_between(vel, yBFit-dely, yBFit+dely, color="#ABABAB",
        #        label='3-$\sigma$ uncertainty band')
        #if cfg_par['gFit']['modName'] !='g1':
        #    comps = result.eval_components()
        #    for i in xrange(0,len(lineInfo['ID'])):
                
        #        ax1.plot(velPlot, comps['g1ln'+str(i)+'_'], 'g--')
            
        #        if cfg_par['gFit']['modName'] =='g2':
        #            ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
            
        #        elif cfg_par['gFit']['modName'] !='g2':
        #            ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
        #            ax1.plot(velPlot, comps['g3ln'+str(i)+'_'], 'c--')    

        outPlotDir = cfg_par['general']['bptDir']+cfg_par['gFit']['modName']+'/'

        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)

        outPlot = outPlotDir+modString[i]+'OIII.png'
        plt.savefig(outPlot,
                    format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        plt.show()
        plt.close()
           
        return 0

    def bptSII(self,cfg_par):


        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data
        lineRatios = hdul['LineRatios_'+cfg_par['gFit']['modName']].data

        # initialize figure
        params = self.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure(figsize =(10,8))
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])

        if cfg_par['gFit']['modName'] == 'g1':
            modString = ['G1']
        elif cfg_par['gFit']['modName'] == 'g2':
            modString = ['G1','G2','ToT']
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']


        for i in xrange (0, len(modString)):
            y = np.log10(lineRatios[modString[i]+'-OIII5006/Hb4861'])
            x = np.log10(lineRatios[modString[i]+'-SII6716/Ha6562'])
            k = lineBPT[modString[i]+'-BPT_SII']
        #ax.set_xticks([])

        ax1.set_xlabel(r'log([SII] 6716,6730/H$_\alpha$ 6562)')
        ax1.set_ylabel(r'log([OIII]5006/H$_\beta$4861)')


        # Calculate axis limits and aspect ratio
        xMin = -2.
        xMax = 2.
        yMin = -2.
        yMax = 4

        # Set axis limits
        ax1.set_xlim(xMin, xMax)
        ax1.set_ylim(yMin, yMax)

        ax1.tick_params(axis='both', which='major', pad=5)
        #ax1.xaxis.set_minor_locator()
        #ax1.yaxis.set_minor_locator()      
        
        idxLin = np.where(k==2.)
        idxSey = np.where(k==1.)
        idxKew = np.where(k==0.)


        ax1.scatter(x[idxLin], y[idxLin], c='green', marker='+', s=80, linewidths=4)
        ax1.scatter(x[idxKew], y[idxKew], c='blue', marker='+', s=80, linewidths=4)
        ax1.scatter(x[idxSey], y[idxSey], c='red', marker='+', s=80, linewidths=4)

        keX = np.log10(np.linspace(np.power(10,-2.),np.power(10,0.5),1e4))
        keY= 0.72 / (keX - 0.32) + 1.30
        ax1.plot(keX, keY, ls=':',c='black', label='Kauffmann et al. 2003')
        
        seyX = np.log10(np.linspace(np.power(10,-0.4),np.power(10,0.5),1e4))
        seyLine = 1.89*seyX + 0.76
        ax1.plot(seyX, seyLine, ls=':',c='black', label='Seyfert-LINER')
    
        #ax1.text(xText, y1_max*0.90, r'BIN ID:\t'+str(singleVorBinInfo['BIN_ID'][0]), {'color': 'k', 'fontsize': 20})
        #ax1.text(xText, y1_max*0.94, r'X,Y:\t'+str(xx)+','+str(yy), {'color': 'k', 'fontsize': 20})

        #ax1.text(xText, x_max*0.85, r'Success:\t'+successStr, {'color': 'b'})
        #ax1.text(xText, y1_max*0.88, r'$\tilde{\chi}^2$:\t'+redchiStr, {'color': 'k', 'fontsize': 20})
        #ax1.text(xText, y1_max*0.82, r'aic:\t'+aicStr, {'color': 'k', 'fontsize': 20})
        #ax1.text(xText, y1_max*0.76, r'bic:\t'+bicStr, {'color': 'k', 'fontsize': 20})

        #ax1.fill_between(vel, yBFit-dely, yBFit+dely, color="#ABABAB",
        #        label='3-$\sigma$ uncertainty band')
        #if cfg_par['gFit']['modName'] !='g1':
        #    comps = result.eval_components()
        #    for i in xrange(0,len(lineInfo['ID'])):
                
        #        ax1.plot(velPlot, comps['g1ln'+str(i)+'_'], 'g--')
            
        #        if cfg_par['gFit']['modName'] =='g2':
        #            ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
            
        #        elif cfg_par['gFit']['modName'] !='g2':
        #            ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
        #            ax1.plot(velPlot, comps['g3ln'+str(i)+'_'], 'c--')    

        outPlotDir = cfg_par['general']['bptDir']+cfg_par['gFit']['modName']+'/'

        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)

        outPlot = outPlotDir+modString[i]+'SII.png'
        plt.savefig(outPlot,
                    format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        plt.show()
        plt.close()
           
        return 0

    def bptOI(self,cfg_par):


        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data
        lineRatios = hdul['LineRatios_'+cfg_par['gFit']['modName']].data

        # initialize figure
        params = self.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure(figsize =(10,8))
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])

        if cfg_par['gFit']['modName'] == 'g1':
            modString = ['G1']
        elif cfg_par['gFit']['modName'] == 'g2':
            modString = ['G1','G2','ToT']
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']


        for i in xrange (0, len(modString)):
            y = np.log10(lineRatios[modString[i]+'-OIII5006/Hb4861'])
            x = np.log10(lineRatios[modString[i]+'-OI6300/Ha6562'])
            k = lineBPT[modString[i]+'-BPT_OI']
        #ax.set_xticks([])

        ax1.set_xlabel(r'log([OI] 6300/H$_\alpha$ 6562)')
        ax1.set_ylabel(r'log([OIII]5006/H$_\beta$4861)')


        # Calculate axis limits and aspect ratio
        xMin = -4.
        xMax = 6.
        yMin = -2.
        yMax = 12

        # Set axis limits
        ax1.set_xlim(xMin, xMax)
        ax1.set_ylim(yMin, yMax)

        ax1.tick_params(axis='both', which='major', pad=5)
        #ax1.xaxis.set_minor_locator()
        #ax1.yaxis.set_minor_locator()      
        
        idxLin = np.where(k==2.)
        idxSey = np.where(k==1.)
        idxKew = np.where(k==0.)


        ax1.scatter(x[idxLin], y[idxLin], c='green', marker='+', s=80, linewidths=4)
        ax1.scatter(x[idxKew], y[idxKew], c='blue', marker='+', s=80, linewidths=4)
        ax1.scatter(x[idxSey], y[idxSey], c='red', marker='+', s=80, linewidths=4)


        kaX = np.log10(np.linspace(np.power(10,-3.),np.power(10,-0.595),1e4))
        kaY = 0.73 / (kaX + 0.59) + 1.33
        ax1.plot(kaX, kaY, ls='--',c='black', label='Kewley et al. 2001')

        seyX = np.log10(np.linspace(np.power(10,-1.),np.power(10,0.0),1e4))
        seyLine = 1.18*seyX + 1.30
        ax1.plot(seyX, seyLine, ls=':',c='black', label='Seyfert-LINER')

    
        #ax1.text(xText, y1_max*0.90, r'BIN ID:\t'+str(singleVorBinInfo['BIN_ID'][0]), {'color': 'k', 'fontsize': 20})
        #ax1.text(xText, y1_max*0.94, r'X,Y:\t'+str(xx)+','+str(yy), {'color': 'k', 'fontsize': 20})

        #ax1.text(xText, x_max*0.85, r'Success:\t'+successStr, {'color': 'b'})
        #ax1.text(xText, y1_max*0.88, r'$\tilde{\chi}^2$:\t'+redchiStr, {'color': 'k', 'fontsize': 20})
        #ax1.text(xText, y1_max*0.82, r'aic:\t'+aicStr, {'color': 'k', 'fontsize': 20})
        #ax1.text(xText, y1_max*0.76, r'bic:\t'+bicStr, {'color': 'k', 'fontsize': 20})

        #ax1.fill_between(vel, yBFit-dely, yBFit+dely, color="#ABABAB",
        #        label='3-$\sigma$ uncertainty band')
        #if cfg_par['gFit']['modName'] !='g1':
        #    comps = result.eval_components()
        #    for i in xrange(0,len(lineInfo['ID'])):
                
        #        ax1.plot(velPlot, comps['g1ln'+str(i)+'_'], 'g--')
            
        #        if cfg_par['gFit']['modName'] =='g2':
        #            ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
            
        #        elif cfg_par['gFit']['modName'] !='g2':
        #            ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
        #            ax1.plot(velPlot, comps['g3ln'+str(i)+'_'], 'c--')    

        outPlotDir = cfg_par['general']['bptDir']+cfg_par['gFit']['modName']+'/'

        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)

        outPlot = outPlotDir+modString[i]+'OI.png'
        plt.savefig(outPlot,
                    format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        plt.show()
        plt.close()
           
        return 0