#!/usr/bin/env python3.6
'''

Set of tools for generating different versions of the kinematical plot (sigma vs. centroid)
and compare the values of sigma and centroid of each fitted line with the expectations of Cold Chaotic Accretion.

Requirements
------------
Fit solutions must have been previously generated. Sigma and centroid of the integrated line must have been computed
  - `gPlay.py` or `gPlayMp.py` : fit solutions
  - `ancels.py`: measure sigma and centroid
'''

import os, sys
from astropy.io import fits
from astropy.table import Table, Column

import numpy as np

import pandas as pd

from matplotlib.patches import Ellipse

from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator, FixedLocator
from matplotlib import transforms as mtransforms
from matplotlib.ticker import LogFormatter 
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
import matplotlib.axes as maxes

import matplotlib.ticker as mticker


import cvPlay, tPlay, momPlot
from scavengers import util as ut

cvP = cvPlay.convert()
tP = tPlay.tplay()
mPl = momPlot.MOMplot()

class ancelsplot():
    '''This class draws different versions of the kinematical plot (centroid vs sigma), 
    for ionised gas (MUSE observations) and cold gas (interferometric observations)

    '''

    def sigmaCentroid(self,cfg_par,outPlotDir=None):
        '''Draws the k-plot from the log_dispIntrXX and log_centroidXX columns of the subtable ancels_YY of ionised gas lines.

            - XX = cfg_par['kinematicalAnalysis']['Name'] : ionised gas emission line name
            - YY = cfg_par['modName']['modName'] : fit ID (g1, g2 or BF)
        
        Other variables are set in section `kinematicalAnalysis` of the parameter file.

        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        outPlotDir: str, optional
            output directory of plot

        '''


        lineInfo = tP.openLineList(cfg_par)

        lineThresh = float(lineInfo['SNThresh'][0])


        hdul = fits.open(cfg_par['general']['outTableName'])
        ancels = hdul['Ancels'+cfg_par['gFit']['modName']].data
        if cfg_par['gFit']['modName'] == 'BF':
                cfg_par['gFit']['modName'] = 'g2'
        lines = hdul['LineRes_'+cfg_par['gFit']['modName']].data
        residuals = hdul['Residuals_'+cfg_par['gFit']['modName']].data
        
        SN = residuals['SN_NII6583']
        BinID = np.argsort(residuals['Bin_ID'])
        SNSort = SN[BinID]

        g1Sort = lines['g1_sigIntr_NII6583'][BinID]

        if cfg_par['gFit']['modName'] == 'g1':
            modString = ['G1']
            x = [np.log10(lines['g1_Centre_NII5683'])]
            y = [np.log10(lines['g1_sigIntr_NII5683'])]

        elif cfg_par['gFit']['modName'] == 'g2':
            modString = ['G1','G2','ToT']

            binCode = [np.argsort(lines['BIN_ID']),np.argsort(ancels['BIN_ID'])]

            x = [np.log10(lines['g1_Centre_NII6583'][binCode[0]]),np.log10(lines['g2_Centre_NII6583'][binCode[0]]),
            ancels['logCentroid_NII6583'][binCode[1]]]
            
            y = [np.log10(lines['g1_sigIntr_NII6583'][binCode[0]]),np.log10(lines['g2_sigIntr_NII6583'][binCode[0]]),
            ancels['logDispIntr_NII6583'][binCode[1]]]
            CCASca=ancels['CCAIN'][binCode[1]]
            rotSca = ancels['RotMod'][binCode[1]]
        
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']
            x = [np.log10(lines['g1_Centre_NII6583']),np.log10(lines['g2_Centre_NII6583']),np.log10(lines['g3_Centre_NII6583']),
            ancels['logCentroid_NII6583']]

            binCode = [lines['BIN_ID'],lines['BIN_ID'],lines['BIN_ID'],ancels['BIN_ID']]
        
            y = [np.log10(lines['g1_sigIntr_NII6583']),np.log10(lines['g2_sigIntr_NII6583']),np.log10(lines['g3_sigIntr_NII6583']),
            ancels['logDispIntr_NII6583']]
        
        cfg_par['gFit']['modName'] = 'BF'
        for i in range (0, len(modString)):
            
            # initialize figure
            params = ut.loadRcParams()
            plt.rcParams.update(params)
            fig = plt.figure(figsize =(10,8),constrained_layout=False)
            fig.set_tight_layout(False)
            fig.subplots_adjust(hspace=0.0)
            gs = gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0])

        
            #ax.set_xticks([])
            
            ax1.set_xlabel(r'log($|v_{\rm los}-v_{\rm sys}|$)\,\, [km s$^{-1}$]')
            ax1.set_ylabel(r'log($\sigma_{\rm los}$)\,\, [km s$^{-1}$]')
        
            # Calculate axis limits and aspect ratio
            xMin = 0.
            xMax = 2.81
            yMin = 0.2
            yMax = 3.12
            
            ax1.minorticks_on()
            ax1.xaxis.set_major_locator(MultipleLocator(0.5))
            ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
            ax1.yaxis.set_major_locator(MultipleLocator(0.5))
            ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
            ax1.tick_params(axis='both', which='major', pad=5)
            ax1.tick_params(axis='both', which='minor', pad=2)

            if cfg_par['kinematicalAnalysis']['ancillaryInfo']['theoreticalCCA'] == 'Ensemble':
                Mean_sigmav = 2.13 
                RMS_sigmav  = 0.13 
                Mean_vshift = 1.59 
                RMS_vshift  = 0.37
                theta = -179.88 #covariance angle (contours inclined)
                ellColor='purple'
            if cfg_par['kinematicalAnalysis']['ancillaryInfo']['theoreticalCCA'] == 'Pencil':
                Mean_sigmav = 1.65 
                RMS_sigmav  = 0.41
                Mean_vshift = 2.0
                RMS_vshift  = 0.47
                theta = 165.16 #covariance angle (contours inclined)
                ellColor='purple'

            idx = SNSort>=lineThresh

            #x[i][idx] = np.nan
            #y[i][idx] = np.nan
            if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotTheoreticalCCA'] == True:
                rmsToFWHM = 2.*np.sqrt(2.*np.log(2))
                ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=rmsToFWHM*RMS_vshift, height=rmsToFWHM*RMS_sigmav, angle=theta,
                    label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['theoreticalCCA']+' beam')     
                ellSigma1.set_clip_box(ax1.bbox)
                ellSigma1.set_alpha(0.2)
                ellSigma1.set_facecolor(ellColor)
                ax1.add_artist(ellSigma1)

                ellip, ellip_lbl = ax1.get_legend_handles_labels()

                ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=2*rmsToFWHM*RMS_vshift, height=2*rmsToFWHM*RMS_sigmav, angle=theta)     
                ellSigma1.set_clip_box(ax1.bbox)
                ellSigma1.set_alpha(0.3)
                ellSigma1.set_facecolor(ellColor)
                ax1.add_artist(ellSigma1)
                
                #ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=3*rmsToFWHM*RMS_vshift, height=3*rmsToFWHM*RMS_sigmav, angle=theta)     
                #ellSigma1.set_clip_box(ax1.bbox)
                #ellSigma1.set_alpha(0.2)
                #ellSigma1.set_facecolor(ellColor)
                #ax1.add_artist(ellSigma1)
            
            #ax1.scatter(x[i][idx], y[i][idx], c='red', marker='.', s=20, linewidths=None,edgecolors='red',
            #   label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel'][i])

            if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotRotation'] == True:

                indexNoRot = np.logical_and(rotSca ==0.,SNSort>=lineThresh)
                #print(indexRot)
                ax1.scatter(x[i][indexNoRot], y[i][indexNoRot], c='seagreen', marker='.', s=20, linewidths=None, alpha=0.3,facecolors='seagreen',edgecolors='green',
                    label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel'][i])

                indexRot = np.logical_and(rotSca == 1.,SNSort>=lineThresh)
                #print(indexRot)
                ax1.scatter(x[i][indexRot], y[i][indexRot], c='blue', marker='.', s=20, linewidths=None, alpha=0.3,facecolors='blue',edgecolors='blue',
                    label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel'][i]+' in rotation')



            if cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCAanalysis'] == True and (cfg_par['kinematicalAnalysis']['ancillaryInfo']['rotOutCCA'] ==True or cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotElse'] ==True):
                
                if cfg_par['kinematicalAnalysis']['ancillaryInfo']['rotOutCCA'] ==True:
                    indexCCA = np.logical_and(np.logical_and(np.logical_and(rotSca!=1.,np.logical_and(CCASca!=np.nan,SNSort>=lineThresh)),
                        CCASca!=-1.),CCASca <= 1.)
                    indexElse = np.logical_and(np.logical_and(rotSca!=1.,np.logical_and(CCASca!=np.nan,SNSort>=lineThresh)),CCASca == -1.)

                    ax1.scatter(x[i][indexCCA], y[i][indexCCA], c='seagreen', marker='.', s=20, linewidths=None, alpha=0.5,facecolors='seagreen',edgecolors=None,
                        label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel'][i]+' CCA within '+
                        str(round(cfg_par['kinematicalAnalysis']['ancillaryInfo']['sigmaInCCA'],1))+r'$\sigma$')                
                else:   
                    indexCCA = np.logical_and(CCASca <= 1.,np.logical_and(np.logical_and(CCASca!=np.nan,SNSort>=lineThresh),CCASca!=-1.))
                    indexElse = np.logical_and(CCASca == -1.,np.logical_and(CCASca!=np.nan,SNSort>=lineThresh))

                    ax1.scatter(x[i][indexCCA], y[i][indexCCA], c='seagreen', marker='.', s=20, linewidths=None, alpha=0.5,facecolors='seagreen',edgecolors=None,
                        label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel'][i]+' CCA within '+
                        str(round(cfg_par['kinematicalAnalysis']['ancillaryInfo']['sigmaInCCA'],1))+r'$\sigma$')     
                
                if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotElse'] ==True:
                    ax1.scatter(x[i][indexElse], y[i][indexElse], c='darkgray', marker='.', s=20, linewidths=None, alpha=0.2,facecolors='darkgray',edgecolors=None,
                        label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel'][i]+' unidentified')
                #ax1.scatter(x[indexCCA], y[indexCCA], c='red', marker='.', s=20, linewidths=None, alpha=0.5,
                #    label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel'][i]+' CCA outliers')
            
            # Set axis limits
            ax1.set_xlim(xMin, xMax)
            ax1.set_ylim(yMin, yMax)
            lgnd = plt.legend(loc=3,prop={'size':18},handletextpad=-0.3,borderaxespad=0.)
            
            lgnd.legendHandles[0]._sizes = [200]
            lgnd.legendHandles[1]._sizes = [200]

            
            #change the marker size manually for both lines
            #ax1.legend._legmarker.set_markersize(6)
            #ax1.legend(handletextpad=0.1,markersize(6))
            if outPlotDir==None:
                outPlotDir = cfg_par['general']['plotMomModDir']

            if not os.path.exists(outPlotDir):
                os.mkdir(outPlotDir)

            if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotRotation'] == True and cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCAanalysis'] == True:
                outPlot = outPlotDir+'sigmaCentroid-'+modString[i]+'-'+cfg_par['kinematicalAnalysis']['Name']+'RotCCA-'+str(int(cfg_par['kinematicalAnalysis']['ancillaryInfo']['sigmaInCCA']))+'sigma.png'
            elif cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotRotation'] == False and cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCAanalysis'] == True:
                outPlot = outPlotDir+'sigmaCentroid-'+modString[i]+'-'+cfg_par['kinematicalAnalysis']['Name']+'CCA-'+str(int(cfg_par['kinematicalAnalysis']['ancillaryInfo']['sigmaInCCA']))+'sigma.png'
            else:
                outPlot = outPlotDir+'sigmaCentroid-'+modString[i]+'-'+cfg_par['kinematicalAnalysis']['Name']+'.pdf'

            plt.savefig(outPlot,format=cfg_par['lineRatios']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=300)#,
            print(outPlot)
                    # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
            #plt.show()
            plt.close()
               
        return 0

    def sigmaCentroidColdGas(self,cfg_par):
        '''Draws the k-plot from the log_dispIntrXX and log_centroidXX columns of the subtable ancels_YY of cold gas lines
        (interferometric observations)

            - XX = cfg_par['otherGasKinAnalysis']['Name'] : cold gas emission line name
            - YY = cfg_par['modName']['modName'] : fit ID (g1, g2 or BF)
        
        Other variables are set in section `otherGasKinAnalysis` of the parameter file.

        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        outPlotDir: str, optional
            output directory of plot

        '''



        hdul = fits.open(cfg_par['general']['outTableName'])

        ancels = hdul['Ancels'+cfg_par['gFit']['modName']].data

        x = ancels['logCentroid_'+cfg_par['otherGasKinAnalysis']['Name']]
        y = ancels['logSigma_'+cfg_par['otherGasKinAnalysis']['Name']]
            
        CCASca=ancels['CCAIN']
        rotSca = ancels['RotMod']
        # initialize figure
        
        params = ut.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure(figsize =(10,8),constrained_layout=False)
        fig.set_tight_layout(False)
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])

        #ax.set_xticks([])
        
        ax1.set_xlabel(r'$\log(|v_{\rm los}-v_{\rm sys}|)$\,\, [km s$^{-1}$]')
        ax1.set_ylabel(r'$\log(\sigma_{\rm los}$)\,\, [km s$^{-1}$]')
    
        # Calculate axis limits and aspect ratio
        xMin = 0.
        xMax = 2.9
        yMin = 0.2
        yMax = 3.2
        
        ax1.minorticks_on()
        ax1.xaxis.set_major_locator(MultipleLocator(0.5))
        ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.yaxis.set_major_locator(MultipleLocator(0.5))
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.tick_params(axis='both', which='major', pad=5)
        ax1.tick_params(axis='both', which='minor', pad=2)

        if cfg_par['otherGasKinAnalysis']['ancillaryInfo']['theoreticalCCA'] == 'Ensemble':
            Mean_sigmav = 2.13 
            RMS_sigmav  = 0.13 
            Mean_vshift = 1.59 
            RMS_vshift  = 0.37
            theta = -179.88 #covariance angle (contours inclined)
            ellColor='purple'
        if cfg_par['otherGasKinAnalysis']['ancillaryInfo']['theoreticalCCA'] == 'Pencil':
            Mean_sigmav = 1.65 
            RMS_sigmav  = 0.41
            Mean_vshift = 2.0
            RMS_vshift  = 0.47
            theta = 165.16 #covariance angle (contours inclined)
            ellColor='darkseagreen'
        
        if cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotTheoreticalCCA'] == True:
            rmsToFWHM = 2.*np.sqrt(2.*np.log(2))
            ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=rmsToFWHM*RMS_vshift, height=rmsToFWHM*RMS_sigmav, angle=theta,
                label=cfg_par['otherGasKinAnalysis']['ancillaryInfo']['theoreticalCCA']+' beam')     
            ellSigma1.set_clip_box(ax1.bbox)
            ellSigma1.set_alpha(0.4)
            ellSigma1.set_facecolor(ellColor)
            ax1.add_artist(ellSigma1)

            ellip, ellip_lbl = ax1.get_legend_handles_labels()

            ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=2*rmsToFWHM*RMS_vshift, height=2*rmsToFWHM*RMS_sigmav, angle=theta)     
            ellSigma1.set_clip_box(ax1.bbox)
            ellSigma1.set_alpha(0.3)
            ellSigma1.set_facecolor(ellColor)
            ax1.add_artist(ellSigma1)
            
            ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=3*rmsToFWHM*RMS_vshift, height=3*rmsToFWHM*RMS_sigmav, angle=theta)     
            ellSigma1.set_clip_box(ax1.bbox)
            ellSigma1.set_alpha(0.2)
            ellSigma1.set_facecolor(ellColor)
            ax1.add_artist(ellSigma1)
        
        #print((idxAGN),(idxKew),(idxKauf),(idxBad))
        if cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotRotation'] == True:

            indexNoRot = rotSca ==0.
            #print(indexRot)
            ax1.scatter(x[indexNoRot], y[indexNoRot], c='seagreen', marker='.', s=50, linewidths=None, alpha=0.8,facecolors='seagreen',edgecolors='green',
                label=cfg_par['otherGasKinAnalysis']['ancillaryInfo']['CCALabel'])

            indexRot = rotSca == 1.
            #print(indexRot8
            ax1.scatter(x[indexRot], y[indexRot], c='blue', marker='.', s=50, linewidths=None, alpha=0.8,facecolors='blue',edgecolors='blue',
                label=cfg_par['otherGasKinAnalysis']['ancillaryInfo']['CCALabel']+' in rotation')



        if cfg_par['otherGasKinAnalysis']['ancillaryInfo']['CCAanalysis'] == True and (cfg_par['otherGasKinAnalysis']['ancillaryInfo']['rotOutCCA'] ==True or cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotElse'] ==True):
                
            
            if cfg_par['otherGasKinAnalysis']['ancillaryInfo']['rotOutCCA'] ==True:
                
                indexCCA = np.logical_and(np.logical_and(np.logical_and(rotSca!=1.,CCASca!=np.nan),
                    CCASca!=-1.),CCASca <= 1.)
                indexElse = np.logical_and(np.logical_and(rotSca!=1.,CCASca!=np.nan),CCASca == -1.)

                ax1.scatter(x[indexCCA], y[indexCCA], c='seagreen', marker='.', s=50, linewidths=None, alpha=0.8,facecolors='seagreen',edgecolors=None,
                    label=cfg_par['otherGasKinAnalysis']['ancillaryInfo']['CCALabel']+' CCA within '+
                    str(int(cfg_par['otherGasKinAnalysis']['ancillaryInfo']['sigmaInCCA']))+r'$\sigma$')
            
            else:   
                indexCCA = np.logical_and(CCASca <= 1.,np.logical_and(CCASca!=np.nan,CCASca!=-1.))
                indexElse = np.logical_and(CCASca == -1.,CCASca!=np.nan)
                
                ax1.scatter(x[indexCCA], y[indexCCA], c='seagreen', marker='.', s=50, linewidths=None, alpha=0.8,facecolors='seagreen',edgecolors=None,
                    label=cfg_par['otherGasKinAnalysis']['ancillaryInfo']['CCALabel']+' CCA within '+
                    str(int(cfg_par['otherGasKinAnalysis']['ancillaryInfo']['sigmaInCCA']))+r'$\sigma$')
                
            if cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotElse'] ==True:
                ax1.scatter(x[indexElse], y[indexElse], c='darkgray', marker='.', s=20, linewidths=None, alpha=0.8,facecolors='darkgray',edgecolors=None,
                    label=cfg_par['otherGasKinAnalysis']['ancillaryInfo']['CCALabel']+' unidentified')



        # Set axis limits
        ax1.set_xlim(xMin, xMax)
        ax1.set_ylim(yMin, yMax)
        ax1.legend = plt.legend(loc=1)
        ax1.legend.get_frame().set_edgecolor('black')
        ax1.legend.get_frame().set_facecolor('white')


        outPlotDir = cfg_par['general']['plotMomModDir']

        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)

        outPlot = outPlotDir+'sigmaCentroid-g1-'+cfg_par['otherGasKinAnalysis']['Name']+'.pdf'
        plt.savefig(outPlot,format=cfg_par['lineRatios']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=300)#,
                # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        #plt.show()
        #plt.close()
               
        return 0   

    def sigmaCentroidMultiple(self,cfg_par):
        '''Draws the k-plot from the log_dispIntrXX and log_centroidXX columns from multiple tables given in cfg_par['multipleRegions']['tableNames']
        of the parameter file.
        
        Other variables are set in section `otherGasKinAnalysis` of the parameter file.

        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        '''

        tableNames = cfg_par['multipleRegions']['tableNames']
        colorNames = cfg_par['multipleRegions']['colors']
        regionNames = cfg_par['multipleRegions']['regions']


        params = ut.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure(figsize=(9,7.24409),constrained_layout=False)
        fig.set_tight_layout(False)
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])

        #ax1.set_xticks([])
        
        ax1.set_xlabel(r'$\log(|v_{\rm los}-v_{\rm sys}|)$\,\, [km s$^{-1}$]')
        ax1.set_ylabel(r'$\log(\sigma_{\rm los}$)\,\, [km s$^{-1}$]')
        
        #ax1.xaxis.set_ticklabels([])
        # Calculate axis limits and aspect ratio
        xMin = 0.
        xMax = 2.9
        yMin = 0.2
        yMax = 3.2
        
        ax1.minorticks_on()
        ax1.xaxis.set_major_locator(MultipleLocator(0.5))
        ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.yaxis.set_major_locator(MultipleLocator(0.5))
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.tick_params(axis='both', which='major', pad=5)
        ax1.tick_params(axis='both', which='minor', pad=2)

        ax1.tick_params(top=True)
        ax1.tick_params(right=True)
        ax1.tick_params(top=True, which='minor')
        ax1.tick_params(right=True, which='minor')

        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['theoreticalCCA'] == 'Ensemble':
            Mean_sigmav = 2.13 
            RMS_sigmav  = 0.13 
            Mean_vshift = 1.59 
            RMS_vshift  = 0.37
            theta = -179.88 #covariance angle (contours inclined)
            ellColor='purple'
        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['theoreticalCCA'] == 'Pencil':
            Mean_sigmav = 1.65 
            RMS_sigmav  = 0.41
            Mean_vshift = 2.0
            RMS_vshift  = 0.47
            theta = 165.16 #covariance angle (contours inclined)
            ellColor='purple'
        



        for i in range(0,len(tableNames)):
            
            hdul = fits.open(cfg_par['general']['runNameDir']+tableNames[i])
            ancels = hdul[1].data

            x = ancels['logCentroid_'+cfg_par['multipleRegions']['Name']]
            if cfg_par['multipleRegions']['Name'] == 'HI' or cfg_par['multipleRegions']['Name'] == 'CO':

                y = ancels['logSigma_'+cfg_par['multipleRegions']['Name']]
            else:
                y = ancels['logDispIntr_'+cfg_par['multipleRegions']['Name']]
            # initialize figure
            #print((idxAGN),(idxKew),(idxKauf),(idxBad))
            ax1.scatter(x, y, c=colorNames[i],facecolors=colorNames[i],edgecolors=colorNames[i],
                marker='.', s=20, linewidths=None, alpha=0.3,label=regionNames[i])
            hdul.close()

        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotTheoreticalCCA'] == True or cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotTheoreticalCCA'] == True:
            rmsToFWHM = 2.*np.sqrt(2.*np.log(2))
            ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=rmsToFWHM*RMS_vshift, height=rmsToFWHM*RMS_sigmav, angle=theta,
                label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['theoreticalCCA']+' beam',fill=False)     
            ellSigma1.set_clip_box(ax1.bbox)
            ellSigma1.set_alpha(1)
            # ellSigma1.set_facecolor('white')

            ellSigma1.set_edgecolor(ellColor)            
            ellSigma1.set_linewidth(2)            

            ax1.add_artist(ellSigma1)

            ellip, ellip_lbl = ax1.get_legend_handles_labels()

            ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=2*rmsToFWHM*RMS_vshift, height=2*rmsToFWHM*RMS_sigmav, 
                angle=theta,fill=False)     
            ellSigma1.set_clip_box(ax1.bbox)
            ellSigma1.set_alpha(1)
            # ellSigma1.set_facecolor('white')
            ellSigma1.set_linewidth(2)            

            ellSigma1.set_edgecolor(ellColor)
            ax1.add_artist(ellSigma1)
            
            #ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=3*rmsToFWHM*RMS_vshift, height=3*rmsToFWHM*RMS_sigmav, angle=theta)     
            #ellSigma1.set_clip_box(ax1.bbox)
            #ellSigma1.set_alpha(0.2)
            #ellSigma1.set_facecolor(ellColor)
            #ax1.add_artist(ellSigma1)


        # Set axis limits
        ax1.set_xlim(xMin, xMax)
        ax1.set_ylim(yMin, yMax)
        ax1.legend = plt.legend(frameon=True,handlelength=0, handletextpad=1,loc=3,borderaxespad=0.,markerscale=2,prop={'size': 18})
         
        for lh in ax1.legend.legendHandles: 
            lh.set_alpha(1)        
            #lh.set_sizes([50.0])
            #lh._legmarker.set_markersize(18)
        ax1.legend.get_frame().set_edgecolor('black')
        ax1.legend.get_frame().set_facecolor('white')
        # for legend_handle in  ax1.legend.legendHandles:
        #     legend_handle._legmarker.set_markersize(9)

        outPlotDir = cfg_par['general']['plotMomModDir']

        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)

        ax1.set_autoscale_on(False)    

        outPlot = outPlotDir+'K-multiple'+cfg_par['multipleRegions']['Name']+'.png'
        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotTheoreticalCCA'] == True or cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotTheoreticalCCA'] == True:        
            outPlot = outPlotDir+'K-multiple'+cfg_par['multipleRegions']['Name']+'ell.png'
        print(outPlot)
        plt.savefig(outPlot,format=cfg_par['lineRatios']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=300)#,
                # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        #plt.show()
        #plt.close()
               
        return 0  

    def taylorNumberRadius(self,cfg_par):
        '''Draws the taylorNumber vs radius for multi-phase gas values of sigma_v and radius given in inTables.


        Parameters
        ----------

        cfg_par['CCAanalysis']: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.
            

        inTables: list, str
            List with tablenames, one for each phase of the gas.
        
        maxRadius: float
            X-limit of plot in kpc

        gasColors: list, str
            List of colors for the different tables

        vRot: np.array()
            Rotational velocity for each phase of the gas. Radial dependency changes size of array.

        
        Requirements:
        -------------

        inTables must have been binned based on radius, binnedTableMean and binnedTableMeanErr extensions are considered
        NGOOD column must exist in BinnedTableMeanErr extension.

        Notes:
        ------
        Taylor Number: T_at=v_rot/sigma_v


        '''

        tableNames = cfg_par['CCAanalysis']['taylorNumber']['inTables']
        colorNames = cfg_par['CCAanalysis']['taylorNumber']['gasColors']

        params = ut.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure(figsize =(10,8),constrained_layout=False)
        fig.set_tight_layout(False)
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])

        ax1.set_xscale('symlog')
        ax1.set_yscale('symlog')
        
        #ax1.set_xlabel(r'$r$\,\, [kpc]')
        ax1.set_ylabel(r'${\rm T}_{{\rm a}_t} \simeq \frac{v_{\rm rot}}{\sigma_{\rm los}}$')
        # Calculate axis limits and aspect ratio
        xMin = 0.
        xMax = cfg_par['CCAanalysis']['taylorNumber']['maxRadius']+cfg_par['CCAanalysis']['taylorNumber']['maxRadius']/100.*2. 
        yMin = -1.5

        ax1.annotate(r'${\rm T}_{{\rm a}_t}>3$: rotation dominated',
            xy=(0.05, np.log10(3.1)), xycoords='data')

        ax1.annotate(r'$1<{\rm T}_{{\rm a}_t}<3$: turbulence $\&$ rotation',
            xy=(0.05, np.log10(1.5)), xycoords='data')


        ax1.annotate(r'${\rm T}_{{\rm a}_t}<1$: turbulence dominated',
            xy=(0.05, np.log10(0.75)), xycoords='data')        

        #ax1.ticklabel_format(axis='both', style='plain')
        
        ax1.minorticks_on()
        #ax1.xaxis.set_major_locator(MultipleLocator(1))
        #ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax1.xaxis.set_major_formatter(mticker.ScalarFormatter())

        ax1.tick_params(axis='both', which='major', pad=5)
        ax1.tick_params(axis='both', which='minor', pad=2)      
        #ax1.xaxis.set_minor_formatter(mticker.ScalarFormatter())
        #ax1.yaxis.set_minor_formatter(mticker.ScalarFormatter())

        majors = [np.log10(1e-1), np.log10(1), np.log10(3), np.log10(1e1)]
        ax1.yaxis.set_major_locator(FixedLocator(majors))
        majors = [1e-1, 1, 3, 1e1]
        #minors = [7.5e-2, 0.25, 0.75]
        #ax1.yaxis.set_minor_locator(FixedLocator(minors))
        ax1.yaxis.set_ticklabels(["$%.2f$" % y for y in majors]) 

        majors = [0, 1, 2, 3,4,5]
        ax1.xaxis.set_major_locator(FixedLocator(majors))
        #ax1.yaxis.set_minor_locator(ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())(base=10,numticks=5))
        #ax1.yaxis.set_minor_locator(AutoMinorLocator())

        ax1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                     bottom=True, top=True, left=True, right=True)

        #ax1.xaxis.set_ticklabels([])

        yMax=[]
        #yMin=[]
        for i in range(0,len(tableNames)):
            yColName = cfg_par['CCAanalysis']['taylorNumber']['yColumns'][i]
            gasName = cfg_par['CCAanalysis']['taylorNumber']['gasNames'][i]
            colorName = cfg_par['CCAanalysis']['taylorNumber']['gasColors'][i]

            hdul = fits.open(cfg_par['general']['runNameDir']+tableNames[i])
            ancels=hdul['BinnedTableMean'].data
            ancelsErr=hdul['BinnedTableMeanErr'].data

            x = ancels['r']
            Tat =np.divide(cfg_par['CCAanalysis']['taylorNumber']['vRot'],np.sqrt(3.)*ancels[yColName])
            print(Tat)
            y = np.log10(Tat)

            yErr = np.divide(cfg_par['CCAanalysis']['taylorNumber']['vRot'],np.sqrt(3.)*ancelsErr[yColName])
            yErrGood = ancelsErr['NGOOD']
            yErr= np.divide(yErr,np.sqrt(yErrGood))
            yErr = np.divide(yErr,Tat)

            index = np.where(x<=cfg_par['CCAanalysis']['taylorNumber']['maxRadius'])

            xPlot = x[index]
            yPlot = y[index]
            

            yErrPlot = yErr[index]/2.
            
            xErrPlot=np.diff(xPlot)/2.
            xErrPlot=np.append((xPlot[1]-xErrPlot[1]-xMin)/2.,xErrPlot)  
            # initialize figure
            #print((idxAGN),(idxKew),(idxKauf),(idxBad))
            ax1.errorbar(xPlot, yPlot, yerr=yErrPlot, xerr=xErrPlot, c=colorName,
                marker='o',markersize=8, alpha=0.8,label=gasName,ls='',capsize=3,lw=1.5)
            yMax.append(np.nanmax(yPlot+yErrPlot))
            #yMin.append(np.nanmin(yPlot-yErrPlot))

        yMax =np.nanmax(yMax)+ np.nanmax(yMax)/100.*2.
        #if np.nanmin(yMin)<0:
        #     yMin =np.nanmin(yMin)- np.nanmin(yMin)/100.*10.
        # else:
        #     yMin=0.

        #ax1.fill_between([xMin,xMax],y1=np.log10(1),y2=np.log10(3),facecolor='purple', alpha=0.2)
        ax1.axhline(y=np.log10(3),xmin=xMin,xmax=xMax,color='black', alpha=1,lw=2,ls=':')
        ax1.axhline(y=np.log10(1),xmin=xMin,xmax=xMax,color='black', alpha=1,lw=2,ls=':')        
        xhoriz = np.logspace(np.log10(1),np.log10(xMax+1),8,base=10.)
        #xhoriz=np.log10(xhoriz)
        yhoriz=np.zeros(8)+np.log10(1e-1)
        
        #upperlimits = np.array([1, 0]*5)
        #lowerlimits = np.array([0, 1]*5)
        ax1.errorbar(xhoriz, yhoriz,marker=r'$\downarrow$',markersize=30,color='red', lw=0,alpha=0.5,
            label='[NII]6583, HI, CO out of disk')
        
        # Set axis limits
        ax1.set_xlim(xMin, xMax)
        ax1.set_ylim(yMin, yMax)

        # get handles
        handles, labels = ax1.get_legend_handles_labels()
        # remove the errorbars
        handles = [h[0] for h in handles]
        # use them in the legend
        ax1.legend=plt.legend(handles, labels, loc='lower left',numpoints=1,handlelength=0, handletextpad=1,borderaxespad=0,markerscale=1,prop={'size': 18})
         
        for lh in ax1.legend.legendHandles: 
            lh.set_alpha(1)        
            #lh.set_sizes([50.0])
            #lh._legmarker.set_markersize(18)
        
        ax1.legend.get_frame().set_edgecolor('black')
        ax1.legend.get_frame().set_facecolor('white')
        # for legend_handle in  ax1.legend.legendHandles:
        #     legend_handle._legmarker.set_markersize(9)



        ccaDir =cfg_par['general']['runNameDir']+'CCA/'
        if not os.path.exists(ccaDir):
            os.mkdir(ccaDir)
        outPlotDir =cfg_par['general']['runNameDir']+'CCA/plots/'

        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)

        #ax1.set_autoscale_on(False)    

        outPlot = outPlotDir+'T_at-multiple.png'

        plt.savefig(outPlot,format=cfg_par['lineRatios']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=300)#

    def cRatioPlot(self,cfg_par):
        '''Plots the t_cool/t_eddy vs radius for different tables binned in radius.

        Parameters
        ----------

        cfg_par['CCAanalysis']['cRatio']: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.
        
        Requirements:
        -------------

        tables must contain a column with velocity dispersion, one with radius and one with t_cool.
        Columnnames are given in the parameter file

        inTables must have been binned based on radius, binnedTableMean and binnedTableMeanErr extensions are considered
        NGOOD column must exist in BinnedTableMeanErr extension.

        Notes:
        ------
        t_eddy: t_eddy=2pi r/sigma_v(r)
        t_cool: t_cool=3/2(ne+ni)kbT/(ne*ni*Lambda(T))

        '''
        
        tableNames = cfg_par['CCAanalysis']['cRatio']['inTables']
        colorNames = cfg_par['CCAanalysis']['cRatio']['gasColors']

        params = ut.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure(figsize =(10,8),constrained_layout=False)
        fig.set_tight_layout(False)
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])

        ax1.set_xscale('symlog')
        ax1.set_yscale('symlog')

        ax1.set_xlabel(r'$r$\,\, [kpc]')
        ax1.set_ylabel(r'$C \simeq \frac{t_{\rm cool}}{t_{\rm eddy}}$')
    
        # Calculate axis limits and aspect ratio
        xMin = 1
        xMax = cfg_par['CCAanalysis']['cRatio']['maxRadius']+cfg_par['CCAanalysis']['cRatio']['maxRadius']/100.*2. 
        yMin = -1.
        
        ax1.minorticks_on()
        #ax1.ticklabel_format(axis='both', style='plain')

        #ax1.xaxis.set_major_locator(MultipleLocator(1))
        #ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
        #ax1.yaxis.set_major_locator(MultipleLocator(3))
        #ax1.yaxis.set_minor_locator(MultipleLocator(1))

        ax1.xaxis.set_major_formatter(mticker.ScalarFormatter())

        ax1.tick_params(axis='both', which='major', pad=5)
        ax1.tick_params(axis='both', which='minor', pad=2)

        majors = [np.log(0.5),np.log10(1e-1), np.log10(1), np.log10(3), np.log10(1e1),np.log10(2e1)]
        ax1.yaxis.set_major_locator(FixedLocator(majors))
        majors = [0.5,1e-1, 1, 3, 1e1,2e1]
        #minors = [7.5e-2, 0.25, 0.75]
        #ax1.yaxis.set_minor_locator(FixedLocator(minors))
        ax1.yaxis.set_ticklabels(["$%.2f$" % y for y in majors]) 
        
        #minors = [7.5e-2, 0.25, 0.75,1.5,]
        
        #ax1.yaxis.set_minor_locator(FixedLocator(minors))
        
        #ax1.yaxis.set_ticklabels(["$%.2f$" % y for y in majors]) 

        majors = [1, 2, 3,4,5]
        ax1.xaxis.set_major_locator(FixedLocator(majors))
        #ax1.xaxis.set_major_formatter(mticker.ScalarFormatter())
        #ax1.yaxis.set_minor_formatter(mticker.ScalarFormatter())
        ax1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                     bottom=True, top=True, left=True, right=True)


        yMax=[]
        #yMin=[]
        #ax1.set_yscale('log',nonposy='clip')

        for i in range(0,len(tableNames)):
            tEddyColumns = cfg_par['CCAanalysis']['cRatio']['tEddy'][i]
            tEddyErrColumnsUp = cfg_par['CCAanalysis']['cRatio']['tEddy'][i]+'Up'
            tEddyErrColumnsDown = cfg_par['CCAanalysis']['cRatio']['tEddy'][i]+'Down'

            rColumns = cfg_par['CCAanalysis']['cRatio']['rColumns'][i]
            tCoolColumns = cfg_par['CCAanalysis']['cRatio']['tCool'][i]

            hdul = fits.open(cfg_par['general']['runNameDir']+tableNames[i])
            ancels=hdul['BinnedTableMean'].data
            ancelsErr=hdul['BinnedTableMeanErr'].data


            Cratio=ancels[tCoolColumns]/ancels[tEddyColumns]
            CratioUp = ancels[tCoolColumns]/ancelsErr[tEddyErrColumnsDown]
            CratioDown = ancels[tCoolColumns]/ancelsErr[tEddyErrColumnsUp]
            print(CratioUp,CratioDown)
            CratioErr=((CratioUp-CratioDown))/2.

            gasName = cfg_par['CCAanalysis']['cRatio']['gasNames'][i]
            colorName = cfg_par['CCAanalysis']['cRatio']['gasColors'][i]

            hdul = fits.open(cfg_par['general']['runNameDir']+tableNames[i])
            ancels=hdul['BinnedTableMean'].data
            ancelsErr=hdul['BinnedTableMeanErr'].data

            x = ancels['r']
            y = Cratio
            print(gasName)
            print('Cratio= ',y)

            #yErr = ancelsErr[tEddyColumns]

            #yErr = CratioErr*yErr
            yErr = CratioErr


            index = np.where(x<=cfg_par['CCAanalysis']['cRatio']['maxRadius'])

            xPlot = x[index]
            yPlot = y[index]
            yErrPlot = yErr[index]/2.
            
            xErrPlot=np.diff(xPlot)/2.
            xErrPlot=np.append((xPlot[1]-xErrPlot[1]-xMin)/2.,xErrPlot)  
            # initialize figure
            #print((idxAGN),(idxKew),(idxKauf),(idxBad))
            ax1.errorbar(xPlot, np.log10(yPlot), yerr=yErrPlot, xerr=xErrPlot, c=colorName,
                marker='o',markersize=8, alpha=0.8,label=gasName,ls='',capsize=3,lw=1.5)
            yMax.append(np.nanmax(np.log10(yPlot)+yErrPlot))
            #yMin.append(np.nanmin(yPlot-yErrPlot))

        yMax =np.nanmax(yMax)+ np.nanmax(yMax)/100.*0.5
        #yMax=np.log10(3e1)
        #if np.nanmin(yMin)<0:
            #yMin =np.nanmin(yMin)- np.nanmin(yMin)/100.*10.
        #    yMin=-1.            
        #else:
        #yMin=yMin

        ax1.fill_between([xMin,xMax],y1=np.log10(0.3),y2=np.log10(2.6),facecolor='purple', alpha=0.2)
        #xhoriz = np.linspace(xMin-1,xMax+1,8)
        #xhoriz=np.log10(xhoriz)
        #yhoriz=np.zeros(8)+np.log10(3)
        # Set axis limits
        #ax1.errorbar(xhoriz, yhoriz,marker=r'$\downarrow$',markersize=30,color='purple', ls='--',alpha=0.5)


        ax1.set_xlim(xMin, xMax)
        ax1.set_ylim(yMin, yMax)

        # get handles
        handles, labels = ax1.get_legend_handles_labels()
        # remove the errorbars
        handles = [h[0] for h in handles]
        # use them in the legend
        ax1.legend=plt.legend(handles, labels, loc='lower left',numpoints=1,
            handlelength=0, handletextpad=1,borderaxespad=0,markerscale=1,prop={'size': 18})
         
        for lh in ax1.legend.legendHandles: 
            lh.set_alpha(1)        
            #lh.set_sizes([50.0])
            #lh._legmarker.set_markersize(18)
        
        ax1.legend.get_frame().set_edgecolor('black')
        ax1.legend.get_frame().set_facecolor('white')

        # for legend_handle in  ax1.legend.legendHandles:
        #     legend_handle._legmarker.set_markersize(9)

        ccaDir =cfg_par['general']['runNameDir']+'CCA/'
        if not os.path.exists(ccaDir):
            os.mkdir(ccaDir)
        outPlotDir =cfg_par['general']['runNameDir']+'CCA/plots/'

        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)

        ax1.set_autoscale_on(False)    


        outPlot = outPlotDir+'Cratio-multiple.png'

        plt.savefig(outPlot,format=cfg_par['lineRatios']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=300)#




    def inCCARegionTable(self,cfg_par):
        '''Computes for each fitted line the difference with the expected values of sigma and centroid predicted by
        Cold Chaotic Accretion.

        Saves in .fits files the following maps:
            - rotating gas (if XXX has been previously run)
            - gas within 3 sigma from CCA

        Calls `momPlay()` to plot all maps. 

        Variables are set in section `kinematicalAnalysis` (for the ionised gas) and in `otherGasKinAnalysis`
        (for the cold gas)

        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        ''' 
        lineInfo = tP.openLineList(cfg_par)
        lineThresh = float(lineInfo['SNThresh'][0])


        f = fits.open(cfg_par['general']['dataCubeName'])
        momHead = f[0].header
        f.close()
        print(cfg_par['general']['outTableName'])

        hdul = fits.open(cfg_par['general']['outTableName'])
 
        if cfg_par['gFit']['modName'] == 'BF' and cfg_par['otherGasKinAnalysis']['enable'] == False:
            modName = 'g2'
            residuals = hdul['Residuals_'+modName].data
            anc = hdul['ancels'+cfg_par['gFit']['modName']].data
        elif cfg_par['gFit']['modName'] == 'BF' and cfg_par['otherGasKinAnalysis']['enable'] == True:
            anc = hdul['ancelsBF'].data

        bins = hdul['BININFO'].data
 
        BinAnc = anc['BIN_ID']
        BinBins = bins['BIN_ID']
        pixXSrt = bins['PixX']
        pixYSrt = bins['PixY']

        
        if not cfg_par['otherGasKinAnalysis']['enable']== True:
            SN = residuals['SN_NII6583']
            
            linesG1 = hdul['LineRes_G1'].data
            
            if cfg_par['gFit']['modName'] == 'BF':
                modName='g2'
            
            residuals = hdul['Residuals_'+modName].data
            x=anc['logCentroid_'+cfg_par['kinematicalAnalysis']['Name']]
            y=anc['logDispIntr_'+cfg_par['kinematicalAnalysis']['Name']]

            sigmaCCA = cfg_par['kinematicalAnalysis']['ancillaryInfo']['sigmaInCCA']

        else:
            x=anc['logCentroid_'+cfg_par['otherGasKinAnalysis']['Name']]
            y=anc['logSigma_'+cfg_par['otherGasKinAnalysis']['Name']]
            sigmaCCA = cfg_par['otherGasKinAnalysis']['ancillaryInfo']['sigmaInCCA']

        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['theoreticalCCA'] == 'Ensemble':
            Mean_sigmav = 2.13 
            RMS_sigmav  = 0.13 
            Mean_vshift = 1.59 
            RMS_vshift  = 0.37
            theta = -179.88 #covariance angle (contours inclined)
            ellColor='purple'
        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['theoreticalCCA'] == 'Pencil':
            Mean_sigmav = 1.65 
            RMS_sigmav  = 0.41
            Mean_vshift = 2.0
            RMS_vshift  = 0.47
            theta = 165.16 #covariance angle (contours inclined)
            ellColor='darkseagreen'

        rmsToFWHM = 2.*np.sqrt(2.*np.log(2))
        ellWidth = rmsToFWHM*RMS_vshift*sigmaCCA
        ellHeight = rmsToFWHM*RMS_sigmav*sigmaCCA
        cos_angle = np.cos(np.radians(180.-theta))
        sin_angle = np.sin(np.radians(180.-theta))

        xc = x - Mean_vshift
        yc = y - Mean_sigmav

        xct = xc * cos_angle - yc * sin_angle
        yct = xc * sin_angle + yc * cos_angle 

        rad_cc = (xct**2/(ellWidth/2.)**2) + (yct**2/(ellHeight/2.)**2)

        CCAvec = np.empty(len(anc['BIN_ID']))*np.nan
        CCAMap = np.zeros([momHead['NAXIS2'],momHead['NAXIS1']])*np.nan
        CCAOnly = np.zeros([momHead['NAXIS2'],momHead['NAXIS1']])*np.nan
        CCAund = np.zeros([momHead['NAXIS2'],momHead['NAXIS1']])*np.nan
        
        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotRotation'] == True or cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotRotation']==True:
            
            rotModCol=anc['RotMod']
            RotMap = np.zeros([momHead['NAXIS2'],momHead['NAXIS1']])*np.nan
            RotMapOnly = np.zeros([momHead['NAXIS2'],momHead['NAXIS1']])*np.nan

            RotMapCCA = np.zeros([momHead['NAXIS2'],momHead['NAXIS1']])*np.nan

        for i in range(0,len(rad_cc)):
            
            match_bin = np.where(BinBins==BinAnc[i])[0]

            for index in match_bin:

                if not cfg_par['otherGasKinAnalysis']['enable']== True:

                    thresHold = SN[index]
                    #if thresHold >= lineThresh and sigmaThresh < cfg_par['moments']['sigmaThresh']:
                    if thresHold >= lineThresh:
                    
                        if rad_cc[i] <= 1.:
                            # point in ellipse
                            CCAvec[i] = rad_cc[i]
                            CCAOnly[int(pixYSrt[index]),int(pixXSrt[index])] = rad_cc[i]
                        else:
                            # point not in ellipse
                            CCAvec[i] = -1.

                        CCAMap[int(pixYSrt[index]),int(pixXSrt[index])] = CCAvec[i]
                    
                        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotRotation'] == True:
                            if rotModCol[i]==1.:
                                RotMap[int(pixYSrt[index]),int(pixXSrt[index])] = 1.

                            else:
                                RotMap[int(pixYSrt[index]),int(pixXSrt[index])] = -1.
                    
                        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['rotOutCCA'] == True:

                            if rotModCol[i] ==1.:
                                RotMapCCA[int(pixYSrt[index]),int(pixXSrt[index])] = 2.
                                RotMapOnly[int(pixYSrt[index]),int(pixXSrt[index])] = 1.                            
                            
                            elif rotModCol[i] !=1. and CCAvec[i]!=-1.:
                            
                                RotMapCCA[int(pixYSrt[index]),int(pixXSrt[index])] = CCAvec[i]

                            else:
                                RotMapCCA[int(pixYSrt[index]),int(pixXSrt[index])] = -1

                                CCAund[int(pixYSrt[index]),int(pixXSrt[index])] = 1.

                else:

                    if not np.isnan(x[i]) and not np.isnan(y[i]):
                        if rad_cc[i] <= 1.:
                            # point in ellipse
                            CCAvec[i] = rad_cc[i]
                            CCAOnly[int(pixYSrt[index]),int(pixXSrt[index])] = rad_cc[i]
                        else:
                            # point not in ellipse
                            CCAvec[i] = -1.

                        CCAMap[int(pixYSrt[index]),int(pixXSrt[index])] = CCAvec[i]
                    
                        if cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotRotation'] == True:
                            if rotModCol[i]==1.:
                                RotMap[int(pixYSrt[index]),int(pixXSrt[index])] = 1.

                            else:
                                RotMap[int(pixYSrt[index]),int(pixXSrt[index])] = -1.
                    
                        if cfg_par['otherGasKinAnalysis']['ancillaryInfo']['rotOutCCA'] == True:

                            if rotModCol[i] ==1.:
                                RotMapCCA[int(pixYSrt[index]),int(pixXSrt[index])] = 2.
                                RotMapOnly[int(pixYSrt[index]),int(pixXSrt[index])] = 1.                            
                            
                            elif rotModCol[i] !=1. and CCAvec[i]!=-1.:
                            
                                RotMapCCA[int(pixYSrt[index]),int(pixXSrt[index])] = CCAvec[i]

                            else:
                                RotMapCCA[int(pixYSrt[index]),int(pixXSrt[index])] = -1

                                CCAund[int(pixYSrt[index]),int(pixXSrt[index])] = 1.                  

        if 'CUNIT3' in momHead:
            del momHead['CUNIT3']
        if 'CTYPE3' in momHead:
            del momHead['CTYPE3']
        if 'CDELT3' in momHead:
            del momHead['CDELT3']
        if 'CRVAL3' in momHead:  
            del momHead['CRVAL3']
        if 'CRPIX3' in momHead:
            del momHead['CRPIX3'] 
        if 'NAXIS3' in momHead:
            del momHead['NAXIS3']
        if 'CRDER3' in momHead:
            del momHead['CRDER3']

        Head = momHead.copy()
        Head['WCSAXES'] = 2
        Head['SPECSYS'] = 'topocent'
        Head['BUNIT'] = 'Jy'

        if cfg_par['kinematicalAnalysis']['enable']== True:
            nameGas =  cfg_par['kinematicalAnalysis']['Name']
            Labels= [nameGas+' unidentified',
            nameGas+' CCA within '+str(cfg_par['kinematicalAnalysis']['ancillaryInfo']['sigmaInCCA'])+r'$\sigma$']
        elif cfg_par['otherGasKinAnalysis']['enable']== True:
            nameGas =  cfg_par['otherGasKinAnalysis']['Name']            
            Labels= [nameGas+' unidentified',
            nameGas+' CCA within '+str(cfg_par['otherGasKinAnalysis']['ancillaryInfo']['sigmaInCCA'])+r'$\sigma$']

        CCAOnlyName = cfg_par['general']['momModDir']+'ccaMapOnly-'+str(sigmaCCA)+'-'+nameGas+'.fits'
        CCAMapName  = cfg_par['general']['momModDir']+'ccaMap-'+str(sigmaCCA)+'-'+nameGas+'.fits'
        CCAundName  = cfg_par['general']['momModDir']+'ccaMapUnd-'+str(sigmaCCA)+'-'+nameGas+'.fits'

        fits.writeto(CCAundName,CCAund,Head,overwrite=True)
        fits.writeto(CCAOnlyName,CCAOnly,Head,overwrite=True)
        fits.writeto(CCAMapName,CCAMap,Head,overwrite=True)
        

        colorList= ['gray','seagreen']
        CustomCmap = ListedColormap(colorList)
        

        mPl.momAncPlot(cfg_par, CCAMapName, Labels,CustomCmap,colorList)

        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotRotation'] == True or cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotRotation'] == True:
            RoTName= cfg_par['general']['momModDir']+'ccaMapRot-'+str(sigmaCCA)+'-'+nameGas+'.fits'

            fits.writeto(RoTName,RotMap,Head,overwrite=True)
            Labels= [nameGas,
            nameGas+' in rotation']
            colorList = ['seagreen','blue']
            CustomCmap = ListedColormap(colorList)

            mPl.momAncPlot(cfg_par, RoTName,Labels,CustomCmap,colorList)

        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['rotOutCCA'] == True or cfg_par['otherGasKinAnalysis']['ancillaryInfo']['rotOutCCA'] == True:
            CCARotName = cfg_par['general']['momModDir']+'ccaMap-'+nameGas+'ccaMapCCARot-'+str(sigmaCCA)+'-'+nameGas+'.fits'
            fits.writeto(CCARotName,RotMapCCA,Head,overwrite=True)
            Labels= [nameGas+' unidentified',
            nameGas+' CCA within '+str(cfg_par['kinematicalAnalysis']['ancillaryInfo']['sigmaInCCA'])+r'$\sigma$',nameGas+' in rotation']
            colorList = ['gray','seagreen','blue']
            CustomCmap = ListedColormap(colorList)           

            RoTNameOnly= cfg_par['general']['momModDir']+'ccaMapRotOnly-'+str(sigmaCCA)+'-'+nameGas+'.fits'
            fits.writeto(RoTNameOnly,RotMapOnly,Head,overwrite=True)
            
            mPl.momAncPlot(cfg_par, CCARotName,Labels,CustomCmap,colorList)
            
            CustomCmapOne = ListedColormap('gray')           
            CustomCmapThree = ListedColormap('blue')           

            mPl.momAncPlotOver(cfg_par, CCAundName, CCAOnlyName, RoTNameOnly, Labels,CustomCmapOne,CustomCmapThree,colorList)


        


        t=Table(anc)
 
        if 'CCAIN' not in anc.dtype.names: 
            t.add_column(Column(CCAvec,name='CCAIN'))
        else:
            t.replace_column('CCAIN',Column(CCAvec,name='CCAIN'))        
        
        hdul['Ancels'+cfg_par['gFit']['modName']] = fits.BinTableHDU(t.as_array(),name='Ancels'+cfg_par['gFit']['modName'])

        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)

        return 0
