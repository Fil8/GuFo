#!/usr/bin/env python3.6

import os, sys
from astropy.io import fits
from astropy.table import Table, Column

import numpy as np

import pandas as pd

from matplotlib.patches import Ellipse

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

import cvPlay, tPlay, momPlot

cvP = cvPlay.convert()
tP = tPlay.tplay()
mPl = momPlot.MOMplot()

class ancelsplot(object):

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
          #'text.latex.unicode'  : True
           }
        
        return params

    def sigmaCentroid(self,cfg_par,xyRange=None,outPlotDir=None):

        lineInfo = tP.openLineList(cfg_par)

        lineThresh = float(lineInfo['SNThresh'][0])


        hdul = fits.open(cfg_par['general']['outTableName'])
        ancels = hdul['Ancels'+cfg_par['gFit']['modName']].data
        print(ancels.dtype.names)
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
            colorSca=ancels['CCAIN'][binCode[1]]
            rotSca = ancels['RotMod'][binCode[1]]
        elif cfg_par['gFit']['modName'] == 'g3':
            modString = ['G1','G2','G3','ToT']
            x = [np.log10(lines['g1_Centre_NII6583']),np.log10(lines['g2_Centre_NII6583']),np.log10(lines['g3_Centre_NII6583']),
            ancels['loCentroid_NII6583']]

            binCode = [lines['BIN_ID'],lines['BIN_ID'],lines['BIN_ID'],ancels['BIN_ID']]
        
            y = [np.log10(lines['g1_sigIntr_NII6583']),np.log10(lines['g2_sigIntr_NII6583']),np.log10(lines['g3_sigIntr_NII6583']),
            ancels['logDispIntr_NII6583']]
        
        cfg_par['gFit']['modName'] = 'BF'
        for i in range (0, len(modString)):
            
            # initialize figure
            params = self.loadRcParams()
            plt.rcParams.update(params)
            fig = plt.figure(figsize =(10,8))
            fig.subplots_adjust(hspace=0.0)
            gs = gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0])

        
            #ax.set_xticks([])
            
            ax1.set_xlabel(r'log($|\overline{v}_{\rm los}|$)\,\, [km s$^{-1}$]')
            ax1.set_ylabel(r'log($\sigma_{\rm los}$)\,\, [km s$^{-1}$]')
        
            # Calculate axis limits and aspect ratio
            xMin = 0.
            xMax = 2.9
            yMin = 0.52
            yMax = 3.2
            
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
                ellColor='darkseagreen'

            sigmaThresh = lines['g1_SigIntr_NII6583']
            
            idx = ((SNSort<lineThresh))
            x[i][idx] = np.nan
            y[i][idx] = np.nan

            rmsToFWHM = 2.*np.sqrt(2.*np.log(2))
            ellSigma1 = Ellipse(xy=(Mean_vshift,Mean_sigmav), width=rmsToFWHM*RMS_vshift, height=rmsToFWHM*RMS_sigmav, angle=theta,
                label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['theoreticalCCA']+' beam')     
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

            #print((idxAGN),(idxKew),(idxKauf),(idxBad))
            if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotRotation'] == True:
                indexRot = rotSca==1.
                indexElse = rotSca==0.
                #print(indexRot)
                ax1.scatter(x[i][indexElse], y[i][indexElse], facecolors='red', marker='.', s=20, linewidths=None,edgecolors='red', alpha=0.5,
                    label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel'])
                ax1.scatter(x[i][indexRot], y[i][indexRot], facecolors=None, marker='.', s=20, linewidths=None,edgecolors='blue', alpha=0.5,
                    label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel']+' in rotation')

                #
            else:    
                ax1.scatter(x[i], y[i], c='red', marker='.', s=20, linewidths=None,edgecolors='red', 
                label=cfg_par['kinematicalAnalysis']['ancillaryInfo']['CCALabel'])

            #idxAGN = np.where(colorSca==1)
            #ax1.scatter(x[i][idxAGN], y[i][idxAGN], c='black', marker='.', s=20, linewidths=None,edgecolors=None,alpha=0.1)

            # Set axis limits
            ax1.set_xlim(xMin, xMax)
            ax1.set_ylim(yMin, yMax)
            ax1.legend = plt.legend(loc=1, prop={'size': 12})
            ax1.legend.get_frame().set_edgecolor('black')
            ax1.legend.get_frame().set_facecolor('white')

            if outPlotDir==None:
                outPlotDir = cfg_par['general']['plotMomModDir']

            if not os.path.exists(outPlotDir):
                os.mkdir(outPlotDir)

            outPlot = outPlotDir+'sigmaCentroid-'+modString[i]+'-NII6583.png'
           
            plt.savefig(outPlot,format=cfg_par['lineRatios']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
                    # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
            plt.show()
            plt.close()
               
        return 0

    def sigmaCentroidColdGas(self,cfg_par,xyRange=None):

        lineInfo = tP.openLineList(cfg_par)

        lineThresh = float(lineInfo['SNThresh'][0])


        hdul = fits.open(cfg_par['general']['outTableName'])
        ancels = hdul['AncelsG1'].data

        x = ancels['logCentroid_'+cfg_par['otherGasKinAnalysis']['Name']]
        y = ancels['logSigma_'+cfg_par['otherGasKinAnalysis']['Name']]

        print(x,y)
        # initialize figure
        params = self.loadRcParams()
        plt.rcParams.update(params)
        fig = plt.figure(figsize =(10,8))
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0])

        #ax.set_xticks([])
        
        ax1.set_xlabel(r'log($|\overline{v}_{\rm los}|$)\,\, [km s$^{-1}$]')
        ax1.set_ylabel(r'log($\sigma_{\rm los}$)\,\, [km s$^{-1}$]')
    
        # Calculate axis limits and aspect ratio
        xMin = 0.
        xMax = 2.9
        yMin = 0.52
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
        if cfg_par['ancillary']['plotRotation'] == True:
            indexRot = np.where(ancels['RotMod']==1.)[0]
            indexElse = np.where(ancels['RotMod']==0.)[0]

            ax1.scatter(x[indexRot], y[indexRot], c='blue', marker='.', s=20, linewidths=None,edgecolors='red', 
                label=cfg_par['otherGasKinAnalysis']['ancillaryInfo']['CCALabel']+'in rotation')

            ax1.scatter(x[indexElse], y[indexElse], c='blue', marker='.', s=20, linewidths=None,edgecolors='red', 
                label=cfg_par['otherGasKinAnalysis']['ancillaryInfo']['CCALabel'])
        else:    
            ax1.scatter(x, y, c='red', marker='.', s=20, linewidths=None,edgecolors='red', 
            label=cfg_par['otherGasKinAnalysis']['ancillaryInfo']['CCALabel'])

        # Set axis limits
        ax1.set_xlim(xMin, xMax)
        ax1.set_ylim(yMin, yMax)
        ax1.legend = plt.legend(loc=1, prop={'size': 12})
        ax1.legend.get_frame().set_edgecolor('black')
        ax1.legend.get_frame().set_facecolor('white')


        outPlotDir = cfg_par['general']['plotMomModDir']

        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)

        outPlot = outPlotDir+'sigmaCentroid-g1-'+cfg_par['otherGasKinAnalysis']['Name']+'.png'
        print(outPlot)
        plt.savefig(outPlot,format=cfg_par['lineRatios']['plotFormat'], bbox_inches = "tight",overwrite=True,dpi=100)#,
                # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        #plt.show()
        #plt.close()
               
        return 0     

    def inCCARegionTable(self,cfg_par):
 
        lineInfo = tP.openLineList(cfg_par)
        lineThresh = float(lineInfo['SNThresh'][0])


        f = fits.open(cfg_par['general']['dataCubeName'])
        momHead = f[0].header
        f.close()

        hdul = fits.open(cfg_par['general']['outTableName'])

        anc = hdul['ancels'+cfg_par['gFit']['modName']].data
        bins = hdul['BININFO'].data
        
        if not cfg_par['otherGasKinAnalysis']['enable']== True:
            
            linesG1 = hdul['LineRes_G1'].data
            if cfg_par['gFit']['modName'] == 'BF':
                cfg_par['gFit']['modName'] = 'g2'
            residuals = hdul['Residuals_'+cfg_par['gFit']['modName']].data
            cfg_par['gFit']['modName'] = 'BF'
            x=anc['logCentroid_'+cfg_par['kinematicalAnalysis']['Name']]
            y=anc['logDispIntr_'+cfg_par['kinematicalAnalysis']['Name']]
            
        else:
            x=anc['logCentroid_'+cfg_par['otherGasKinAnalysis']['coldGas']['Name']]
            y=anc['logSigma_'+cfg_par['otherGasKinAnalysis']['coldGas']['Name']]

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
        ellWidth = rmsToFWHM*RMS_vshift
        ellHeight = rmsToFWHM*RMS_sigmav
        cos_angle = np.cos(np.radians(180.-theta))
        sin_angle = np.sin(np.radians(180.-theta))

        xc = x - Mean_vshift
        yc = y - Mean_sigmav

        xct = xc * cos_angle - yc * sin_angle
        yct = xc * sin_angle + yc * cos_angle 

        rad_cc = (xct**2/(ellWidth/2.)**2) + (yct**2/(ellHeight/2.)**2)

        CCAvec = np.empty(len(anc['BIN_ID']))*np.nan
        CCAMap = np.zeros([momHead['NAXIS2'],momHead['NAXIS1']])*np.nan
        
        if cfg_par['kinematicalAnalysis']['ancillaryInfo']['plotRotation'] == True or cfg_par['otherGasKinAnalysis']['ancillaryInfo']['plotRotation']:
            rotModCol=anc['RotMod']

            RotMap = np.zeros([momHead['NAXIS2'],momHead['NAXIS1']])*np.nan
            RotMapCCA = np.zeros([momHead['NAXIS2'],momHead['NAXIS1']])*np.nan

        for i in range(0,len(rad_cc)):
            
            match_bin = np.where(bins['BIN_ID']==anc['BIN_ID'][i])[0]

            for index in match_bin:

                if not cfg_par['otherGasKinAnalysis']['enable']== True:

                    thresHold = residuals['SN_NII6583'][index]
                    sigmaThresh = linesG1['g1_SigIntr_NII6583'][index]

                    #if thresHold >= lineThresh and sigmaThresh < cfg_par['moments']['sigmaThresh']:
                    if thresHold >= lineThresh:
                    
                        if rad_cc[i] <= cfg_par['kinematicalAnalysis']['ancillaryInfo']['sigmaInCCA']:
                            # point in ellipse
                            CCAvec[i] = rad_cc
                        else:
                            # point not in ellipse
                            CCAvec[i] = 0 

                    CCAMap[int(bins['PixY'][index]),int(bins['PixX'][index])] = CCAvec[i]
                    RotMap[int(bins['PixY'][index]),int(bins['PixX'][index])] = rotModCol[index]
                else:

                    if not np.isnan(x[i]) and not np.isnan(y[i]):
                        if rad_cc[i] <= cfg_par['otherGasKinAnalysis']['ancillaryInfo']['sigmaInCCA']:
                            # point in ellipse
                            CCAvec[i] = rad_cc 
                    
                        else:
                            # point not in ellipse
                            CCAvec[i] = 0 

                    CCAMap[int(bins['PixY'][index]),int(bins['PixX'][index])] = CCAvec[i]
                    RotMap[int(bins['PixY'][index]),int(bins['PixX'][index])] = rotModCol[index]


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

        if not cfg_par['otherGasKinAnalysis']['enable']== True:
            nameGas =  cfg_par['kinematicalAnalysis']['Name']
        else:
            nameGas =  cfg_par['otherGasKinAnalysis']['Name']            


        fits.writeto(cfg_par['general']['momModDir']+'ccaMap-'+nameGas+'.fits',CCAMap,Head,overwrite=True)
        if cfg_par['ancillary']['plotRotation'] == True:
            fits.writeto(cfg_par['general']['momModDir']+'ccaMap-'+nameGas+'rot.fits',RotMap,Head,overwrite=True)
            fits.writeto(cfg_par['general']['momModDir']+'ccaMap-'+nameGas+'rotCCA.fits',RotMapCCA,Head,overwrite=True)

        mPl.momAncPlot(cfg_par, cfg_par['general']['momModDir']+'ccaMap-'+nameGas+'.fits',
            nameGas,nameGas,nameGas,'ancillary')

        mPl.momAncPlot(cfg_par, cfg_par['general']['momModDir']+'ccaMap-'+nameGas+'rot.fits',
            nameGas,nameGas,nameGas,'ancillary')

        mPl.momAncPlot(cfg_par, cfg_par['general']['momModDir']+'ccaMap-'+nameGas+'rotCCA.fits',
            nameGas,nameGas,nameGas,'ancillary')


        t=Table(anc)
 
        if 'CCAIN' not in anc.dtype.names: 
            t.add_column(Column(CCAvec,name='CCAIN'))
        else:
            t.replace_column('CCAIN',Column(CCAvec,name='CCAIN'))        
        
        hdul['Ancels'+cfg_par['gFit']['modName']] = fits.BinTableHDU(t.as_array(),name='Ancels'+cfg_par['gFit']['modName'])

        hdul.writeto(cfg_par['general']['outTableName'],overwrite=True)

        return
