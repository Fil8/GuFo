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

class specplot:

#----------------------#
# rc param initialize
#----------------------#
    def loadRcParamsBig(self):
    
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

    def loadRcParamsZoom(self):
    
        params = {
          'font.family'         :' serif',
          'font.serif'          :'times',
          'font.style'          : 'normal',
          'font.weight'         : 'book',
          'font.size'           : 12,
          'axes.linewidth'      : 2,
          'lines.linewidth'     : 1.5,
          'xtick.labelsize'     : 9,
          'ytick.labelsize'     : 9, 
          'xtick.direction'     :'in',
          'ytick.direction'     :'in',
          'xtick.major.size'    : 4,
          'xtick.major.width'   : 1.5,
          'xtick.minor.size'    : 2,
          'xtick.minor.width'   : 0.75,
          'ytick.major.size'    : 4,
          'ytick.major.width'   : 1.5,
          'ytick.minor.size'    : 2,
          'ytick.minor.width'   : 0.75, 
          'text.usetex'         : True,
          'text.latex.unicode'  : True,
          'legend.fontsize'     : 10
           }
        
        return params


    def plotSpecFit(self,cfg_par,vel,y,result,noise,xx,yy,lineInfo,singleVorBinInfo):
        
        velPlot = np.exp(vel)
        yBFit = result.best_fit
        yRes = result.residual
        yInFit = result.init_fit
        key = 'general'
        

        outPlot = cfg_par['general']['outPlotDir']+str(xx)+'_'+str(yy)+'_'+cfg_par['gFit']['modName']+'.png'       
        

        #dely = result.eval_uncertainty(sigma=3)
            
         # initialize figure
        params = self.loadRcParamsBig()
        plt.rcParams.update(params)
        fig = plt.figure(figsize =(10,8))
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        plt.rc('xtick')

        # Initialize subplots
        ax1 = fig.add_subplot(gs[0])

        divider = make_axes_locatable(ax1)
        ax2 = divider.append_axes("bottom", size='15%', pad=0)
        ax1.figure.add_axes(ax2)
        

        #ax.set_xticks([])

        ax1.set_xlabel(r'Wavelength [$\AA$]')
        ax1.set_ylabel(r'Flux [-]')


        # Calculate axis limits and aspect ratio
        x_min = np.min(velPlot)
        x_max = np.max(velPlot)
        if cfg_par['gPlot']['fixed_scale']:
            y1_min = np.min(y)*1.2
            y1_max = np.max(y)*1.2
        else:
            y1_min = np.nanmin(y)*1.1
            y1_max = np.nanmax(y)*1.1

        # Set axis limits
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(y1_min, y1_max)

        ax1.tick_params(axis='both', which='major', pad=5)
        #ax1.xaxis.set_minor_locator()
        #ax1.yaxis.set_minor_locator()      
        
        ax1.step(velPlot, y, where='mid', color='black', linestyle='-')
        ax1.plot(velPlot, yBFit, 'r-', label='best fit')
        #ax1.step(vel, yInFit, 'b-', label='init fit')
        aicStr = str(int(result.aic))
        bicStr = str(int(result.bic))
        redchiStr = str(int(result.redchi))
        successStr = str(result.success)       
        xText = cfg_par['gFit']['lambdaMin']+50


        ax1.text(xText, y1_max*0.90, r'BIN ID:\t'+str(singleVorBinInfo['BIN_ID'][0]), {'color': 'k', 'fontsize': 20})
        ax1.text(xText, y1_max*0.94, r'X,Y:\t'+str(xx)+','+str(yy), {'color': 'k', 'fontsize': 20})

        #ax1.text(xText, x_max*0.85, r'Success:\t'+successStr, {'color': 'b'})
        ax1.text(xText, y1_max*0.88, r'$\tilde{\chi}^2$:\t'+redchiStr, {'color': 'k', 'fontsize': 20})
        ax1.text(xText, y1_max*0.82, r'aic:\t'+aicStr, {'color': 'k', 'fontsize': 20})
        ax1.text(xText, y1_max*0.76, r'bic:\t'+bicStr, {'color': 'k', 'fontsize': 20})

        #ax1.fill_between(vel, yBFit-dely, yBFit+dely, color="#ABABAB",
        #        label='3-$\sigma$ uncertainty band')
        if cfg_par['gFit']['modName'] !='g1':
            comps = result.eval_components()
            for i in xrange(0,len(lineInfo['ID'])):
                
                ax1.plot(velPlot, comps['g1ln'+str(i)+'_'], 'g--')
            
                if cfg_par['gFit']['modName'] =='g2':
                    ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
            
                elif cfg_par['gFit']['modName'] !='g2':
                    ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
                    ax1.plot(velPlot, comps['g3ln'+str(i)+'_'], 'c--')    


        # Calculate axis limits and aspect ratio
        x_min = np.min(velPlot)
        x_max = np.max(velPlot)
        if cfg_par['gPlot']['Res-fixed_scale']:
            y1_min = np.nanmin([-200.,np.nanmin(-noise)*1.1,np.nanmin(-yRes)*1.1])
            y1_max = np.nanmax([+200.,np.nanmax(+noise)*1.1,np.nanmax(+yRes)*1.1])
        else:
            y1_min = np.nanmin(yRes)*1.1
            y1_max = np.nanmax(yRes)*1.1  

        # Set axis limits
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y1_min, y1_max) 

        #ax2.plot(vel, amp(x, p(x,m,n)))
        ax2.step(velPlot, yRes, 'g-', label='residuals')
        ax2.axhline(color='k', linestyle=':', zorder=0)                           
        ax2.fill_between(velPlot, -noise, noise,
                         facecolor='grey', alpha=0.5,step='mid')

        #ax1.legend(loc='best')



        plt.savefig(outPlot,
                    format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        plt.show()
        plt.close()
           
        return 0

    def plotLineZoom(self,cfg_par,vel,y,result,noise,i,j,lineInfo,singleVorBinInfo):

        velExp = np.exp(vel)
        yBFit = result.best_fit
        yRes = result.residual
        print result.residual
        yInFit = result.init_fit
        key = 'general'
        
        outPlotDir = cfg_par['general']['outPlotDir']
        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)
        
        outPlot = outPlotDir+str(singleVorBinInfo['BIN_ID'][0])+'_'+cfg_par['gFit']['modName']+'.png'       

        params = self.loadRcParamsZoom()
        plt.rcParams.update(params)
        # add one row for the plot with full channel width
        n_plots = len(lineInfo['Wave'])
        n_rows = int(np.ceil(n_plots/3.))
        #fig, ax = plt.subplots(squeeze=False,
        #    ncols=3, nrows=n_rows, figsize=(8.25, 11.67))


        fig = plt.figure(figsize=(8.25, 11.67), constrained_layout=False)
        fig.subplots_adjust(hspace=0.)

        gs_top = plt.GridSpec(nrows=n_rows+1, ncols=3,  figure=fig, top=0.95)
        gs_base = plt.GridSpec(nrows=n_rows+1, ncols=3,  figure=fig, hspace=0,top=0.91)


        #gs = fig.add_gridspec(nrows=n_rows+1, ncols=3, left=0.05, figsize=(8.25, 11.67))
        #gs = gridspec.GridSpec(nrows=n_rows+1, ncols=3,  figure=fig)
        
        wave_ax = self.addFullSubplot(cfg_par,fig,gs_top,vel,y,result,noise,i,j,lineInfo,singleVorBinInfo)
        
 
        #for plot_count in range(n_plots):
        k=0
        for i in xrange(0,len(lineInfo['Wave'])):


            if i == 0:
                j = 0
            elif i % 3 == 0:
                j +=1 
                k =0

            waveMin =  np.log(lineInfo['Wave'][i] - lineInfo['lineRangeAng'][i])
            waveMax =  np.log(lineInfo['Wave'][i] + lineInfo['lineRangeAng'][i])
            idxMin = int(np.where(abs(vel-waveMin)==abs(vel-waveMin).min())[0]) 
            idxMax = int(np.where(abs(vel-waveMax)==abs(vel-waveMax).min())[0] )
        
            x_data_plot = velExp[idxMin:idxMax]
            x_data_plot = cvP.lambdaVRad(x_data_plot,lineInfo['Wave'][i])
            y_data_plot = y[idxMin:idxMax]
            y_BFit_plot = result.best_fit[idxMin:idxMax]
            y_Res_plot = result.residual[idxMin:idxMax]
            y_sigma_plot = noise[idxMin:idxMax]

            # Calculate axis limits and aspect ratio
            x_min = np.min(x_data_plot)
            x_max = np.max(x_data_plot)
            if cfg_par['gPlot']['fixed_scale']:
                y1_min = -75
                y1_max = np.nanmax(y_data_plot)*1.1
            else:
                y1_min = np.nanmin(y_data_plot)*1.1
                y1_max = np.nanmax(y_data_plot)*1.1

            minLoc = np.floor(np.min(x_data_plot)/np.float(cfg_par['gPlot']['deltaVlabel']))
            maxLoc = np.ceil(np.max(x_data_plot)/np.float(cfg_par['gPlot']['deltaVlabel']))
            
            xTicks = np.arange(minLoc*cfg_par['gPlot']['deltaVlabel'],(maxLoc+1)*cfg_par['gPlot']['deltaVlabel'],
                cfg_par['gPlot']['deltaVlabel'])
            
            if np.abs(xTicks[0]) > xTicks[-1]:
                    xTicks=np.delete(xTicks,0)
            elif np.abs(xTicks[0]) < xTicks[-1]:
                    xTicks = np.delete(xTicks,-1)


            xTicksStr = [str(hh) for hh in xTicks]
        
            ax = fig.add_subplot(gs_base[j+1,k])
            

            ax.set_xticks(xTicks)
            ax.set_xticklabels([])

            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y1_min, y1_max)
            ax.xaxis.labelpad = 6
            ax.yaxis.labelpad = 10
            ax.minorticks_on()
            ax.tick_params(axis='both', bottom='on', top='on',
                                   left='on', right='on', which='major', direction='in')
            ax.tick_params(axis='both', bottom='on', top='on',
                                   left='on', right='on', which='minor', direction='in')
            if k==0:
                ylabh = ax.set_ylabel(
                    r'Flux [-]')
                ylabh.set_verticalalignment('center')

            ax.tick_params(axis='both', which='major', pad=5)
            #ax1.xaxis.set_minor_locator()
            #ax1.yaxis.set_minor_locator()      
            if lineInfo['Name'][i] == 'Hb':
                lineInfoName =r'H$_\beta$'
            elif lineInfo['Name'][i] == 'Ha':
                lineInfoName =r'H$_\alpha$'
            else:
                lineInfoName =lineInfo['Name'][i]

            ax.step(x_data_plot, y_data_plot, where='mid', color='black', linestyle='-')
            ax.plot(x_data_plot, y_BFit_plot, 'r-', label=lineInfoName+str(int(lineInfo['Wave'][i])))
            #ax2.fill_between(0, y_BFit_plot, y_sigma_plot,
            #                 facecolor='grey', alpha=0.5,step='mid')
            #ax1.step(vel, yInFit, 'b-', label='init fit')

            #ax1.fill_between(vel, yBFit-dely, yBFit+dely, color="#ABABAB",
            #        label='3-$\sigma$ uncertainty band')
            if cfg_par['gFit']['modName'] !='g1':
                comps = result.eval_components()
                for ii in xrange(0,len(lineInfo['ID'])):

                    ax.plot(x_data_plot, comps['g1ln'+str(ii)+'_'][idxMin:idxMax], 'g--')
                
                    if cfg_par['gFit']['modName'] =='g2':
                        ax.plot(x_data_plot, comps['g2ln'+str(ii)+'_'][idxMin:idxMax], 'm--')    
                
                    elif cfg_par['gFit']['modName'] !='g2':
                        ax.plot(x_data_plot, comps['g2ln'+str(ii)+'_'][idxMin:idxMax], 'm--')    
                        ax.plot(x_data_plot, comps['g3ln'+str(ii)+'_'][idxMin:idxMax], 'c--')    


            ax.axvline(color='k', linestyle=':', zorder=0)                           
            legend = ax.legend(loc='best',handlelength=0.0, handletextpad=0.0,frameon=False)
            legend.get_frame().set_facecolor('none')

 
            divider = make_axes_locatable(ax)
            ax2 = divider.append_axes("bottom", size='20%',pad=0)
            ax2.minorticks_on()

            ax.figure.add_axes(ax2)
            # Calculate axis limits
            x_min = np.nanmin(x_data_plot)
            x_max = np.nanmax(x_data_plot)
            if cfg_par['gPlot']['Res-fixed_scale']:
                y1_min = np.nanmin([-200.,np.nanmin(-y_sigma_plot)*1.5,np.nanmin(-y_Res_plot)*1.5])
                y1_max = np.nanmax([+200.,np.nanmax(+y_sigma_plot)*1.5,np.nanmax(+y_Res_plot)*1.5])
            else:
                y1_min = np.nanmin(y_Res_plot)*1.1
                y1_max = np.nanmax(y_Res_plot)*1.1

            ax2.set_xticks(xTicks)

            #ax2.set_yticks([-150,0,150])
            #ax2.set_yticklabels([])

            # Set axis limits
            ax2.set_xlim(x_min, x_max)
            ax2.set_ylim(y1_min, y1_max) 

            #ax2.plot(vel, amp(x, p(x,m,n)))
            ax2.step(x_data_plot, y_Res_plot, 'g-', label='residuals')
            ax2.axhline(color='k', linestyle=':', zorder=0)                           
            ax2.axvline(color='k', linestyle=':', zorder=0)                           
            ax2.fill_between(x_data_plot, -y_sigma_plot, y_sigma_plot,
                             facecolor='grey', alpha=0.5,step='mid')

           # for the last plot add the x-axis label
            if i >= len(lineInfo['Wave'])-3:                
                ax2.set_xticks(xTicks)
                ax2.set_xlabel(
                        r'v$_\mathrm{radio}$\,[$\mathrm{km}\,\mathrm{s}^{-1}$]', labelpad=2)
            else:
                ax2.set_xticklabels([])


            k+=1
        #delete unused subplots
        #i+=1
        #while not i % 3 ==0:   
        #    fig.delaxes(ax)

        #    gs[j-1][k].axes.get_xaxis().set_ticklabels(xTicksStr)
        #    gs[j-1][k].axes.get_xaxis().set_ticklabels(xTicksStr)
        #    k +=1
        #    i +=1


        plt.savefig(outPlot,dpi=200,bbox_inches='tight',
                    format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        #plt.show()
        plt.close()  

    def addFullSubplot(self,cfg_par,fig,gs,vel,y,result,noise,xx,yy,lineInfo,singleVorBinInfo):
        
        velPlot = np.exp(vel)
        yBFit = result.best_fit
        yRes = result.residual
        yInFit = result.init_fit
        key = 'general'
        
        #outPlotDir = cfg_par[key]['runNameDir']+'/'+cfg_par['gPlot']['outPlotDirName']
        #if not os.path.exists(outPlotDir):
        #    os.mkdir(outPlotDir)


        #outPlot = outPlotDir+str(xx)+'_'+str(yy)+'_'+cfg_par['gFit']['modName']+'.png'       
        

        #dely = result.eval_uncertainty(sigma=3)
            
         # initialize figure
        #params = self.loadRcParamsBig()
        #plt.rcParams.update(params)
        #fig = plt.figure(figsize =(10,8))
        #fig.subplots_adjust(hspace=0.0)
        #gs = gridspec.GridSpec(1, 1)
        #plt.rc('xtick')

        # Initialize subplots
        ax1 = fig.add_subplot(gs[0,:])

        divider = make_axes_locatable(ax1)
        ax2 = divider.append_axes("bottom", size='15%', pad=0)
        ax1.figure.add_axes(ax2)
        

        #ax.set_xticks([])

        ax1.set_xlabel(r'Wavelength [$\AA$]',labelpad=15)
        ax1.set_ylabel(r'Flux [-]')


        # Calculate axis limits and aspect ratio
        x_min = np.min(velPlot)
        x_max = np.max(velPlot)
        if cfg_par['gPlot']['fixed_scale']:
            y1_min = np.min(y)*1.2
            y1_max = np.max(y)*1.2
        else:
            y1_min = np.nanmin(y)*1.1
            y1_max = np.nanmax(y)*1.1

        # Set axis limits
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(y1_min, y1_max)

        ax1.tick_params(axis='both', which='major', pad=5)
        #ax1.xaxis.set_minor_locator()
        #ax1.yaxis.set_minor_locator()      
        
        ax1.step(velPlot, y, where='mid', color='black', linestyle='-')
        ax1.plot(velPlot, yBFit, 'r-', label='best fit')
        #ax1.step(vel, yInFit, 'b-', label='init fit')
        aicStr = str(int(result.aic))
        bicStr = str(int(result.bic))
        redchiStr = str(int(result.redchi))
        successStr = str(result.success)       
        xText = cfg_par['gFit']['lambdaMin']+50


        ax1.text(xText, y1_max*0.90, r'BIN ID:\t'+str(singleVorBinInfo['BIN_ID'][0]), {'color': 'k', 'fontsize': 8})
        ax1.text(xText, y1_max*0.80, r'X,Y:\t'+str(xx)+','+str(yy), {'color': 'k', 'fontsize': 8})

        #ax1.text(xText, x_max*0.85, r'Success:\t'+successStr, {'color': 'b'})
        ax1.text(xText, y1_max*0.70, r'$\tilde{\chi}^2$:\t'+redchiStr, {'color': 'k', 'fontsize': 8})
        ax1.text(xText, y1_max*0.60, r'aic:\t'+aicStr, {'color': 'k', 'fontsize': 8})
        ax1.text(xText, y1_max*0.50, r'bic:\t'+bicStr, {'color': 'k', 'fontsize': 8})

        #ax1.fill_between(vel, yBFit-dely, yBFit+dely, color="#ABABAB",
        #        label='3-$\sigma$ uncertainty band')
        if cfg_par['gFit']['modName'] !='g1':
            comps = result.eval_components()
            for i in xrange(0,len(lineInfo['ID'])):
                
                ax1.plot(velPlot, comps['g1ln'+str(i)+'_'], 'g--')
            
                if cfg_par['gFit']['modName'] =='g2':
                    ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
            
                elif cfg_par['gFit']['modName'] !='g2':
                    ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
                    ax1.plot(velPlot, comps['g3ln'+str(i)+'_'], 'c--')    



        # Calculate axis limits and aspect ratio
        x_min = np.min(velPlot)
        x_max = np.max(velPlot)
        print noise
        print yRes
        if cfg_par['gPlot']['Res-fixed_scale']:
            y1_min = np.nanmin([-200.,np.nanmin(-noise)*1.1,np.nanmin(-yRes)*1.1])
            y1_max = np.nanmax([+200.,np.nanmax(+noise)*1.1,np.nanmax(+yRes)*1.1])
        else:
            y1_min = np.nanmin(yRes)*1.1
            y1_max = np.nanmax(yRes)*1.1  

        # Set axis limits
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y1_min, y1_max) 

        #ax2.plot(vel, amp(x, p(x,m,n)))
        ax2.step(velPlot, yRes, 'g-', label='residuals')
        ax2.axhline(color='k', linestyle=':', zorder=0)                           
        ax2.fill_between(velPlot, -noise, noise,
                         facecolor='grey', alpha=0.5,step='mid')

        #ax1.legend(loc='best')


        #plt.savefig(outPlot,
        #            format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        #plt.show()
        #plt.close()
           
        return ax1