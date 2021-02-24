#!/usr/bin/env python3.6
'''

Set of tools for computing line ratios and make BPT plots, electron density maps

Requirements
------------
Moment maps and fit results must exist. 

'''

import os, sys
import yaml


from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy.wcs import WCS

import numpy as np
from scavengers import util as ut
from scavengers import tPlay,cvPlay,momPlot

tP = tPlay.tplay()
cvP = cvPlay.convert()
mPl = momPlot.MOMplot()
class lineplay:
    '''Modules to create different line ratio maps
    
    - makeElectronDensity
        builds electron density map from the [SII] doublet ratio, creates n_e column in lineRatios subtable
    - plotElectronDensity
        plots electron density map produced by makeElectronDensity()     
    '''

    def electronDensityMap(self,cfg_par):
        '''Computes the electron density from the [SII]6716/[SII]6731 ratio of the flux density maps of each gaussian component,
        plots the maps using mPl.momTriplet()
        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        outPlotDir: str, optional
            output directory of plot

        Notes
        ------

        Empirical formula for electron density taken from : 10.1051/0004-6361/201322581
        previously published in Osterbork & Farland 2003, Astrophysics of Gaseous Nebulae and Active Galactic Nuclei

        '''
    
        neDir =cfg_par['general']['bptDir']+'electronDensity/'
        if not os.path.exists(neDir):
            os.mkdir(neDir)
        cfg_par['general']['bptDir'] = neDir

        modName=cfg_par['gFit']['modName']
        momModDir = cfg_par['general']['momDir']+modName+'/'

        SII6716='SII6716'
        SII6730='SII6730'

        if modName=='g1':    
            mom0Names6716=[momModDir+'mom0_g1-'+SII6716+'.fits']
            mom0Names6730=[momModDir+'mom0_g1-'+SII6730+'.fits']
        else:    
            mom0Names6716=[momModDir+'mom0_g1-'+SII6716+'.fits',momModDir+'mom0_g2-'+SII6716+'.fits',momModDir+'mom0_tot-'+SII6716+'.fits']
            mom0Names6730=[momModDir+'mom0_g1-'+SII6730+'.fits',momModDir+'mom0_g2-'+SII6730+'.fits',momModDir+'mom0_tot-'+SII6730+'.fits']

        nEName=[]
        
        for i in range(0,len(mom0Names6730)):

            mom06716=fits.getdata(mom0Names6716[i])
            mom06730=fits.getdata(mom0Names6730[i])

            R=np.divide(mom06716,mom06730)
            nEHead=fits.getheader(mom0Names6730[i])

            R2=np.power(R,2)
            R3=np.power(R,3)

            C1=0.0543*np.tan(-3.0553*R+2.8506)+6.98
            C2=-10.6905*R
            C3=+9.9186*R2
            C4=-3.5442*R3

            nE=C1+C2+C3+C4

            if i==0:
                iStr='g1'
            elif i==1:
                iStr='g2'
            else:
                iStr='tot'

            nE=np.power(10,nE)

            idxMin = nE<40.
            nE[idxMin] = np.nan
            idxMax = nE>1e4
            nE[idxMax] = np.nan

            nEName.append(neDir+'nE_'+iStr+'.fits')

            fits.writeto(nEName[i],np.log10(nE),nEHead,overwrite=True)

        if modName=='g1':
            mPl.mom0plot(cfg_par)
        else:
            mPl.momTriplet(cfg_par,nEName[0],nEName[1],nEName[2],interpolation='nearest',
    shareScale=True,kind='electronDensity')



    def dustExtinctionMap(self,cfg_par):
        '''Computes the total dust extinction in V band AV from the measured Balmer decrement (Halfa/Hbeta flux ratio)
        of each gaussian component, plots the maps using mPl.momTriplet()
        

        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        outPlotDir: str, optional
            output directory of plot

        Notes
        ------

        Av and AHalpha are computed as in Dominguez A, et al. 2012 
        Momcheva, et al. 2013
        https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145
        https://iopscience.iop.org/article/10.1088/0004-6256/145/2/47
        '''
    
        neDir =cfg_par['general']['bptDir']+'dustExtinction/'
        if not os.path.exists(neDir):
            os.mkdir(neDir)
        cfg_par['general']['bptDir'] = neDir

        neDirPlot =cfg_par['general']['bptDir']+'plots/'
        if not os.path.exists(neDirPlot):
            os.mkdir(neDirPlot)
        cfg_par['general']['bptPlotDir'] = neDirPlot

        modName=cfg_par['gFit']['modName']
        momModDir = cfg_par['general']['momDir']+modName+'/'

        Hbeta='Hb4861'
        Halfa='Ha6562'

        if modName=='g1':    
            mom0NamesHalfa=[momModDir+'mom0_g1-'+Halfa+'.fits']
            mom0NamesHbeta=[momModDir+'mom0_g1-'+Hbeta+'.fits']
        else:    
            mom0NamesHalfa=[momModDir+'mom0_g1-'+Halfa+'.fits',momModDir+'mom0_g2-'+Halfa+'.fits',momModDir+'mom0_tot-'+Halfa+'.fits']
            mom0NamesHbeta=[momModDir+'mom0_g1-'+Hbeta+'.fits',momModDir+'mom0_g2-'+Hbeta+'.fits',momModDir+'mom0_tot-'+Hbeta+'.fits']

        AHalphaName=[]
        AvName=[]
        for i in range(0,len(mom0NamesHalfa)):

            mom0Halfa=fits.getdata(mom0NamesHalfa[i])
            mom0Hbeta=fits.getdata(mom0NamesHbeta[i])

            R=np.divide(mom0Halfa,mom0Hbeta)
            nEHead=fits.getheader(mom0NamesHalfa[i])
            nEHead['BUNIT'] = 'mag'
            EBV = 1.97*np.log10(R/2.86)

            if i==0:
                iStr='g1'
            elif i==1:
                iStr='g2'
            else:
                iStr='tot'

            AHalpha = 3.33*EBV
            Av = 4.05*EBV

            index = np.power(10,0.4*AHalpha)<0.

            AHalpha[index] = 1.

            index = np.power(10,0.4*AHalpha)>5.

            AHalpha[index] = 5.

            AHalphaName.append(neDir+'AHalpha_'+iStr+'.fits')
            
            AvName.append(neDir+'Av_'+iStr+'.fits')

            fits.writeto(AHalphaName[i],np.power(10,0.4*AHalpha),nEHead,overwrite=True)
            fits.writeto(AvName[i],Av,nEHead,overwrite=True)



        if modName=='g1':
            mPl.mom0plot(cfg_par)
        else:
            mPl.momTriplet(cfg_par,AHalphaName[0],AHalphaName[1],AHalphaName[2],interpolation='nearest',
    shareScale=True,kind='dustExtinction')
            #mPl.momTriplet(cfg_par,AVName[0],AVName[1],AHalphaName[2],interpolation='nearest',
    #shareScale=True,kind='dustExtinction')


    # def plotElectronDensity(self,cfg_par,im1,im2,im3,inMom0=None,
    # cRange=None,imLevels=None,imContColors=['black','black','black'],beamCoords=None,
    # contName=None,contValues=None,contColors=None,interpolation=None,title=False,
    # shareScale=True,kind='mom1'):
    #     '''Plots the electron density map for each gaussian component and the for the total fitted line

    #     Parameters
    #     ----------

    #     cfg_par: OrderedDictionary
    #         Dictionary with alla parameters or gufo. 
    #         Parameters are given from terminal, or in a `.yml' file.


    #     Returns
    #     ----------
    #     outFig: str
    #         full path to output figure
    #     '''

    #     objCoordsRA = cfg_par['moments']['centreRA']
    #     objCoordsDec = cfg_par['moments']['centreDec']
        
    #     centre = SkyCoord(ra=objCoordsRA*u.degree, dec=objCoordsDec*u.degree, frame='fk5')
    #     size = u.Quantity((cfg_par['moments']['sizePlots'],cfg_par['moments']['sizePlots']), u.arcmin)
  

    #     params = ut.loadRcParams()
    #     plt.rcParams.update(params)


    #     fig = plt.figure( constrained_layout=False)
    #     fig.set_tight_layout(False)

        
    #     gs = plt.GridSpec(nrows=1, ncols=2,  figure=fig,wspace=0.0,hspace=0.0)
    
    #     ax1 = fig.add_subplot(gs[0,0],projection=hduImCut1.wcs)






