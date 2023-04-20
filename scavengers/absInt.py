#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
from astropy import wcs
from astropy import units as u

from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib import rc
from scavengers import util as ut


import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap

from scavengers import cvPlay, headPlay, fitsPlay
from scavengers import util as ut
from scavengers import plotUtil as plUt


cvP = cvPlay.convert()
hP = headPlay.headplay()
fP = fitsPlay.fitsplay()

class absint:


    def __init__(self):

        self.rootdir = os.getcwd()+'/'
        #self.filename = sys.argv[1]
        self.T = 100.
        self.nhi = 1.8216E18
        self.mp=1.67492728E-24   #g
        self.msun=1.98855e33      #g
        self.pc = 3.08577758149137e18
    def zaxis(self,cubename):

        cubefile = fits.open(cubename)  # read input
        hdr = cubefile[0].header
                
        freq = (np.linspace(1, hdr['NAXIS3'], hdr['NAXIS3']) - hdr['CRPIX3']) * hdr['CDELT3'] + hdr['CRVAL3']

        return freq 

    def measContFlux(self,cfg_par,contIm):
        
        contData = fits.getdata(contIm)
        contHead = hP.cleanHead(contIm,writeFile=True)
        contHead = fits.getheader(contIm)

        w = wcs.WCS(contHead)    
        objCoordsRA = cfg_par['HIabs']['centreRA']
        objCoordsDec = cfg_par['HIabs']['centreDec']
        ra=cvP.hms2deg(objCoordsRA)
        dec=cvP.dms2deg(objCoordsDec)
        objCoordsRA = ra
        objCoordsDec = dec
        print(dec)
        #convert coordinates in pixels
        #cen_x,cen_y=w.wcs_world2pix(ra,dec,0)
        cen_x, cen_y = w.wcs_world2pix(ra, dec, 1)
        print(cen_x,cen_y)
        contFlux = np.nanmax(contData[int(np.round(cen_y)-4):int(np.round(cen_y)+4),int(np.round(cen_x-4)):int(np.round(cen_x+4))])
        contFluxScience="{:.2e}".format(contFlux*1e3)
        if 'contFlux' in cfg_par['HIabs'] :
            contFlux = cfg_par['HIabs']['contFlux']
            contFluxScience="{:.2e}".format(cfg_par['HIabs']['contFlux']*1e3)
        return contFlux, contFluxScience 

    def optical_depth(self,sabs,scont):
        '''
        Estimates the optical depth of an absorption line
        INPUT
            scont: continuum flux (float)
            sabs: flux of the absorbed component (negative float)
        OUTPUT
            hi.tau: optical depth of the line   
        '''

        tau=np.log(1.-(-np.divide(sabs,(float(scont)))))
        if tau.size == 1:
            print('Optical depth = '+str(round(tau,3)))

        return tau

    def nhiAbs(self,tau, dv):
        '''
        Estimates the column density of the absorption line
        Parameters:
            tau: optical depth of the line (float)
            dv: width of the line in km/s
        
        Returns:
            hi.nhi_abs: column density of the absorption line in cm-2
        '''

        nhiabs = self.nhi*self.T*tau*dv
        #print('N(HI) = '+str(np.round(nhiabs, 3))+' cm-2')

        return nhiabs

    def mhi_abs(self,nhi_abs, bMaj,bMin, dL,z):
        '''Estimates the mass of the absorbed HI
        Parameters:
            nhi_abs: column density of the absorption line in cm-2
            area: area over which the column density is integrated in pc^2 (output of hi.beam_area)
        
        Returns:
            hi.mhi_abs: hi mass inferred by the absorption line in Msun
        '''
        
        #r = ang * dl / (RAD2DEG * 3600 * (1+z)**2) # Mpc
        bMaj_pc = (bMaj.arcsecond*dL/(180./np.pi*3600.*np.power(1+z,2)))*1e6
        bMin_pc = (bMaj.arcsecond*dL/(180./np.pi*3600.*np.power(1+z,2)))*1e6

        area = bMaj_pc*bMin_pc#pc

        area*=np.power(self.pc,2)#cm-2
        mhiabs = area*self.mp*nhi_abs/self.msun

        #print('M(HI) = '+str(round(mhiabs, 6)/1e8)+' x10^8 mSun')

        return mhiabs

        
    def fwzi(self,cfg_par,vel,flux):


        peakFlux,idxPeak = self.findPeak(cfg_par,vel,flux)
        idxPeak = np.abs(flux-peakFlux).argmin()
        fluxLeft=flux[0:idxPeak][::-1]
        zeroLeft=np.argmax(fluxLeft>0)

        fluxRight=flux[idxPeak:]
        zeroRight=np.argmax(fluxRight>0)+1
        width=-(zeroRight+zeroLeft)*np.nanmean(np.diff(vel))

        return width,zeroLeft,zeroRight

    def findPeak(self,cfg_par,vel,flux):
        velLineMin=cfg_par['galaxy']['vsys']-cfg_par['HIabs']['stats']['velRange']
        velLineMax=cfg_par['galaxy']['vsys']+cfg_par['HIabs']['stats']['velRange']
        
        if vel[0]>vel[-1]:
            idxMax = np.abs(vel-velLineMin).argmin()
            idxMin = np.abs(vel-velLineMax).argmin()
        else:
            idxMin = np.abs(vel-velLineMin).argmin()
            idxMax = np.abs(vel-velLineMax).argmin()            
        peakFlux=np.nanmin(flux[idxMin:idxMax])

        idxPeak = np.abs(flux-peakFlux).argmin()

        return peakFlux, idxPeak

    def abSex(self,cubeName,specName,fluxCont,ra,dec,raN,decN,vunit='m/s'):
        '''
        Extract spectra from all l.o.s. exctracted using a catalog of sources or a source finder
        WARNING:
            source finder, or extraction of sources from a catalog (cont_src module) must be run first
        INPUT:
            dictionary of parameter file
        OUTPUT:
            spectra are saved in ascii format in cfg_par['general']['workdir']+/spec/
            spectra have the following columns: 
            frequency, flux, noise (MADFM), optical depth, optical depth noise, mean noise
        
        OPTIONS:
            chromatic aberration correction
            continuum subtraction
            hanning smoothing
        '''
            
        #verb = cfg_par['general']['verbose']

        #cubename = cfg_par['general'].get('cubename',None)
        cubefile = fits.open(cubeName)  # read input
        hdr = cubefile[0].header
        sci = cubefile[0].data 
        sci = sci.squeeze()
        x = hdr['NAXIS1']
        y = hdr['NAXIS2']
        z = hdr['NAXIS3']
        cen_imx = hdr['CRPIX1']
        cen_imy = hdr['CRPIX2']
        freq0 = hdr['CRVAL3']
        freq_del = hdr['CDELT3']

        #key = 'source_catalog'
        #if cfg_par['source_catalog'].get('enable',False) == True:

        #    catalog_table = str(cfg_par['general'].get('absdir')) + 'cat_src_sharpener.txt'
           
        #    tab = ascii.read(catalog_table)
        
        #    if cfg_par[key].get('catalog', 'NVSS') == 'NVSS':

        #        J2000_name = tab['NVSS']
        #        ra = tab['RAJ2000']
        #        dec = tab['DEJ2000']
        #        flux_cont = tab['S1.4']*1e-3

        #elif cfg_par['source_finder'].get('enable',False) == True:

        #    src_list_csv = cfg_par['general']['absdir']+'mir_src_sharpener.csv'
            # open file
        #    src_list_vec = ascii.read(src_list_csv)
        #    J2000_name = np.array(src_list_vec['J2000'],dtype=str)
        #    ra = np.array(src_list_vec['ra'],dtype=str)
        #    dec = np.array(src_list_vec['dec'],dtype=str)
        #    flux_cont = np.array(src_list_vec['peak'],dtype=float)
        

        pixels = fP.coordToPix(cubeName,ra,dec, verbose=False)

        pixelsN = fP.coordToPix(cubeName,raN,decN, verbose=False)
        print(pixelsN)

        #key = 'spec_ex'

        #src_id = np.arange(0,ra.size+1,1)
        #src_id = src_list_vec['ID']
        freq = self.zaxis(cubeName)

        abs_mean_rms = np.zeros(pixels.shape)
        abs_los_rms = np.zeros(pixels.shape)
        tau_los_rms = np.zeros(pixels.shape)
        outNames = []
        count_thresh =0
        count_fov = 0
        count_blanks = 0
        average_noise = []
        #for i in xrange(0,pixels.shape[0]):

            # extract spectrum from each line of sight
        flux = np.zeros(freq.shape[0])
        madfm = np.zeros(freq.shape[0])

        # better to use the numpy function
        #if str(pixels[i,0]) == 'nan' or str(pixels[i,1]) == 'nan':
        #if np.isnan(pixels[i, 0]) or np.isnan(pixels[i, 1]):
        #    count_thresh +=1
        #    pass

        #elif (0 < int(pixels[i,0]) < x and
        #        0 < int(pixels[i,1]) < y): 
                
        pix_x_or = int(pixels[0])
        pix_y_or = int(pixels[1])

        pix_x_orN = int(pixelsN[0])
        pix_y_orN = int(pixelsN[1])


        for j in range(0, z):
            #chrom_aber = cfg_par[key].get('chrom_aberration', False)
            #correct for chromatic aberration
            #if chrom_aber == True:

            #    if (cfg_par[key].get('zunit','Hz') == 'm/s'):
            #        freq_real= freq* 1e2
            #        freq_real = (kk.C*kk.HI) /  (freq_real + kk.C)
            #        freq_real0 = (kk.C*kk.HI) /  (hdr['CRVAL3']*1e2 + kk.C)
            #        freq_del = (freq_real0 - freq_real[-1] )/ len(freq_real)
                #depending if the cube is in velocity or frequency ?
            #        scale = (freq_real0 - j*freq_del) / freq_real0

            #    pix_x = (pix_x_or - hdr['CRPIX1']) * scale + hdr['CRPIX1']
            #    pix_y = (pix_y_or - hdr['CRPIX2']) * scale + hdr['CRPIX2']
                #print('before rounding: x={0:.3f}, y={1:.3f}'.format(pix_x, pix_y))
            #    pix_x = int(round(pix_x,0))
            #    pix_y = int(round(pix_y,0))
            #else:
            pix_x = pix_x_or
            pix_y = pix_y_or
            
            pix_xN = pix_x_orN
            pix_yN = pix_y_orN
            

            if  (0 < pix_x < x and
                 0 < pix_y < y): 
                flux[j] = sci[j, pix_y, pix_x]
            else:
                flux[j] = 0.0
            
            #print('x={0:d}, y={1:d}, flux={2:.5f}'.format(pix_x, pix_y, flux[j]))


            # determine the noise of the spectrum [Whiting 2012 et al.] in each channel
            # MADMF: median absolute deviation from the median
            # extract a region were to determine the noise: A BOX around the l.o.s.
            if (pix_xN+20 < hdr['NAXIS1'] and  pix_xN-20 > 0 and
               pix_yN+20 < hdr['NAXIS2'] and pix_yN - 20 > 0):
                    rms = np.nanmedian(sci[j, pix_xN -20:pix_xN + 20, pix_yN - 20:pix_yN + 20])
                    #print pix_x +10:pix_x + 10, pix_y - 5:pix_y + 5]
                    #print rms
                    if rms != 0.0:
                        med2 = np.abs(sci[j, pix_yN, pix_xN] - rms)
                        madfm[j] = np.nanmedian(med2) / 0.6744888
                    else:
                        madfm[j] = 0.0
            else:
                madfm[j] = 0.0


            abs_mean_rms = np.nanmean(madfm) 

        if np.nansum(flux) == 0.:
            count_blanks +=1
            #if verbose == True:
                #print '# Blank spectrum:\t'+str(src_id[i])+' '+J2000_name[i]+' #'
            #    print '# Blank spectrum:\t'

            #continue

        # measure noise in the spectrum outside of the line
        end_spec = float(sci.shape[0])
        end_spec_th = int(end_spec/3.)
        end_spec = int(end_spec)
        mean_rms = (np.std(flux[0:end_spec_th]) +
                            np.std(flux[end_spec-end_spec_th:end_spec])) / 2.
        mean_rms_arr = np.zeros(sci.shape[0])+mean_rms
        
        average_noise.append(mean_rms)
        tau = self.optical_depth(flux, fluxCont)
        if vunit == 'm/s':
            freq_del/=1e3
        nhi = self.nhiAbs(tau,freq_del)
        if np.nansum(madfm)!= 0.0:
            tau_noise = self.optical_depth(madfm, fluxCont)
            nhi_noise = self.nhiAbs(tau_noise,freq_del)
        else:
            tau_noise = np.zeros(sci.shape[0])
            nhi_noise = np.zeros(sci.shape[0])

        #write spectrum
        #out_spec = str(cfg_par['general']['specdir']+str(src_id[i])+'_J'+J2000_name[i])+'.txt'
        #out_spec = "{0:s}{1:02d}_J{2:s}.txt".format(
        #                   cfg_par['general']['specdir'], src_id[i], J2000_name[i])
        out_spec = specName+'.txt'
        outNames.append(out_spec)

        #flag_chans = cfg_par[key].get('flag_chans', None)
        #if flag_chans != None:
        #    index_flags_l = (np.abs(freq - flag_chans[0])).argmin()
        #    for k in xrange(1,len(flag_chans)):
        #        index_flags = (np.abs(freq - flag_chans[k])).argmin()
        #        flux[index_flags_l:index_flags] = np.nan
        #        index_flags_l = index_flags

        #if cfg_par[key].get('zunit','Hz') == 'm/s':
        #    xcol = 'Velocity [m/s]'
        #elif cfg_par[key].get('zunit','Hz') == 'km/s':
        xcol = 'Velocity [km/s]'
        #elif cfg_par[key].get('zunit','Hz') == 'MHz':
        #    xcol = 'Frequency [MHz]'
        #else:
        #    xcol = 'Frequency [Hz]'

        t = Table([freq, flux, madfm, tau, tau_noise, nhi, nhi_noise, mean_rms_arr], 
            names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'NHI', 'noise_NHI', 'Mean noise [Jy]'),
            meta={'name': 'Spectrum'})
        ascii.write(t,out_spec,overwrite=True)
        #if verb==True:
            #print '# Extracted spectrum: \t' +str(src_id[i])+' '+J2000_name[i]+' #'
        #    print '# Extracted spectrum: \t'
            #polysub = cfg_par['polynomial_subtraction'].get('enable', False) 
            #if polysub == True:

            #    deg = cfg_par['polynomial_subtraction'].get('degree_pol',3)
            #    sub_flux = poly_sub(cfg_par,freq,flux,deg)
            #    sub_madfm = madfm.copy()
            #    sub_od = tau.copy()
            #    sub_noise_od = tau_noise.copy()

            #    out_spec_polysub = string.split(out_spec,'.')[0]
            #    out_spec_polysub= out_spec_polysub+'_psub.txt'

            #    t = Table([freq, sub_flux, sub_madfm, sub_od, sub_noise_od, mean_rms_arr], 
            #        names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
            #        meta={'name': 'Spectrum'})
            #    ascii.write(t,out_spec_polysub,overwrite=True)

            #dohan = cfg_par['hanning'].get('enable', False)
            #if dohan == True and polysub == False:
                
                #window = cfg_par['hanning'].get('window', 1)
                #han_flux  = hanning_spec(flux)
                #han_madfm = hanning_spec(madfm)
                #han_od = hanning_spec(tau)
                #han_noise_od = hanning_spec(tau_noise)
                
                #out_spec_han = string.split(out_spec,'.')[0]
                #out_spec_han= out_spec_han+'_han.txt'

                #t = Table([freq, han_flux, han_madfm, han_od, han_noise_od, mean_rms_arr], 
                #    names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
                #    meta={'name': 'Spectrum'})
                #ascii.write(t,out_spec_han,overwrite=True)

            #elif polysub == True and dohan == True:

            #    han_sub_flux  = hanning_spec(sub_flux)
            #    han_sub_madfm = hanning_spec(sub_madfm)
            #    han_sub_od = hanning_spec(sub_od)
            #    han_sub_noise_od = hanning_spec(sub_noise_od)
                
            #    out_spec_han = string.split(out_spec,'.')[0]
            #    out_spec_han= out_spec_han+'psub_han.txt'


            #    t = Table([freq, han_sub_flux, han_sub_madfm, han_sub_od, han_sub_noise_od, mean_rms_arr], 
            #        names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
            #        meta={'name': 'Spectrum'})
            #    ascii.write(t,out_spec_han,overwrite=True)

        # close fits file
        cubefile.close()
        
        #print '\n# Total number of sources: \t'+str(pixels.shape[0])
        #print '# Sources flagged: \t\t'+str(count_thresh)
        #print '# Blank spectra:\t\t'+str(count_blanks)
        #print '# Total number of spectra: \t'+str(pixels.shape[0]-count_thresh-count_fov-count_blanks)
        #print '# Average noise in spectra: \t'+str(round(np.nanmean(average_noise)*1e3,1))+' mJy/beam'

        return out_spec

    def absPlotInt(self,cfg_par,specName,contFlux=None,xLims=None,yLims=None,vSys=None,y_sigma=None,pxBeam=None,source=None,opt='abs'):
        '''Plots an integrated line profile. If contFlux is given, the y-axis is in optical depth.

        Parameters
        ----------

        cfg_par: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        specName: str, 
            _default=None_, full path to spectrum. 2-columns txt file: velocity vs flux
        
        contFlux: float, optional
            _default=None_, integrated continuum flux in Jy. Used to convert the spectrum in optical depth.

        xLims: list, optional
            _default=None_, x-axis range: [x_min, x_max]

        yLims: list, optional
            _default=None_, y-axis range: [y_min, y_max]

        vSys: float, optional
            _default=None_, systemic velocity of source in km/s

        y_sigma: np.array, optional
            _default=None_, noise array (must have the same lenght as spectrum)

        opt: str, optional
            _default='abs'_, kind of spectrum to plot: 'abs': integrated absorption, 'emSof': emission spectrum from SoFiA,
            'em': generic emission spectrum

        Returns
        ----------
        outFig: str
            full path to output figure

        Notes
        ----------
        Conversion to optical depth using self.optical_depth()

        '''

        params = plUt.loadRcParams('fwR')
        plt.rcParams.update(params)


        if opt=='abs':
            spec=np.genfromtxt(specName)

            vel=spec[:,0]
            x_data=vel
            flux=spec[:,1]
            y_data=flux
            
        elif opt =='emSof':
            spec=ascii.read(specName)

            if 'freq' in spec.dtype.names:
                x_data = spec['freq']
                xLabel = r'Freq. [Hz]'
            elif 'vel' in spec.dtype.names:
                x_data = spec['vel']/1e3
            else:
                x_data = spec['felo']/1e3
            flux = spec['f_sum']
            n_pix=spec['n_pix']
    
            y_data=flux*1e3
         
            y_sigma *=1e3
            y_sigma *=np.sqrt(n_pix/pxBeam)
        elif opt=='em':
            spec=ascii.read(specName)

            vel=spec['Velocity [km/s]']
            x_data=vel
            flux=spec['Flux [Jy]']
            y_data=flux*1e3

        if contFlux is not None:
            tau = self.optical_depth(flux,contFlux)
            y_data=-tau

        fig = plt.figure(figsize=(7.24409,4.074800625),constrained_layout=False)
        fig.set_tight_layout(False)
        
        gs = plt.GridSpec(nrows=1, ncols=1,  figure=fig,wspace=0.0,hspace=0.0)

        ax1 = fig.add_subplot(gs[0,0])

        # Initialize subplots
        #ax1 = fig.add_subplot(gs[0])
        ax1.set_xlabel('')
        ax1.set_ylabel('')
        ax1.tick_params(axis='both', bottom='on', top='on',
                        left='on', right='on', which='major', direction='in')
        ax1.tick_params(axis='both', bottom='on', top='on',
                        left='on', right='on', which='minor', direction='in')


        if contFlux is None and opt!='emSof':

            ylabh = ax1.set_ylabel(
                r'S\,$[\mathrm{mJy}\,\,\mathrm{beam}^{-1}]$')
        if contFlux is None and opt=='emSof':
            ylabh = ax1.set_ylabel(
                r'Flux density $[\mathrm{mJy}]$')
        else:
            ylabh = ax1.set_ylabel(
                r'$\tau$',labelpad=-10)            
            ylabh.set_verticalalignment('center')        



        if cfg_par['galaxy']['vsys']> 0:
            if vSys==None:
                ax1.axvline(cfg_par['galaxy']['vsys'],color='k',linestyle='-.',linewidth=0.8)
                ax1.set_xlabel(r'$cz\,[\mathrm{km}\,\mathrm{s}^{-1}]$')

            else:
                ax1.axvline(0,color='k',linestyle='-.',linewidth=0.8)
                x_data-=vSys
                ax1.set_xlabel(r'velocity $[\mathrm{km}\,\mathrm{s}^{-1}]$')


            # Plot spectra
    #                if self.abs_ex_plot_linestyle == 'step':
            # ax1.plot(x_data, y_data, color='black', linestyle='-')

            #if flag_chans != None:
            #    flag_chans = np.array(flag_chans)
            #    if cfg_par['spec_ex'].get('zunit') == 'm/s':
            #        flag_chans = np.divide(flag_chans, 1e3)
            #    index_flags_l = (np.abs(x_data - flag_chans[0])).argmin()
            #    for k in xrange(1, len(flag_chans)):
            #        index_flags = (np.abs(x_data - flag_chans[k])).argmin()
            #        # y_data[index_flags] = 0.0
            #    y_data[index_flags_l:index_flags] = 0.0
            
        ax1.step(x_data, y_data, where='mid', color='black', linestyle='-')
        # ax1.plot(x_data, y_data,  color='black', linestyle='-')

        if y_sigma is not None:
            ax1.fill_between(x_data, y_data-y_sigma, y_data+y_sigma,
                         facecolor='grey', alpha=0.4,step='mid', joinstyle="miter")
            # ax1.errorbar(x_data, y_data, elinewidth=0.75,
            #                 yerr=y_sigma, capsize=1, color='black', linestyle='-')
        # Calculate axis limits and aspect ratio
        if xLims==None:
            x_min = np.min(x_data)
            x_max = np.max(x_data)        
        else:
            x_min=xLims[0]
            x_max = xLims[1]
        y1_array = y_data[np.where((x_data > x_min) & (x_data < x_max))]
        #if cfg_par[key]['fixed_scale']:
        #    y1_min = -50
        #    y1_max = 50
        #else:

        if opt=='emSof':
            if yLims==None:
                peak = np.nanmax([np.nanmin(y_data),np.nanmax(y_data)])
                y1_min = (-peak/10.)*1.1
                #y1_max = np.nanmax(peak)*1.1
                y1_max = peak*1.1
            else:
                y1_min=yLims[0]
                y1_max = yLims[1]
           
            if source is not None:
                ax1.axvline(x_data[source['z_min']],color='grey', alpha=0.25,linestyle='-',linewidth=0.8,drawstyle='steps-pre')
                ax1.axvline(x_data[source['z_max']],color='grey', alpha=0.25,linestyle='-',linewidth=0.8,drawstyle='steps-post')


        else:
            peak = np.nanmax([-np.nanmin(y_data),np.nanmax(y_data)])
            y1_min = np.nanmin(-peak)*1.1
            #y1_max = np.nanmax(peak)*1.1
            y1_max = -y1_min/2.
            # Set axis limits
        ax1.set_ylim(y1_min, y1_max)

        ax1.set_xlim(x_min, x_max)
        ax1.xaxis.labelpad = 6
        ax1.yaxis.labelpad = 10
        
        ax1.axhline(color='k', linestyle='--', zorder=0)
   

        # Plot noise


        # Add minor tick marks
        ax1.minorticks_on()
        ax1.set_autoscale_on(False)    

        # Save figure to file
        # name of plot is combination of beam number and name of extracted source
        #outplot = os.path.basename(outPlot)
        # changed this to account for the different name convention
        # outplot = string.split(outplot,'.')[0]
        # outplot = cfg_par['general']['plotdir']+outplot+'.png'

        #outplot = "{0:s}{1:s}_{2:s}".format(
        #    cfg_par['general']['plotdir'], cfg_par['general']['workdir'].split("/")[-2], outplot.replace('.txt', '_compact.png'))
        #if cfg_par['general']['plot_format'] == "pdf":
        #plt.savefig(outplot.replace('.png', ".pdf"),
        #            overwrite=True, bbox_inches='tight')
        #else:
        outPlot = cfg_par['HIabs']['absDir']+cfg_par['HIabs']['title']+'.'+cfg_par['HIabs']['spectrum']['plotFormat']

        plt.savefig(outPlot,
                        overwrite=True, bbox_inches='tight', dpi=100)
        plt.show()
        #plt.close("all")

        return(outPlot)

    def absPlot(self,specName,outPlot,sVel,detPlot,yunit='flux',zunit='m/s'):
        '''
        Plots spectra of all radio sources found by find_src_imsad
        saved in basedir/beam/abs/spec.
        Plots are stored in basedir/beam/abs/plot
        IN
                Spectra extracted by spec_ex
        IN cfga
                abs_ex_plot_xaxis= ' '      #: X-axis units ['velocity','frequency']
                abs_ex_plot_yaxis= ' '      #: Y axis units ['flux','optical depth']
                #: plots line at redshift of source in spectrum redshift must be stored in table of load_src_csv
                abs_ex_plot_redsrc= True
                abs_ex_plot_title= True     #: plot title: J2000 name of radio source
                abs_ex_plot_format= ' '     #: format of plot ['.pdf','.jpeg','.png']
        OUT
                For each source outputs have the following name:
                J2000_xaxis-unit_yaxis-unit.plot_format = J220919.87+180920.17_vel_flux.pdf
        '''


        params = ut.loadRcParams()
        plt.rcParams.update(params)



        if os.path.isfile(specName) == True:



            spec_vec = ascii.read(specName)
            x_data = np.array(spec_vec[spec_vec.colnames[0]], dtype=float)
            n_channels = np.size(x_data)

            # Set plot specs
            font_size = 16
            plt.ioff()
            plt.rc('xtick', labelsize=font_size-2)
            plt.rc('ytick', labelsize=font_size-2)

            fig, ax1 = plt.subplots(figsize=(10, 5))
            #fig = plt.figure(figsize=(9, 6))
            # fig.subplots_adjust(hspace=0.0)
            #gs = gridspec.GridSpec(1, 1)

            # Initialize subplots
            #ax1 = fig.add_subplot(gs[0])
            ax1.set_xlabel('')
            ax1.set_ylabel('')
            ax1.tick_params(axis='both', bottom='on', top='on',
                            left='on', right='on', which='major', direction='in')
            ax1.tick_params(axis='both', bottom='on', top='on',
                            left='on', right='on', which='minor', direction='in')

            #flag_chans = cfg_par['spec_ex'].get('flag_chans', None)
            #flag_chans1 = cfg_par['spec_ex'].get('flag_chans', None)

            if zunit == 'm/s':
                x_data /= 1e3
            ax1.set_xlabel(r'$cz\,(\mathrm{km}\,\mathrm{s}^{-1})$', fontsize=font_size)
            if yunit=='flux':
                y_data = np.array(spec_vec[spec_vec.colnames[1]], dtype=float)*1e3
                y_sigma = np.array(spec_vec[spec_vec.colnames[2]])*1e3
                ylabh = ax1.set_ylabel(r'S\,$[\mathrm{mJy}\,\mathrm{beam}^{-1}]$', fontsize=font_size)
            elif yunit=='nhi':
                y_data = np.array(-spec_vec[spec_vec.colnames[5]], dtype=float)/1e18
                y_sigma = np.array(spec_vec[spec_vec.colnames[6]])/1e18
                
                ylabh = ax1.set_ylabel(r'N$_\mathrm{HI}$ [$\times 10^{18} \mathrm{cm}\,^{-2}$]', fontsize=font_size)
            elif yunit=='tau':

                y_data = np.array(-spec_vec[spec_vec.colnames[3]], dtype=float)
                y_sigma = np.array(spec_vec[spec_vec.colnames[4]])
                
                ylabh = ax1.set_ylabel(r'$\tau$', fontsize=font_size)
            #if cfg_par['spec_ex'].get('zunit') == 'MHz':
            #    x_data /= 1e6
            #    ax1.set_xlabel(r'Frequency [MHz]', fontsize=font_size)


            ylabh.set_verticalalignment('center')

            

            if sVel> 0:
                ax1.axvline(sVel,color='k',linestyle='-.',linewidth=1)
            # Plot spectra
            ax1.axvline(0,color='k',linestyle='-.',linewidth=1)


    #                if self.abs_ex_plot_linestyle == 'step':
            # ax1.plot(x_data, y_data, color='black', linestyle='-')

            #if flag_chans != None:
            #    flag_chans = np.array(flag_chans)
            #    if cfg_par['spec_ex'].get('zunit') == 'm/s':
            #        flag_chans = np.divide(flag_chans, 1e3)
            #    index_flags_l = (np.abs(x_data - flag_chans[0])).argmin()
            #    for k in xrange(1, len(flag_chans)):
            #        index_flags = (np.abs(x_data - flag_chans[k])).argmin()
            #        # y_data[index_flags] = 0.0
            #    y_data[index_flags_l:index_flags] = 0.0
            
            ax1.step(x_data-sVel, y_data, where='mid', color='black', linestyle='-')
            # Calculate axis limits and aspect ratio
            print(x_data,sVel)
            x_min = np.nanmin(x_data-sVel)
            x_max = np.nanmax(x_data-sVel)
            y1_array = y_data[np.where((x_data > x_min) & (x_data < x_max))]
            #if cfg_par[key]['fixed_scale']:
            #    y1_min = -50
            #    y1_max = 50
            #else:
            peak = np.nanmax([-np.nanmin(y_data),np.nanmax(y_data)])
            y1_min = np.nanmin(-peak)*1.1
            y1_max = np.nanmax(peak)*1.1

            # Set axis limits
            ax1.set_xlim(-500, 500)
            ax1.set_ylim(y1_min, 1)
            ax1.xaxis.labelpad = 6
            ax1.yaxis.labelpad = 10

            #if flag_chans1 != None:
            #    flag_chans = np.array(flag_chans1)
            #    if cfg_par['spec_ex'].get('zunit') == 'm/s':
            #        flag_chans = np.divide(flag_chans, 1e3)
            #    index_flags_l = (np.abs(x_data - flag_chans[0])).argmin()
            #    for k in xrange(1, len(flag_chans)):
            #        index_flags = (np.abs(x_data - flag_chans[k])).argmin()
            #        ax1.fill_between([x_data[index_flags_l], x_data[index_flags]],
            #                         y1_min, y1_max, facecolor='grey', alpha=0.3)

            # Plot noise
           #print y_sigma
            ax1.fill_between(x_data-sVel, -y_sigma, y_sigma,
                             facecolor='grey',step='mid',alpha=0.5)
            # Plot stuff
            ax1.axhline(color='k', linestyle=':', zorder=0)

            #redshifts = cfg_par[key].get('redshift_sources', None)
            #if len(redshifts) == 2:
            #    ax1.fill_between([redshifts[0], redshifts[1]], y1_min,
            #                     y1_max, facecolor='red', alpha=0.1)

            #if cfg_par[key]['title'] == True:
            #    ax1.set_title("{0:s} (\#{1:d}): {2:s}".format(cfg_par['general']['workdir'].split("/")[-2], int(os.path.basename(spec_name).split('_')[0]), os.path.basename(spec_name).replace(
            #        '.txt', '').split('_')[-1]), fontsize=font_size+2)
                # if self.abs_ex_plot_title == True:

            # Add title
            #print(outPlot)
            #aaa = outPlot.split('.png')[0]
            #print(aaa)
            #tt=aaa.split('/')[-1]
            #ax1.set_title("{0:s}".format(tt), fontsize=font_size+2)
            #ax1.axes.titlepad = 8

            # Add minor tick marks
            ax1.minorticks_on()

            # Save figure to file
            # name of plot is combination of beam number and name of extracted source
            #outplot = os.path.basename(outPlot)
            # changed this to account for the different name convention
            # outplot = string.split(outplot,'.')[0]
            # outplot = cfg_par['general']['plotdir']+outplot+'.png'

            #outplot = "{0:s}{1:s}_{2:s}".format(
            #    cfg_par['general']['plotdir'], cfg_par['general']['workdir'].split("/")[-2], outplot.replace('.txt', '_compact.png'))
            #if cfg_par['general']['plot_format'] == "pdf":
            #plt.savefig(outplot.replace('.png', ".pdf"),
            #            overwrite=True, bbox_inches='tight')
            #else:
            plt.savefig(outPlot,
                            overwrite=True, bbox_inches='tight', dpi=100)
            #plt.show()
            plt.close("all")

            # also create multi-plot spectra
            if detPlot > 0:
            #     # print(n_channels)

            #     # number of channels at which to split
                n_channel_per_plot = detPlot

            #     # get the number of plots
                n_plots = int(np.ceil(float(n_channels)/float(n_channel_per_plot)))

            #     # print(n_plots)

                n_rows = n_plots

            #     # add one row for the plot with full channel width
                fig, ax = plt.subplots(ncols=1, nrows=n_rows, figsize=(10, 2*n_rows),frameon=False,sharey=True)
                fig.subplots_adjust(hspace=0.2)
                fig.text(0.04, 0.5, r'S\,$[\rm{mJy}\,\rm{beam}^{-1}]$', fontsize=font_size, va='center', rotation='vertical')

            #      # ax1.annotate("Full spectrum", xy=(
            #      #     0.05, 0.95), xycoords='axes fraction', ha='left')

            #     # ax[1].annotate("Detailed spectrum", xy=(
            #     #     0.05, 0.95), xycoords='axes fraction', ha='left')
                
            #      if cfg_par[key]['title'] == True:
            #          ax[0].set_title("{0:s} (\#{1:d}): {2:s}".format(cfg_par['general']['workdir'].split("/")[-2], int(os.path.basename(spec_name).split('_')[0]), os.path.basename(spec_name).replace(
            #              '.txt', '').split('_')[-1]), fontsize=font_size+2)

            #     # go through the rest of the plots and create them
                for plot_count in range(n_rows):

            #         # the chunk of data corresponding to the plot
                    data_indices_min = plot_count * n_channel_per_plot
                    data_indices_max = (plot_count+1) * n_channel_per_plot
            #         # print(data_indices_min)
            #         # print(data_indices_max)
                    x_data_plot = x_data[data_indices_min:data_indices_max]
                    y_data_plot = y_data[data_indices_min:data_indices_max]
                    y_sigma_plot = y_sigma[data_indices_min:data_indices_max]*1e3
            

            #         # set the plot limits (only the x-axis needs to be adjusted)
                    ax[plot_count].set_ylim(y1_min, y1_max)
                    ax[plot_count].xaxis.labelpad = 6
                    ax[plot_count].yaxis.labelpad = 10
                    ax[plot_count].minorticks_on()
                    ax[plot_count].tick_params(axis='both', bottom='on', top='on',
                                               left='on', right='on', which='major', direction='in')
                    ax[plot_count].tick_params(axis='both', bottom='on', top='on',
                                               left='on', right='on', which='minor', direction='in')
                    
            #         # adjust the plot range of the last plot to match the others if the number
            #         # of channels cannot be divided by the number of channels per plot without rest
                    if plot_count == n_plots-1 and float(n_channels) % float(n_channel_per_plot) != 0:
            #             # get channel spacing assuming the spacing is linear
                        channel_width = np.diff(x_data_plot)[0]

            #             # number of missing channels
                        n_missing_channels = n_channel_per_plot - len(x_data_plot)

                        x_data_plot_min = np.min(
                            x_data_plot) + channel_width * n_missing_channels

            #             # fill array up to the number of channels
            #             # for k in range(n_missing_channels):
            #             #     x_data_plot = np.append(
            #             #         x_data_plot, x_data_plot[-1]+channel_width)
            #             #     y_data_plot = np.append(y_data_plot, 0.)
                    else:
                        x_data_plot_min = np.min(x_data_plot)

                    x_data_plot_max = np.max(x_data_plot)

                    if (sVel> 0 and sVel>x_data_plot_min and sVel< x_data_plot_max):
                        ax[plot_count].axvline(sVel,color='k',linestyle='-.',linewidth=1)
            #         #print(x_data_plot_max - x_data_plot_min)
                    ax[plot_count].step(x_data_plot, y_data_plot, where='mid', color='black', linestyle='-')
                    ax[plot_count].fill_between(x_data_plot, -y_sigma_plot, y_sigma_plot,facecolor='grey',step='mid',alpha=0.5)
                    #ax[plot_count].plot(x_data_plot, y_data_plot)

                    ax[plot_count].set_xlim(x_data_plot_min, x_data_plot_max)

            #         # ax[plot_count].fill_between(x_data_plot, -y_sigma_plot, y_sigma_plot,
            #         #                             facecolor='grey', alpha=0.5)
                    ax[plot_count].axhline(color='k', linestyle=':', zorder=0)

            #         # for the last plot add the x-axis label
                    if plot_count == n_plots-1:

            #            if cfg_par['spec_ex'].get('zunit') == 'm/s':
                            x_data /= 1e3
                            ax[plot_count].set_xlabel(r'$cz\,(\mathrm{km}\,\mathrm{s}^{-1})$', fontsize=font_size)
            #            elif cfg_par['spec_ex'].get('zunit') == 'MHz':
            #                 x_data /= 1e6
            #                 ax[plot_count].set_xlabel(
            #                     r'Frequency [MHz]', fontsize=font_size)

            #             # in case the last plot contains less channels than the others
            #             # if n_channels % n_channel_per_plot == 0:
            #             #     n_channels_subplot = np.size(x_data_plot)
            #             #     figwidth = fig.get_figwidth()
            #             #     ax[plot_count].set_figwidth(
            #             #         figwidth * n_channels_subplot/n_channel_per_plot)
                    if plot_count == 0:
                        ax[plot_count].set_title("{0:s}".format(tt+' - zoom'), fontsize=font_size+2)

            #     # name of plot
            #     outplot = os.path.basename(spec_name)
            #     outplot = "{0:s}{1:s}_{2:s}".format(
            #         cfg_par['general']['plotdir'], cfg_par['general']['workdir'].split("/")[-2], outplot.replace('.txt', '_detailed.png'))
            #     if cfg_par['general']['plot_format'] == "pdf":
            #         plt.savefig(outplot.replace('.png', ".pdf"),
            #                     overwrite=True, bbox_inches='tight')
            #     else:
                #ylabh = ax[plot_count].set_ylabel(
                #        r'S\,$[\mathrm{mJy}\,\mathrm{beam}^{-1}]$', fontsize=font_size)
                #ylabh.set_verticalalignment('center')

                outPlotDet = outPlot.replace('.png','_Zoom.png')
                plt.savefig(outPlotDet,overwrite=True, bbox_inches='tight', dpi=100)
                #plt.show()
                plt.close("all")
            #     plt.close("all")
            #if verb == True:
            #    print '# Plotted spectrum of source ' + os.path.basename(specName)+'. #'
        else:
            print('# Missing spectrum of source ' + os.path.basename(specName)+'. #')

 