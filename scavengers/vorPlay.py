#!/usr/bin/env python3.6
import os, sys
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from vorbin.voronoi_2d_binning import voronoi_2d_binning
import functools

from scavengers import tPlay, cvPlay, starSub

tP = tPlay.tplay()
cvP = cvPlay.convert()
ss = starSub.starsub()

class vorplay(object):
    """
    COPYRIGHT
        This class is deliberately taken from the GIST pipeline: 
        A multi-purpose tool for the analysis and visualisation of (integral-field) spectroscopic data
        (abittner.gitlab.io/thegistpipeline/index.html)
        (ui.adsabs.harvard.edu/abs/2019A%26A...628A.117B/abstract)
    PURPOSE: 
        This file contains a collection of functions necessary to Voronoi-bin the data. 
        The Voronoi-binning makes use of the algorithm from Cappellari & Copin 2003
        (ui.adsabs.harvard.edu/?#abs/2003MNRAS.342..345C). 
    """
    def vorBin(self,cfg_par):

        
        key = 'general'

        workDir = cfg_par[key]['workdir']
        
        #open line lineList
        lineInfo = tP.openLineList(cfg_par)
        #diffusion = 1e-5
        
        #open table for bins
        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo  = tP.openTablesPPXF(cfg_par,workDir+cfg_par[key]['outVorTableName'],
            workDir+cfg_par[key]['tableSpecName'])
        
        cvel      = 299792.458
        velscale  = (wave[1]-wave[0])*cvel/np.mean(wave)
        #sys.exit(0)
        #open datacube
        f = fits.open(cfg_par['general']['cubeDir']+cfg_par['general']['dataCubeName'])
        hh = f[0].header
        dd = f[0].data
        #s     = np.shape(dd)
        #spec  = np.reshape(dd,[s[0],s[1]*s[2]])
        
        #estimate the noise around the Hbeta line
        velRangeInf = cvP.vRadLambda(15000,
                lineInfo['Wave'][0])-lineInfo['Wave'][0] 
        velRangeSup = cvP.vRadLambda(4000,
                lineInfo['Wave'][0])-lineInfo['Wave'][0] 

        waveLeftSup = np.log(lineInfo['Wave'][0]-velRangeSup)
        #print(waveLeftSup, lineInfo['Wave'][0], velRangeSup)
        idxWaveLeftSup = int(np.where(abs(wave-waveLeftSup)==abs(wave-waveLeftSup).min())[0])

        waveLeftInf = np.log(lineInfo['Wave'][0]-velRangeInf)
        idxWaveLeftInf = int(np.where(abs(wave-waveLeftInf)==abs(wave-waveLeftInf).min())[0])

        waveRightSup = np.log(lineInfo['Wave'][0]+velRangeInf)
        idxWaveRightSup = int(np.where(abs(wave-waveRightSup)==abs(wave-waveRightSup).min())[0])

        waveRightInf = np.log(lineInfo['Wave'][0]+velRangeSup)
        idxWaveRightInf = int(np.where(abs(wave-waveRightInf)==abs(wave-waveRightInf).min())[0])

        #find peak of HbetaLine
        waveAmpIn1Min = np.log(lineInfo['Wave'][0]-lineInfo['cenRangeAng'][0])
        indexMin = int(np.where(abs(wave-waveAmpIn1Min)==abs(wave-waveAmpIn1Min).min())[0]) 

        waveAmpIn1Max = np.log(lineInfo['Wave'][0]+lineInfo['cenRangeAng'][0])
        indexMax = int(np.where(abs(wave-waveAmpIn1Max)==abs(wave-waveAmpIn1Max).min())[0])
        
        peak = np.empty([dd.shape[1],dd.shape[2]])
        stdLeft = np.empty([dd.shape[1],dd.shape[2]])
        stdRight = np.empty([dd.shape[1],dd.shape[2]])
        noise = np.empty([dd.shape[1],dd.shape[2]])

        for j in range(0,dd.shape[1]):
            for i in range(0,dd.shape[2]):
                peak[j,i] = np.nanmax(dd[indexMin:indexMax,j,i])
                stdLeft[j,i] = np.nanstd(dd[idxWaveLeftInf:idxWaveLeftSup,j,i])
                stdRight[j,i] = np.nanstd(dd[idxWaveRightInf:idxWaveRightSup,j,i])
                noise[j,i] = np.divide(np.nansum([stdLeft[j,i], stdRight[j,i]]),2.)       

        snr = np.divide(peak,noise)
        #

        snr = np.reshape(snr,[dd.shape[1]*dd.shape[2]])
        noise = np.reshape(noise,[dd.shape[1]*dd.shape[2]])
        spec = np.reshape(dd[indexMin:indexMax],[indexMax-indexMin,dd.shape[1]*dd.shape[2]])
        specFull = np.reshape(dd,[dd.shape[0],dd.shape[1]*dd.shape[2]])
        x, y  = np.meshgrid(xAxis,yAxis)
        x     = np.reshape(x,[dd.shape[1]*dd.shape[2]])
        y     = np.reshape(y,[dd.shape[1]*dd.shape[2]])

        # Removing obviously defective pixels: Remove spaxel with any nan or negative values
        idx_good = np.where( np.median(spec, axis=0) > 0.0 )[0]
        spec     = spec[:,idx_good]
        specFull     = specFull[:,idx_good]

        noise    = noise[idx_good]
        x        = x[idx_good]
        y        = y[idx_good]
        signal   = np.nanmax(spec,axis=0)
        
        #sys.exit(0)       
        #noise = np.nanmean([np.sum(np.nanstd(spec[idxWaveLeftInf:idxWaveLeftSup,:]),np.nanstd(spec[idxWaveRightInf:idxWaveRightSup,:]))])
        #estimate the errors with the der_snr algorithm
        #espec = np.zeros( spec.shape )
        #for i in range( 0, spec.shape[1] ):
        #    espec[:,i] = der_snr.der_snr( spec[:,i] )
                
        # Storing eveything into a structure
        #cube = {'x':xAxis, 'y':yAxis, 'wave':wave, 'spec':spec, 'error':espec, 'snr':snr,
        #'signal':signal, 'noise':noise, 'velscale':velscale, 'pixelsize':pxSize}    
        #print(np.nanmean(signal),np.nanmin(signal),np.nanmin(noise))
        binNum = self.define_voronoi_bins(cfg_par, x, y, signal,noise, pxSize,
            snr, cfg_par['vorBin']['snr'], cfg_par['vorBin']['covarNoise'])

        self.apply_voronoi_bins(binNum, specFull, noise,velscale, wave, 'lin')

        ss.makeCubesVorLine(gPar.cfg_par)


        return

    def sn_func(self,index, signal=None, noise=None, covar_vor=0.00 ):
        """
           This function is passed to the Voronoi binning routine of Cappellari &
           Copin 2003 (ui.adsabs.harvard.edu/?#abs/2003MNRAS.342..345C) and used to
           estimate the noise in the bin from the noise in the spaxels. This
           implementation is identical to the default one, but accounts for spatial
           correlations in the noise by applying an empirical equation (see e.g.
           Garcia-Benito et al. 2015;
           ui.adsabs.harvard.edu/?#abs/2015A&A...576A.135G) together with the
           parameter defined in the Config-file. 
        """

        # Add the noise in the spaxels to obtain the noise in the bin
        sn = np.sum(signal[index])/np.sqrt(np.sum(noise[index]**2))

        # Account for spatial correlations in the noise by applying an empirical
        # equation (see e.g. Garcia-Benito et al. 2015;
        # ui.adsabs.harvard.edu/?#abs/2015A&A...576A.135G)
        sn /= 1 + covar_vor * np.log10(index.size)

        return(sn)


    def define_voronoi_bins(self, cfg_par, x, y, signal, noise, pixelsize, snr, target_snr, covar_vor):
        """
        This function applies the Voronoi-binning algorithm of Cappellari & Copin
        2003 (ui.adsabs.harvard.edu/?#abs/2003MNRAS.342..345C) to the data. It can
        be accounted for spatial correlations in the noise (see function sn_func()).
        A BIN_ID is assigned to every spaxel. Spaxels which do not satisfy the
        minimum SNR threshold are excluded from the Voronoi-binning, but are
        assigned a negative BIN_ID, with the absolute value of the BIN_ID
        corresponding to the nearest Voronoi-bin that satisfies the minimum SNR
        threshold.  All results are saved in a dedicated table to provide easy means
        of matching spaxels and bins. 
        """
        # Pass a function for the SNR calculation to the Voronoi-binning algorithm,
        # in order to account for spatial correlations in the noise
        sn_func_covariances = functools.partial(self.sn_func, covar_vor=covar_vor )

        try:
            print('\n\t *********** --- GuFo: VorBinning --- ***********\n')
            # Do the Voronoi binning
            binNum, xNode, yNode, xBar, yBar, sn, nPixels, _ = voronoi_2d_binning(x, y,
                signal, noise, target_snr, plot=False, quiet=True, pixelsize=pixelsize,
                sn_func=sn_func_covariances )
        
            # Handle common exceptions
        except ValueError as e:

            # Sufficient SNR and no binning needed
            if str(e) == 'All pixels have enough S/N and binning is not needed': 

                print("             "+"The Voronoi-binning routine of Cappellari & Copin (2003) returned the following error:")
                print("             "+str(e)+"\n")
                #print("             "+str(len(idx_inside))+" spaxels will be treated as Voronoi-bins.")
                
                logging.warning("Defining the Voronoi bins failed. The Voronoi-binning routine of Cappellari & Copin "+\
                                "(2003) returned the following error: \n"+str(e))
                #logging.info("Analysis will continue without Voronoi-binning! "+str(len(idx_inside))+" spaxels will be treated as Voronoi-bins.")

                binNum, xNode, yNode, sn, nPixels = self.noBinning( x, y, snr, idx_inside )

            # Any uncaught exceptions
            else:
                print("The Voronoi-binning routine of Cappellari & Copin (2003) returned the following error: \n"+str(e))
                return( True )

        # Find the nearest Voronoi bin for the pixels outside the Voronoi region
        #binNum_outside = find_nearest_voronoibin( x, y, idx_outside, xNode, yNode )

        # Generate extended binNum-list: 
        #   Positive binNum indicate the Voronoi bin of the spaxel (for spaxels inside the Voronoi region)
        #   Negative binNum indicate the nearest Voronoi bin of the spaxel (for spaxels outside of the Voronoi region)
        #ubins = np.unique(binNum)
        #nbins = len(ubins)
        #binNum_long = np.zeros( len(x) ); binNum_long[:] = np.nan
        #binNum_long[idx_inside]  = binNum
        #binNum_long[idx_outside] = -1 * binNum_outside

        # Save bintable: data for *ALL* spectra inside and outside of the Voronoi region!
        self.save_table(cfg_par, x, y, signal, snr, binNum, np.unique(binNum), xNode, yNode, sn, nPixels, pixelsize)

        return(binNum)

    def find_nearest_voronoibin(self,x, y, idx_outside, xNode, yNode):
        """
        This function determines the nearest Voronoi-bin for all spaxels which do
        not satisfy the minimum SNR threshold. 
        """
        x = x[idx_outside]
        y = y[idx_outside]
        pix_coords = np.concatenate( (x.reshape((len(x),1)),         y.reshape((len(y),1))),         axis=1 )
        bin_coords = np.concatenate( (xNode.reshape((len(xNode),1)), yNode.reshape((len(yNode),1))), axis=1 )

        dists = dist.cdist( pix_coords, bin_coords, 'euclidean' ) 
        closest = np.argmin( dists, axis=1 )

        return(closest)

    def noBinning(self, x, y, snr, idx_inside):
        """ 
        In case no Voronoi-binning is required/possible, treat spaxels in the input
        data as Voronoi bins, in order to continue the analysis. 
        """
        binNum  = np.arange( 0, len(idx_inside) )
        xNode   = x[idx_inside]
        yNode   = y[idx_inside]
        sn      = snr[idx_inside]
        nPixels = np.ones( len(idx_inside) )

        return( binNum, xNode, yNode, sn, nPixels )

    def save_table(self,cfg_par, x, y, signal, snr, binNum_new, ubins, xNode, yNode, sn, nPixels, pixelsize):
        """ 
        Save all relevant information about the Voronoi binning to disk. In
        particular, this allows to later match spaxels and their corresponding bins. 
        """
        outfits_table = cfg_par['general']['workdir']+cfg_par['general']['outVorLineName']

        # Expand data to spaxel level
        xNode_new = np.zeros( len(x) )
        yNode_new = np.zeros( len(x) )
        sn_new = np.zeros( len(x) )
        nPixels_new = np.zeros( len(x) )
        for i in range( len(ubins) ):
            idx = np.where( ubins[i] == np.abs(binNum_new) )[0]
            xNode_new[idx] = xNode[i]
            yNode_new[idx] = yNode[i]
            sn_new[idx] = sn[i]
            nPixels_new[idx] = nPixels[i]

        cols = []
        print(np.arange(len(x)),binNum_new)
        cols.append(fits.Column(name='ID',        format='J',   array=np.arange(len(x)) ))
        cols.append(fits.Column(name='BIN_ID',    format='J',   array=binNum_new        ))
        cols.append(fits.Column(name='X',         format='D',   array=x                 ))
        cols.append(fits.Column(name='Y',         format='D',   array=y                 ))
        cols.append(fits.Column(name='FLUX',      format='D',   array=signal            ))
        cols.append(fits.Column(name='SNR',       format='D',   array=snr               ))
        cols.append(fits.Column(name='XBIN',      format='D',   array=xNode_new         ))
        cols.append(fits.Column(name='YBIN',      format='D',   array=yNode_new         ))
        cols.append(fits.Column(name='SNRBIN',    format='D',   array=sn_new            ))
        cols.append(fits.Column(name='NSPAX',     format='J',   array=nPixels_new       ))
        tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        tbhdu.writeto(outfits_table, overwrite=True)
        fits.setval(outfits_table, "PIXSIZE", value=pixelsize)

        print("Writing: "+outfits_table)


    def apply_voronoi_bins(self,cfg_par, binNum, spec, espec, velscale, wave):
        """
        The constructed Voronoi-binning is applied to the underlying spectra. The
        resulting Voronoi-binned spectra are saved to disk. 
        """
        # Apply Voronoi bins
        print("Applying the Voronoi bins to data")
        idx_inside = np.where( binNum == ubins[i] )[0]
        bin_data, bin_error, bin_flux = self.voronoi_binning( binNum, spec, espec )
        # Save Voronoi binned spectra
        self.save_vorspectra(cfg_par, bin_data, bin_error, velscale, wave)
        return(None)


    def voronoi_binning(self, binNum, spec, error ):
        """ Spectra belonging to the same Voronoi-bin are added. """
        ubins     = np.unique(binNum[1])
        nbins     = len(ubins)
        npix      = spec.shape[0]
        print('Npix,NBins')
        print(npix,nbins)
        print(binNum.shape)
        bin_data  = np.zeros([npix,nbins])
        bin_error = np.zeros([npix,nbins])
        bin_flux  = np.zeros(nbins)

        for i in range(nbins):
            k = np.where( binNum == ubins[i] )[0]
            valbin = len(k)
            if valbin == 1:
               av_spec     = spec[:,k]
               av_err_spec = error[:,k]
            else:
               av_spec     = np.nansum(spec[:,k],axis=1)
               av_err_spec = np.sqrt(np.sum(error[:,k],axis=1))
        
            bin_data[:,i]  = np.ravel(av_spec)
            bin_error[:,i] = np.ravel(av_err_spec)
            bin_flux[i]    = np.mean(av_spec,axis=0)

        return(bin_data, bin_error, bin_flux)

    def save_vorspectra(self, cfg_par, log_spec, log_error, velscale, logLam):
        """ Voronoi-binned spectra and error spectra are saved to disk. """
        #if flag == 'log':
        #    outfits_spectra  = outdir+rootname+'_LineVorSpectra.fits'
        #elif flag == 'lin':
        outfits_spectra  = cfg_par['general']['outVorSpectra']

        npix = len(log_spec)

        # Create primary HDU
        priHDU = fits.PrimaryHDU()

        # Table HDU for spectra
        cols = []
        cols.append( fits.Column(name='SPEC',  format=str(npix)+'D', array=log_spec.T  ))
        cols.append( fits.Column(name='ESPEC', format=str(npix)+'D', array=log_error.T ))
        dataHDU = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        dataHDU.name = 'VOR_SPECTRA'

        # Table HDU for LOGLAM
        cols = []
        cols.append( fits.Column(name='LOGLAM', format='D', array=logLam ))
        loglamHDU = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        loglamHDU.name = 'LOGLAM'

        # Create HDU list and save to file
        HDUList = fits.HDUList([priHDU, dataHDU, loglamHDU])
        HDUList.writeto(outfits_spectra, overwrite=True)

        # Set header values
        fits.setval(outfits_spectra,'VELSCALE',value=velscale)
        fits.setval(outfits_spectra,'CRPIX1',  value=1.0)
        fits.setval(outfits_spectra,'CRVAL1',  value=logLam[0])
        fits.setval(outfits_spectra,'CDELT1',  value=logLam[1]-logLam[0])

        #if flag == 'log':
        #    pipeline.prettyOutput_Done("Writing: "+rootname+'_VorSpectra.fits')
        #elif flag == 'lin':
        #    pipeline.prettyOutput_Done("Writing: "+rootname+'_VorSpectra_linear.fits')
        #logging.info("Wrote: "+outfits_spectra)
