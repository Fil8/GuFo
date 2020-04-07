#!/usr/bin/env python3.6
 
from astropy.io import fits
from astropy import wcs
import os
import numpy as np
from astropy.table import Table, Column

from ppxf.ppxf_util import log_rebin

from scavengers import cvPlay, tPlay
cvP = cvPlay.convert()
tP = tPlay.tplay()

class starsub(object):
    '''
    Modules to prepare for gaussian fitting of optical lines: 
    - subtract the stellar continuum from PPXF tables run with GISTPipeline ()
        Stellar subtractions is performed by subtracting the stellar and the nuclear emission from the original 
        unbinned data cube, by rescaling the fitted continuum – constant within each Voronoi bin – to the original 
        unbinned observed continuum flux of each spaxel.
    - Create a datacube of the stellar subtracted data by 
    - Perform Voronoi binning of the stellar subtracted datacube
    - Create the voronoi binned datacube
    '''

###these functions are likely not used
    # def run_logrebinning(self, bin_data, velscale, nbins, wave ):
    #     """
    #     Calls the log-rebinning routine of pPXF (see Cappellari & Emsellem 2004;
    #     ui.adsabs.harvard.edu/?#abs/2004PASP..116..138C;
    #     ui.adsabs.harvard.edu/?#abs/2017MNRAS.466..798C).
    #     """
    #     # Setup arrays
    #     lamRange = np.array([np.amin(wave),np.amax(wave)])
    #     sspNew, logLam, _ = log_rebin(lamRange, bin_data[:,0], velscale=velscale)
    #     log_bin_data = np.zeros([len(logLam),nbins])

    #     # Do log-rebinning 
    #     for i in range(0, nbins):
    #         log_bin_data[:,i] = self.corefunc_logrebin(lamRange, bin_data[:,i], velscale, len(logLam), i, nbins)

    #     return(log_bin_data)


    # def corefunc_logrebin(self, lamRange, bin_data, velscale, npix, iterate, nbins):
    #     """
    #     Calls the log-rebinning routine of pPXF (see Cappellari & Emsellem 2004;
    #     ui.adsabs.harvard.edu/?#abs/2004PASP..116..138C;
    #     ui.adsabs.harvard.edu/?#abs/2017MNRAS.466..798C). 

    #     TODO: Should probably be merged with run_logrebinning. 
    #     """
    #     try:
    #         sspNew, logLam, _ = log_rebin(lamRange, bin_data, velscale=velscale)
    #         pipeline.printProgress(iterate+1, nbins, barLength = 50)
    #         return(sspNew)

    #     except:
    #         out = np.zeros(npix); out[:] = np.nan
    #         return(out)


    def makeCubesVorBin(self,cfg_par):
        '''
        Creates voronoi binned data cubes from output of PPXF
        
        Parameters:
            - cfg_par: configuration file

        Uses:
            - output tables of PPXF run using the GISTpipeline 
            - GISTpipeline table of stellar continuum emission   
        Returns:
            voronoi binned datacube 
            voronoi binned noisecube
            star subtracted voronoi binned datacube
            voronoi binned datacube of stellar continuum emission

        Options:
            cfg_par['starSub']['outputs'] defines the kind of output datacubes
            all, stars, data, lines
        '''
        workDir = cfg_par['general']['workdir']

        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec,dataStar = tP.openPPXFforSubtraction(cfg_par,workDir+cfg_par['general']['tableBinName'],
            workDir+cfg_par['general']['tableSpecName'],workDir+cfg_par['general']['tableStarName'])

        data=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        Stars=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        Lines=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        noiseCube=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        
        header = self.makeHeader(cfg_par, wave, pxSize)

        hdu = fits.PrimaryHDU(data=data,header=header)

        diffusion = 1e-5
        del header['LATPOLE']
        del header['LONPOLE']
        xxVec = []
        yyVec = []

        for i in range(0,vorBinInfo['ID'].shape[0]):

            indexX = ((xAxis <= (np.round(vorBinInfo['X'][i],4)+diffusion)) & 
                      ((np.round(vorBinInfo['X'][i],4)-diffusion) < xAxis))
            
            indexY = ((yAxis <= (np.round(vorBinInfo['Y'][i],4)+diffusion)) & 
                      ((np.round(vorBinInfo['Y'][i],4)-diffusion) < yAxis))       
            
            xx = np.where(indexX)[0]
            yy = np.where(indexY)[0]

            xxVec.append(xx[0])
            yyVec.append(yy[0])
            indexBin =  vorBinInfo['BIN_ID'][i]
            
            if indexBin>0 and xx and yy: 
                   
                tmpD = np.array(dataSpec[indexBin][0][:])
                tmp = tmpD.tolist()
                data[:,yy[0],xx[0]] = tmp
                
                tmpN = np.array(dataSpec[indexBin][1][:])
                tmpN = tmpD/tmpN
                tmpN = tmpN.tolist()                
                noiseCube[:,yy[0],xx[0]] = tmpN

                tmp = np.array(dataStar[indexBin][1][:])
                tmp = tmp.tolist()
                Stars[:,yy[0],xx[0]] = tmp

            else:
                pass

        xxVecArr= Column(np.array(xxVec), name='PixX')
        yyVecArr= Column(np.array(yyVec), name='PixY')


        t = Table(vorBinInfo)
        t.add_column(xxVecArr,index=0)
        t.add_column(yyVecArr,index=0) 

        tab = fits.open(workDir+cfg_par['general']['tableBinName'])
        head = tab[0].header

        empty_primary = fits.PrimaryHDU(header=head)           

        #t2 = fits.BinTableHDU.from_columns(t,name='vorBinInfo')
        hdul = fits.HDUList([empty_primary])      

        hdul.append(fits.BinTableHDU(t.as_array(), name='vorBinInfo'))

        hdul.writeto(workDir+cfg_par['general']['outVorTableName'],overwrite=True)

        Lines = np.subtract(data,Stars)

        outputs = cfg_par['starSub']['outputs']
        
        if outputs == 'lines':
            fits.writeto(cfg_par['general']['outLines'],Lines,header,overwrite=True)
            fits.writeto(cfg_par['general']['outNoise'],noiseCube,header,overwrite=True)
            print('\n\t         +++\t\tLine Cube saved\t\t +++')
            return 
        else: 
            if outputs == 'stars':
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                print('\n\t         +++\t\tStar Cube saved\t\t +++')
                return 
            elif outputs == 'data':
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('\n\t         +++\t\tData Cube saved\t\t +++')
                return
            elif outputs=='all':
                fits.writeto(cfg_par['general']['outNoise'],noiseCube,header,overwrite=True)
                fits.writeto(cfg_par['general']['outLines'],Lines,header,overwrite=True)
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('\n\t         +++\t\tAll Cube saved\t\t +++')
                return

    def makeCubesPix(self,cfg_par):
        '''
        Creates data cubes from output of PPXF run using the GISTpipeline
        with the pixel size of the input data
        
        Parameters:
            - cfg_par: configuration file
        Uses:
            - output tables of voronoi binning 
                - cfg_par['general']['tableAllSpecName'], 
            - GISTpipeline table of stellar continuum emission   
            - module of tPlay: openPPXFforSubtraction
        Returns:
            - datacube 
            - noisecube
            - line cube (star subtracted datacube)
            - datacube of stellar continuum emission

        Options:
            cfg_par['starSub']['outputs'] defines the kind of output datacubes
            all, stars, data, lines

        Notes:
            the stellar emission is subtracted from the original unbinned datacube
            by rescaling the fitted continuum – constant within each Voronoi bin – to the original 
            unbinned observed continuum flux of each spaxel. 
        '''

        workDir = cfg_par['general']['workdir']

        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec,dataStar = tP.openPPXFforSubtraction(cfg_par,workDir+cfg_par['general']['tableBinName'],
            workDir+cfg_par['general']['tableAllSpecName'],workDir+cfg_par['general']['tableStarName'])
        

        lineInfo = tP.openLineList(cfg_par)

        dataSpec = np.array(dataSpec['SPEC'][:])
        
        dataSub=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        data=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])

        Stars=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        noiseCube=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        
        header = self.makeHeader(cfg_par, wave, pxSize)

        hdu = fits.PrimaryHDU(data=dataSub,header=header)

        f = fits.open(workDir+cfg_par['general']['inputCube'])
        hh = f[0].header
        dd = f[0].data


        waveCube = hh['CRVAL3']+(np.arange(dd.shape[0]))*hh['CD3_3']
        waveCube = waveCube / (1+cfg_par['general']['redshift'])     

        idx = np.where( np.logical_and( waveCube >= cfg_par['starSub']['waveMin'], waveCube <= cfg_par['starSub']['waveMax'] ) )[0]
   
        dd = dd[idx,:,:]


        waveCube = waveCube[idx]
        velscale  = (wave[1]-wave[0])*(cfg_par['general']['C']*1e-3)/np.mean(wave)
        diffusion = 1e-5
        del header['LATPOLE']
        del header['LONPOLE']
        xxVec = []
        yyVec = []
        
        #create AllSpectra datacube
        for i in range(0,vorBinInfo['ID'].shape[0]):

            indexX = ((xAxis <= (np.round(vorBinInfo['X'][i],4)+diffusion)) & 
                      ((np.round(vorBinInfo['X'][i],4)-diffusion) < xAxis))
            
            indexY = ((yAxis <= (np.round(vorBinInfo['Y'][i],4)+diffusion)) & 
                      ((np.round(vorBinInfo['Y'][i],4)-diffusion) < yAxis))       
            
            xx = np.where(indexX)[0]
            yy = np.where(indexY)[0]

            xxVec.append(xx[0])
            yyVec.append(yy[0])
            indexBin =  vorBinInfo['BIN_ID'][i]

            if indexBin>0 and xx and yy: 

                tmpS = np.array(dataStar[indexBin][1][:])
                tmpD = np.array(dataSpec[i,:])

                data[:,yy[0],xx[0]] = tmpD
                if cfg_par['starSub'].get('scaleHow',None) == 'mean':
                    starsScale=np.multiply(tmpS,(np.divide(np.nanmean(tmpD),np.nanmean(tmpS))))                
                elif cfg_par['starSub'].get('scaleHow',None) == 'median':
                    starsScale=np.multiply(tmpS,(np.divide(np.nanmedian(tmpD),np.nanmedian(tmpS))))

                dataSub[:,yy[0],xx[0]] = np.subtract(tmpD,starsScale)
                Stars[:,yy[0],xx[0]] = starsScale

            else:
                pass

        xxVecArr= Column(np.array(xxVec), name='PixX')
        yyVecArr= Column(np.array(yyVec), name='PixY')



        t = Table(vorBinInfo)
        t.add_column(xxVecArr,index=0)
        t.add_column(yyVecArr,index=0) 

        tab = fits.open(workDir+cfg_par['general']['tableBinName'])
        head = tab[0].header
   

        empty_primary = fits.PrimaryHDU(header=head)           

        hdul = fits.HDUList([empty_primary])      

        hdul.append(fits.BinTableHDU(t.as_array(), name='vorBinInfo'))


        hdul.writeto(cfg_par['general']['outVorTableName'],overwrite=True)
        outputs = cfg_par['starSub']['outputs']
        
        if outputs == 'lines':
            fits.writeto(cfg_par['general']['outLines'],dataSub,header,overwrite=True)
            print('\n\t         +++\t\tLine Cube saved\t\t +++')
            return 
        else: 
            if outputs == 'stars':
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                print('\n\t         +++\t\tStar Cube saved\t\t +++')    
                return 
            elif outputs == 'data':
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('\n\t         +++\t\tData Cube saved\t\t +++')                 
            elif outputs=='all':
                fits.writeto(cfg_par['general']['outNoise'],noiseCube,header,overwrite=True)
                fits.writeto(cfg_par['general']['outLines'],Lines,header,overwrite=True)
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('\n\t         +++\t\tAll Cube saved\t\t +++')    
                return

    def makeCubesVorLine(self,cfg_par):
        '''
        Creates line datacube from results of Voronoi binning
        
        Parameters:
            - cfg_par: configuration file
        Uses:
            - output tables of voronoi binning (created by vorPlay and stored in workdir/runName/tables) 
                - cfg_par['general']['outVorLineTableName'], cfg_par['general']['outVorSpectra']
            - module of tPlay: openVorLineOutput
        
        Returns:
            - voronoi binned datacube 
            - voronoi binned noisecube

        Notes:
            the output table of voronoi binning is updated with the pixel coordinates of each bin

        '''


        key = 'general'
        workDir = cfg_par['general']['workdir']
        
        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec = tP.openVorLineOutput(cfg_par,cfg_par['general']['outVorLineTableName'],
            cfg_par['general']['outVorSpectra'])

        data=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])*np.nan
        noiseCube=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])*np.nan
        
        header = self.makeHeader(cfg_par, wave, pxSize)

        hdu = fits.PrimaryHDU(data=data,header=header)

        diffusion = 1e-5
        del header['LATPOLE']
        del header['LONPOLE']
        xxVec = []
        yyVec = []
        #create AllSpectra datacube

        for i in range(0,vorBinInfo['ID'].shape[0]):
            #print xAxis
            #print yAxis
            indexBin =  np.abs(vorBinInfo['BIN_ID'][i])

            indexX = ((xAxis <= (np.round(vorBinInfo['X'][i],4)+diffusion)) & 
                      ((np.round(vorBinInfo['X'][i],4)-diffusion) < xAxis))
            
            indexY = ((yAxis <= (np.round(vorBinInfo['Y'][i],4)+diffusion)) & 
                      ((np.round(vorBinInfo['Y'][i],4)-diffusion) < yAxis))       
            
            xx = np.where(indexX)[0]
            yy = np.where(indexY)[0]
            #print(indexBin,xx,yy,vorBinInfo['X'][i],vorBinInfo['Y'][i],indexX)

            xxVec.append(xx[0])
            yyVec.append(yy[0])

            if  xx and yy and vorBinInfo['FLUX'][i]!=0.0:
                
                tmpD = np.array(dataSpec[indexBin][0][:])
                tmp = tmpD.tolist()
                data[:,yy[0],xx[0]] = tmp
                
                tmpN = np.array(dataSpec[indexBin][1][:])
                #tmpN = tmpD/tmpN
                tmpN = tmpN.tolist()                
                noiseCube[:,yy[0],xx[0]] = tmpN

                #tmp = np.array(dataStar[indexBin][1][:])
                #tmp = tmp.tolist()
                #Stars[:,yy[0],xx[0]] = tmp

            else:
                pass
        
        #idx = np.where(data==0.0)[0]
        #data[idx] = np.nan
        xxVecArr= Column(np.array(xxVec), name='PixX')
        yyVecArr= Column(np.array(yyVec), name='PixY')

        #yyVecArr=np.array(yyVec,dtype={'names':('PixY')})
        #print(vorBinInfo.shape)
        #print(xxVecArr.shape)

        t = Table(vorBinInfo)
        t.add_column(xxVecArr,index=0)
        t.add_column(yyVecArr,index=0) 

        tab = fits.open(workDir+cfg_par['general']['tableBinName'])
        head = tab[0].header
  

        empty_primary = fits.PrimaryHDU(header=head)           

        hdul = fits.HDUList([empty_primary])      

        hdul.append(fits.BinTableHDU(t.as_array(), name='vorBinInfo'))


        hdul.writeto(cfg_par['general']['outVorLineTableName'],overwrite=True)       

        fits.writeto(cfg_par['general']['outVorLines'],data,header,overwrite=True)
        fits.writeto(cfg_par['general']['outVorNoise'],noiseCube,header,overwrite=True)
        print('''\t+---------+\n\t Line Cube saved\n\t+---------+''')     
        return 
    


    def makeHeader(self,cfg_par,wave,pxSize):
        '''
        Defines header of output datacube
        Puts wcs coordinates in datacube from information provided in the configuration file
        
        Parameters:
            - cfg_par: configuration file
            - wave: x-axis of spectra in log(lamdba) units
            - pxSize: pixel size of output cubes
        
        Returns:
            - header of datacubes
        '''       
        crPix3 =cfg_par['starSub']['pixZ']
        crPix1=cfg_par['starSub']['pixX']
        crPix2=cfg_par['starSub']['pixY']

        crVal3 = np.min(wave)
        deltaLambda = (np.max(wave)-np.min(wave))/len(wave)

        crpix3 = 1
        crVal1=cfg_par['starSub']['ra']
        crVal2=cfg_par['starSub']['dec']

        ra = cvP.hms2deg(crVal1)
        dec = cvP.dms2deg(crVal2)

        w = wcs.WCS(naxis=3)

        w.wcs.crpix = [crPix1, crPix2, crPix3]
        w.wcs.cdelt = np.array([-pxSize,pxSize,deltaLambda])
        w.wcs.crval = [ra, dec, crVal3]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN","AWAV"]

        header = w.to_header()

        header['EQUINOX'] = 2000
        header['CRDER3'] = 0.026
        header['CTYPE3'] = "AWAV"
        header['CUNIT3'] = "Angstrom"

        return header
