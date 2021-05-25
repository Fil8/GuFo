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

        #pxSize*=3600.

        header = self.makeHeader(cfg_par, wave, pxSize)
        cvel      = 299792.458
        velscale  = (wave[1]-wave[0])*cvel/np.mean(wave)

        lineInfo = tP.openLineList(cfg_par)

        dataSpec = np.array(dataSpec['SPEC'][:])
        
        dataSub=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        data=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])

        Stars=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        noiseCube=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        

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
        noiseVec = []
        specVec = []
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
                tmpN = np.array(noiseBin[i,:])

                data[:,yy[0],xx[0]] = tmpD
                if cfg_par['starSub'].get('scaleHow',None) == 'mean':
                    starsScale=np.multiply(tmpS,(np.divide(np.nanmean(tmpD),np.nanmean(tmpS))))                
                elif cfg_par['starSub'].get('scaleHow',None) == 'median':
                    starsScale=np.multiply(tmpS,(np.divide(np.nanmedian(tmpD),np.nanmedian(tmpS))))

                dataSub[:,yy[0],xx[0]] = np.subtract(tmpD,starsScale)
                Stars[:,yy[0],xx[0]] = starsScale

                noiseVec.append([tmpN])
                specVec.append([np.subtract(tmpD,starsScale)])

            else:
                noiseVec.append([np.zeros(len(noiseBin[i,:]))])
                specVec.append([np.zeros(len(dataSpec[i,:]))])



        xxVecArr= Column(np.array(xxVec), name='PixX')
        yyVecArr= Column(np.array(yyVec), name='PixY')


        if cfg_par['gFit']['method'] == 'pixel':
            tab = fits.open(workDir+cfg_par['general']['tableAllSpecName'])
            priHDU = fits.PrimaryHDU()
            # Table HDU for spectra
            cols = []

            #print(len(specVec),len(noiseVec))
            #sys.exit(0)
            cols.append( fits.Column(name='SPEC',  format=str(len(np.array(specVec).T))+'D', array=np.array(specVec) ))
            cols.append( fits.Column(name='ESPEC', format=str(len(np.array(specVec).T))+'D', array=np.array(noiseVec) ))
            #print(cols)
            print(len(cols))
            #sys.exit(0)
            dataHDU = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
            dataHDU.name = 'SPECTRA'
            
            #loglamHDU = tab[2].data
            #oglamHDU.name = 'LOGLAM'

            hdl = fits.HDUList([priHDU,dataHDU,tab[2]])

            hdl.writeto(cfg_par['general']['outPixSpectra'],overwrite=True)

            fits.setval(cfg_par['general']['outPixSpectra'],'VELSCALE',value=velscale)
            fits.setval(cfg_par['general']['outPixSpectra'],'CRPIX1',  value=1.0)
            fits.setval(cfg_par['general']['outPixSpectra'],'CRVAL1',  value=tab[2].data[0][0])
            fits.setval(cfg_par['general']['outPixSpectra'],'CDELT1',  value=tab[2].data[1][0]-tab[2].data[0][0])


        # Table HDU for spectra
        # cols = []
        # npix = len(xxVec)*len(yyVec)
        # cols.append( fits.Column(name='SPEC',  format=str(npix)+'D', array=log_spec.T  ))
        # cols.append( fits.Column(name='ESPEC', format=str(npix)+'D', array=log_error.T ))
        # dataHDU = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        # dataHDU.name = 'VOR_SPECTRA'


        t = Table(vorBinInfo)
        t.add_column(xxVecArr,index=0)
        t.add_column(yyVecArr,index=0) 
        if cfg_par['gFit']['method'] == 'pixel':
            t['BIN_ID'] = t['ID']        

        tab = fits.open(workDir+cfg_par['general']['tableBinName'])
        head = tab[0].header
   

        empty_primary = fits.PrimaryHDU(header=head)           

        hdul = fits.HDUList([empty_primary])      

        hdul.append(fits.BinTableHDU(t.as_array(), name='vorBinInfo'))


        hdul.writeto(cfg_par['general']['outVorTableName'],overwrite=True)
        outputs = cfg_par['starSub']['outputs']
        
        if outputs == 'lines':
            fits.writeto(cfg_par['general']['outLines'],dataSub,header,overwrite=True)
            print(header)
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
                fits.writeto(cfg_par['general']['outLines'],dataSub,header,overwrite=True)
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

    def makeHeaderMoms(self,cfg_par,pxSize):
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
        crPix1=cfg_par['starSub']['pixX']
        crPix2=cfg_par['starSub']['pixY']

        crVal1=cfg_par['starSub']['ra']
        crVal2=cfg_par['starSub']['dec']

        ra = cvP.hms2deg(crVal1)
        dec = cvP.dms2deg(crVal2)

        w = wcs.WCS(naxis=2)

        w.wcs.crpix = [crPix1, crPix2]
        w.wcs.cdelt = np.array([-pxSize,pxSize])
        w.wcs.crval = [ra, dec]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        header = w.to_header()

        header['EQUINOX'] = 2000

        return header

    def makePPFXmaps(self,cfg_par):
        '''
        Denerates the moment maps of stars from the PPFX output
        
        Parameters:
            - cfg_par: configuration file
               - tableBinName = _table of ppxf output
               - tableStarKin = _ppxf of ppxf output
        
        Returns:
            - output names of [flux,v,sigma] maps
        '''               

        print(cfg_par)
        crPix1=cfg_par['starSub']['pixX']
        crPix2=cfg_par['starSub']['pixY']
      
        tab = fits.open(cfg_par['general']['workdir']+cfg_par['general']['tableBinName'])
        head = tab[0].header
        headTab = tab[1].header
        dataTab = tab[1].data    

        tab = fits.open(cfg_par['general']['workdir']+cfg_par['general']['tableStarKin'])
        dataStar = tab[1].data    

        head['CRPIX1'] = crPix1
        head['CRPIX2'] = crPix2 
        
        xMin = np.min(dataTab['X'])
        xMax = np.max(dataTab['X'])

        shapeX = (xMax-xMin)/head['PIXSIZE']

        yMin = np.min(dataTab['Y'])
        yMax = np.max(dataTab['Y'])

        shapeY = (yMax-yMin)/head['PIXSIZE']

        xAxis = np.arange(xMin, xMax,head['PIXSIZE'])
        yAxis = np.arange(yMin, yMax,head['PIXSIZE'])
        
        pxSize = head['PIXSIZE']/3600.

        head = self.makeHeaderMoms(cfg_par,pxSize)

        diffusion = 1e-5
        del head['LATPOLE']
        del head['LONPOLE']
        xxVec = []
        yyVec = []
        #create AllSpectra datacube
        mom0=np.empty([yAxis.shape[0],xAxis.shape[0]])*np.nan
        mom1=np.empty([yAxis.shape[0],xAxis.shape[0]])*np.nan
        mom2=np.empty([yAxis.shape[0],xAxis.shape[0]])*np.nan

        for i in range(0,dataTab['ID'].shape[0]):
            #print xAxis
            #print yAxis
            match_bin = np.where(dataStar['BIN_ID']==(dataTab['BIN_ID'][i]))[0]
            
            indexX = ((xAxis <= (np.round(dataTab['X'][i],4)+diffusion)) & 
                      ((np.round(dataTab['X'][i],4)-diffusion) < xAxis))
            
            indexY = ((yAxis <= (np.round(dataTab['Y'][i],4)+diffusion)) & 
                      ((np.round(dataTab['Y'][i],4)-diffusion) < yAxis))       
            
            xx = np.where(indexX)[0]
            yy = np.where(indexY)[0]
            #print(indexBin,xx,yy,vorBinInfo['X'][i],vorBinInfo['Y'][i],indexX)

            xxVec.append(xx[0])
            yyVec.append(yy[0])

            if  xx and yy and dataTab['FLUX'][i]!=0.0 and len(match_bin)==1:
                
                    mom0[yy[0],xx[0]] = dataTab['Flux'][i]
                                    
                    mom1[yy[0],xx[0]] = dataStar['V'][match_bin]
                    mom2[yy[0],xx[0]] = dataStar['SIGMA'][match_bin]
                

                #tmp = np.array(dataStar[indexBin][1][:])
                #tmp = tmp.tolist()
                #Stars[:,yy[0],xx[0]] = tmp

            else:
                pass        

        print(np.nanmedian(mom1))

        mom1-=np.nanmedian(mom1)
        mom1+=2035
        head['BUNIT']='10-20 erg s^{-1} cm^{-2} AA^{-1}'
        outMom0=cfg_par['general']['momDir']+'mom0Stars.fits'
        fits.writeto(outMom0,mom0,head,overwrite=True)
        head['BUNIT'] = 'km/s'
        head['SPECSYS'] = 'topocent'
        outMom1=cfg_par['general']['momDir']+'mom1Stars.fits'
        fits.writeto(cfg_par['general']['momDir']+'mom1Stars.fits',mom1,head,overwrite=True)
        outMom2=cfg_par['general']['momDir']+'mom2Stars.fits'
        fits.writeto(cfg_par['general']['momDir']+'mom2Stars.fits',mom2,head,overwrite=True)

        return [mom0,mom1]
