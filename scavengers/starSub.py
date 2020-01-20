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

    def run_logrebinning(self, bin_data, velscale, nbins, wave ):
        """
        Calls the log-rebinning routine of pPXF (see Cappellari & Emsellem 2004;
        ui.adsabs.harvard.edu/?#abs/2004PASP..116..138C;
        ui.adsabs.harvard.edu/?#abs/2017MNRAS.466..798C).
        """
        # Setup arrays
        lamRange = np.array([np.amin(wave),np.amax(wave)])
        sspNew, logLam, _ = log_rebin(lamRange, bin_data[:,0], velscale=velscale)
        log_bin_data = np.zeros([len(logLam),nbins])

        # Do log-rebinning 
        for i in range(0, nbins):
            log_bin_data[:,i] = self.corefunc_logrebin(lamRange, bin_data[:,i], velscale, len(logLam), i, nbins)

        return(log_bin_data)


    def corefunc_logrebin(self, lamRange, bin_data, velscale, npix, iterate, nbins):
        """
        Calls the log-rebinning routine of pPXF (see Cappellari & Emsellem 2004;
        ui.adsabs.harvard.edu/?#abs/2004PASP..116..138C;
        ui.adsabs.harvard.edu/?#abs/2017MNRAS.466..798C). 

        TODO: Should probably be merged with run_logrebinning. 
        """
        try:
            sspNew, logLam, _ = log_rebin(lamRange, bin_data, velscale=velscale)
            pipeline.printProgress(iterate+1, nbins, barLength = 50)
            return(sspNew)

        except:
            out = np.zeros(npix); out[:] = np.nan
            return(out)


    def makeCubesVorBin(self,cfg_par):

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
        #create AllSpectra datacube
        for i in range(0,vorBinInfo['ID'].shape[0]):
            #print xAxis
            #print yAxis
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

        #yyVecArr=np.array(yyVec,dtype={'names':('PixY')})
        #print(vorBinInfo.shape)
        #print(xxVecArr.shape)

        t = Table(vorBinInfo)
        t.add_column(xxVecArr,index=0)
        t.add_column(yyVecArr,index=0) 
        #vorBinInfo = np.column_stack((vorBinInfo,xxVec,yyVec))
        #vorBinInfo = np.vstack([vorBinInfo,yyVecArr])
        tab = fits.open(workDir+cfg_par['general']['tableBinName'])
        head = tab[0].header
        #data = tab[0].data
        #tab[1] = vorBinInfo    

        empty_primary = fits.PrimaryHDU(header=head)           

        #t2 = fits.BinTableHDU.from_columns(t,name='vorBinInfo')
        hdul = fits.HDUList([empty_primary])      

        hdul.append(fits.BinTableHDU(t.as_array(), name='vorBinInfo'))


        #hdul.append(t2)  
        hdul.writeto(workDir+cfg_par['general']['outVorTableName'],overwrite=True)

        Lines = np.subtract(data,Stars)

        outputs = cfg_par['starSub']['outputs']
        
        if outputs == 'lines':
            fits.writeto(cfg_par['general']['outLines'],Lines,header,overwrite=True)
            fits.writeto(cfg_par['general']['outNoise'],noiseCube,header,overwrite=True)
            print('''\t+---------+\n\t Line Cube saved\n\t+---------+''')     
            return 
        else: 
            if outputs == 'stars':
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                print('''\t+---------+\n\t Star Cube saved\n\t+---------+''')    
                return 
            elif outputs == 'data':
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('''\t+---------+\n\t Data Cube saved\n\t+---------+''')                 
            elif outputs=='all':
                fits.writeto(cfg_par['general']['outNoise'],noiseCube,header,overwrite=True)
                fits.writeto(cfg_par['general']['outLines'],Lines,header,overwrite=True)
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('''\t+---------+\n\t All Cubes saved\n\t+---------+''')    
                return

    def makeCubesPix(self,cfg_par):

        workDir = cfg_par['general']['workdir']

        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec,dataStar = tP.openPPXFforSubtraction(cfg_par,workDir+cfg_par['general']['tableBinName'],
            workDir+cfg_par['general']['tableAllSpecName'],workDir+cfg_par['general']['tableStarName'])
        

        lineInfo = tP.openLineList(cfg_par)

        dataSpec = np.array(dataSpec['SPEC'][:])
        #dataSpec = np.reshape(dataSpec.T,[dataSpec.shape[1],yAxis.shape[0],xAxis.shape[0]])
        
        dataSub=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        data=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])

        Stars=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        noiseCube=np.empty([len(wave),yAxis.shape[0],xAxis.shape[0]])
        
        header = self.makeHeader(cfg_par, wave, pxSize)

        hdu = fits.PrimaryHDU(data=dataSub,header=header)

        f = fits.open(cfg_par['general']['inputCube'])
        hh = f[0].header
        dd = f[0].data


        waveCube = hh['CRVAL3']+(np.arange(dd.shape[0]))*hh['CD3_3']
        waveCube = waveCube / (1+cfg_par['general']['redshift'])     

        idx = np.where( np.logical_and( waveCube >= cfg_par['starSub']['waveMin'], waveCube <= cfg_par['starSub']['waveMax'] ) )[0]
   
        dd = dd[idx,:,:]

        #spec = np.reshape(dd,[dd.shape[0],dd.shape[1]*dd.shape[2]])
        #idx_good = np.where( np.median(spec, axis=0) > 0.0 )[0]
        #spec     = spec[:,idx_good]

        waveCube = waveCube[idx]
        velscale  = (wave[1]-wave[0])*(cfg_par['general']['C']*1e-3)/np.mean(wave)
        #log_spec, logLam = self.run_logrebinning(spec, velscale,dd.shape[1]*dd.shape[2], waveCube )
        diffusion = 1e-5
        del header['LATPOLE']
        del header['LONPOLE']
        xxVec = []
        yyVec = []
        #create AllSpectra datacube
        for i in range(0,vorBinInfo['ID'].shape[0]):
            #print xAxis
            #print yAxis
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
                #tmpS = tmpS.tolist()
                #print(xx[0],yy[0])
                
                #for j in range(0,len(yy)):
                #   for k in range(0,len(xx)):

                tmpD = np.array(dataSpec[i,:])
                #tmp = tmpD.tolist()
                data[:,yy[0],xx[0]] = tmpD
                if cfg_par['starSub'].get('scaleHow',None) == 'mean':
                    starsScale=np.multiply(tmpS,(np.divide(np.nanmean(tmpD),np.nanmean(tmpS))))                
                elif cfg_par['starSub'].get('scaleHow',None) == 'median':
                    starsScale=np.multiply(tmpS,(np.divide(np.nanmedian(tmpD),np.nanmedian(tmpS))))

                dataSub[:,yy[0],xx[0]] = np.subtract(tmpD,starsScale)
                Stars[:,yy[0],xx[0]] = starsScale

                        #tmpN = np.array(dataSpec[indexBin][1][:])
                        #tmpN = tmpD/tmpN
                        #tmpN = tmpN.tolist()                
                        #noiseCube[:,yy[0],xx[0]] = tmpN
            else:
                pass

        xxVecArr= Column(np.array(xxVec), name='PixX')
        yyVecArr= Column(np.array(yyVec), name='PixY')

        #yyVecArr=np.array(yyVec,dtype={'names':('PixY')})
        #print(vorBinInfo.shape)
        #print(xxVecArr.shape)

        t = Table(vorBinInfo)
        t.add_column(xxVecArr,index=0)
        t.add_column(yyVecArr,index=0) 
        #vorBinInfo = np.column_stack((vorBinInfo,xxVec,yyVec))
        #vorBinInfo = np.vstack([vorBinInfo,yyVecArr])
        tab = fits.open(workDir+cfg_par['general']['tableBinName'])
        head = tab[0].header
        #data = tab[0].data
        #tab[1] = vorBinInfo    

        empty_primary = fits.PrimaryHDU(header=head)           

        #t2 = fits.BinTableHDU.from_columns(t,name='vorBinInfo')
        hdul = fits.HDUList([empty_primary])      

        hdul.append(fits.BinTableHDU(t.as_array(), name='vorBinInfo'))


        #hdul.append(t2)  
        hdul.writeto(workDir+cfg_par['general']['outVorTableName'],overwrite=True)

        #Lines = np.subtract(data,Stars)

        outputs = cfg_par['starSub']['outputs']
        
        if outputs == 'lines':
            fits.writeto(cfg_par['general']['outLines'],dataSub,header,overwrite=True)
            #fits.writeto(cfg_par['general']['outNoise'],noiseCube,header,overwrite=True)
            print('''\t+---------+\n\t Line Cube saved\n\t+---------+''')     
            return 
        else: 
            if outputs == 'stars':
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                print('''\t+---------+\n\t Star Cube saved\n\t+---------+''')    
                return 
            elif outputs == 'data':
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('''\t+---------+\n\t Data Cube saved\n\t+---------+''')                 
            elif outputs=='all':
                fits.writeto(cfg_par['general']['outNoise'],noiseCube,header,overwrite=True)
                fits.writeto(cfg_par['general']['outLines'],Lines,header,overwrite=True)
                fits.writeto(cfg_par['general']['outStars'],Stars,header,overwrite=True)
                fits.writeto(cfg_par['general']['outCube'],data,header,overwrite=True)
                print('''\t+---------+\n\t All Cubes saved\n\t+---------+''')    
                return

    def makeCubesVorLine(self,cfg_par):

        key = 'general'
        workDir = cfg_par['general']['workdir']
        
        #wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec  = tP.openTablesPPXFforSubtraction(cfg_par,workDir+cfg_par['general']['outVorLineName'],
        #    cfg_par['general']['outVorSpectra'])

        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo,dataSpec,dataStar = tP.openPPXFforSubtraction(cfg_par,cfg_par['general']['outVorLineName'],
            cfg_par['general']['outVorSpectra'],workDir+cfg_par['general']['tableStarName']

        print(xAxis.shape,len(wave))
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
        #create AllSpectra datacube
        for i in range(0,vorBinInfo['ID'].shape[0]):
            #print xAxis
            #print yAxis
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

                #tmp = np.array(dataStar[indexBin][1][:])
                #tmp = tmp.tolist()
                #Stars[:,yy[0],xx[0]] = tmp

            else:
                pass

        #xxVecArr= Column(np.array(xxVec), name='PixX')
        #yyVecArr= Column(np.array(yyVec), name='PixY')

        #yyVecArr=np.array(yyVec,dtype={'names':('PixY')})
        #print(vorBinInfo.shape)
        #print(xxVecArr.shape)

        #t = Table(vorBinInfo)
        #t.add_column(xxVecArr,index=0)
        #t.add_column(yyVecArr,index=0) 
        
        #vorBinInfo = np.column_stack((vorBinInfo,xxVec,yyVec))
        #vorBinInfo = np.vstack([vorBinInfo,yyVecArr])
        #tab = fits.open(workDir+cfg_par['general']['tableBinName'])
        #head = tab[0].header
        #data = tab[0].data
        #tab[1] = vorBinInfo    

        #empty_primary = fits.PrimaryHDU(header=head)           

        #t2 = fits.BinTableHDU.from_columns(t,name='vorBinInfo')
        #hdul = fits.HDUList([empty_primary])      

        #hdul.append(fits.BinTableHDU(t.as_array(), name='vorBinInfo'))

        #hdul.append(t2)  
        #hdul.writeto(workDir+cfg_par['general']['outVorLineTableName'],overwrite=True)

        fits.writeto(cfg_par['general']['outVorLines'],data,header,overwrite=True)
            #fits.writeto(cfg_par['general']['outNoise'],noiseCube,header,overwrite=True)
        print('''\t+---------+\n\t Line Cube saved\n\t+---------+''')     
        return 
    


    def makeHeader(self,cfg_par,wave,pxSize):
        
        # these are arbitrary but must correspond, I choose the centre of the galaxy(usually offset with fov)
        crPix3 =cfg_par['starSub']['pixZ']
        crPix1=cfg_par['starSub']['pixX']
        crPix2=cfg_par['starSub']['pixY']

        crVal3 = np.min(wave)
        deltaLambda = (np.max(wave)-np.min(wave))/len(wave)

        #only for Fornax E
        #crVal3 = np.exp(dataSpecStraight[0])
        #crVal3 =4750.2734375
        crpix3 = 1
        crVal1=cfg_par['starSub']['ra']
        crVal2=cfg_par['starSub']['dec']

        ra = cvP.hms2deg(crVal1)
        dec = cvP.dms2deg(crVal2)

        w = wcs.WCS(naxis=3)


        # Set up an "Airy's zenithal" projection
        # Vector properties may be set with Python lists, or Numpy arrays
        w.wcs.crpix = [crPix1, crPix2, crPix3]
        w.wcs.cdelt = np.array([-pxSize,pxSize,deltaLambda])
        w.wcs.crval = [ra, dec, crVal3]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN","AWAV"]

        header = w.to_header()

        #w.wcs.set_pv([(2, 1, 45.0)])
        header['EQUINOX'] = 2000
        header['CRDER3'] = 0.026
        header['CTYPE3'] = "AWAV"
        header['CUNIT3'] = "Angstrom"
        #header['CRPIX1'] = crPix1
        #head['CRPIX2'] = crPix2

        return header
