#!/usr/bin/env python3.6
import os, sys
import yaml


from scavengers import gPlay, tPlay, specPlot

import pkg_resources
try:
    __version__ = pkg_resources.require("gufo")[0].version
except pkg_resources.DistributionNotFound:
    __version__ = "dev"



####################################################################################################


class gufo(object):

    def __init__(self,file=None):
        '''
        Set paths and constants
        Load config file
        Call self.set_dirs()
        '''

        #self.rootdir = os.getcwd()+'/'
        self.C = 2.99792458e8

        # get directories
        GFIT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        sys.path.append(os.path.join(GFIT_PATH, 'gufo'))
        
        GFIT_DIR = GFIT_PATH+'/gufo/'

        file_default = GFIT_DIR + 'gufo_default.yaml'

        if file != None:
            cfg = open(file)
        else:
            cfg = open(file_default)

        self.cfg_par = yaml.load(cfg, yaml.SafeLoader)
        if self.cfg_par['general']['verbose'] == True:
            print(yaml.dump(self.cfg_par))
        self.cfg_par['general']['gfitPath'] = GFIT_DIR
        self.cfg_par['general']['C'] = self.C

        self.set_dirs()

        cfg.close()

        return 

    def set_dirs(self):

        runDir = self.cfg_par['general']['workdir']+self.cfg_par['general']['runName']+'/'
        if not os.path.exists(runDir):
            os.mkdir(runDir)
        self.cfg_par['general']['runNameDir'] = runDir

        cubeDir = self.cfg_par['general']['runNameDir']+'cubes/'
        if not os.path.exists(cubeDir):
            os.mkdir(cubeDir)
        self.cfg_par['general']['cubeDir'] = cubeDir
 
        cubeletsDir = self.cfg_par['general']['cubeDir']+'cubelets/'
        if not os.path.exists(cubeletsDir):
            os.mkdir(cubeletsDir)
        self.cfg_par['general']['cubeletsDir'] = cubeletsDir

        pvDir = self.cfg_par['general']['cubeletsDir']+'pv/'
        if not os.path.exists(pvDir):
            os.mkdir(pvDir )
        self.cfg_par['general']['pvDir'] = pvDir

        tableDir = self.cfg_par['general']['runNameDir']+'tables/'
        if not os.path.exists(tableDir):
            os.mkdir(tableDir)
        self.cfg_par['general']['tableDir'] = tableDir

        if self.cfg_par['starSub'].get('enable',False) == True:
            self.cfg_par['general']['outVorTableName'] = self.cfg_par['general']['tableDir']+'GuFo_LineTable.fits'
            if self.cfg_par['starSub'].get('scaleFlux',False) == True:
            
                if self.cfg_par['starSub'].get('scaleHow',None) == 'mean':
                    nameEnd = 'PixMean.fits'
                elif self.cfg_par['starSub'].get('scaleHow',None) == 'median':
                    nameEnd = 'PixMed.fits'
                self.cfg_par['general']['outCube'] = cubeDir+'CubePix'+nameEnd
                self.cfg_par['general']['outStars'] = cubeDir+'StarCube'+nameEnd
                self.cfg_par['general']['outLines'] = cubeDir+'LineCube'+nameEnd
                self.cfg_par['general']['outNoise'] = cubeDir+'noiseCube'+nameEnd
            
            elif self.cfg_par['starSub'].get('scaleFlux',False) == False:
                self.cfg_par['general']['outCube'] = cubeDir+'CubeCubeVor.fits'
                self.cfg_par['general']['outStars'] = cubeDir+'StarCubeVor.fits'
                self.cfg_par['general']['outLines'] = cubeDir+'LineCubeVor.fits'
                self.cfg_par['general']['outNoise'] = cubeDir+'noiseCubeVor.fits'
        else: 
            self.cfg_par['general']['outVorTableName'] = self.cfg_par['general']['tableDir']+'GuFo_LineTable.fits'
            if self.cfg_par['starSub'].get('scaleFlux',False) == True:
                if self.cfg_par['starSub'].get('scaleHow',None) == 'mean':
                    nameEnd = 'PixMean.fits'
                elif self.cfg_par['starSub'].get('scaleHow',None) == 'median':
                    nameEnd = 'PixMed.fits'
                self.cfg_par['general']['outLines'] = cubeDir+'LineCube'+nameEnd
                self.cfg_par['general']['outNoise'] = cubeDir+'noiseCube'+nameEnd
            elif self.cfg_par['starSub'].get('scaleFlux',False) == False:
                self.cfg_par['general']['outLines'] = cubeDir+'LineCubeVor.fits'
                self.cfg_par['general']['outNoise'] = cubeDir+'noiseCubeVor.fits'

        
        self.cfg_par['general']['outVorSpectra'] = self.cfg_par['general']['tableDir']+'GuFo_LineVorSpectra.fits'
        self.cfg_par['general']['outVorLineTableName'] = self.cfg_par['general']['tableDir']+'GuFo_LineVorTable.fits'
        self.cfg_par['general']['outVorLines'] = cubeDir+'LineVor.fits'
        self.cfg_par['general']['outVorNoise'] = cubeDir+'NoiseVor.fits'

        if not self.cfg_par['general'].get('dataCubeName',None):
           self.cfg_par['general']['dataCubeName'] =  self.cfg_par['general']['outVorLines'] 
        if not self.cfg_par['general'].get('noiseCubeName',None):
           self.cfg_par['general']['noiseCubeName'] =  self.cfg_par['general']['outVorNoise'] 
        if self.cfg_par['gFit']['method'] == 'pixel':
            self.cfg_par['general']['dataCubeName'] = cubeDir+self.cfg_par['general']['outCube']
            self.cfg_par['general']['outPixSpectra'] = self.cfg_par['general']['tableDir']+'GuFo_LinePixSpectra.fits'
            self.cfg_par['general']['outVorSpectra'] = self.cfg_par['general']['outPixSpectra']
            self.cfg_par['general']['outVorLineTableName'] =self.cfg_par['general']['outVorTableName']
            self.cfg_par['general']['outVorLines'] = self.cfg_par['general']['outLines']

        outTableName = self.cfg_par['general']['runNameDir']+self.cfg_par['general']['outTable']
        outTableName2 = self.cfg_par['general']['runNameDir']+self.cfg_par['general']['outTable2']

        self.cfg_par['general']['outTableName'] = outTableName
        self.cfg_par['general']['outTableName2'] = outTableName2

        outPlotDir = self.cfg_par['general']['runNameDir']+'spectra/'
        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)
        
        self.cfg_par['general']['outPlotDir'] = outPlotDir

        ccaDir =self.cfg_par['general']['runNameDir']+'CCA/'
        self.cfg_par['general']['CCAanalysis'] = ccaDir
        ccaPlotDir =self.cfg_par['general']['runNameDir']+'CCA/plots'

        self.cfg_par['general']['CCAanalysisPlotDir'] = ccaPlotDir

        momDir =self.cfg_par['general']['runNameDir']+'moments/'
        if not os.path.exists(momDir):
            os.mkdir(momDir)
        self.cfg_par['general']['momDir'] = momDir

        momModDir =self.cfg_par['general']['momDir']+self.cfg_par['gFit']['modName']+'/'
        if not os.path.exists(momModDir):
            os.mkdir(momModDir)
        self.cfg_par['general']['momModDir'] = momModDir

        plotMomModDir =self.cfg_par['general']['momModDir']+'plots/'
        if not os.path.exists(plotMomModDir):
            os.mkdir(plotMomModDir)
        self.cfg_par['general']['plotMomModDir'] = plotMomModDir

        bptDir =self.cfg_par['general']['runNameDir']+'lineRatios/'
        if not os.path.exists(bptDir):
            os.mkdir(bptDir)
        self.cfg_par['general']['bptDir'] = bptDir

        bptPlotDir =self.cfg_par['general']['bptDir']+'plots/'
        if not os.path.exists(bptPlotDir):
            os.mkdir(bptPlotDir)
        self.cfg_par['general']['bptPlotDir'] = bptPlotDir

        bptModDir =self.cfg_par['general']['runNameDir']+'bpt/'+self.cfg_par['gFit']['modName']+'/'
        if not os.path.exists(bptModDir):
            os.mkdir(bptModDir)
        self.cfg_par['general']['bptDir'] = bptModDir

        modNameDir = self.cfg_par['general']['runNameDir']+'models/'
        if not os.path.exists(modNameDir):
            os.mkdir(modNameDir)
        self.cfg_par['general']['modNameDir'] = modNameDir

        resDir =self.cfg_par['general']['runNameDir']+'residuals/'
        if not os.path.exists(resDir):
            os.mkdir(resDir)
        self.cfg_par['general']['resDir'] = resDir

        resModDir =self.cfg_par['general']['resDir']+self.cfg_par['gFit']['modName']+'/'
        if not os.path.exists(resModDir):
            os.mkdir(resModDir)
        self.cfg_par['general']['resModDir'] = resModDir

        noiseDir =self.cfg_par['general']['runNameDir']+'residuals/noise/'
        if not os.path.exists(noiseDir):
            os.mkdir(noiseDir)
        self.cfg_par['general']['noiseDir'] = noiseDir

        return

