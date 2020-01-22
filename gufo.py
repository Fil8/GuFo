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

        #self.rootdir = os.getcwd()+'/'
        self.C = 2.99792458e8

        # get directories
        GFIT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        sys.path.append(os.path.join(GFIT_PATH, 'gufo'))
        
        GFIT_DIR = GFIT_PATH+'/gufo/'

        file_default = GFIT_DIR + 'gPlay_default.yaml'

        if file != None:
            cfg = open(file)
        else:
            cfg = open(file_default)

        self.cfg_par = yaml.load(cfg)
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
   
        if self.cfg_par['starSub'].get('enable',False) == True:
            self.cfg_par['general']['outVorTableName'] = self.cfg_par['general']['workdir']+'GuFo_LineTable.fits'
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
            self.cfg_par['general']['outVorTableName'] = self.cfg_par['general']['workdir']+'GuFo_LineTable.fits'
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

        
        self.cfg_par['general']['outVorSpectra'] = self.cfg_par['general']['workdir']+'GuFo_LineVorSpectra.fits'
        self.cfg_par['general']['outVorLineTableName'] = self.cfg_par['general']['workdir']+'GuFo_LineVorTable.fits'
        self.cfg_par['general']['outVorLines'] = cubeDir+'LineVor.fits'
        self.cfg_par['general']['outVorNoise'] = cubeDir+'NoiseVor.fits'

        outTableName = self.cfg_par['general']['runNameDir']+'gPlayOut.fits'

        self.cfg_par['general']['outTableName'] = outTableName

        outPlotDir = self.cfg_par['general']['runNameDir']+'spectra/'
        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)
        
        self.cfg_par['general']['outPlotDir'] = outPlotDir

        momDir =self. cfg_par['general']['runNameDir']+'moments/'
        if not os.path.exists(momDir):
            os.mkdir(momDir)
        self.cfg_par['general']['momDir'] = momDir

        bptDir =self. cfg_par['general']['runNameDir']+'bpt/'
        if not os.path.exists(bptDir):
            os.mkdir(bptDir)
        self.cfg_par['general']['bptDir'] = bptDir

        modNameDir = self.cfg_par['general']['runNameDir']+'models/'
        if not os.path.exists(modNameDir):
            os.mkdir(modNameDir)
        self.cfg_par['general']['modNameDir'] = modNameDir

        return

