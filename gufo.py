#!/usr/bin/env python
import os, sys
import yaml


from scavengers import gPlay, tPlay, specPlot

import pkg_resources
try:
    __version__ = pkg_resources.require("gufo")[0].version
except pkg_resources.DistributionNotFound:
    __version__ = "dev"



####################################################################################################


class gufo:

    def __init__(self,file=None):

        #self.rootdir = os.getcwd()+'/'
        self.C = 2.99792458e8

        # get directories
        GFIT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        GFIT_DIR = GFIT_PATH+'/scavengers/'
        sys.path.append(os.path.join(GFIT_PATH, 'gufo'))
        
        file_default = GFIT_DIR + 'gPlay_default.yaml'

        if file != None:
            cfg = open(file)
        else:
            cfg = open(file_default)

        self.cfg_par = yaml.load(cfg)
        if self.cfg_par['general']['verbose'] == True:
            print yaml.dump(self.cfg_par)

        self.set_dirs()

        cfg.close()

        return

    def set_dirs(self):

        runDir = self.cfg_par['general']['workdir']+self.cfg_par['general']['runName']+'/'
        if not os.path.exists(runDir):
            os.mkdir(runDir)
        self.cfg_par['general']['runNameDir'] = runDir

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


        modNameDir = self.cfg_par[key]['runNameDir']+'myMods/'
        cfg_par[key]['modNameDir'] = modNameDir
        if not os.path.exists(modNameDir):
            os.mkdir(modNameDir)
        self.cfg_par['general']['modNameDir'] = momDir

        return

