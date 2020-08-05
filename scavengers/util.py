import sys, os
import yaml

from tqdm import tqdm


def loadCfg(file=None):
    #self.rootdir = os.getcwd()+'/'
    C = 2.99792458e8

    # get directories
    GFIT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(os.path.join(GFIT_PATH, 'gufo'))

    GFIT_DIR = GFIT_PATH+'/gufo/'

    file_default = GFIT_DIR + 'gufo_default.yaml'

    if file != None:
        cfg = open(file)
    else:
        cfg = open(file_default)

    cfg_par = yaml.load(cfg, yaml.SafeLoader)
    if cfg_par['general']['verbose'] == True:
        print(yaml.dump(cfg_par))
    cfg_par['general']['gfitPath'] = GFIT_DIR
    cfg_par['general']['C'] = C

    #set_dirs()

    cfg.close()

    return cfg_par


def set_dirs(cfg_par):

    runDir = cfg_par['general']['workdir']+cfg_par['general']['runName']+'/'
    if not os.path.exists(runDir):
        os.mkdir(runDir)
    cfg_par['general']['runNameDir'] = runDir

    cubeDir = cfg_par['general']['runNameDir']+'cubes/'
    if not os.path.exists(cubeDir):
        os.mkdir(cubeDir)
    cfg_par['general']['cubeDir'] = cubeDir

    tableDir = cfg_par['general']['runNameDir']+'tables/'
    if not os.path.exists(tableDir):
        os.mkdir(tableDir)
    cfg_par['general']['tableDir'] = tableDir


    if cfg_par['starSub'].get('enable',False) == True:
        cfg_par['general']['outVorTableName'] = cfg_par['general']['tableDir']+'GuFo_LineTable.fits'
        if cfg_par['starSub'].get('scaleFlux',False) == True:
        
            if cfg_par['starSub'].get('scaleHow',None) == 'mean':
                nameEnd = 'PixMean.fits'
            elif cfg_par['starSub'].get('scaleHow',None) == 'median':
                nameEnd = 'PixMed.fits'
            cfg_par['general']['outCube'] = cubeDir+'CubePix'+nameEnd
            cfg_par['general']['outStars'] = cubeDir+'StarCube'+nameEnd
            cfg_par['general']['outLines'] = cubeDir+'LineCube'+nameEnd
            cfg_par['general']['outNoise'] = cubeDir+'noiseCube'+nameEnd
        
        elif cfg_par['starSub'].get('scaleFlux',False) == False:
            cfg_par['general']['outCube'] = cubeDir+'CubeCubeVor.fits'
            cfg_par['general']['outStars'] = cubeDir+'StarCubeVor.fits'
            cfg_par['general']['outLines'] = cubeDir+'LineCubeVor.fits'
            cfg_par['general']['outNoise'] = cubeDir+'noiseCubeVor.fits'
    else: 
        cfg_par['general']['outVorTableName'] = cfg_par['general']['tableDir']+'GuFo_LineTable.fits'
        if cfg_par['starSub'].get('scaleFlux',False) == True:
            if cfg_par['starSub'].get('scaleHow',None) == 'mean':
                nameEnd = 'PixMean.fits'
            elif cfg_par['starSub'].get('scaleHow',None) == 'median':
                nameEnd = 'PixMed.fits'
            cfg_par['general']['outLines'] = cubeDir+'LineCube'+nameEnd
            cfg_par['general']['outNoise'] = cubeDir+'noiseCube'+nameEnd
        elif cfg_par['starSub'].get('scaleFlux',False) == False:
            cfg_par['general']['outLines'] = cubeDir+'LineCubeVor.fits'
            cfg_par['general']['outNoise'] = cubeDir+'noiseCubeVor.fits'

    
    cfg_par['general']['outVorSpectra'] = cfg_par['general']['tableDir']+'GuFo_LineVorSpectra.fits'
    cfg_par['general']['outVorLineTableName'] = cfg_par['general']['tableDir']+'GuFo_LineVorTable.fits'
    cfg_par['general']['outVorLines'] = cubeDir+'LineVor.fits'
    cfg_par['general']['outVorNoise'] = cubeDir+'NoiseVor.fits'

    if not cfg_par['general'].get('dataCubeName',None):
       cfg_par['general']['dataCubeName'] =  cfg_par['general']['outVorLines'] 
    if not cfg_par['general'].get('noiseCubeName',None):
       cfg_par['general']['noiseCubeName'] =  cfg_par['general']['outVorNoise'] 

    outTableName = cfg_par['general']['runNameDir']+'gPlayOut.fits'

    cfg_par['general']['outTableName'] = outTableName

    outPlotDir = cfg_par['general']['runNameDir']+'spectra/'
    if not os.path.exists(outPlotDir):
        os.mkdir(outPlotDir)
    
    cfg_par['general']['outPlotDir'] = outPlotDir

    momDir =cfg_par['general']['runNameDir']+'moments/'
    if not os.path.exists(momDir):
        os.mkdir(momDir)
    cfg_par['general']['momDir'] = momDir

    momModDir =cfg_par['general']['momDir']+cfg_par['gFit']['modName']+'/'
    if not os.path.exists(momModDir):
        os.mkdir(momModDir)
    cfg_par['general']['momModDir'] = momModDir

    bptDir =cfg_par['general']['runNameDir']+'bpt/'
    if not os.path.exists(bptDir):
        os.mkdir(bptDir)
    cfg_par['general']['bptDir'] = bptDir

    modNameDir = cfg_par['general']['runNameDir']+'models/'
    if not os.path.exists(modNameDir):
        os.mkdir(modNameDir)
    cfg_par['general']['modNameDir'] = modNameDir

    return cfg_par



# ==============================================================================
#                          P R E T T Y   O U T P U T
# ==============================================================================
""" A collection of functions to generate the pretty output in stdout. """
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


def printProgress(iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 80, color = 'g'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        color       - Optional  : color identifier (Str)
    """
    if   color == 'y': color = '\033[43m'
    elif color == 'k': color = '\033[40m'
    elif color == 'r': color = '\033[41m'
    elif color == 'g': color = '\033[42m'
    elif color == 'b': color = '\033[44m'
    elif color == 'm': color = '\033[45m'
    elif color == 'c': color = '\033[46m'

    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = color + ' '*filledLength + '\033[49m' + ' '*(barLength - filledLength -1)
    sys.stdout.write('\r%s |%s| %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()

def prettyOutput_Running(outputlabel):
    sys.stdout.write("\r [ "+'\033[0;37m'+"RUNNING "+'\033[0;39m'+"] "+outputlabel)
    sys.stdout.flush(); print("")

def prettyOutput_Done(outputlabel, progressbar=False):
    if progressbar == True:
        sys.stdout.write("\033[K")
    sys.stdout.write("\033[F"); sys.stdout.write("\033[K")
    sys.stdout.write("\r\r [ "+'\033[0;32m'+"DONE    "+'\033[0;39m'+"] "+outputlabel)
    sys.stdout.flush(); print("")

def prettyOutput_Warning(outputlabel, progressbar=False):
    if progressbar == True:
        sys.stdout.write("\033[K")
    sys.stdout.write("\033[F"); sys.stdout.write("\033[K")
    sys.stdout.write("\r\r [ "+'\033[0;33m'+"WARNING "+'\033[0;39m'+"] "+outputlabel)
    sys.stdout.flush(); print("")

def prettyOutput_Failed(outputlabel, progressbar=False):
    if progressbar == True:
        sys.stdout.write("\033[K")
    sys.stdout.write("\033[F"); sys.stdout.write("\033[K")
    sys.stdout.write("\r\r [ "+'\033[0;31m'+"FAILED  "+'\033[0;39m'+"] "+outputlabel)
    sys.stdout.flush(); print("")

def prettyOutput_DonePrefix():
    return(" [ "+'\033[0;32m'+"DONE    "+'\033[0;39m'+"] ")

def prettyOutput_WarningPrefix():
    return(" [ "+'\033[0;33m'+"WARNING "+'\033[0;39m'+"] ")

def prettyOutput_FailedPrefix():
    return(" [ "+'\033[0;31m'+"FAILED  "+'\033[0;39m'+"] ")
