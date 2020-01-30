import sys
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

def update(*a):
    pbar.update()

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
