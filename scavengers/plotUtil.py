#!/usr/bin/env python3.6

'''

These define the `rcparameters` of GaNGiaLF plots.
Plot specifics are thought for publication. 
Two default sets of parameters are possible: single-column and fullwidth Figures in A&A.
Default fontsize for a fullwidth Figure is 16pt. Default A&A fontsize is 11pt.
Single-column figures will be automatically halfed in size by latTeX, hence fontsizes double: 32
'''

import sys, os

def loadRcParams(option):
    '''Loads the standard rc parameters for uniform plots.

    
        Parameters
        ----------
        option: specify figure size
                - single-column (sc), square
                - fullwidth (fw), square
                - fullwidth (fwR), rectangular (16:9)

    Returns
    -------
    dict()

    '''
    if option='fw'
        figSize= '7.24409,7.24409'
        font=16
    elif option='sc'
        figSize= '3.54331,3.54331'
        font=16
    elif option='fwR'
        figSize= '7.24409,4.074800625'
        font=16

    params = {'figure.figsize'      : figSize,
        'figure.autolayout' : True,
        'font.family'         :'serif',
        'pdf.fonttype'        : 3,
        'font.serif'          :'times',
        'font.style'          : 'normal',
        'font.weight'         : 'book',
        'font.size'           : font,
        'axes.linewidth'      : 1.5,
        'lines.linewidth'     : 1.,
        'xtick.labelsize'     : font,
        'ytick.labelsize'     : font,
        'legend.fontsize'     : font, 
        'xtick.direction'     :'in',
        'ytick.direction'     :'in',
        'xtick.major.size'    : 3,
        'xtick.major.width'   : 1.5,
        'xtick.minor.size'    : 2.5,
        'xtick.minor.width'   : 1.,
        'ytick.major.size'    : 3,
        'ytick.major.width'   : 1.5,
        'ytick.minor.size'    : 2.5,
        'ytick.minor.width'   : 1., 
        'text.usetex'         : True,
        'text.latex.preamble' : r'\usepackage{amsmath}',
        'text.latex.preamble' : r'\usepackage{lmodern}',    # latin modern, recommended to replace computer modern sans serif
        'text.latex.preamble' : r'\usepackage{helvet}',    # set the normal font here
        'text.latex.preamble' : r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        'text.latex.preamble' : r'\sansmath' 
        #'text.latex.unicode'  : True
         }
  
    return params