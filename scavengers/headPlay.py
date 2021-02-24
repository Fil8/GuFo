#!/usr/bin/env python

'''

Set of tools for head-play with fits files.

'''


import sys, os, string
import numpy as np
from astropy.io import fits

class headplay():
    '''This class plays modifies the header of .fits files in useful various ways.
    
    '''

    def cleanHead(self,inFile,writeFile=False):
        '''Cleans the header of 2D images so that wcs libraries do not have issues.
            Squeezes the shapes of data to 2-dimensions accordingly.

            Parameters
            ----------
            inFile: str
                full path of input .fits file

            writeFile: bool, optional
                if set to True inFile is overwritten

            Returns
            -------
            dd: np.array
                2D-array with data of fits file
            hh: header
                modified fits header

        '''
        base = fits.open(inFile)
        heads = base[0].header
        
        datas = base[0].data
        datas = np.squeeze(datas)

        if 'NAXIS3' in heads:
            del heads['NAXIS3']
        if 'CRVAL3' in heads:
            heads['FREQ'] = heads['CRVAL3']
            del heads['CRVAL3']
        if 'CDELT3' in heads:
            del heads['CDELT3']
        if 'CRPIX3' in heads: 
            del heads['CRPIX3']
        if 'CTYPE3' in heads:        
            del heads['CTYPE3']  
        if 'CROTA3' in heads:
            del heads['CROTA3']
        if 'CUNIT3' in heads:
            del heads['CUNIT3']        
            
        if 'NAXIS4' in heads:
            del heads['NAXIS4']     
        if 'CRVAL4' in heads:        
            del heads['CRVAL4']
        if 'CDELT4' in heads:    
            del heads['CDELT4']
        if 'CRPIX4' in heads:    
            del heads['CRPIX4']
        if 'CTYPE4' in heads:    
            del heads['CTYPE4'] 
        if 'CROTA4' in heads:
            del heads['CROTA4']  
        if 'CUNIT4' in heads:
            del heads['CUNIT4']

        if 'DATAMIN' in heads:
            heads.set('DATAMIN', heads['DATAMIN'])
        if 'DATAMAC' in heads:    
            heads.set('DATAMAX', heads['DATAMAX'])
        
        if 'CELLSCAL' in heads:        
            heads.set('CELLSCAL', heads['CELLSCAL'])

        if 'RESTFREQ' in heads:
            heads.set('RESTFREQ', heads['RESTFREQ'])
        
        if 'SPECSYS3' in heads:
            heads.set('SPECSYS3', heads['SPECSYS3'])

        
        if 'WCSAXES' in heads:
            heads['WCSAXES'] = 2
        
        if 'PC1_1' in heads:   
            del heads['PC1_1']            
        if 'PC2_1' in heads:   
            del heads['PC2_1']           
        if 'PC3_1' in heads:   
            del heads['PC3_1']
        if 'PC4_1' in heads:   
            del heads['PC4_1']    
        if 'PC1_2' in heads:
            del heads['PC1_2']
        if 'PC2_2' in heads:
            del heads['PC2_2']
        if 'PC3_2' in heads:
            del heads['PC3_2']
        if 'PC4_2' in heads:
            del heads['PC4_2']
        if 'PC1_3' in heads:
            del heads['PC1_3']
        if 'PC2_3' in heads:
            del heads['PC2_3']
        if 'PC3_3' in heads:    
            del heads['PC3_3']
        if 'PC4_3' in heads:    
            del heads['PC4_3']
        if 'PC1_4' in heads:
            del heads['PC1_4']            
        if 'PC2_4' in heads:
            del heads['PC2_4']
        if 'PC3_4' in heads:
            del heads['PC3_4']
        if 'PC4_4' in heads:
            del heads['PC4_4']
        if 'PV2_2' in heads:
            del heads['PV2_1']
        if 'PV2_2' in heads:
            del heads['PV2_2']
        if 'CD1_1' in heads:
            heads['CDELT1'] = heads['CD1_1']
            del heads['CD1_1']

        if 'CD2_2' in heads:
            heads['CDELT2'] = heads['CD2_2']
            del heads['CD2_2']

        if 'CD1_2' in heads:
            del heads['CD1_2']

        if 'CD2_1' in heads:
            del heads['CD2_1']

        heads['NAXIS'] = 2

        if writeFile == True:
            fits.writeto(inFile,datas,heads,overwrite=True)

        return heads, datas

    def putHead(self,inFile,key,value,outFile=False):
        '''Puts a key header and value into a fits file.

            Parameters
            ----------
            inFile: str
                full path of input .fits file

            key: str
                name of the key to set

            value: str
                value of the key

            outFile: str, optional
                ful path of output .fits file, if given file is not overwritten

            Returns
            -------
            output: str
                equal to inFile, if output is not given

        '''
        files=fits.open(inFile)

        datas=files[0].data
        heads=files[0].header

        heads[key] = value
        # if 'DATAMIN' in heads:
        #     del heads['DATAMIN']
        # if 'DATAMAX' in heads:
        #     del heads['DATAMAX']

        if outFile==False:
            outFile=inFile

        fits.writeto(outFile,datas,heads,overwrite=True)

        return outFile

    def delHead(self,inFile,key,outFile=False):
        '''Deletes a key header.

            Parameters
            ----------
            inFile: str
                full path of input .fits file

            key: str
                name of the key to set

            outFile: str, optional
                ful path of output .fits file, if given file is not overwritten

            Returns
            -------
            output: str
                equal to inFile, if output is not given

        '''
        files=fits.open(inFile)

        datas=files[0].data
        heads=files[0].header

        del heads[key]        

        if outFile==False:
            outFile=inFile

        fits.writeto(outFile,datas,heads,overwrite=True)

        return outFile