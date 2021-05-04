#!/usr/bin/env python3.6

import os, sys, shutil
import yaml
import cubePlay
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import patches as mpatches

import mpl_toolkits.axes_grid1.axes_grid as axes_grid
#from mpl_toolkits.axes_grid.colorbar import colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


import numpy as np

from astropy.io import fits

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from reproject import reproject_interp

cP = cubePlay.cubeplay()

class Velo(object):
    def __init__(self, header):
        wcs = WCS(header)
        self.wcs_vel = wcs.sub([3])

    def to_vel(self, p):
        v = self.wcs_vel.wcs_pix2sky([[p]], 0)
        return v[0][0]

class cubeplot:
    '''Modules to plot channel maps of different datacubes
    - setub_axess
        module to set the axes of the datacube given the header
    - chanMaps
        module to plot the channel maps of a datacube
    '''  

    def setup_axes(self, header,nrows):

        #gh = pywcsgrid2.GridHelper(wcs=header)
        #gh.locator_params(nbins=3)
        
        wcsIm = WCS(header)        
        wcs2Dim=wcsIm.slice(1,2)

        if nrows>5 and nrows<=7:
            ncols=5
        elif nrows>7:
            ncols=6
        else:
            ncols=4
        fig = plt.figure(figsize=(ncols*2.5,nrows*2.5),constrained_layout=False)
        fig.subplots_adjust(hspace=0.,wspace=0.)
        gs = plt.GridSpec(nrows=nrows, ncols=ncols,  figure=fig, top=0.95)
       
        return fig,gs,wcs2Dim,ncols


    def chanMaps(self,cfg_par):
        '''Plots the channel maps of a datacube with contours

        Parameters
        ----------

        cfg_par['cubePlay']['chanMaps']: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        Returns:
        --------
        
        outFig: str
            full path to output figure
            
        Requirements:
        -------------

        datacube must exist, CUNIT3 is assumed to be km/s. 
        Axes range must be known and set in the configuration file.
        The step with which plot the channel maps must also be specified.
        '''
        
        inIm=cfg_par['cubePlay']['chanMaps']['inCube']
        header=fits.getheader(inIm)
        fits_cube=fits.getdata(inIm)

        xyzRange = np.array(cfg_par['cubePlay']['chanMaps']['xyzRange'],dtype=int)
        if xyzRange[0,1]==0:
           xyzRange[0,1]= header['NAXIS1']
        if xyzRange[1,1]==0:
           xyzRange[1,1]= header['NAXIS2']


        vel = ((np.linspace(1, fits_cube.shape[0], fits_cube.shape[0]) - header['CRPIX3']) 
            * header['CDELT3'] + header['CRVAL3'])/1e3
        vsys = cfg_par['cubePlay']['chanMaps']['vsys']
        vel -= vsys
        idxMax = np.where(abs(vel-xyzRange[2,0])==abs(vel-xyzRange[2,0]).min())[0][0]
        idxMin = np.where(abs(vel-xyzRange[2,1])==abs(vel-xyzRange[2,1]).min())[0][0]
        step=cfg_par['cubePlay']['chanMaps']['step']/(-header['CDELT3']/1e3)
        nrows= int(round((idxMax-idxMin)/(step*4)))
        fig,gs, wcsIm, ncols = self.setup_axes(header,nrows)


        i = idxMin
        chanStop=idxMax
        ext_xmin=xyzRange[0,0]
        ext_xmax=xyzRange[0,1]
        ext_ymin=xyzRange[1,0]
        ext_ymax=xyzRange[1,1]

        c = SkyCoord(cfg_par['cubePlay']['chanMaps']['raSep'],cfg_par['cubePlay']['chanMaps']['decSep'],unit=(u.hourangle,u.deg))

        cmap = plt.cm.Blues
        norm = mcolors.Normalize()
        images = []
        start_channel = 2

                #for plot_count in range(n_plots):
        k=0
        counter=0
        figs=((chanStop-idxMin)/int(round(step)))
        for i in range(idxMin,chanStop,int(round(step))):

            if counter == 0:
                j = 0
            elif counter % ncols == 0:
                j +=1 
                k =0
            
            channel_number = chanStop - int(round(step))*counter
            #v0=header['CRVAL3']/1e3-((header['CRPIX3']-chanStop)*header['CDELT3'])/1e3
            v = int(vel[channel_number]+vsys)

            channel = fits_cube[channel_number,:,:]
            ax = fig.add_subplot(gs[j,k],projection=wcsIm)

            ax.tick_params(axis='both', 
                                     which='major', direction='in')

        #     if k!=0:
        #         ax.coords[1].set_ticklabel_visible(False)
        #     else:
        #         ax.coords[1].set_ticklabel_visible(True)

        #     if j==4:       
        #         ax.coords[0].set_ticklabel_visible(True)        
            cRange=np.array(cfg_par['cubePlay']['chanMaps']['cRange'],dtype=float)
            im = ax.imshow(channel[ext_ymin:ext_ymax,ext_xmin:ext_xmax], origin="lower", 
                           norm=norm, cmap=str(cfg_par['cubePlay']['chanMaps']['colorMap']),aspect='auto',
                           extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax],
                              vmin=cRange[0],vmax=cRange[1])

            if 'cubeContours' in cfg_par['cubePlay']['chanMaps']:
                imPos = np.nonzero(cfg_par['cubePlay']['chanMaps']['cubeContours'][0,:,0])
                
                ax.contour(channel[ext_ymin:ext_ymax,ext_xmin:ext_xmax],levels=cfg_par['cubePlay']['chanMaps']['cubeContours'][0,imPos,0][0],
                    origin="lower",aspect='auto',linewidths=0.5,
                    colors=cfg_par['cubePlay']['chanMaps']['contourColor'][0],extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])
                
                imNeg = np.nonzero(cfg_par['cubePlay']['chanMaps']['cubeContours'][0,:,1])
                if np.nansum(imNeg) != np.nan:
                    ax.contour(channel[ext_ymin:ext_ymax,ext_xmin:ext_xmax],levels=cfg_par['cubePlay']['chanMaps']['cubeContours'][0,imNeg,1][0],
                        origin="lower",aspect='auto',linestyles = 'dashed',linewidths=0.5,
                        colors=cfg_par['cubePlay']['chanMaps']['contourColor'][1],extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])

            if 'overCubes' in cfg_par['cubePlay']['chanMaps']:
                header2dIm = header.copy()
                cP.delete_3rd_axis(header2dIm)
                wcsIm2d=WCS(header2dIm)
                for h in range(0,len(cfg_par['cubePlay']['chanMaps']['overCubes'])):
                    hduCont = fits.open(cfg_par['cubePlay']['chanMaps']['overCubes'][h])[0]

                    vell = ((np.linspace(1, hduCont.data.shape[0], hduCont.data.shape[0]) - hduCont.header['CRPIX3']) 
                        * hduCont.header['CDELT3'] + hduCont.header['CRVAL3'])/1e3
                    header2d= hduCont.header.copy()
                    cP.delete_3rd_axis(header2d)

                    wcsCont = WCS(header2d)


                    vell -= vsys
                    idx = np.where(abs(vell-vel[channel_number])==abs(vell-vel[channel_number]).min())[0][0]

                    # size= (ext_ymaxCont-ext_yminCont,ext_xmaxCont-ext_xminCont)
                    # centre = (ext_yminCont+size[1]/2.,ext_xminCont+size[0]/2.)
                    # print(size,centre)
                    # hduContCut = Cutout2D(hduCont.data, centre, size)    
                    # print(hduContCut.shape)
                    #hduContCut = Cutout2D(hduCont.data, [vsys,], size)    
                    array, footprint = reproject_interp((hduCont.data[idx,:,:], wcsCont) ,
                                                wcsIm2d, shape_out=[channel.shape[0],channel.shape[1]])

                    #ax1.contour(array.data,levels=imLevels[i,:,0], colors=imColors[i])
                    imPos = np.nonzero(cfg_par['cubePlay']['chanMaps']['otherContours'][h,:,0])
                    ax.contour(array.data,levels=cfg_par['cubePlay']['chanMaps']['otherContours'][h,imPos,0][0],
                        origin="lower",aspect='auto',linewidths=0.5,
                        colors=cfg_par['cubePlay']['chanMaps']['otherColors'][h],extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])
                   
                    imNeg = np.nonzero(cfg_par['cubePlay']['chanMaps']['otherContours'][h,:,1])
                    if np.nansum(imNeg) != np.nan:
                        ax.contour(array.data,levels=cfg_par['cubePlay']['chanMaps']['otherContours'][h,imNeg,1][0],
                            origin="lower",aspect='auto',linestyles = 'dashed',linewidths=0.5,
                            colors=cfg_par['cubePlay']['chanMaps']['otherColors'][h],extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])
           
            patches = [ mpatches.Patch(edgecolor=None, facecolor=None, linewidth=0, linestyle=None,
                                       color='black',  label=str(v)+r' km s$^{-1}$' )]
            # put those patched as legend-handles into the legend
            legend = ax.legend(handles=patches, loc=3, borderaxespad=0.,frameon=False,
                               handlelength=0,handletextpad=0)

            ax.coords[0].set_ticks(spacing=c.ra.degree*u.degree)
            #ax.coords[1].set_ticks(spacing=c.ra.degree*u.degree)

            ax.coords[0].set_axislabel(r'RA (J2000)') 
            ax.coords[1].set_axislabel(r'Dec (J2000)') 
            
            ax.coords[1].set_major_formatter('dd:mm')
            ax.coords[0].set_major_formatter('hh:mm')

            ax.coords[1].set_ticklabel_visible(False)   
            
            # if j==4:    
            #     ax.coords[0].set_ticklabel_visible(True)  
            # else:
            #     ax.coords[0].set_ticklabel_visible(False)  
            if k==0:    
                ax.coords[1].set_ticklabel_visible(True)  
            else:
                ax.coords[1].set_ticklabel_visible(False)          
            if counter >= (figs-ncols):                
                ax.coords[0].set_ticklabel_visible(True)  
            else:
                ax.coords[0].set_ticklabel_visible(False)  
            
            k+=1
            counter+=1

        #     hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        #     array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
        #                                         hduImCut3.wcs, shape_out=hduImCut3.shape)

        #     cs = ax3.contour(array.data,levels=contValues[i], colors=contColors[i])    

        # label with velocities
        # use_path_effect = True
        # try:
        #     from matplotlib.patheffects import withStroke
        # except ImportError:
        #     use_path_effect = False

        # for i, ax in enumerate(g):
        #     channel_number = start_channel + i
        #     v = vel.to_vel(channel_number) / 1.e3
        #     t = ax.add_inner_title(r"$v=%4.1f$ km s$^{-1}$" % (v), loc=2, frameon=False)
        #     if use_path_effect:
        #         t.txt._text.set_path_effects([withStroke(foreground="w",
        #                                                  linewidth=3)])


        # make colorbar
        # cb = plt.colorbar(im, cax=cax)
        # cb.set_label("T [K]")
        # cb.set_ticks([0, 1, 2, 3])

        # # adjust norm
        # norm.vmin = -0.1
        # norm.vmax = 3.5
        # for im in images:
        #     im.changed()

        # if 0:
        #     plt.savefig("co_channel_maps.eps", dpi=70, bbox_inches="tight")
        outFigDir = cfg_par['general']['workdir']+'chanMaps/'
        if not os.path.exists(outFigDir):
            os.mkdir(outFigDir)
        cfg_par['cubePlay']['chanMaps']['outDir'] = outFigDir
        outFigName=outFigDir+cfg_par['cubePlay']['chanMaps']['title']+'.'+cfg_par['cubePlay']['chanMaps']['plotFormat']
        plt.savefig(outFigName,bbox_inches='tight',dpi=300,format=cfg_par['cubePlay']['chanMaps']['plotFormat']) 
        plt.close()# if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        return outFigName

    def chanContoursOverOpt(self,cfg_par):
        '''Plots channel maps contours over an optical image

        Parameters
        ----------

        cfg_par['cubePlay']['chanMaps']: OrderedDictionary
            Dictionary with alla parameters or gufo. 
            Parameters are given from terminal, or in a `.yml' file.

        Returns:
        --------
        
        outFig: str
            full path to output figure
            
        Requirements:
        -------------

        datacube must exist, CUNIT3 is assumed to be km/s. 
        Axes range must be known and set in the configuration file.
        The step with which plot the channel maps must also be specified.
        '''
        
        inIm=cfg_par['cubePlay']['chanMaps']['inOpt']
        imHead = fits.getheader(inIm)

        inCube1=cfg_par['cubePlay']['chanMaps']['overCubes'][0]        
        header=fits.getheader(inCube1)
        fits_cube=fits.getdata(inCube1)



        xyzRange = np.array(cfg_par['cubePlay']['chanMaps']['xyzRange'],dtype=int)
        if xyzRange[0,1]==0:
           xyzRange[0,1]= imHead['NAXIS1']
        if xyzRange[1,1]==0:
           xyzRange[1,1]= imHead['NAXIS2']


        vel = ((np.linspace(1, fits_cube.shape[0], fits_cube.shape[0]) - header['CRPIX3']) 
            * header['CDELT3'] + header['CRVAL3'])/1e3
        vsys = cfg_par['cubePlay']['chanMaps']['vsys']
        vel -= vsys
        idxMax = np.where(abs(vel-xyzRange[2,0])==abs(vel-xyzRange[2,0]).min())[0][0]
        idxMin = np.where(abs(vel-xyzRange[2,1])==abs(vel-xyzRange[2,1]).min())[0][0]
        step=cfg_par['cubePlay']['chanMaps']['step']/(-header['CDELT3']/1e3)
        nrows= int(round((idxMax-idxMin)/(step*4)))
        fig,gs, wcsIm, ncols = self.setup_axes(header,nrows)


        i = idxMin
        chanStop=idxMax
        ext_xmin=xyzRange[0,0]
        ext_xmax=xyzRange[0,1]
        ext_ymin=xyzRange[1,0]
        ext_ymax=xyzRange[1,1]




        c = SkyCoord(cfg_par['cubePlay']['chanMaps']['raSep'],cfg_par['cubePlay']['chanMaps']['decSep'],unit=(u.hourangle,u.deg))

        cmap = plt.cm.Blues
        norm = mcolors.Normalize()
        images = []
        start_channel = 2

                #for plot_count in range(n_plots):
        k=0
        counter=0
        figs=((chanStop-idxMin)/int(round(step)))
        for i in range(idxMin,chanStop,int(round(step))):

            if counter == 0:
                j = 0
            elif counter % ncols == 0:
                j +=1 
                k =0
            
            channel_number = chanStop - int(round(step))*counter
            #v0=header['CRVAL3']/1e3-((header['CRPIX3']-chanStop)*header['CDELT3'])/1e3
            v = int(vel[channel_number]+vsys)

            ax = fig.add_subplot(gs[j,k],projection=wcsIm)

            ax.tick_params(axis='both', 
                                     which='major', direction='in')

        #     if k!=0:
        #         ax.coords[1].set_ticklabel_visible(False)
        #     else:
        #         ax.coords[1].set_ticklabel_visible(True)

        #     if j==4:       
        #         ax.coords[0].set_ticklabel_visible(True)        
            cRange=np.array(cfg_par['cubePlay']['chanMaps']['cRange'],dtype=float)
            im = ax.imshow(imHead[ext_ymin:ext_ymax,ext_xmin:ext_xmax], origin="lower", 
                       norm=norm, cmap=str(cfg_par['cubePlay']['chanMaps']['colorMap']),aspect='auto',
                       extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax],
                          vmin=cRange[0],vmax=cRange[1])

            for h in range(0,len(cfg_par['cubePlay']['chanMaps']['overCubes'])):
                    hduCont = fits.open(cfg_par['cubePlay']['chanMaps']['overCubes'][h])[0]
                    
                    header2dIm = imHead.copy()
                    cP.delete_3rd_axis(header2dIm)
                    wcsIm2d=WCS(header2dIm)
                    
                    vell = ((np.linspace(1, hduCont.data.shape[0], hduCont.data.shape[0]) - hduCont.header['CRPIX3']) 
                        * hduCont.header['CDELT3'] + hduCont.header['CRVAL3'])/1e3
                    header2d= hduCont.header.copy()
                    cP.delete_3rd_axis(header2d)

                    wcsCont = WCS(header2d)


                    vell -= vsys
                    idx = np.where(abs(vell-vel[channel_number])==abs(vell-vel[channel_number]).min())[0][0]

                    # size= (ext_ymaxCont-ext_yminCont,ext_xmaxCont-ext_xminCont)
                    # centre = (ext_yminCont+size[1]/2.,ext_xminCont+size[0]/2.)
                    # print(size,centre)
                    # hduContCut = Cutout2D(hduCont.data, centre, size)    
                    # print(hduContCut.shape)
                    #hduContCut = Cutout2D(hduCont.data, [vsys,], size)    
                    array, footprint = reproject_interp((hduCont.data[idx,:,:], wcsCont) ,
                                                wcsIm2d, shape_out=[channel.shape[0],channel.shape[1]])


                    #ax1.contour(array.data,levels=imLevels[i,:,0], colors=imColors[i])
                    imPos = np.nonzero(cfg_par['cubePlay']['chanMaps']['otherContours'][h,:,0])
                    ax.contour(array.data,levels=cfg_par['cubePlay']['chanMaps']['otherContours'][h,imPos,0][0],
                        origin="lower",aspect='auto',linewidths=0.5,
                        colors=cfg_par['cubePlay']['chanMaps']['otherColors'][h],extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])
                   
                    imNeg = np.nonzero(cfg_par['cubePlay']['chanMaps']['otherContours'][h,:,1])
                    if np.nansum(imNeg) != np.nan:
                        ax.contour(array.data,levels=cfg_par['cubePlay']['chanMaps']['otherContours'][h,imNeg,1][0],
                            origin="lower",aspect='auto',linestyles = 'dashed',linewidths=0.5,
                            colors=cfg_par['cubePlay']['chanMaps']['otherColors'][h],extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax])
           
            patches = [ mpatches.Patch(edgecolor=None, facecolor=None, linewidth=0, linestyle=None,
                                       color='black',  label=str(v)+r' km s$^{-1}$' )]
            # put those patched as legend-handles into the legend
            legend = ax.legend(handles=patches, loc=3, borderaxespad=0.,frameon=False,
                               handlelength=0,handletextpad=0)

            ax.coords[0].set_ticks(spacing=c.ra.degree*u.degree)
            #ax.coords[1].set_ticks(spacing=c.ra.degree*u.degree)

            ax.coords[0].set_axislabel(r'RA (J2000)') 
            ax.coords[1].set_axislabel(r'Dec (J2000)') 
            
            ax.coords[1].set_major_formatter('dd:mm')
            ax.coords[0].set_major_formatter('hh:mm')

            ax.coords[1].set_ticklabel_visible(False)   
            
            # if j==4:    
            #     ax.coords[0].set_ticklabel_visible(True)  
            # else:
            #     ax.coords[0].set_ticklabel_visible(False)  
            if k==0:    
                ax.coords[1].set_ticklabel_visible(True)  
            else:
                ax.coords[1].set_ticklabel_visible(False)          
            if counter >= (figs-ncols):                
                ax.coords[0].set_ticklabel_visible(True)  
            else:
                ax.coords[0].set_ticklabel_visible(False)  
            
            k+=1
            counter+=1

        #     hduContCut = Cutout2D(hduCont.data, centre, size, wcs=wcsCont)    
        #     array, footprint = reproject_interp((hduContCut.data, hduContCut.wcs) ,
        #                                         hduImCut3.wcs, shape_out=hduImCut3.shape)

        #     cs = ax3.contour(array.data,levels=contValues[i], colors=contColors[i])    

        # label with velocities
        # use_path_effect = True
        # try:
        #     from matplotlib.patheffects import withStroke
        # except ImportError:
        #     use_path_effect = False

        # for i, ax in enumerate(g):
        #     channel_number = start_channel + i
        #     v = vel.to_vel(channel_number) / 1.e3
        #     t = ax.add_inner_title(r"$v=%4.1f$ km s$^{-1}$" % (v), loc=2, frameon=False)
        #     if use_path_effect:
        #         t.txt._text.set_path_effects([withStroke(foreground="w",
        #                                                  linewidth=3)])


        # make colorbar
        # cb = plt.colorbar(im, cax=cax)
        # cb.set_label("T [K]")
        # cb.set_ticks([0, 1, 2, 3])

        # # adjust norm
        # norm.vmin = -0.1
        # norm.vmax = 3.5
        # for im in images:
        #     im.changed()

        # if 0:
        #     plt.savefig("co_channel_maps.eps", dpi=70, bbox_inches="tight")
        outFigDir = cfg_par['general']['workdir']+'chanMaps/'
        if not os.path.exists(outFigDir):
            os.mkdir(outFigDir)
        cfg_par['cubePlay']['chanMaps']['outDir'] = outFigDir
        outFigName=outFigDir+cfg_par['cubePlay']['chanMaps']['title']+'.'+cfg_par['cubePlay']['chanMaps']['plotFormat']
        plt.savefig(outFigName,bbox_inches='tight',dpi=300,format=cfg_par['cubePlay']['chanMaps']['plotFormat']) 
        plt.close()# if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        return outFigName







