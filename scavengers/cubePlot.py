#!/usr/bin/env python3.6

import os, sys, shutil
import yaml

import matplotlib.pyplot as plt
import pywcsgrid2
import pywcs

from astropy.io import fits

import mpl_toolkits.axes_grid1.axes_grid as axes_grid
#from mpl_toolkits.axes_grid.colorbar import colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.colors as mcolors

import numpy as np




class Velo(object):
    def __init__(self, header):
        wcs = pywcs.WCS(header)
        self.wcs_vel = wcs.sub([3])

    def to_vel(self, p):
        v = self.wcs_vel.wcs_pix2sky([[p]], 0)
        return v[0][0]

class cubeplot:
    '''Modules to create cubelets of real fitted lines and residuals
    - makeHeader
        make header of line in vrad velocity frame 
    - makeLineCube
        make cubelets for each line marked in lineList.txt
    '''

    def setup_axes( header):

        #gh = pywcsgrid2.GridHelper(wcs=header)
        #gh.locator_params(nbins=3)
        
        wcsIm = WCS(header)        
        wcs2Dim=wcsIm.slice(1,2)
        fig = plt.figure(figsize=(10,10),constrained_layout=False)
        fig.subplots_adjust(hspace=0.,wspace=0.)
        gs = plt.GridSpec(nrows=5, ncols=4,  figure=fig, top=0.95)
       
        return fig,gs,wcs2Dim


    def chanMaps(self):

        fig,gs, wcsIm = setup_axes(header)

        # draw images

        i = 0
        chanStop=20
        ext_xmin=190
        ext_xmax=1160
        ext_ymin=245
        ext_ymax=1120

        c = SkyCoord('00:00:08.0','00:00:10.0',unit=(u.hourangle,u.deg))

        cmap = plt.cm.Blues
        import matplotlib.colors as mcolors
        norm = mcolors.Normalize()
        images = []
        start_channel = 2
                #for plot_count in range(n_plots):
        k=0
        for i in range(0,chanStop):

            if i == 0:
                j = 0
            elif i % 4 == 0:
                j +=1 
                k =0
            
            
            channel_number = chanStop-i
            v0=header['CRVAL3']/1e3-((header['CRPIX3']-chanStop-1)*header['CDELT3'])/1e3
            v = int(v0-((header['CDELT3']*i)/1e3))

            channel = fits_cube[0].data[channel_number,:,:]
            ax = fig.add_subplot(gs[j,k],projection=wcsIm)



            ax.tick_params(axis='both', 
                                     which='major', direction='in')

        #     if k!=0:
        #         ax.coords[1].set_ticklabel_visible(False)
        #     else:
        #         ax.coords[1].set_ticklabel_visible(True)

        #     if j==4:       
        #         ax.coords[0].set_ticklabel_visible(True)        
            
            im = ax.imshow(channel[ext_ymin:ext_ymax,ext_xmin:ext_xmax], origin="lower", 
                           norm=norm, cmap=cmap,aspect='auto',
                           extent=[ext_xmin,ext_xmax,ext_ymin,ext_ymax],
                              vmin=-30.,vmax=55.)

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
            ax.coords[0].set_major_formatter('hh:mm:ss')

            ax.coords[1].set_ticklabel_visible(False)   
            if j==4:    
                ax.coords[0].set_ticklabel_visible(True)  
            else:
                ax.coords[0].set_ticklabel_visible(False)  
            if k==0:    
                ax.coords[1].set_ticklabel_visible(True)  
            else:
                ax.coords[1].set_ticklabel_visible(False)          

            k+=1


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
        plt.savefig(outFigName,bbox_inches='tight',dpi=300,format='png') 
        plt.close()# if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)










