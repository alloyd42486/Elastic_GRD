#
# Collection of routines based on Al-Attar et al. 2023 to make geographic plots
#

import numpy as np
import pyshtools as pysh
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


########################### Plot Global Field  #######################
#
# function to plot geographic data on GL grid.
#
def plot(fun,**kwargs):

    # deal with optional arguments
    cstring = kwargs.get('cmap',None)
    ofile = kwargs.get('ofile',None)
    contour = kwargs.get('contour',False)
    ncont = kwargs.get('ncont',6)
    label = kwargs.get('label','')
    marker = kwargs.get('marker',None)
    clim = kwargs.get('clim',None)
    clim_sym = kwargs.get('clim_sym',True)
    clim_pos = kwargs.get('clim_pos',False)
    clim_neg = kwargs.get('clim_neg',False)
    clim_scale = kwargs.get('clim_scale',None)
    xlim =  kwargs.get('xlim',None)
    ylim =  kwargs.get('ylim',None)

    if(cstring == None):
        if(clim_pos):
            cstring = "Blues"
        elif(clim_neg):
            cstring = "Reds_r"
        else:
            cstring = "RdBu"
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

    #Cartopy and options within it and dependecies sometimes are limited to -180 to 180 long
    #pyshtools is limited to 0 to 360 lon
    #here we convert from 0 to 360 lon to -180 to 180 lon
    par = np.zeros([fun.nlat,fun.nlon-1]) #dropping wrapp-around at 0 & 360
    lon = np.zeros(fun.nlon-1)

    imid = np.floor_divide(fun.nlon,2)
    lon[0:imid] = fun.lons()[imid:fun.nlon] - 360.0
    lon[imid:fun.nlon-1] = fun.lons()[1:imid]
    par[:,0:imid] = fun.data[:,imid:fun.nlon]
    par[:,imid:fun.nlon-1] = fun.data[:,1:imid]

    if(contour):
        plt.contourf(lon,fun.lats(),par,cmap=cstring,levels = ncont)
    else:
        plt.pcolormesh(lon,fun.lats(),par,shading='gouraud',cmap=cstring)
    cbar = plt.colorbar(location = "bottom",fraction=0.072, pad=0.04)
    ax.coastlines()
    cbar.set_label(label,labelpad = 10)
    if(marker != None):
        lat = marker[0]
        lon = marker[1]
        plt.plot([lon],[lat],marker='o', markersize=10, color="green")
    cm = np.nanmax(np.abs(fun.data))
    if(clim != None):
        plt.clim(clim)
    else:
        if(clim_sym):
            clim = [-cm,cm]
        if(clim_pos):
            clim = [0,cm]
        if(clim_neg):
            clim = [-cm,0]
        if(clim_scale != None):
            c1 = clim_scale*clim[0]
            c2 = clim_scale*clim[1]
            clim = [c1,c2]
        plt.clim(clim)

    if(xlim != None):
        plt.xlim(xlim)

    if(ylim != None):
        plt.ylim(ylim)

    if(ofile == None):
        plt.show()
    else:
        plt.savefig(ofile,bbox_inches='tight')
        plt.close()
    return

