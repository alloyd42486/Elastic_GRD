#
# A collection of subroutines for handing io for different surface fields that are to be
# represented on a GLQ grid in pyshtools
#
#

import numpy as np
import pyshtools as pysh
from scipy.interpolate import RegularGridInterpolator
import os
from pathlib import Path

########################### Generic Surface Field  #######################
#
# Generic spherical field input
# 3 columns longitude, latitude, parameter
#

def load_generic_shell_field(file,**kwargs):
  # we can specify the column of the parameter with par_column
  # lon and lat are assumed to be columns 1,2 respectively (i.e. index 0,1)

  ## optional argument
  par_column = kwargs.get('par_column',3)

  ## load input data file
  mf1=Path(file+'.npy')
  mf2=Path(file)

  if mf1.is_file():
    data = np.load(mf1)

  elif mf2.is_file():
    data = np.loadtxt(mf2) #lon, lat, thickness change

    #numpy .npy binary format example
    np.save(mf1,data)

  else:
    print("Can not find file: "+mf2)


  ## reshape
  # unique returns values in  assending order
  lon = np.unique(data[:,0])
  lat = np.unique(data[:,1])

  num1 = lon.size #columns
  num2 = lat.size # row
  ip = par_column - 1


  # TODO - code other cases
  # latitude changes rapidly in assending order
  if (data[0,0] == data[1,0] and data[0,1] < data[1,1]):
    # longitude in assending order
    if (data[0,0] < data[-1,0]):
      par = np.flipud(np.transpose(np.reshape(data[:,ip],(num1,num2))))
      lat = np.flip(lat)
    # longitude in decending order
    else:
      print("finish storage code")
      exit()

  # longitude changes rapidly in assending order
  elif (data[0,0] < data[1,0] and data[0,1] == data[1,1]):
    # latitude in assending order
    if (data[0,1] < data[-1,1]):
      par = np.flipud(np.reshape(data[:,ip],(num2,num1)))
      lat = np.flip(lat)
    # latitude in decending order
    else:
      par = np.reshape(data[:,ip],(num2,num1))
      lat = np.flip(lat)
  else:
    print("Error with format of lat and lon in "+file)
    exit()

  ## Correction for -180, 180 to 0,360 
  if (not(any(x > 180.0 for x in lon))):
    lon = np.flip(np.mod(lon+360.0,360.0))
    par = np.fliplr(par)

  ## correction for insufficient grid coverage
  if ( (lat[0] < 90.0 and lat[0] > 0.0) or (lat[0] > -90.0 and lat[0] < 0.0) ):
    print('Expanding latitude grid in 0 position: '+file) 
    lat = np.append([np.sign(lat[0])*90.0],lat)
    par = np.vstack((np.average(par[0,:])*np.ones(num1),par))   

  if ( (lat[-1] < 90.0 and lat[-1] > 0.0) or (lat[-1] > -90.0 and lat[-1] < 0.0) ):
    print('Expanding latitude grid in -1 position: '+file)
    lat = np.append(lat,[np.sign(lat[-1])*90.0])
    par = np.vstack((par,np.average(par[-1,:])*np.ones(num1)))

  if ( (0.0 < lon[0] and 0.0 < lon[-1]) or ( 360.0 > lon[0] and 360.0 > lon[-1]) ):
    print('Expanding longitude grid: '+file)
    # calculating difference between points 0 and -1; likely only works for 0 to 360
    lon_0 = lon[0] - (lon[0]+180)%360.0 + (lon[-1]+180)%360.0
    lon_n = lon[-1] - (lon[-1]+180)%360.0 + (lon[0]+180)%360.0
    lon = np.concatenate(([lon_0],lon,[lon_n]))
    par = np.column_stack((par[:,-1],par,par[:,0]))


  return lon,lat,par


########################### Initial States ##############################
#
# function to load sea level (i.e. -Topography)
# and store on needed GLQ grid
#

def get_sl(file,L):

  ## load data
  lon,lat,sl = load_generic_shell_field(file)

  ## build interpolation
  sl_interp = RegularGridInterpolator((lat,lon), sl, method = 'linear', bounds_error=True)
 
  ## Perform interpolation
  # initialize grid
  sl_glq = pysh.SHGrid.from_zeros(lmax=L,grid = 'GLQ')
  lon2_glq,lat2_glq = np.meshgrid(sl_glq.lons(),sl_glq.lats())

  sl_glq.data[:,:] = sl_interp((lat2_glq,lon2_glq))

  return sl_glq

def get_topo(file,L):

  ## load data
  lon,lat,topo = load_generic_shell_field(file)

  ## build interpolation
  topo_interp = RegularGridInterpolator((lat,lon), topo, method = 'linear', bounds_error=True)

  ## Perform interpolation
  # initialize grid
  topo_glq = pysh.SHGrid.from_zeros(lmax=L,grid = 'GLQ')
  lon2_glq,lat2_glq = np.meshgrid(topo_glq.lons(),topo_glq.lats())

  topo_glq.data[:,:] = topo_interp((lat2_glq,lon2_glq))

  return topo_glq


########################### surface loads  ##############################
#
# function to load surface load returning thickness
#

def get_generic_shell_field(file,L):

  ## load data
  lon,lat,f = load_generic_shell_field(file)

  ## Perform interpolation
  # initialize grid
  f_glq = pysh.SHGrid.from_zeros(lmax=L,grid = 'GLQ')
  lon2_glq,lat2_glq = np.meshgrid(f_glq.lons(),f_glq.lats())

  ## build interpolation
  f_interp = RegularGridInterpolator((lat,lon), f, method = 'linear', bounds_error=True)

  # interpolate onto new grid
  f_glq.data[:,:] = f_interp((lat2_glq,lon2_glq))

  return f_glq


#
# function to load surface load returning mass
#


########################### Love Numbers  ###############################
#
# Load love numbers obtained following Al-Attar and Trump 2014
#
def love_numbers(L):

  # load input data file
  dirname = os.path.dirname(__file__)

  f1 = os.path.join(dirname, 'data/love.npy')
  f2 = os.path.join(dirname, 'data/love.dat')

  mf1=Path(f1)
  mf2=Path(f2)

  if mf1.is_file():
    data = np.load(mf1)

  elif mf2.is_file():
    data = np.loadtxt(mf2)

    #numpy .npy binary format example
    np.save(mf1,data)

  else:
    print("Can not find file: "+mf2)

  ## store data
  h = data[:L+1,1] + data[:L+1,3]
  k = data[:L+1,2] + data[:L+1,4]
  ht = data[2,5]
  kt = data[2,6]

  return h,k,ht,kt
