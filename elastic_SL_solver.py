#
# Collection of routines based on Al-Attar et al. 2023 to solve the elastic sea level finger print
# This code is part of a bigger package to perform sea level fingerprints for NASA's ModelE Earth system model
#

import pyshtools as pysh
import numpy as np
from numpy import pi as pi
from scipy.interpolate import RegularGridInterpolator
import my_io as io
import fp_util as fpu
sentinel = object()


if __name__ == "__main__":
    pass


############################ Physical Constants ##########################
b    = 6368000.          # mean radius of the solid Earth 
g    = 9.825652323       # mean surface gravity
G    = 6.6723e-11        # gravitational constant
rhofw = 1000.0           # density of fresh water
rhoow = 1030.0           # density of ocean water
rhoi =  917.0            # density of ice
rhos = 2600.0            # surface density of solid Earth
CC = 8.038e37            # polar moment of inertia 
AA = 8.012e37            # equatorial moment of inertia
Om = 2.*pi/(24.*3600.)   # mean rotation rate
Me = 5.974e24            # mass of the Earth
#########################################################################


########################## Numerical Constants ##########################
ep = 1.e-8             # tolerance for iterations
#########################################################################


############################ Ocean Function #############################
def ocean_function(sl0,ice0):
  # function that returns the ocean function and ocean area

  # data for Caspian sea
  lat_C = 42.2
  lon_C = 50.7+180.0
  th_C = (90-lat_C)*pi/180
  ph_C= (180+lon_C)*pi/180
  delta_C = 7.5*pi/180

  # data for Lake Eyre
  lat_E = -29.2
  lon_E = 137.3+180.0
  th_E = (90-lat_E)*pi/180
  ph_E= (180+lon_E)*pi/180
  delta_E = 2.0*pi/180
        
  C = sl0.copy()        
  for ilat,lat in enumerate(C.lats()):
    for ilon,lon in enumerate(C.lons()):
      sll  = sl0.data[ilat,ilon] 
      icel = ice0.data[ilat,ilon]    
      if (rhoow*sll - rhoi*icel >= 0.): 
        C.data[ilat,ilon] = 1.
      else:
        C.data[ilat,ilon] = 0.
    
      # remove the Caspian sea
      th = (90-lat)*pi/180
      ph = lon*pi/180
      delta = np.cos(th)*np.cos(th_C) + np.sin(th)*np.sin(th_C)*np.cos(ph-ph_C)
      delta = np.arccos(delta)
#      if (delta < delta_C):
#        C.data[ilat,ilon] = 0.0

      # remove the Lake Eyre
      delta = np.cos(th)*np.cos(th_E) + np.sin(th)*np.sin(th_E)*np.cos(ph-ph_E)
      delta = np.arccos(delta)
#      if (delta < delta_E):
#        C.data[ilat,ilon] = 0.0

  return C


############################### Rotation ################################
def inertia_tensor_perturbation(phi_lm):
  # returns the jx = Om*i_{zx} and jy = Om*i_{zy} components of the
  # torque perturbation associated with the gravitational potential

  j = Om*np.sqrt(5./(12*pi))*(b**3/G)*phi_lm.coeffs[:,2,1]
    
  return j


def rotation_vector_perturbation(j):
  # returns the rotation vector perturbations given those for the
  # torque vector, j.

  om = j/(CC-AA)
  
  return om
        
        
def centrifugal_perturbation_coefficients(om):
  # returns the centrifugal potential perturbation in spherical harmonic
  # domain given the rotation vector perturbation

  psi_2m = np.zeros([2,3])
  psi_2m[:,1] = b**2*Om*np.sqrt((4*pi)/15.)*om
  
  return psi_2m

############################# Fingerprints ##############################
#
# function to solve the fingerprint problem for a given direct load
#

def fingerprint(C,zeta,rotation=True):
  
  ## Gather needed information
  L = C.lmax # get the maximum degree
  A = fpu.surface_integral(C) # Ocean area
  h,k,ht,kt = io.love_numbers(L) # get the love numbers

  ## Calculate the average change in sea level
  slu = -fpu.surface_integral(zeta)/(rhoow*A)
  onegrid = pysh.SHGrid.from_zeros(lmax=L,grid = 'GLQ')
  onegrid.data[:,:] = 1.
  sl = slu*onegrid

  ## Initialise displacement and potential perturbations
  u_lm   = pysh.SHCoeffs.from_zeros(lmax=sl.lmax,normalization = 'ortho')
  phi_lm = u_lm.copy()
  psi_lm = u_lm.copy()

  ## Prepare to solve iteratively
  sl0 = sl.copy()
  err = 1.
  it = -1

  ## loop to solve
  while (err > ep):
    # Compute the current loads
    sigma = rhoow*C*sl + zeta # ocean + surface load
    sigma_lm  = sigma.expand(normalization = 'ortho')

    # Determine the response to the loading
    for l in range(L+1):
      u_lm.coeffs[:,l,:]   =  h[l]*sigma_lm.coeffs[:,l,:]
      phi_lm.coeffs[:,l,:] =  k[l]*sigma_lm.coeffs[:,l,:]
   
    # Add in the centrifugal contribution
    u_lm.coeffs[:,2,:]   +=  ht*psi_lm.coeffs[:,2,:]
    phi_lm.coeffs[:,2,:] +=  kt*psi_lm.coeffs[:,2,:]

    # get the centrifugal potential perturbation
    if (rotation):
      j  = inertia_tensor_perturbation(phi_lm)
      om = rotation_vector_perturbation(j)
      psi_2m = centrifugal_perturbation_coefficients(om)
      psi_lm.coeffs[:,2,:3] = psi_2m

    # get the spatial fields
    u   =   u_lm.expand(grid = 'GLQ')
    phi = phi_lm.expand(grid = 'GLQ')
    psi = psi_lm.expand(grid = 'GLQ')


    ## Update the sea level
    fac = fpu.ocean_integral(C,g*u + phi + psi)/(g*A) + slu
    sl = -1.*(g*u + phi + psi)/g + fac*onegrid

    it = it+1
    if (it == 0):
      slnorm = np.max(np.abs(sl.data))
    else:
      err = np.max(np.abs(sl.data - sl0.data))/np.abs(slnorm)
      print('iteration = ',it,'relative change = ',err)

    # store the most recent solution
    sl0 = sl.copy()

  if (not rotation):
    om = np.zeros(2)


  return sl,u,phi,om,psi

