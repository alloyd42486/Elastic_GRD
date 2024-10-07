#
# A collection of functions for the use of adjoint methods in the elastic sea level
# fingerprints
#

import numpy as np
import pyshtools as pysh
from numpy import pi as pi
import elastic_SL_solver as ess
sentinel = object()


############################# Generic Loads  ###########################

def point_load(L,lats,lons,grid = 'GLQ',angle = 0.,w=[sentinel]):
  # returns a point load at a given geographic location with optional
  # inverse Laplacian smoothing

  if (len(lats)!=len(lons)):
    raise SystemExit('lats and lons are different sizes!')

  if (w[0]==sentinel):
    w = np.ones(len(lats))

  th = 0.4*angle*pi/180
  t  = th*th

  ths = 90.- lats # colatitude
  phs = lons

  pl_lm = pysh.SHCoeffs.from_zeros(lmax=L,normalization = 'ortho')

  for isource in range(len(lats)):

    ylm = pysh.expand.spharm(pl_lm.lmax,ths[isource],phs[isource],normalization = 'ortho')

    for l in range(0,pl_lm.lmax+1):
      fac = np.exp(-l*(l+1)*t)
      pl_lm.coeffs[0,l,0] +=  w[isource]*ylm[0,l,0]*fac
      
      for m in range(1,l+1):
        pl_lm.coeffs[0,l,m] += w[isource]*ylm[0,l,m]*fac
        pl_lm.coeffs[1,l,m] += w[isource]*ylm[1,l,m]*fac
    
  pl_lm = (1/ess.b**2)*pl_lm
  pl = pl_lm.expand(grid = 'GLQ')

  return pl



############################# Adjoint Loads  ###########################
#
# 
#

def sea_level_load(L,lat,lon,grid = 'GLQ',angle = 1.):
  # returns the adjoint loads for a sea level measurement at 

  lats = np.full((1),lat)
  lons = np.full((1),lon)

  zeta     = point_load(L,lats,lons,angle = angle,grid = grid) # sea level load
  zeta_u   = pysh.SHGrid.from_zeros(lmax=L,grid = grid) # solid Earth load
  zeta_phi = pysh.SHGrid.from_zeros(lmax=L,grid = grid) # geoid load
  kk       = np.zeros(2)

  return zeta,zeta_u,zeta_phi,kk
