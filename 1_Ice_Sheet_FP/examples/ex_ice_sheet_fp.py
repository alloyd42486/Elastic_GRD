#
# An example of an elastic SL FP for Greenland comparing rotation on and offs
#

# import libraries
import numpy as np
import pyshtools as pysh

import sys
sys.path.append('../../')
import my_io as io
import plot as plt
import elastic_SL_solver as ess
import fp_util as fpu

fpath = 'Figures'

# set the truncation degree
L=512


##### Get inputs and calculate some needed quantities
# read in the present day topography and determine SL
topo_glq = io.get_topo('../../data/etopo_bed.llz',L)
sl0_glq = -1.0*topo_glq

# read in initial ice thickness
I0_glq = io.get_generic_shell_field('../../data/IT_ANT-BM2_GR-BM5.llt',L)

# calculate dI
dI_glq = I0_glq.copy()
dI_glq.data[:,:] = 0.0

for i,lat in enumerate(dI_glq.lats()):
  for j,lon in enumerate(dI_glq.lons()):
    if (I0_glq.data[i,j] > 0.0):
      if ( 290.0 < lon < 350.0 and 58.0 < lat < 90.0):
        dI_glq.data[i,j] = -1.0

# plot inputs
plt.plot(sl0_glq,label = 'Initial Sea Level (m)',ofile = fpath+'/sl0.png')
plt.plot(I0_glq,label = 'Initial Ice thickness (m)',ofile = fpath+'/I0.png')
plt.plot(dI_glq,label = 'Change in Ice thickness (m)',ofile = fpath+'/dI.png')

# compute ocean function and load
C = ess.ocean_function(sl0_glq,I0_glq)
zeta = ess.rhoi*dI_glq*(1-C)

plt.plot(C,label = 'Ocean Function (1 or 0)',ofile = fpath+'/ocean_func.png',clim = [0,1])
plt.plot(zeta,label = 'load Chage (kg/m^2)',ofile = fpath+'/load_change.png')

# Compute fingerprint
sl,u,phi,om,psi = ess.fingerprint(C,zeta,rotation=True)

# Plot fingerprint
plt.plot(sl,label = 'RSL Change (m)',clim=[-0.01,0.01],ofile = fpath+'/rsl_change.png')
plt.plot(u,label = 'u (m)',clim=[-0.01,0.01],ofile = fpath+'/u.png')
plt.plot(phi,label = 'phi Change (m^2/s^2)',clim=[-0.01,0.01],ofile = fpath+'/phi_change.png')



#### No Rotation
# Compute fingerprint
sl_nr,u_nr,phi_nr,om_nr,psi_nr = ess.fingerprint(C,zeta,rotation=False)

# Plot fingerprint
plt.plot(sl_nr,label = 'NR: RSL Change (m)',clim=[-0.01,0.01],ofile = fpath+'/nr_rsl_change.png')
plt.plot(u_nr,label = 'NR: u (m)',clim=[-0.01,0.01],ofile = fpath+'/nr_u.png')
plt.plot(phi_nr,label = 'NR: phi Change (m^2/s^2)',clim=[-0.01,0.01],ofile = fpath+'/nr_phi_change.png')

# Plot fingerprint
plt.plot(sl-sl_nr,label = 'RMNR: RSL Change (m)',clim=[-0.001,0.001],ofile = fpath+'/rMnr_rsl_change.png')
plt.plot(u-u_nr,label = 'RMNR: u (m)',clim=[-0.001,0.001],ofile = fpath+'/rMnr_u.png')
plt.plot(phi-phi_nr,label = 'RMNR: phi Change (m^2/s^2)',clim=[-0.001,0.001],ofile = fpath+'/rMnr_phi_change.png')
