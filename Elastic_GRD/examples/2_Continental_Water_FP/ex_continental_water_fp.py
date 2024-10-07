#
# An example of an elastic SL FP for lakes in NA being increased by 1 meter
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

# calculate dL



# plot inputs
plt.plot(sl0_glq,label = 'Initial Sea Level (m)',ofile = fpath+'/sl0.png')
plt.plot(I0_glq,label = 'Initial Ice thickness (m)',ofile = fpath+'/I0.png')
plt.plot(dI_glq,label = 'Change in Ice thickness (m)',ofile = fpath+'/dI.png')

# compute ocean function and load
C = ess.ocean_function(sl0_glq,I0_glq)
zeta = 

plt.plot(C,label = 'Ocean Function (1 or 0)',ofile = fpath+'/ocean_func.png',clim = [0,1])
plt.plot(zeta,label = 'load Chage (kg/m^2)',ofile = fpath+'/load_change.png')

# Compute fingerprint
sl,u,phi,om,psi = ess.fingerprint(C,zeta,rotation=True)

# Plot fingerprint
plt.plot(sl,label = 'RSL Change (m)',clim=[-0.01,0.01],ofile = fpath+'/rsl_change.png')
plt.plot(u,label = 'u (m)',clim=[-0.01,0.01],ofile = fpath+'/u.png')
plt.plot(phi,label = 'phi Change (m^2/s^2)',clim=[-0.01,0.01],ofile = fpath+'/phi_change.png')

