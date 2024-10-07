#
# Compute benchmark for WAIS against SLcode from J. Austermann
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

fpath = 'Figures_WAIS'

# set the truncation degree
L=512


##### Get inputs and calculate some needed quantities
# read in the present day topography and determine SL
topo_glq = io.get_topo('data/topography_matlab_benchmark.txt',L)
sl0_glq = -1.0*topo_glq

# read in initial ice thickness
I0_glq = io.get_generic_shell_field('data/ice_0_nointerp.txt',L)

# calculate dI
dI_glq = -1.0*io.get_generic_shell_field('data/WAIS_nointerp_del_ice.txt',L)

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
plt.plot(sl,label = 'RSL Change (m)',ofile = fpath+'/rsl_change.png')
plt.plot(u,label = 'u (m)',ofile = fpath+'/u.png')
plt.plot(phi,label = 'phi Change (m^2/s^2)',ofile = fpath+'/phi_change.png')



# Load SLcode result
sl_SLcode = io.get_generic_shell_field('data/WAIS_dsl_512.txt',L)

plt.plot(sl_SLcode,label = 'RSL Change (m)',ofile = fpath+'/SL_code_rsl_change.png')
plt.plot(sl-sl_SLcode,label = 'Benchmark Residual RSL (m)',ofile = fpath+'/BM_residual_rsl_change.png')
