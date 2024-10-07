#
# A collection of function to aid in calculating SL Fingerprints
#
#

import numpy as np
import pyshtools as pysh
from numpy import pi as pi
import elastic_SL_solver as ess

########################### General Utilities  #########################
#
# 
#

def surface_integral(fun):
  # integrates a function over the surface

  b = ess.b

  fun_lm = fun.expand(lmax_calc = 0,normalization = 'ortho')
  int = np.sqrt(4*pi)*b*b*fun_lm.coeffs[0,0,0]

  return int

def ocean_integral(C,fun):
  # integrates a function over the oceans

  b = ess.b

  tmp = C*fun
  tmp_lm = tmp.expand(lmax_calc = 0,normalization = 'ortho')
  int = np.sqrt(4*pi)*b*b*tmp_lm.coeffs[0,0,0]
    
  return int
