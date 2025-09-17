import numpy as np

# Any initialization step can be done here
# This method is called during CldFraction::initialize_impl
def init ():
    pass

#########################################################
def main (ice_threshold, ice_4out_threshold,
          qi, liq_cld_frac,
          ice_cld_frac, tot_cld_frac,
          ice_cld_frac_4out, tot_cld_frac_4out):
#########################################################

    ice_cld_frac[:] = 0
    ice_cld_frac_4out[:] = 0
    ice_cld_frac[qi > ice_threshold] = 1
    ice_cld_frac_4out[qi > ice_4out_threshold] = 1

    np.maximum(ice_cld_frac,liq_cld_frac, out=tot_cld_frac)
    np.maximum(ice_cld_frac_4out,liq_cld_frac,out=tot_cld_frac_4out)
