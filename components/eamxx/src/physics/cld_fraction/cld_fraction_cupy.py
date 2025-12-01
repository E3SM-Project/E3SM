import cupy as cp

# Any initialization step can be done here
# This method is called during CldFraction::initialize_impl
def init ():
    pass

#########################################################
def get_cu_array(np_arr):
#########################################################
    shape   = np_arr.shape
    dtype   = np_arr.dtype
    ptr     = np_arr.__array_interface__['data'][0]
    strides = np_arr.strides

    # The exact size here does not really matter, as we are just creating an
    # unmanaged mem block, of which we then simply grab the start address.
    # Still, use the correct size for code clarity
    size = shape[0]*strides[0]
    mem = cp.cuda.UnownedMemory(ptr=ptr,owner=None,size=size)
    memptr = cp.cuda.MemoryPointer(mem, 0)

    return cp.ndarray(shape=shape,dtype=dtype,memptr=memptr,strides=strides)

#########################################################
def main (ice_threshold, ice_4out_threshold,
          qi, liq_cld_frac,
          ice_cld_frac, tot_cld_frac,
          ice_cld_frac_4out, tot_cld_frac_4out):
#########################################################

    cu_qi = get_cu_array(qi)
    cu_liq_cld_frac = get_cu_array(liq_cld_frac)
    cu_ice_cld_frac = get_cu_array(ice_cld_frac)
    cu_tot_cld_frac = get_cu_array(tot_cld_frac)
    cu_ice_cld_frac_4out = get_cu_array(ice_cld_frac_4out)
    cu_tot_cld_frac_4out = get_cu_array(tot_cld_frac_4out)

    cu_ice_cld_frac[:] = 0
    cu_ice_cld_frac_4out[:] = 0
    cu_ice_cld_frac[cu_qi > ice_threshold] = 1
    cu_ice_cld_frac_4out[cu_qi > ice_4out_threshold] = 1

    cp.maximum(cu_ice_cld_frac,cu_liq_cld_frac, out=cu_tot_cld_frac)
    cp.maximum(cu_ice_cld_frac_4out,cu_liq_cld_frac,out=cu_tot_cld_frac_4out)
