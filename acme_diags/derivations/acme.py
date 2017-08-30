from __future__ import print_function

import copy
from collections import OrderedDict
from numbers import Number
import cdms2
from genutil import udunits
import MV2
import numpy as np


def process_derived_var(var_key, derived_vars_dict, nc_file, parameter):
    """Given a key (var_key) to the derived_vars_dict dict, compute and return
     whatever is described in derived_vars_dict[var_key] for the nc_file"""
    if hasattr(parameter, 'derived_variables'):
        _add_user_derived_vars(derived_vars_dict, parameter)
    if var_key in derived_vars_dict:
        return _compute_derived_var(var_key, derived_vars_dict, nc_file)
    else:
        raise RuntimeError(
            'The variable %s was not in the derived variables dictionary' % var_key)


def _add_user_derived_vars(derived_vars_dict, parameter):
    """Append parameter.derived_variables to the correct part of derived_vars_dict"""
    for key, user_derived_vars in parameter.derived_variables.iteritems():
        # append the user-defined vars to the already defined ones
        # add to an existing entry, otherwise create a new one
        if key in derived_vars_dict:
            new_dict = OrderedDict(user_derived_vars)  # has user-defined derived vars first
            # add all of the default derived vars to the end of new_dict
            for k in derived_vars_dict[key]:
                if k in new_dict:  # don't overwrite the user-defined var with a default derived var
                    continue
                new_dict[k] = derived_vars_dict[key][k]
            derived_vars_dict[key] = new_dict
        else:
            derived_vars_dict[key] = user_derived_vars


def _compute_derived_var(var_key, derived_vars_dict, nc_file):
    """Call the first valid derivation from the derived_vars_dict dict."""
    derived_vars = derived_vars_dict[var_key]
    # store a list of all inputs visited, so if Exception, we get a good msg.
    derived_var_inputs = []

    # get the first function and inputs from the derived_vars_dict dict
    for inputs, func in derived_vars.iteritems():
        derived_var_inputs.append(inputs)
        # tuples with a single string [ex: ('pr')] become just a string ['pr']
        # are all of the variables (inputs) in the nc_file?
        if isinstance(inputs, str) and inputs in nc_file.variables:
            args = [nc_file(inputs)(squeeze=1)]
            return func(*args)

        elif isinstance(inputs, tuple) and set(inputs).issubset(nc_file.variables):
            args = [nc_file(var)(squeeze=1) for var in inputs]
            return func(*args)

    # When nc_file is obs, there is var_key in nc_file, i.e. nc_file(var_key)
    if var_key in nc_file.variables:
        return nc_file(var_key)(squeeze=1)

    raise RuntimeError('None of the variables (%s) are in the file: %s' % (
        derived_var_inputs, nc_file.id))


def rename(new_name):
    """Given the new name, just return it."""
    return new_name


def aplusb(var1, var2, target_units=None):
    """Returns var1 + var2. If both of their units are not the same,
    it tries to convert both of their units to target_units"""

    if target_units is not None:
        var1 = convert_units(var1, target_units)
        var2 = convert_units(var2, target_units)

    return var1 + var2


def convert_units(var, target_units):
    """Converts units of var to target_units.
    var is a cdms.TransientVariable."""

    if not hasattr(var, 'units') and var.id == 'SST':
        var.units = target_units
    elif not hasattr(var, 'units') and var.id =='ICEFRAC':
        var.units = target_units
        var = 100.0 * var
    elif var.units == 'fraction':
        var = 100.0 * var
        var.units = target_units
    elif var.units == 'mb':
        var.units = target_units
    elif var.units == 'gpm':  # geopotential meter
        var = var / 9.8 / 100  # convert to hecto meter
        var.units = target_units
    elif var.units == 'Pa/s':
         var = var /100.0*24*3600
         var.units = target_units
    elif var.units == 'mb/day':
         var = var 
         var.units = target_units
    else:
        temp = udunits(1.0, var.units)
        coeff, offset = temp.how(target_units)
        var = coeff * var + offset
        var.units = target_units

    return var


def mask_by(input_var, maskvar, low_limit=None, high_limit=None):
    """masks a variable var to be missing except where maskvar>=low_limit and maskvar<=high_limit. 
    None means to omit the constrint, i.e. low_limit = -infinity or high_limit = infinity. var is changed and returned; we don't make a new variable.
    var and maskvar: dimensioned the same variables.
    low_limit and high_limit: scalars.
    """
    var = copy.deepcopy(input_var)
    if low_limit is None and high_limit is None:
        return var
    if low_limit is None and high_limit is not None:
        maskvarmask = maskvar > high_limit
    elif low_limit is not None and high_limit is None:
        maskvarmask = maskvar < low_limit
    else:
        maskvarmask = (maskvar < low_limit) | (maskvar > high_limit)
    if var.mask is False:
        newmask = maskvarmask
    else:
        newmask = var.mask | maskvarmask
    var.mask = newmask
    return var


def qflxconvert_units(var):
    print(var.units)
    if var.units == 'kg/m2/s' or var.units == 'kg m-2 s-1':
        # need to find a solution for units not included in udunits
        # var = convert_units( var, 'kg/m2/s' )
        var = var * 3600.0 * 24  # convert to mm/day
        var.units = 'mm/day'
    return var

def prect(precc, precl):
    """Total precipitation flux = convective + large-scale"""
    var = precc + precl
    var = convert_units(var, "mm/day")
    var.long_name = "Total precipitation rate (convective + large-scale)"
    return var


def albedo(rsdt, rsut):
    """TOA (top-of-atmosphere) albedo, rsut / rsdt, unit is nondimension"""
    var = rsut / rsdt
    var.units = "dimensionless"
    var.long_name = "TOA albedo"
    return var

def albedoc(rsdt, rsutcs):
    """TOA (top-of-atmosphere) albedo clear-sky, rsutcs / rsdt, unit is nondimension"""
    var = rsutcs / rsdt
    var.units = "dimensionless"
    var.long_name = "TOA albedo clear-sky"
    return var

def albedo_srf(rsds, rsus):
    """Surface albedo, rsus / rsds, unit is nondimension"""
    var = rsus / rsds
    var.units = "dimensionless"
    var.long_name = "Surface albedo"
    return var

def rst(rsdt, rsut):
    """TOA (top-of-atmosphere) net shortwave flux"""
    var = rsdt - rsut
    var.long_name = "TOA net shortwave flux"
    return var

def rstcs(rsdt, rsutcs):
    """TOA (top-of-atmosphere) net shortwave flux clear-sky"""
    var = rsdt - rsutcs
    var.long_name = "TOA net shortwave flux clear-sky"
    return var

def swcfsrf(fsns, fsnsc):
    """Surface shortwave cloud forcing """
    var = fsns - fsnsc
    var.long_name = "Surface shortwave cloud forcing"
    return var

def lwcfsrf(flns, flnsc):
    """Surface longwave cloud forcing, for ACME model, upward is postitive for LW , for ceres, downward is postive for both LW and SW"""
    var = -(flns - flnsc)
    var.long_name = "Surface longwave cloud forcing"
    return var

def swcf(fsntoa, fsntoac):
    """TOA shortwave cloud forcing """
    var = fsntoa - fsntoac
    var.long_name = "TOA shortwave cloud forcing"
    return var

def lwcf(flntoa, flntoac):
    """TOA longwave cloud forcing """
    var = flntoa - flntoac
    var.long_name = "TOA longwave cloud forcing"
    return var

def fldsc(ts, flnsc):
    """Clearsky Surf LW downwelling flux"""
    var = 5.67e-8 * ts**4 - flnsc
    var.units = "W/m2"
    var.long_name = "Clearsky Surf LW downwelling flux"
    return var

def restom(fsnt, flnt):
    """TOM(top of model) Radiative flux"""
    var = fsnt - flnt
    var.long_name = "TOM(top of model) Radiative flux"
    return var

def restoa(fsnt, flnt):
    """TOA(top of atmosphere) Radiative flux"""
    var = fsnt - flnt
    var.long_name = "TOA(top of atmosphere) Radiative flux"
    return var

def cosp_bin_sum(cld, prs_low, prs_high, tau_low, tau_high):
    """sum of cosp bins to calculate cloud fraction in specified cloud top pressure and cloud thickness bins, input variable has dimention (cosp_prs,cosp_tau,lat,lon)"""
    prs = cld.getAxis(0)
    tau = cld.getAxis(1)
    if prs_low is None:  prs_low = prs[0]
    if prs_high is None: prs_high = prs[-1]
    if tau_low is None:  tau_low = tau[0]
    if tau_high is None: tau_high = tau[-1]

    if cld.id == 'FISCCP1_COSP':   #ISCCP model
        cld_bin = cld(cosp_prs = (prs_low, prs_high), cosp_tau = (tau_low, tau_high))
    if cld.id == 'CLISCCP':        #ISCCP obs
        cld_bin = cld(isccp_prs = (prs_low, prs_high), isccp_tau = (tau_low, tau_high))

    if cld.id == 'CLMODIS':   #MODIS 
        try: 
            cld_bin = cld(cosp_prs = (prs_low, prs_high), cosp_tau_modis = (tau_low, tau_high)) #MODIS model
        except:
            cld_bin = cld(modis_prs = (prs_low, prs_high), modis_tau = (tau_low, tau_high))  #MODIS obs

    if cld.id == 'CLD_MISR':        #MISR model
        cld_bin = cld(cosp_htmisr = (prs_low, prs_high), cosp_tau = (tau_low, tau_high))
    if cld.id == 'CLMISR':        #MISR obs
        cld_bin = cld(misr_cth = (prs_low, prs_high), misr_tau = (tau_low, tau_high))

    cld_bin_sum = MV2.sum(MV2.sum(cld_bin,axis=1),axis=0)
    return cld_bin_sum

def cosp_histogram_standardize(cld):
    """standarize cloud top pressure and cloud thickness bins to dimensions that suitable for plotting, input variable has dimention (cosp_prs,cosp_tau)"""
    prs = cld.getAxis(0)
    tau = cld.getAxis(1)

    prs_low = prs[0]
    prs_high = prs[-1]
    tau_low = tau[0]
    tau_high = tau[-1]

    prs_bounds = getattr(prs,'bounds')
    if prs_bounds is None:
        cloud_prs_bounds = np.array([[1000.,800.],[800.,680.],[680.,560.],[560.,440.],[440.,310.],[310.,180.],[180.,0.]])  # length 7
        prs.setBounds( np.array(cloud_prs_bounds,dtype=np.float32) )

    tau_bounds =  getattr(tau,'bounds')
    if tau_bounds is None:
        cloud_tau_bounds = np.array([[0.3,1.3],[1.3,3.6],[3.6,9.4],[9.4,23],[23,60],[60,379]]) # length 6
        tau.setBounds( np.array(cloud_tau_bounds,dtype=np.float32) )

    if cld.id == 'FISCCP1_COSP':   #ISCCP model
        cld_hist = cld(cosp_tau = (0.3, tau_high))
    if cld.id == 'CLISCCP':        #ISCCP obs
        cld_hist = cld(isccp_tau = (0.3, tau_high))

    if cld.id == 'CLMODIS':   #MODIS 
        try: 
            cld_hist = cld( cosp_tau_modis = (0.3, tau_high)) #MODIS model
        except:
            cld_hist = cld( modis_tau = (0.3, tau_high))  #MODIS obs

    if cld.id == 'CLD_MISR':        #MISR model
        cld_hist = cld(cosp_tau = (0.3, tau_high), cosp_htmisr=(0,prs_high))
    if cld.id == 'CLMISR':        #MISR obs
        cld_hist = cld(misr_tau = (0.3, tau_high),misr_cth =(0,prs_high))

    return cld_hist
       

# derived_variables is a dictionary to accomodate user specified derived variables. 
# The driver search for available variable keys and functions to calculate derived variable.
# For example 
# In derived_variable there is an entry for 'PRECT':
# If 'PRECT' is not available, but both 'PRECC' and 'PRECT' are available in the netcdf variable keys,
# PRECT is calculated using fuction prect() with precc and precl as inputs.

derived_variables = {
    'PRECT': OrderedDict([
        (('pr'), lambda pr:  qflxconvert_units(rename(pr))),
        (('PRECC', 'PRECL'), lambda precc, precl: prect(precc, precl))
    ]),
    'SST': OrderedDict([
        (('sst'),rename),# lambda sst: convert_units(rename(sst),target_units="degC")),
        (('SST'), lambda sst: convert_units(sst, target_units="degC")),
        (('TS', 'OCNFRAC'), lambda ts, ocnfrac: mask_by(
            convert_units(ts, target_units="degC"), ocnfrac, low_limit=0.9))
    ]),
    'PREH2O': OrderedDict([
        (('TMQ'), rename),
        (('prw'), rename)
    ]),
    'SOLIN': OrderedDict([
        (('rsdt'), rename)
    ]),
    'ALBEDO': OrderedDict([
        (('ALBEDO'), rename),
        (('SOLIN', 'FSNTOA'), lambda solin, fsntoa: albedo(solin, solin-fsntoa)),
        (('rsdt', 'rsut'), lambda rsdt, rsut: albedo(rsdt, rsut))
    ]),
    'ALBEDOC': OrderedDict([
        (('ALBEDOC'), rename),
        (('SOLIN', 'FSNTOAC'), lambda solin, fsntoac: albedoc(solin, solin-fsntoac)),
        (('rsdt', 'rsutcs'), lambda rsdt, rsutcs: albedoc(rsdt, rsutcs))
    ]),
    'ALBEDO_SRF': OrderedDict([
        (('ALBEDO_SRF'), rename),
        (('rsds', 'rsus'), lambda rsds, rsus: albedo_srf(rsds, rsus)),
        (('FSDS', 'FSNS'), lambda fsds, fsns: albedo_srf(fsds, fsds-fsns))
    ]),
    #Pay attention to the positive direction of SW and LW fluxes
    'SWCF': OrderedDict([
        (('SWCF'), rename),
        (('toa_net_sw_all_mon','toa_net_sw_clr_mon'), lambda net_all,net_clr: swcf(net_all,net_clr)),
        (('toa_cre_sw_mon'), rename),
        (('FSNTOA', 'FSNTOAC'), lambda fsntoa, fsntoac: swcf(fsntoa, fsntoac))
    ]),
    'SWCFSRF': OrderedDict([
        (('SWCFSRF'), rename),
        (('sfc_net_sw_all_mon','sfc_net_sw_clr_mon'), lambda net_all,net_clr: swcfsrf(net_all,net_clr)),
        (('sfc_cre_net_sw_mon'), rename),
        (('FSNS', 'FSNSC'), lambda fsns, fsnsc: swcfsrf(fsns, fsnsc))
    ]),
    'LWCF': OrderedDict([
        (('LWCF'), rename),
        (('toa_net_lw_all_mon','toa_net_lw_clr_mon'), lambda net_all,net_clr: lwcf(net_clr,net_all)),
        (('toa_cre_lw_mon'), rename),
        (('FLNTOA', 'FLNTOAC'), lambda flntoa, flntoac: lwcf(flntoa, flntoac))
    ]),
    'LWCFSRF': OrderedDict([
        (('LWCFSRF'), rename),
        (('sfc_net_lw_all_mon','sfc_net_lw_clr_mon'), lambda net_all,net_clr: lwcfsrf(net_clr,net_all)),
        (('sfc_cre_net_lw_mon'), rename),
        (('FLNS', 'FLNSC'), lambda flns, flnsc: lwcfsrf(flns, flnsc))
    ]),
    'FLNS': OrderedDict([
        (('sfc_net_lw_all_mon'), lambda sfc_net_lw_all_mon: -sfc_net_lw_all_mon)
    ]),
    'FLNSC': OrderedDict([
        (('sfc_net_lw_clr_mon'), lambda sfc_net_lw_clr_mon: -sfc_net_lw_clr_mon)
    ]),
    'FLDS': OrderedDict([
        (('rlds'), rename)
    ]),
    'FLDSC': OrderedDict([
        (('rldscs'), rename),
        (('TS', 'FLNSC'), lambda ts, flnsc: fldsc(ts, flnsc))
    ]),
    'FSNS': OrderedDict([
        (('sfc_net_sw_all_mon'), rename)
    ]),
    'FSNSC': OrderedDict([
        (('sfc_net_sw_clr_mon'), rename)
    ]),
    'FSDS': OrderedDict([
        (('rsds'), rename)
    ]),
    'FSDSC': OrderedDict([
        (('rsdscs'), rename),
        (('rsdsc'), rename)
    ]),
    'FLUT': OrderedDict([
        (('rlut'), rename)
    ]),
    'FLNT': OrderedDict([
        (('FLNT'), rename)
    ]),
    'FLUTC': OrderedDict([
        (('rlutcs'), rename)
    ]),
    'FSNTOA': OrderedDict([
        (('FSNTOA'), rename),
        (('rsdt', 'rsut'), lambda rsdt, rsut: rst(rsdt, rsut))
    ]),
    'FSNTOAC': OrderedDict([
        # Note: CERES_EBAF data in amwg obs sets misspells "units" as "lunits"
        (('FSNTOAC'), rename),
        (('rsdt', 'rsutcs'), lambda rsdt, rsutcs: rstcs(rsdt, rsutcs))
    ]),
    'RESTOM': OrderedDict([
        (('RESTOA'), rename),
        (('toa_net_all_mon'), rename),
        (('FSNT', 'FLNT'), lambda fsnt, flnt: restom(fsnt, flnt))
    ]),
    'RESTOA': OrderedDict([
        (('RESTOM'), rename),
        (('toa_net_all_mon'), rename),
        (('FSNT', 'FLNT'), lambda fsnt, flnt: restoa(fsnt, flnt))
    ]),
    'TREFHT_LAND': OrderedDict([
        (('TREFHT_LAND'), rename),
        (('TREFHT', 'LANDFRAC'), lambda trefht, landfrac: mask_by(
            convert_units(trefht, target_units="K"), landfrac, low_limit=0.65))
    ]),
    'PRECT_LAND': OrderedDict([
        (('PRECIP_LAND'), rename),
        # 0.5 just to match amwg
        (('PRECC', 'PRECL', 'LANDFRAC'), lambda precc, precl, landfrac: mask_by(
            prect(precc, precl), landfrac, low_limit=0.5))
    ]),
    'Z3': OrderedDict([
        (('zg'), lambda zg: convert_units(rename(zg), target_units="hectometer")),
        (('Z3'), lambda z3: convert_units(z3, target_units="hectometer"))
    ]),
    'PSL': OrderedDict([
        (('PSL'), lambda psl: convert_units(psl, target_units="mbar"))
    ]),
    'T': OrderedDict([
        (('ta'), rename),
        (('T'), lambda t: convert_units(t, target_units="K"))
    ]),
    'U': OrderedDict([
        (('ua'), rename),
        (('U'), lambda u: convert_units(u, target_units="m/s"))
    ]),
    'V': OrderedDict([
        (('va'), rename),
        (('V'), lambda u: convert_units(u, target_units="m/s"))
    ]),
    'TREFHT': OrderedDict([
        (('TREFHT'), lambda t: convert_units(t, target_units="K"))
    ]),
    'TREFHT': OrderedDict([
        (('TREFHT'), lambda t: convert_units(t, target_units="K"))
    ]),
    'QFLX': OrderedDict([
        (('QFLX'), lambda qflx: qflxconvert_units(qflx))
    ]),
    'LHFLX': OrderedDict([
        (('hfls'), rename)
    ]),
    'SHFLX': OrderedDict([
        (('SHFLX'), rename)
    ]),
    'TGCLDLWP_OCN': OrderedDict([
        (('TGCLDLWP_OCEAN'), lambda x: convert_units(x, target_units='g/m^2')),
        (('TGCLDLWP', 'OCNFRAC'), lambda tgcldlwp, ocnfrac: mask_by(convert_units(tgcldlwp, target_units="g/m^2"), ocnfrac, low_limit=0.65))
    ]),
    'PRECT_OCN': OrderedDict([
        (('PRECT_OCEAN'), lambda x: convert_units(x, target_units='mm/day')),
        (('PRECC', 'PRECL', 'OCNFRAC'), lambda a, b, ocnfrac: mask_by(aplusb(a, b, target_units="mm/day"), ocnfrac, low_limit=0.65))
    ]),
    'PREH2O_OCN': OrderedDict([
        (('PREH2O_OCEAN'), lambda x: convert_units(x, target_units='mm')),
        (('TMQ', 'OCNFRAC'), lambda preh2o, ocnfrac: mask_by(preh2o, ocnfrac, low_limit=0.65))
    ]),
    'CLDHGH': OrderedDict([
        (('CLDHGH'), lambda cldhgh: convert_units(cldhgh, target_units="%"))
    ]),
    'CLDLOW': OrderedDict([
        (('CLDLOW'), lambda cldlow: convert_units(cldlow, target_units="%"))
    ]),
    'CLDMED': OrderedDict([
        (('CLDMED'), lambda cldmed: convert_units(cldmed, target_units="%"))
    ]),
    'CLDTOT': OrderedDict([
        (('CLDTOT'), lambda cldtot: convert_units(cldtot, target_units="%"))
    ]),
#below for COSP output
    #CLIPSO
    'CLDHGH_CAL': OrderedDict([
        (('CLDHGH_CAL'), lambda cldhgh: convert_units(cldhgh, target_units="%"))
    ]),
    'CLDLOW_CAL': OrderedDict([
        (('CLDLOW_CAL'), lambda cldlow: convert_units(cldlow, target_units="%"))
    ]),
    'CLDMED_CAL': OrderedDict([
        (('CLDMED_CAL'), lambda cldmed: convert_units(cldmed, target_units="%"))
    ]),
    'CLDTOT_CAL': OrderedDict([
        (('CLDTOT_CAL'), lambda cldtot: convert_units(cldtot, target_units="%"))
    ]),
    #ISCCP
    'CLDTOT_TAU1.3_ISCCP': OrderedDict([
        (('FISCCP1_COSP'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, None), target_units="%")),
        (('CLISCCP'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, None), target_units="%"))
    ]),
    'CLDTOT_TAU1.3_9.4_ISCCP': OrderedDict([
        (('FISCCP1_COSP'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, 9.4), target_units="%")),
        (('CLISCCP'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, 9.4), target_units="%"))
    ]),
    'CLDTOT_TAU9.4_ISCCP': OrderedDict([
        (('FISCCP1_COSP'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 9.4, None), target_units="%")),
        (('CLISCCP'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 9.4, None), target_units="%"))
    ]),
    #MODIS
    'CLDTOT_TAU1.3_MODIS': OrderedDict([
        (('CLMODIS'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, None), target_units="%")),
    ]),
    'CLDTOT_TAU1.3_9.4_MODIS': OrderedDict([
        (('CLMODIS'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, 9.4), target_units="%")),
    ]),
    'CLDTOT_TAU9.4_MODIS': OrderedDict([
        (('CLMODIS'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 9.4, None), target_units="%")),
    ]),
    'CLDHGH_TAU1.3_MODIS': OrderedDict([
        (('CLMODIS'), lambda cld: convert_units(cosp_bin_sum(cld,0, 440, 1.3, None), target_units="%")),
    ]),
    'CLDHGH_TAU1.3_9.4_MODIS': OrderedDict([
        (('CLMODIS'), lambda cld: convert_units(cosp_bin_sum(cld,0, 440, 1.3, 9.4), target_units="%")),
    ]),
    'CLDHGH_TAU9.4_MODIS': OrderedDict([
        (('CLMODIS'), lambda cld: convert_units(cosp_bin_sum(cld,0, 440, 9.4, None), target_units="%")),
    ]),
    #MISR
    'CLDTOT_TAU1.3_MISR': OrderedDict([
        (('CLD_MISR'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, None), target_units="%")),
        (('CLMISR'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, None), target_units="%"))
    ]),
    'CLDTOT_TAU1.3_9.4_MISR': OrderedDict([
        (('CLD_MISR'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, 9.4), target_units="%")),
        (('CLMISR'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 1.3, 9.4), target_units="%"))
    ]),
    'CLDTOT_TAU9.4_MISR': OrderedDict([
        (('CLD_MISR'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 9.4, None), target_units="%")),
        (('CLMISR'), lambda cld: convert_units(cosp_bin_sum(cld,None,None, 9.4, None), target_units="%"))
    ]),
    'CLDLOW_TAU1.3_MISR': OrderedDict([
        (('CLD_MISR'), lambda cld: convert_units(cosp_bin_sum(cld,0, 3, 1.3, None), target_units="%")),
        (('CLMISR'), lambda cld: convert_units(cosp_bin_sum(cld,0, 3, 1.3, None), target_units="%"))
    ]),
    'CLDLOW_TAU1.3_9.4_MISR': OrderedDict([
        (('CLD_MISR'), lambda cld: convert_units(cosp_bin_sum(cld, 0, 3, 1.3, 9.4), target_units="%")),
        (('CLMISR'), lambda cld: convert_units(cosp_bin_sum(cld, 0, 3, 1.3, 9.4), target_units="%"))
    ]),
    'CLDLOW_TAU9.4_MISR': OrderedDict([
        (('CLD_MISR'), lambda cld: convert_units(cosp_bin_sum(cld, 0, 3, 9.4, None), target_units="%")),
        (('CLMISR'), lambda cld: convert_units(cosp_bin_sum(cld, 0, 3, 9.4, None), target_units="%"))
    ]),
# COSP cloud fraction joint histogram     
    'COSP_HISTOGRAM_MISR': OrderedDict([
        (('CLD_MISR'),lambda cld: cosp_histogram_standardize(rename(cld))),
        (('CLMISR'),lambda cld: cosp_histogram_standardize(rename(cld)))
    ]),
    'COSP_HISTOGRAM_MODIS': OrderedDict([
        (('CLMODIS'),lambda cld: cosp_histogram_standardize(rename(cld))),
    ]),
    'COSP_HISTOGRAM_ISCCP': OrderedDict([
        (('FISCCP1_COSP'),lambda cld: cosp_histogram_standardize(rename(cld))),
        (('CLISCCP'),lambda cld: cosp_histogram_standardize(rename(cld)))
    ]),
    'ICEFRAC': OrderedDict([
        (('ICEFRAC'), lambda icefrac: convert_units(icefrac, target_units="%"))
    ]),
    'RELHUM': OrderedDict([
        (('hur'), lambda hur: convert_units(hur, target_units="%")),
        (('RELHUM'), lambda relhum: convert_units(relhum, target_units="%"))
#        (('RELHUM'), rename)
    ]),
    'OMEGA': OrderedDict([
        (('wap'), lambda wap: convert_units(wap, target_units="mbar/day")),
        (('OMEGA'), lambda omega: convert_units(omega, target_units="mbar/day"))
    ]),
    'SHUM': OrderedDict([
        (('Q'), lambda q: convert_units(rename(q), target_units="g/kg")),
        (('SHUM'), lambda shum: convert_units(shum, target_units="g/kg"))
    ])
}
