#!/usr/bin/env python
"""
Resolution-based adjustment of PS for E3SM initial conditions.

This script adjusts surface pressure to match topography at the target
resolution and adds it to the vertical coordinate file for use in NCO
vertical remapping.

Usage:
    python adjust_ps_m2m.py --vertfile VERTFILE --rawicfile RAWICFILE --topofile TOPOFILE --mapped-phis MAPPED_PHISFILE [--phis-var PHIS_VAR]
    
Arguments:
    --vertfile: Path to vertical coordinate file
    --rawicfile: Path to raw initial condition file (remapped, with PS and T)
    --topofile: Path to target grid topography file
    --mapped-phis: Path to mapped PHIS file
    --phis-var: Variable name for PHIS (default: 'PHIS')
"""

import argparse
import numpy as np
import netCDF4 as nc


def ana_ini_psreset(ps, t0, phis1, phis2, tsair_in, sigma=None):
    """
    Adjust surface pressure to match topography changes.
    
    Based on IDL code ini3psreset.pro for CAM3/E3SM.
    Refer to physics/cam1/cpslec.F90
    
    Parameters
    ----------
    ps : array_like
        Original surface pressure (Pa)
    t0 : array_like
        Temperature at lowest model level (K)
    phis1 : array_like
        Directly interpolated geopotential height (m^2/s^2)
    phis2 : array_like
        Adjusted/filtered geopotential height for target resolution (m^2/s^2)
    tsair_in : array_like
        Initial surface air temperature (K) - will be adjusted
    sigma : float, optional
        Sigma level at bottom of model (hybm at lowest level).
        If not provided, assumes t0 is 0-30mb above ground.
    
    Returns
    -------
    psl : ndarray
        Adjusted sea level pressure (Pa)
    tsair : ndarray
        Adjusted surface air temperature (K)
    """
    
    # Convert to numpy arrays and flatten
    t = np.atleast_1d(np.array(t0)).flatten()
    ps = np.atleast_1d(np.array(ps)).flatten()
    phis1 = np.atleast_1d(np.array(phis1)).flatten()
    phis2 = np.atleast_1d(np.array(phis2)).flatten()
    tsair = np.atleast_1d(np.array(tsair_in)).flatten().copy()
    
    if t.ndim > 1:
        raise ValueError("Only a 1D t array needed (lowest level)")
    
    # Calculate difference between directly interpolated and filtered PHIS
    phis = phis1 - phis2
    
    # Physical constants
    gravit = 9.80616  # m/s^2
    rair = 287.04     # J/kg/K
    xlapse = 6.5e-3   # K/m
    alpha = rair * xlapse / gravit
    
    # Determine pressure at model level
    if sigma is None:
        pmid = ps - 30 * 100.0  # 30 hPa = 3000 Pa
    else:
        pmid = ps * sigma
    
    # Initialize output arrays
    ncol = len(ps)
    psl = np.zeros(ncol, dtype=np.float64)
    
    print(f"\nPHIS difference (phis1 - phis2) stats:")
    print(f"  Min: {np.min(phis):.6f} m²/s²")
    print(f"  Max: {np.max(phis):.6f} m²/s²")
    print(f"  Mean: {np.mean(phis):.6f} m²/s²")
    
    # Adjust pressure for each column
    for i in range(ncol):
        if abs(phis[i] / gravit) < 1.e-4:
            psl[i] = float(ps[i])
        else:
            Tstar = float(t[i]) * (1.0 + alpha * (float(ps[i]) / float(pmid[i]) - 1.0))
            TT0 = Tstar + xlapse * float(phis[i]) / gravit
            tsair[i] = float(tsair[i]) + xlapse * float(phis[i]) / gravit
            
            if Tstar <= 290.5 and TT0 > 290.5:
                alph = rair / float(phis[i]) * (290.5 - Tstar)
            elif Tstar > 290.5 and TT0 > 290.5:
                alph = 0.0
                Tstar = 0.5 * (290.5 + Tstar)
            else:
                alph = alpha
                if Tstar < 255.0:
                    Tstar = 0.5 * (255.0 + Tstar)
            
            beta = float(phis[i]) / (rair * Tstar)
            psl[i] = float(ps[i]) * np.exp(beta * (1.0 - alph * beta / 2.0 + 
                                             (alph * beta)**2 / 3.0))
    
    return psl, tsair


def read_netcdf_vars(filename, varnames):
    """
    Read variables from a NetCDF file.
    
    Parameters
    ----------
    filename : str
        Path to NetCDF file
    varnames : list of str
        Variable names to read. Can use 'var1:var2' syntax to rename var1 to var2
    
    Returns
    -------
    dict
        Dictionary with variable names as keys and data as values
    """
    data = {}
    with nc.Dataset(filename, 'r') as ds:
        for varname in varnames:
            if ':' in varname:
                source_name, target_name = varname.split(':')
            else:
                source_name = target_name = varname
            
            if source_name in ds.variables:
                ncvar = ds.variables[source_name]
                var_data = ncvar[:]
                
                print(f"    Reading {source_name}: shape={var_data.shape}, dims={ncvar.dimensions}")
                
                if np.ma.is_masked(var_data):
                    var_data = var_data.filled(np.nan)
                
                if source_name in ['hyam', 'hybm', 'hyai', 'hybi', 'lat', 'lon']:
                    var_data = np.squeeze(var_data)
                
                data[target_name.lower()] = var_data
            else:
                print(f"Warning: Variable {source_name} not found in {filename}")
    
    return data


def main(vertfile, rawicfile, topofile, mapped_phis_file, phis_var='PHIS'):
    """Main function to adjust PS based on topography."""
    
    print(f"Vert file: {vertfile}")
    print(f"Raw IC file: {rawicfile}")
    print(f"Topo file: {topofile}")
    print(f"Mapped PHIS file: {mapped_phis_file}")
    print(f"PHIS varname: {phis_var}")
    
    # Read topography data
    print("\nReading topography file...")
    dm = read_netcdf_vars(topofile, [f'{phis_var}:PHIS', 'lat', 'lon'])
    print(f"  PHIS shape: {dm['phis'].shape}")
    
    # Read initial condition data
    print("\nReading initial condition file...")
    ds = read_netcdf_vars(rawicfile, ['PS', 'hyam', 'hybm', 'hyai', 'hybi', 
                                       'T', 'lat', 'lon'])
    print(f"  PS shape: {ds['ps'].shape}")
    print(f"  T shape: {ds['t'].shape}")
    print(f"  hybm shape: {ds['hybm'].shape}")
    
    # Read mapped PHIS
    print("\nReading mapped PHIS file...")
    dh = read_netcdf_vars(mapped_phis_file, [f'{phis_var}:PHIS'])
    print(f"  Mapped PHIS shape: {dh['phis'].shape}")
    
    ds['phis'] = dh['phis']
    
    nlev = len(ds['hybm'])
    print(f"\nNumber of vertical levels: {nlev}")
    
    # Extract PS
    ps_data = ds['ps']
    print(f"Original PS shape: {ps_data.shape}")
    
    if ps_data.ndim == 1:
        ps_data = ps_data
    elif ps_data.ndim == 2:
        if ps_data.shape[0] == 1:
            ps_data = ps_data[0, :]
        elif ps_data.shape[1] == 1:
            ps_data = ps_data[:, 0]
        else:
            ps_data = ps_data[0, :]
    elif ps_data.ndim == 3:
        ps_data = ps_data[0, :, 0] if ps_data.shape[0] == 1 else ps_data[0, 0, :]
    
    ps_data = np.atleast_1d(ps_data).flatten()
    print(f"Processed PS shape: {ps_data.shape}")
    
    # Extract temperature
    t_data = ds['t']
    print(f"Original T shape: {t_data.shape}")
    
    if t_data.ndim == 1:
        tbot = t_data
    elif t_data.ndim == 2:
        if t_data.shape[1] == nlev:
            tbot = t_data[:, -1]
        elif t_data.shape[0] == nlev:
            tbot = t_data[-1, :]
        else:
            raise ValueError(f"Cannot determine T dimension order: shape={t_data.shape}, nlev={nlev}")
    elif t_data.ndim == 3:
        lev_dim = None
        for i, dim_size in enumerate(t_data.shape):
            if dim_size == nlev:
                lev_dim = i
                break
        
        if lev_dim is None:
            raise ValueError(f"Cannot find lev dimension in T: shape={t_data.shape}, nlev={nlev}")
        
        if lev_dim == 0:
            tbot = t_data[-1, 0, :] if t_data.shape[1] < t_data.shape[2] else t_data[-1, :, 0]
        elif lev_dim == 1:
            tbot = t_data[0, -1, :] if t_data.shape[0] < t_data.shape[2] else t_data[:, -1, 0]
        else:
            tbot = t_data[0, :, -1] if t_data.shape[0] < t_data.shape[1] else t_data[:, 0, -1]
    elif t_data.ndim == 4:
        lev_dim = None
        for i, dim_size in enumerate(t_data.shape):
            if dim_size == nlev:
                lev_dim = i
                break
        if lev_dim == 1:
            tbot = t_data[0, -1, :, :].flatten()
        else:
            raise ValueError(f"Unexpected 4D T array structure: shape={t_data.shape}")
    else:
        raise ValueError(f"Unexpected T array dimensions: {t_data.shape}")
    
    tbot = np.atleast_1d(tbot).flatten()
    print(f"Processed T shape: {tbot.shape}")
    
    phis_src = np.atleast_1d(ds['phis']).flatten()
    phis_tgt = np.atleast_1d(dm['phis']).flatten()
    
    print(f"\nFinal data shapes:")
    print(f"  PS: {ps_data.shape}")
    print(f"  T: {tbot.shape}")
    print(f"  PHIS_src: {phis_src.shape}")
    print(f"  PHIS_tgt: {phis_tgt.shape}")
    
    if not (len(ps_data) == len(tbot) == len(phis_src) == len(phis_tgt)):
        raise ValueError(f"Dimension mismatch: PS={len(ps_data)}, T={len(tbot)}, "
                        f"PHIS_src={len(phis_src)}, PHIS_tgt={len(phis_tgt)}")
    
    # Physical constants
    gravit = 9.80616
    rair = 287.04
    xlapse = 6.5e-3
    alpha = rair * xlapse / gravit
    
    pmid_bot = ps_data * ds['hybm'][-1]
    tsair = tbot * (1.0 + alpha * (ps_data / pmid_bot - 1.0))
    
    print(f"\nPhysical constants:")
    print(f"  gravit: {gravit}")
    print(f"  rair: {rair}")
    print(f"  xlapse: {xlapse}")
    print(f"  alpha: {alpha}")
    print(f"  Bottom sigma (hybm[-1]): {ds['hybm'][-1]}")
    
    print(f"\nBottom level temperature (tbot) stats:")
    print(f"  Min: {np.min(tbot):.6f} K")
    print(f"  Max: {np.max(tbot):.6f} K")
    print(f"  Mean: {np.mean(tbot):.6f} K")
    
    print(f"\nOriginal PS stats:")
    print(f"  Min: {np.min(ps_data):.6f} Pa")
    print(f"  Max: {np.max(ps_data):.6f} Pa")
    print(f"  Mean: {np.mean(ps_data):.6f} Pa")
    
    print(f"\nBottom level pressure (pmid_bot) stats:")
    print(f"  Min: {np.min(pmid_bot):.6f} Pa")
    print(f"  Max: {np.max(pmid_bot):.6f} Pa")
    print(f"  Mean: {np.mean(pmid_bot):.6f} Pa")
    
    print(f"\nDerived surface air temperature (tsair) stats:")
    print(f"  Min: {np.min(tsair):.6f} K")
    print(f"  Max: {np.max(tsair):.6f} K")
    print(f"  Mean: {np.mean(tsair):.6f} K")
    
    print(f"\nPHIS_src (mapped) stats:")
    print(f"  Min: {np.min(phis_src):.6f} m²/s²")
    print(f"  Max: {np.max(phis_src):.6f} m²/s²")
    print(f"  Mean: {np.mean(phis_src):.6f} m²/s²")
    
    print(f"\nPHIS_tgt (topo) stats:")
    print(f"  Min: {np.min(phis_tgt):.6f} m²/s²")
    print(f"  Max: {np.max(phis_tgt):.6f} m²/s²")
    print(f"  Mean: {np.mean(phis_tgt):.6f} m²/s²")
    
    # Adjust surface pressure
    psnew, tsair_adjusted = ana_ini_psreset(
        ps_data, tbot, phis_src, phis_tgt, tsair,
        sigma=float(ds['hybm'][-1])
    )
    
    print(f"\nAdjusted PS stats:")
    print(f"  Min: {np.min(psnew):.6f} Pa")
    print(f"  Max: {np.max(psnew):.6f} Pa")
    print(f"  Mean: {np.mean(psnew):.6f} Pa")
    
    print(f"\nPS change stats:")
    print(f"  Min: {np.min(psnew - ps_data):.6f} Pa")
    print(f"  Max: {np.max(psnew - ps_data):.6f} Pa")
    print(f"  Mean: {np.mean(psnew - ps_data):.6f} Pa")
    print(f"  Std: {np.std(psnew - ps_data):.6f} Pa")
    
    ncol = len(phis_tgt)
    print(f"\nnlev and ncol: {nlev}, {ncol}")
    
    # Write adjusted PS to vertical coordinate file
    with nc.Dataset(vertfile, 'a') as ncfile:
        print(f"\nExisting dimensions in {vertfile}:")
        for dim_name, dim in ncfile.dimensions.items():
            print(f"  {dim_name}: {len(dim)} {'(unlimited)' if dim.isunlimited() else ''}")
        
        if 'ncol' not in ncfile.dimensions:
            print(f"Creating ncol dimension with size {ncol}")
            ncfile.createDimension('ncol', ncol)
        else:
            print(f"ncol dimension already exists with size {len(ncfile.dimensions['ncol'])}")
        
        if 'time' not in ncfile.dimensions:
            print("Creating time dimension (unlimited)")
            ncfile.createDimension('time', None)
        else:
            print(f"time dimension already exists")
        
        if 'PS' in ncfile.variables:
            existing_ps = ncfile.variables['PS']
            print(f"Existing PS variable dimensions: {existing_ps.dimensions}")
            print(f"Existing PS variable shape: {existing_ps.shape}")
            ps_var = existing_ps
        else:
            dim_order = ('time', 'ncol')
            print(f"Creating new PS variable with dimensions: {dim_order}")
            ps_var = ncfile.createVariable('PS', 'f8', dim_order)
            ps_var.long_name = 'Surface pressure'
            ps_var.units = 'Pa'
        
        if ps_var.dimensions == ('ncol', 'time'):
            print(f"Writing PS with shape ({len(psnew)}, 1) to dimensions (ncol, time)")
            ps_var[:, 0] = psnew
        elif ps_var.dimensions == ('time', 'ncol'):
            print(f"Writing PS with shape (1, {len(psnew)}) to dimensions (time, ncol)")
            ps_var[0, :] = psnew
        else:
            raise ValueError(f"Unexpected PS dimension order: {ps_var.dimensions}")
    
    print(f"Successfully wrote adjusted PS to {vertfile}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Adjust surface pressure for E3SM initial conditions based on topography'
    )
    parser.add_argument('--vertfile', required=True, 
                       help='Path to vertical coordinate file')
    parser.add_argument('--rawicfile', required=True,
                       help='Path to raw initial condition file (remapped, with PS and T)')
    parser.add_argument('--topofile', required=True,
                       help='Path to target grid topography file')
    parser.add_argument('--mapped-phis', dest='mapped_phis', required=True,
                       help='Path to mapped PHIS file')
    parser.add_argument('--phis-var', dest='phis_var', default='PHIS',
                       help='Variable name for PHIS (default: PHIS)')
    
    args = parser.parse_args()
    
    main(args.vertfile, args.rawicfile, args.topofile, args.mapped_phis, args.phis_var)
