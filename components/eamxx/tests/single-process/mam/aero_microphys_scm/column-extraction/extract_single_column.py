#!/usr/bin/env python3
"""
Extract a single column from an EAMxx NetCDF output file.

This script extracts all data for a specified column index from an EAMxx
output file and writes it to a new NetCDF file suitable for use as initial
conditions in single-column tests.

Usage:
    python extract_single_column.py input.nc output.nc --col-idx 0
    python extract_single_column.py input.nc output.nc --col-idx 50 --time-idx 0
"""

import argparse
import sys
import netCDF4 as nc
import numpy as np
from pathlib import Path


def extract_column(input_file, output_file, col_idx, time_idx=None, verbose=False):
    """
    Extract a single column from an EAMxx NetCDF file.
    
    Args:
        input_file: Path to input NetCDF file
        output_file: Path to output NetCDF file
        col_idx: Column index to extract (0-based)
        time_idx: Time index to extract (0-based). If None, extracts last time step.
        verbose: Print detailed information
    """
    
    if verbose:
        print(f"Reading from: {input_file}")
        print(f"Writing to: {output_file}")
        print(f"Column index: {col_idx}")
    
    # Open input file
    with nc.Dataset(input_file, 'r') as src:
        
        # Get dimensions
        if 'ncol' not in src.dimensions:
            print("ERROR: 'ncol' dimension not found in input file")
            print(f"Available dimensions: {list(src.dimensions.keys())}")
            return False
        
        ncol = len(src.dimensions['ncol'])
        
        if col_idx >= ncol:
            print(f"ERROR: Column index {col_idx} out of range [0, {ncol-1}]")
            return False
        
        # Determine time index
        if 'time' in src.dimensions:
            ntime = len(src.dimensions['time'])
            if time_idx is None:
                time_idx = ntime - 1  # Last time step
            elif time_idx >= ntime:
                print(f"ERROR: Time index {time_idx} out of range [0, {ntime-1}]")
                return False
        else:
            time_idx = 0
            ntime = 1
        
        if verbose:
            print(f"Input has {ncol} columns and {ntime} time steps")
            print(f"Extracting time index: {time_idx}")
        
        # Create output file
        with nc.Dataset(output_file, 'w', format='NETCDF4') as dst:
            
            # Copy global attributes
            dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
            
            # Add provenance
            dst.setncattr('column_extraction_source', str(input_file))
            dst.setncattr('column_extraction_col_idx', col_idx)
            dst.setncattr('column_extraction_time_idx', time_idx)
            
            # Create dimensions with ncol=1
            for dim_name, dimension in src.dimensions.items():
                if dim_name == 'ncol':
                    dst.createDimension(dim_name, 1)
                elif dim_name == 'time':
                    dst.createDimension(dim_name, None if dimension.isunlimited() else 1)
                else:
                    dst.createDimension(
                        dim_name, 
                        len(dimension) if not dimension.isunlimited() else None
                    )
            
            if verbose:
                print(f"\nProcessing {len(src.variables)} variables:")
            
            # Copy variables
            for var_name, variable in src.variables.items():
                
                if verbose:
                    print(f"  {var_name:30s} {variable.dimensions} {variable.shape}")
                
                # Create variable in output
                dst_var = dst.createVariable(
                    var_name,
                    variable.datatype,
                    variable.dimensions,
                    fill_value=variable._FillValue if hasattr(variable, '_FillValue') else None,
                    zlib=True,
                    complevel=1
                )
                
                # Copy variable attributes
                dst_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs() 
                                   if k != '_FillValue'})
                
                # Copy data, slicing along ncol and/or time dimensions
                dims = variable.dimensions
                
                if 'ncol' not in dims and 'time' not in dims:
                    # No slicing needed
                    dst_var[:] = variable[:]
                    
                elif 'ncol' in dims and 'time' not in dims:
                    # Slice along ncol only
                    col_axis = dims.index('ncol')
                    slices = [slice(None)] * len(dims)
                    slices[col_axis] = col_idx
                    data = variable[tuple(slices)]
                    
                    # Reshape to add dimension back
                    new_shape = list(data.shape)
                    new_shape.insert(col_axis, 1)
                    data = data.reshape(new_shape)
                    dst_var[:] = data
                    
                elif 'time' in dims and 'ncol' not in dims:
                    # Slice along time only
                    time_axis = dims.index('time')
                    slices = [slice(None)] * len(dims)
                    slices[time_axis] = time_idx
                    data = variable[tuple(slices)]
                    
                    # Reshape to add dimension back
                    new_shape = list(data.shape)
                    new_shape.insert(time_axis, 1)
                    data = data.reshape(new_shape)
                    dst_var[:] = data
                    
                else:
                    # Slice along both ncol and time
                    col_axis = dims.index('ncol')
                    time_axis = dims.index('time')
                    slices = [slice(None)] * len(dims)
                    slices[col_axis] = col_idx
                    slices[time_axis] = time_idx
                    data = variable[tuple(slices)]
                    
                    # Reshape to add dimensions back
                    new_shape = []
                    old_idx = 0
                    for i, dim in enumerate(dims):
                        if dim == 'ncol':
                            new_shape.append(1)
                        elif dim == 'time':
                            new_shape.append(1)
                        else:
                            new_shape.append(data.shape[old_idx] if old_idx < len(data.shape) else 1)
                            old_idx += 1
                    
                    # Handle scalar case
                    if len(data.shape) == 0:
                        data = np.array([[[data]]])
                    else:
                        data = data.reshape(new_shape)
                    
                    dst_var[:] = data
    
    if verbose:
        print(f"\nSuccessfully extracted column {col_idx} to {output_file}")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Extract a single column from an EAMxx NetCDF output file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract column 0 from last time step
  %(prog)s input.nc output.nc --col-idx 0
  
  # Extract column 50 from time step 0
  %(prog)s input.nc output.nc --col-idx 50 --time-idx 0
  
  # Extract with verbose output
  %(prog)s input.nc output.nc --col-idx 10 -v
        """
    )
    
    parser.add_argument('input', type=Path,
                        help='Input NetCDF file')
    parser.add_argument('output', type=Path,
                        help='Output NetCDF file')
    parser.add_argument('--col-idx', type=int, required=True,
                        help='Column index to extract (0-based)')
    parser.add_argument('--time-idx', type=int, default=None,
                        help='Time index to extract (0-based). Default: last time step')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print detailed information')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}")
        return 1
    
    # Check output file doesn't exist or confirm overwrite
    if args.output.exists():
        response = input(f"Output file {args.output} exists. Overwrite? [y/N] ")
        if response.lower() not in ['y', 'yes']:
            print("Aborted.")
            return 1
    
    # Extract column
    success = extract_column(
        args.input,
        args.output,
        args.col_idx,
        args.time_idx,
        args.verbose
    )
    
    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())
