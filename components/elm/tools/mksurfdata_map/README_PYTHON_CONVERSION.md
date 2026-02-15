# Python Conversion of mksurfdata.pl

## Overview

This directory now contains a Python version of the `mksurfdata.pl` script called `mksurfdata.py`. The Python script provides the same functionality as the original Perl script for creating namelist files and running the mksurfdata_map tool.

## Usage

The Python script supports the same command-line arguments as the Perl version:

```bash
./mksurfdata.py -r 0.5x0.5 -y 1850 -d -l /global/cfs/cdirs/e3sm/inputdata
```

This is equivalent to the Perl command:

```bash
./mksurfdata.pl -res 0.5x0.5 -years 1850 -d -dinlc /global/cfs/cdirs/e3sm/inputdata
```

## Command-Line Arguments

### Basic Options
- `-r, --res`: Resolution(s) to use (default: 'all')
- `-y, --years`: Simulation year(s) or year range (default: '1850,2000')
- `-l, --dinlc`: Directory location for inputdata (default: '/compyfs/inputdata')
- `-d, --debug`: Debug mode - print commands without executing
- `-c, --rcp`: Representative concentration pathway (default: '-999.9' for historical)

### Feature Flags
- `--crop`: Add crop datasets
- `--hirespft`: Use high-resolution PFT dataset
- `--merge_gis`: Merge Greenland Ice Sheet data from CISM
- `--inlandwet`: Allow inland wetlands
- `--mv`: Move files to correct location in inputdata after creation
- `--allownofile`: Allow script to run even if input files don't exist

### Advanced Options
- `--exedir`: Directory where mksurfdata_map program is located
- `--glc_nec`: Number of glacier elevation classes (default: 0)
- `--dynpft`: Dynamic PFT/harvesting file to use
- `--usrname`: CLM user data name to find grid file

### User-Specified Resolutions
- `-usr_gname`: User resolution name (for -res usrspec)
- `-usr_gdate`: User map date (for -res usrspec)
- `--usr_mapdir`: Directory for user-supplied mapping files

### Override Options
- `--pft_frc`: Comma-delimited list of PFT fractions
- `--pft_idx`: Comma-delimited list of PFT indices
- `--soil_cly`: Percentage of soil that is clay
- `--soil_snd`: Percentage of soil that is sand
- `--soil_col`: Soil color (1=light to 20=dark)
- `--soil_fmx`: Soil maximum saturated fraction (0-1)

## Examples

### Basic Usage
```bash
# Generate surface data for 0.5x0.5 resolution for year 1850
./mksurfdata.py -r 0.5x0.5 -y 1850 -l /global/cfs/cdirs/e3sm/inputdata

# Debug mode (don't actually run, just show what would happen)
./mksurfdata.py -r 0.5x0.5 -y 1850 -d -l /global/cfs/cdirs/e3sm/inputdata

# Generate for multiple years
./mksurfdata.py -r 0.5x0.5 -y 1850,2000 -l /global/cfs/cdirs/e3sm/inputdata

# Generate for year range (transient)
./mksurfdata.py -r 0.5x0.5 -y 1850-2000 -l /global/cfs/cdirs/e3sm/inputdata
```

### With Additional Features
```bash
# Include crop datasets
./mksurfdata.py -r 0.5x0.5 -y 2000 --crop -l /global/cfs/cdirs/e3sm/inputdata

# Use high-resolution PFT data
./mksurfdata.py -r 0.5x0.5 -y 2000 --hirespft -l /global/cfs/cdirs/e3sm/inputdata

# Allow inland wetlands
./mksurfdata.py -r 0.5x0.5 -y 1850 --inlandwet -l /global/cfs/cdirs/e3sm/inputdata

# Move files to inputdata directory after creation
./mksurfdata.py -r 0.5x0.5 -y 1850 --mv -l /global/cfs/cdirs/e3sm/inputdata
```

### Multiple Resolutions and Years
```bash
# Multiple resolutions
./mksurfdata.py -r 0.5x0.5,1x1 -y 1850 -l /global/cfs/cdirs/e3sm/inputdata

# Multiple years and RCP scenarios
./mksurfdata.py -r 0.5x0.5 -y 2000,2050,2100 -c 4.5 -l /global/cfs/cdirs/e3sm/inputdata
```

## Key Differences from Perl Version

1. **Command-line arguments**: 
   - Python uses `-r` instead of `-res`
   - Python uses `-y` instead of `-years`
   - Python uses `-l` instead of `-dinlc`
   - Long form arguments now use `--` prefix (e.g., `--debug` instead of `-debug`)

2. **Dependencies**:
   - Requires Python 3.6+
   - No external Python packages needed (uses only standard library)
   - Still requires the same Perl scripts for querying XML defaults (`queryDefaultNamelist.pl`)

3. **XML Parsing**:
   - The script still relies on the Perl `queryDefaultNamelist.pl` script for XML queries
   - The 'all' resolution option requires XML parsing and is not fully implemented
   - Specify explicit resolutions instead of using 'all'

4. **Error Handling**:
   - Python version has more explicit error messages
   - Better validation of input parameters

## Implementation Notes

- The Python script maintains the same logic flow as the Perl version
- All validation functions (`check_soil`, `check_pft`, etc.) are preserved
- Namelist generation follows the same format
- The script still calls the Fortran `mksurfdata_map` executable
- Query operations still use `queryDefaultNamelist.pl` via subprocess calls

## Requirements

- Python 3.6 or higher
- Perl (for queryDefaultNamelist.pl)
- Access to E3SM inputdata directory
- Compiled `mksurfdata_map` executable

## Known Limitations

1. **Resolution 'all' option**: Not fully implemented due to XML parsing requirements. Use explicit resolution lists instead.

2. **XML Dependencies**: The script still depends on Perl scripts for XML queries. A future enhancement would be to implement direct XML parsing in Python.

3. **Validation**: Some validation that would be performed by querying XML definitions is simplified or skipped.

## Troubleshooting

### Script fails with "queryDefaultNamelist.pl not found"
- Ensure the Perl script exists at `../../bld/queryDefaultNamelist.pl`
- Check that the script is executable

### "Resolution not supported" error
- Do not use 'all' - specify explicit resolutions
- Check that the resolution exists in the XML configuration

### Files not created
- Check that mksurfdata_map executable exists and is compiled
- Verify inputdata paths are correct
- Run with `-d` (debug) flag to see commands without executing

## Migration Guide

To migrate from Perl to Python:

```bash
# Old Perl command:
./mksurfdata.pl -res 0.5x0.5 -years 1850 -d -dinlc /path/to/inputdata

# New Python command (equivalent):
./mksurfdata.py -r 0.5x0.5 -y 1850 -d -l /path/to/inputdata
```

Both scripts can coexist in the same directory, so you can test the Python version while keeping the Perl version as a backup.

## Contributing

If you find bugs or have improvements, please report them to the E3SM development team.
