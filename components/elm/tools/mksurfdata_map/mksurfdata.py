#!/usr/bin/env python3
"""
Python version of mksurfdata.pl
Converts Perl script to Python for generating surface datasets.
"""

import argparse
import os
import subprocess
import sys
from datetime import datetime
from validation import validate_resolution, validate_years, get_valid_resolutions
from xml_query import get_namelist_definition

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Make surface datasets for all resolutions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  For supported resolutions:
    %(prog)s -res <res> [OPTIONS]
  
  For unsupported, user-specified resolutions:
    %(prog)s -res usrspec -usr_gname <name> -usr_gdate <date> [OPTIONS]
        """
    )
    
    # Main resolution options
    parser.add_argument('-r', '-res', '--res', dest='hgrid', default='all',
                        help='Resolution(s) to use for files (default: all)')
    parser.add_argument('-usr_gname', '--usr_gname',
                        help='User resolution name (only with -res usrspec)')
    parser.add_argument('-usr_gdate', '--usr_gdate',
                        help='User map date (only with -res usrspec)')
    parser.add_argument('-usr_mapdir', '--usr_mapdir', 
                        default='../../shared/mkmapdata',
                        help='Directory with user-supplied mapping files')
    
    # Year and scenario options
    parser.add_argument('-y', '-years', '--years', dest='years', default='1850,2000',
                        help='Simulation year(s) to run over (default: 1850,2000)')
    parser.add_argument('-c', '-rcp', '--rcp', dest='rcp', default='-999.9',
                        help='Representative concentration pathway(s) (default: -999.9=historical)')
    
    # Data options
    parser.add_argument('-l', '-dinlc', '--dinlc', dest='csmdata',
                        default=None,
                        help='Directory location for inputdata (DIN_LOC_ROOT)')
    parser.add_argument('-crop', '--crop', action='store_true',
                        help='Add in crop datasets')
    parser.add_argument('-hirespft', '--hirespft', action='store_true',
                        help='Use high-resolution pft dataset (3min vs half-degree)')
    parser.add_argument('-glc_nec', '--glc_nec', type=int, default=0,
                        help='Number of glacier elevation classes (default: 0)')
    parser.add_argument('-merge_gis', '--merge_gis', action='store_true',
                        help='Merge Greenland Ice Sheet data from CISM')
    parser.add_argument('-inlandwet', '--inlandwet', action='store_true',
                        help='Allow inland wetlands')
    parser.add_argument('-dynpft', '--dynpft',
                        help='Dynamic PFT/harvesting file to use')
    
    # Execution options
    parser.add_argument('-d', '-debug', '--debug', dest='debug', action='store_true',
                        help='Debug mode: print actions, do not run mksurfdata_map')
    parser.add_argument('-exedir', '--exedir',
                        help='Directory where mksurfdata_map program is')
    parser.add_argument('-allownofile', '--allownofile', action='store_true',
                        help='Allow script to run even if input files do not exist')
    parser.add_argument('-mv', '--mv', action='store_true',
                        help='Move files after creation to inputdata location')
    
    # User-specific options
    parser.add_argument('-usrname', '--usrname', default='',
                        help='CLM user data name to find grid file with')
    
    # Override options
    parser.add_argument('-pft_frc', '--pft_frc',
                        help='Comma-delimited list of percentages for veg types')
    parser.add_argument('-pft_idx', '--pft_idx',
                        help='Comma-delimited veg index for each fraction')
    parser.add_argument('-soil_cly', '--soil_cly', type=float,
                        help='Percentage of soil that is clay')
    parser.add_argument('-soil_snd', '--soil_snd', type=float,
                        help='Percentage of soil that is sand')
    parser.add_argument('-soil_col', '--soil_col', type=int,
                        help='Soil color (1=light to 20=dark)')
    parser.add_argument('-soil_fmx', '--soil_fmx', type=float,
                        help='Soil maximum saturated fraction (0-1)')
    
    return parser.parse_args()


def trim(s):
    """Trim whitespace from string."""
    return s.strip() if s else ""


def query_namelist(scrdir, csmdata, namelist, var, options=None, res=None, usrname=None):
    """
    Query the Perl queryDefaultNamelist.pl script.
    
    Args:
        scrdir: Script directory
        csmdata: CSMDATA directory
        namelist: Namelist name (e.g., 'default_settings', 'elmexp')
        var: Variable to query
        options: Dictionary of options (e.g., {'type': 'veg', 'hgrid': '0.9x1.25'})
        res: Resolution
        usrname: User resolution name
    
    Returns:
        String value from query
    """
    query_script = f"{scrdir}/../../bld/queryDefaultNamelist.pl"
    
    cmd = [query_script, '-csmdata', csmdata, '-silent', '-justvalue', 
           '-namelist', namelist]
    
    if res and res != 'all' and res != 'usrspec':
        cmd.extend(['-res', res])
    
    if usrname:
        cmd.extend(['-usrname', usrname])
    
    if options:
        opts_str = ','.join([f"{k}={v}" for k, v in options.items()])
        cmd.extend(['-options', opts_str])
    
    cmd.extend(['-var', var])
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return trim(result.stdout)
    except subprocess.CalledProcessError as e:
        # Only show detailed error in debug mode or for specific cases
        stderr_msg = e.stderr.strip() if e.stderr else ''
        if 'invalid resolution' in stderr_msg.lower():
            # Don't print warning here - let the caller handle it with better context
            return ""
        else:
            print(f"\nWarning: Query failed for variable '{var}'")
            if options:
                print(f"         Options: {options}")
            print(f"         Error: {stderr_msg}")
            return ""


def check_soil(opts):
    """Check that soil options are set correctly."""
    for soil_type in ['soil_cly', 'soil_snd']:
        if opts.get(soil_type) is None:
            sys.exit(f"ERROR: Soil variables were set, but {soil_type} was NOT set")
    
    texsum = opts['soil_cly'] + opts['soil_snd']
    loam = 100.0 - texsum
    if texsum < 0.0 or texsum > 100.0:
        sys.exit(f"ERROR: Soil textures are out of range: clay={opts['soil_cly']} "
                 f"sand={opts['soil_snd']} loam={loam}")


def check_soil_col_fmx(opts):
    """Check that soil color or soil fmax option is set correctly."""
    if opts.get('soil_col') is not None:
        if opts['soil_col'] < 0 or opts['soil_col'] > 20:
            sys.exit(f"ERROR: Soil color is out of range = {opts['soil_col']}")
    
    if opts.get('soil_fmx') is not None:
        if opts['soil_fmx'] < 0.0 or opts['soil_fmx'] > 1.0:
            sys.exit(f"ERROR: Soil fmax is out of range = {opts['soil_fmx']}")


def check_pft(opts, numpft):
    """Check that pft options are set correctly."""
    # Eliminate starting and ending square brackets
    pft_idx = opts['pft_idx'].strip('[]')
    pft_frc = opts['pft_frc'].strip('[]')
    
    for pft_type in ['pft_idx', 'pft_frc']:
        if opts.get(pft_type) is None:
            sys.exit(f"ERROR: PFT variables were set, but {pft_type} was NOT set")
    
    pft_idx_list = [int(x) for x in pft_idx.split(',')]
    pft_frc_list = [float(x) for x in pft_frc.split(',')]
    
    if len(pft_idx_list) != len(pft_frc_list):
        sys.exit("ERROR: PFT arrays are different sizes: pft_idx and pft_frc")
    
    sumfrc = 0.0
    for i, idx in enumerate(pft_idx_list):
        # Check index in range
        if idx < 0 or idx > numpft:
            sys.exit(f"ERROR: pft_idx out of range = {opts['pft_idx']}")
        
        # Make sure there are no duplicates
        for j in range(i):
            if idx == pft_idx_list[j]:
                sys.exit(f"ERROR: pft_idx has duplicates = {opts['pft_idx']}")
        
        # Check fraction in range
        frc = pft_frc_list[i]
        if frc <= 0.0 or frc > 100.0:
            sys.exit(f"ERROR: pft_frc out of range (>0.0 and <=100.0) = {opts['pft_frc']}")
        
        sumfrc += frc
    
    # Check that fraction sums up to 100%
    if abs(sumfrc - 100.0) > 1.0e-6:
        sys.exit(f"ERROR: pft_frc does NOT add up to 100% = {opts['pft_frc']}")


def query_data_files(res, opts, scrdir):
    """
    Query for all necessary mapping and data files.
    
    Args:
        res: Resolution string
        opts: Options dictionary
        scrdir: Script directory
    
    Returns:
        Dictionary with 'map' and 'datfil' subdictionaries
    """
    # Data types to query
    # Note: 'veg' is queried separately later with sim_year option
    data_types = ['lak', 'veg', 'voc', 'top', 'tex', 'col', 'ord',
                  'fmx', 'lai', 'urb', 'org', 'glc', 'utp', 'wet',
                  'gdp', 'peat', 'abm', 'topostats', 'vic', 'ch4',
                  'pho', 'grvl', 'slp10', 'ero', 'fert', 'toprad']
    
    # Data file types (excludes veg since mksrf_fvegtyp needs sim_year)
    data_file_types = [t for t in data_types if t != 'veg']
    
    hirespft = 'on' if opts['hirespft'] else 'off'
    merge_gis = 'on' if opts['merge_gis'] else 'off'
    csmdata = opts['csmdata']
    mkcrop = 'on' if opts['crop'] else 'off'
    usrname = opts['usrname'] if opts['usrname'] else None
    
    map_files = {}
    data_files = {}
    hgrid_map = {}
    lmask_map = {}
    
    for typ in data_types:
        # Query for lmask, hgrid, and filename variable for this type
        lmask = query_namelist(scrdir, csmdata, 'default_settings', 'lmask',
                              options={'type': typ, 'mergeGIS': merge_gis, 'hirespft': hirespft},
                              res=res, usrname=usrname)
        
        if not lmask:
            sys.exit(f"\nERROR: Could not determine land mask for resolution '{res}' and type '{typ}'.\n"
                    f"       This usually means the resolution is not supported.\n"
                    f"       Check that '{res}' is a valid resolution in the namelist definition files.")
        
        hgrid = query_namelist(scrdir, csmdata, 'default_settings', 'hgrid',
                              options={'type': typ, 'hirespft': hirespft},
                              res=res, usrname=usrname)
        
        if not hgrid:
            sys.exit(f"\nERROR: Could not determine source grid for resolution '{res}' and type '{typ}'.\n"
                    f"       This usually means the resolution is not supported.\n"
                    f"       Check that '{res}' is a valid resolution in the namelist definition files.")
        
        hgrid_map[typ] = hgrid
        lmask_map[typ] = lmask
        
        # Query for mapping file
        if opts['hgrid'] == 'usrspec':
            mapdate = opts['usr_gdate']
            map_files[typ] = f"{opts['usr_mapdir']}/map_{hgrid}_{lmask}_to_{res}_nomask_aave_da_c{mapdate}.nc"
        else:
            map_files[typ] = query_namelist(scrdir, csmdata, 'elmexp', 'map',
                                           options={'frm_hgrid': hgrid, 'frm_lmask': lmask,
                                                   'to_hgrid': res, 'to_lmask': 'nomask'},
                                           usrname=usrname)
        
        # Check mapping file
        if not map_files[typ] or not map_files[typ].strip():
            sys.exit(f"\nERROR: Could not find mapping file for resolution '{res}' (type '{typ}').\n"
                    f"       Source grid: {hgrid}, Land mask: {lmask}\n"
                    f"       Target resolution: {res}")
        
        if not opts['allownofile'] and not os.path.isfile(map_files[typ]):
            sys.exit(f"\nERROR: Mapping file does not exist:\n       {map_files[typ]}\n"
                    f"       Resolution: {res}, Type: {typ}")
        
        # Query for data file (skip veg - will be queried later with sim_year)
        if typ in data_file_types:
            filnm = query_namelist(scrdir, csmdata, 'default_settings', 'mksrf_filename',
                                  options={'type': typ},
                                  res=res, usrname=usrname)
            
            data_files[typ] = query_namelist(scrdir, csmdata, 'elmexp', filnm,
                                            options={'hgrid': hgrid, 'lmask': lmask,
                                                    'mergeGIS': merge_gis, 'crop': mkcrop},
                                            usrname=usrname)
            
            # Check data file
            if not data_files[typ] or not data_files[typ].strip():
                sys.exit(f"\nERROR: Could not find data file '{filnm}' for type '{typ}'.\n"
                        f"       Source grid: {hgrid}, Land mask: {lmask}")
            
            if not opts['allownofile'] and not os.path.isfile(data_files[typ]):
                sys.exit(f"\nERROR: Data file does not exist:\n       {data_files[typ]}\n"
                        f"       Type: {typ}")
    
    # Get grid data file from the pft map file or grid if not found
    griddata = map_files.get('veg', '')
    if not griddata:
        griddata = query_namelist(scrdir, csmdata, 'default_settings', 'fatmgrid',
                                 res=res, usrname=usrname)
        if not griddata:
            sys.exit(f"ERROR: could NOT find a grid data file for this resolution: {res}.")
    
    return {
        'map': map_files,
        'datfil': data_files,
        'hgrd': hgrid_map,
        'lmsk': lmask_map,
        'griddata': griddata
    }


def main():
    """Main execution function."""
    args = parse_args()
    
    # Convert args to dict for easier handling
    opts = vars(args)
    
    # Determine DIN_LOC_ROOT
    csmdata = opts['csmdata'] or os.environ.get('DIN_LOC_ROOT')
    if not csmdata:
        sys.exit('ERROR: DIN_LOC_ROOT must be set via --dinlc or environment variable.')
    opts['csmdata'] = csmdata
    
    # Get basic options
    hgrid = opts['hgrid']
    years = opts['years']
    rcp = opts['rcp']
    debug = opts['debug']
    
    # Set numpft based on crop option
    numpft = 50 if opts['crop'] else 16
    
    # Validate usrspec options
    if hgrid == 'usrspec':
        if not opts['usr_gname'] or not opts['usr_gdate']:
            sys.exit("ERROR: -usr_gname and -usr_gdate required when -res usrspec")
    
    # Check soil options
    soil_override = False
    if opts['soil_cly'] is not None or opts['soil_snd'] is not None:
        check_soil(opts)
        soil_override = True
    check_soil_col_fmx(opts)
    
    # Check pft options
    pft_override = False
    if opts['pft_frc'] or opts['pft_idx']:
        check_pft(opts, numpft)
        pft_override = True
    
    # Check dynpft file exists
    if opts['dynpft'] and not os.path.isfile(opts['dynpft']):
        sys.exit(f"ERROR: Dynamic PFT file does NOT exist: {opts['dynpft']}")
    
    # Generate date string
    sdate = 'c' + datetime.now().strftime('%y%m%d')
    
    # Set up basic namelist filename
    nl = 'namelist'
    
    # Get resolution list
    if hgrid == 'all':
        # Query XML for all valid resolutions
        hresols = get_valid_resolutions()
        if not hresols:
            print("WARNING: Could not query valid resolutions from XML")
            hresols = ['0.9x1.25']  # Fallback
    elif hgrid == 'usrspec':
        hresols = [opts['usr_gname']]
    else:
        hresols = hgrid.split(',')
        # Validate non-usrspec resolutions unless usrname is specified
        if not opts['usrname']:
            valid_resolutions = get_valid_resolutions()
            for res in hresols:
                if res not in valid_resolutions:
                    print(f"\nERROR: Invalid resolution: '{res}'")
                    print(f"\nValid resolutions are:")
                    for vr in sorted(valid_resolutions):
                        print(f"  {vr}")
                    sys.exit(1)
    
    # Parse years list
    years_list = years.split(',')
    
    # Parse rcp list  
    rcp_list = rcp.split(',')
    
    # Determine inland wetlands setting
    no_inlandwet = '.false.' if opts['inlandwet'] else '.true.'
    
    # Get script directory for queries
    scrdir = os.path.dirname(os.path.abspath(__file__))
    
    # Process each resolution
    for res in hresols:
        print(f"\nProcessing resolution: {res}")
        
        # Query for data files for this resolution
        print(f"  Querying data files for resolution {res}...")
        file_info = query_data_files(res, opts, scrdir)
        map_files = file_info['map']
        data_files = file_info['datfil']
        griddata = file_info['griddata']
        
        # Check for all-urban single point datasets
        all_urb_list = ['1x1_camdenNJ', '1x1_vancouverCAN', 
                        '1x1_mexicocityMEX', '1x1_urbanc_alpha']
        all_urb = '.true.' if res in all_urb_list else '.false.'
        urb_pt = 0
        if res in all_urb_list and res != '1x1_camdenNJ':
            urb_pt = 1
        
        # Always run at double precision for output
        double = '.true.'
        
        # Process each RCP
        for rcp_val in rcp_list:
            # Process each simulation year
            for sim_year in years_list:
                # Skip urban points unless sim_year=2000
                if urb_pt and sim_year != '2000':
                    print(f"  For urban -- skip simulation year = {sim_year}")
                    continue
                
                print(f"  rcp={rcp_val}, sim_year={sim_year}")
                
                # Parse year range if applicable
                if '-' in sim_year:
                    sim_yr0, sim_yrn = sim_year.split('-')
                    # Special case: 1850-2000 becomes 1850-2005
                    if sim_year == '1850-2000':
                        sim_year = '1850-2005'
                        sim_yr0, sim_yrn = '1850', '2005'
                        print(f"    Note: Converting 1850-2000 to 1850-2005")
                else:
                    sim_yr0 = sim_yrn = sim_year
                
                # Generate landuse timeseries file if needed
                landuse_ts_file = None
                if sim_year != sim_yr0:
                    if float(rcp_val) == -999.9:
                        desc = f"hist_simyr{int(sim_yr0):04d}-{int(sim_yrn):04d}"
                    else:
                        desc = f"rcp{float(rcp_val):.1f}_simyr{int(sim_yr0):04d}-{int(sim_yrn):04d}"
                    
                    # Only generate if not using override and dynpft not specified
                    if not pft_override and not opts['dynpft']:
                        landuse_ts_file = f"landuse_timeseries_{desc}.txt"
                        print(f"  Generating landuse timeseries file: {landuse_ts_file}")
                        
                        # This is a stub - real implementation would query for vegtyp 
                        # files for each year in the range
                        with open(landuse_ts_file, 'w') as luf:
                            for yr in range(int(sim_yr0), int(sim_yrn) + 1):
                                # Placeholder - would query actual vegtyp file
                                vegtyp_yr = f"{csmdata}/lnd/clm2/rawdata/vegtyp_{yr}.nc"
                                luf.write(f"{vegtyp_yr:<195} {yr:4d}\n")
                                if yr % 100 == 0:
                                    print(f"    Processing year: {yr}")
                        
                        print(f"  Done generating landuse timeseries file")
                    elif opts['dynpft']:
                        landuse_ts_file = opts['dynpft']
                        print(f"  Using user-provided dynpft file: {landuse_ts_file}")
                    elif pft_override and opts['dynpft']:
                        landuse_ts_file = opts['dynpft']
                    else:
                        # Generate override version
                        landuse_ts_file = f"landuse_timeseries_override_{desc}.txt"
                        print(f"  Generating override landuse timeseries file: {landuse_ts_file}")
                        
                        frstpft = (f"<pft_f>{opts['pft_frc']}</pft_f>"
                                   f"<pft_i>{opts['pft_idx']}</pft_i>"
                                   "<harv>0,0,0,0,0</harv><graz>0</graz>")
                        
                        with open(landuse_ts_file, 'w') as luf:
                            luf.write(f"{frstpft:<195} {int(sim_yr0):4d}\n")
                
                # Generate namelist content with all mapping and data files
                with open(nl, 'w') as f:
                    f.write('&elmexp\n')
                    f.write(f' nglcec            = {opts["glc_nec"]}\n')
                    f.write(f" mksrf_fgrid       = '{griddata}'\n")
                    f.write(f" map_fpft          = '{map_files['veg']}'\n")
                    f.write(f" map_fglacier      = '{map_files['glc']}'\n")
                    f.write(f" map_fsoicol       = '{map_files['col']}'\n")
                    f.write(f" map_fsoiord       = '{map_files['ord']}'\n")
                    f.write(f" map_furban        = '{map_files['urb']}'\n")
                    f.write(f" map_fmax          = '{map_files['fmx']}'\n")
                    f.write(f" map_forganic      = '{map_files['org']}'\n")
                    f.write(f" map_flai          = '{map_files['lai']}'\n")
                    f.write(f" map_fharvest      = '{map_files['lai']}'\n")
                    f.write(f" map_flakwat       = '{map_files['lak']}'\n")
                    f.write(f" map_fwetlnd       = '{map_files['wet']}'\n")
                    f.write(f" map_fvocef        = '{map_files['voc']}'\n")
                    f.write(f" map_fsoitex       = '{map_files['tex']}'\n")
                    f.write(f" map_furbtopo      = '{map_files['utp']}'\n")
                    f.write(f" map_flndtopo      = '{map_files['top']}'\n")
                    f.write(f" map_fgdp          = '{map_files['gdp']}'\n")
                    f.write(f" map_fpeat         = '{map_files['peat']}'\n")
                    f.write(f" map_fabm          = '{map_files['abm']}'\n")
                    f.write(f" map_ftopostats    = '{map_files['topostats']}'\n")
                    f.write(f" map_fvic          = '{map_files['vic']}'\n")
                    f.write(f" map_fch4          = '{map_files['ch4']}'\n")
                    f.write(f" map_fphosphorus   = '{map_files['pho']}'\n")
                    f.write(f" map_fgrvl         = '{map_files['grvl']}'\n")
                    f.write(f" map_fslp10        = '{map_files['slp10']}'\n")
                    f.write(f" map_fero          = '{map_files['ero']}'\n")
                    f.write(f" map_ffert         = '{map_files['fert']}'\n")
                    f.write(f" map_ftoprad       = '{map_files['toprad']}'\n")
                    f.write(f" mksrf_fsoitex     = '{data_files['tex']}'\n")
                    f.write(f" mksrf_forganic    = '{data_files['org']}'\n")
                    f.write(f" mksrf_flakwat     = '{data_files['lak']}'\n")
                    f.write(f" mksrf_fwetlnd     = '{data_files['wet']}'\n")
                    f.write(f" mksrf_fmax        = '{data_files['fmx']}'\n")
                    f.write(f" mksrf_fglacier    = '{data_files['glc']}'\n")
                    f.write(f" mksrf_fvocef      = '{data_files['voc']}'\n")
                    f.write(f" mksrf_furbtopo    = '{data_files['utp']}'\n")
                    f.write(f" mksrf_flndtopo    = '{data_files['top']}'\n")
                    f.write(f" mksrf_fgdp        = '{data_files['gdp']}'\n")
                    f.write(f" mksrf_fpeat       = '{data_files['peat']}'\n")
                    f.write(f" mksrf_fabm        = '{data_files['abm']}'\n")
                    f.write(f" mksrf_ftopostats  = '{data_files['topostats']}'\n")
                    f.write(f" mksrf_fvic        = '{data_files['vic']}'\n")
                    f.write(f" mksrf_fch4        = '{data_files['ch4']}'\n")
                    f.write(f" outnc_double      = {double}\n")
                    f.write(f' all_urban         = {all_urb}\n')
                    f.write(f' no_inlandwet      = {no_inlandwet}\n')
                    f.write(f" mksrf_furban      = '{data_files['urb']}'\n")
                    f.write(f" mksrf_fphosphorus = '{data_files['pho']}'\n")
                    f.write(f" mksrf_fgrvl       = '{data_files['grvl']}'\n")
                    f.write(f" mksrf_fslp10      = '{data_files['slp10']}'\n")
                    f.write(f" mksrf_fero        = '{data_files['ero']}'\n")
                    f.write(f" mksrf_ffert       = '{data_files['fert']}'\n")
                    f.write(f" mksrf_ftoprad     = '{data_files['toprad']}'\n")
                    
                    # Query for vegetation type file
                    resol_veg = file_info['hgrd']['veg']
                    query_opts = {
                        'sim_year': sim_yr0,
                        'crop': 'on' if opts['crop'] else 'off',
                        'rcp': rcp_val
                    }
                    
                    vegtyp = query_namelist(scrdir, csmdata, 'elmexp', 'mksrf_fvegtyp',
                                           options=query_opts,
                                           res=resol_veg, usrname=opts['usrname'])
                    
                    if not vegtyp:
                        sys.exit(f"ERROR: could not query vegtyp file for sim_year={sim_yr0}")
                    
                    # Add soil overrides if specified
                    if soil_override:
                        f.write(f' soil_clay     = {opts["soil_cly"]}\n')
                        f.write(f' soil_sand     = {opts["soil_snd"]}\n')
                    
                    # Add pft overrides if specified
                    if pft_override:
                        f.write(f' pft_frc      = {opts["pft_frc"]}\n')
                        f.write(f' pft_idx      = {opts["pft_idx"]}\n')
                    
                    # Add vegetation type and other specified files
                    f.write(f" mksrf_fvegtyp  = '{vegtyp}'\n")
                    f.write(f" mksrf_fsoicol  = '{data_files['col']}'\n")
                    f.write(f" mksrf_fsoiord  = '{data_files['ord']}'\n")
                    f.write(f" mksrf_flai     = '{data_files['lai']}'\n")
                    
                    # Generate output filenames
                    if float(rcp_val) == -999.9:
                        desc_yr0 = f"simyr{int(sim_yr0):04d}"
                        if sim_year != sim_yr0:
                            desc = f"hist_simyr{int(sim_yr0):04d}-{int(sim_yrn):04d}"
                        else:
                            desc = desc_yr0
                    else:
                        desc_yr0 = f"rcp{float(rcp_val):.1f}_simyr{int(sim_yr0):04d}"
                        if sim_year != sim_yr0:
                            desc = f"rcp{float(rcp_val):.1f}_simyr{int(sim_yr0):04d}-{int(sim_yrn):04d}"
                        else:
                            desc = desc_yr0
                    
                    crpdes = "mp50_" if opts['crop'] else ""
                    ofile = f"surfdata_{res}_{crpdes}{desc_yr0}_{sdate}"
                    f.write(f" fsurdat        = '{ofile}.nc'\n")
                    f.write(f" fsurlog        = '{ofile}.log'\n")
                    
                    # Handle transient case
                    if sim_year != sim_yr0 and landuse_ts_file:
                        f.write(f" mksrf_fdynuse  = '{landuse_ts_file}'\n")
                        ofile_ts = f"landuse.timeseries_{res}_{desc}_{sdate}"
                        f.write(f" fdyndat        = '{ofile_ts}.nc'\n")
                    else:
                        f.write(" mksrf_fdynuse  = ' '\n")
                        f.write(" fdyndat        = ' '\n")
                    
                    # Set crop numpft if applicable
                    if opts['crop']:
                        f.write(f' numpft = {numpft}\n')
                    
                    f.write('/\n')
                
                print(f"  Generated namelist: {nl}")
                
                # Print namelist for debugging
                if debug or True:  # Always print for now
                    print("  Namelist contents:")
                    with open(nl, 'r') as f:
                        for line in f:
                            print(f"    {line}", end='')
                
                # Run mksurfdata_map or debug
                if debug:
                    print(f"  Debug mode: would run mksurfdata_map < {nl}")
                    # Create dummy output files
                    open(f'{ofile}.nc', 'a').close()
                    open(f'{ofile}.log', 'a').close()
                else:
                    # Determine executable directory
                    exedir = opts['exedir'] if opts['exedir'] else '.'
                    cmd = f"{exedir}/mksurfdata_map < {nl}"
                    print(f"  Running: {cmd}")
                    ret = subprocess.call(cmd, shell=True)
                    if ret != 0:
                        sys.exit(f'ERROR: mksurfdata_map failed with exit code {ret}')
                    
                    # Check that output files were created
                    if not os.path.isfile(f'{ofile}.nc'):
                        sys.exit(f"ERROR: Output file {ofile}.nc was NOT created")
                    if not os.path.isfile(f'{ofile}.log'):
                        sys.exit(f"ERROR: Log file {ofile}.log was NOT created")
                    
                    # Handle file movement if requested
                    if opts['mv']:
                        surfdir = "lnd/clm2/surfdata"
                        outdir = f"{csmdata}/{surfdir}"
                        
                        if not os.path.isdir(outdir):
                            print(f"  WARNING: Output directory does not exist: {outdir}")
                        else:
                            # Move netcdf file
                            src_nc = f'{ofile}.nc'
                            dst_nc = f"{outdir}/{ofile}.nc"
                            print(f"  Moving {src_nc} to {dst_nc}")
                            os.rename(src_nc, dst_nc)
                            os.chmod(dst_nc, 0o444)
                            
                            # Move log file
                            src_log = f'{ofile}.log'
                            dst_log = f"{outdir}/{ofile}.log"
                            print(f"  Moving {src_log} to {dst_log}")
                            os.rename(src_log, dst_log)
                            os.chmod(dst_log, 0o444)
                            
                            # Move transient file if exists
                            if sim_year != sim_yr0 and os.path.isfile(f'{ofile_ts}.nc'):
                                src_ts = f'{ofile_ts}.nc'
                                dst_ts = f"{outdir}/{ofile_ts}.nc"
                                print(f"  Moving {src_ts} to {dst_ts}")
                                os.rename(src_ts, dst_ts)
                                os.chmod(dst_ts, 0o444)
                
                print("  Done with this case.\n")
    
    print("Successfully completed all cases.")

if __name__ == '__main__':
    main()
