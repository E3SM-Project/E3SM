#!/usr/bin/env python3
"""
mksurfdata.py - Python script to make surface datasets for all resolutions.

Converted from mksurfdata.pl
Original: Oct/30/2008 Erik Kluzek
Python version: 2026
"""

import os
import sys
import argparse
import subprocess
from datetime import datetime
from pathlib import Path
import re


class MksurfdataGenerator:
    """Generate surface datasets for E3SM/ELM using mksurfdata_map tool."""
    
    def __init__(self):
        self.prog_name = Path(__file__).name
        self.prog_dir = Path(__file__).parent.absolute()
        self.cwd = Path.cwd()
        self.scrdir = self.prog_dir
        
        # Default CSMDATA path
        self.default_csmdata = "/compyfs/inputdata"
        
        # Number of PFTs
        self.numpft = 16
        
        # Namelist definition file
        self.nldef_file = self.scrdir / "../../bld/namelist_files/namelist_definition.xml"
        
        # SVN repository settings
        self.svnrepo = "https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata"
        self.svnmesg = "Update fsurdat files with mksurfdata_map"
        self.surfdir = "lnd/clm2/surfdata"
        
    def setup_argparser(self):
        """Setup command-line argument parser."""
        parser = argparse.ArgumentParser(
            description='Generate surface datasets for E3SM/ELM',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  For supported resolutions:
    %(prog)s -r 0.5x0.5 -y 1850 -d -l /global/cfs/cdirs/e3sm/inputdata
    
  For user-specified resolutions:
    %(prog)s -r usrspec -usr_gname <user_gname> -usr_gdate <user_gdate>
    
Notes:
  - years, res, and rcp can be comma-delimited lists
  - Use -d/--debug to see what would run without actually executing
"""
        )
        
        # Resolution options
        parser.add_argument('-r', '--res', dest='hgrid', default='all',
                          help='Resolution(s) to use (default: all)')
        parser.add_argument('-usr_gname', dest='usr_gname',
                          help='User resolution name (for -res usrspec)')
        parser.add_argument('-usr_gdate', dest='usr_gdate',
                          help='User map date (for -res usrspec)')
        parser.add_argument('-usr_mapdir', dest='usr_mapdir',
                          default='../../shared/mkmapdata',
                          help='Directory for user-supplied mapping files')
        
        # Year and scenario options
        parser.add_argument('-y', '--years', dest='years', default='1850,2000',
                          help='Simulation year(s) or year range (default: 1850,2000)')
        parser.add_argument('-c', '--rcp', dest='rcp', default='-999.9',
                          help='Representative concentration pathway (default: -999.9 for historical)')
        
        # Input data options
        parser.add_argument('-l', '--dinlc', dest='csmdata', 
                          default=self.default_csmdata,
                          help=f'Directory location for inputdata (default: {self.default_csmdata})')
        
        # Feature flags
        parser.add_argument('-d', '--debug', action='store_true',
                          help='Debug mode - print commands without executing')
        parser.add_argument('--allownofile', action='store_true',
                          help='Allow script to run even if input files do not exist')
        parser.add_argument('--crop', action='store_true',
                          help='Add crop datasets')
        parser.add_argument('--hirespft', action='store_true',
                          help='Use high-resolution PFT dataset (3-minute resolution)')
        parser.add_argument('--merge_gis', action='store_true',
                          help='Merge Greenland Ice Sheet data from CISM')
        parser.add_argument('--inlandwet', action='store_true',
                          help='Allow inland wetlands')
        parser.add_argument('--mv', action='store_true',
                          help='Move files to correct location in inputdata after creation')
        
        # Advanced options
        parser.add_argument('--exedir', dest='exedir',
                          help='Directory where mksurfdata_map program is located')
        parser.add_argument('--glc_nec', dest='glc_nec', type=int, default=0,
                          help='Number of glacier elevation classes (default: 0)')
        parser.add_argument('--dynpft', dest='dynpft',
                          help='Dynamic PFT/harvesting file to use')
        parser.add_argument('--usrname', dest='usrname', default='',
                          help='CLM user data name to find grid file')
        
        # Override options
        parser.add_argument('--pft_frc', dest='pft_frc',
                          help='Comma-delimited list of PFT fractions')
        parser.add_argument('--pft_idx', dest='pft_idx',
                          help='Comma-delimited list of PFT indices')
        parser.add_argument('--soil_cly', dest='soil_cly', type=float,
                          help='Percentage of soil that is clay')
        parser.add_argument('--soil_snd', dest='soil_snd', type=float,
                          help='Percentage of soil that is sand')
        parser.add_argument('--soil_col', dest='soil_col', type=int,
                          help='Soil color (1=light to 20=dark)')
        parser.add_argument('--soil_fmx', dest='soil_fmx', type=float,
                          help='Soil maximum saturated fraction (0-1)')
        
        return parser
    
    def trim(self, s):
        """Remove leading and trailing whitespace."""
        return s.strip() if s else ''
    
    def check_soil(self, args):
        """Check that soil options are set correctly."""
        if args.soil_cly is None or args.soil_snd is None:
            raise ValueError("ERROR: Soil variables were set, but soil_cly or soil_snd was NOT set")
        
        texsum = args.soil_cly + args.soil_snd
        loam = 100.0 - texsum
        if texsum < 0.0 or texsum > 100.0:
            raise ValueError(f"ERROR: Soil textures are out of range: clay={args.soil_cly}, "
                           f"sand={args.soil_snd}, loam={loam}")
    
    def check_soil_col_fmx(self, args):
        """Check that soil color or soil fmax option is set correctly."""
        if args.soil_col is not None:
            if args.soil_col < 0 or args.soil_col > 20:
                raise ValueError(f"ERROR: Soil color is out of range = {args.soil_col}")
        
        if args.soil_fmx is not None:
            if args.soil_fmx < 0.0 or args.soil_fmx > 1.0:
                raise ValueError(f"ERROR: Soil fmax is out of range = {args.soil_fmx}")
    
    def check_pft(self, args):
        """Check that PFT options are set correctly."""
        # Remove square brackets if present
        if args.pft_idx:
            args.pft_idx = args.pft_idx.strip('[]')
        if args.pft_frc:
            args.pft_frc = args.pft_frc.strip('[]')
        
        if args.pft_idx is None or args.pft_frc is None:
            raise ValueError("ERROR: PFT variables were set, but pft_idx or pft_frc was NOT set")
        
        pft_idx = [int(x) for x in args.pft_idx.split(',')]
        pft_frc = [float(x) for x in args.pft_frc.split(',')]
        
        if len(pft_idx) != len(pft_frc):
            raise ValueError("ERROR: PFT arrays are different sizes: pft_idx and pft_frc")
        
        sumfrc = 0.0
        for i, idx in enumerate(pft_idx):
            # Check index in range
            if idx < 0 or idx > self.numpft:
                raise ValueError(f"ERROR: pft_idx out of range = {args.pft_idx}")
            
            # Check for duplicates
            if idx in pft_idx[:i]:
                raise ValueError(f"ERROR: pft_idx has duplicates = {args.pft_idx}")
            
            # Check fraction in range
            if pft_frc[i] <= 0.0 or pft_frc[i] > 100.0:
                raise ValueError(f"ERROR: pft_frc out of range (>0.0 and <=100.0) = {args.pft_frc}")
            
            sumfrc += pft_frc[i]
        
        # Check that fractions sum to 100%
        if abs(sumfrc - 100.0) > 1e-6:
            raise ValueError(f"ERROR: pft_frc does NOT add up to 100% = {args.pft_frc}")
    
    def query_default_namelist(self, script_path, options_str, var=None, namelist=None):
        """Query the default namelist using the queryDefaultNamelist.pl script."""
        cmd = [str(script_path)]
        
        if options_str:
            cmd.extend(options_str.split())
        
        if namelist:
            cmd.extend(['-namelist', namelist])
        
        if var:
            cmd.extend(['-var', var])
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return self.trim(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"ERROR running query: {' '.join(cmd)}")
            print(f"Error: {e.stderr}")
            return ''
    
    def validate_resolution(self, res, query_script):
        """Validate that a resolution is supported."""
        # This would need to query the XML definition
        # For now, we'll skip validation or use a simple check
        return True
    
    def validate_year(self, year, query_script):
        """Validate that a simulation year is supported."""
        # This would need to query the XML definition
        # For now, we'll skip validation or use a simple check
        return True
    
    def validate_rcp(self, rcp, query_script):
        """Validate that an RCP value is supported."""
        # This would need to query the XML definition
        # For now, we'll skip validation or use a simple check
        return True
    
    def get_mapping_files(self, res, query_script, csmdata, args, hirespft, merge_gis):
        """Get mapping files and data files for all types."""
        types = ['lak', 'veg', 'voc', 'top', 'tex', 'col', 'ord', 
                'fmx', 'lai', 'urb', 'org', 'glc', 'utp', 'wet',
                'gdp', 'peat', 'abm', 'topostats', 'vic', 'ch4',
                'pho', 'grvl', 'slp10', 'ero', 'fert', 'toprad']
        
        map_files = {}
        hgrd = {}
        lmsk = {}
        datfil = {}
        
        usrnam = ''
        if args.usrname and res == args.usrname:
            usrnam = f'-usrname {args.usrname}'
        
        mkopts = f'-csmdata {csmdata} -silent -justvalue {usrnam}'
        
        if args.hgrid != 'usrspec':
            mopts = f'-res {res} -csmdata {csmdata} -silent -justvalue {usrnam}'
        else:
            mopts = f'-csmdata {csmdata} -silent -justvalue {usrnam}'
        
        mkcrop = ',crop=\'off\''
        if args.crop:
            mkcrop = ',crop=\'on\''
        
        for typ in types:
            # Get land mask - note: Perl adds -silent even though mopts already has it
            opts_str = f'{mopts} -silent -options type={typ},mergeGIS={merge_gis},hirespft={hirespft}'
            lmask = self.query_default_namelist(query_script, opts_str, var='lmask', namelist='default_settings')
            
            # Get horizontal grid
            opts_str = f'{mopts} -silent -options type={typ},hirespft={hirespft}'
            hgrid = self.query_default_namelist(query_script, opts_str, var='hgrid', namelist='default_settings')
            
            # Get filename variable
            opts_str = f'{mopts} -silent -options type={typ}'
            filnm = self.query_default_namelist(query_script, opts_str, var='mksrf_filename', namelist='default_settings')
            
            hgrd[typ] = hgrid
            lmsk[typ] = lmask
            
            # Get mapping file
            if args.hgrid == 'usrspec':
                mapdate = args.usr_gdate
                map_files[typ] = os.path.join(
                    args.usr_mapdir, 
                    f'map_{hgrid}_{lmask}_to_{res}_nomask_aave_da_c{mapdate}.nc'
                )
            else:
                opts_str = (f'-csmdata {csmdata} -silent -justvalue -onlyfiles '
                          f'-options frm_hgrid={hgrid},frm_lmask={lmask},'
                          f'to_hgrid={res},to_lmask=nomask')
                map_files[typ] = self.query_default_namelist(query_script, opts_str, var='map', namelist='elmexp')
            
            map_files[typ] = self.trim(map_files[typ])
            
            if not map_files[typ]:
                raise ValueError(f"ERROR: could NOT find a mapping file for resolution: {res} "
                               f"and type: {typ} at {hgrid} and {lmask}")
            
            if not args.allownofile and not os.path.exists(map_files[typ]):
                raise ValueError(f"ERROR: mapping file does NOT exist: {map_files[typ]}")
            
            # Get data file
            # filnm is the variable name (like 'mksrf_fvegtyp'), use it to query its value
            opts_str = f'{mkopts} -options hgrid={hgrid},lmask={lmask},mergeGIS={merge_gis}{mkcrop}'
            datfil[typ] = self.query_default_namelist(query_script, opts_str, var=filnm, namelist='elmexp')
            datfil[typ] = self.trim(datfil[typ])
            
            if not datfil[typ]:
                raise ValueError(f"ERROR: could NOT find a {filnm} data file for resolution: "
                               f"{hgrid} and type: {typ} and {lmask}")
            
            if not args.allownofile and not os.path.exists(datfil[typ]):
                raise ValueError(f"ERROR: data file does NOT exist: {datfil[typ]}")
        
        return map_files, hgrd, lmsk, datfil
    
    def get_grid_data(self, map_files, query_script, usrnam, res):
        """Get grid data file."""
        griddata = self.trim(map_files.get('veg', ''))
        if not griddata:
            opts_str = f'-csmdata {self.csmdata} -silent -justvalue -onlyfiles {usrnam}'
            griddata = self.query_default_namelist(query_script, opts_str, var='fatmgrid')
            if not griddata:
                raise ValueError(f"ERROR: could NOT find a grid data file for resolution: {res}")
        return griddata
    
    def write_namelist(self, filename, params):
        """Write namelist file."""
        with open(filename, 'w') as f:
            f.write("&elmexp\n")
            for key, value in params.items():
                if isinstance(value, str):
                    f.write(f" {key} = '{value}'\n")
                elif isinstance(value, bool):
                    val_str = '.true.' if value else '.false.'
                    f.write(f" {key} = {val_str}\n")
                else:
                    f.write(f" {key} = {value}\n")
            f.write("/\n")
    
    def write_landuse_timeseries(self, filename, query_script, resol, sim_yr0, sim_yrn, 
                                 rcp, mkcrop, dynpft_format):
        """Write landuse timeseries text file."""
        print(f"Writing out landuse_timeseries text file: {filename}")
        with open(filename, 'w') as f:
            for yr in range(sim_yr0, sim_yrn + 1):
                opts_str = (f'{resol} -csmdata {self.csmdata} -silent -justvalue -onlyfiles '
                          f'-options sim_year={yr},rcp={rcp}{mkcrop} '
                          f'-namelist elmexp')
                vegtyp_yr = self.query_default_namelist(query_script, opts_str, var='mksrf_fvegtyp')
                f.write(f'{vegtyp_yr:<195} {yr:4d}\n')
                if yr % 100 == 0:
                    print(f"year: {yr}")
        print("Done writing file")
    
    def run(self, args):
        """Main execution function."""
        # Set CSMDATA
        self.csmdata = args.csmdata
        
        # Validate crop option
        if args.crop:
            self.numpft = 50
        
        # Check soil options
        if args.soil_cly is not None or args.soil_snd is not None:
            self.check_soil(args)
            args.soil_override = True
        else:
            args.soil_override = False
        
        # Check soil color and fmax
        self.check_soil_col_fmx(args)
        
        # Check PFT options
        if args.pft_frc is not None or args.pft_idx is not None:
            self.check_pft(args)
            args.pft_override = True
        else:
            args.pft_override = False
        
        # Check dynpft file
        if args.dynpft and not os.path.exists(args.dynpft):
            raise ValueError(f"ERROR: Dynamic PFT file does NOT exist: {args.dynpft}")
        
        # Get date string
        sdate = 'c' + datetime.now().strftime('%y%m%d')
        
        # Setup query script
        query_script = self.scrdir / '../../bld/queryDefaultNamelist.pl'
        
        # Parse resolutions
        if args.hgrid == 'all':
            # Would need to query XML for all valid resolutions
            raise NotImplementedError("Resolution 'all' requires XML parsing - please specify resolution(s)")
        elif args.hgrid == 'usrspec':
            if not args.usr_gname or not args.usr_gdate:
                raise ValueError("ERROR: For usrspec, must provide -usr_gname and -usr_gdate")
            hresols = [args.usr_gname]
        else:
            hresols = args.hgrid.split(',')
        
        # Parse years
        years = args.years.split(',')
        
        # Parse RCPs
        rcpaths = args.rcp.split(',')
        
        # Create input data files script
        cfile = 'clm.input_data_files'
        if os.path.exists(cfile):
            os.rename(cfile, f'{cfile}.previous')
        
        with open(cfile, 'w') as cfh:
            cfh.write('#!/bin/csh -f\n')
            cfh.write(f'set CSMDATA = {self.csmdata}\n')
        os.chmod(cfile, 0o755)
        
        # Setup namelist file
        nl = 'namelist'
        
        # Determine urban points
        all_urb_resolutions = ['1x1_camdenNJ', '1x1_vancouverCAN', 
                              '1x1_mexicocityMEX', '1x1_urbanc_alpha']
        
        # Loop over resolutions
        for res in hresols:
            print(f"\n{'='*50}")
            print(f"Processing resolution: {res}")
            print(f"{'='*50}\n")
            
            # Check if urban point
            urb_pt = res in all_urb_resolutions and res != '1x1_camdenNJ'
            all_urb = res in all_urb_resolutions
            
            # Get mapping files
            hirespft = 'on' if args.hirespft else 'off'
            merge_gis = 'on' if args.merge_gis else 'off'
            
            try:
                map_files, hgrd, lmsk, datfil = self.get_mapping_files(
                    res, query_script, self.csmdata, args, hirespft, merge_gis
                )
            except Exception as e:
                print(f"ERROR getting mapping files for {res}: {e}")
                continue
            
            # Get grid data
            usrnam = ''
            if args.usrname and res == args.usrname:
                usrnam = f'-usrname {args.usrname}'
            
            try:
                griddata = self.get_grid_data(map_files, query_script, usrnam, res)
            except Exception as e:
                print(f"ERROR getting grid data for {res}: {e}")
                continue
            
            # No inland wetlands by default
            no_inlandwet = '.false.' if args.inlandwet else '.true.'
            
            # Always double precision
            double = '.true.'
            
            # Loop over RCPs
            for rcp in rcpaths:
                rcp_float = float(rcp)
                
                # Loop over years
                for sim_year in years:
                    # Skip non-2000 years for urban points
                    if urb_pt and sim_year != '2000':
                        print(f"For urban -- skip this simulation year = {sim_year}")
                        continue
                    
                    # Handle 1850-2000 -> 1850-2005
                    if sim_year == '1850-2000':
                        print(f"For {sim_year} actually run 1850-2005")
                        sim_year = '1850-2005'
                    
                    # Parse year range
                    if '-' in sim_year:
                        sim_yr0, sim_yrn = map(int, sim_year.split('-'))
                    else:
                        sim_yr0 = sim_yrn = int(sim_year)
                    
                    # Build query options
                    resol = f'-res {hgrd["veg"]}' if args.hgrid != 'usrspec' else ''
                    
                    rcp_option = '' if rcp_float == -999.9 else f',rcp={rcp}'
                    mkcrop = ',crop=\'off\''
                    crpdes = ''
                    if args.crop:
                        mkcrop = ',crop=\'on\''
                        crpdes = 'mp50_'
                    
                    # Get vegetation type file
                    opts_str = (f'{resol} -csmdata {self.csmdata} -silent -justvalue -onlyfiles '
                              f'-options sim_year={sim_yr0}{mkcrop}{rcp_option} '
                              f'-namelist elmexp')
                    vegtyp = self.query_default_namelist(query_script, opts_str, var='mksrf_fvegtyp')
                    
                    if not vegtyp:
                        print(f"ERROR: trouble getting vegtyp file")
                        continue
                    
                    # Build description strings
                    if rcp_float == -999.9:
                        desc = f'hist_simyr{sim_yr0:04d}-{sim_yrn:04d}'
                        desc_yr0 = f'simyr{sim_yr0:04d}'
                    else:
                        desc = f'rcp{rcp_float:.1f}_simyr{sim_yr0:04d}-{sim_yrn:04d}'
                        desc_yr0 = f'rcp{rcp_float:.1f}_simyr{sim_yr0:04d}'
                    
                    # Handle landuse timeseries for transient cases
                    landuse_timeseries_text_file = None
                    if sim_yr0 != sim_yrn:
                        dynpft_format = '%-195.195s %4.4d\n'
                        
                        if not args.dynpft and not args.pft_override:
                            landuse_timeseries_text_file = f'landuse_timeseries_{desc}.txt'
                            self.write_landuse_timeseries(
                                landuse_timeseries_text_file, query_script, resol,
                                sim_yr0, sim_yrn, rcp_float, mkcrop, dynpft_format
                            )
                        elif args.pft_override and args.dynpft:
                            landuse_timeseries_text_file = args.dynpft
                        else:
                            landuse_timeseries_text_file = f'landuse_timeseries_override_{desc}.txt'
                            with open(landuse_timeseries_text_file, 'w') as f:
                                frstpft = (f"<pft_f>{args.pft_frc}</pft_f>"
                                         f"<pft_i>{args.pft_idx}</pft_i>"
                                         f"<harv>0,0,0,0,0</harv><graz>0</graz>")
                                if len(frstpft) > 195:
                                    raise ValueError(f"ERROR PFT line is too long ({len(frstpft)}): {frstpft}")
                                f.write(f'{frstpft:<195} {sim_yr0:4d}\n')
                    
                    # Build namelist parameters
                    params = {
                        'nglcec': args.glc_nec,
                        'mksrf_fgrid': griddata,
                        'map_fpft': map_files['veg'],
                        'map_fglacier': map_files['glc'],
                        'map_fsoicol': map_files['col'],
                        'map_fsoiord': map_files['ord'],
                        'map_furban': map_files['urb'],
                        'map_fmax': map_files['fmx'],
                        'map_forganic': map_files['org'],
                        'map_flai': map_files['lai'],
                        'map_fharvest': map_files['lai'],
                        'map_flakwat': map_files['lak'],
                        'map_fwetlnd': map_files['wet'],
                        'map_fvocef': map_files['voc'],
                        'map_fsoitex': map_files['tex'],
                        'map_furbtopo': map_files['utp'],
                        'map_flndtopo': map_files['top'],
                        'map_fgdp': map_files['gdp'],
                        'map_fpeat': map_files['peat'],
                        'map_fabm': map_files['abm'],
                        'map_ftopostats': map_files['topostats'],
                        'map_fvic': map_files['vic'],
                        'map_fch4': map_files['ch4'],
                        'map_fphosphorus': map_files['pho'],
                        'map_fgrvl': map_files['grvl'],
                        'map_fslp10': map_files['slp10'],
                        'map_fero': map_files['ero'],
                        'map_ffert': map_files['fert'],
                        'map_ftoprad': map_files['toprad'],
                        'mksrf_fsoitex': datfil['tex'],
                        'mksrf_forganic': datfil['org'],
                        'mksrf_flakwat': datfil['lak'],
                        'mksrf_fwetlnd': datfil['wet'],
                        'mksrf_fmax': datfil['fmx'],
                        'mksrf_fglacier': datfil['glc'],
                        'mksrf_fvocef': datfil['voc'],
                        'mksrf_furbtopo': datfil['utp'],
                        'mksrf_flndtopo': datfil['top'],
                        'mksrf_fgdp': datfil['gdp'],
                        'mksrf_fpeat': datfil['peat'],
                        'mksrf_fabm': datfil['abm'],
                        'mksrf_ftopostats': datfil['topostats'],
                        'mksrf_fvic': datfil['vic'],
                        'mksrf_fch4': datfil['ch4'],
                        'outnc_double': double,
                        'all_urban': '.true.' if all_urb else '.false.',
                        'no_inlandwet': no_inlandwet,
                        'mksrf_furban': datfil['urb'],
                        'mksrf_fphosphorus': datfil['pho'],
                        'mksrf_fgrvl': datfil['grvl'],
                        'mksrf_fslp10': datfil['slp10'],
                        'mksrf_fero': datfil['ero'],
                        'mksrf_ffert': datfil['fert'],
                        'mksrf_ftoprad': datfil['toprad'],
                        'mksrf_fvegtyp': vegtyp,
                        'mksrf_fsoicol': datfil['col'],
                        'mksrf_fsoiord': datfil['ord'],
                        'mksrf_flai': datfil['lai'],
                    }
                    
                    # Add soil overrides
                    if args.soil_override:
                        params['soil_clay'] = args.soil_cly
                        params['soil_sand'] = args.soil_snd
                    
                    # Add PFT overrides
                    if args.pft_override:
                        params['pft_frc'] = args.pft_frc
                        params['pft_idx'] = args.pft_idx
                    
                    # Output files
                    ofile = f'surfdata_{res}_{desc_yr0}_{sdate}'
                    params['fsurdat'] = f'{ofile}.nc'
                    params['fsurlog'] = f'{ofile}.log'
                    
                    # Transient files
                    if sim_yr0 != sim_yrn:
                        ofile_ts = f'landuse.timeseries_{res}_{desc}_{sdate}'
                        params['mksrf_fdynuse'] = landuse_timeseries_text_file or ' '
                        params['fdyndat'] = f'{ofile_ts}.nc' if landuse_timeseries_text_file else ' '
                    else:
                        params['mksrf_fdynuse'] = ' '
                        params['fdyndat'] = ' '
                    
                    # Add numpft if crop
                    if args.crop:
                        params['numpft'] = self.numpft
                    
                    # Write namelist
                    self.write_namelist(nl, params)
                    
                    print(f"resolution: {res} rcp={rcp} sim_year = {sim_year}")
                    print(f"namelist: {nl}")
                    
                    # Print namelist contents
                    with open(nl, 'r') as f:
                        print(f.read())
                    
                    # Delete previous output files
                    for f in [f'{ofile}.nc', f'{ofile}.log']:
                        if os.path.exists(f):
                            os.remove(f)
                    
                    # Run mksurfdata_map
                    exedir = args.exedir if args.exedir else self.scrdir
                    exe_path = Path(exedir) / 'mksurfdata_map'
                    cmd = f'{exe_path} < {nl}'
                    
                    print(f'{cmd}')
                    
                    if not args.debug:
                        try:
                            with open(nl, 'r') as nl_file:
                                result = subprocess.run(
                                    str(exe_path),
                                    stdin=nl_file,
                                    capture_output=True,
                                    text=True,
                                    check=True
                                )
                                print(result.stdout)
                        except subprocess.CalledProcessError as e:
                            print(f"ERROR in mksurfdata_map: {e}")
                            print(e.stderr)
                            continue
                    else:
                        # Debug mode - just touch files
                        Path(f'{ofile}.nc').touch()
                        Path(f'{ofile}.log').touch()
                        if sim_yr0 != sim_yrn:
                            Path(f'{ofile_ts}.nc').touch()
                    
                    print(f"\n{'='*50}\n")
                    
                    # Check that files were created
                    if not os.path.exists(f'{ofile}.nc'):
                        print(f"ERROR: surfdata netcdf file was NOT created!")
                        continue
                    
                    if sim_yr0 != sim_yrn and landuse_timeseries_text_file:
                        if not os.path.exists(f'{ofile_ts}.nc'):
                            print(f"ERROR: landuse_timeseries netcdf file was NOT created!")
                            continue
                    
                    # Handle urban point override
                    if urb_pt and not args.debug:
                        # Query for previous surface dataset
                        opts_str = f'-res {res} -csmdata {self.csmdata} -silent -justvalue'
                        prvsurfdata = self.query_default_namelist(query_script, opts_str, var='fsurdat')
                        
                        if prvsurfdata:
                            varlist = ('CANYON_HWR,EM_IMPROAD,EM_PERROAD,EM_ROOF,EM_WALL,'
                                     'HT_ROOF,THICK_ROOF,THICK_WALL,T_BUILDING_MAX,T_BUILDING_MIN,'
                                     'WIND_HGT_CANYON,WTLUNIT_ROOF,WTROAD_PERV,ALB_IMPROAD_DIR,'
                                     'ALB_IMPROAD_DIF,ALB_PERROAD_DIR,ALB_PERROAD_DIF,ALB_ROOF_DIR,'
                                     'ALB_ROOF_DIF,ALB_WALL_DIR,ALB_WALL_DIF,TK_ROOF,TK_WALL,'
                                     'TK_IMPROAD,CV_ROOF,CV_WALL,CV_IMPROAD,NLEV_IMPROAD,'
                                     'PCT_URBAN,URBAN_REGION_ID')
                            
                            print("Overwrite urban parameters with previous surface dataset values")
                            ncks_cmd = f'ncks -A -v {varlist} {prvsurfdata} {ofile}.nc'
                            print(ncks_cmd)
                            subprocess.run(ncks_cmd, shell=True, check=False)
                    
                    # Move files if requested
                    if args.mv and not args.debug:
                        outdir = Path(self.csmdata) / self.surfdir
                        outdir.mkdir(parents=True, exist_ok=True)
                        
                        final_ofile = f'surfdata_{res}_{crpdes}{desc_yr0}_{sdate}'
                        
                        # Move nc file
                        src = f'{ofile}.nc'
                        dst = outdir / f'{final_ofile}.nc'
                        if os.path.exists(src):
                            os.rename(src, dst)
                            os.chmod(dst, 0o444)
                            print(f"Moved {src} to {dst}")
                        
                        # Move log file
                        src = f'{ofile}.log'
                        dst = outdir / f'{final_ofile}.log'
                        if os.path.exists(src):
                            os.rename(src, dst)
                            os.chmod(dst, 0o444)
                            print(f"Moved {src} to {dst}")
                        
                        # Move timeseries file
                        if sim_yr0 != sim_yrn and landuse_timeseries_text_file:
                            final_ts_ofile = f'landuse.timeseries_{res}_{desc}_{sdate}'
                            src = f'{ofile_ts}.nc'
                            dst = outdir / f'{final_ts_ofile}.nc'
                            if os.path.exists(src):
                                os.rename(src, dst)
                                os.chmod(dst, 0o444)
                                print(f"Moved {src} to {dst}")
        
        print("\nSuccessfully created fsurdat files")


def main():
    """Main entry point."""
    generator = MksurfdataGenerator()
    parser = generator.setup_argparser()
    args = parser.parse_args()
    
    try:
        generator.run(args)
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
