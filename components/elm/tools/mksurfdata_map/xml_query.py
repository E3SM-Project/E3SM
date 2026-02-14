"""
XML query utilities for mksurfdata.py
Python equivalent of queryDefaultNamelist.pl functionality
"""

import xml.etree.ElementTree as ET
import os
import sys


class NamelistDefinition:
    """Handle queries against namelist definition XML files."""
    
    def __init__(self, xml_path):
        """Initialize with path to namelist definition XML."""
        self.xml_path = xml_path
        if not os.path.isfile(xml_path):
            sys.exit(f"ERROR: Cannot find namelist definition file: {xml_path}")
        
        try:
            self.tree = ET.parse(xml_path)
            self.root = self.tree.getroot()
        except ET.ParseError as e:
            sys.exit(f"ERROR: Cannot parse XML file {xml_path}: {e}")
    
    def get_valid_values(self, entry_id):
        """Get valid values for a given entry ID."""
        for entry in self.root.findall('.//entry'):
            if entry.attrib.get('id') == entry_id:
                valid = entry.attrib.get('valid_values', '')
                if valid:
                    return [v.strip().strip("'") for v in valid.split(',') if v.strip()]
        return []
    
    def is_valid_value(self, entry_id, value):
        """Check if a value is valid for a given entry."""
        valid_values = self.get_valid_values(entry_id)
        # Remove quotes from value for comparison
        value_clean = value.strip().strip("'\"")
        return value_clean in valid_values


class QueryDefaultNamelist:
    """Query default values from XML database files."""
    
    def __init__(self, csmdata, scrdir):
        """Initialize query system."""
        self.csmdata = csmdata
        self.scrdir = scrdir
        self.bld_dir = os.path.join(scrdir, '../../bld')
    
    def query(self, namelist, var, res=None, options=None):
        """
        Query for a variable value with given options.
        
        Args:
            namelist: Namelist group (e.g., 'elmexp', 'default_settings')
            var: Variable name to query
            res: Resolution (optional)
            options: Dictionary of options (e.g., {'type': 'veg', 'hgrid': '0.9x1.25'})
        
        Returns:
            String value or empty string if not found
        """
        # Map common variable names to XML entries
        var_map = {
            'lmask': 'lmask',
            'hgrid': 'hgrid',
            'mksrf_filename': 'mksrf_filename',
            'map': 'map',
            'fatmgrid': 'fatmgrid',
            'mksrf_fvegtyp': 'mksrf_fvegtyp',
            'fsurdat': 'fsurdat',
        }
        
        # This is a simplified stub - full implementation would need to
        # parse the actual namelist XML files and match against options
        
        # For now, return placeholder values to allow testing
        if var == 'lmask':
            return 'nomask'
        elif var == 'hgrid':
            if options and 'type' in options:
                # Return appropriate grid based on type
                return '0.9x1.25'
            return res if res else '0.9x1.25'
        elif var == 'mksrf_filename':
            return f"mksrf_f{options.get('type', '')}"
        
        return ''
    
    def get_mapping_file(self, frm_hgrid, frm_lmask, to_hgrid, to_lmask='nomask'):
        """Get path to mapping file."""
        # In real implementation, this would query the XML database
        # For now, construct expected filename
        mapdir = f"{self.csmdata}/lnd/clm2/mappingdata/maps"
        mapfile = f"{mapdir}/map_{frm_hgrid}_{frm_lmask}_to_{to_hgrid}_{to_lmask}_aave_da.nc"
        return mapfile
    
    def get_data_file(self, filetype, hgrid, lmask, **kwargs):
        """Get path to data file."""
        # Stub implementation
        # Real version would query XML database for actual file paths
        datadir = f"{self.csmdata}/lnd/clm2/rawdata"
        return f"{datadir}/{filetype}_{hgrid}_{lmask}.nc"


def get_scrdir():
    """Get script directory."""
    return os.path.dirname(os.path.abspath(__file__))


def get_namelist_definition():
    """Get the namelist definition object."""
    scrdir = get_scrdir()
    nldef_file = os.path.join(scrdir, '../../bld/namelist_files/namelist_definition.xml')
    return NamelistDefinition(nldef_file)
