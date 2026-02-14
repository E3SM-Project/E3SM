
import xml.etree.ElementTree as ET
import os
import re

NAMELIST_DEF_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    '../../bld/namelist_files/namelist_definition.xml'))

def parse_xml_with_regex(xml_path, entry_id):
    """Parse XML using regex as fallback when ET fails."""
    try:
        with open(xml_path, 'r') as f:
            content = f.read()
        # Find entry with matching id
        pattern = rf'<entry\s+id="{entry_id}"[^>]*valid_values\s*=\s*["\']([^"\']*)["\']'
        match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)
        if match:
            valid_values = match.group(1)
            return [v.strip() for v in valid_values.split(',') if v.strip()]
    except Exception as e:
        print(f"Warning: Could not parse XML with regex: {e}")
    return []

def get_valid_resolutions(xml_path=NAMELIST_DEF_PATH):
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        for entry in root.findall('entry'):
            if entry.attrib.get('id') == 'res':
                valid = entry.attrib.get('valid_values', '')
                return [v.strip() for v in valid.split(',') if v.strip()]
    except ET.ParseError:
        # Silently fall back to regex parsing (XML file has malformed tags)
        return parse_xml_with_regex(xml_path, 'res')
    return []

def get_valid_years(xml_path=NAMELIST_DEF_PATH):
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        for entry in root.findall('entry'):
            if entry.attrib.get('id') == 'sim_year':
                valid = entry.attrib.get('valid_values', '')
                return [v.strip() for v in valid.split(',') if v.strip()]
    except ET.ParseError:
        # Silently fall back to regex parsing (XML file has malformed tags)
        return parse_xml_with_regex(xml_path, 'sim_year')
    return []

def get_valid_year_ranges(xml_path=NAMELIST_DEF_PATH):
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        for entry in root.findall('entry'):
            if entry.attrib.get('id') == 'sim_year_range':
                valid = entry.attrib.get('valid_values', '')
                return [v.strip() for v in valid.split(',') if v.strip()]
    except ET.ParseError:
        # Silently fall back to regex parsing (XML file has malformed tags)
        return parse_xml_with_regex(xml_path, 'sim_year_range')
    return []

def validate_resolution(res):
    valid_res = get_valid_resolutions()
    if res not in valid_res:
        print(f"ERROR: Invalid resolution: {res}")
        print(f"Valid options: {', '.join(valid_res)}")
        exit(1)

def validate_years(years):
    valid_years = get_valid_years()
    valid_ranges = get_valid_year_ranges()
    for y in years.split(","):
        if "-" in y:
            if y not in valid_ranges:
                print(f"ERROR: Invalid year range: {y}")
                print(f"Valid ranges: {', '.join(valid_ranges)}")
                exit(1)
        else:
            if y not in valid_years:
                print(f"ERROR: Invalid year: {y}")
                print(f"Valid years: {', '.join(valid_years)}")
                exit(1)
