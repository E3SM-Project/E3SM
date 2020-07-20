"""
Base class for archive files.  This class inherits from generic_xml.py
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML

logger = logging.getLogger(__name__)

class ArchiveBase(GenericXML):

    def get_entry(self, compname):
        """
        Returns an xml node corresponding to compname in comp_archive_spec
        """
        return self.scan_optional_child('comp_archive_spec',
                                        attributes={"compname":compname})

    def _get_file_node_text(self, attnames, archive_entry):
        """
        get the xml text associated with each of the attnames
        based at root archive_entry
        returns a list of text entries or
        an empty list if no entries are found
        """
        nodes = []
        textvals = []
        for attname in attnames:
            nodes.extend(self.get_children(attname, root=archive_entry))
        for node in nodes:
            textvals.append(self.text(node))
        return textvals

    def get_rest_file_extensions(self, archive_entry):
        """
        get the xml text associated with each of the rest_file_extensions
        based at root archive_entry (root is based on component name)
        returns a list of text entries or
        an empty list if no entries are found
        """
        return self._get_file_node_text(['rest_file_extension'],archive_entry)

    def get_hist_file_extensions(self, archive_entry):
        """
        get the xml text associated with each of the hist_file_extensions
        based at root archive_entry (root is based on component name)
        returns a list of text entries or
        an empty list if no entries are found
        """
        return self._get_file_node_text(['hist_file_extension'],archive_entry)

    def get_entry_value(self, name, archive_entry):
        """
        get the xml text associated with name under root archive_entry
        returns None if no entry is found, expects only one entry
        """
        node = self.get_optional_child(name, root=archive_entry)
        if node is not None:
            return self.text(node)
        return None

    def get_latest_hist_files(self, casename, model, from_dir, suffix="", ref_case=None):
        """
        get the most recent history files in directory from_dir with suffix if provided
        """
        test_hists = self.get_all_hist_files(casename, model, from_dir, suffix=suffix, ref_case=ref_case)
        latest_files = {}
        histlist = []
        for hist in test_hists:
            ext = _get_extension(model, hist)
            latest_files[ext] = hist

        for key in latest_files.keys():
            histlist.append(latest_files[key])
        return histlist

    def get_all_hist_files(self, casename, model, from_dir, suffix="", ref_case=None):
        """
        gets all history files in directory from_dir with suffix (if provided)
        ignores files with ref_case in the name if ref_case is provided
        """
        dmodel = model
        if model == "cpl":
            dmodel = "drv"
        # remove when component name is changed
        if model == "fv3gfs":
            model = "fv3"
        hist_files = []
        extensions = self.get_hist_file_extensions(self.get_entry(dmodel))
        if suffix and len(suffix) > 0:
            has_suffix = True
        else:
            has_suffix = False

        # Strip any trailing $ if suffix is present and add it back after the suffix
        for ext in extensions:
            if ext.endswith('$') and has_suffix:
                ext = ext[:-1]
            string = model+r'\d?_?(\d{4})?\.'+ext
            if has_suffix:
                string += '.'+suffix+'$'


            logger.debug ("Regex is {}".format(string))
            pfile = re.compile(string)
            hist_files.extend([f for f in os.listdir(from_dir) if pfile.search(f) and ( (f.startswith(casename) or f.startswith(model)) and not f.endswith("cprnc.out") )])

        if ref_case:
            expect(ref_case not in casename,"ERROR: ref_case name {} conflicts with casename {}".format(ref_case,casename))
            hist_files = [h for h in hist_files if not (ref_case in os.path.basename(h))]

        hist_files = list(set(hist_files))
        hist_files.sort()
        logger.debug("get_all_hist_files returns {} for model {}".format(hist_files, model))

        return hist_files

def _get_extension(model, filepath):
    r"""
    For a hist file for the given model, return what we call the "extension"

    model - The component model
    filepath - The path of the hist file
    >>> _get_extension("cpl", "cpl.hi.nc")
    'hi'
    >>> _get_extension("cpl", "cpl.h.nc")
    'h'
    >>> _get_extension("cpl", "cpl.h1.nc.base")
    'h1'
    >>> _get_extension("cpl", "TESTRUNDIFF.cpl.hi.0.nc.base")
    'hi'
    >>> _get_extension("cpl", "TESTRUNDIFF_Mmpi-serial.f19_g16_rx1.A.melvin_gnu.C.fake_testing_only_20160816_164150-20160816_164240.cpl.h.nc")
    'h'
    >>> _get_extension("clm","clm2_0002.h0.1850-01-06-00000.nc")
    '0002.h0'
    >>> _get_extension("pop","PFS.f09_g16.B1850.cheyenne_intel.allactive-default.GC.c2_0_b1f2_int.pop.h.ecosys.nday1.0001-01-02.nc")
    'h'
    >>> _get_extension("mom", "ga0xnw.mom6.frc._0001_001.nc")
    'frc'
    >>> _get_extension("mom", "ga0xnw.mom6.sfc.day._0001_001.nc")
    'sfc.day'
    >>> _get_extension("mom", "bixmc5.mom6.prog._0001_01_05_84600.nc")
    'prog'
    >>> _get_extension("mom", "bixmc5.mom6.hm._0001_01_03_42300.nc")
    'hm'
    >>> _get_extension("mom", "bixmc5.mom6.hmz._0001_01_03_42300.nc")
    'hmz'
    >>> _get_extension("pop", "casename.pop.dd.0001-01-02-00000")
    'dd'
    """
    # Remove with component namechange
    if model == "fv3gfs":
        model = "fv3"
    basename = os.path.basename(filepath)
    m = None
    ext_regexes = []

    # First add any model-specific extension regexes; these will be checked before the
    # general regex
    if model == "mom":
        # Need to check 'sfc.day' specially: the embedded '.' messes up the
        # general-purpose regex
        ext_regexes.append(r'sfc\.day')

    # Now add the general-purpose extension regex
    ext_regexes.append(r'\w+')

    for ext_regex in ext_regexes:
        full_regex_str = model+r'\d?_?(\d{4})?\.('+ext_regex+r')[-\w\.]*'
        full_regex = re.compile(full_regex_str)
        m = full_regex.search(basename)
        if m is not None:
            if m.group(1) is not None:
                result = m.group(1)+'.'+m.group(2)
            else:
                result = m.group(2)
            return result

    expect(m, "Failed to get extension for file '{}'".format(filepath))

    return result
