"""
Base class for archive files.  This class inherits from generic_xml.py
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.hist_utils import get_extension

logger = logging.getLogger(__name__)

class ArchiveBase(GenericXML):

    def get_entry(self, compname):
        return self.scan_optional_child('comp_archive_spec',
                                        attributes={"compname":compname})

    def get_file_node_text(self, attnames, archive_entry):
        nodes = []
        textvals = []
        for attname in attnames:
            nodes.extend(self.get_children(attname, root=archive_entry))
        for node in nodes:
            textvals.append(self.text(node))
        return textvals

    def get_rest_file_extensions(self, archive_entry):
        return self.get_file_node_text(['rest_file_extension'],archive_entry)

    def get_rest_file_regex(self, archive_entry):
        return self.get_file_node_text(['rest_file_regex'],archive_entry)

    def get_hist_file_extensions(self, archive_entry):
        return self.get_file_node_text(['hist_file_extension'],archive_entry)

    def get_hist_file_regex(self, archive_entry):
        return self.get_file_node_text(['hist_file_regex'],archive_entry)

    def get_entry_value(self, name, archive_entry):
        node = self.get_optional_child(name, root=archive_entry)
        if node is not None:
            return self.text(node)
        return None

    def get_latest_hist_files(self, model, from_dir, suffix="", ref_case=None):

        test_hists = self.get_all_hist_files(model, from_dir, suffix=suffix, ref_case=ref_case)
        latest_files = {}
        histlist = []
        regex = self.get_hist_file_regex(self.get_entry(model))
        for hist in test_hists:
            ext = get_extension(model, hist, regex=regex)
            latest_files[ext] = hist

        for key in latest_files.keys():
            histlist.append(os.path.join(from_dir,latest_files[key]))
            # Special case for fv3gfs which outputs in cubed sphere tiles
            if "tile[1-6].nc" in key:
                for i in range(1,5):
                    new_file = latest_files[key].replace("tile6.nc","tile{}.nc".format(i))
                    histlist.append(os.path.join(from_dir, new_file))

        return histlist

    def get_all_hist_files(self, model, from_dir, suffix="", ref_case=None):
        dmodel = model
        if model == "cpl":
            dmodel = "drv"
        hist_files = []
        extensions = self.get_hist_file_extensions(self.get_entry(dmodel))
        regex = self.get_hist_file_regex(self.get_entry(dmodel))

        for ext in extensions:
            if 'initial' in ext:
                continue
            if ext.endswith('$') and len(suffix)>0:
                ext = ext[:-1]
            string = model+r'\d?_?(\d{4})?\.'+ext
            if suffix and len(suffix)>0:
                string += suffix+'$'

            logger.debug ("Regex is {}".format(string))

            pfile = re.compile(string)
            hist_files.extend([f for f in os.listdir(from_dir) if pfile.search(f)])

        for match in regex:
            pfile = re.compile(match)
            hist_files.extend([f for f in os.listdir(from_dir) if pfile.search(f)])

        if ref_case:
            hist_files = [h for h in hist_files if not (ref_case in os.path.basename(h))]

        hist_files = list(set(hist_files))
        hist_files.sort()
        logger.debug("get_all_hist_files returns {} for model {}".format(hist_files, model))
        return hist_files

    def get_all_rest_files(self, model, from_dir, suffix="", ref_case=None):
        if model == "cpl":
            model = "drv"
        rest_files = []
        extensions = self.get_rest_file_extensions(self.get_entry(model))
        regex = self.get_rest_file_regex(self.get_entry(model))
        for ext in extensions:
            if 'initial' in ext:
                continue
            if ext.endswith('$'):
                ext = ext[:-1]
            string = model+r'\d?_?(\d{4})?\.'+ext+suffix+'$'
            logger.debug ("Regex is {}".format(string))
            pfile = re.compile(string)
            rest_files.extend([f for f in os.listdir(from_dir) if pfile.search(f)])
        for match in regex:
            pfile = re.compile(match)
            rest_files.extend([f for f in os.listdir(from_dir) if pfile.search(f)])

        if ref_case:
            rest_files = [h for h in rest_files if not (ref_case in os.path.basename(h))]

        rest_files = list(set(rest_files))
        rest_files.sort()
        logger.debug("get_all_rest_files returns {} for model {}".format(rest_files, model))
        return rest_files
