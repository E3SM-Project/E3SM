"""
Interface to the env_archive.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.archive_base import ArchiveBase
from CIME.XML.env_base import EnvBase
from CIME.hist_utils import get_extension

logger = logging.getLogger(__name__)
# pylint: disable=super-init-not-called
class EnvArchive(ArchiveBase,EnvBase):

    def __init__(self, case_root=None, infile="env_archive.xml", read_only=False):
        """
        initialize an object interface to file env_archive.xml in the case directory
        """
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_archive.xsd")
        EnvBase.__init__(self, case_root, infile, schema=schema, read_only=read_only)

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
        hist_files = []
        extensions = self.get_hist_file_extensions(self.get_entry(model))
        regex = self.get_hist_file_regex(self.get_entry(model))
        for ext in extensions:
            if 'initial' in ext:
                continue
            if extension.endswith('$'):
                extension = extension[:-1]
            string = model+r'\d?_?(\d{4})?\.'+extension+suffix+'$'
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
        rest_files = []
        extensions = self.get_rest_file_extensions(self.get_entry(model))
        regex = self.get_rest_file_regex(self.get_entry(model))
        for ext in extensions:
            if 'initial' in ext:
                continue
            if extension.endswith('$'):
                extension = extension[:-1]
            string = model+r'\d?_?(\d{4})?\.'+extension+suffix+'$'
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

    def get_entries(self):
        return self.get_children('comp_archive_spec')

    def get_entry_info(self, archive_entry):
        compname = self.get(archive_entry, 'compname')
        compclass = self.get(archive_entry, 'compclass')
        return compname,compclass

    def get_rpointer_contents(self, archive_entry):
        rpointer_items = []
        rpointer_nodes = self.get_children('rpointer', root=archive_entry)
        for rpointer_node in rpointer_nodes:
            file_node = self.get_child('rpointer_file', root=rpointer_node)
            content_node = self.get_child('rpointer_content', root=rpointer_node)
            rpointer_items.append([self.text(file_node),self.text(content_node)])
        return rpointer_items

    def get_type_info(self, vid):
        return "char"
