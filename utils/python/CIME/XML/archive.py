"""
Interface to the archive.xml file.  This class inherits from GenericXML.py

"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.utils import expect, get_cime_root, get_model

logger = logging.getLogger(__name__)

class Archive(GenericXML):

    def __init__(self, infile=None, files=None):
        """
        initialize an object
        """
        if files is None:
            files = Files()
        schema = files.get_schema("ARCHIVE_SPEC_FILE")

        GenericXML.__init__(self, infile, schema)

    def setup(self, env_archive, components, files=None):
        if files is None:
            files = Files()

        if 'cpl' not in components:
            components.append('cpl')
        for comp in components:
            infile = files.get_value("ARCHIVE_SPEC_FILE", {"component":comp})

            if infile is not None and os.path.isfile(infile):
                arch = Archive(infile=infile, files=files)
                specs = arch.get_node("comp_archive_spec", {"compname":comp})
            else:
                if infile is None:
                    logger.debug("No archive file defined for component %s"%comp)
                else:
                    logger.debug("Archive file %s for component %s not found"%(infile,comp))

                specs = self.get_optional_node("comp_archive_spec", attributes={"compname":comp})
            if specs is None:
                logger.debug("No archive specs found for component %s"%comp)
            else:
                logger.debug("adding archive spec for %s"%comp)
                env_archive.add_child(specs)

