"""
Interface to the archive.xml file.  This class inherits from GenericXML.py
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.archive_base import ArchiveBase
from CIME.XML.files import Files
from CIME.utils import expect, get_model
from copy import deepcopy

logger = logging.getLogger(__name__)

class Archive(ArchiveBase):

    def __init__(self, infile=None, files=None):
        """
        initialize an object
        """
        if files is None:
            files = Files()
        schema = files.get_schema("ARCHIVE_SPEC_FILE")
        super(Archive, self).__init__(infile, schema)

    def setup(self, env_archive, components, files=None):
        if files is None:
            files = Files()

        components_node = env_archive.make_child("components", attributes={"version":"2.0"})
        arch_components = deepcopy(components)
        model = get_model()
        if 'drv' not in arch_components and model != 'ufs':
            arch_components.append('drv')
        if 'dart' not in arch_components and model == 'cesm':
            arch_components.append('dart')

        for comp in arch_components:
            infile = files.get_value("ARCHIVE_SPEC_FILE", {"component":comp})

            if infile is not None and os.path.isfile(infile):
                arch = Archive(infile=infile, files=files)
                specs = arch.get_optional_child(name="comp_archive_spec", attributes={"compname":comp})
            else:
                if infile is None:
                    logger.debug("No archive file defined for component {}".format(comp))
                else:
                    logger.debug("Archive file {} for component {} not found".format(infile,comp))

                specs = self.get_optional_child(name="comp_archive_spec", attributes={"compname":comp})

            if specs is None:
                logger.debug("No archive specs found for component {}".format(comp))
            else:
                logger.debug("adding archive spec for {}".format(comp))
                env_archive.add_child(specs, root=components_node)

    def get_all_config_archive_files(self, files):
        """
        Returns the list of ARCHIVE_SPEC_FILES that exist on disk as defined in config_files.xml
        """
        archive_spec_node = files.get_child("entry", {"id" : "ARCHIVE_SPEC_FILE"})
        component_nodes = files.get_children("value", root=files.get_child("values", root=archive_spec_node))
        config_archive_files = []
        for comp in component_nodes:
            attr = self.get(comp,"component")
            if attr:
                compval = files.get_value("ARCHIVE_SPEC_FILE", attribute={"component":attr})
            else:
                compval = self.get_resolved_value(self.text(comp))

            if os.path.isfile(compval):
                config_archive_files.append(compval)

        config_archive_files = list(set(config_archive_files))
        return config_archive_files
