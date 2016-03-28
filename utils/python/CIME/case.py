"""
Wrapper around all env XML for a case.

All interaction with and between the module files in XML/ takes place
through the Case module.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, run_cmd
from CIME.XML.machines          import Machines
from CIME.XML.files             import Files
from CIME.XML.component         import Component
from CIME.XML.compsets          import Compsets
from CIME.XML.grids             import Grids

from CIME.XML.env_test          import EnvTest
from CIME.XML.env_mach_specific import EnvMachSpecific
from CIME.XML.env_case          import EnvCase
from CIME.XML.env_mach_pes      import EnvMachPes
from CIME.XML.env_build         import EnvBuild
from CIME.XML.env_run           import EnvRun
from CIME.XML.env_archive       import EnvArchive
from CIME.XML.env_batch         import EnvBatch

logger = logging.getLogger(__name__)

class Case(object):

    def __init__(self, case_root=os.getcwd()):
        expect(os.path.isdir(case_root),
               "Case root directory '%s' does not exist" % case_root)

        self._env_files_that_need_rewrite = set()

        self._env_entryid_files = []
        self._env_entryid_files.append(EnvRun(case_root))
        self._env_entryid_files.append(EnvBuild(case_root))
        self._env_entryid_files.append(EnvMachPes(case_root))
        self._env_entryid_files.append(EnvCase(case_root))
        self._env_entryid_files.append(EnvBatch(case_root))
        if os.path.isfile(os.path.join(case_root,"env_test.xml")):
            self._env_entryid_files.append(EnvTest(case_root))
        self._env_generic_files = []
        self._env_generic_files.append(EnvMachSpecific(case_root))
        self._env_generic_files.append(EnvArchive(case_root))

        self._case_root = case_root
        self._target_component = None
        self._compset = None
        self._grid = None
        self._components = []
        self._component_classes = []
        self._component_grids = []
        self._gridmaps = []

    def __del__(self):
        self.flush()

    def flush(self):
        for env_file in self._env_files_that_need_rewrite:
            env_file.write()

        self._env_files_that_need_rewrite = set()

    def get_value(self, item, attribute={}, resolved=True, subgroup=None):
        result = None
        for env_file in self._env_entryid_files:
            # Wait and resolve in self rather than in env_file
            result = env_file.get_value(item, attribute, resolved=False, subgroup=subgroup)
            logging.debug("CASE %s %s"%(item,result))
            if result is not None:
                if resolved and type(result) is str:
                    return self.get_resolved_value(result)
                return result

        logging.info("Not able to retreive value for item '%s'" % item)

    def get_type_info(self, item):
        result = None
        for env_file in self._env_entryid_files:
            result = env_file.get_type_info(item)
            if result is not None:
                return result

        logging.info("Not able to retreive type for item '%s'" % item)

    def get_resolved_value(self, item, recurse=0):
        num_unresolved = item.count("$")
        recurse_limit = 10
        if (num_unresolved > 0 and recurse < recurse_limit ):
            for env_file in self._env_entryid_files:
                result = env_file.get_resolved_value(item)
                item = result
            if ("$" not in item):
                return item
            else:
                self.get_resolved_value(item,recurse=recurse+1)

        if(recurse >= recurse_limit):
            logging.warning("Not able to fully resolve item '%s'" % item)

        return item

    def set_value(self, item, value, subgroup=None, ignore_type=False):
        """
        If a file has not been defined, set an id/value pair in the
        case dictionary, this will be used later. Note that in
        create_newcase, when this is called and are setting the
       command line options none of these files have been defined
        If a file has been defined, and the variable is in the file,
        then that value will be set in the file object and the file
        name is returned
        """
        result = None;
        for env_file in self._env_entryid_files:
            result = env_file.set_value(item, value, subgroup, ignore_type)
            if (result is not None):
                self._env_files_that_need_rewrite.add(env_file)
                return result

    def _get_compset_longname(self, compset_name):
        """
        Find the compset longname and the target component that sets
        the compset, given the input compset name (which could be a
        long name or an alias) component whose compset file contained
        an entry that matched the compset longname
        """
        files = Files()
        components = files.get_components("COMPSETS_SPEC_FILE")
        logger.info(" Possible components for COMPSETS_SPEC_FILE are %s" % components)
        
        # Loop through all of the files listed in COMPSETS_SPEC_FILE and find the file
        # that has a match for either the alias or the longname in that order
        for component in components:
            
            # Determine the compsets file for the possible target component
            file = files.get_value("COMPSETS_SPEC_FILE", {"component":component})
            
            # If the file exists, read it and see if there is a match for the compset alias or longname
            if (os.path.isfile(file)):
                compsets = Compsets(file)
                match = compsets.get_compset_match(name=compset_name)
                if match is not None:
                    self._compset = match
                    self._target_component = component
                    logger.debug("Successful compset match %s found in file %s " %(self._compset, file))
                    return 

        logger.debug("Could not find a compset match for either alias or longname in %s" %file)

    def _get_grid_info(self, grid_name):
        files = Files()
        gridfile = files.get_value("GRIDS_SPEC_FILE")
        logger.debug(" Grid specification file is %s" % gridfile)
        
        grids = Grids(gridfile)
        self._grid, self._component_grids, self._domains, self._gridmaps \
            = grids.get_grid_match(name=grid_name, compset=self._compset)

        self.set_value("GRID", self._grid)
        for idx, grid in enumerate(self._component_grids):
            self.set_value(grid[0].upper(), grid[1])
            print "DEBUG1: name is %s and set value is %s" %(grid[0].upper(),grid[1])
        for idx, map in enumerate(self._gridmaps):
            self.set_value(map[0],map[1])

        expect ((self._grid is not None),
                "No valid grid found for input grid %s and compset %s" %(self._grid, self._compset))

    def _get_compset_components(self):
        elements = self._compset.split('_')
        for element in elements:
            # ignore the initial date in the compset longname
            if re.search(r'^\d+$',element):
                pass
            else:
                element_component = element.split('%')[0].lower()
                element_component = re.sub(r'[0-9]*',"",element_component)
                self._components.append(element_component)
                
    def __iter__(self):
        for entryid_file in self._env_entryid_files:
            for key, val in entryid_file:
                if type(val) is str and '$' in val:
                    yield key, self.get_resolved_value(val)
                else:
                    yield key, val

    def _get_component_config_data(self):

        # attributes used for multi valued defaults ($attlist is a hash reference)
        attlist = {"component":self._target_component, "compset":self._compset, "grid":self._grid}

        for env_file in self._env_entryid_files:
            env_file.write()

        # Determine list of component classes that this coupler/driver knows how
        # to deal with. This list follows the same order as compset longnames follow.
        files = Files()
        drv_config_file = files.get_value("CONFIG_DRV_FILE")
        drv_comp = Component(drv_config_file)
        for env_file in self._env_entryid_files:
            print " DEBUG: filename is ",os.path.basename(env_file.filename)
            #env_file.write() # works here
            nodes = env_file.add_elements_by_group(drv_comp, attlist, os.path.basename(env_file.filename));
            print "DEBUG: number of nodes is ",len(nodes)
            for node in nodes:
                env_file.root.append(node)
            env_file.write() # does not work here
            

        # # loop over all elements of both component_classes and components - and get config_component_file for
        # # for each component
        # self._component_classes = drv_comp.get_valid_model_components()
        # for i in xrange(1,len(self._component_classes)):
        #     comp_class = self._component_classes[i]
        #     comp_name  = self._components[i-1]
	#     node_name = 'CONFIG_' + comp_class + '_FILE';
        #     comp_config_file = files.get_value(node_name, {"component":comp_name}, resolved=True) 
        #     comp_comp = Component(comp_config_file)
        #     print "DEBUG: ",i,comp_class,comp_name,node_name,comp_config_file
        #     for env_file in self._env_entryid_files:
        #         env_file.add_elements_by_group(comp_comp, attlist, os.path.basename(env_file.filename));


    def configure(self, compset_name, grid_name):

        self._get_compset_longname(compset_name)
        self._get_compset_components()
        self._get_grid_info(grid_name)
        self._get_component_config_data()

        logger.info(" Component that sets compsets is: %s" %self._target_component)
        logger.info(" Compset is: %s " %self._compset)
        logger.info(" Grid is: %s " %self._grid )
        logger.info(" Components in compset are: %s " %self._components)
        logger.info(" Component grids in compset are: %s " %self._component_grids)
>>>>>>> first pass at create_newcase
