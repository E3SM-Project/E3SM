"""
Wrapper around all env XML for a case.

All interaction with and between the module files in XML/ takes place
through the Case module.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils                     import expect, run_cmd, get_cime_root, convert_to_type
from CIME.XML.machines              import Machines
from CIME.XML.pes                   import Pes
from CIME.XML.files                 import Files
from CIME.XML.component             import Component
from CIME.XML.compsets              import Compsets
from CIME.XML.grids                 import Grids

from CIME.XML.env_test              import EnvTest
from CIME.XML.env_mach_specific     import EnvMachSpecific
from CIME.XML.env_case              import EnvCase
from CIME.XML.env_mach_pes          import EnvMachPes
from CIME.XML.env_build             import EnvBuild
from CIME.XML.env_run               import EnvRun
from CIME.XML.env_archive           import EnvArchive
from CIME.XML.env_batch             import EnvBatch

from CIME.XML.generic_xml           import GenericXML

logger = logging.getLogger(__name__)

class Case(object):

    def __init__(self, case_root=os.getcwd()):
        #expect(os.path.isdir(case_root),
        #       "Case root directory '%s' does not exist" % case_root)

        #FIXME - define a method for when case_root does not exist

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

        # Hold arbitary values. In create_newcase we may set values
        # for xml files that haven't been created yet. We need a place
        # to store them until we are ready to create the file. At file
        # creation we get the values for those fields from this lookup
        # table and then remove the entry. This was what I came up
        # with in the perl anyway and I think that we still need it here.
        self.lookups = {}
        self.lookups['CIMEROOT'] = get_cime_root()

        self._compsetname = None
        self._target_component = None
        self._gridname = None
        self._components = []
        self._component_config_files = []

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

        result = None
        if item in self.lookups.keys():
            result = self.lookups[item]

        if result is None:
            logger.debug("No value available for item '%s'" % item)
        elif resolved:
            result = self.get_resolved_value(result)
            # Return value as right type
            type_str = self._get_type_info(node)
            return convert_to_type(result, type_str, item)

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
        if result is None:
            self.lookups[item] = value

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
                    self._compsetname = match
                    self._target_component = component
                    logger.debug("Successful compset match %s found in file %s " %(self._compsetname, file))
                    return

        logger.debug("Could not find a compset match for either alias or longname in %s" %file)

    def get_compset_components(self):
        elements = self._compsetname.split('_')
        for element in elements:
            # ignore the initial date in the compset longname
            if re.search(r'^\d+$',element):
                pass
            else:
                element_component = element.split('%')[0].lower()
                element_component = re.sub(r'[0-9]*',"",element_component)
                self._components.append(element_component)
        return self._components

    def __iter__(self):
        for entryid_file in self._env_entryid_files:
            for key, val in entryid_file:
                if type(val) is str and '$' in val:
                    yield key, self.get_resolved_value(val)
                else:
                    yield key, val

    def _get_component_config_data(self):
        # attributes used for multi valued defaults ($attlist is a hash reference)
        attlist = {"component":self._target_component, "compset":self._compsetname, "grid":self._gridname}

        # Determine list of component classes that this coupler/driver knows how
        # to deal with. This list follows the same order as compset longnames follow.
        files = Files()
        drv_config_file = files.get_value("CONFIG_DRV_FILE")
        drv_comp = Component(drv_config_file)
        for env_file in self._env_entryid_files:
            nodes = env_file.add_elements_by_group(drv_comp, attlist, os.path.basename(env_file.filename));

        # loop over all elements of both component_classes and components - and get config_component_file for
        # for each component
        component_classes =drv_comp.get_valid_model_components()
        for i in xrange(1,len(component_classes)):
            comp_class = component_classes[i]
            comp_name  = self._components[i-1]
	    node_name = 'CONFIG_' + comp_class + '_FILE';
            comp_config_file = files.get_value(node_name, {"component":comp_name}, resolved=True)
            component = Component(comp_config_file)
            for env_file in self._env_entryid_files:
                env_file.add_elements_by_group(component, attlist, os.path.basename(env_file.filename));
            self._component_config_files.append((node_name,comp_config_file))

        # Add the group and elements for the config_files.xml
        for env_file in self._env_entryid_files:
            env_file.add_elements_by_group(files, attlist, os.path.basename(env_file.filename));

        for key,value in self.lookups.items():
            result = self.set_value(key,value)
            if result is not None:
                del self.lookups[key]

    def configure(self, compset_name, grid_name, machine_name=None, compiler=None, mpilib=None):
        self._get_compset_longname(compset_name)
        self.get_compset_components()

        files = Files()

        # get grid info and overwrite files with that data
        if self.get_value("GRIDS_SPEC_FILE") is not None:
            gridfile = self.get_value("GRIDS_SPEC_FILE")
        else:
            gridfile = files.get_value("GRIDS_SPEC_FILE")
        grids = Grids(gridfile)
        gridinfo = grids.get_grid_info(name=grid_name, compset=self._compsetname)
        self._gridname = gridinfo["GRID"]
        logger.debug(" Grid specification file is %s" % gridfile)

        self._get_component_config_data()

        for key,value in gridinfo.iteritems():
            self.set_value(key,value)

        # Add the group and elements for the config_files.xml
        for idx, config_file in enumerate(self._component_config_files):
            self.set_value(config_file[0],config_file[1])

        #--------------------------------------------
        # machine
        #--------------------------------------------
        # set machine values in env_xxx files
        machobj = Machines(machine=machine_name)
        nodenames = machobj.get_node_names()

        if "COMPILER" in nodenames: nodenames.remove("COMPILER")
        if "MPILIB" in nodenames: nodenames.remove("MPILIB")
        nodenames =  [x for x in nodenames if
                      '_system' not in x and '_variables' not in x and 'mpirun' not in x]

        for nodename in nodenames:
            value = machobj.get_value(nodename)
            type_str = self.get_type_info(nodename)
            if type_str is not None:
                self.set_value(nodename, convert_to_type(value, type_str, nodename))

        if compiler is None:
            compiler = machobj.get_default_compiler()
        else:
            expect(machobj.is_valid_compiler(compiler),
                   "compiler %s is not supported on machine %s" %(compiler, machine_name))
        self.set_value("COMPILER",compiler)

        if mpilib is None:
            mpilib = machobj.get_default_MPIlib()
        else:
            expect(machobj.is_valid_MPIlib(mpilib),
                   "MPIlib %s is not supported on machine %s" %(mpilib, machine_name))
        self.set_value("MPILIB",mpilib)

        machdir = machobj.get_machines_dir()
        self.set_value("MACHDIR", machdir)

        # the following go into the env_mach_specific file
        vars = ("module_system", "environment_variables", "batch_system", "mpirun")
        for var in vars:
            nodes = machobj.get_first_child_nodes(var)
            for node in nodes:
                self._env_generic_files[0].add_child(node)

        #--------------------------------------------
        # batch system
        #--------------------------------------------
        batch_system = machobj.get_batch_system_type()
        #FIXME - this needs to be filled in - this is just a place holeder


        #--------------------------------------------
        # pe payout
        #--------------------------------------------
        pesobj = Pes(self._target_component)

        #FIXME - add pesize_opts as optional argument below
        pes_ntasks, pes_nthrds, pes_rootpe = pesobj.find_pes_layout(self._gridname, self._compsetname, machine_name)

        for key, value in pes_ntasks.iteritems():
            self.set_value(key,int(value))
        for key, value in pes_nthrds.iteritems():
            self.set_value(key,int(value))
        for key, value in pes_rootpe.iteritems():
            self.set_value(key,int(value))

        self.set_value("COMPSET",self._compsetname)

        logger.info(" Component that sets compsets is: %s" %self._target_component)
        logger.info(" Compset is: %s " %self._compsetname)
        logger.info(" Grid is: %s " %self._gridname )
        logger.info(" Components in compset are: %s " %self._components)
