"""
Wrapper around all env XML for a case.

All interaction with and between the module files in XML/ takes place
through the Case module.
"""
from copy   import deepcopy
import glob, os, shutil, traceback
from CIME.XML.standard_module_setup import *

from CIME.utils                     import expect, get_cime_root, append_status
from CIME.utils                     import convert_to_type, get_model, get_project
from CIME.XML.build                 import Build
from CIME.XML.machines              import Machines
from CIME.XML.pes                   import Pes
from CIME.XML.files                 import Files
from CIME.XML.component             import Component
from CIME.XML.compsets              import Compsets
from CIME.XML.grids                 import Grids
from CIME.XML.batch                 import Batch
from CIME.XML.pio                   import PIO

from CIME.XML.env_test              import EnvTest
from CIME.XML.env_mach_specific     import EnvMachSpecific
from CIME.XML.env_case              import EnvCase
from CIME.XML.env_mach_pes          import EnvMachPes
from CIME.XML.env_build             import EnvBuild
from CIME.XML.env_run               import EnvRun
from CIME.XML.env_archive           import EnvArchive
from CIME.XML.env_batch             import EnvBatch

from CIME.user_mod_support          import apply_user_mods
from CIME.case_setup import case_setup

logger = logging.getLogger(__name__)

class Case(object):
    """
    https://github.com/ESMCI/cime/wiki/Developers-Introduction
    The Case class is the heart of the CIME Case Control system.  All
    interactions with a Case take part through this class.  All of the
    variables used to create and manipulate a case are defined in xml
    files and for every xml file there is a python class to interact
    with that file.

    XML files which are part of the CIME distribution and are meant to
    be readonly with respect to a case are typically named
    config_something.xml and the corresponding python Class is
    Something and can be found in file CIME.XML.something.py.  I'll
    refer to these as the CIME config classes.

    XML files which are part of a case and thus are read/write to a
    case are typically named env_whatever.xml and the cooresponding
    python modules are CIME.XML.env_whatever.py and classes are
    EnvWhatever.  I'll refer to these as the Case env classes.

    The Case Class includes an array of the Case env classes, in the
    configure function and it's supporting functions defined below
    the case object creates and manipulates the Case env classes
    by reading and interpreting the CIME config classes.

    """
    def __init__(self, case_root=None, read_only=True):

        if case_root is None:
            case_root = os.getcwd()
        self._caseroot = case_root
        logger.debug("Initializing Case.")
        self._env_files_that_need_rewrite = set()
        self._read_only_mode = True
        self._force_read_only = read_only

        self._env_entryid_files = []
        self._env_generic_files = []
        self._files = []

        self.read_xml()

        # Hold arbitary values. In create_newcase we may set values
        # for xml files that haven't been created yet. We need a place
        # to store them until we are ready to create the file. At file
        # creation we get the values for those fields from this lookup
        # table and then remove the entry.
        self.lookups = {}
        self.set_lookup_value('CIMEROOT',os.path.abspath(get_cime_root()))

        self._compsetname = None
        self._gridname = None
        self._compsetsfile = None
        self._pesfile = None
        self._gridfile = None
        self._components = []
        self._component_classes = []

    # Define __enter__ and __exit__ so that we can use this as a context manager
    # and force a flush on exit.
    def __enter__(self):
        if not self._force_read_only:
            self._read_only_mode = False
        return self

    def __exit__(self, *_):
        self.flush()
        self._read_only_mode = True
        return False

    def schedule_rewrite(self, env_file):
        assert not self._read_only_mode, \
            "case.py scripts error: attempted to modify an env file while in " \
            "read-only mode"
        self._env_files_that_need_rewrite.add(env_file)

    def read_xml(self):
        if(len(self._env_files_that_need_rewrite)>0):
            files = ""
            for env_file in self._env_files_that_need_rewrite:
                files += " "+env_file.filename
            expect(False,"Object(s) %s seem to have newer data than the corresponding case file"%files)

        self._env_entryid_files = []
        self._env_entryid_files.append(EnvRun(self._caseroot))
        self._env_entryid_files.append(EnvBuild(self._caseroot))
        self._env_entryid_files.append(EnvMachPes(self._caseroot))
        self._env_entryid_files.append(EnvCase(self._caseroot))
        self._env_entryid_files.append(EnvBatch(self._caseroot))
        if os.path.isfile(os.path.join(self._caseroot,"env_test.xml")):
            self._env_entryid_files.append(EnvTest(self._caseroot))
        self._env_generic_files = []
        self._env_generic_files.append(EnvMachSpecific(self._caseroot))
        self._env_generic_files.append(EnvArchive(self._caseroot))
        self._files = self._env_entryid_files + self._env_generic_files

    def get_env(self, short_name):
        full_name = "env_%s.xml" % (short_name)
        for env_file in self._files:
            if os.path.basename(env_file.filename) == full_name:
                return env_file

        expect(False, "Could not find object for %s in case"%full_name)

    def copy(self, newcasename, newcaseroot, newcimeroot=None, newsrcroot=None):
        newcase = deepcopy(self)
        for env_file in newcase._files: # pylint: disable=protected-access
            basename = os.path.basename(env_file.filename)
            env_file.filename = os.path.join(newcaseroot,basename)

        if newcimeroot is not None:
            newcase.set_value("CIMEROOT", newcimeroot)

        if newsrcroot is not None:
            newcase.set_value("SRCROOT", newsrcroot)

        newcase.set_value("CASE",newcasename)
        newcase.set_value("CASEROOT",newcaseroot)
        newcase.set_value("CONTINUE_RUN","FALSE")
        newcase.set_value("RESUBMIT",0)
        return newcase

    def flush(self, flushall=False):
        if not os.path.isdir(self._caseroot):
            # do not flush if caseroot wasnt created
            return
        if flushall:
            for env_file in self._files:
                self.schedule_rewrite(env_file)
        for env_file in self._env_files_that_need_rewrite:
            env_file.write()
        self._env_files_that_need_rewrite = set()

    def get_value(self, item, attribute=None, resolved=True, subgroup=None):
        result = None
        for env_file in self._env_entryid_files:
            # Wait and resolve in self rather than in env_file

            result = env_file.get_value(item, attribute, resolved=False, subgroup=subgroup)

            if result is not None:
                if resolved and type(result) is str:
                    result = self.get_resolved_value(result)
                    vtype = env_file.get_type_info(item)
                    result = convert_to_type(result, vtype, item)
                return result

        for env_file in self._env_generic_files:

            result = env_file.get_value(item, attribute, resolved=False, subgroup=subgroup)

            if result is not None:
                if resolved and type(result) is str:
                    return self.get_resolved_value(result)
                return result

        # Return empty result
        return result


    def get_values(self, item=None, attribute=None, resolved=True, subgroup=None):

        """
        Return info object for given item, return all info for all item if item is empty.
        """

        logger.debug("(get_values) Input values: %s , %s , %s , %s , %s" , self.__class__.__name__ , item, attribute, resolved, subgroup)

        # Empty result list
        results = []

        for env_file in self._env_entryid_files:
            # Wait and resolve in self rather than in env_file
            logger.debug("Searching in %s" , env_file.__class__.__name__)
            result = None

            try:
                # env_batch has its own implementation of get_values otherwise in entry_id
                result = env_file.get_values(item, attribute, resolved=False, subgroup=subgroup)
                # Method exists, and was used.
            except AttributeError:
                # Method does not exist.  What now?
                traceback.print_exc()
                logger.debug("No get_values method for class %s (%s)" , env_file.__class__.__name__ , AttributeError)

            if result is not None and (len(result) >= 1):

                if resolved :
                    for r in result :
                        if type(r['value']) is str:
                            logger.debug("Resolving %s" , r['value'])
                            r['value'] = self.get_resolved_value(r['value'])

                if subgroup :
                    found = []
                    for r in result :
                        if r['group'] == subgroup :
                            found.append(r)
                    results += found
                else:
                    results = results + result

        logger.debug("(get_values) Return value:  %s" , results )
        return results

    def get_type_info(self, item):
        result = None
        for env_file in self._env_entryid_files:
            result = env_file.get_type_info(item)
            if result is not None:
                return result

        logging.debug("Not able to retreive type for item '%s'" % item)

    def get_resolved_value(self, item, recurse=0):
        num_unresolved = item.count("$")
        recurse_limit = 10
        if (num_unresolved > 0 and recurse < recurse_limit ):
            for env_file in self._env_entryid_files:
                item = env_file.get_resolved_value(item)
            if ("$" not in item):
                return item
            else:
                item = self.get_resolved_value(item,recurse=recurse+1)

        if recurse >= 2*recurse_limit:
            logging.warning("Not able to fully resolve item '%s'" % item)
        elif recurse >= recurse_limit:
            #try env_batch first
            env_batch = self.get_env("batch")
            item = env_batch.get_resolved_value(item)
            logger.debug("item is %s, checking env_batch"%item)
            if item is not None:
                if ("$" not in item):
                    return item
                else:
                    item = self.get_resolved_value(item,recurse=recurse+1)
            else:
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
        if item == "CASEROOT":
            self._caseroot = value
        result = None
        for env_file in self._env_entryid_files:
            result = env_file.set_value(item, value, subgroup, ignore_type)
            if (result is not None):
                logger.debug("Will rewrite file %s %s",env_file.filename, item)
                self._env_files_that_need_rewrite.add(env_file)
                return result

    def set_lookup_value(self, item, value):
        if item in self.lookups.keys() and self.lookups[item] is not None:
            logger.warn("Item %s already in lookups with value %s"%(item,self.lookups[item]))
        else:
            self.lookups[item] = value


    def _set_compset_and_pesfile(self, compset_name, user_compset=False, pesfile=None):
        """
        Loop through all the compset files and find the compset
        specifation file that matches either the input 'compset_name'.
        Note that the input compset name (i.e. compset_name) can be
        either a longname or an alias.  This will also set the
        compsets and pes specfication files.
        """
        files = Files()
        components = files.get_components("COMPSETS_SPEC_FILE")
        logger.debug(" Possible components for COMPSETS_SPEC_FILE are %s" % components)

        # Loop through all of the files listed in COMPSETS_SPEC_FILE and find the file
        # that has a match for either the alias or the longname in that order
        for component in components:

            # Determine the compsets file for this component
            compsets_filename = files.get_value("COMPSETS_SPEC_FILE", {"component":component})

            # If the file exists, read it and see if there is a match for the compset alias or longname
            if (os.path.isfile(compsets_filename)):
                compsets = Compsets(compsets_filename)
                match = compsets.get_compset_match(name=compset_name)
                pesfile = files.get_value("PES_SPEC_FILE"     , {"component":component})
                if match is not None:
                    self._pesfile = pesfile
                    self._compsetsfile = compsets_filename
                    self._compsetname = match
                    tests_filename    = files.get_value("TESTS_SPEC_FILE"   , {"component":component}, resolved=False)
                    tests_mods_dir    = files.get_value("TESTS_MODS_DIR"    , {"component":component}, resolved=False)
                    user_mods_dir     = files.get_value("USER_MODS_DIR"     , {"component":component}, resolved=False)
                    self.set_lookup_value("COMPSETS_SPEC_FILE" ,
                                   files.get_value("COMPSETS_SPEC_FILE", {"component":component}, resolved=False))
                    self.set_lookup_value("TESTS_SPEC_FILE"    , tests_filename)
                    self.set_lookup_value("TESTS_MODS_DIR"     , tests_mods_dir)
                    self.set_lookup_value("USER_MODS_DIR"      , user_mods_dir)
                    self.set_lookup_value("PES_SPEC_FILE"      ,
                                   files.get_value("PES_SPEC_FILE"     , {"component":component}, resolved=False))
                    logger.info("Compset longname is %s " %(match))
                    logger.info("Compset specification file is %s" %(compsets_filename))
                    logger.info("Pes     specification file is %s" %(pesfile))
                    return

        if user_compset is True:
            #Do not error out for user_compset
            logger.warn("Could not find a compset match for either alias or longname in %s" %(compset_name))
            self._compsetname = compset_name
            self._pesfile = pesfile
            self.set_lookup_value("PES_SPEC_FILE", pesfile)
        else:
            expect(False,
                   "Could not find a compset match for either alias or longname in %s" %(compset_name))


    def get_compset_components(self):
        #If are doing a create_clone then, self._compsetname is not set yet
        components = []
        compset = self.get_value("COMPSET")
        if compset is None:
            compset = self._compsetname
        expect(compset is not None,
               "ERROR: compset is not set")
        # the first element is always the date operator - skip it
        elements = compset.split('_')[1:] # pylint: disable=maybe-no-member
        for element in elements:
            # ignore the possible BGC or TEST modifier
            if element.startswith("BGC%") or element.startswith("TEST"):
                continue
            else:
                element_component = element.split('%')[0].lower()
                element_component = re.sub(r'[0-9]*',"",element_component)
                components.append(element_component)
        return components


    def __iter__(self):
        for entryid_file in self._env_entryid_files:
            for key, val in entryid_file:
                if type(val) is str and '$' in val:
                    yield key, self.get_resolved_value(val)
                else:
                    yield key, val


    def _get_component_config_data(self):
        # attributes used for multi valued defaults ($attlist is a hash reference)
        attlist = {"compset":self._compsetname, "grid":self._gridname}

        # Determine list of component classes that this coupler/driver knows how
        # to deal with. This list follows the same order as compset longnames follow.
        files = Files()
        # Add the group and elements for the config_files.xml
        for env_file in self._env_entryid_files:
            env_file.add_elements_by_group(files, attlist)

        drv_config_file = files.get_value("CONFIG_DRV_FILE")
        drv_comp = Component(drv_config_file)
        for env_file in self._env_entryid_files:
            env_file.add_elements_by_group(drv_comp, attributes=attlist)

        # loop over all elements of both component_classes and components - and get config_component_file for
        # for each component
        self._component_classes =drv_comp.get_valid_model_components()
        if len(self._component_classes) > len(self._components):
            self._components.append('sesp')

        for i in xrange(1,len(self._component_classes)):
            comp_class = self._component_classes[i]
            comp_name  = self._components[i-1]
            node_name = 'CONFIG_' + comp_class + '_FILE'
            # Add the group and elements for the config_files.xml
            comp_config_file = files.get_value(node_name, {"component":comp_name}, resolved=False)
            self.set_value(node_name, comp_config_file)
            comp_config_file = self.get_resolved_value(comp_config_file)
            expect(comp_config_file is not None,"No config file for component %s"%comp_name)
            compobj = Component(comp_config_file)
            for env_file in self._env_entryid_files:
                env_file.add_elements_by_group(compobj, attributes=attlist)


        for key,value in self.lookups.items():
            result = self.set_value(key,value)
            if result is not None:
                del self.lookups[key]

    def get_components(self):
        """
        return dictionary of the form [component_class:component],
        e.g. [atm:cam], for all compset components
        """

        files = Files()
        drv_comp = Component(files.get_value("CONFIG_DRV_FILE"))

        # Determine list of component classes that this coupler/driver knows how
        # to deal with. This list follows the same order as compset longnames follow.
        component_classes = drv_comp.get_valid_model_components()
        components = self.get_compset_components()

        # Note that component classes can have a bigger range than
        # compents since stub esp (sesp) is an optional component - so
        # need to take the min of the two below
        comp_dict = {}
        for i in xrange(0,len(components)):
            comp_name  = components[i]
            comp_class = component_classes[i+1]
            comp_dict[comp_class] = comp_name
        return comp_dict

    def configure(self, compset_name, grid_name, machine_name=None,
                  project=None, pecount=None, compiler=None, mpilib=None,
                  user_compset=False, pesfile=None,
                  user_grid=False, gridfile=None, ninst=1, test=False,
                  walltime=None, queue=None):

        #--------------------------------------------
        # compset, pesfile, and compset components
        #--------------------------------------------
        self._set_compset_and_pesfile(compset_name, user_compset=user_compset, pesfile=pesfile)

        self._components = self.get_compset_components()
        #FIXME - if --user-compset is True then need to determine that
        #all of the compset settings are valid

        #--------------------------------------------
        # grid
        #--------------------------------------------
        if user_grid is True and gridfile is not None:
            self.set_value("GRIDS_SPEC_FILE", gridfile)
        grids = Grids(gridfile)

        gridinfo = grids.get_grid_info(name=grid_name, compset=self._compsetname)

        self._gridname = gridinfo["GRID"]
        for key,value in gridinfo.items():
            logger.debug("Set grid %s %s"%(key,value))
            self.set_lookup_value(key,value)

        #--------------------------------------------
        # component config data
        #--------------------------------------------
        self._get_component_config_data()

        self.get_compset_var_settings()

        #--------------------------------------------
        # machine
        #--------------------------------------------
        # set machine values in env_xxx files
        machobj = Machines(machine=machine_name)
        machine_name = machobj.get_machine_name()
        self.set_value("MACH",machine_name)
        nodenames = machobj.get_node_names()
        nodenames =  [x for x in nodenames if
                      '_system' not in x and '_variables' not in x and 'mpirun' not in x and\
                      'COMPILER' not in x and 'MPILIB' not in x]

        for nodename in nodenames:
            value = machobj.get_value(nodename, resolved=False)
            type_str = self.get_type_info(nodename)
            if type_str is not None:
                logger.debug("machine nodname %s value %s"%(nodename, value))
                self.set_value(nodename, convert_to_type(value, type_str, nodename))

        if compiler is None:
            compiler = machobj.get_default_compiler()
        else:
            expect(machobj.is_valid_compiler(compiler),
                   "compiler %s is not supported on machine %s" %(compiler, machine_name))

        self.set_value("COMPILER",compiler)

        if mpilib is None:
            mpilib = machobj.get_default_MPIlib({"compiler":compiler})
        else:
            expect(machobj.is_valid_MPIlib(mpilib, {"compiler":compiler}),
                   "MPIlib %s is not supported on machine %s" %(mpilib, machine_name))
        self.set_value("MPILIB",mpilib)

        machdir = machobj.get_machines_dir()
        self.set_value("MACHDIR", machdir)

        # Create env_mach_specific settings from machine info.
        env_mach_specific_obj = self.get_env("mach_specific")
        env_mach_specific_obj.populate(machobj)
        self.schedule_rewrite(env_mach_specific_obj)

        #--------------------------------------------
        # pe payout
        #--------------------------------------------
        match1 = re.match('([0-9]+)x([0-9]+)', "" if pecount is None else pecount)
        match2 = re.match('([0-9]+)', "" if pecount is None else pecount)
        pes_ntasks = {}
        pes_nthrds = {}
        pes_rootpe = {}
        if match1:
            opti_tasks = match1.group(1)
            opti_thrds = match1.group(2)
        elif match2:
            opti_tasks = match2.group(1)
            opti_thrds = 1

        other = {}
        if match1 or match2:
            for component_class in self._component_classes:
                if component_class == "DRV":
                    component_class = "CPL"
                string = "NTASKS_" + component_class
                pes_ntasks[string] = opti_tasks
                string = "NTHRDS_" + component_class
                pes_nthrds[string] = opti_thrds
                string = "ROOTPE_" + component_class
                pes_rootpe[string] = 0
        else:
            pesobj = Pes(self._pesfile)

            pes_ntasks, pes_nthrds, pes_rootpe, other = pesobj.find_pes_layout(self._gridname, self._compsetname,
                                                                    machine_name, pesize_opts=pecount)

        mach_pes_obj = self.get_env("mach_pes")
        totaltasks = {}
        # Since other items may include PES_PER_NODE we need to do this first
        # we can get rid of this code when all of the perl is removed
        for key, value in other.items():
            self.set_value(key, value)
        pes_per_node = self.get_value("PES_PER_NODE")
        for key, value in pes_ntasks.items():
            totaltasks[key[-3:]] = int(value)
            mach_pes_obj.set_value(key,int(value), pes_per_node=pes_per_node)
        for key, value in pes_rootpe.items():
            totaltasks[key[-3:]] += int(value)
            mach_pes_obj.set_value(key,int(value), pes_per_node=pes_per_node)
        for key, value in pes_nthrds.items():
            totaltasks[key[-3:]] *= int(value)
            mach_pes_obj.set_value(key,int(value), pes_per_node=pes_per_node)

        maxval = 1
        if mpilib != "mpi-serial":
            for key, val in totaltasks.items():
                if val < 0:
                    val = -1*val*pes_per_node
                if val > maxval:
                    maxval = val

        # Make sure that every component has been accounted for
        # set, nthrds and ntasks to 1 otherwise. Also set the ninst values here.
        for compclass in self._component_classes:
            if compclass == "DRV":
                continue
            key = "NINST_%s"%compclass
            mach_pes_obj.set_value(key, ninst)
            key = "NTASKS_%s"%compclass
            if key not in pes_ntasks.keys():
                mach_pes_obj.set_value(key,1)
            key = "NTHRDS_%s"%compclass
            if compclass not in pes_nthrds.keys():
                mach_pes_obj.set_value(compclass,1)

        # FIXME - this is a short term fix for dealing with the restriction that
        # CISM1 cannot run on multiple cores
        if "CISM1" in self._compsetname:
            mach_pes_obj.set_value("NTASKS_GLC",1)
            mach_pes_obj.set_value("NTHRDS_GLC",1)

        #--------------------------------------------
        # batch system
        #--------------------------------------------
        batch_system_type = machobj.get_value("BATCH_SYSTEM")
        batch = Batch(batch_system=batch_system_type, machine=machine_name)
        bjobs = batch.get_batch_jobs()
        env_batch = self.get_env("batch")
        env_batch.set_batch_system(batch, batch_system_type=batch_system_type)
        env_batch.create_job_groups(bjobs)
        env_batch.set_job_defaults(bjobs, pesize=maxval, walltime=walltime, force_queue=queue)
        self.schedule_rewrite(env_batch)

        self.set_value("COMPSET",self._compsetname)

        self._set_pio_xml()
        logger.info(" Compset is: %s " %self._compsetname)
        logger.info(" Grid is: %s " %self._gridname )
        logger.info(" Components in compset are: %s " %self._components)

        # Set project id
        if project is None:
            project = get_project(machobj)
        if project is not None:
            self.set_value("PROJECT", project)
        elif machobj.get_value("PROJECT_REQUIRED"):
            expect(project is not None, "PROJECT_REQUIRED is true but no project found")

        # Overwriting an existing exeroot or rundir can cause problems
        exeroot = self.get_value("EXEROOT")
        rundir = self.get_value("RUNDIR")
        for wdir in (exeroot, rundir):
            logging.debug("wdir is %s"%wdir)
            if os.path.exists(wdir):
                expect(not test, "Directory %s already exists, aborting test"% wdir)
                response = raw_input("\nDirectory %s already exists, (r)eplace, (a)bort, or (u)se existing?"% wdir)
                if response.startswith("r"):
                    shutil.rmtree(wdir)
                else:
                    expect(response.startswith("u"), "Aborting by user request")

        # miscellaneous settings
        if self.get_value("RUN_TYPE") == 'hybrid':
            self.set_value("GET_REFCASE", True)

        # Turn on short term archiving as cesm default setting
        model = get_model()
        if model == "cesm" and not test:
            self.set_value("DOUT_S",True)



    def get_compset_var_settings(self):
        compset_obj = Compsets(infile=self.get_value("COMPSETS_SPEC_FILE"))
        matches = compset_obj.get_compset_var_settings(self._compsetname, self._gridname)
        for name, value in matches:
            if len(value) > 0:
                logger.debug("Compset specific settings: name is %s and value is %s"%(name,value))
                self.set_value(name, value)

    def set_initial_test_values(self):
        testobj = self.get_env("test")
        testobj.set_initial_values(self)

    def get_batch_jobs(self):
        batchobj = self.get_env("batch")
        return batchobj.get_jobs()

    def _set_pio_xml(self):
        pioobj = PIO()
        grid = self.get_value("GRID")
        compiler = self.get_value("COMPILER")
        mach = self.get_value("MACH")
        compset = self.get_value("COMPSET")
        mpilib = self.get_value("MPILIB")
        defaults = pioobj.get_defaults(grid=grid,compset=compset,mach=mach,compiler=compiler, mpilib=mpilib)
        for vid, value in defaults.items():
            self.set_value(vid,value)

    def _create_caseroot_tools(self):
        machines_dir = os.path.abspath(self.get_value("MACHDIR"))
        toolsdir = os.path.join(self.get_value("CIMEROOT"),"scripts","Tools")
        # setup executable files in caseroot/
        exefiles = (os.path.join(toolsdir, "case.setup"),
                    os.path.join(toolsdir, "case.build"),
                    os.path.join(toolsdir, "case.submit"),
                    os.path.join(toolsdir, "preview_namelists"),
                    os.path.join(toolsdir, "check_input_data"),
                    os.path.join(toolsdir, "check_case"),
                    os.path.join(toolsdir, "taskmaker"),
                    os.path.join(toolsdir, "archive_metadata.sh"),
                    os.path.join(toolsdir, "xmlchange"),
                    os.path.join(toolsdir, "xmlquery"))
        try:
            for exefile in exefiles:
                destfile = os.path.join(self._caseroot,os.path.basename(exefile))
                os.symlink(exefile, destfile)
        except Exception as e:
            logger.warning("FAILED to set up exefiles: %s" % str(e))

        # set up utility files in caseroot/Tools/
        toolfiles = (os.path.join(toolsdir, "check_lockedfiles"),
                     os.path.join(toolsdir, "lt_archive.sh"),
                     os.path.join(toolsdir, "getTiming"),
                     os.path.join(toolsdir, "save_provenance"),
                     os.path.join(machines_dir,"Makefile"),
                     os.path.join(machines_dir,"mkSrcfiles"),
                     os.path.join(machines_dir,"mkDepends"))

        for toolfile in toolfiles:
            destfile = os.path.join(self._caseroot,"Tools",os.path.basename(toolfile))
            expect(os.path.isfile(toolfile)," File %s does not exist"%toolfile)
            try:
                os.symlink(toolfile, destfile)
            except Exception as e:
                logger.warning("FAILED to set up toolfiles: %s %s %s" % (str(e), toolfile, destfile))

        # Create Macros file.
        machine = self.get_value("MACH")
        files = Files()
        # Use config_build if the environment variable is set, or if there is no
        # config_compilers file.
        if os.getenv("CIME_USE_CONFIG_BUILD") == "TRUE" or \
           files.get_value("COMPILERS_SPEC_FILE") is None:
            build_file = files.get_value("BUILD_SPEC_FILE")
            machobj = Machines(machine=machine, files=files)
            macro_maker = Build(machobj)
            macros_path = os.path.join(self._caseroot, "Macros")
            with open(macros_path, "w") as macros_file:
                macro_maker.write_macros('Makefile', build_file, macros_file)

        # Copy any system or compiler Depends files to the case.
        compiler = self.get_value("COMPILER")
        for dep in (machine, compiler):
            dfile = "Depends.%s"%dep
            if os.path.isfile(os.path.join(machines_dir,dfile)):
                shutil.copyfile(os.path.join(machines_dir,dfile), os.path.join(self._caseroot,dfile))
        dfile = "Depends.%s.%s"%(machine,compiler)
        if os.path.isfile(os.path.join(machines_dir,dfile)):
            shutil.copyfile(os.path.join(machines_dir,dfile), os.path.join(self._caseroot, dfile))
            # set up infon files
            # infofiles = os.path.join(os.path.join(toolsdir, README.post_process")
            #FIXME - the following does not work
            # print "DEBUG: infofiles are ",infofiles
            #    try:
            #        for infofile in infofiles:
            #            print "DEBUG: infofile is %s, %s"  %(infofile, os.path.basename(infofile))
            #            dst_file = caseroot + "/" + os.path.basename(infofile)
            #            shutil.copyfile(infofile, dst_file)
            #            os.chmod(dst_file, os.stat(dst_file).st_mode | stat.S_IXUSR | stat.S_IXGRP)
            #    except Exception as e:
            #        logger.warning("FAILED to set up infofiles: %s" % str(e))

    def _create_caseroot_sourcemods(self):
        components = self.get_compset_components()
        for component in components:
            directory = os.path.join(self._caseroot,"SourceMods","src.%s"%component)
            if not os.path.exists(directory):
                os.makedirs(directory)

        directory = os.path.join(self._caseroot, "SourceMods", "src.share")
        if not os.path.exists(directory):
            os.makedirs(directory)

        directory = os.path.join(self._caseroot,"SourceMods","src.drv")
        if not os.path.exists(directory):
            os.makedirs(directory)

        if get_model() == "cesm":
        # Note: this is CESM specific, given that we are referencing cism explitly
            if "cism" in components:
                directory = os.path.join(self._caseroot, "SourceMods", "src.cism", "glimmer-cism")
                if not os.path.exists(directory):
                    os.makedirs(directory)
                readme_file = os.path.join(directory, "README")

                str_to_write = """
                Put source mods for the glimmer-cism library in the glimmer-cism subdirectory
                This includes any files that are in the glimmer-cism subdirectory of $cimeroot/../components/cism
                Anything else (e.g., mods to source_glc or drivers) goes in this directory, NOT in glimmer-cism/"""

                with open(readme_file, "w") as fd:
                    fd.write(str_to_write)

    def create_caseroot(self, clone=False):
        if not os.path.exists(self._caseroot):
            # Make the case directory
            logger.info(" Creating Case directory %s" %self._caseroot)
            os.makedirs(self._caseroot)
        os.chdir(self._caseroot)

        # Create relevant directories in $self._caseroot
        if clone:
            newdirs = ("LockedFiles", "Tools")
        else:
            newdirs = ("SourceMods", "LockedFiles", "Buildconf", "Tools")
        for newdir in newdirs:
            os.makedirs(newdir)
        # Open a new README.case file in $self._caseroot

        append_status(" ".join(sys.argv), caseroot=self._caseroot, sfile="README.case")
        append_status("Compset longname is %s"%self.get_value("COMPSET"),
                      caseroot=self._caseroot, sfile="README.case")
        append_status("Compset specification file is %s" %
                      (self.get_value("COMPSETS_SPEC_FILE")),
                      caseroot=self._caseroot, sfile="README.case")
        append_status("Pes     specification file is %s" %
                      (self.get_value("PES_SPEC_FILE")),
                      caseroot=self._caseroot, sfile="README.case")
        for component_class in self._component_classes:
            if component_class == "DRV":
                continue
            comp_grid = "%s_GRID"%component_class
            append_status("%s is %s"%(comp_grid,self.get_value(comp_grid)),
                          caseroot=self._caseroot, sfile="README.case")
        if not clone:
            self._create_caseroot_sourcemods()
        self._create_caseroot_tools()

    def apply_user_mods(self, user_mods_dir=None):
        if user_mods_dir is not None:
            if os.path.isabs(user_mods_dir):
                user_mods_path = user_mods_dir
            else:
                user_mods_path = self.get_value('USER_MODS_DIR')
                user_mods_path = os.path.join(user_mods_path, user_mods_dir)
            self.set_value("USER_MODS_FULLPATH",user_mods_path)
            ninst_vals = {}
            for i in xrange(1,len(self._component_classes)):
                comp_class = self._component_classes[i]
                comp_name  = self._components[i-1]
                if comp_class == "DRV":
                    continue
                ninst_comp = self.get_value("NINST_%s"%comp_class)
                if ninst_comp > 1:
                    ninst_vals[comp_name] = ninst_comp
            apply_user_mods(self._caseroot, user_mods_path, ninst_vals)

    def create_clone(self, newcase, keepexe=False, mach_dir=None, project=None):

        newcaseroot = os.path.abspath(newcase)
        expect(not os.path.isdir(newcaseroot),
               "New caseroot directory %s already exists" % newcaseroot)
        newcasename = os.path.basename(newcaseroot)
        newcase_cimeroot = os.path.abspath(get_cime_root())

        # create clone from self to case
        clone_cimeroot = self.get_value("CIMEROOT")
        if newcase_cimeroot != clone_cimeroot:
            logger.warning(" case  CIMEROOT is %s " %newcase_cimeroot)
            logger.warning(" clone CIMEROOT is %s " %clone_cimeroot)
            logger.warning(" It is NOT recommended to clone cases from different versions of CIME.")


        # *** create case object as deepcopy of clone object ***
        srcroot = os.path.join(newcase_cimeroot,"..")
        newcase = self.copy(newcasename, newcaseroot, newsrcroot=srcroot)
        newcase.set_value("CIMEROOT", newcase_cimeroot)

        # determine if will use clone executable or not
        if keepexe:
            orig_exeroot = self.get_value("EXEROOT")
            newcase.set_value("EXEROOT", orig_exeroot)
            newcase.set_value("BUILD_COMPLETE","TRUE")
        else:
            newcase.set_value("BUILD_COMPLETE","FALSE")

        # set machdir
        if mach_dir is not None:
            newcase.set_value("MACHDIR", mach_dir)

        # Set project id
        # Note: we do not just copy this from the clone because it seems likely that
        # users will want to change this sometimes, especially when cloning another
        # user's case. However, note that, if a project is not given, the fallback will
        # be to copy it from the clone, just like other xml variables are copied.
        if project is None:
            project = self.get_value("PROJECT", subgroup="case.run")
        if project is not None:
            newcase.set_value("PROJECT", project)

        # create caseroot
        newcase.create_caseroot(clone=True)
        newcase.flush(flushall=True)

        # copy user_nl_files
        cloneroot = self._caseroot
        files = glob.glob(cloneroot + '/user_nl_*')
        for item in files:
            shutil.copy(item, newcaseroot)

        # copy SourceMod and Buildconf files
        for casesub in ("SourceMods", "Buildconf"):
            shutil.copytree(os.path.join(cloneroot, casesub), os.path.join(newcaseroot, casesub))

        # copy env_case.xml to LockedFiles
        shutil.copy(os.path.join(newcaseroot,"env_case.xml"), os.path.join(newcaseroot,"LockedFiles"))

        # Update README.case
        fclone   = open(cloneroot + "/README.case", "r")
        fnewcase = open(newcaseroot  + "/README.case", "a")
        fnewcase.write("\n    *** original clone README follows ****")
        fnewcase.write("\n " +  fclone.read())

        clonename = self.get_value("CASE")
        logger.info(" Successfully created new case %s from clone case %s " %(newcasename, clonename))

        case_setup(newcase, clean=False, test_mode=False)

        return newcase

    def submit_jobs(self, no_batch=False, job=None):
        env_batch = self.get_env('batch')
        env_batch.submit_jobs(self, no_batch=no_batch, job=job)
