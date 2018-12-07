"""
Wrapper around all env XML for a case.

All interaction with and between the module files in XML/ takes place
through the Case module.
"""
from copy import deepcopy
import glob, os, shutil, math, six
from CIME.XML.standard_module_setup import *
#pylint: disable=import-error,redefined-builtin
from six.moves import input
from CIME.utils                     import expect, get_cime_root, append_status
from CIME.utils                     import convert_to_type, get_model
from CIME.utils                     import get_project, get_charge_account, check_name
from CIME.utils                     import get_current_commit, safe_copy
from CIME.locked_files              import LOCKED_DIR, lock_file
from CIME.XML.machines              import Machines
from CIME.XML.pes                   import Pes
from CIME.XML.files                 import Files
from CIME.XML.testlist              import Testlist
from CIME.XML.component             import Component
from CIME.XML.compsets              import Compsets
from CIME.XML.grids                 import Grids
from CIME.XML.batch                 import Batch
from CIME.XML.pio                   import PIO
from CIME.XML.archive               import Archive
from CIME.XML.env_test              import EnvTest
from CIME.XML.env_mach_specific     import EnvMachSpecific
from CIME.XML.env_case              import EnvCase
from CIME.XML.env_mach_pes          import EnvMachPes
from CIME.XML.env_build             import EnvBuild
from CIME.XML.env_run               import EnvRun
from CIME.XML.env_archive           import EnvArchive
from CIME.XML.env_batch             import EnvBatch
from CIME.XML.generic_xml           import GenericXML
from CIME.user_mod_support          import apply_user_mods
from CIME.aprun import get_aprun_cmd_for_case

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

    This class extends across multiple files, class members external to this file
    are listed in the following imports
    """
    from CIME.case.case_setup import case_setup
    from CIME.case.case_clone import create_clone, _copy_user_modified_to_clone
    from CIME.case.case_test  import case_test
    from CIME.case.case_submit import check_DA_settings, check_case, submit
    from CIME.case.case_st_archive import case_st_archive, restore_from_archive, \
        archive_last_restarts, test_st_archive, test_env_archive
    from CIME.case.case_run import case_run
    from CIME.case.case_cmpgen_namelists import case_cmpgen_namelists
    from CIME.case.check_lockedfiles import check_lockedfile, check_lockedfiles, check_pelayouts_require_rebuild
    from CIME.case.preview_namelists import create_dirs, create_namelists
    from CIME.case.check_input_data import check_all_input_data, stage_refcase, check_input_data

    def __init__(self, case_root=None, read_only=True):

        if case_root is None:
            case_root = os.getcwd()
        self._caseroot = case_root
        logger.debug("Initializing Case.")
        self._read_only_mode = True
        self._force_read_only = read_only
        self._primary_component = None

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
        self._cime_model = get_model()
        self.set_lookup_value('MODEL', self._cime_model)
        self._compsetname = None
        self._gridname = None
        self._pesfile = None
        self._gridfile = None
        self._components = []
        self._component_classes = []
        self._component_description = {}
        self._is_env_loaded = False

        # these are user_mods as defined in the compset
        # Command Line user_mods are handled seperately

        # Derived attributes
        self.thread_count = None
        self.total_tasks = None
        self.tasks_per_node = None
        self.num_nodes = None
        self.spare_nodes = None
        self.tasks_per_numa = None
        self.cores_per_task = None
        self.srun_binding = None

        # check if case has been configured and if so initialize derived
        if self.get_value("CASEROOT") is not None:
            self.initialize_derived_attributes()

    def check_if_comp_var(self, vid):
        for env_file in self._env_entryid_files:
            new_vid, new_comp, iscompvar = env_file.check_if_comp_var(vid)
            if iscompvar:
                return new_vid, new_comp, iscompvar

        return vid, None, False

    def initialize_derived_attributes(self):
        """
        These are derived variables which can be used in the config_* files
        for variable substitution using the {{ var }} syntax
        """
        env_mach_pes  = self.get_env("mach_pes")
        env_mach_spec = self.get_env('mach_specific')
        comp_classes  = self.get_values("COMP_CLASSES")
        max_mpitasks_per_node  = self.get_value("MAX_MPITASKS_PER_NODE")

        self.thread_count = env_mach_pes.get_max_thread_count(comp_classes)

        mpi_attribs = {
            "compiler" : self.get_value("COMPILER"),
            "mpilib"   : self.get_value("MPILIB"),
            "threaded" : self.get_build_threaded(),
            }

        job = self.get_primary_job()
        executable = env_mach_spec.get_mpirun(self, mpi_attribs, job, exe_only=True)[0]
        if executable is not None and "aprun" in executable:
            _, self.num_nodes, self.total_tasks, self.tasks_per_node, self.thread_count = get_aprun_cmd_for_case(self, "e3sm.exe")
            self.spare_nodes = env_mach_pes.get_spare_nodes(self.num_nodes)
            self.num_nodes += self.spare_nodes
        else:
            self.total_tasks = env_mach_pes.get_total_tasks(comp_classes)
            self.tasks_per_node = env_mach_pes.get_tasks_per_node(self.total_tasks, self.thread_count)

            self.num_nodes, self.spare_nodes = env_mach_pes.get_total_nodes(self.total_tasks, self.thread_count)
            self.num_nodes += self.spare_nodes

        logger.debug("total_tasks {} thread_count {}".format(self.total_tasks, self.thread_count))

        self.tasks_per_numa = int(math.ceil(self.tasks_per_node / 2.0))
        smt_factor = max(1,int(self.get_value("MAX_TASKS_PER_NODE") / max_mpitasks_per_node))

        threads_per_node = self.tasks_per_node * self.thread_count
        threads_per_core = 1 if (threads_per_node <= max_mpitasks_per_node) else smt_factor
        self.cores_per_task = self.thread_count / threads_per_core

        os.environ["OMP_NUM_THREADS"] = str(self.thread_count)

        self.srun_binding = smt_factor*max_mpitasks_per_node / self.tasks_per_node

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

    def read_xml(self):
        self._env_entryid_files = []
        self._env_entryid_files.append(EnvCase(self._caseroot, components=None, read_only=self._force_read_only))
        components = self._env_entryid_files[0].get_values("COMP_CLASSES")
        self._env_entryid_files.append(EnvRun(self._caseroot, components=components, read_only=self._force_read_only))
        self._env_entryid_files.append(EnvBuild(self._caseroot, components=components, read_only=self._force_read_only))
        self._env_entryid_files.append(EnvMachPes(self._caseroot, components=components, read_only=self._force_read_only))
        self._env_entryid_files.append(EnvBatch(self._caseroot, read_only=self._force_read_only))
        if os.path.isfile(os.path.join(self._caseroot,"env_test.xml")):
            self._env_entryid_files.append(EnvTest(self._caseroot, components=components, read_only=self._force_read_only))
        self._env_generic_files = []
        self._env_generic_files.append(EnvMachSpecific(self._caseroot, read_only=self._force_read_only))
        self._env_generic_files.append(EnvArchive(self._caseroot, read_only=self._force_read_only))
        self._files = self._env_entryid_files + self._env_generic_files

    def get_case_root(self):
        """Returns the root directory for this case."""
        return self._caseroot

    def get_env(self, short_name, allow_missing=False):
        full_name = "env_{}.xml".format(short_name)
        for env_file in self._files:
            if os.path.basename(env_file.filename) == full_name:
                return env_file
        if allow_missing:
            return None
        expect(False,"Could not find object for {} in case".format(full_name))

    def copy(self, newcasename, newcaseroot, newcimeroot=None, newsrcroot=None):
        newcase = deepcopy(self)
        for env_file in newcase._files: # pylint: disable=protected-access
            basename = os.path.basename(env_file.filename)
            newfile = os.path.join(newcaseroot, basename)
            env_file.change_file(newfile, copy=True)

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
        for env_file in self._files:
            env_file.write(force_write=flushall)

    def get_values(self, item, attribute=None, resolved=True, subgroup=None):
        for env_file in self._files:
            # Wait and resolve in self rather than in env_file
            results = env_file.get_values(item, attribute, resolved=False, subgroup=subgroup)
            if len(results) > 0:
                new_results = []
                if resolved:
                    for result in results:
                        if isinstance(result, six.string_types):
                            result = self.get_resolved_value(result)
                            vtype = env_file.get_type_info(item)
                            if vtype is not None or vtype != "char":
                                result = convert_to_type(result, vtype, item)

                            new_results.append(result)

                        else:
                            new_results.append(result)

                else:
                    new_results = results

                return new_results

        # Return empty result
        return []

    def get_value(self, item, attribute=None, resolved=True, subgroup=None):
        result = None
        for env_file in self._files:
            # Wait and resolve in self rather than in env_file
            result = env_file.get_value(item, attribute, resolved=False, subgroup=subgroup)

            if result is not None:
                if resolved and isinstance(result, six.string_types):
                    result = self.get_resolved_value(result)
                    vtype = env_file.get_type_info(item)
                    if vtype is not None and vtype != "char":
                        result = convert_to_type(result, vtype, item)

                return result

        # Return empty result
        return result

    def get_record_fields(self, variable, field):
        """ get_record_fields gets individual requested field from an entry_id file
        this routine is used only by xmlquery """
        # Empty result
        result = []

        for env_file in self._env_entryid_files:
            # Wait and resolve in self rather than in env_file
            logger.debug("(get_record_field) Searching in {}".format(env_file.__class__.__name__))
            if field == "varid":
                roots = env_file.scan_children("entry")
            else:
                roots = env_file.get_nodes_by_id(variable)

            for root in roots:
                if root is not None:
                    if field == "raw":
                        result.append(env_file.get_raw_record(root))
                    elif field == "desc":
                        result.append(env_file.get_description(root))
                    elif field == "varid":
                        result.append(env_file.get(root, "id"))
                    elif field == "group":
                        result.extend(env_file.get_groups(root))
                    elif field == "valid_values":
                        # pylint: disable=protected-access
                        vv = env_file._get_valid_values(root)
                        if vv:
                            result.extend(vv)
                    elif field == "file":
                        result.append(env_file.filename)

        if not result:
            for env_file in self._env_generic_files:
                roots = env_file.scan_children(variable)
                for root in roots:
                    if root is not None:
                        if field == "raw":
                            result.append(env_file.get_raw_record(root))
                        elif field == "group":
                            result.extend(env_file.get_groups(root))
                        elif field == "file":
                            result.append(env_file.filename)

        return list(set(result))

    def get_type_info(self, item):
        result = None
        for env_file in self._env_entryid_files:
            result = env_file.get_type_info(item)
            if result is not None:
                return result

        return result

    def get_resolved_value(self, item, recurse=0, allow_unresolved_envvars=False):
        num_unresolved = item.count("$") if item else 0
        recurse_limit = 10
        if (num_unresolved > 0 and recurse < recurse_limit ):
            for env_file in self._env_entryid_files:
                item = env_file.get_resolved_value(item,
                                                   allow_unresolved_envvars=allow_unresolved_envvars)
            if ("$" not in item):
                return item
            else:
                item = self.get_resolved_value(item, recurse=recurse+1,
                                               allow_unresolved_envvars=allow_unresolved_envvars)

        return item

    def set_value(self, item, value, subgroup=None, ignore_type=False, allow_undefined=False, return_file=False):
        """
        If a file has been defined, and the variable is in the file,
        then that value will be set in the file object and the resovled value
        is returned unless return_file is True, in which case (resolved_value, filename)
        is returned where filename is the name of the modified file.
        """
        if item == "CASEROOT":
            self._caseroot = value
        result = None

        for env_file in self._files:
            result = env_file.set_value(item, value, subgroup, ignore_type)
            if (result is not None):
                logger.debug("Will rewrite file {} {}".format(env_file.filename, item))
                return (result, env_file.filename) if return_file else result

        if len(self._files) == 1:
            expect(allow_undefined or result is not None,
                   "No variable {} found in file {}".format(item, self._files[0].filename))
        else:
            expect(allow_undefined or result is not None,
                   "No variable {} found in case".format(item))

    def set_valid_values(self, item, valid_values):
        """
        Update or create a valid_values entry for item and populate it
        """
        result = None
        for env_file in self._env_entryid_files:
            result = env_file.set_valid_values(item, valid_values)
            if (result is not None):
                logger.debug("Will rewrite file {} {}".format(env_file.filename, item))
                return result

    def set_lookup_value(self, item, value):
        if item in self.lookups and self.lookups[item] is not None:
            logger.warning("Item {} already in lookups with value {}".format(item,self.lookups[item]))
        else:
            logger.debug("Setting in lookups: item {}, value {}".format(item,value))
            self.lookups[item] = value

    def clean_up_lookups(self, allow_undefined=False):
        # put anything in the lookups table into existing env objects
        for key,value in list(self.lookups.items()):
            logger.debug("lookup key {} value {}".format(key, value))
            result = self.set_value(key,value, allow_undefined=allow_undefined)
            if result is not None:
                del self.lookups[key]

    def _set_compset(self, compset_name, files, driver="mct"):
        """
        Loop through all the compset files and find the compset
        specifation file that matches either the input 'compset_name'.
        Note that the input compset name (i.e. compset_name) can be
        either a longname or an alias. This will set various compset-related
        info.

        Returns a tuple: (compset_alias, science_support, component_defining_compset)
        (For a user-defined compset - i.e., a compset without an alias - these
        return values will be None, [], None.)
        """
        science_support = []
        compset_alias = None
        components = files.get_components("COMPSETS_SPEC_FILE")
        logger.debug(" Possible components for COMPSETS_SPEC_FILE are {}".format(components))

        self.set_lookup_value("COMP_INTERFACE", driver)
        if self._cime_model == 'cesm':
            comp_root_dir_cpl = files.get_value("COMP_ROOT_DIR_CPL")
            self.set_lookup_value("COMP_ROOT_DIR_CPL",comp_root_dir_cpl)

        # Loop through all of the files listed in COMPSETS_SPEC_FILE and find the file
        # that has a match for either the alias or the longname in that order
        for component in components:

            # Determine the compsets file for this component
            compsets_filename = files.get_value("COMPSETS_SPEC_FILE", {"component":component})

            # If the file exists, read it and see if there is a match for the compset alias or longname
            if (os.path.isfile(compsets_filename)):
                compsets = Compsets(compsets_filename)
                match, compset_alias, science_support = compsets.get_compset_match(name=compset_name)
                if match is not None:
                    self._compsetname = match
                    logger.info("Compset longname is {}".format(match))
                    logger.info("Compset specification file is {}".format(compsets_filename))
                    return compset_alias, science_support

        if compset_alias is None:
            logger.info("Did not find an alias or longname compset match for {} ".format(compset_name))
            self._compsetname = compset_name
            # if this is a valiid compset longname there will be at least 7 components.
            components = self.get_compset_components()
            expect(len(components) > 6, "No compset alias {} found and this does not appear to be a compset longname.".format(compset_name))

        return None, science_support

    def get_primary_component(self):
        if self._primary_component is None:
            self._primary_component = self._find_primary_component()
        return self._primary_component

    def _find_primary_component(self):
        """
        try to glean the primary component based on compset name
        """
        progcomps = {}
        spec = {}
        primary_component = None

        for comp in self._component_classes:

            if comp == "CPL":
                continue
            spec[comp] = self.get_value("COMP_{}".format(comp))
            notprogcomps = ("D{}".format(comp),"X{}".format(comp),"S{}".format(comp))
            if spec[comp].upper() in notprogcomps:
                progcomps[comp] = False
            else:
                progcomps[comp] = True
        expect("ATM" in progcomps and "LND" in progcomps and "OCN" in progcomps and \
               "ICE" in progcomps, " Not finding expected components in {}".format(self._component_classes))
        if progcomps["ATM"] and progcomps["LND"] and progcomps["OCN"] and \
           progcomps["ICE"]:
            primary_component = "allactive"
        elif progcomps["LND"] and progcomps["OCN"] and progcomps["ICE"]:
            # this is a "J" compset
            primary_component = "allactive"
        elif progcomps["ATM"]:
            if "DOCN%SOM" in self._compsetname and progcomps["LND"]:
                # This is an "E" compset
                primary_component = "allactive"
            else:
                # This is an "F" or "Q" compset
                primary_component = spec["ATM"]
        elif progcomps["LND"]:
            # This is an "I" compset
            primary_component = spec["LND"]
        elif progcomps["OCN"]:
            # This is a "C" or "G" compset
            primary_component = spec["OCN"]
        elif progcomps["ICE"]:
            # This is a "D" compset
            primary_component = spec["ICE"]
        elif "GLC" in progcomps and progcomps["GLC"]:
            # This is a "TG" compset
            primary_component = spec["GLC"]
        else:
            # This is "A", "X" or "S"
            primary_component = "drv"

        return primary_component


    def _set_info_from_primary_component(self, files, pesfile=None):
        """
        Sets file and directory paths that depend on the primary component of
        this compset.

        Assumes that self._primary_component has already been set.
        """
        component = self.get_primary_component()

        compset_spec_file = files.get_value("COMPSETS_SPEC_FILE",
                                            {"component":component}, resolved=False)

        self.set_lookup_value("COMPSETS_SPEC_FILE" ,compset_spec_file)
        if pesfile is None:
            self._pesfile = files.get_value("PES_SPEC_FILE", {"component":component})
            pesfile_unresolved = files.get_value("PES_SPEC_FILE", {"component":component}, resolved=False)
            logger.info("Pes     specification file is {}".format(self._pesfile))
        else:
            self._pesfile = pesfile
            pesfile_unresolved = pesfile
        expect(self._pesfile is not None,"No pesfile found for component {}".format(component))

        self.set_lookup_value("PES_SPEC_FILE", pesfile_unresolved)

        tests_filename = files.get_value("TESTS_SPEC_FILE", {"component":component}, resolved=False)
        tests_mods_dir = files.get_value("TESTS_MODS_DIR" , {"component":component}, resolved=False)
        user_mods_dir  = files.get_value("USER_MODS_DIR"  , {"component":component}, resolved=False)
        self.set_lookup_value("TESTS_SPEC_FILE", tests_filename)
        self.set_lookup_value("TESTS_MODS_DIR" , tests_mods_dir)
        self.set_lookup_value("USER_MODS_DIR"  , user_mods_dir)


    def get_compset_components(self):
        #If are doing a create_clone then, self._compsetname is not set yet
        components = []
        compset = self.get_value("COMPSET")
        if compset is None:
            compset = self._compsetname
        expect(compset is not None,
               "compset is not set")
        # the first element is always the date operator - skip it
        elements = compset.split('_')[1:] # pylint: disable=maybe-no-member
        for element in elements:
            # ignore the possible BGC or TEST modifier
            if element.startswith("BGC%") or element.startswith("TEST"):
                continue
            else:
                element_component = element.split('%')[0].lower()
                if "ww" not in element_component:
                    element_component = re.sub(r'[0-9]*',"",element_component)
                components.append(element_component)
        return components


    def __iter__(self):
        for entryid_file in self._env_entryid_files:
            for key, val in entryid_file:
                if isinstance(val, six.string_types) and '$' in val:
                    yield key, self.get_resolved_value(val)
                else:
                    yield key, val

    def set_comp_classes(self, comp_classes):
        self._component_classes = comp_classes
        for env_file in self._env_entryid_files:
            env_file.set_components(comp_classes)

    def _get_component_config_data(self, files, driver=None):
        # attributes used for multi valued defaults
        # attlist is a dictionary used to determine the value element that has the most matches
        attlist = {"compset":self._compsetname, "grid":self._gridname, "cime_model":self._cime_model}

        # Determine list of component classes that this coupler/driver knows how
        # to deal with. This list follows the same order as compset longnames follow.

        # Add the group and elements for the config_files.xml
        for env_file in self._env_entryid_files:
            env_file.add_elements_by_group(files, attlist)

        drv_config_file = files.get_value("CONFIG_CPL_FILE")
        drv_comp = Component(drv_config_file, "CPL")
        for env_file in self._env_entryid_files:
            env_file.add_elements_by_group(drv_comp, attributes=attlist)

        drv_config_file_model_specific = files.get_value("CONFIG_CPL_FILE_MODEL_SPECIFIC")
        drv_comp_model_specific = Component(drv_config_file_model_specific, 'CPL')

        self._component_description["forcing"] = drv_comp_model_specific.get_forcing_description(self._compsetname)
        logger.info("Compset forcing is {}".format(self._component_description["forcing"]))
        self._component_description["CPL"] = drv_comp_model_specific.get_description(self._compsetname)
        if len(self._component_description["CPL"]) > 0:
            logger.info("Com forcing is {}".format(self._component_description["CPL"]))
        for env_file in self._env_entryid_files:
            env_file.add_elements_by_group(drv_comp_model_specific, attributes=attlist)

        self.clean_up_lookups(allow_undefined=True)

        # loop over all elements of both component_classes and components - and get config_component_file for
        # for each component
        self.set_comp_classes(drv_comp.get_valid_model_components())

        if len(self._component_classes) > len(self._components):
            self._components.append('sesp')

        # will need a change here for new cpl components
        root_dir_node_name = 'COMP_ROOT_DIR_CPL'
        comp_root_dir = files.get_value(root_dir_node_name, {"component":driver}, resolved=False)

        if comp_root_dir is not None:
            self.set_value(root_dir_node_name, comp_root_dir)

        for i in range(1,len(self._component_classes)):
            comp_class = self._component_classes[i]
            comp_name  = self._components[i-1]
            root_dir_node_name = 'COMP_ROOT_DIR_' + comp_class
            node_name = 'CONFIG_' + comp_class + '_FILE'
            comp_root_dir = files.get_value(root_dir_node_name, {"component":comp_name}, resolved=False)
            if comp_root_dir is not None:
                self.set_value(root_dir_node_name, comp_root_dir)

            compatt = {"component":comp_name}
            # Add the group and elements for the config_files.xml
            comp_config_file = files.get_value(node_name, compatt, resolved=False)
            expect(comp_config_file is not None,"No component {} found for class {}".format(comp_name, comp_class))
            self.set_value(node_name, comp_config_file)
            comp_config_file =  files.get_value(node_name, compatt)
            expect(comp_config_file is not None and os.path.isfile(comp_config_file),
                   "Config file {} for component {} not found.".format(comp_config_file, comp_name))
            compobj = Component(comp_config_file, comp_class)
            # For files following version 3 schema this also checks the compsetname validity

            self._component_description[comp_class] = compobj.get_description(self._compsetname)
            expect(self._component_description[comp_class] is not None,
                   "No description found in file {} for component {} in comp_class {}".format(comp_config_file, comp_name, comp_class))
            logger.info("{} component is {}".format(comp_class, self._component_description[comp_class]))
            for env_file in self._env_entryid_files:
                env_file.add_elements_by_group(compobj, attributes=attlist)

        self.clean_up_lookups()

    def _setup_mach_pes(self, pecount, multi_driver, ninst, machine_name, mpilib):
        #--------------------------------------------
        # pe layout
        #--------------------------------------------
        mach_pes_obj = None
        # self._pesfile may already be env_mach_pes.xml if so we can just return
        gfile = GenericXML(infile=self._pesfile)
        ftype = gfile.get_id()
        expect(ftype == "env_mach_pes.xml" or ftype == "config_pes", " Do not recognize {} as a valid CIME pes file {}".format(self._pesfile, ftype))
        if ftype == "env_mach_pes.xml":
            new_mach_pes_obj = EnvMachPes(infile=self._pesfile, components=self._component_classes)
            self.update_env(new_mach_pes_obj, "mach_pes")
            return new_mach_pes_obj.get_value("TOTALPES")

        pesobj = Pes(self._pesfile)

        match1 = re.match('(.+)x([0-9]+)', "" if pecount is None else pecount)
        match2 = re.match('([0-9]+)', "" if pecount is None else pecount)

        pes_ntasks = {}
        pes_nthrds = {}
        pes_rootpe = {}
        pes_pstrid = {}
        other      = {}
        comment = None
        force_tasks = None
        force_thrds = None

        if match1:
            opti_tasks = match1.group(1)
            if opti_tasks.isdigit():
                force_tasks = int(opti_tasks)
            else:
                pes_ntasks = pesobj.find_pes_layout(self._gridname, self._compsetname, machine_name,
                                                    pesize_opts=opti_tasks, mpilib=mpilib)[0]
            force_thrds = int(match1.group(2))
        elif match2:
            force_tasks = int(match2.group(1))
            pes_nthrds = pesobj.find_pes_layout(self._gridname, self._compsetname, machine_name, mpilib=mpilib)[1]
        else:
            pes_ntasks, pes_nthrds, pes_rootpe, pes_pstrid, other, comment = pesobj.find_pes_layout(self._gridname, self._compsetname,
                                                                               machine_name, pesize_opts=pecount, mpilib=mpilib)

        if match1 or match2:
            for component_class in self._component_classes:
                if force_tasks is not None:
                    string_ = "NTASKS_" + component_class
                    pes_ntasks[string_] = force_tasks

                if force_thrds is not None:
                    string_ = "NTHRDS_" + component_class
                    pes_nthrds[string_] = force_thrds

                # Always default to zero rootpe if user forced procs and or threads
                string_ = "ROOTPE_" + component_class
                pes_rootpe[string_] = 0

        mach_pes_obj = self.get_env("mach_pes")
        mach_pes_obj.add_comment(comment)

        if other is not None:
            for key, value in other.items():
                self.set_value(key, value)

        totaltasks = []
        for comp_class in self._component_classes:
            ntasks_str = "NTASKS_{}".format(comp_class)
            nthrds_str = "NTHRDS_{}".format(comp_class)
            rootpe_str = "ROOTPE_{}".format(comp_class)
            pstrid_str = "PSTRID_{}".format(comp_class)

            ntasks = pes_ntasks[ntasks_str] if ntasks_str in pes_ntasks else 1
            nthrds = pes_nthrds[nthrds_str] if nthrds_str in pes_nthrds else 1
            rootpe = pes_rootpe[rootpe_str] if rootpe_str in pes_rootpe else 0
            pstrid = pes_pstrid[pstrid_str] if pstrid_str in pes_pstrid else 1

            totaltasks.append( (ntasks + rootpe) * nthrds )

            mach_pes_obj.set_value(ntasks_str, ntasks)
            mach_pes_obj.set_value(nthrds_str, nthrds)
            mach_pes_obj.set_value(rootpe_str, rootpe)
            mach_pes_obj.set_value(pstrid_str, pstrid)

        if multi_driver:
            mach_pes_obj.set_value("MULTI_DRIVER", True)

        # Make sure that every component has been accounted for
        # set, nthrds and ntasks to 1 otherwise. Also set the ninst values here.
        for compclass in self._component_classes:
            key = "NINST_{}".format(compclass)
            if compclass == "CPL":
                continue
            mach_pes_obj.set_value(key, ninst)

            key = "NTASKS_{}".format(compclass)
            if key not in pes_ntasks:
                mach_pes_obj.set_value(key,1)
            key = "NTHRDS_{}".format(compclass)
            if compclass not in pes_nthrds:
                mach_pes_obj.set_value(compclass,1)

    def configure(self, compset_name, grid_name, machine_name=None,
                  project=None, pecount=None, compiler=None, mpilib=None,
                  pesfile=None, gridfile=None,
                  multi_driver=False, ninst=1, test=False,
                  walltime=None, queue=None, output_root=None,
                  run_unsupported=False, answer=None,
                  input_dir=None, driver=None, non_local=False):

        expect(check_name(compset_name, additional_chars='.'), "Invalid compset name {}".format(compset_name))

        #--------------------------------------------
        # compset, pesfile, and compset components
        #--------------------------------------------
        files = Files(comp_interface=driver)

        compset_alias, science_support = self._set_compset(compset_name, files, driver)

        self._components = self.get_compset_components()

        #--------------------------------------------
        # grid
        #--------------------------------------------
        grids = Grids(gridfile)

        gridinfo = grids.get_grid_info(name=grid_name, compset=self._compsetname, driver=driver)

        self._gridname = gridinfo["GRID"]
        for key,value in gridinfo.items():
            logger.debug("Set grid {} {}".format(key,value))
            self.set_lookup_value(key,value)

        #--------------------------------------------
        # component config data
        #--------------------------------------------
        self._get_component_config_data(files, driver=driver)

        # This needs to be called after self.set_comp_classes, which is called
        # from self._get_component_config_data
        self._primary_component = self.get_primary_component()

        self._set_info_from_primary_component(files, pesfile=pesfile)

        self.clean_up_lookups(allow_undefined=True)

        self.get_compset_var_settings()

        self.clean_up_lookups(allow_undefined=True)

        #--------------------------------------------
        # machine
        #--------------------------------------------
        # set machine values in env_xxx files
        machobj = Machines(machine=machine_name)
        probed_machine = machobj.probe_machine_name()
        machine_name = machobj.get_machine_name()
        self.set_value("MACH", machine_name)
        if probed_machine != machine_name and probed_machine is not None:
            logger.warning("WARNING: User-selected machine '{}' does not match probed machine '{}'".format(machine_name, probed_machine))
        else:
            logger.info("Machine is {}".format(machine_name))

        nodenames = machobj.get_node_names()
        nodenames = [x for x in nodenames if
                     '_system' not in x and '_variables' not in x and 'mpirun' not in x and\
                     'COMPILER' not in x and 'MPILIB' not in x]

        for nodename in nodenames:
            value = machobj.get_value(nodename, resolved=False)
            type_str = self.get_type_info(nodename)
            if type_str is not None:
                logger.debug("machine nodname {} value {}".format(nodename, value))
                self.set_value(nodename, convert_to_type(value, type_str, nodename))

        if compiler is None:
            compiler = machobj.get_default_compiler()
        else:
            expect(machobj.is_valid_compiler(compiler),
                   "compiler {} is not supported on machine {}".format(compiler, machine_name))

        self.set_value("COMPILER",compiler)

        if mpilib is None:
            mpilib = machobj.get_default_MPIlib({"compiler":compiler})
        else:
            expect(machobj.is_valid_MPIlib(mpilib, {"compiler":compiler}),
                   "MPIlib {} is not supported on machine {}".format(mpilib, machine_name))
        self.set_value("MPILIB",mpilib)

        machdir = machobj.get_machines_dir()
        self.set_value("MACHDIR", machdir)

        # Create env_mach_specific settings from machine info.
        env_mach_specific_obj = self.get_env("mach_specific")
        env_mach_specific_obj.populate(machobj)

        self._setup_mach_pes(pecount, multi_driver, ninst, machine_name, mpilib)

        if multi_driver and ninst>1:
            logger.info(" Driver/Coupler has %s instances" % ninst)

        #--------------------------------------------
        # archiving system
        #--------------------------------------------
        env_archive = self.get_env("archive")
        infile_node = files.get_child("entry", {"id":"ARCHIVE_SPEC_FILE"})
        infile = files.get_default_value(infile_node)
        infile = self.get_resolved_value(infile)
        logger.debug("archive defaults located in {}".format(infile))
        archive = Archive(infile=infile, files=files)
        archive.setup(env_archive, self._components, files=files)

        self.set_value("COMPSET",self._compsetname)

        self._set_pio_xml()
        logger.info(" Compset is: {} ".format(self._compsetname))
        logger.info(" Grid is: {} ".format(self._gridname ))
        logger.info(" Components in compset are: {} ".format(self._components))

        if not test and not run_unsupported and self._cime_model == "cesm":
            if grid_name in science_support:
                logger.info("\nThis is a CESM scientifically supported compset at this resolution.\n")
            else:
                self._check_testlists(compset_alias, grid_name, files)

        self.set_value("REALUSER", os.environ["USER"])

        # Set project id
        if project is None:
            project = get_project(machobj)
        if project is not None:
            self.set_value("PROJECT", project)
        elif machobj.get_value("PROJECT_REQUIRED"):
            expect(project is not None, "PROJECT_REQUIRED is true but no project found")
        # Get charge_account id if it exists
        charge_account = get_charge_account(machobj, project)
        if charge_account is not None:
            self.set_value("CHARGE_ACCOUNT", charge_account)

        # Resolve the CIME_OUTPUT_ROOT variable, other than this
        # we don't want to resolve variables until we need them
        if output_root is None:
            output_root = self.get_value("CIME_OUTPUT_ROOT")
        self.set_value("CIME_OUTPUT_ROOT", output_root)
        if non_local:
            self.set_value("EXEROOT", os.path.join(output_root, self.get_value("CASE"), "bld"))
            self.set_value("RUNDIR", os.path.join(output_root, self.get_value("CASE"), "run"))
            self.set_value("NONLOCAL", True)

        # Overwriting an existing exeroot or rundir can cause problems
        exeroot = self.get_value("EXEROOT")
        rundir = self.get_value("RUNDIR")
        for wdir in (exeroot, rundir):
            logging.debug("wdir is {}".format(wdir))
            if os.path.exists(wdir):
                expect(not test, "Directory {} already exists, aborting test".format(wdir))
                if answer is None:
                    response = input("\nDirectory {} already exists, (r)eplace, (a)bort, or (u)se existing?".format(wdir))
                else:
                    response = answer

                if response.startswith("r"):
                    shutil.rmtree(wdir)
                else:
                    expect(response.startswith("u"), "Aborting by user request")

        # miscellaneous settings
        if self.get_value("RUN_TYPE") == 'hybrid':
            self.set_value("GET_REFCASE", True)

        # Turn on short term archiving as cesm default setting
        model = get_model()
        self.set_model_version(model)
        if model == "cesm" and not test:
            self.set_value("DOUT_S",True)
            self.set_value("TIMER_LEVEL", 4)

        if test:
            self.set_value("TEST",True)

        self.initialize_derived_attributes()

        #--------------------------------------------
        # batch system (must come after initialize_derived_attributes)
        #--------------------------------------------
        env_batch = self.get_env("batch")

        batch_system_type = machobj.get_value("BATCH_SYSTEM")
        logger.info("Batch_system_type is {}".format(batch_system_type))
        batch = Batch(batch_system=batch_system_type, machine=machine_name)
        bjobs = batch.get_batch_jobs()

        env_batch.set_batch_system(batch, batch_system_type=batch_system_type)
        env_batch.create_job_groups(bjobs, test)

        if walltime:
            self.set_value("USER_REQUESTED_WALLTIME", walltime, subgroup=self.get_primary_job())
        if queue:
            self.set_value("USER_REQUESTED_QUEUE", queue, subgroup=self.get_primary_job())

        env_batch.set_job_defaults(bjobs, self)

        # Make sure that parallel IO is not specified if total_tasks==1
        if self.total_tasks == 1:
            for compclass in self._component_classes:
                key = "PIO_TYPENAME_{}".format(compclass)
                pio_typename = self.get_value(key)
                if pio_typename in ("pnetcdf", "netcdf4p"):
                    self.set_value(key, "netcdf")

        if input_dir is not None:
            self.set_value("DIN_LOC_ROOT", os.path.abspath(input_dir))

    def get_compset_var_settings(self):
        compset_obj = Compsets(infile=self.get_value("COMPSETS_SPEC_FILE"))
        matches = compset_obj.get_compset_var_settings(self._compsetname, self._gridname)
        for name, value in matches:
            if len(value) > 0:
                logger.info("Compset specific settings: name is {} and value is {}".format(name, value))
                self.set_lookup_value(name, value)

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
        machine = self.get_value("MACH")
        toolsdir = os.path.join(self.get_value("CIMEROOT"),"scripts","Tools")
        casetools = os.path.join(self._caseroot, "Tools")
        # setup executable files in caseroot/
        exefiles = (os.path.join(toolsdir, "case.setup"),
                    os.path.join(toolsdir, "case.build"),
                    os.path.join(toolsdir, "case.submit"),
                    os.path.join(toolsdir, "case.qstatus"),
                    os.path.join(toolsdir, "case.cmpgen_namelists"),
                    os.path.join(toolsdir, "preview_namelists"),
                    os.path.join(toolsdir, "preview_run"),
                    os.path.join(toolsdir, "check_input_data"),
                    os.path.join(toolsdir, "check_case"),
                    os.path.join(toolsdir, "xmlchange"),
                    os.path.join(toolsdir, "xmlquery"),
                    os.path.join(toolsdir, "pelayout"))
        try:
            for exefile in exefiles:
                destfile = os.path.join(self._caseroot,os.path.basename(exefile))
                os.symlink(exefile, destfile)
        except Exception as e:
            logger.warning("FAILED to set up exefiles: {}".format(str(e)))

        # set up utility files in caseroot/Tools/
        toolfiles = [os.path.join(toolsdir, "check_lockedfiles"),
                     os.path.join(toolsdir, "getTiming"),
                     os.path.join(toolsdir, "save_provenance"),
                     os.path.join(toolsdir,"Makefile"),
                     os.path.join(toolsdir,"mkSrcfiles"),
                     os.path.join(toolsdir,"mkDepends")]

        # used on Titan
        if os.path.isfile( os.path.join(toolsdir,"mdiag_reduce.csh") ):
            toolfiles.append( os.path.join(toolsdir,"mdiag_reduce.csh") )
            toolfiles.append( os.path.join(toolsdir,"mdiag_reduce.pl") )

        for toolfile in toolfiles:
            destfile = os.path.join(casetools, os.path.basename(toolfile))
            expect(os.path.isfile(toolfile)," File {} does not exist".format(toolfile))
            try:
                os.symlink(toolfile, destfile)
            except Exception as e:
                logger.warning("FAILED to set up toolfiles: {} {} {}".format(str(e), toolfile, destfile))

        if get_model() == "e3sm":
            if os.path.exists(os.path.join(machines_dir, "syslog.{}".format(machine))):
                safe_copy(os.path.join(machines_dir, "syslog.{}".format(machine)), os.path.join(casetools, "mach_syslog"))
            else:
                safe_copy(os.path.join(machines_dir, "syslog.noop"), os.path.join(casetools, "mach_syslog"))

        # add archive_metadata to the CASEROOT but only for CESM
        if get_model() == "cesm":
            try:
                exefile = os.path.join(toolsdir, "archive_metadata")
                destfile = os.path.join(self._caseroot,os.path.basename(exefile))
                os.symlink(exefile, destfile)
            except Exception as e:
                logger.warning("FAILED to set up exefiles: {}".format(str(e)))


    def _create_caseroot_sourcemods(self):
        components = self.get_compset_components()
        components.extend(['share', 'drv'])
        readme_message = """Put source mods for the {component} library in this directory.

WARNING: SourceMods are not kept under version control, and can easily
become out of date if changes are made to the source code on which they
are based. We only recommend using SourceMods for small, short-term
changes that just apply to one or two cases. For larger or longer-term
changes, including gradual, incremental changes towards a final
solution, we highly recommend making changes in the main source tree,
leveraging version control (git or svn).
"""

        for component in components:
            directory = os.path.join(self._caseroot,"SourceMods","src.{}".format(component))
            if not os.path.exists(directory):
                os.makedirs(directory)
                # Besides giving some information on SourceMods, this
                # README file serves one other important purpose: By
                # putting a file inside each SourceMods subdirectory, we
                # prevent aggressive scrubbers from scrubbing these
                # directories due to being empty (which can cause builds
                # to fail).
                readme_file = os.path.join(directory, "README")
                with open(readme_file, "w") as fd:
                    fd.write(readme_message.format(component=component))

        if get_model() == "cesm":
        # Note: this is CESM specific, given that we are referencing cism explitly
            if "cism" in components:
                directory = os.path.join(self._caseroot, "SourceMods", "src.cism", "source_cism")
                if not os.path.exists(directory):
                    os.makedirs(directory)
                    readme_file = os.path.join(directory, "README")
                    str_to_write = """Put source mods for the source_cism library in this subdirectory.
This includes any files from $COMP_ROOT_DIR_GLC/source_cism. Anything
else (e.g., mods to source_glc or drivers) goes in the src.cism
directory, NOT in this subdirectory."""

                    with open(readme_file, "w") as fd:
                        fd.write(str_to_write)

    def create_caseroot(self, clone=False):
        if not os.path.exists(self._caseroot):
            # Make the case directory
            logger.info(" Creating Case directory {}".format(self._caseroot))
            os.makedirs(self._caseroot)
        os.chdir(self._caseroot)

        # Create relevant directories in $self._caseroot
        if clone:
            newdirs = (LOCKED_DIR, "Tools")
        else:
            newdirs = ("SourceMods", LOCKED_DIR, "Buildconf", "Tools")
        for newdir in newdirs:
            os.makedirs(newdir)

        # Open a new README.case file in $self._caseroot
        append_status(" ".join(sys.argv), "README.case", caseroot=self._caseroot)
        compset_info = "Compset longname is {}".format(self.get_value("COMPSET"))
        append_status(compset_info,
                      "README.case", caseroot=self._caseroot)
        append_status("Compset specification file is {}".format(self.get_value("COMPSETS_SPEC_FILE")),
                      "README.case", caseroot=self._caseroot)
        append_status("Pes     specification file is {}".format(self.get_value("PES_SPEC_FILE")),
                      "README.case", caseroot=self._caseroot)
        if "forcing" in self._component_description:
            append_status("Forcing is {}".format(self._component_description["forcing"])
                      ,"README.case", caseroot=self._caseroot)
        for component_class in self._component_classes:
            if component_class in self._component_description and \
               len(self._component_description[component_class])>0:
                append_status("Component {} is {}".format(component_class, self._component_description[component_class]),"README.case", caseroot=self._caseroot)
            if component_class == "CPL":
                append_status("Using %s coupler instances" %
                              (self.get_value("NINST_CPL")),
                              "README.case", caseroot=self._caseroot)
                continue
            comp_grid = "{}_GRID".format(component_class)

            append_status("{} is {}".format(comp_grid,self.get_value(comp_grid)),
                          "README.case", caseroot=self._caseroot)
            comp = str(self.get_value("COMP_{}".format(component_class)))
            user_mods = self._get_comp_user_mods(comp)
            if user_mods is not None:
                note = "This component includes user_mods {}".format(user_mods)
                append_status(note, "README.case", caseroot=self._caseroot)
                logger.info(note)
        if not clone:
            self._create_caseroot_sourcemods()
        self._create_caseroot_tools()

    def apply_user_mods(self, user_mods_dir=None):
        """
        User mods can be specified on the create_newcase command line (usually when called from create test)
        or they can be in the compset definition, or both.
        """
        all_user_mods = []
        for comp in self._component_classes:
            component = str(self.get_value("COMP_{}".format(comp)))
            if component == self._primary_component:
                continue
            comp_user_mods = self._get_comp_user_mods(component)
            if comp_user_mods is not None:
                all_user_mods.append(comp_user_mods)
        # get the primary last so that it takes precidence over other components
        comp_user_mods = self._get_comp_user_mods(self._primary_component)
        if comp_user_mods is not None:
            all_user_mods.append(comp_user_mods)
        if user_mods_dir is not None:
            all_user_mods.append(user_mods_dir)

        # This looping order will lead to the specified user_mods_dir taking
        # precedence over self._user_mods, if there are any conflicts.
        for user_mods in all_user_mods:
            if os.path.isabs(user_mods):
                user_mods_path = user_mods
            else:
                user_mods_path = self.get_value('USER_MODS_DIR')
                user_mods_path = os.path.join(user_mods_path, user_mods)
            apply_user_mods(self._caseroot, user_mods_path)

        # User mods may have modified underlying XML files
        if all_user_mods:
            self.read_xml()

    def _get_comp_user_mods(self, component):
        """
        For a component 'foo', gets the value of FOO_USER_MODS.

        Returns None if no value was found, or if the value is an empty string.
        """
        comp_user_mods = self.get_value("{}_USER_MODS".format(component.upper()))
        #pylint: disable=no-member
        if comp_user_mods is None or comp_user_mods == "" or comp_user_mods.isspace():
            return None
        else:
            return comp_user_mods

    def submit_jobs(self, no_batch=False, job=None, skip_pnl=None, prereq=None, allow_fail=False,
                    resubmit_immediate=False, mail_user=None, mail_type=None, batch_args=None,
                    dry_run=False):
        env_batch = self.get_env('batch')
        result =  env_batch.submit_jobs(self, no_batch=no_batch, skip_pnl=skip_pnl,
                                        job=job, user_prereq=prereq, allow_fail=allow_fail,
                                        resubmit_immediate=resubmit_immediate,
                                        mail_user=mail_user, mail_type=mail_type,
                                        batch_args=batch_args, dry_run=dry_run)
        return result

    def get_job_info(self):
        """
        Get information on batch jobs associated with this case
        """
        xml_job_ids = self.get_value("JOB_IDS")
        if not xml_job_ids:
            return {}
        else:
            result = {}
            job_infos = xml_job_ids.split(", ") # pylint: disable=no-member
            for job_info in job_infos:
                jobname, jobid = job_info.split(":")
                result[jobname] = jobid

            return result

    def report_job_status(self):
        jobmap = self.get_job_info()
        if not jobmap:
            logger.info("No job ids associated with this case. Either case.submit was not run or was run with no-batch")
        else:
            for jobname, jobid in jobmap.items():
                status = self.get_env("batch").get_status(jobid)
                if status:
                    logger.info("{}: {}".format(jobname, status))
                else:
                    logger.info("{}: Unable to get status. Job may be complete already.".format(jobname))

    def cancel_batch_jobs(self, jobids):
        env_batch = self.get_env('batch')
        for jobid in jobids:
            success = env_batch.cancel_job(jobid)
            if not success:
                logger.warning("Failed to kill {}".format(jobid))

    def get_mpirun_cmd(self, job=None, allow_unresolved_envvars=True):
        if job is None:
            job = self.get_primary_job()

        env_mach_specific = self.get_env('mach_specific')
        run_exe = env_mach_specific.get_value("run_exe")
        run_misc_suffix = env_mach_specific.get_value("run_misc_suffix")
        run_misc_suffix = "" if run_misc_suffix is None else run_misc_suffix
        run_suffix = run_exe + run_misc_suffix

        mpirun_cmd_override = self.get_value("MPI_RUN_COMMAND")
        if mpirun_cmd_override not in ["", None, "UNSET"]:
            return self.get_resolved_value(mpirun_cmd_override + " " + run_exe + " " + run_misc_suffix)

        # Things that will have to be matched against mpirun element attributes
        mpi_attribs = {
            "compiler" : self.get_value("COMPILER"),
            "mpilib"   : self.get_value("MPILIB"),
            "threaded" : self.get_build_threaded(),
            "queue" : self.get_value("JOB_QUEUE", subgroup=job),
            "unit_testing" : False
            }

        executable, mpi_arg_list = env_mach_specific.get_mpirun(self, mpi_attribs, job)

        # special case for aprun
        if executable is not None and "aprun" in executable and not "theta" in self.get_value("MACH"):
            aprun_args, num_nodes = get_aprun_cmd_for_case(self, run_exe)[0:2]
            expect( (num_nodes + self.spare_nodes) == self.num_nodes, "Not using optimized num nodes")
            return self.get_resolved_value(executable + aprun_args + " " + run_misc_suffix, allow_unresolved_envvars=allow_unresolved_envvars)

        else:
            mpi_arg_string = " ".join(mpi_arg_list)

        if self.get_value("BATCH_SYSTEM") == "cobalt":
            mpi_arg_string += " : "

        return self.get_resolved_value("{} {} {}".format(executable if executable is not None else "", mpi_arg_string, run_suffix), allow_unresolved_envvars=allow_unresolved_envvars)

    def set_model_version(self, model):
        version = "unknown"
        srcroot = self.get_value("SRCROOT")
        version = get_current_commit(True, srcroot, tag=(model=="cesm"))

        self.set_value("MODEL_VERSION", version)

        if version != "unknown":
            logger.info("{} model version found: {}".format(model, version))
        else:
            logger.warning("WARNING: No {} Model version found.".format(model))

    def load_env(self, reset=False, job=None, verbose=False):
        if not self._is_env_loaded or reset:
            if job is None:
                job = self.get_primary_job()
            os.environ["OMP_NUM_THREADS"] = str(self.thread_count)
            env_module = self.get_env("mach_specific")
            env_module.load_env(self, job=job, verbose=verbose)
            self._is_env_loaded = True

    def get_build_threaded(self):
        """
        Returns True if current settings require a threaded build/run.
        """
        force_threaded = self.get_value("FORCE_BUILD_SMP")
        smp_present = bool(force_threaded) or self.thread_count > 1
        return smp_present

    def _check_testlists(self, compset_alias, grid_name, files):
        """
        CESM only: check the testlist file for tests of this compset grid combination

        compset_alias should be None for a user-defined compset (i.e., a compset
        without an alias)
        """
        if "TESTS_SPEC_FILE" in self.lookups:
            tests_spec_file = self.get_resolved_value(self.lookups["TESTS_SPEC_FILE"])
        else:
            tests_spec_file = self.get_value("TESTS_SPEC_FILE")

        testcnt = 0
        if os.path.isfile(tests_spec_file) and compset_alias is not None:
            # It's important that we not try to find matching tests if
            # compset_alias is None, since compset=None tells get_tests to find
            # tests of all compsets!
            # Only collect supported tests as this _check_testlists is only
            #   called if run_unsupported is False.
            tests = Testlist(tests_spec_file, files)
            testlist = tests.get_tests(compset=compset_alias, grid=grid_name, supported_only=True)
            for test in testlist:
                if test["category"] == "prealpha" or test["category"] == "prebeta" or "aux_" in test["category"]:
                    testcnt += 1
        if testcnt > 0:
            logger.warning("\n*********************************************************************************************************************************")
            logger.warning("This compset and grid combination is not scientifically supported, however it is used in {:d} tests.".format(testcnt))
            logger.warning("*********************************************************************************************************************************\n")
        else:
            expect(False, "\nThis compset and grid combination is untested in CESM.  "
                   "Override this warning with the --run-unsupported option to create_newcase.",
                   error_prefix="STOP: ")

    def set_file(self, xmlfile):
        """
        force the case object to consider only xmlfile
        """
        expect(os.path.isfile(xmlfile), "Could not find file {}".format(xmlfile))

        if not self._read_only_mode:
            self.flush(flushall=True)

        gfile = GenericXML(infile=xmlfile)
        ftype = gfile.get_id()

        logger.warning("setting case file to {}".format(xmlfile))
        components = self.get_value("COMP_CLASSES")
        new_env_file = None
        for env_file in self._env_entryid_files:
            if os.path.basename(env_file.filename) == ftype:
                if ftype == "env_run.xml":
                    new_env_file = EnvRun(infile=xmlfile, components=components)
                elif ftype == "env_build.xml":
                    new_env_file = EnvBuild(infile=xmlfile, components=components)
                elif ftype == "env_case.xml":
                    new_env_file = EnvCase(infile=xmlfile, components=components)
                elif ftype == "env_mach_pes.xml":
                    new_env_file = EnvMachPes(infile=xmlfile, components=components)
                elif ftype == "env_batch.xml":
                    new_env_file = EnvBatch(infile=xmlfile)
                elif ftype == "env_test.xml":
                    new_env_file = EnvTest(infile=xmlfile)
            if new_env_file is not None:
                self._env_entryid_files = []
                self._env_generic_files = []
                self._env_entryid_files.append(new_env_file)
                break
        if new_env_file is None:
            for env_file in self._env_generic_files:
                if os.path.basename(env_file.filename) == ftype:
                    if ftype == "env_archive.xml":
                        new_env_file = EnvArchive(infile=xmlfile)
                    elif ftype == "env_mach_specific.xml":
                        new_env_file = EnvMachSpecific(infile=xmlfile)
                    else:
                        expect(False, "No match found for file type {}".format(ftype))
                if new_env_file is not None:
                    self._env_entryid_files = []
                    self._env_generic_files = []
                    self._env_generic_files.append(new_env_file)
                    break

        self._files = self._env_entryid_files + self._env_generic_files

    def update_env(self, new_object, env_file):
        """
        Replace a case env object file
        """
        old_object = self.get_env(env_file)
        new_object.filename = old_object.filename
        if old_object in self._env_entryid_files:
            self._env_entryid_files.remove(old_object)
            self._env_entryid_files.append(new_object)
        elif old_object in self._env_generic_files:
            self._env_generic_files.remove(old_object)
            self._env_generic_files.append(new_object)
        self._files.remove(old_object)
        self._files.append(new_object)

    def get_latest_cpl_log(self, coupler_log_path=None, cplname="cpl"):
        """
        find and return the latest cpl log file in the
        coupler_log_path directory
        """
        if coupler_log_path is None:
            coupler_log_path = self.get_value("RUNDIR")
        cpllog = None
        cpllogs = glob.glob(os.path.join(coupler_log_path, '{}.log.*'.format(cplname)))
        if cpllogs:
            cpllog = max(cpllogs, key=os.path.getctime)
            return cpllog
        else:
            return None

    def create(self, casename, srcroot, compset_name, grid_name,
               user_mods_dir=None, machine_name=None,
               project=None, pecount=None, compiler=None, mpilib=None,
               pesfile=None, gridfile=None,
               multi_driver=False, ninst=1, test=False,
               walltime=None, queue=None, output_root=None,
               run_unsupported=False, answer=None,
               input_dir=None, driver=None, non_local=False):
        try:
            # Set values for env_case.xml
            self.set_lookup_value("CASE", os.path.basename(casename))
            self.set_lookup_value("CASEROOT", self._caseroot)
            self.set_lookup_value("SRCROOT", srcroot)

            # Configure the Case
            self.configure(compset_name, grid_name, machine_name=machine_name,
                           project=project,
                           pecount=pecount, compiler=compiler, mpilib=mpilib,
                           pesfile=pesfile, gridfile=gridfile,
                           multi_driver=multi_driver, ninst=ninst, test=test,
                           walltime=walltime, queue=queue,
                           output_root=output_root,
                           run_unsupported=run_unsupported, answer=answer,
                           input_dir=input_dir, driver=driver, non_local=non_local)

            self.create_caseroot()

            # Write out the case files
            self.flush(flushall=True)
            self.apply_user_mods(user_mods_dir)

            # Lock env_case.xml
            lock_file("env_case.xml", self._caseroot)
        except:
            if os.path.exists(self._caseroot):
                if not logger.isEnabledFor(logging.DEBUG) and not test:
                    logger.warning("Failed to setup case, removing {}\nUse --debug to force me to keep caseroot".format(self._caseroot))
                    shutil.rmtree(self._caseroot)
                else:
                    logger.warning("Leaving broken case dir {}".format(self._caseroot))

            raise

    def is_save_timing_dir_project(self,project):
        """
        Check whether the project is permitted to archive performance data in the location
        specified for the current machine
        """
        save_timing_dir_projects = self.get_value("SAVE_TIMING_DIR_PROJECTS")
        if not save_timing_dir_projects:
            return False
        else:
            save_timing_dir_projects = save_timing_dir_projects.split(",") # pylint: disable=no-member
            for save_timing_dir_project in save_timing_dir_projects:
                regex = re.compile(save_timing_dir_project)
                if regex.match(project):
                    return True

            return False

    def get_primary_job(self):
        return "case.test" if self.get_value("TEST") else "case.run"
