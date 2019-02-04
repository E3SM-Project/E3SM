"""
short term archiving
case_st_archive, restore_from_archive, archive_last_restarts
are members of class Case from file case.py
"""

import shutil, glob, re, os

from CIME.XML.standard_module_setup import *
from CIME.utils                     import run_and_log_case_status, ls_sorted_by_mtime, symlink_force, safe_copy, find_files
from CIME.date                      import get_file_date
from CIME.XML.archive       import Archive
from CIME.XML.files            import Files
from os.path                        import isdir, join

logger = logging.getLogger(__name__)

###############################################################################
def _get_archive_file_fn(copy_only):
###############################################################################
    """
    Returns the function to use for archiving some files
    """
    return safe_copy if copy_only else shutil.move

###############################################################################
def _get_datenames(casename, rundir):
###############################################################################
    """
    Returns the date objects specifying the times of each file
    Note we are assuming that the coupler restart files exist and are consistent with other component datenames
    Not doc-testable due to filesystem dependence
    """
    expect(isdir(rundir), 'Cannot open directory {} '.format(rundir))

    files = sorted(glob.glob(os.path.join(rundir, casename +      '.cpl.r.*.nc')))
    if not files:
        files = sorted(glob.glob(os.path.join(rundir, casename + '.cpl_0001.r.*.nc')))

    logger.debug("  cpl files : {} ".format(files))

    if not files:
        logger.warning('Cannot find a {}.cpl*.r.*.nc file in directory {} '.format(casename, rundir))

    datenames = []
    for filename in files:
        file_date = get_file_date(filename)
        datenames.append(file_date)
    return datenames

def _datetime_str(_date):
    """
    Returns the standard format associated with filenames.

    >>> _datetime_str(date(5, 8, 22))
    '0005-08-22-00000'
    >>> _datetime_str(get_file_date("0011-12-09-00435"))
    '0011-12-09-00435'
    """

    format_string = "{year:04d}-{month:02d}-{day:02d}-{seconds:05d}"
    return format_string.format(year = _date.year(),
                                month = _date.month(),
                                day = _date.day(),
                                seconds = _date.second_of_day())

def _datetime_str_mpas(_date):
    """
    Returns the mpas format associated with filenames.

    >>> _datetime_str_mpas(date(5, 8, 22))
    '0005-08-22_00:00:00'
    >>> _datetime_str_mpas(get_file_date("0011-12-09-00435"))
    '0011-12-09_00:07:15'
    """

    format_string = "{year:04d}-{month:02d}-{day:02d}_{hours:02d}:{minutes:02d}:{seconds:02d}"
    return format_string.format(year = _date.year(),
                                month = _date.month(),
                                day = _date.day(),
                                hours = _date.hour(),
                                minutes = _date.minute(),
                                seconds = _date.second())

###############################################################################
def _get_ninst_info(case, compclass):
###############################################################################
    """
    Returns the number of instances used by a component and suffix strings for filenames
    Not doc-testable due to case dependence
    """

    ninst = case.get_value('NINST_' + compclass.upper())
    ninst_strings = []
    if ninst is None:
        ninst = 1
    for i in range(1,ninst+1):
        if ninst > 1:
            ninst_strings.append('_' + '{:04d}'.format(i))

    logger.debug("ninst and ninst_strings are: {} and {} for {}".format(ninst, ninst_strings, compclass))
    return ninst, ninst_strings

###############################################################################
def _get_component_archive_entries(components, archive):
###############################################################################
    """
    Each time this generator function is called, it yields a tuple
    (archive_entry, compname, compclass) for one component in this
    case's compset components.
    """
    for compname in components:
        logger.debug("compname is {} ".format(compname))
        archive_entry = archive.get_entry(compname)
        if archive_entry is None:
            logger.debug("No entry found for {}".format(compname))
            compclass = None
        else:
            compclass = archive.get(archive_entry, "compclass")
        yield(archive_entry, compname, compclass)


###############################################################################
def _archive_rpointer_files(casename, ninst_strings, rundir, save_interim_restart_files, archive,
                            archive_entry, archive_restdir, datename, datename_is_last):
###############################################################################

    if datename_is_last:
        # Copy of all rpointer files for latest restart date
        rpointers = glob.glob(os.path.join(rundir, 'rpointer.*'))
        for rpointer in rpointers:
            safe_copy(rpointer, os.path.join(archive_restdir, os.path.basename(rpointer)))
    else:
        # Generate rpointer file(s) for interim restarts for the one datename and each
        # possible value of ninst_strings
        if save_interim_restart_files:

            # parse env_archive.xml to determine the rpointer files
            # and contents for the given archive_entry tag
            rpointer_items = archive.get_rpointer_contents(archive_entry)

            # loop through the possible rpointer files and contents
            for rpointer_file, rpointer_content in rpointer_items:
                temp_rpointer_file = rpointer_file
                temp_rpointer_content = rpointer_content

                # put in a temporary setting for ninst_strings if they are empty
                # in order to have just one loop over ninst_strings below
                if rpointer_content != 'unset':
                    if not ninst_strings:
                        ninst_strings = ["empty"]

                    for ninst_string in ninst_strings:
                        rpointer_file = temp_rpointer_file
                        rpointer_content = temp_rpointer_content
                        if ninst_string == 'empty':
                            ninst_string = ""
                        for key, value in [('$CASE', casename),
                                           ('$DATENAME', _datetime_str(datename)),
                                           ('$MPAS_DATENAME', _datetime_str_mpas(datename)),
                                           ('$NINST_STRING', ninst_string)]:
                            rpointer_file = rpointer_file.replace(key, value)
                            rpointer_content = rpointer_content.replace(key, value)

                        # write out the respective files with the correct contents
                        rpointer_file = os.path.join(archive_restdir, rpointer_file)
                        logger.info("writing rpointer_file {}".format(rpointer_file))
                        f = open(rpointer_file, 'w')
                        for output in rpointer_content.split(','):
                            f.write("{} \n".format(output))
                        f.close()
                else:
                    logger.info("rpointer_content unset, not creating rpointer file {}".format(rpointer_file))

###############################################################################
def _archive_log_files(dout_s_root, rundir, archive_incomplete, archive_file_fn):
###############################################################################
    """
    Find all completed log files, or all log files if archive_incomplete is True, and archive them.
    Each log file is required to have ".log." in its name, and completed ones will end with ".gz"
    Not doc-testable due to file system dependence
    """
    archive_logdir = os.path.join(dout_s_root, 'logs')
    if not os.path.exists(archive_logdir):
        os.makedirs(archive_logdir)
        logger.debug("created directory {} ".format(archive_logdir))

    if archive_incomplete == False:
        log_search = '*.log.*.gz'
    else:
        log_search = '*.log.*'

    logfiles = glob.glob(os.path.join(rundir, log_search))
    for logfile in logfiles:
        srcfile = join(rundir, os.path.basename(logfile))
        destfile = join(archive_logdir, os.path.basename(logfile))
        archive_file_fn(srcfile, destfile)
        logger.info("moving {} to {}".format(srcfile, destfile))

###############################################################################
def _archive_history_files(archive, archive_entry,
                           compclass, compname, histfiles_savein_rundir,
                           last_date, archive_file_fn, dout_s_root, casename, rundir):
###############################################################################
    """
    perform short term archiving on history files in rundir

    Not doc-testable due to case and file system dependence
    """

    # determine history archive directory (create if it does not exist)

    archive_histdir = os.path.join(dout_s_root, compclass, 'hist')
    if not os.path.exists(archive_histdir):
        os.makedirs(archive_histdir)
        logger.debug("created directory {}".format(archive_histdir))
    # the compname is drv but the files are named cpl
    if compname == 'drv':
        compname = 'cpl'

    if compname == 'clm':
        compname = r'clm2?'

    # determine ninst and ninst_string

    # archive history files - the only history files that kept in the
    # run directory are those that are needed for restarts

    for suffix in archive.get_hist_file_extensions(archive_entry):
        if compname.find('mpas') == 0 or compname == 'mali':
            newsuffix =                    compname + r'\d*'
        else:
            newsuffix = casename + r'\.' + compname + r'_?' + r'\d*'
        newsuffix += r'\.' + suffix
        if not suffix.endswith('$'):
            newsuffix += r'\.'

        logger.debug("short term archiving suffix is {} ".format(newsuffix))

        pfile = re.compile(newsuffix)
        histfiles = [f for f in os.listdir(rundir) if pfile.search(f)]
        logger.debug("histfiles = {} ".format(histfiles))

        if histfiles:
            for histfile in histfiles:
                file_date = get_file_date(os.path.basename(histfile))
                if last_date is None or file_date is None or file_date <= last_date:
                    srcfile = join(rundir, histfile)
                    expect(os.path.isfile(srcfile),
                           "history file {} does not exist ".format(srcfile))
                    destfile = join(archive_histdir, histfile)
                    if histfile in histfiles_savein_rundir:
                        logger.info("copying {} to {} ".format(srcfile, destfile))
                        safe_copy(srcfile, destfile)
                    else:
                        logger.info("moving {} to {} ".format(srcfile, destfile))
                        archive_file_fn(srcfile, destfile)

###############################################################################
def get_histfiles_for_restarts(rundir, archive, archive_entry, restfile, testonly=False):
###############################################################################
    """
    query restart files to determine history files that are needed for restarts

    Not doc-testable due to filesystem dependence
    """

    # Make certain histfiles is a set so we don't repeat
    histfiles = set()
    rest_hist_varname = archive.get_entry_value('rest_history_varname', archive_entry)
    if rest_hist_varname != 'unset':
        cmd = "ncdump -v {} {} ".format(rest_hist_varname, os.path.join(rundir, restfile))
        if testonly:
            out = "{} =".format(rest_hist_varname)
        else:
            rc, out, error = run_cmd(cmd)
            if rc != 0:
                logger.info(" WARNING: {} failed rc={:d}\n    out={}\n    err={}".format(cmd, rc, out, error))
        logger.debug(" get_histfiles_for_restarts: \n    out={}".format(out))

        searchname = "{} =".format(rest_hist_varname)
        if searchname in out:
            offset = out.index(searchname)
            items = out[offset:].split(",")
            for item in items:
                # the following match has an option of having any number of '.'s and '/'s
                # at the beginning of the history filename
                matchobj = re.search(r"\"\S+\s*\"", item)
                if matchobj:
                    histfile = matchobj.group(0).strip('" ')
                    histfile = os.path.basename(histfile)
                    # append histfile to the list ONLY if it exists in rundir before the archiving
                    if histfile in histfiles:
                        logger.warning("WARNING, tried to add a duplicate file to histfiles")
                    if os.path.isfile(os.path.join(rundir,histfile)):
                        histfiles.add(histfile)
                    else:
                        logger.debug(" get_histfiles_for_restarts: histfile {} does not exist ".format(histfile))
    return histfiles

###############################################################################
def _archive_restarts_date(case, casename, rundir, archive,
                           datename, datename_is_last, last_date,
                           archive_restdir, archive_file_fn, components=None,
                           link_to_last_restart_files=False, testonly=False):
###############################################################################
    """
    Archive restart files for a single date

    Returns a dictionary of histfiles that need saving in the run
    directory, indexed by compname
    """
    logger.info('-------------------------------------------')
    logger.info('Archiving restarts for date {}'.format(datename))
    logger.info('-------------------------------------------')
    logger.debug("last date: {}".format(last_date))

    if components is None:
        components = case.get_compset_components()
        components.append('drv')
        components.append('dart')

    histfiles_savein_rundir_by_compname = {}

    for (archive_entry, compname, compclass) in _get_component_archive_entries(components, archive):
        if compclass:
            logger.info('Archiving restarts for {} ({})'.format(compname, compclass))

            # archive restarts
            histfiles_savein_rundir = _archive_restarts_date_comp(case, casename, rundir,
                                                                  archive, archive_entry,
                                                                  compclass, compname,
                                                                  datename, datename_is_last,
                                                                  last_date, archive_restdir,
                                                                  archive_file_fn,
                                                                  link_to_last_restart_files=
                                                                  link_to_last_restart_files,
                                                                  testonly=testonly)
            histfiles_savein_rundir_by_compname[compname] = histfiles_savein_rundir

    return histfiles_savein_rundir_by_compname

###############################################################################
def _archive_restarts_date_comp(case, casename, rundir, archive, archive_entry,
                                compclass, compname, datename, datename_is_last,
                                last_date, archive_restdir, archive_file_fn,
                                link_to_last_restart_files=False, testonly=False):
###############################################################################
    """
    Archive restart files for a single date and single component

    If link_to_last_restart_files is True, then make a symlink to the
    last set of restart files (i.e., the set with datename_is_last
    True); if False (the default), copy them. (This has no effect on the
    history files that are associated with these restart files.)
    """
    datename_str = _datetime_str(datename)

    if datename_is_last or case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES'):
        if not os.path.exists(archive_restdir):
            os.makedirs(archive_restdir)

    # archive the rpointer file(s) for this datename and all possible ninst_strings
    _archive_rpointer_files(casename, _get_ninst_info(case, compclass)[1], rundir,
                            case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES'),
                            archive, archive_entry, archive_restdir, datename, datename_is_last)

    # move all but latest restart files into the archive restart directory
    # copy latest restart files to archive restart directory
    histfiles_savein_rundir = []

    # determine function to use for last set of restart files
    if link_to_last_restart_files:
        last_restart_file_fn = symlink_force
        last_restart_file_fn_msg = "linking"
    else:
        last_restart_file_fn = safe_copy
        last_restart_file_fn_msg = "copying"

    # the compname is drv but the files are named cpl
    if compname == 'drv':
        compname = 'cpl'

    # get file_extension suffixes
    for suffix in archive.get_rest_file_extensions(archive_entry):
#        logger.debug("suffix is {} ninst {}".format(suffix, ninst))
        restfiles = ""
        if compname.find('mpas') == 0 or compname == 'mali':
            pattern = compname + r'\.' + suffix + r'\.' + '_'.join(datename_str.rsplit('-', 1))
            pfile = re.compile(pattern)
            restfiles = [f for f in os.listdir(rundir) if pfile.search(f)]
        else:
            pattern = r"^{}\.{}[\d_]*\.".format(casename, compname)
            pfile = re.compile(pattern)
            files = [f for f in os.listdir(rundir) if pfile.search(f)]
            pattern =  r'_?' + r'\d*' + r'\.' + suffix + r'\.' + r'[^\.]*' + r'\.?' + datename_str
            pfile = re.compile(pattern)
            restfiles = [f for f in files if pfile.search(f)]
            logger.debug("pattern is {} restfiles {}".format(pattern, restfiles))
        for restfile in restfiles:
            restfile = os.path.basename(restfile)

            file_date = get_file_date(restfile)
            if last_date is not None and file_date > last_date:
                # Skip this file
                continue

            if not os.path.exists(archive_restdir):
                os.makedirs(archive_restdir)

            # obtain array of history files for restarts
            # need to do this before archiving restart files
            histfiles_for_restart = get_histfiles_for_restarts(rundir, archive,
                                                               archive_entry, restfile,
                                                               testonly=testonly)

            if datename_is_last and histfiles_for_restart:
                for histfile in histfiles_for_restart:
                    if histfile not in histfiles_savein_rundir:
                        histfiles_savein_rundir.append(histfile)

            # archive restart files and all history files that are needed for restart
            # Note that the latest file should be copied and not moved
            if datename_is_last:
                srcfile = os.path.join(rundir, restfile)
                destfile = os.path.join(archive_restdir, restfile)
                last_restart_file_fn(srcfile, destfile)
                logger.info("{} file {} to {}".format(last_restart_file_fn_msg, srcfile, destfile))
                for histfile in histfiles_for_restart:
                    srcfile = os.path.join(rundir, histfile)
                    destfile = os.path.join(archive_restdir, histfile)
                    expect(os.path.isfile(srcfile),
                           "history restart file {} for last date does not exist ".format(srcfile))
                    logger.info("Copying {} to {}".format(srcfile, destfile))
                    safe_copy(srcfile, destfile)
                    logger.debug("datename_is_last + histfiles_for_restart copying \n  {} to \n  {}".format(srcfile, destfile))
            else:
                # Only archive intermediate restarts if requested - otherwise remove them
                if case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES'):
                    srcfile = os.path.join(rundir, restfile)
                    destfile = os.path.join(archive_restdir, restfile)
                    expect(os.path.isfile(srcfile),
                           "restart file {} does not exist ".format(srcfile))
                    archive_file_fn(srcfile, destfile)
                    logger.info("moving file {} to {}".format(srcfile, destfile))

                    # need to copy the history files needed for interim restarts - since
                    # have not archived all of the history files yet
                    for histfile in histfiles_for_restart:
                        srcfile = os.path.join(rundir, histfile)
                        destfile = os.path.join(archive_restdir, histfile)
                        expect(os.path.isfile(srcfile),
                               "hist file {} does not exist ".format(srcfile))
                        logger.info("copying {} to {}".format(srcfile, destfile))
                        safe_copy(srcfile, destfile)
                else:
                    srcfile = os.path.join(rundir, restfile)
                    logger.info("removing interim restart file {}".format(srcfile))
                    if (os.path.isfile(srcfile)):
                        try:
                            os.remove(srcfile)
                        except OSError:
                            logger.warning("unable to remove interim restart file {}".format(srcfile))
                    else:
                        logger.warning("interim restart file {} does not exist".format(srcfile))

    return histfiles_savein_rundir

###############################################################################
def _archive_process(case, archive, last_date, archive_incomplete_logs, copy_only,
                     components=None,dout_s_root=None, casename=None, rundir=None, testonly=False):
###############################################################################
    """
    Parse config_archive.xml and perform short term archiving
    """

    logger.debug('In archive_process...')

    if dout_s_root is None:
        dout_s_root = case.get_value("DOUT_S_ROOT")
    if rundir is None:
        rundir = case.get_value("RUNDIR")
    if casename is None:
        casename = case.get_value("CASE")
    if components is None:
        components = case.get_compset_components()
        components.append('drv')
        components.append('dart')

    archive_file_fn = _get_archive_file_fn(copy_only)

    # archive log files
    _archive_log_files(dout_s_root, rundir,
                       archive_incomplete_logs, archive_file_fn)

    # archive restarts and all necessary associated files (e.g. rpointer files)
    datenames = _get_datenames(casename, rundir)
    logger.debug("datenames {} ".format(datenames))
    histfiles_savein_rundir_by_compname = {}
    for datename in datenames:
        datename_is_last = False
        if datename == datenames[-1]:
            datename_is_last = True

        logger.debug("datename {} last_date {}".format(datename,last_date))
        if last_date is None or datename <= last_date:
            archive_restdir = join(dout_s_root, 'rest', _datetime_str(datename))

            histfiles_savein_rundir_by_compname_this_date = _archive_restarts_date(
                case, casename, rundir, archive, datename, datename_is_last,
                last_date, archive_restdir, archive_file_fn, components, testonly=testonly)
            if datename_is_last:
                histfiles_savein_rundir_by_compname = histfiles_savein_rundir_by_compname_this_date

    # archive history files

    for (archive_entry, compname, compclass) in _get_component_archive_entries(components, archive):
        if compclass:
            logger.info('Archiving history files for {} ({})'.format(compname, compclass))
            histfiles_savein_rundir = histfiles_savein_rundir_by_compname.get(compname, [])
            logger.debug("_archive_process: histfiles_savein_rundir {} ".format(histfiles_savein_rundir))
            _archive_history_files(archive, archive_entry,
                                   compclass, compname, histfiles_savein_rundir,
                                   last_date, archive_file_fn,
                                   dout_s_root, casename, rundir)

###############################################################################
def restore_from_archive(self, rest_dir=None, dout_s_root=None, rundir=None):
###############################################################################
    """
    Take archived restart files and load them into current case.  Use rest_dir if provided otherwise use most recent
    restore_from_archive is a member of Class Case
    """
    if dout_s_root is None:
        dout_s_root = self.get_value("DOUT_S_ROOT")
    if rundir is None:
        rundir = self.get_value("RUNDIR")
    if rest_dir is not None:
        if not os.path.isabs(rest_dir):
            rest_dir = os.path.join(dout_s_root, "rest", rest_dir)
    else:
        rest_dir = os.path.join(dout_s_root, "rest", ls_sorted_by_mtime(os.path.join(dout_s_root, "rest"))[-1])

    logger.info("Restoring restart from {}".format(rest_dir))

    for item in glob.glob("{}/*".format(rest_dir)):
        base = os.path.basename(item)
        dst = os.path.join(rundir, base)
        if os.path.exists(dst):
            os.remove(dst)
        logger.info("Restoring {} from {} to {}".format(item, rest_dir, rundir))

        safe_copy(item, rundir)

###############################################################################
def archive_last_restarts(self, archive_restdir, rundir, last_date=None, link_to_restart_files=False):
###############################################################################
    """
    Convenience function for archiving just the last set of restart
    files to a given directory. This also saves files attached to the
    restart set, such as rpointer files and necessary history
    files. However, it does not save other files that are typically
    archived (e.g., history files, log files).

    Files are copied to the directory given by archive_restdir.

    If link_to_restart_files is True, then symlinks rather than copies
    are done for the restart files. (This has no effect on the history
    files that are associated with these restart files.)
    """
    archive = self.get_env('archive')
    casename = self.get_value("CASE")
    datenames = _get_datenames(casename, rundir)
    expect(len(datenames) >= 1, "No restart dates found")
    last_datename = datenames[-1]

    # Not currently used for anything if we're only archiving the last
    # set of restart files, but needed to satisfy the following interface
    archive_file_fn = _get_archive_file_fn(copy_only=False)

    _ = _archive_restarts_date(case=self,
                               casename=casename,
                               rundir=rundir,
                               archive=archive,
                               datename=last_datename,
                               datename_is_last=True,
                               last_date=last_date,
                               archive_restdir=archive_restdir,
                               archive_file_fn=archive_file_fn,
                               link_to_last_restart_files=link_to_restart_files)

###############################################################################
def case_st_archive(self, last_date_str=None, archive_incomplete_logs=True, copy_only=False, resubmit=True):
###############################################################################
    """
    Create archive object and perform short term archiving
    """
    caseroot = self.get_value("CASEROOT")
    self.load_env(job="case.st_archive")
    if last_date_str is not None:
        try:
            last_date = get_file_date(last_date_str)
        except ValueError:
            expect(False, 'Could not parse the last date to archive')
    else:
        last_date = None

    dout_s_root = self.get_value('DOUT_S_ROOT')
    if dout_s_root is None or dout_s_root == 'UNSET':
        expect(False,
               'XML variable DOUT_S_ROOT is required for short-term achiver')
    if not isdir(dout_s_root):
        os.makedirs(dout_s_root)

    dout_s_save_interim = self.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES')
    if dout_s_save_interim == 'FALSE' or dout_s_save_interim == 'UNSET':
        rest_n = self.get_value('REST_N')
        stop_n = self.get_value('STOP_N')
        if rest_n < stop_n:
            logger.warning('Restart files from end of run will be saved'
                        'interim restart files will be deleted')

    logger.info("st_archive starting")

    archive = self.get_env('archive')
    functor = lambda: _archive_process(self, archive, last_date, archive_incomplete_logs, copy_only)
    run_and_log_case_status(functor, "st_archive", caseroot=caseroot)

    logger.info("st_archive completed")

    # resubmit case if appropriate
    resubmit_cnt = self.get_value("RESUBMIT")
    logger.debug("resubmit_cnt {} resubmit {}".format(resubmit_cnt, resubmit))
    if resubmit_cnt > 0 and resubmit:
        logger.info("resubmitting from st_archive, resubmit={:d}".format(resubmit_cnt))
        if self.get_value("MACH") == "mira":
            expect(os.path.isfile(".original_host"), "ERROR alcf host file not found")
            with open(".original_host", "r") as fd:
                sshhost = fd.read()
            run_cmd("ssh cooleylogin1 ssh {} '{case}/case.submit {case} --resubmit' "\
                        .format(sshhost, case=caseroot), verbose=True)
        else:
            self.submit(resubmit=True)

    return True

def test_st_archive(self, testdir="st_archive_test"):
    archive = Archive()
    files = Files()
    components = []
#    expect(not self.get_value("MULTI_DRIVER"),"Test not configured for multi-driver cases")

    config_archive_files = archive.get_all_config_archive_files(files)
    # create the run directory testdir and populate it with rest_file and hist_file from
    # config_archive.xml test_file_names
    if os.path.exists(testdir):
        logger.info("Removing existing test directory {}".format(testdir))
        shutil.rmtree(testdir)
    dout_s_root=os.path.join(testdir,"archive")
    archive = Archive()
    schema = files.get_schema("ARCHIVE_SPEC_FILE")
    for config_archive_file in config_archive_files:
        archive.read(config_archive_file, schema)
    comp_archive_specs = archive.get_children("comp_archive_spec")
    for comp_archive_spec in comp_archive_specs:
        components.append(archive.get(comp_archive_spec, 'compname'))
        test_file_names = archive.get_optional_child("test_file_names", root=comp_archive_spec)
        if test_file_names is not None:
            if not os.path.exists(testdir):
                os.makedirs(os.path.join(testdir,"archive"))

            for file_node in archive.get_children("tfile", root=test_file_names):
                fname = os.path.join(testdir,archive.text(file_node))
                disposition = archive.get(file_node, "disposition")
                logger.info("Create file {} with disposition {}".
                            format(fname, disposition))
                with open(fname, 'w') as fd:
                    fd.write(disposition+"\n")

    logger.info("testing components: {} ".format(list(set(components))))
    _archive_process(self, archive, None, False, False,components=list(set(components)),
                     dout_s_root=dout_s_root,
                     casename="casename", rundir=testdir, testonly=True)

    _check_disposition(testdir)

    # Now test the restore capability
    testdir2 = os.path.join(testdir,"run2")
    os.makedirs(testdir2)

    restore_from_archive(self, rundir=testdir2, dout_s_root=dout_s_root)

    restfiles = [f for f in os.listdir(os.path.join(testdir,"archive","rest","1976-01-01-00000"))]
    for _file in restfiles:
        expect(os.path.isfile(os.path.join(testdir2,_file)), "Expected file {} to be restored from rest dir".format(_file))

    return True

def test_env_archive(self, testdir="env_archive_test"):
    components = self.get_values("COMP_CLASSES")
    comps_in_case = []
    # create the run directory testdir and populate it with rest_file and hist_file from
    # config_archive.xml test_file_names
    if os.path.exists(testdir):
        logger.info("Removing existing test directory {}".format(testdir))
        shutil.rmtree(testdir)
    dout_s_root=os.path.join(testdir,"archive")
    archive = self.get_env('archive')
    comp_archive_specs = archive.scan_children("comp_archive_spec")

    # ignore stub and dead components
    for comp in list(components):
        compname = self.get_value("COMP_{}".format(comp))
        if (compname == 's'+comp.lower() or compname == 'x'+comp.lower()) and comp != 'ESP':
            logger.info("Not testing component {}".format(comp))
            components.remove(comp)
        elif comp == 'ESP' and self.get_value('MODEL') == 'e3sm':
            components.remove(comp)
        else:
            if compname == 'cpl':
                compname = 'drv'
            comps_in_case.append(compname)

    for comp_archive_spec in comp_archive_specs:
        comp_expected = archive.get(comp_archive_spec, 'compname')
        if comp_expected == "ww3":
            comp_expected = "ww"
        comp_class = archive.get(comp_archive_spec, 'compclass').upper()
        if comp_class in components:
            components.remove(comp_class)
        else:
            expect(False,"Error finding comp_class {} in components".format(comp_class))
        if comp_expected == 'cpl':
            comp_expected = 'drv'
        if comp_expected != 'dart':
            expect(comp_expected in comps_in_case, "env_archive defines component {} not defined in case".format(comp_expected))

        test_file_names = archive.get_optional_child("test_file_names", root=comp_archive_spec)
        if test_file_names is not None:
            if not os.path.exists(testdir):
                os.makedirs(os.path.join(testdir,"archive"))

            for file_node in archive.get_children("tfile", root=test_file_names):
                fname = os.path.join(testdir,archive.text(file_node))
                disposition = archive.get(file_node, "disposition")
                logger.info("Create file {} with disposition {}".
                            format(fname, disposition))
                with open(fname, 'w') as fd:
                    fd.write(disposition+"\n")

    expect(not components, "No archive entry found for components: {}".format(components))
    if 'dart' not in comps_in_case:
        comps_in_case.append('dart')
    logger.info("testing components: {} ".format(comps_in_case))
    _archive_process(self, archive, None, False, False,components=comps_in_case,
                     dout_s_root=dout_s_root,
                     casename="casename", rundir=testdir, testonly=True)

    _check_disposition(testdir)

    # Now test the restore capability
    testdir2 = os.path.join(testdir,"run2")
    os.makedirs(testdir2)

    restore_from_archive(self, rundir=testdir2, dout_s_root=dout_s_root)

    restfiles = [f for f in os.listdir(os.path.join(testdir,"archive","rest","1976-01-01-00000"))]
    for _file in restfiles:
        expect(os.path.isfile(os.path.join(testdir2,_file)), "Expected file {} to be restored from rest dir".format(_file))

    return True

def _check_disposition(testdir):
    copyfilelist = []
    for root, _, files in os.walk(testdir):
        for _file in files:
            with open(os.path.join(root, _file), "r") as fd:
                disposition = fd.readline()
            logger.info("Checking testfile {} with disposition {}".format(_file, disposition))
            if root == testdir:
                if "move" in disposition:
                    if find_files(os.path.join(testdir, "archive"), _file):
                        expect(False,
                               "Copied file {} to archive with disposition move".format(_file))
                    else:
                        expect(False,
                               "Failed to move file {} to archive".format(_file))
                if "copy" in disposition:
                    copyfilelist.append(_file)
            elif "ignore" in disposition:
                expect(False, "Moved file {} with dispostion ignore to directory {}".format(_file, root))
            elif "copy" in disposition:
                expect(_file in copyfilelist, "File {} with disposition copy was moved to directory {}"
                       .format(_file, root))
    for _file in copyfilelist:
        expect(find_files(os.path.join(testdir,"archive"), _file) != [],
               "File {} was not copied to archive.".format(_file))
