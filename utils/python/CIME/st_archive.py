"""
functions for performing short term archiving
"""

from XML.standard_module_setup import *
from CIME.case import Case
from CIME.utils import expect, appendStatus
from CIME.XML.env_archive import EnvArchive
from os.path import isfile, isdir, join
import shutil, glob, re

logger = logging.getLogger(__name__)

def check_run(case):
    logger.debug('In check_run...')
    dout_s_root = case.get_value('DOUT_S_ROOT')
    if dout_s_root is None or dout_s_root == 'UNSET':
        expect(False, 'XML variable DOUT_S_ROOT is required for root location of short-term achiver')
    if not isdir(dout_s_root):
        os.makedirs(dout_s_root)
    dout_s_save_interim = case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES')
    if dout_s_save_interim == 'FALSE' or dout_s_save_interim == 'UNSET':
        rest_n = int(case.get_value('REST_N'))
        stop_n = int(case.get_value('STOP_N'))
        if rest_n < stop_n:
            logger.warn('Restart files from end of run will be saved, interim restart files will be deleted')
    statusFile = 'CaseStatus'
    runComplete = False
    caseroot = case.get_value('CASEROOT')
    if isfile(join(caseroot, statusFile)):
        if 'Run SUCCESSFUL' in open(join(caseroot, statusFile)).read():
            runComplete = True
    return runComplete


def list_xml(case, archive):
    logger.debug('In list_xml...')
    for archive_spec_node in archive.get_nodes('comp_archive_spec'):
        comp = archive_spec_node.attrib['name']
        rootdir_node = archive.get_node('rootdir', root=archive_spec_node)
        rootdir = rootdir_node.text
        ninst = case.get_value('NINST_' + rootdir.upper())
        multi = ninst > 1
        logger.info('\n============================================================================')
        logger.info('component name = %s ' % comp)
        logger.info('rootdir = %s' % rootdir)
        logger.info('multiple-instance support = %s ' % multi)
        casename = case.get_value('CASENAME')
        dout_s_root = case.get_value('DOUT_S_ROOT')
        for file_extension in archive.get_nodes('file_extension', root=archive_spec_node):
            suffix = file_extension.attrib['regex_suffix']
            subdir = archive.get_node('subdir', root=file_extension).text
            keep_last_in_rundir = archive.get_node('keep_last_in_rundir', root=file_extension).text
            logger.info('\n  ***** File extension specification')
            logger.info('  regex_suffix = %s ' % suffix)
            logger.info('  subdir = %s ' % join(dout_s_root, rootdir, subdir))
            logger.info('  keep_last_in_rundir %s = ' % keep_last_in_rundir)


def list_archive(dirname):
    logger.debug('In list_archivoe: %s ...' % dirname)
    logger.info('%s ' % dirname)
    if isdir(dirname):
        list_dirs = os.walk(dirname)
        for root, dirs, files in list_dirs:
            for f in files:
                print join(root, f)


def get_datenames(case):
    logger.debug('In get_datename...')
    rundir = case.get_value('RUNDIR')
    expect(isdir(rundir), 'Cannot open directory %s ' % rundir)
    casename = case.get_value("CASE")
    files = sorted(glob.glob(os.path.join(rundir,casename + '.cpl.r*.nc')))
    if not files:
        expect(False, 'Cannot find a %s.cpl.r.*.nc file in directory %s ' % (casename,rundir))
    datenames = []
    for filename in files:
        names = filename.split('.')
        datename = names[-2]
        datenames.append(datename)
        logger.debug('cpl dateName: %s ' % datename)
    return datenames


def archive_move_files(files, case, keep_last_in_rundir, archive_dir):
    logger.debug('In archive_move_files...')
    rundir = case.get_value('RUNDIR')
    if keep_last_in_rundir:
        # move all but the last file in the rundir
        keepfile = files[-1]
        for filename in files:
            if filename != keepfile:
                srcfile = join(rundir, filename)
                destfile = join(archive_dir, filename)
                shutil.move(srcfile, destfile)
        # now simply copy the last file in the rundir
        srcfile = join(rundir, keepfile)
        destfile = join(archive_dir, keepfile)
        print "DEBUG: srcfile %s and destfile %s " %(srcfile,destfile) 
        shutil.copy(srcfile, destfile)
        logger.debug('keepfile = %s ' % keepfile)
    else:
        # move all the files in the rundir
        for filename in files:
            srcfile = join(rundir, filename)
            destfile = join(archive_dir, filename)
            shutil.move(srcfile, destfile)

def get_ninst_info(case, compclass):
        ninst = case.get_value('NINST_' + compclass.upper())
        if ninst is None:
            ninst = 1

        ninst_strings = []
        for i in range(ninst):
            if ninst > 1:
                ninst_strings.append('_' + '%04d' % i)
            else:
                ninst_strings.append('')
        return ninst, ninst_strings
        logger.debug("ninst and ninst_strings are: %s and %s for %s" %(ninst,ninst_strings,compclass))

           
def archive_rpointer(case, archive, archive_spec_node, archive_dir, 
                     datename, datename_is_last):

    casename = case.get_value("CASE")
    compclass = archive.get_node('rootdir', root=archive_spec_node).text
    ninst, ninst_strings = get_ninst_info(case, compclass)

    if (datename_is_last):
        # Copy of all rpointer files for latest restart date
        rundir = case.get_value("RUNDIR")
        rpointers = glob.glob(os.path.join(rundir,'rpointer.*'))
        for rpointer in rpointers:
            shutil.copy(rpointer,os.path.join(archive_dir,os.path.basename(rpointer)))
    else:
        # Generate rpointer file(s) for interim restarts for the one datename and each 
        # possible value of ninst_strings
        if case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES'):
            rpointer_content = archive.get_node('rpointer_content', root=archive_spec_node).text
            if rpointer_content is not 'unset':
                for ninst_string in ninst_strings:
                    subs = dict()
                    subs['$CASE'] = casename
                    subs['$DATENAME'] = datename
                    subs['$NINST_STRING'] = ninst_string
                    for key in subs.keys():
                        rpointer_content = rpointer_content.replace(key, subs[key])

                    rpointer_name = 'rpointer' + ninst_string + '.' + compclass
                    rpointer_file = os.path.join(archive_dir, rpointer_name)
                    logger.debug("writing rpointer_file %s" %rpointer_file)
                    f = open(rpointer_file,'w')
                    for output in rpointer_content.split(','):
                        f.write("%s \n" %output)
                    f.close()


def archive_log_files(case):
    dout_s_root = case.get_value("DOUT_S_ROOT")
    rundir = case.get_value("RUNDIR")
    archive_dir = os.path.join(dout_s_root, 'logs')
    if not os.path.exists(archive_dir):
        os.makedirs(archive_dir)
        logger.debug("created directory archive_dir")
                          
    logfiles = glob.glob(os.path.join(rundir,'*.log.*'))
    for logfile in logfiles:
        srcfile = join(rundir, os.path.basename(logfile))
        destfile = join(archive_dir, os.path.basename(logfile))
        shutil.move(srcfile, destfile)

def archive_history_files(case, archive, archive_spec_node, compclass, compname):
    casename = case.get_value("CASE")
    dout_s_root = case.get_value("DOUT_S_ROOT")
    rundir = case.get_value("RUNDIR")

    # determine ninst and ninst_string
    ninst, ninst_string = get_ninst_info(case, compclass)

    for file_extension in archive.get_nodes('file_extension', root=archive_spec_node):
        subdir = archive.get_node('subdir', root=file_extension).text
        if subdir == 'rest':
            # handle restart files outside of this loop
            continue
        keep_last_in_rundir = archive.get_node('keep_last_in_rundir', root=file_extension).text

        archive_dir = os.path.join(dout_s_root, compclass, subdir)
        if not os.path.exists(archive_dir):
            os.makedirs(archive_dir)
            print "DEBUG: created directory %s" %archive_dir
            logger.debug("created directory %s" %archive_dir)

        suffix = file_extension.attrib['regex_suffix']
        for i in range(ninst):
            newsuffix = compname + ".*" + ninst_string[i] + suffix
            logger.debug("short term archiving suffix is %s " %newsuffix)
            p = re.compile(newsuffix) 
            files = [ f for f in os.listdir(rundir) if p.search(f) ]
            if files:
                logger.debug("hist files are %s " %files)
                print "DEBUG: hist files are %s " %files
                archive_move_files(files, case, keep_last_in_rundir, archive_dir)


def archive_history_files_for_restarts(case, archive, archive_spec_node, archive_dir, restfile):
    restart_hist_varname = archive.get_node("restart_history_varname",root=archive_spec_node).text
    if restart_hist_varname != 'unset':
        rundir = case.get_value("RUNDIR")
        cmd = "ncdump -v %s %s " %(restart_hist_varname, os.path.join(rundir,restfile))
        rc,out,error = run_cmd(cmd, ok_to_fail=True)
        searchname = "%s =" %restart_hist_varname
        if searchname in out:
            offset = out.index(searchname)
            items  = out[offset:].split(",")
            for item in items:
                # the following match has an option of having a './' at the beginning of
                # the history filename
                matchobj = re.search("\"(\.*\/*\w.*)\s?\"",item)
                if matchobj:
                    histfile = matchobj.group(1).strip()
                    if "./" in histfile:
                        histfile.replace("./","")
                    srcfile = os.path.join(rundir,histfile)
                    destfile = os.path.join(archive_dir,histfile)
                    if os.path.isfile(srcfile):
                        shutil.copy(srcfile, destfile)
                    else:
                        logger.info("WARNING: file %s does not exist, will not be archived to restart dir " %srcfile)


def archive_restarts(case, archive, archive_spec_node, compclass, compname, datename, datename_is_last):
    # determine directory for archiving restarts based on datename
    dout_s_root = case.get_value("DOUT_S_ROOT")
    rundir = case.get_value("RUNDIR")
    casename = case.get_value("CASE")


    archive_dir = join(dout_s_root, 'rest', datename)
    if not os.path.exists(archive_dir):
        os.makedirs(archive_dir)
    os.listdir(archive_dir) #DEBUG

    print "DEBUG: datename %s and datename_is_last %s " %(datename,datename_is_last)
    # determine ninst and ninst_string
    ninst, ninst_strings = get_ninst_info(case, compclass)

    # move all but latest restart files into the archive restart directory
    # copy latest restart files to archive restart directory
    # loop over each <file_extension></file_extension> group and find all files that have
    # the regex_suffix
    for file_extension in archive.get_nodes('file_extension', root=archive_spec_node):
        subdir = archive.get_node('subdir', root=file_extension).text
        if subdir != 'rest':
            continue
        regex_suffix = file_extension.attrib['regex_suffix']
        for i in range(ninst):
            rundir = case.get_value("RUNDIR")
            pattern = compname + ninst_strings[i] + regex_suffix + datename 
            p = re.compile(pattern) 
            files = [ f for f in os.listdir(rundir) if p.search(f) ]
            for filename in files:
                restfile = os.path.basename(filename)

                # archive all history files that are needed for restart
                # need to do this before archiving restart files - since
                # are querying the restart files in rundir below
                archive_history_files_for_restarts(case, archive, archive_spec_node, restfile, archive_dir)

                # archive restart files
                # Note that the latest file should be copied and not moved
                srcfile = os.path.join(rundir,filename)
                destfile = os.path.join(archive_dir,restfile)
                if (datename_is_last):
                    shutil.copy(srcfile, destfile)
                else:
                    # Only archive intermediate restarts if requested - otherwise remove them
                    if case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES'):
                        print "DEBUG: moving %s to %s" %(srcfile,destfile)
                        shutil.move(srcfile, destfile)
                    else:
                        os.remove(srcfile)

    # archive the rpointer file(s) for this datename and all possible
    # ninst_strings
    archive_rpointer(case, archive, archive_spec_node, archive_dir, 
                     datename, datename_is_last)

def archive_process(case, archive):
    logger.debug('In archive_process...')
    dout_s_root = case.get_value('DOUT_S_ROOT')
    rundir = case.get_value('RUNDIR')
    compset_comps = case.get_compset_components()
    compset_comps.append('cpl')

    for archive_spec_node in archive.get_nodes('comp_archive_spec'):
        compname = archive_spec_node.attrib['name']
        compname.split('[')[0]
        if compname not in compset_comps:
            continue

        # determine component class
        compclass = archive.get_node('rootdir', root=archive_spec_node).text
        logger.info('doing short term archiving for %s (%s)' % (compname, compclass))

        # archive log files
        archive_log_files(case)

        # archive history files 
        archive_history_files(case, archive, archive_spec_node, compclass, compname)

        # archive restarts and all necessary associated fields (e.g. rpointer files)
        datenames = get_datenames(case)
        for datename in datenames:
            datename_is_last = False
            if datename == datenames[-1]:
                datename_is_last = True
            # archive restarts
            archive_restarts(case, archive, archive_spec_node, 
                             compclass, compname, datename, datename_is_last)


def short_term_archive(input_flag, output_flag, undo_flag):
    case = Case()
    caseroot = case.get_value('CASEROOT')
    archive = EnvArchive(infile=os.path.join(caseroot, 'env_archive.xml'))
    rundir = case.get_value('RUNDIR')
    runComplete = check_run(case)
    if input_flag:
        list_xml(case, archive)

    elif output_flag:
        if runComplete:
            dout_s_root = case.get_value('DOUT_S_ROOT')
            logger.info('Short-term archive listing of %s ' % dout_s_root)
            list_archive(dout_s_root)
        else:
            expect(False, 'st_archive: run is not complete')

    elif undo_flag:
        undoArchive(case)

    elif runComplete:
        # perform short term archiving
        caseroot = case.get_value('CASEROOT')
        appendStatus('st_archiving starting', caseroot=caseroot, sfile='CaseStatus')
        archive_process(case, archive)
        appendStatus('st_archiving completed', caseroot=caseroot, sfile='CaseStatus')
        logger.info('short term archiving is complete.')

    else:
        expect(False, 'based on CaseStatus output, run is either not complete or was not successful')
