"""
functions for performing short term archiving
"""

from XML.standard_module_setup import *
from CIME.case import Case
from CIME.utils import expect, appendStatus
from CIME.XML.env_archive import EnvArchive
from os.path import isfile, isdir, join, basename
import shutil, glob, re

logger = logging.getLogger(__name__)

def checkRun(case):
    logger.debug('In checkRun...')
    dout_s_root = case.get_value('DOUT_S_ROOT')
    if dout_s_root is None or dout_s_root == 'UNSET':
        expect(False, 'XML variable DOUT_S_ROOT is required for root location of short-term achiver')
    if not isdir(dout_s_root):
        print "DEBUG: dout_s_root is ",dout_s_root
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


def listXMLin(case, archive):
    logger.debug('In listXMLin...')
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


def listArchive(dirname):
    logger.debug('In listArchive: %s ...' % dirname)
    logger.info('%s ' % dirname)
    if isdir(dirname):
        list_dirs = os.walk(dirname)
        for root, dirs, files in list_dirs:
            for f in files:
                print join(root, f)


def getDatename(case):
    logger.debug('In getDatename...')
    rundir = case.get_value('RUNDIR')
    expect(isdir(rundir), 'Cannot open directory %s ' % rundir)
    files = sorted(glob.glob(rundir + '/*cpl.r*.nc'))
    rfile = files[-1]
    names = rfile.split('.')
    if names:
        datename = names[-2]
        logger.debug('Cpl dateName: %s ' % datename)
        return datename
    expect(False, 'Cannot find a cpl.r.*.nc file in directory %s ' % rundir)


def moveFiles(case, keep_last_in_rundir, suffix, destdir, restdir, datename, runfiles):
    logger.debug('In moveFiles...')
    rundir = case.get_value('RUNDIR')
    p = re.compile(suffix)
    files = [ f for f in os.listdir(rundir) if p.match(f) ]
    numfiles = len(files)
    if numfiles > 0:
        if keep_last_in_rundir:
            keepfile = files[-1]
            shutil.copy(join(rundir, keepfile), join(destdir, keepfile))
            logger.debug('keepfile = %s ' % keepfile)
            for filename in files:
                if filename != keepfile and isfile(filename):
                    shutil.move(join(rundir, filename), join(destdir, filename))

        else:
            for filename in files:
                shutil.move(join(rundir, filename), join(destdir, filename))

    for filename in files:
        runfiles.remove(filename)

    counter = 0
    for filename in files:
        if datename in filename:
            counter = counter + 1
            if counter > 1:
                expect(False, 'Multiple restart files found for suffix %s ' % suffix)
            restfile = basename(filename)
            if isfile(join(destdir, restfile)):
                shutil.copy(join(destdir, restfile), join(restdir, restfile))


def archiveProcess(case, archive, datename, runfiles):
    logger.debug('In archiveProcess...')
    dout_s_root = case.get_value('DOUT_S_ROOT')
    rundir = case.get_value('RUNDIR')
    restdir = join(dout_s_root, 'rest', datename)
    if not os.path.exists(restdir):
        os.makedirs(restdir)
    items_to_copy = [ item for item in glob.glob(rundir + '/rpointer.*') ]
    for item in items_to_copy:
        if isfile(item):
            shutil.copy(item, join(restdir, basename(item)))

    for archive_spec_node in archive.get_nodes('comp_archive_spec'):
        compset_comps = case.get_compset_components()
        compset_comps.append('cpl')
        comp = archive_spec_node.attrib['name']
        comp.split('[')[0]
        if comp not in compset_comps:
            continue
        rootdir_node = archive.get_node('rootdir', root=archive_spec_node)
        rootdir = rootdir_node.text
        logger.info('doing short term archiving for %s (%s)' % (comp, rootdir))
        ninst = case.get_value('NINST_' + rootdir.upper())
        if ninst is None and comp == 'cpl':
            ninst = 1
        for file_extension in archive.get_nodes('file_extension', root=archive_spec_node):
            suffix = file_extension.attrib['regex_suffix']
            subdir = archive.get_node('subdir', root=file_extension).text
            keep_last_in_rundir = archive.get_node('keep_last_in_rundir', root=file_extension).text
            destdir = join(dout_s_root, rootdir, subdir)
            if not os.path.exists(destdir):
                os.makedirs(destdir)
            for i in range(ninst):
                ninst_suffix = ''
                if ninst > 1:
                    ninst_suffix = '_' + '%04d' % i
                print 'DEBUG: comp %s with subdir %s' % (comp, subdir)
                newSuffix = suffix
                if suffix[0:1] == '.':
                    if subdir == 'logs':
                        newSuffix = rootdir + ninst_suffix + suffix
                        print 'DEBUG: comp %s log newSuffix %s ' % (comp, newSuffix)
                    else:
                        casename = case.get_value('CASE')
                        newSuffix = casename + '.' + comp + ninst_suffix + suffix
                moveFiles(case, keep_last_in_rundir.lower, newSuffix, destdir, restdir, datename, runfiles)


def short_term_archive(input_flag, output_flag, undo_flag):
    case = Case()
    caseroot = case.get_value('CASEROOT')
    archive = EnvArchive(infile=os.path.join(caseroot, 'env_archive.xml'))
    rundir = case.get_value('RUNDIR')
    runComplete = checkRun(case)
    if input_flag:
        listXMLin(case, archive)
    elif output_flag:
        if runComplete:
            dout_s_root = case.get_value('DOUT_S_ROOT')
            logger.info('Short-term archive listing of %s ' % dout_s_root)
            listArchive(dout_s_root)
        else:
            expect(False, 'st_archive: run is not complete')
    elif undo_flag:
        undoArchive(case)
    elif runComplete:
        datename = getDatename(case)
        runfiles = [ f for f in os.listdir(rundir) if isfile(join(rundir, f)) ]
        caseroot = case.get_value('CASEROOT')
        appendStatus('st_archiving starting', caseroot=caseroot, sfile='CaseStatus')
        archiveProcess(case, archive, datename, runfiles)
        appendStatus('st_archiving completed', caseroot=caseroot, sfile='CaseStatus')
        dout_save_interim = case.get_value('DOUT_S_SAVE_INTERIM_RESTART_FILES')
        logger.info('short term archiving is complete.')
    else:
        expect(False, 'based on CaseStatus output, run is either not complete or was not successful')
