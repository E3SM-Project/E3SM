#!/usr/bin/env python
"""
This script creates coupling files needed to run an MPAS-Ocean and
MPAS-seaice within E3SM.

Load the lastest e3sm-unified conda package.
"""
# import modules # {{{
import os
import shutil
import subprocess
import configparser
import argparse
import numpy as np
import glob
from datetime import datetime
import traceback
import sys
# }}}


def main():
# {{{

    print("****** Creating E3SM coupling files ******")
    # obtain configuration settings
    config = configparser.ConfigParser()
    config.read("config_E3SM_coupling_files.ini")

    function_list = [initial_condition_ocean,
                     graph_partition_ocean,
                     initial_condition_seaice,
                     scrip,
                     transects_and_regions,
                     mapping_CORE_Gcase,
                     mapping_JRA_Gcase,
                     mapping_ne30,
                     domain_CORE_Gcase,
                     domain_JRA_Gcase,
                     domain_ne30,
                     mapping_runoff,
                     salinity_restoring]

    # clean: Delete all directories
    parser = argparse.ArgumentParser()
    parser.add_argument('--clean', action='store_true')
    parser.add_argument('--ice_shelf_cavities', action='store_true')
    if parser.parse_args().clean:
        print('****** clean out directories ******')
        for function in function_list:
            function_name = function.__name__
            try:
                shutil.rmtree(function_name)
                print('removed directory: ', function_name)
            except:
                print('failed to remove directory: ', function_name)
        return

    if parser.parse_args().ice_shelf_cavities:
        ice_shelf_cavities = True
    else:
        ice_shelf_cavities = False
    print("- ice_shelf_cavities set to", ice_shelf_cavities)

    # determine mesh name
    mesh_name = config.get('main', 'mesh_name')
    currentDir = os.getcwd()
    if mesh_name == 'autodetect':
        path = currentDir.split('/')
        try:
            index = path.index('global_ocean') + 1
            mesh_name = 'o' + path[index]
            print("- mesh name autodetected from path: " + mesh_name)
        except BaseException:
            print("mesh name not found in path. Please specify " +
                  "the mesh_name in config_E3SM_coupling_files.ini. Stopping.")
            return
    else:
        print("- mesh name specified in config file: ", mesh_name)

    # determine date string
    date_string = config.get('main', 'date_string')
    currentDir = os.getcwd()
    if date_string == 'autodetect':
        now = datetime.now()
        date_string = now.strftime("%y%m%d")
        print("- date string autodetected from today's date:", date_string)
    else:
        print("- date string specified in config file:", date_string)

    # create inputdata directories
    make_dir('assembled_files_for_upload/inputdata/ocn/mpas-o/' + mesh_name)
    make_dir('assembled_files_for_upload/inputdata/ice/mpas-cice/' + mesh_name)
    make_dir('assembled_files_for_upload/inputdata/cpl/cpl6')
    make_dir('assembled_files_for_upload/inputdata/share/domains')

    success = True
    print()
    for function in function_list:
        function_name = function.__name__
        print("****** " + function_name + " ******")

        if config.get(function_name, 'enable').lower() == 'false':
            print("Disabled in .ini file")
        else:
            make_dir(function_name)
            os.chdir(function_name)

            try:
                function(config, mesh_name, date_string, ice_shelf_cavities)
                print('SUCCESS')
            except BaseException:
                print('!!! FAILURE !!!')
                traceback.print_exc(file=sys.stdout)
                success = False
            os.chdir(currentDir)
        print(" ")

    if success:
        print("****** SUCCESS for all enabled steps ******")
    else:
        print("!!!!!! FAILURE: One or more steps failed. See output above !!!!!!")
# }}}


def initial_condition_ocean(config, mesh_name, date_string, ice_shelf_cavities):
# {{{

    # create links
    make_link('../init.nc', mesh_name + '.nc')

    # command line execution
    args = ['ncks', '-x', '-v', 'xtime', '-O',
            mesh_name + '.nc',
            mesh_name + '_no_xtime.nc'
            ]
    run_command(args)

    # create link to output directory
    os.chdir('../assembled_files_for_upload/inputdata/ocn/mpas-o/' + mesh_name)
    make_link(
        '../../../../../initial_condition_ocean/' +
        mesh_name + '_no_xtime.nc',
        mesh_name + '_no_xtime.nc')
# }}}

def graph_partition_ocean(config, mesh_name, date_string, ice_shelf_cavities):
# {{{

    # create links
    make_link('../graph.info', 'mpas-o.graph.info.' + date_string)

    # command line execution
    nCells = sum(1 for l in open('../graph.info'))
    min_graph_size = int(nCells / 6000)
    max_graph_size = int(nCells / 100)
    print(
        "Creating graph files between ",
        min_graph_size,
        " and ",
        max_graph_size)
    n_power2 = 2**np.arange(1, 20 + 1)
    n_multiples12 = 12 * np.arange(1, 8 + 1)

    n = n_power2
    for power10 in range(3):
        n = np.concatenate([n, 10**power10 * n_multiples12])

    for j in range(len(n)):
        if n[j] >= min_graph_size and n[j] <= max_graph_size:
            args = ['gpmetis', 'mpas-o.graph.info.' + date_string, str(n[j])]
            run_command(args)

    # create link to output directory
    files = glob.glob('mpas-o.graph.info.*')
    os.chdir('../assembled_files_for_upload/inputdata/ocn/mpas-o/' + mesh_name)
    for file in files:
        make_link('../../../../../graph_partition_ocean/' + file, './' + file)
# }}}


def initial_condition_seaice(config, mesh_name, date_string, ice_shelf_cavities):
# {{{

    # create links
    make_link('../init.nc', mesh_name + '.nc')

    # command line execution
    args = ['ncks', '-x', '-v',
            'bottomDepth,refBottomDepth,restingThickness, temperature,salinity,\
             temperatureSurfaceValue,salinitySurfaceValue,surfaceVelocityZonal,\
             surfaceVelocityMeridional,SSHGradientZonal,SSHGradientMeridional,\
             vertNonLocalFluxTemp,normalVelocity,layerThickness,normalBarotropicVelocity,\
             vertCoordMovementWeights,boundaryLayerDepth,seaIcePressure,\
             atmosphericPressure,filteredSSHGradientZonal,filteredSSHGradientMeridional',
            '-O',  # Overwrite existing file
            mesh_name + '.nc',
            'seaice.' + mesh_name + '.nc'
            ]
    run_command(args)

    # make links to output directory
    os.chdir(
        '../assembled_files_for_upload/inputdata/ice/mpas-cice/' +
        mesh_name)
    make_link('../../../../../initial_condition_seaice/seaice.' +
              mesh_name + '.nc', 'seaice.' + mesh_name + '.nc')
# }}}


def scrip(config, mesh_name, date_string, ice_shelf_cavities):
# {{{

    if ice_shelf_cavities:
        nomaskStr='.nomask'
    else:
        nomaskStr=''

    # create links
    make_link('../init.nc', mesh_name + '.nc')

    # command line execution
    scrip_file = ('ocean.' + mesh_name + nomaskStr + '.scrip.' + date_string + '.nc')

    args = ['create_SCRIP_file_from_MPAS_mesh.py',
        '-m', mesh_name + '.nc',
        '-s', scrip_file]
    run_command(args)

    if ice_shelf_cavities:
        scrip_file_mask = ('ocean.' + mesh_name + '.mask.scrip.' + date_string + '.nc')
        args = ['create_SCRIP_file_from_MPAS_mesh.py',
            '-m', mesh_name + '.nc',
            '-s', scrip_file_mask,
            '--landice']
        run_command(args)

    # make links to output directories
    os.chdir('../assembled_files_for_upload/inputdata/ocn/mpas-o/' + mesh_name)
    make_link('../../../../../scrip/' + scrip_file, scrip_file)
    if ice_shelf_cavities:
        make_link('../../../../../scrip/' + scrip_file_mask, scrip_file_mask)
# }}}


def transects_and_regions(config, mesh_name, date_string, ice_shelf_cavities):
# {{{

    # obtain configuration settings
    transect_region_geojson = config.get(
        'transects_and_regions',
        'transect_region_geojson')

    maskfile = 'masks_SingleRegionAtlanticWTransportTransects.' + mesh_name + '.nc'

    # make links
    make_link('../init.nc', mesh_name + '.nc')
    make_link(transect_region_geojson,
              'SingleRegionAtlanticWTransportTransects.geojson')

    # check: how does this call MpasMaskCreator?
    args = ['MpasMaskCreator.x', mesh_name + '.nc', maskfile,
            '-f', 'SingleRegionAtlanticWTransportTransects.geojson']
    run_command(args)

    # make links in output directory
    os.chdir('../assembled_files_for_upload/inputdata/ocn/mpas-o/' + mesh_name)
    make_link(
        '../../../../../transects_and_regions/masks_SingleRegionAtlanticWTransportTransects.' +
        mesh_name + '.nc',
        'masks_SingleRegionAtlanticWTransportTransects.' + mesh_name + '.nc')
# }}}


def mapping_CORE_Gcase(config, mesh_name, date_string, ice_shelf_cavities):
# {{{
    atm_scrip_tag = config.get('mapping_CORE_Gcase', 'atm_scrip_tag')
    mapping(config, mesh_name, date_string, ice_shelf_cavities, atm_scrip_tag)

    # make links in output directory
    files = glob.glob('map_*')
    os.chdir('../assembled_files_for_upload/inputdata/cpl/cpl6')
    for file in files:
        make_link('../../../../mapping_CORE_Gcase/' + file, './' + file)
# }}}


def mapping_JRA_Gcase(config, mesh_name, date_string, ice_shelf_cavities):
# {{{
    atm_scrip_tag = config.get('mapping_JRA_Gcase', 'atm_scrip_tag')
    mapping(config, mesh_name, date_string, ice_shelf_cavities, atm_scrip_tag)

    # make links in output directory
    files = glob.glob('map_*')
    os.chdir('../assembled_files_for_upload/inputdata/cpl/cpl6')
    for file in files:
        make_link('../../../../mapping_JRA_Gcase/' + file, './' + file)
# }}}


def mapping_ne30(config, mesh_name, date_string, ice_shelf_cavities):
# {{{
    atm_scrip_tag = config.get('mapping_ne30', 'atm_scrip_tag')
    mapping(config, mesh_name, date_string, ice_shelf_cavities, atm_scrip_tag)

    # make links in output directory
    files = glob.glob('map_*')
    os.chdir('../assembled_files_for_upload/inputdata/cpl/cpl6')
    for file in files:
        make_link('../../../../mapping_CORE_Gcase/' + file, './' + file)
# }}}


def mapping(config, mesh_name, date_string, ice_shelf_cavities, atm_scrip_tag):
# {{{

    # obtain configuration settings
    nprocs = config.get('main', 'nprocs')
    if ice_shelf_cavities:
        nomaskStr='.nomask'
    else:
        nomaskStr=''
    atm_scrip_path = config.get('main', 'atm_scrip_path')

    # make links
    ocn_scrip_file = 'ocean.' + mesh_name + nomaskStr + '.scrip.' + date_string + '.nc'
    make_link('../scrip/' + ocn_scrip_file, ocn_scrip_file)
    atm_scrip_file = atm_scrip_tag + '.nc'
    make_link(atm_scrip_path + '/' + atm_scrip_file, atm_scrip_file)

    # execute commands
    try:
        for method, short in [['conserve', 'aave'], ['bilinear', 'blin'], ['patch', 'patc']]:

            # Ocean to atmosphere
            args = ['mpirun', '-n', nprocs, 'ESMF_RegridWeightGen',
                    '--method', method,
                    '--source', ocn_scrip_file,
                    '--destination', atm_scrip_file,
                    '--weight', 'map_' + mesh_name + nomaskStr + '_TO_' + atm_scrip_tag + '_'
                       + short + '.' + date_string + '.nc',
                    '--ignore_unmapped']
            run_command(args)

            # Atmosphere to ocean
            args = ['mpirun', '-n', nprocs, 'ESMF_RegridWeightGen',
                    '--method', method,
                    '--source', atm_scrip_file,
                    '--destination', ocn_scrip_file,
                    '--weight', 'map_' + atm_scrip_tag + '_TO_' + mesh_name + nomaskStr + '_'
                       + short + '.' + date_string + '.nc',
                    '--ignore_unmapped']
            run_command(args)

    except OSError:
        print('mapping must be run on one compute node')

    if ice_shelf_cavities:
        print("\n Mapping files with masks for ice shelf cavities")
        # make links
        ocn_scrip_file = 'ocean.' + mesh_name + '.mask.scrip.' + date_string + '.nc'
        make_link('../scrip/' + ocn_scrip_file, ocn_scrip_file)

        # execute commands
        try:
            for method, short in [['conserve', 'aave'], ['bilinear', 'blin'], ['patch', 'patc']]:

                # Ocean to atmosphere
                args = ['mpirun', '-n', nprocs, 'ESMF_RegridWeightGen',
                        '--method', method,
                        '--source', ocn_scrip_file,
                        '--destination', atm_scrip_file,
                        '--weight', 'map_' + mesh_name + '.mask_TO_' + atm_scrip_tag + '_'
                           + short + '.' + date_string + '.nc',
                        '--ignore_unmapped']
                run_command(args)

                # Atmosphere to ocean
                args = ['mpirun', '-n', nprocs, 'ESMF_RegridWeightGen',
                        '--method', method,
                        '--source', atm_scrip_file,
                        '--destination', ocn_scrip_file,
                        '--weight', 'map_' + atm_scrip_tag + '_TO_' + mesh_name + '.mask_'
                           + short + '.' + date_string + '.nc',
                        '--ignore_unmapped']
                run_command(args)

        except OSError:
            print('mapping_CORE_Gcase must be run on one compute node')
# }}}


def domain_CORE_Gcase(config, mesh_name, date_string, ice_shelf_cavities):
# {{{
    # obtain configuration settings
    domain_exe = config.get('main', 'domain_exe')
    atm_scrip_tag = config.get('mapping_CORE_Gcase', 'atm_scrip_tag')
    if ice_shelf_cavities:
        nomaskStr='.nomask'
    else:
        nomaskStr=''

    # make links
    make_link(domain_exe, 'domain_exe')
    mapping_file = 'map_' + mesh_name + nomaskStr + '_TO_' + \
        atm_scrip_tag + '_aave.' + date_string + '.nc'
    make_link('../mapping_CORE_Gcase/' + mapping_file, mapping_file)

    # execute commands
    args = ['./domain_exe', '-m', mapping_file, '-o', mesh_name, '-l', 'T62']
    run_command(args)

    # make links in output directories
    files = glob.glob('domain*.nc')
    os.chdir('../assembled_files_for_upload/inputdata/share/domains')
    for file in files:
        make_link('../../../../domain/' + file, './' + file)
# }}}


def domain_JRA_Gcase(config, mesh_name, date_string, ice_shelf_cavities):
# {{{
    # obtain configuration settings
    domain_exe = config.get('main', 'domain_exe')
    atm_scrip_tag = config.get('mapping_JRA_Gcase', 'atm_scrip_tag')
    if ice_shelf_cavities:
        nomaskStr='.nomask'
    else:
        nomaskStr=''

    # make links
    make_link(domain_exe, 'domain_exe')
    mapping_file = 'map_' + mesh_name + nomaskStr + '_TO_' + \
        atm_scrip_tag + '_aave.' + date_string + '.nc'
    make_link('../mapping_JRA_Gcase/' + mapping_file, mapping_file)

    # execute commands
    args = ['./domain_exe', '-m', mapping_file, '-o', mesh_name, '-l', 'T62']
    run_command(args)

    # make links in output directories
    files = glob.glob('domain*.nc')
    os.chdir('../assembled_files_for_upload/inputdata/share/domains')
    for file in files:
        make_link('../../../../domain/' + file, './' + file)
# }}}


def domain_ne30(config, mesh_name, date_string, ice_shelf_cavities):
# {{{
    # obtain configuration settings
    domain_exe = config.get('main', 'domain_exe')
    atm_scrip_tag = config.get('mapping_ne30', 'atm_scrip_tag')
    if ice_shelf_cavities:
        nomaskStr='.nomask'
    else:
        nomaskStr=''

    # make links
    make_link(domain_exe, 'domain_exe')
    mapping_file = 'map_' + mesh_name + nomaskStr + '_TO_' + \
        atm_scrip_tag + '_aave.' + date_string + '.nc'
    make_link('../mapping_ne30/' + mapping_file, mapping_file)

    # execute commands
    args = ['./domain_exe', '-m', mapping_file, '-o', mesh_name, '-l', 'T62']
    run_command(args)

    # make links in output directories
    files = glob.glob('domain*.nc')
    os.chdir('../assembled_files_for_upload/inputdata/share/domains')
    for file in files:
        make_link('../../../../domain/' + file, './' + file)
# }}}


def mapping_runoff(config, mesh_name, date_string, ice_shelf_cavities):
# {{{

    print("WARNING: This works, but uses a version of runoff_map in cime at")
    print("    cime/tools/mapping/gen_mapping_files/runoff_to_ocn")
    print("    This needs to be replaced with a newer version")
    print("    -- Mark Petersen Dec 2019")

    # obtain configuration settings
    runoff_map_exe = config.get('mapping_runoff', 'runoff_map_exe')
    runoff_map_lnd_file = config.get('mapping_runoff', 'runoff_map_lnd_file')
    if ice_shelf_cavities:
        nomaskStr='.nomask'
    else:
        nomaskStr=''

    # make links
    make_link(runoff_map_exe, 'runoff_map_exe')
    make_link(runoff_map_lnd_file, 'runoff.daitren.annual.nc')
    ocn_scrip_file = ('ocean.' + mesh_name + nomaskStr + '.scrip.' + date_string + '.nc')
    make_link('../scrip/' + ocn_scrip_file, ocn_scrip_file)

    # write namelist file
    # could put eFold and rMax in ini file
    # not sure if 'coastal mask' flag applies.
    f = open("runoff_map.nml", "w+")
    f.write("&input_nml\n")
    f.write("   gridtype     = 'obs'\n")
    f.write("   file_roff    = 'runoff.daitren.annual.nc'\n")
    f.write("   file_ocn     = '" + ocn_scrip_file + "'\n")
    #f.write("   file_ocn_coastal_mask = '" + ocn_scrip_file  + "'\n")
    f.write("   file_nn      = 'map_rx1_to_" + mesh_name + "_coast_nearestdtos_" + date_string + ".nc'\n")
    f.write("   file_smooth  = 'map_" + mesh_name + "_coast_to_" + mesh_name + "_sm_e1000r300_" + date_string + ".nc'\n")
    f.write("   file_new     = 'map_rx1_to_" + mesh_name + "_nnsm_e1000r300_" + date_string + ".nc'\n")
    f.write("   title        = 'runoff map: rx1 -> " + mesh_name + ", nearest neighbor and smoothed'\n")
    f.write("   eFold        = 1000000.0\n")
    f.write("   rMax         =  300000.0\n")
    f.write("   restrict_smooth_src_to_nn_dest = .true.\n")
    f.write("   step1 = .true.\n")
    f.write("   step2 = .true.\n")
    f.write("   step3 = .true.\n")
    f.write("/\n")
    f.close()

    # execute commands
    args = ['./runoff_map_exe']
    run_command(args)

    # Alter runoff mapping so runoff does not go under ice shelves
    # WARNING: this is not hooked up yet. I need to know which mapping files this applies to.
    # Also, this is pointing to the correct -w and -n flags, but it only works if I
    # switch those files.
    if ice_shelf_cavities:
        make_link('../copy_cell_indices_ISC.py', 'copy_cell_indices_ISC.py')
        make_link('../init.nc', 'init.nc')
        make_link('../no_ISC_culled_mesh.nc', 'no_ISC_culled_mesh.nc.nc')
        args = ['./copy_cell_indices_ISC.py',
            '-i', 'map_oQU240wISC_coast_to_oQU240wISC_sm_e1000r300_200202.nc',
            '-o', 'map_output.nc',
            '-w', 'init.nc',
            '-n', 'no_ISC_culled_mesh.nc.nc'
            ]
        run_command(args)

    # make links in output directories
    files = glob.glob('map*.nc')
    os.chdir('../assembled_files_for_upload/inputdata/cpl/cpl6')
    for file in files:
        make_link('../../../../mapping_runoff/' + file, './' + file)
# }}}


def salinity_restoring(config, mesh_name, date_string, ice_shelf_cavities):
# {{{

    # obtain configuration settings
    grid_Levitus_1x1_scrip_file = config.get(
        'salinity_restoring', 'grid_Levitus_1x1_scrip_file')
    salinity_restoring_input_file = config.get(
        'salinity_restoring', 'salinity_restoring_input_file')
    nprocs = config.get('main', 'nprocs')
    if ice_shelf_cavities:
        nomaskStr='.nomask'
    else:
        nomaskStr=''

    # make links
    make_link(grid_Levitus_1x1_scrip_file, 'grid_Levitus_1x1_scrip_file.nc')
    make_link(
        salinity_restoring_input_file,
        'salinity_restoring_input_file.nc')
    ocn_scrip_file = 'ocean.' + mesh_name + nomaskStr + '.scrip.' + date_string + '.nc'
    make_link('../scrip/' + ocn_scrip_file, ocn_scrip_file)

    # execute commands
    salinity_restoring_output_file = 'sss.PHC2_monthlyClimatology.' + \
        mesh_name + '.' + date_string + '.nc'
    try:
        for method, short in [['bilinear', 'blin']]:

            # mapping file, 1x1 to ocean mesh
            map_Levitus_file = 'map_' + 'Levitus_1x1' + '_TO_' + \
                mesh_name + '_' + short + '.' + date_string + '.nc'
            args = ['mpirun', '-n', nprocs, 'ESMF_RegridWeightGen',
                    '--method', method,
                    '--source', 'grid_Levitus_1x1_scrip_file.nc',
                    '--destination', ocn_scrip_file,
                    '--weight', map_Levitus_file,
                    '--ignore_unmapped']
            run_command(args)

    except OSError:
        print('salinity_restoring must be run on compute node')

    # remap from 1x1 to model grid
    args = ['ncremap',
        '-i', 'salinity_restoring_input_file.nc',
        '-o', 'intermediate_file.nc',
        '-m', map_Levitus_file,
        '--no_cll_msr', '--no_frm_trm', '--no_stg_grd']
    run_command(args)

    # Remove all bounds attributes.  This is necessary to remove the lon, lat
    # of vertices.
    args = ['ncatted', '-a', 'bounds,,d,,', 'intermediate_file.nc']
    run_command(args)

    # Remove lon, lat of vertices since they take a lot of space.
    args = ['ncks', '-O', '-x', '-v', 'lat_vertices,lon_vertices',
            'intermediate_file.nc', salinity_restoring_output_file]
    run_command(args)

    # Rename ncol and SALT to be consistent with what MPAS-Ocean expects.
    args = [
        'ncrename',
        '-d', 'ncol,nCells',
        '-v', 'SALT,surfaceSalinityMonthlyClimatologyValue',
        salinity_restoring_output_file]
    run_command(args)

    # make links in output directories
    os.chdir('../assembled_files_for_upload/inputdata/ocn/mpas-o/' + mesh_name)
    make_link(
        '../../../../../salinity_restoring/' + salinity_restoring_output_file,
        './' + salinity_restoring_output_file)
# }}}


def make_dir(dirName):
# {{{
    try:
        os.makedirs(dirName)
    except OSError:
        pass
# }}}


def make_link(source, linkName):
# {{{
    try:
        if os.path.exists(linkName):
            os.remove(linkName)
        os.symlink(source, linkName)
    except OSError:
        pass
# }}}


def write_command_history(text):
# {{{
    try:
        print(text)
        with open('command_history', 'a') as outstream:
            outstream.write(text + '\n')
    except OSError:
        pass
# }}}


def run_command(args):
# {{{
    try:
        write_command_history(' '.join(args))
        with open('log.out', 'a') as outstream:
            outstream.write('Command: ' + ' '.join(args) + '\n')
            subprocess.check_call(args, stdout=outstream, stderr=outstream)
            outstream.write('\n')
    except OSError:
        pass
# }}}


if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
