#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import getpass
import datetime
import importlib
import traceback
import cdp.cdp_run
from acme_diags.acme_parser import ACMEParser
from acme_diags.acme_viewer import create_viewer
from acme_diags.driver.utils import get_set_name
import json
import collections
import csv


def create_lat_lon_metrics_table(root_dir, parameters):

    metrics_info = collections.OrderedDict()
    for parameter in parameters:
        for set_num in parameter.sets:
            set_num =  get_set_name(set_num)
            print(set_num)
            if set_num in ['lat_lon','5']:
                
                # ref_name-variable-season-region
                # or
                # ref_name-variable-plev'mb'-season-region
                ref_name = parameter.ref_name
                
                for var in parameter.variables:
                    for season in parameter.seasons:
                        for region in parameter.regions:
                            # because some parameters have plevs, there might be
                            # more than one row_name, fnm pair
                            row_name_and_fnm = []
    
                            if parameter.plevs == []:  # 2d variables
                                row_name = '{} {} {}'.format(var, region, ref_name)
                                fnm = '{}-{}-{}-{}'.format(ref_name,
                                                           var, season, region)
                                row_name_and_fnm.append((row_name, fnm))
                            else:  # 3d variables
                                for plev in parameter.plevs:
                                    row_name = '{} {} {} {}'.format(
                                        var, str(int(plev)) + ' mb', region, ref_name)
                                    fnm = '{}-{}-{}-{}-{}'.format(
                                        ref_name, var, int(plev), season, region)
                                    row_name_and_fnm.append((row_name, fnm))
                            print(row_name_and_fnm)
                            metrics_path = os.path.join(
                                #'..', '{}'.format(set_num), parameter.case_id, fnm)
                               parameter.results_dir, '{}'.format(set_num), parameter.case_id, fnm)
                            try:
                                metrics_dic = json.load(open(metrics_path + '.json'))
                            except Exception as e:
                                print(e)
                                continue
                     
                        if season not in metrics_info:
                            metrics_info[season] = collections.OrderedDict()
                        if row_name not in metrics_info[season]:
                            metrics_info[season][row_name] = collections.OrderedDict()
                        metrics_info[season][row_name]['metrics'] = metrics_dic
                            
                            

                # save metrics information in .csv table
                header = ['Variables','Unit','Model mean','Obs mean','Mean Bias','RMSE','correlation']
                for season in parameter.seasons:
                    table_path = os.path.abspath(os.path.join(
                                #'..', '{}'.format(set_num), parameter.case_id, fnm)
                               parameter.results_dir, '{}'.format(set_num)))
                    print(table_path)
                    with open(table_path + '/' + season + '_metrics_table.csv','w') as f1:
                        writer=csv.writer(f1, delimiter=',',lineterminator='\n', quoting=csv.QUOTE_NONE)
                        writer.writerow(header)
                        for key, metrics_dic in metrics_info[season].items():
                            metrics = metrics_dic['metrics']
                            row = [key, metrics['unit'], round(metrics['test_regrid']['mean'],3),round(metrics['ref_regrid']['mean'],3), round(metrics['test_regrid']['mean'] - metrics['ref_regrid']['mean'],3), round(metrics['misc']['rmse'],3), round(metrics['misc']['corr'],3)]
                            writer.writerow(row)

                    # convert csv to html

                    read_csv = csv.reader(open(table_path + '/' + season + '_metrics_table.csv'))
                    htmlfile = open(table_path + '/' + season + '_metrics_table.html','w+')
                    htmlfile.write('<p><th><b>'+ season + ' Mean' + '</b></th></p>')
                    # initialize rownum variable
                    rownum = 0
                    # write <table> tag
                    htmlfile.write('<table>')
                    # generate table contents
                
                    for row in read_csv: # Read a single row from the CSV file
                
                     # write header row. assumes first row in csv contains header
                         if rownum == 0:
                            htmlfile.write('<tr>') # write <tr> tag
                            for column in row:
                                htmlfile.write('<th>' + column +'</th>')
                            htmlfile.write('</tr>')
                
                      #  write all other rows 
                         else:
                            htmlfile.write('<tr><div style="width: 50px" >')
                            #htmlfile.write('<tr>')    
                            for column in row:
                                htmlfile.write('<td>' + column +'</td>')
                            #htmlfile.write('</tr>')
                            htmlfile.write('</div></tr>')
                         #increment row count 
                         rownum += 1
                   # write </table> tag
                    htmlfile.write('</table>') 
                    



def _get_default_diags(set_num, dataset):
    """Returns the path from the json corresponding to set_num"""
    set_num = get_set_name(set_num)
    folder = '{}'.format(set_num)
    fnm = '{}_{}.json'.format(set_num, dataset)
    pth = os.path.join(sys.prefix, 'share', 'acme_diags', folder, fnm)
    print('Using {} for {}'.format(pth, set_num))
    if not os.path.exists(pth):
        raise RuntimeError(
            "Plotting via set '{}' not supported, file {} not installed".format(set_num, fnm))
    return pth


def _collaspse_results(parameters):
    """When using cdp_run, parameters is a list of lists: [[Parameters], ...].
       Make this just a list: [Parameters]."""
    output_parameters = []

    for p1 in parameters:
        if isinstance(p1, list):
            for p2 in p1:
                output_parameters.append(p2)
        else:
            output_parameters.append(p1)

    return output_parameters


def run_diag(parameters):
    """For a single set of parameters, run the corresponding diags."""
    results = []
    for pset in parameters.sets:
        set_name = get_set_name(pset)

        parameters.current_set = set_name
        mod_str = 'acme_diags.driver.{}_driver'.format(set_name)
        try:
            module = importlib.import_module(mod_str)
            print('Starting to run ACME Diagnostics.')
            single_result = module.run_diag(parameters)
            results.append(single_result)
        except Exception as e:
            print('Error in {}'.format(mod_str))
            traceback.print_exc()
            if parameters.debug:
                sys.exit()

    return results


if __name__ == '__main__':
    parser = ACMEParser()
    args = parser.view_args()

    if args.parameters and not args.other_parameters:  # -p only
        original_parameter = parser.get_orig_parameters(default_vars=False)

        if not hasattr(original_parameter, 'sets'):
            original_parameter.sets = [
                'lat_lon', 'polar', 'zonal_mean_xy', 'zonal_mean_2d', 'cosp_histogram']

        # load the default jsons
        default_jsons_paths = []

        for set_num in original_parameter.sets:
            datasets = ['ACME']
            if hasattr(original_parameter, 'datasets'):
                datasets = original_parameter.datasets
            for ds in datasets:
                default_jsons_paths.append(_get_default_diags(set_num, ds))
        other_parameters = parser.get_other_parameters(
            files_to_open=default_jsons_paths, check_values=False)
        # Don't put the sets from the Python parameters to the default.
        # Ex. if sets=[5, 7] in the Python parameters, don't change sets in the
        # default jsons like lat_lon_AMWG_default.json
        vars_to_ignore = ['sets']
        parameters = parser.get_parameters(
            orig_parameters=original_parameter, other_parameters=other_parameters, vars_to_ignore=vars_to_ignore)

    elif not args.parameters and args.other_parameters:  # -d only
        other_parameters = parser.get_other_parameters(check_values=True)
        parameters = parser.get_parameters(other_parameters=other_parameters)

    elif args.parameters and args.other_parameters:  # -p and -d
        original_parameter = parser.get_orig_parameters(default_vars=False)
        other_parameters = parser.get_other_parameters(check_values=False)
        parameters = parser.get_parameters(
            orig_parameters=original_parameter, other_parameters=other_parameters)

    else:
        raise RuntimeError('You tried running the diags without -p and/or -d')

    dt = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    for p in parameters:
        # needed for distributed running
        # chown of all generated files to the user who ran the diags
        p.user = getpass.getuser()

        if not hasattr(p, 'results_dir'):
            p.results_dir = '{}-{}'.format('acme_diags_results', dt)

    if not os.path.exists(parameters[0].results_dir):
        os.makedirs(parameters[0].results_dir, 0o775)
    if parameters[0].multiprocessing:
        parameters = cdp.cdp_run.multiprocess(run_diag, parameters)
    elif parameters[0].distributed:
        parameters = cdp.cdp_run.distribute(run_diag, parameters)
    else:
        parameters = cdp.cdp_run.serial(run_diag, parameters)

    parameters = _collaspse_results(parameters)

    if parameters:
        pth = os.path.join(parameters[0].results_dir, 'viewer')
        if not os.path.exists(pth):
            os.makedirs(pth)

        create_lat_lon_metrics_table(pth, parameters)

        create_viewer(pth, parameters, parameters[0].output_format[0])
    else:
        print('There was not a single valid diagnostics run, no viewer created')
