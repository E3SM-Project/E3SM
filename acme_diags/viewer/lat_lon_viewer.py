"""
Functionality to create the Taylor diagrams and
metrics table for the Latitude-Longitude set.
"""

import os
import csv
import json
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import acme_diags
from acme_diags.plot.cartopy.taylor_diagram import TaylorDiagram
from . import utils


def generate_lat_lon_metrics_table(lat_lon_table_info, seasons, viewer, root_dir, parameters):
    """
    For each season in lat_lon_table_info, create a csv,
    convert it to an html and append that html to the viewer.
    """
    table_dir = os.path.join(root_dir, 'table-data')  # output_dir/viewer/table-data

    if not os.path.exists(table_dir):
        os.mkdir(table_dir)

    for season in lat_lon_table_info:
        test_name = parameters[0].short_test_name if parameters[0].short_test_name else parameters[0].test_name
        if parameters[0].run_type == 'model_vs_obs':
            ref_name = 'Observation and Reanalysis'
        else:
            ref_name = parameters[0].short_ref_name if parameters[0].short_ref_name else parameters[0].ref_name            
        csv_path = _create_csv_from_dict(lat_lon_table_info, table_dir, season, test_name, parameters[0].run_type)
        html_path = _cvs_to_html(csv_path, season, test_name, ref_name)

        # Ex: change this: /Users/zshaheen/output_dir/viewer/table-data/ANN_metrics_table.html
        # to this: viewer/table-data/ANN_metrics_table.html
        html_path = '/'.join(html_path.split('/')[-3:])

        lat_lon_table_info[season]['html_path'] = html_path

    url = _create_lat_lon_table_index(lat_lon_table_info, seasons, viewer, root_dir)
    utils.add_header(root_dir, os.path.join(root_dir, url), parameters)
    _edit_table_html(lat_lon_table_info, seasons, root_dir)

    return 'Table', url


def _create_csv_from_dict(lat_lon_table_info, output_dir, season, test_name, run_type):
    """
    Create a csv for a season in lat_lon_table_info
    in output_dir and return the path to it.
    """
    table_path = os.path.join(output_dir, season + '_metrics_table.csv')

    col_names = ['Variables', 'Unit', 'Test_mean',
                'Ref._mean', 'Mean_Bias', 'Test_STD',
                'Ref._STD', 'RMSE', 'Correlation']

    with open(table_path, 'w') as table_csv:
        writer = csv.writer(table_csv, delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE)
        writer.writerow(col_names)
        for key, metrics_dic in list(lat_lon_table_info[season].items()):
            metrics = metrics_dic['metrics']
            if run_type == 'model_vs_model':
                key = key.split()[0] + ' ' + key.split()[1]
            row = [key, metrics['unit'], round(metrics['test_regrid']['mean'], 3),
                    round(metrics['ref_regrid']['mean'], 3),
                    round(metrics['test_regrid']['mean'] - metrics['ref_regrid']['mean'], 3),
                    round(metrics['test_regrid']['std'], 3), round(metrics['ref_regrid']['std'], 3),
                    round(metrics['misc']['rmse'], 3), round(metrics['misc']['corr'], 3)]
            writer.writerow(row)

    return table_path


def _cvs_to_html(csv_path, season, test_name, ref_name):
    """
    Convert the csv for a season located at csv_path
    to an HTML, returning the path to the HTML.
    """
    html_path = csv_path.replace('csv', 'html')

    with open(html_path, 'w') as htmlfile:
        htmlfile.write('<p><b>Test: {}</b><br>'.format(test_name))
        htmlfile.write('<b>Reference: {}</b></p>'.format(ref_name))
        htmlfile.write('<p><th><b>{} Mean </b></th></p>'.format(season))
        htmlfile.write('<table>')

        with open(csv_path) as csv_file:
            read_csv = csv.reader(csv_file)

            # Generate the table's contents.
            for num, row in enumerate(read_csv):
                # Write the header row, assuming the first
                # row in csv contains the header.
                if num == 0:
                    htmlfile.write('<tr>')
                    for column in row:
                        htmlfile.write('<th>{}</th>'.format(column))
                    htmlfile.write('</tr>')

                # Write all other rows.
                else:
                    htmlfile.write('<tr><div style="width: 50px">')
                    for column in row:
                        htmlfile.write('<td>{}</td>'.format(column))
                    htmlfile.write('</div></tr>')

        htmlfile.write('</table>')

    return html_path


def _create_lat_lon_table_index(lat_lon_table_info, seasons, viewer, root_dir):
    """
    Create an index in the viewer that links the
    individual htmls for the lat-lon table.
    """
    viewer.add_page('Table', seasons)
    viewer.add_group('Summary Table')
    viewer.add_row('All variables')

    for s in seasons:
        if s in lat_lon_table_info:
            viewer.add_col(lat_lon_table_info[s]['html_path'], is_file=True, title=s)
        else:
            viewer.add_col('-----', is_file=True, title='-----')

    url = viewer.generate_page()
    return url


def generate_lat_lon_taylor_diag(lat_lon_table_info, seasons, viewer, root_dir, parameters):
    """
    For each season in lat_lon_table_info, create a csv, plot using
    taylor diagrams and append that html to the viewer.
    """
    taylor_diag_dir = os.path.join(root_dir, 'taylor-diagram-data')  # output_dir/viewer/taylor-diagram-data

    if not os.path.exists(taylor_diag_dir):
        os.mkdir(taylor_diag_dir)

    season_to_png = {}
    for season in lat_lon_table_info:
        test_name = parameters[0].short_test_name if parameters[0].short_test_name else parameters[0].test_name
        if parameters[0].run_type == 'model_vs_obs':
            ref_name = 'Observation and Reanalysis'
        else:
            ref_name = parameters[0].short_ref_name if parameters[0].short_ref_name else parameters[0].ref_name            

        csv_path = _create_csv_from_dict_taylor_diag(lat_lon_table_info, taylor_diag_dir,
            season, test_name, parameters[0].run_type, ref_name)
        # Remove any reference to the results_dir when inserting the links into HTML pages.
        # This is because that folder can be renamed.
        csv_path = csv_path.split('viewer')[-1]
        csv_path = 'viewer' + csv_path

        season_to_png[season] = csv_path.replace('csv', 'png')

    url = _create_taylor_index(seasons, viewer, root_dir, season_to_png)
    utils.add_header(root_dir, os.path.join(root_dir, url), parameters)

    return 'Taylor Diagram', url


def _create_csv_from_dict_taylor_diag(lat_lon_table_info, output_dir, season, test_name, run_type, ref_name):
    """
    Create a csv for a season in lat_lon_table_info in
    output_dir and return the path to it.
    
    Since the Taylor diagram uses the same seasons
    as lat_lon_table_info, we can use that.
    """
    taylor_diag_path = os.path.join(output_dir, '{}_metrics_taylor_diag.csv'.format(season))
    control_runs_path = os.path.join(acme_diags.INSTALL_PATH, 'control_runs',
                        '{}_metrics_taylor_diag_B1850_v0.csv'.format(season))

    col_names = ['Variables', 'Test_STD', 'Ref._STD', 'Correlation']

    with open(taylor_diag_path, 'w') as table_csv:
        writer=csv.writer(table_csv, delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE)
        writer.writerow(col_names)

        for key, metrics_dic in list(lat_lon_table_info[season].items()):
            # Only include variables from a certain list in the Taylor diagram.
            var_list_1 = ['PRECT', 'PSL', 'SWCF', 'LWCF', 'TREFHT']
            var_list_2 = ['PRECT_GPCP', 'PSL_ERA-Interim', 'SWCF_ceres', 'LWCF_ceres', 'TREFHT_CRU']
            if run_type == 'model_vs_obs':
                if key.split()[0] in var_list_1 and '_'.join((key.split()[0], key.split()[2].split('_')[0])) in var_list_2:
                    metrics = metrics_dic['metrics']
                    row = [key, round(metrics['test_regrid']['std'], 3), round(metrics['ref_regrid']['std'], 3),
                        round(metrics['misc']['corr'], 3)]
                    writer.writerow(row)
            else:
                if key.split()[0] in var_list_1:
                    metrics = metrics_dic['metrics']
                    row = [key, round(metrics['test_regrid']['std'], 3), round(metrics['ref_regrid']['std'], 3),
                        round(metrics['misc']['corr'], 3)]
                    writer.writerow(row)
    
    with open(taylor_diag_path, 'r') as taylor_csv:
        reader = csv.reader(taylor_csv, delimiter=",")
        data = list(reader)
        row_count = len(data)
   
    # Read the control run data.
    with open(control_runs_path, 'r') as control_runs_taylor_csv:
        reader = csv.reader(control_runs_taylor_csv, delimiter=",")
        control_runs_data = list(reader)

    keys_control_runs = []
    for i in range(0, len(control_runs_data)):
        keys_control_runs.append(control_runs_data[i][0])
   
    # Generate Taylor diagram plot if there is metrics
    # saved for any variable within the list.
    marker = ['o', 'd', '+', 's', '>', '<', 'v', '^', 'x', 'h', 'X', 'H'] 
    color = ['k', 'r', 'g', 'y', 'm']
 
    if row_count > 0:
        matplotlib.rcParams.update({'font.size': 20})
        fig = plt.figure(figsize=(9,8))
        refstd = 1.0
        taylordiag = TaylorDiagram(refstd, fig=fig, rect=111, label="REF")
        ax = taylordiag._ax
        
        # Add samples to Taylor diagram.
        for irow in range(1,row_count):
            std_norm, correlation = float(data[irow][1])/float(data[irow][2]), float(data[irow][3])
            taylordiag.add_sample(std_norm, correlation, marker=marker[irow], c=color[0], ms=10,
                label=data[irow][0], markerfacecolor='None', markeredgecolor=color[0], linestyle='None')

        # Add a legend to the figure.
        fig.legend(taylordiag.samplePoints,
                   [p.get_label() for p in taylordiag.samplePoints],
                   numpoints=1, loc='center right', bbox_to_anchor=(1.0, .5), prop={'size':10})

        # Add samples for baseline simulation.
        if run_type == 'model_vs_obs':
            for irow in range(1,row_count):
                if data[irow][0] in keys_control_runs:
                    control_irow = keys_control_runs.index(data[irow][0])
                    std_norm = float(control_runs_data[control_irow][1])/float(control_runs_data[control_irow][2])
                    correlation = float(control_runs_data[control_irow][3])
                    taylordiag.add_sample(std_norm, correlation, marker=marker[irow], c=color[1],
                        ms=10, label=data[irow][0]+'E3sm_v0 B1850',
                        markerfacecolor='None', markeredgecolor=color[1], linestyle='None')
    
                baseline_text = 'E3SMv0_B1850'
                ax.text(0.6, 0.95, baseline_text, ha='left', va='center', transform=ax.transAxes, color=color[1], fontsize=12)

        model_text = 'Test Model: ' + test_name
        ax.text(0.6, 1, model_text, ha='left', va='center', transform=ax.transAxes, color=color[0], fontsize=12)
        if run_type == 'model_vs_model':
            ax.text(0.6, 0.95, 'Ref. Model: ' + ref_name, ha='left', va='center', transform=ax.transAxes, color='k', fontsize=12)

        plt.title(season + ': Spatial Variability', y=1.08)
        fig.savefig(os.path.join(output_dir, season + '_metrics_taylor_diag.png'))

    return taylor_diag_path


def _create_taylor_index(seasons, viewer, root_dir, season_to_png):
    """
    Create an index in the viewer that links the
    individual htmls for the lat-lon table.
    """
    viewer.add_page('Taylor Diagram', seasons)
    viewer.add_group('Summary Taylor Diagrams')
    viewer.add_row('All variables')

    for s in seasons:
        if s in season_to_png:
            pth = os.path.join('..', season_to_png[s])
            viewer.add_col(pth, is_file=True, title=s)
        else:
            viewer.add_col('-----', is_file=True, title='-----')

    url = viewer.generate_page()
    return url


def _edit_table_html(lat_lon_table_info, seasons, root_dir):
    """
    After the viewer is created, edit the table html to
    insert the custom htmls.
    """
    def _add_html_to_col(season, season_path, html_path):
        """
        Since the output viewer doesn't support html images, do this hack.
        For the col in the html at html_path, insert the link to col_path.
        """
        # Change:
        # <tr class="output-row">
        #  <!-- ... -->
        #  <td colspan="1">
        #   <!-- what needs to be changed -->
        #  </td>
        # <!-- ... -->
        # </tr>
        # to:
        # <tr class="output-row">
        #  <!-- ... -->
        #  <td colspan="1">
        #   <a href="{season_path}"> {season} </a> <!-- this was changed -->
        #  </td>
        # <!-- ... -->
        # </tr>

        soup = BeautifulSoup(open(html_path), "lxml")

        for tr in soup.find_all("tr", {"class": "output-row"}):
            lst = ['All variables'] + seasons
            index = lst.index(season)
            cols = tr.find_all("td")  # the cols are ['All variables', 'ANN', 'DJF', 'MAM', 'JJA', 'SON']
            td = cols[index]  # get the HTML element related to the season

            url = os.path.join('..', '..', season_path)
            a = soup.new_tag("a", href=url)
            a.append(season)

            td.string = ''
            td.append(a)

        html = soup.prettify("utf-8")
        with open(html_path, "wb") as f:
            f.write(html)

    for s in seasons:
        if s in lat_lon_table_info:
            _add_html_to_col(s, lat_lon_table_info[s]['html_path'],
                os.path.join(root_dir, 'table', 'index.html'))
