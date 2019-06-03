import os
import sys
import fnmatch
import datetime
import shutil
import collections
import copy
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from cdp.cdp_viewer import OutputViewer
import acme_diags
from acme_diags.driver.utils.general import get_set_name
from acme_diags.plot.cartopy.taylor_diagram import TaylorDiagram
from acme_diags import container

# Dict of
# {
#   set_num: {
#       case_id: {
#           row_name: {
#               season: filename
#           }
#       }
#   }
# }
# Order is page(set_num) -> group(case_id) -> row(row_name) -> col(season) -> filename
# Used when actually inserting the cols into the viewer
# needed so we can have a cols in order of ANN, DJF, MAM, JJA, SON
ROW_INFO = collections.OrderedDict()

# A similar dict for creating the lat-lon tables
LAT_LON_TABLE_INFO = collections.OrderedDict()

def _copy_acme_logo(root_dir):
    """Copy over e3sm_logo.png to root_dir/viewer"""
    src_pth = os.path.join(acme_diags.INSTALL_PATH, "e3sm_logo.png")
    dst_path = os.path.join(root_dir, "viewer")
    shutil.copy(src_pth, dst_path)


def _get_acme_logo_path(root_dir, html_path):
    """Based of the root dir of the viewer and the current
    dir of the html, get the relative path of the E3SM logo"""
    # when current_dir = myresults-07-11/viewer/index.html, the image is in myresults-07-11/viewer/viewer/e3sm_logo.png
    # so there's no need to move some number of directories up.
    # That's why we have - 3
    relative_dir = html_path.replace(root_dir + '/', '')
    dirs_to_go_up = len(relative_dir.split('/')) - 1
    pth = os.path.join('.')
    for _ in range(0, dirs_to_go_up):
        pth = os.path.join(pth, '..')
    pth = os.path.join(pth, 'viewer', 'e3sm_logo.png')
    return pth


def _add_header(path, version, test, ref, time, logo_path):
    """Add the header to the html located at path"""

    # We're inserting the following in the body under navbar navbar-default
    # <div id="e3sm-header" style="background-color:#dbe6c5; float:left; width:45%">
    # 	<p style="margin-left:5em">
    # 		<b>E3SM Diagnostics Package [VERSION]</b><br>
    # 		Test: [SOMETHING]<br>
    # 		Reference: [SOMETHING]<br>
    # 		Created: [DATE]<br>
    # 	</p>
    # </div>
    # <div id="e3sm-header2" style="background-color:#dbe6c5; float:right; width:55%">
    # 	<img src="e3sm_logo.png" alt="logo" style="width:201px; height:91px; background-color:#dbe6c5">
    # </div>

    soup = BeautifulSoup(open(path), "lxml")
    old_header = soup.find_all("nav", "navbar navbar-default")
    if len(old_header) is not 0:
        old_header[0].decompose()

    header_div = soup.new_tag(
        "div", id="e3sm-header", style="background-color:#dbe6c5; float:left; width:45%")
    p = soup.new_tag("p", style="margin-left:5em")

    bolded_title = soup.new_tag("b")
    bolded_title.append("E3SM Diagnostics Package {}".format(version))
    bolded_title.append(soup.new_tag("br"))
    p.append(bolded_title)

    p.append("Test: {}".format(test))
    p.append(soup.new_tag("br"))

    p.append("Reference: {}".format(ref))
    p.append(soup.new_tag("br"))

    p.append("Created: {}".format(time))

    header_div.append(p)
    soup.body.insert(0, header_div)

    img_div = soup.new_tag("div", id="e3sm-header2",
                           style="background-color:#dbe6c5; float:right; width:55%")
    img = soup.new_tag("img", src=logo_path, alt="logo",
                       style="width:201px; height:91px; background-color:#dbe6c5")
    img_div.append(img)
    soup.body.insert(1, img_div)

    html = soup.prettify("utf-8")
    with open(path, "wb") as f:
        f.write(html)


def h1_to_h3(path):
    """Change any <h1> to <h3> because h1 is just too big"""
    soup = BeautifulSoup(open(path), "lxml")
    h1 = soup.find('h1')
    if h1 is None:
        return
    h1.name = 'h3'

    html = soup.prettify("utf-8")
    with open(path, "wb") as f:
        f.write(html)


def _extras(root_dir, parameters):
    """Add extras (header, etc) to the generated htmls"""
    _copy_acme_logo(root_dir)

    index_files = []
    for root, _, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, 'index.html'):
            pth = os.path.join(root, filename)
            index_files.append(pth)
    dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    os.path.join('viewer', 'e3sm_logo.png')

    for f in index_files:
        path = _get_acme_logo_path(root_dir, f)

        test_name = parameters[0].short_test_name if parameters[0].short_test_name else parameters[0].test_name
        if parameters[0].run_type == 'model_vs_obs':
            ref_name = 'Observation and Reanalysis'
        else:
            ref_name = parameters[0].short_ref_name if parameters[0].short_ref_name else parameters[0].ref_name            
        _add_header(f, acme_diags.__version__, test_name, ref_name, dt, path)
        h1_to_h3(f)

    _edit_table_html(root_dir)
    _add_table_and_taylor_to_viewer_index(root_dir)
    _add_provenance_to_viewer_index(root_dir)


def _add_pages_and_top_row(viewer, parameters):
    """Add the page and columns of the page"""
    set_to_seasons = collections.OrderedDict()  # dict of {set: [seasons]}
    for p in parameters:
        for set_num in p.sets:
            set_num = get_set_name(set_num)
            for ssn in p.seasons:
                if set_num not in set_to_seasons:
                    set_to_seasons[set_num] = []
                set_to_seasons[set_num].append(ssn)

    for set_num, seasons in list(set_to_seasons.items()):
        ROW_INFO[set_num] = collections.OrderedDict()
        col_labels = ['Description']
        for s in ['ANN', 'DJF', 'MAM', 'JJA', 'SON']:
            if s in seasons:
                col_labels.append(s)
        viewer.add_page("{}".format(_better_page_name(set_num)), col_labels)

    # Add the provenance page to the end.
    viewer.add_page('Provenance', [])

def _get_description(var, parameters):
    """Get the description for the variable from the parameters"""
    if hasattr(parameters, 'viewer_descr') and var in parameters.viewer_descr:
        return parameters.viewer_descr[var]
    return var


def _better_page_name(old_name):
    """Use a longer, more descriptive name for the pages."""
    if old_name == 'zonal_mean_xy':
        return 'Zonal mean line plots'
    elif old_name == 'zonal_mean_2d':
        return 'Pressure-Latitude zonal mean contour plots'
    elif old_name == 'lat_lon':
        return 'Latitude-Longitude contour maps'
    elif old_name == 'polar':
        return 'Polar contour maps'
    elif old_name == 'cosp_histogram':
        return 'CloudTopHeight-Tau joint histograms'
    elif old_name == 'meridional_mean_2d':
        return 'Pressure-Longitude meridional mean contour plots'
    else:
        return old_name
    
def _add_to_lat_lon_metrics_table(metrics_path, season, row_name):
    """Add the metrics for the current season and row_name to the lat-lon table"""
    with open(metrics_path + '.json') as json_file:
        metrics_dict = json.load(json_file)

        if season not in LAT_LON_TABLE_INFO:
            LAT_LON_TABLE_INFO[season] = collections.OrderedDict()
        if row_name not in LAT_LON_TABLE_INFO[season]:
            LAT_LON_TABLE_INFO[season][row_name] = collections.OrderedDict()
        LAT_LON_TABLE_INFO[season][row_name]['metrics'] = metrics_dict

def _create_csv_from_dict(output_dir, season, test_name, run_type):
    """Create a csv for a season in LAT_LON_TABLE_INFO in output_dir and return the path to it"""
    table_path = os.path.join(output_dir, season + '_metrics_table.csv')

    col_names = ['Variables', 'Unit', 'Test_mean', 'Ref._mean', 'Mean_Bias', 'Test_STD', 'Ref._STD', 'RMSE', 'Correlation']

    with open(table_path, 'w') as table_csv:
        writer=csv.writer(table_csv, delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE)
        writer.writerow(col_names)
        for key, metrics_dic in list(LAT_LON_TABLE_INFO[season].items()):
            metrics = metrics_dic['metrics']
            if run_type == 'model_vs_model':
                key = key.split()[0] + ' ' +key.split()[1]
            row = [key, metrics['unit'], round(metrics['test_regrid']['mean'],3), round(metrics['ref_regrid']['mean'],3), round(metrics['test_regrid']['mean'] - metrics['ref_regrid']['mean'],3), round(metrics['test_regrid']['std'],3), round(metrics['ref_regrid']['std'],3),round(metrics['misc']['rmse'],3), round(metrics['misc']['corr'],3)]
            writer.writerow(row)

    return table_path

def _create_csv_from_dict_taylor_diag(output_dir, season, test_name, run_type, ref_name):
    """Create a csv for a season in LAT_LON_TABLE_INFO in output_dir and return the path to it.
    Since the Taylor Diagram uses the same seasons as LAT_LON_TABLE_INFO, we can use that."""
    taylor_diag_path = os.path.join(output_dir, '{}_metrics_taylor_diag.csv'.format(season))
    control_runs_path =  os.path.join(acme_diags.INSTALL_PATH, 'control_runs', '{}_metrics_taylor_diag_B1850_v0.csv'.format(season))

    col_names = ['Variables', 'Test_STD', 'Ref._STD', 'Correlation']

    with open(taylor_diag_path, 'w') as table_csv:
        writer=csv.writer(table_csv, delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE)
        writer.writerow(col_names)

        for key, metrics_dic in list(LAT_LON_TABLE_INFO[season].items()):
            # only include variables in a a certain list for taylor diagram
            if run_type == 'model_vs_obs':
                if key.split()[0] in ['PRECT', 'PSL', 'SWCF', 'LWCF', 'TREFHT'] and '_'.join((key.split()[0], key.split()[2].split('_')[0])) in ['PRECT_GPCP','PSL_ERA-Interim','SWCF_ceres','LWCF_ceres', 'TREFHT_CRU']:
                    metrics = metrics_dic['metrics']
                    row = [key, round(metrics['test_regrid']['std'],3), round(metrics['ref_regrid']['std'],3), round(metrics['misc']['corr'],3)]
                    writer.writerow(row)
            else:
                if key.split()[0] in ['PRECT', 'PSL', 'SWCF', 'LWCF', 'TREFHT']:
                    metrics = metrics_dic['metrics']
                    row = [key, round(metrics['test_regrid']['std'],3), round(metrics['ref_regrid']['std'],3), round(metrics['misc']['corr'],3)]
                    writer.writerow(row)
                  
    
    with open(taylor_diag_path, 'r') as taylor_csv:
        reader = csv.reader(taylor_csv, delimiter = ",")
        data =  list(reader)
        row_count=len(data)
   
    # read control run data
    with open(control_runs_path, 'r') as control_runs_taylor_csv:
        reader = csv.reader(control_runs_taylor_csv, delimiter = ",")
        control_runs_data =  list(reader)

    keys_control_runs = []
    for i in range(0, len(control_runs_data)):
        keys_control_runs.append(control_runs_data[i][0])
    #print keys_control_runs
   
    # generate taylor diagram plot if there is metrics saved for any variable within the list.
    marker = ['o', 'd', '+', 's', '>', '<', 'v' , '^', 'x', 'h', 'X', 'H'] 
    color = ['k', 'r', 'g', 'y', 'm']
 
    if row_count > 0:
        matplotlib.rcParams.update({'font.size': 20})
        fig = plt.figure(figsize=(9,8))
        refstd =  1.0
        taylordiag = TaylorDiagram(refstd, fig = fig,rect = 111, label = "REF")
        ax = taylordiag._ax
        
        # Add samples to taylor diagram
        for irow in range(1,row_count):
            std_norm, correlation = float(data[irow][1])/float(data[irow][2]), float(data[irow][3])
            taylordiag.add_sample(std_norm, correlation, marker = marker[irow], c = color[0],ms = 10, label = data[irow][0], markerfacecolor = 'None', markeredgecolor = color[0], linestyle = 'None')

        # Add a figure legend
        fig.legend(taylordiag.samplePoints,
                   [p.get_label() for p in taylordiag.samplePoints],
                   numpoints=1, loc='center right', bbox_to_anchor=(1.0, .5), prop={'size':10})


        # Add samples for baseline simulation:
        if run_type == 'model_vs_obs':
            for irow in range(1,row_count):
                if data[irow][0] in keys_control_runs:
                    control_irow = keys_control_runs.index(data[irow][0])
                    std_norm, correlation = float(control_runs_data[control_irow][1])/float(control_runs_data[control_irow][2]), float(control_runs_data[control_irow][3])
                    taylordiag.add_sample(std_norm, correlation, marker = marker[irow], c = color[1],ms = 10, label = data[irow][0]+'E3sm_v0 B1850', markerfacecolor = 'None', markeredgecolor = color[1], linestyle = 'None')
    
                baseline_text = 'E3SMv0_B1850'
                ax.text(0.6, 0.95, baseline_text, ha='left', va='center', transform=ax.transAxes,color=color[1], fontsize=12)

        model_text = 'Test Model: ' + test_name
        ax.text(0.6, 1, model_text, ha='left', va='center', transform=ax.transAxes,color=color[0], fontsize=12)
        if run_type == 'model_vs_model':
            ax.text(0.6, 0.95, 'Ref. Model: ' + ref_name, ha='left', va='center', transform=ax.transAxes, color='k', fontsize=12)

        plt.title(season + ': Spatial Variability', y = 1.08)
        fig.savefig(os.path.join(output_dir, season + '_metrics_taylor_diag.png'))

    return taylor_diag_path

def _cvs_to_html(csv_path, season, test_name, ref_name):
    """Convert the csv for a season located at csv_path to an HTML, returning the path to the HTML"""
    html_path = csv_path.replace('csv', 'html')

    with open(html_path, 'w') as htmlfile:
        htmlfile.write('<p><b>Test: {}</b><br>'.format(test_name))
        htmlfile.write('<b>Reference: {}</b></p>'.format(ref_name))
        htmlfile.write('<p><th><b>{} Mean </b></th></p>'.format(season))
        htmlfile.write('<table>')

        with open(csv_path) as csv_file:
            read_csv = csv.reader(csv_file)

            # generate table contents
            for num, row in enumerate(read_csv):

                # write the header row, assuming the first row in csv contains the header
                if num == 0:
                    htmlfile.write('<tr>')
                    for column in row:
                        htmlfile.write('<th>{}</th>'.format(column))
                    htmlfile.write('</tr>')

                # write all other rows 
                else:
                    htmlfile.write('<tr><div style="width: 50px">')
                    for column in row:
                        htmlfile.write('<td>{}</td>'.format(column))
                    htmlfile.write('</div></tr>')

        htmlfile.write('</table>')

    return html_path

def _add_html_to_col(season, season_path, html_path):
    """Since the output viewer doesn't support html images, do this hack.
    For the col in the html at html_path, insert the link to col_path."""
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
        index = ['All variables', 'ANN', 'DJF', 'MAM', 'JJA', 'SON'].index(season)
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

def _edit_table_html(root_dir):
    """After the viewer is created, edit the table html to insert the custom htmls"""
    for s in ['ANN', 'DJF', 'MAM', 'JJA', 'SON']:
        if s in LAT_LON_TABLE_INFO:
            _add_html_to_col(s, LAT_LON_TABLE_INFO[s]['html_path'], os.path.join(root_dir, 'table', 'index.html'))

def _create_lat_lon_table_index(viewer, root_dir):
    """Create an index in the viewer that links the individual htmls for the lat-lon table."""
    seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']
    viewer.add_page('Table', seasons)
    viewer.add_group('Summary Table')
    viewer.add_row('All variables')

    for s in seasons:
        if s in LAT_LON_TABLE_INFO:
            viewer.add_col(LAT_LON_TABLE_INFO[s]['html_path'], is_file=True, title=s)
        else:
            viewer.add_col('-----', is_file=True, title='-----')

def _create_taylor_index(viewer, root_dir, season_to_png):
    """Create an index in the viewer that links the individual htmls for the lat-lon table."""
    seasons = ['ANN', 'DJF', 'MAM', 'JJA', 'SON']
    viewer.add_page('Taylor Diagram', seasons)
    viewer.add_group('Summary Taylor Diagrams')
    viewer.add_row('All variables')

    for s in seasons:
        if s in season_to_png:
            pth = os.path.join('..', season_to_png[s])
            viewer.add_col(pth, is_file=True, title=s)
        else:
            viewer.add_col('-----', is_file=True, title='-----')

def _add_table_and_taylor_to_viewer_index(root_dir):
    """Move the link to 'Table' and 'Taylor Diagram' next to the link to Latitude-Longitude contour maps"""
    index_page = os.path.join(root_dir, 'index.html')
    soup = BeautifulSoup(open(index_page), "lxml")

    # append the new tag underneath the old one, so add it to the parent of the old one
    td_to_move_table = None
    td_to_move_taylor = None

    # rows to delete
    delete_this_table = None
    delete_this_taylor = None

    for tr in soup.find_all("tr"):
        for td in tr.find_all("td"):
            for a in td.find_all("a"):
                if 'table' in a['href']:
                    td_to_move_table = copy.deepcopy(td)
                    delete_this_table = tr
                if 'taylor' in a['href']:
                    td_to_move_taylor = copy.deepcopy(td)
                    delete_this_taylor = tr

    if delete_this_table:
        delete_this_table.decompose()
    if delete_this_taylor:
        delete_this_taylor.decompose()

    for tr in soup.find_all("tr"):
        for td in tr.find_all("td"):
            for a in td.find_all("a"):
                if _better_page_name('lat_lon') in a.string and td_to_move_table:
                    td.append(td_to_move_table)
                if _better_page_name('lat_lon') in a.string and td_to_move_taylor:
                    td.append(td_to_move_taylor)    

    html = soup.prettify("utf-8")
    with open(index_page, "wb") as f:
        f.write(html)

def _add_provenance_to_viewer_index(root_dir):
    """
    In the index, link 'Provenance' to the prov folder.
    """
    index_page = os.path.join(root_dir, 'index.html')
    soup = BeautifulSoup(open(index_page), "lxml")

    for tr in soup.find_all("tr"):
        for td in tr.find_all("td"):
            for a in td.find_all("a"):
                if 'Provenance' in a.text:
                    a['href'] = '../prov'

    html = soup.prettify("utf-8")
    with open(index_page, "wb") as f:
        f.write(html)

def generate_lat_lon_metrics_table(viewer, root_dir, parameters):
    """For each season in LAT_LON_TABLE_INFO, create a csv, convert it to an html and append that html to the viewer."""
    table_dir = os.path.join(root_dir, 'table-data')  # output_dir/viewer/table-data

    if not os.path.exists(table_dir):
        os.mkdir(table_dir)

    for season in LAT_LON_TABLE_INFO:
        test_name = parameters[0].short_test_name if parameters[0].short_test_name else parameters[0].test_name
        if parameters[0].run_type == 'model_vs_obs':
            ref_name = 'Observation and Reanalysis'
        else:
            ref_name = parameters[0].short_ref_name if parameters[0].short_ref_name else parameters[0].ref_name            
        csv_path = _create_csv_from_dict(table_dir, season, test_name, parameters[0].run_type)
        html_path = _cvs_to_html(csv_path, season, test_name, ref_name)

        # Ex: change this: /Users/zshaheen/output_dir/viewer/table-data/ANN_metrics_table.html
        # to this: viewer/table-data/ANN_metrics_table.html
        html_path = '/'.join(html_path.split('/')[-3:])

        LAT_LON_TABLE_INFO[season]['html_path'] = html_path

    _create_lat_lon_table_index(viewer, root_dir)

def generate_lat_lon_taylor_diag(viewer, root_dir, parameters):
    """For each season in LAT_LON_TABLE_INFO, create a csv, plot using taylor diagram  and append that html to the viewer."""
    taylor_diag_dir = os.path.join(root_dir, 'taylor-diagram-data')  # output_dir/viewer/taylor-diagram-data

    if not os.path.exists(taylor_diag_dir):
        os.mkdir(taylor_diag_dir)

    season_to_png = {}
    for season in LAT_LON_TABLE_INFO:
        test_name = parameters[0].short_test_name if parameters[0].short_test_name else parameters[0].test_name
        if parameters[0].run_type == 'model_vs_obs':
            ref_name = 'Observation and Reanalysis'
        else:
            ref_name = parameters[0].short_ref_name if parameters[0].short_ref_name else parameters[0].ref_name            

        csv_path = _create_csv_from_dict_taylor_diag(taylor_diag_dir, season, test_name, parameters[0].run_type, ref_name)
        # Remove any reference to the results_dir when inserting the links into HTML pages.
        # This is because that folder can be renamed.
        csv_path = csv_path.split('viewer')[-1]
        csv_path = 'viewer' + csv_path

        season_to_png[season] = csv_path.replace('csv', 'png')

    _create_taylor_index(viewer, root_dir, season_to_png)

def create_metadata(parameter):
    """
    From a set of parameters, extract the metadata.
    """
    metadata = collections.OrderedDict()
    msg = 'Use this command to recreate this image:'
    metadata[msg] = ''
    cmd = 'e3sm_diags --no_viewer '

    from acme_diags.acme_parser import ACMEParser
    parser = ACMEParser()

    args = parser.view_args()
    supported_cmd_args = list(args.__dict__.keys())
    
    if 'other_parameters' in supported_cmd_args:
        supported_cmd_args.remove('other_parameters')
    
    if 'parameters' in supported_cmd_args:
        supported_cmd_args.remove('parameters')

    if container.is_container():
        container.decontainerize_parameter(parameter)

    for param_name in parameter.__dict__:
        param = parameter.__dict__[param_name]
        # we don't want to include blank values
        if not param:
            continue

        if param_name in supported_cmd_args:
            if isinstance(param, list) or isinstance(param, tuple):
                # ex: --diff_levels -7, -6, -5, -4
                cmd += "--{} ".format(param_name)
                for p in param:
                    if isinstance(p, str) and p.isdigit():
                        cmd += " {} ".format(str(p))
                    else:
                        cmd += " '{}' ".format(str(p))
            
            elif isinstance(param, bool):
                # ex: --multiprocessing
                # note there's no value after the parameter, it's just a flag
                if param:  # command is True, so add --command to set it to True
                    cmd += "--{} ".format(param_name)

            elif isinstance(param, str) and param.isdigit():
                cmd += "--{} {} ".format(param_name, param)
            else:
                cmd += "--{} '{}' ".format(param_name, param)
    
    metadata[msg] = cmd

    return metadata

def create_viewer(root_dir, parameters, ext):
    """Based of the parameters, find the files with
    extension ext and create the viewer in root_dir."""

    viewer = OutputViewer(path=root_dir, index_name='E3SM Diagnostics')
    _add_pages_and_top_row(viewer, parameters)

    for parameter in parameters:
        results_dir = parameter.results_dir
        for set_num in parameter.sets:
            set_num = get_set_name(set_num)

            # Filenames are like:
            # ref_name-variable-season-region
            # or
            # ref_name-variable-plev'mb'-season-region
            ref_name = getattr(parameter, 'ref_name', '')
            for var in parameter.variables:
                for season in parameter.seasons:
                    for region in parameter.regions:
                        # because some parameters have plevs, there might be
                        # more than one row_name, fnm pair
                        row_name_and_fnm = []

                        if parameter.plevs == []:  # 2d variables
                            if parameter.run_type == 'model_vs_model':
                                row_name = '{} {}'.format(var, region)
                            else: 
                                row_name = '{} {} {}'.format(var, region, ref_name)
                            fnm = '{}-{}-{}-{}'.format(ref_name,
                                                       var, season, region)
                            row_name_and_fnm.append((row_name, fnm))
                        else:  # 3d variables
                            for plev in parameter.plevs:
                                if parameter.run_type == 'model_vs_model':
                                    row_name = '{}-{} {}'.format(
                                        var, str(int(plev)) + 'mb', region)
                                else:
                                    row_name = '{}-{} {} {}'.format(
                                        var, str(int(plev)) + 'mb', region, ref_name)
                                fnm = '{}-{}-{}-{}-{}'.format(
                                    ref_name, var, int(plev), season, region)
                                row_name_and_fnm.append((row_name, fnm))

                        if set_num in ['lat_lon', '5']:
                            metrics_path = os.path.join(results_dir, '{}'.format(set_num), parameter.case_id, fnm)
                            if os.path.exists(metrics_path + '.json'):
                                _add_to_lat_lon_metrics_table(metrics_path, season, row_name)
                            else:
                                print(('JSON does not exist: {}'.format(metrics_path + '.json')))
                                continue
                        for row_name, fnm in row_name_and_fnm:
                            if parameter.case_id not in ROW_INFO[set_num]:
                                ROW_INFO[set_num][parameter.case_id] = collections.OrderedDict(
                                )
                            if row_name not in ROW_INFO[set_num][parameter.case_id]:
                                ROW_INFO[set_num][parameter.case_id][row_name] = collections.OrderedDict(
                                )
                                ROW_INFO[set_num][parameter.case_id][row_name]['descr'] = _get_description(
                                    var, parameter)
                            # each season has a image_path and metadata linked to it, thus we use a dict
                            ROW_INFO[set_num][parameter.case_id][row_name][season] = {}
                            # format fnm to support relative paths
                            ROW_INFO[set_num][parameter.case_id][row_name][season]['image_path'] = os.path.join(
                                '..', '{}'.format(set_num), parameter.case_id, fnm)
                            # If ran in a container, create_metadata() will modify *_data_path and results_dir to their original value.
                            ROW_INFO[set_num][parameter.case_id][row_name][season]['metadata'] = create_metadata(parameter)

    # add all of the files in from the case_id/ folder in ANN, DJF, MAM, JJA,
    # SON order
    for set_num in ROW_INFO:
        viewer.set_page("{}".format(_better_page_name(set_num)))

        for group in ROW_INFO[set_num]:
            try:
                viewer.set_group(group)
            except RuntimeError:
                viewer.add_group(group)

            for row_name in ROW_INFO[set_num][group]:
                try:
                    viewer.set_row(row_name)
                except RuntimeError:
                    # row of row_name wasn't in the viewer, so add it
                    viewer.add_row(row_name)
                    # the description, currently the var
                    viewer.add_col(ROW_INFO[set_num][group][row_name]['descr'])

                # [1:] is to ignore 'Description' col
                for col_season in viewer.page.columns[1:]:
                    if col_season not in ROW_INFO[set_num][group][row_name]:
                        viewer.add_col('-----', is_file=True, title='-----')
                    else:
                        metadata = ROW_INFO[set_num][group][row_name][col_season]['metadata']
                        fnm = ROW_INFO[set_num][group][row_name][col_season]['image_path']
                        formatted_files = []
                        if parameters[0].save_netcdf:
                            nc_files = [
                                fnm + nc_ext for nc_ext in ['_test.nc', '_ref.nc', '_diff.nc']]
                            formatted_files = [
                                {'url': f, 'title': f} for f in nc_files]
                        viewer.add_col(fnm + '.' + ext, is_file=True,
                                       title=col_season, other_files=formatted_files, meta=metadata)

    generate_lat_lon_metrics_table(viewer, root_dir, parameters)
    generate_lat_lon_taylor_diag(viewer, root_dir, parameters)
    viewer.generate_viewer(prompt_user=False)
    _extras(root_dir, parameters)

