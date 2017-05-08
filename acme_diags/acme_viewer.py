import os
import glob
from cdp.cdp_viewer import OutputViewer
from acme_diags.driver.utils import get_output_dir

# Dict of {set_num: {row_name: {season: filename}}}, follows the order
# of pages in the viewer.
# Used when actually inserting the cols into the viewer
# needed so we can have a cols in order of ANN, DJF, MAM, JJA, SON
ROW_INFO = {}

def add_pages_and_top_row(viewer, parameters):
    """Add the page and columns of the page"""
    set_to_seasons = {}  # dict of {set: [seasons]}
    for p in parameters:
        for set_num in p.sets:
            for ssn in p.seasons:
                if set_num not in set_to_seasons:
                    set_to_seasons[set_num] = []
                set_to_seasons[set_num].append(ssn)

    for set_num, seasons in set_to_seasons.iteritems():
        ROW_INFO[set_num] = {}
        col_labels = ['Description']
        for s in ['ANN', 'DJF', 'MAM', 'JJA', 'SON']:
            if s in seasons:
                col_labels.append(s)
        viewer.add_page("Set-{}".format(set_num), col_labels)

def create_viewer(root_dir, parameters, ext):
    """Based of the parameters, find the files with
    extension ext and create the viewer in root_dir."""

    viewer = OutputViewer(path=root_dir, index_name='ACME Diagnostics')
    add_pages_and_top_row(viewer, parameters)


    for parameter in parameters:
        for set_num in parameter.sets:
            viewer.set_page("Set-{}".format(set_num))
            viewer.add_group(parameter.case_id)

            # Add all of the files with extension ext from the case_id/ folder'
            pth = get_output_dir(set_num, parameter)
            files = glob.glob(pth + '/*.' + ext)  # ex: files[0] = myresults/set5/set5_SST_HadISST/HadISST_CL-SST-SON-global.png

            for ext_fnm in files:
                fnm = ext_fnm.replace('.' + ext, '')
                fnm = fnm.split('/')[-1]  # ex: HadISST_CL-SST-SON-global
                keywords = fnm.split('-')

                # 2d vars, format is [ref_name, var, season, region]
                # 3d vars, format is [ref_name, var, plev, season, region]
                # ref_name and/or var can be something_like_this, so we always use negative indices
                region = keywords[-1]
                season = keywords[-2]
                if keywords[-3].isdigit():  # for when we have a 3d variable, ex: keywords[-3] = 880
                    plev = keywords[-3]
                    var = keywords[-4]
                else:
                    plev = None
                    var = keywords[-3]

                if plev is None:  # 2d variable
                    row_name = '%s %s' % (var, region)
                else:  # 3d variable
                    row_name = '%s %s %s' % (var, plev + ' mb ', region)

                try:
                    viewer.set_row(row_name)
                except RuntimeError:
                    # row of row_name wasn't in the viewer, so add it
                    viewer.add_row(row_name)
                    viewer.add_col(var)  # the description

                if row_name not in ROW_INFO[set_num]:
                    ROW_INFO[set_num][row_name] = {}
                # format fnm to support relative paths
                ROW_INFO[set_num][row_name][season] = os.path.join('set{}'.format(set_num), parameter.case_id, fnm)
    
    # add all of the files in from the case_id/ folder in ANN, DJF, MAM, JJA, SON order
    for set_num in ROW_INFO:
        viewer.set_page("Set-{}".format(set_num))
        for row_name in ROW_INFO[set_num]:
            for col_season in viewer.page.columns[1:]:  # [1:] is to ignore 'Description' col 
                viewer.set_row(row_name)
                if col_season in ROW_INFO[set_num][row_name]:
                    fnm = ROW_INFO[set_num][row_name][col_season]
                    nc_files = [fnm + nc_ext for nc_ext in ['_test.nc', '_ref.nc', '_diff.nc']]
                    formatted_files = [{'url': f, 'title': f} for f in nc_files]
                    viewer.add_col(fnm + '.' + ext, is_file=True, title=col_season, other_files=formatted_files)
                else:
                    # insert a blank value
                    viewer.add_col('-----')

    viewer.generate_viewer()
