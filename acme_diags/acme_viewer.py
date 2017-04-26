import os
import glob
from cdp.cdp_viewer import OutputViewer


def add_page_and_top_row(viewer, parameters):
    ''' Setup for OutputViewer '''
    col_labels = ['Description']
    seasons = []
    for p in parameters:
        for s in p.season:
            if s not in seasons:
                seasons.append(s)
    for s in seasons:
        col_labels.append(s)

    viewer.add_page("Set 5", col_labels)

def create_viewer(root_dir, parameters):
    ''' Based of the parameters, find
    the files and create the viewer '''

    viewer = OutputViewer(path=root_dir, index_name='ACME Diagnostics')
    add_page_and_top_row(viewer, parameters)

    for parameter in parameters:
        viewer.add_group(parameter.case_id)

        # Add all of the .png files from the case_id/ folder'
        pth = os.path.join(parameter.results_dir, parameter.case_id)
        files = glob.glob(pth + '/*.png')
        for png_file in files:
            fnm = png_file.replace('.png', '')
            keywords = fnm.split('_')

            if len(keywords) == 5:
                # 2d vars, format is [set_num, ref_name, var, season, region]
                row_name = '%s %s' % (keywords[2], keywords[4])
            elif len(keywords) == 6:
                # 3d vars, format is [set_num, ref_name, var, plev, season, region]
                row_name = '%s %s %s' % (keywords[2], keywords[3] + ' mb ', keywords[5])
            else:
                raise RuntimeError('Invalid name of plot files')

            try:
                viewer.set_row(row_name)
            except RuntimeError:
                # row of row_name wasn't added, so add it
                viewer.add_row(row_name)
                viewer.add_col(keywords[2])

            nc_files = [fnm + ext for ext in ['_test.nc', '_ref.nc', '_diff.nc']]
            formatted_files = [{'url': f, 'title': f} for f in nc_files]
            viewer.add_col(os.path.join( '..', png_file), is_file=True, title=keywords[-2], other_files=formatted_files)
    viewer.generate_viewer()
