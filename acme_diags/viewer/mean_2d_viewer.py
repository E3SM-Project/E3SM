import os
from .utils import add_header, h1_to_h3
from .default_viewer import create_metadata, seasons_used, SEASONS
from cdp.cdp_viewer import OutputViewer


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for a certain set of diagnostics,
    create a single page.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)

    set_name = parameters[0].sets[0]
    # The name that's displayed on the viewer.
    if set_name == 'zonal_mean_2d':
        display_name = 'Pressure-Latitude zonal mean contour plots'
    elif set_name == 'meridional_mean_2d':
        display_name = 'Pressure-Longitude meridional mean contour plots'

    cols = ['Description'] + seasons_used(parameters)
    viewer.add_page(display_name, short_name=set_name, columns=cols)

    # Sort the parameters so that the viewer is created in the correct order.
    # Using SEASONS.index(), we make sure we get the parameters in
    # ['ANN', 'DJF', ..., 'SON'] order instead of alphabetical.
    parameters.sort(key=lambda x: (x.case_id, x.variables[0], SEASONS.index(x.seasons[0])))

    for param in parameters:
        ref_name = getattr(param, 'ref_name', '')
        for var in param.variables:
            for season in param.seasons:
                for region in param.regions:

                    try:
                        viewer.set_group(param.case_id)
                    except RuntimeError:
                        viewer.add_group(param.case_id)

                    if param.run_type == 'model_vs_model':
                        row_name = '{} {}'.format(var, region)
                    else: 
                        row_name = '{} {} {}'.format(var, region, ref_name)

                    try:
                        viewer.set_row(row_name)
                    except RuntimeError:
                        # A row of row_name wasn't in the viewer, so add it.
                        viewer.add_row(row_name)
                        # Adding the description for the row.
                        viewer.add_col(param.viewer_descr[var])

                    ext = param.output_format[0]
                    fnm = '{}-{}-{}-{}'.format(ref_name, var, season, region)
                    file_name = os.path.join('..', set_name, param.case_id, '{}.{}'.format(fnm, ext))
                    viewer.add_col(file_name, is_file=True, title=season,
                        meta=create_metadata(param))

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
