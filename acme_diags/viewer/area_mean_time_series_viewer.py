import os
from .default_viewer import create_metadata
from cdp.cdp_viewer import OutputViewer


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for a certain set of diagnostics,
    create a single page.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)
    viewer.add_group('Variable')

    # The name that's displayed on the viewer.
    display_name = 'Area Mean Time Series'
    set_name = 'area_mean_time_series'
    cols = ['Description', 'Plot']
    viewer.add_page(display_name, short_name=set_name, columns=cols)

    for param in parameters:
        for var in param.variables:
            viewer.add_row(var)
            # Adding the description for this var to the current row.
            # This was obtained and stored in the driver for this plotset.
            viewer.add_col(param.viewer_descr[var])

            ext = param.output_format[0]
            file_name = os.path.join('..', set_name, '{}.{}'.format(var, ext))
            viewer.add_col(file_name, is_file=True, title='Plot')

