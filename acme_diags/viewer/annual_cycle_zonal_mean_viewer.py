import os

from cdp.cdp_viewer import OutputViewer

from .default_viewer import create_metadata
from .utils import add_header, h1_to_h3


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for a the diff_diags set, create a single webpage.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)
    display_name = "Annual Cycle Zonal Mean Contour Plots"
    set_name = "annual_cycle_zonal_mean"
    cols = ["Description", "Plot"]

    viewer.add_page(display_name, short_name=set_name, columns=cols)
    viewer.add_group("Variable")

    for param in parameters:
        for var in param.variables:
            viewer.add_row("{} {}".format(var, param.ref_name))
            # Adding the description for this var to the current row.
            # This was obtained and stored in the driver for this plotset.
            viewer.add_col(param.viewer_descr[var])

            file_name = "{}-{}-{}.png".format(param.ref_name, var, "Annual-Cycle")
            # We need to make sure we have relative paths, and not absolute ones.
            # This is why we don't use get_output_dir() as in the plotting script
            # to get the file name.
            file_name = os.path.join("..", set_name, param.case_id, file_name)
            viewer.add_col(
                file_name, is_file=True, title="Plot", meta=create_metadata(param)
            )

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
