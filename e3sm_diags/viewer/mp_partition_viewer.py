import os

from cdp.cdp_viewer import OutputViewer

from e3sm_diags.logger import custom_logger

from .default_viewer import create_metadata
from .utils import add_header, h1_to_h3

logger = custom_logger(__name__)


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for the mixed phase partition (mp_partition) set,
    create a single webpage.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)

    # The name that's displayed on the viewer.
    display_name = "Mixed-phase Partition"
    set_name = "mp_partition"
    # The title of the colums on the webpage.
    # Appears in the second and third columns of the bolded rows.
    cols = ["Description", "Plot"]
    viewer.add_page(display_name, short_name=set_name, columns=cols)
    # Appears in the first column of the bolded rows.
    viewer.add_group("Variable")
    for param in parameters:
        # We need to make sure we have relative paths, and not absolute ones.
        # This is why we don't use get_output_dir() as in the plotting script
        # to get the file name.
        ext = param.output_format[0]
        relative_path = os.path.join("..", set_name, param.case_id, param.output_file)
        image_relative_path = "{}.{}".format(relative_path, ext)
        if param.print_statements:
            logger.info("image_relative_path: {}".format(image_relative_path))
        for var in param.variables:
            viewer.add_row("LCF")
            # Adding the description for this var to the current row.
            # This was obtained and stored in the driver for this plotset.
            # Appears in the second column of the non-bolded rows.
            viewer.add_col("Liquid Condensate Fraction ")
            # Link to an html version of the plot png file.
            # Appears in the third column of the non-bolded rows.
            viewer.add_col(
                image_relative_path,
                is_file=True,
                title="Plot",
                meta=create_metadata(param),
            )

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
