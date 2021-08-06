import os

from cdp.cdp_viewer import OutputViewer

from .utils import add_header, h1_to_h3


def create_viewer(root_dir, parameters):
    """
    Given a set of parameter for a certain set of diagnostics,
    create a single page.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)
    # The name that's displayed on the viewer.
    set_name = "tc_analysis"
    display_name = "Diagnostics for Tropical Cyclones"
    # The title of the colums on the webpage.
    # Appears in the second and third columns of the bolded rows.
    cols = ["Description", "Plot"]

    viewer.add_page(display_name, short_name=set_name, columns=cols)

    relative_path = os.path.join("..", set_name)
    ext = parameters[0].output_format[0]

    viewer.add_group("Bar/Line plots")
    viewer.add_row("TC Frequency")
    output_file = "tc-frequency.{}".format(ext)
    image_relative_path = os.path.join(relative_path, output_file)
    viewer.add_col("Frequency of TCs by ocean basins")
    viewer.add_col(image_relative_path, is_file=True, title="Plot")

    viewer.add_row("TC Intensity")
    output_file = "tc-intensity.{}".format(ext)
    image_relative_path = os.path.join(relative_path, output_file)
    viewer.add_col("Intensity of TCs by ocean basins")
    viewer.add_col(image_relative_path, is_file=True, title="Plot")

    viewer.add_row("TC Seasonality")
    output_file = "tc-frequency-annual-cycle.{}".format(ext)
    image_relative_path = os.path.join(relative_path, output_file)
    viewer.add_col("Annual cycle of TC frequency by ocean basins")
    viewer.add_col(image_relative_path, is_file=True, title="Plot")

    viewer.add_row("ACE Distribution")
    output_file = "ace-distribution.{}".format(ext)
    image_relative_path = os.path.join(relative_path, output_file)
    viewer.add_col("Distribution of Accumulated Cyclone Energy (ACE)")
    viewer.add_col(image_relative_path, is_file=True, title="Plot")

    viewer.add_group("Maps")
    viewer.add_row("Cyclone Density Map")
    output_file = "cyclone-density-map.{}".format(ext)
    image_relative_path = os.path.join(relative_path, output_file)
    viewer.add_col("Density map of Tropical Cyclones")
    viewer.add_col(image_relative_path, is_file=True, title="Plot")

    viewer.add_row("AEW Density Map")
    output_file = "aew-density-map.{}".format(ext)
    image_relative_path = os.path.join(relative_path, output_file)
    viewer.add_col("Density map of African Easterly Wave")
    viewer.add_col(image_relative_path, is_file=True, title="Plot")

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
