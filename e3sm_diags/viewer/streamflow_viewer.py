import os

from cdp.cdp_viewer import OutputViewer

from .default_viewer import create_metadata
from .utils import add_header, h1_to_h3


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for the streamflow set,
    create a single webpage.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)

    # The name that's displayed on the viewer.
    display_name = "Streamflow"
    set_name = "streamflow"
    # The title of the columns on the webpage.
    # Appears in the second and third columns of the bolded rows.
    cols = ["Description", "Plot"]
    viewer.add_page(display_name, short_name=set_name, columns=cols)
    for plot_type in ["seasonality_map", "annual_map", "annual_scatter"]:
        if plot_type == "seasonality_map":
            group_name = "Seasonality Map"
        elif plot_type == "annual_map":
            group_name = "Mean Annual Streamflow Map"
        elif plot_type == "annual_scatter":
            group_name = "Mean Annual Streamflow Scatter Plot"
        else:
            raise Exception("Invalid plot_type={}".format(plot_type))
        # Appears in the first column of the bolded rows.
        viewer.add_group(group_name)
        for param in parameters:
            # We need to make sure we have relative paths, and not absolute ones.
            # This is why we don't use get_output_dir() as in the plotting script
            # to get the file name.
            ext = param.output_format[0]
            # TODO: will param.variables ever be longer than one variable?
            #  If so, we'll need to create unique image_relative_paths
            for var in param.variables:
                param.var_id = "{}-{}".format(var, plot_type)
                if plot_type == "seasonality_map":
                    param.viewer_descr[var] = param.main_title_seasonality_map
                    output_file = param.output_file_seasonality_map
                elif plot_type == "annual_map":
                    param.viewer_descr[var] = param.main_title_annual_map
                    output_file = param.output_file_annual_map
                elif plot_type == "annual_scatter":
                    param.viewer_descr[var] = param.main_title_annual_scatter
                    output_file = param.output_file_annual_scatter
                else:
                    raise Exception("Invalid plot_type={}".format(plot_type))
                # param.output_file_{seasonality, annual_map, annual_scatter} is defined in
                # e3sm_diags/parameter/streamflow_parameter.py.
                # relative_path must use param.case_id and output_file
                # to match the file_path determined in
                # e3sm_diags/plot/cartopy/streamflow_plot.py.
                # Otherwise, the plot will not be properly linked from the viewer.
                relative_path = os.path.join("..", set_name, param.case_id, output_file)
                image_relative_path = "{}.{}".format(relative_path, ext)
                if param.print_statements:
                    print("image_relative_path: {}".format(image_relative_path))
                # Appears in the first column of the non-bolded rows.
                # To match the output_dir determined by get_output_dir in e3sm_diags/plot/cartopy/streamflow_plot.py,
                # the row name should be param.case_id.
                # Otherwise, the plot image and the plot HTML file will have URLs with different final directory names.
                # (The viewer will still display everything properly though).
                viewer.add_row("{}-{}".format(param.case_id, plot_type))
                # Adding the description for this var to the current row.
                # Appears in the second column of the non-bolded rows.
                viewer.add_col(param.viewer_descr[var])
                # Link to an html version of the plot png file.
                # Appears in the third column of the non-bolded rows.
                meta = create_metadata(param)
                viewer.add_col(
                    image_relative_path,
                    is_file=True,
                    title="Plot",
                    other_files=[],
                    meta=meta,
                )

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
