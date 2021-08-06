import os

from cdp.cdp_viewer import OutputViewer

from .utils import add_header, h1_to_h3

region_name = {
    "sgp": "Southern Great Plains",
    "nsa": "North Slope Alaska",
    "twpc1": "Manus",
    "twpc2": "Nauru",
    "twpc3": "Darwin",
}


def create_viewer(root_dir, parameters):
    """
    Given a set of parameter for a certain set of diagnostics,
    create a single page.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)
    # The name that's displayed on the viewer.
    set_name = "arm_diags"
    display_name = "Diagnostics at ARM stations"
    # The title of the colums on the webpage.
    # Appears in the second and third columns of the bolded rows.
    cols = ["Description", "Plot", "", "", "", ""]

    viewer.add_page(display_name, short_name=set_name, columns=cols)
    param_dict = {}  # type: ignore

    relative_path = os.path.join("..", set_name)
    for param in parameters:
        key = param.diags_set
        if key not in param_dict.keys():
            param_dict[key] = []
        param_dict[key].append(param)

    for diags_set in param_dict.keys():
        valid_parameters = param_dict[diags_set]
        if diags_set == "annual_cycle":
            viewer.add_group("Annual Cycle")
            for param in valid_parameters:
                ext = param.output_format[0]
                viewer.add_row(
                    "{} at {} ({})".format(
                        param.variables[0],
                        region_name[param.regions[0]],
                        param.regions[0],
                    )
                )
                viewer.add_col("Annual cycles of " + param.var_name)
                image_relative_path = os.path.join(
                    relative_path, "{}.{}".format(param.output_file, ext)
                )
                viewer.add_col(image_relative_path, is_file=True, title="Plot")
        if diags_set == "diurnal_cycle":
            viewer.add_group("Diurnal Cycle")
            for param in valid_parameters:
                ext = param.output_format[0]
                viewer.add_row("PRECT at SGP")
                viewer.add_col("Diurnal cycle of precipitation")
                for season in ["DJF", "MAM", "JJA", "SON"]:
                    output_file = "{}-PRECT-{}-sgp-diurnal-cycle.{}".format(
                        param.ref_name, season, ext
                    )
                    image_relative_path = os.path.join(relative_path, output_file)
                    viewer.add_col(image_relative_path, is_file=True, title=season)
        if diags_set == "convection_onset":
            viewer.add_group("Convection Onset Statistics")
            for param in valid_parameters:
                ext = param.output_format[0]
                viewer.add_row(
                    "{} ({})".format(region_name[param.regions[0]], param.regions[0])
                )
                viewer.add_col("Convection Onset Statistics")
                region = param.regions[0]
                output_file = "{}-convection-onset-{}.{}".format(
                    param.ref_name, region, ext
                )
                image_relative_path = os.path.join(relative_path, output_file)
                viewer.add_col(image_relative_path, is_file=True, title="Plot")
        if diags_set == "diurnal_cycle_zt":
            viewer.add_group("Monthly Diurnal Cycle of Cloud")
            for param in valid_parameters:
                ext = param.output_format[0]
                viewer.add_row(
                    "{} ({})".format(region_name[param.regions[0]], param.regions[0])
                )
                region = param.regions[0]
                viewer.add_col("Monthly Diurnal Cycle of Cloud")
                output_file1 = "{}-CLOUD-ANNUALCYCLE-{}-{}.{}".format(
                    param.ref_name, region, "test", ext
                )
                output_file2 = "{}-CLOUD-ANNUALCYCLE-{}-{}.{}".format(
                    param.ref_name, region, "ref", ext
                )
                image_relative_path = os.path.join(relative_path, output_file1)
                viewer.add_col(image_relative_path, is_file=True, title="Test")
                image_relative_path = os.path.join(relative_path, output_file2)
                viewer.add_col(image_relative_path, is_file=True, title="Reference")

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
