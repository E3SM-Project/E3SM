import os
import shutil

from e3sm_diags.viewer.core_viewer import OutputViewer

from .default_viewer import seasons_used
from .lat_lon_viewer import _cvs_to_html
from .utils import _fix_table_col_links, add_header, h1_to_h3


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for a the diff_diags set,
    create a single webpage.

    Return the title and url for this page.
    """
    # Viewer configurations.
    display_name = "Aerosol Budget Tables"
    set_name = "aerosol_budget"
    seasons = seasons_used(parameters)

    test_name = (
        parameters[0].short_test_name
        if parameters[0].short_test_name
        else parameters[0].test_name
    )
    ref_name = (
        parameters[0].short_ref_name
        if parameters[0].short_ref_name
        else parameters[0].ref_name
    )

    # Directory for the CSVs and Tables.
    # Example: aerosol_budgets/Aerosol_table/aerosol_budget/
    csv_table_dir = os.path.join(root_dir, "..", f"{set_name}")

    # Table data directory -- stores data as CSV and HTML table files, which
    # populate the viewer's main table.
    # Example: aerosol_budgets/viewer/table
    table_dir = os.path.join(root_dir, "table")
    if not os.path.exists(table_dir):
        os.mkdir(table_dir)

    # Create the Viewer object and configure.
    viewer = OutputViewer(path=root_dir)
    viewer.add_page(display_name, columns=seasons)
    viewer.add_group("All Species")
    viewer.add_row("All Species")

    # Populate each column with a link to each HTML table file by season.
    season_html_paths = {}
    for season in seasons:
        # Copy the CSV file to the table-data viewer directory.
        # /aerosol_budgets/Aerosol_table/viewer/../aerosol_budget/v2.LR.historical_0101-ANN-budget-table.csv
        csv_filename = f"{parameters[0].test_name}-{season}-budget-table.csv"
        csv_path = os.path.join(csv_table_dir, csv_filename)
        shutil.copy(csv_path, table_dir)

        # Example: /aerosol_budgets/Aerosol_table/viewer/table/v2.LR.historical_0101-ANN-budget-table.csv
        csv_copy_path = os.path.join(table_dir, csv_filename)

        # Convert the CSV to an HTML table.
        # Example: aerosol_budgets/Aerosol_table/viewer/table/v2.LR.historical_0101-ANN-budget-table.html
        html_table_path = _cvs_to_html(csv_copy_path, season, test_name, ref_name)

        # Create a column using the relative path to the HTML table.
        # Example: ../viewer/table-data/v2.LR.historical_0101-ANN-budget-table.html
        rel_html_path = os.path.join("..", "/".join(html_table_path.split("/")[-3:]))
        viewer.add_col(rel_html_path, is_file=True, title=season)

        season_html_paths[season] = rel_html_path

    url = viewer.generate_page()

    # Example: aerosol_budgets/Aerosol_table/viewer/aerosol/index.html
    index_path = f"{root_dir}/{url}"
    _fix_table_col_links(index_path, "All Species", season_html_paths)

    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
