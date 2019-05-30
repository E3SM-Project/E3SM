import collections
from . import default_viewer

# A mapping of each diagnostics set to the viewer
# that handles creating of the HTML pages.
SET_TO_VIEWER = {
    'lat_lon': default_viewer.create_viewer,
    'polar': default_viewer.create_viewer,
    'zonal_mean_xy': default_viewer.create_viewer,
    'zonal_mean_2d': default_viewer.create_viewer,
    'meridional_mean_2d': default_viewer.create_viewer,
    'cosp_histogram': default_viewer.create_viewer
}

def create_index(root_dir, title_and_url_list):
    """
    Creates the index which joins the individual viewers.
    A list of (title, url) tuples are passed
    in to create the index.
    """
    pass


def create_viewer(root_dir, parameters, extension='png'):
    """
    Based of the parameters, find the files with the
    certain extension and create the viewer in root_dir.
    """
    
    # Group each parameter object based on the `sets` parameter.
    set_to_parameters = collections.defaultdict(list)
    for param in parameters:
        for set_name in param.sets:
            set_to_parameters.append(param)
    
    print(set_to_parameters)

    # A list of (title, url) tuples that each viewer generates.
    # This is used to create the main index.
    title_and_url_list = []
    # Now call the viewers with the parameters as the arguments.
    for set_name, parameters in set_to_parameters.items():
        viewer_function = SET_TO_VIEWER[set_name]
        title, url = viewer_function(parameters)
        title_and_url_list.append((title, url))
    
    index_url = create_index(title_and_url_list)

    return index_url
