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
    
    # A map of the viewer to a list of parameters that are
    # needed to be processed to that viewer.
    # Ex: {default_viewer.create_viewer: [Parameter, Parameter, ...]}
    viewer_to_parameters = collections.defaultdict(list)
    for param in parameters:
        # For the `sets` parameter in each of the parameter objects,
        # group them based on the viewer they are mapped to.
        # Parameters with the same viewer are grouped together.
        for set_name in param.sets:
            viewer = SET_TO_VIEWER[set_name]
            viewer_to_parameters[viewer].append(param)
    
    print(viewer_to_parameters)

    # A list of (title, url) tuples that each viewer generates.
    # This is used to create the main index.
    title_and_url_list = []
    # Now call the viewers with the parameters as the arguments.
    for viewer, parameters in viewer_to_parameters.items():
        titles_and_url_for_viewer = viewer(parameters)
        title_and_url_list.extend(titles_and_url_for_viewer)
    
    index_url = create_index(title_and_url_list)

    return index_url
