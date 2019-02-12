from collections import namedtuple
from acme_diags.driver import utils
from acme_diags.metrics import mean

RefTestMetrics = namedtuple('RefTestMetrics', ['ref', 'test', 'metrics'])

def create_metrics(ref_domain):
    """
    For this plotset, calculate the mean of the
    reference data and return a dict of that.
    """
    return {'mean': mean(ref_domain)}


def run_diag(parameter):
    variables = parameter.variables
    seasons = parameter.seasons
    # TODO: Idk if this is needed:
    ref_name = getattr(parameter, 'ref_name', '')
    regions = parameter.regions

    # Initialize the Dataset class here.
    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)

    # Both input data sets must be time-series files.
    # Raising an error will cause this specific set of
    # diagnostics with these parameters to be skipped.
    #if test_data.is_climo() or ref_data.is_climo():
    #    msg = 'Cannot run the plotset regional_mean_time_series '
    #    msg += 'because both the test and ref data need to be time-series files.'
    #    raise RuntimeError(msg)

    # The data that'll be sent to the plotting function.
    # There are six tuples, each will be plotted like so:
    # [ 0 ]   [ 1 ]
    # [   ]   [   ]
    #
    # [ 2 ]   [ 3 ]
    # [   ]   [   ]
    #
    # [ 4 ]   [ 5 ]
    # [   ]   [   ]
    data = []

    for season in seasons:
        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data, season)
        parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data, season)

        # Get land/ocean fraction for masking.
        # TODO: Is this needed?
        try:
            land_frac = test_data.get_climo_variable('LANDFRAC', season)
            ocean_frac = test_data.get_climo_variable('OCNFRAC', season)
        except:
            mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
            with cdms2.open(mask_path) as f:
                land_frac = f('LANDFRAC')
                ocean_frac = f('OCNFRAC')

        for var in variables:
            print('Variable: {}'.format(var))
            parameter.var_id = var

            test = test_data.get_timeseries_variable(var, season)
            ref = ref_data.get_timeseries_variable(var, season)

            # TODO: Remove this if using the A-Prime viewer.
            parameter.viewer_descr[var] = test.long_name if hasattr(
                test, 'long_name') else 'No long_name attr in test data.'

            for region in regions:
                # The regions that are supported are in acme_diags/derivations/default_regions.py
                # You can add your own if it's not in there.
                print("Selected region: {}".format(region))

                # TODO: Will this work if ref and test are timeseries data,
                # but land_frac and ocean_frac are climo'ed.
                test_domain, ref_domain = utils.general.select_region(
                    region, test, ref, land_frac, ocean_frac, parameter)

                metrics_dict = create_metrics(ref_domain)

                result = RefTestMetrics(test=test_domain, ref=ref_domain, metrics=metrics_dict)
                data.append(result)
            
            # plot(parameter.current_set, data, parameter)
            plot(data, parameter)
            # TODO: How will this work when there are a bunch of plots for each image?
            # utils.general.save_ncfiles(parameter.current_set,
            #                     mv1_domain, mv2_domain, diff, parameter)

def plot(data, parameter):
    # Data is a list of size 6, where each element is a
    # tuple with ref data, test data, and metrics.
    for single_data_set in data:
        ref = single_data_set.ref
        test = single_data_set.test
        metrics = single_data_set.metrics
