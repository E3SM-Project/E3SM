import os
import collections
import cdms2
import acme_diags
from acme_diags.driver import utils
from acme_diags.metrics import mean

RefsTestMetrics = collections.namedtuple('RefsTestMetrics', ['refs', 'test', 'metrics'])

def create_metrics(ref_domain):
    """
    For this plotset, calculate the mean of the
    reference data and return a dict of that.
    """
    return {'mean': mean(ref_domain)}


def run_diag(parameter):
    variables = parameter.variables
    regions = parameter.regions
    ref_names_regional_mean = getattr(parameter, 'ref_names_regional_mean', [])

    # Both input data sets must be time-series files.
    # Raising an error will cause this specific set of
    # diagnostics with these parameters to be skipped.
    #if test_data.is_climo() or ref_data.is_climo():
    #    msg = 'Cannot run the plotset regional_mean_time_series '
    #    msg += 'because both the test and ref data need to be time-series files.'
    #    raise RuntimeError(msg)

    for var in variables:
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
        print('Variable: {}'.format(var))
        parameter.var_id = var

        # Get land/ocean fraction for masking.
        # For now, we're only using the climo data that we saved below.
        # So no time-series LANDFRAC or OCNFRAC from the user is used.
        mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
        with cdms2.open(mask_path) as f:
            land_frac = f('LANDFRAC')
            ocean_frac = f('OCNFRAC')

        for region in regions:
            # The regions that are supported are in acme_diags/derivations/default_regions.py
            # You can add your own if it's not in there.
            print("Selected region: {}".format(region))
            test_data = utils.dataset.Dataset(parameter, test=True)
            test = test_data.get_timeseries_variable(var)

            # TODO: Remove this if using the A-Prime viewer.
            parameter.viewer_descr[var] = test.long_name if hasattr(
                test, 'long_name') else 'No long_name attr in test data.'
            # Get the name of the data, appended with the years averaged.
            parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)

            refs = []

            for ref_name in ref_names_regional_mean:    
                setattr(parameter, 'ref_name', ref_name)
                ref_data = utils.dataset.Dataset(parameter, ref=True)
            
                parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data)

                ref = ref_data.get_timeseries_variable(var)
                # TODO: Will this work if ref and test are timeseries data,
                # but land_frac and ocean_frac are climo'ed.
                test_domain, ref_domain = utils.general.select_region(
                    region, test, ref, land_frac, ocean_frac, parameter)

                refs.append(ref_domain)

            metrics_dict = create_metrics(ref_domain)

            result = RefsTestMetrics(test=test_domain, refs=refs, metrics=metrics_dict)
            data.append(result)
            
        # plot(parameter.current_set, data, parameter)
        plot(data, parameter)
        # TODO: How will this work when there are a bunch of plots for each image?
        # Yes, these files should be saved.
        # utils.general.save_ncfiles(parameter.current_set,
        #                     mv1_domain, mv2_domain, diff, parameter)

def plot(data, parameter):
    # Data is a list based on the len of the regions parameter.
    # Each element is a tuple with ref data,
    # test data, and metrics for that region.
    for single_data_set in data:
        refs = single_data_set.refs
        test = single_data_set.test
        # You can have multiple reference data, so we
        # make a list of all the data to plot.s
        data_to_plot = refs + [test]
        metrics = single_data_set.metrics
    
        for data in data_to_plot:
            # Plot each of these data sets.
            pass
