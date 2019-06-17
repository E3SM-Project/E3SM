import os
import collections
import cdms2
import cdutil
import acme_diags
from acme_diags.driver import utils
from acme_diags.metrics import mean
from acme_diags.plot.cartopy import area_mean_time_series_plot

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
    ref_names = parameter.ref_names

    # Both input data sets must be time-series files.
    # Raising an error will cause this specific set of
    # diagnostics with these parameters to be skipped.
    #if test_data.is_climo() or ref_data.is_climo():
    #    msg = 'Cannot run the plotset regional_mean_time_series '
    #    msg += 'because both the test and ref data need to be time-series files.'
    #    raise RuntimeError(msg)

    for var in variables:
        # The data that'll be sent to the plotting function.
        # There are eight tuples, each will be plotted like so:
        # [ 0 ]   [ 1 ]   [ 2 ]
        # [   ]   [   ]   [   ]
        #
        # [ 3 ]   [ 4 ]   [ 5 ]
        # [   ]   [   ]   [   ]
        #
        # [ 6 ]   [ 7 ]
        # [   ]   [   ]
        regions_to_data = collections.OrderedDict()
        print('Variable: {}'.format(var))

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
            
            # Make sure data have correct montly Bounds
            cdutil.setTimeBoundsMonthly(test)
            print('test shape',test.shape, test.units)

            parameter.viewer_descr[var] = getattr(test, 'long_name', var)
            # Get the name of the data, appended with the years averaged.
            parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)

            refs = []

            for ref_name in ref_names:    
                setattr(parameter, 'ref_name', ref_name)
                ref_data = utils.dataset.Dataset(parameter, ref=True)
            
                parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data)

                ref = ref_data.get_timeseries_variable(var)

                cdutil.setTimeBoundsMonthly(ref)
 
                
                # TODO: Will this work if ref and test are timeseries data,
                # but land_frac and ocean_frac are climo'ed.
                test_domain, ref_domain = utils.general.select_region(
                    region, test, ref, land_frac, ocean_frac, parameter)

                # Average over selected region, and average
                # over months to get the yearly mean.
                test_domain = cdutil.averager(test_domain,axis = 'xy')
                test_domain_year = cdutil.YEAR(test_domain)
                ref_domain = cdutil.averager(ref_domain,axis = 'xy')
                ref_domain_year = cdutil.YEAR(ref_domain)
                ref_domain_year.ref_name = ref_name

                refs.append(ref_domain_year)

            # metrics_dict = create_metrics(ref_domain)
            metrics_dict = ref_domain_year.mean()
            # print(test_domain_year.getTime().asComponentTime())
            # print(test.getTime().asComponentTime())

            result = RefsTestMetrics(test=test_domain_year, refs=refs, metrics=metrics_dict)
            regions_to_data[region] = result
 
        area_mean_time_series_plot.plot(var, regions_to_data, parameter)
        # TODO: How will this work when there are a bunch of plots for each image?
        # Yes, these files should be saved.
        # utils.general.save_ncfiles(parameter.current_set,
        #                     mv1_domain, mv2_domain, diff, parameter)
    return parameter

