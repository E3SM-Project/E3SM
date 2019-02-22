import os
import collections
import cdms2
import acme_diags
from acme_diags.driver import utils
from acme_diags.metrics import mean
import cdutil

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np


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
        regions_to_data = collections.OrderedDict()
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
            # Make sure data have correct montly Bounds
            cdutil.setTimeBoundsMonthly(test)
            print('test shape',test.shape)

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

                cdutil.setTimeBoundsMonthly(ref)
 
                
                # TODO: Will this work if ref and test are timeseries data,
                # but land_frac and ocean_frac are climo'ed.
                test_domain, ref_domain = utils.general.select_region(
                    region, test, ref, land_frac, ocean_frac, parameter)

                # average over selected region, average over months to get yearly mean
                test_domain = cdutil.averager(test_domain,axis = 'xy')
                test_domain_year = cdutil.YEAR(test_domain)
                ref_domain = cdutil.averager(ref_domain,axis = 'xy')
                ref_domain_year = cdutil.YEAR(ref_domain)

                refs.append(ref_domain_year)

            #metrics_dict = create_metrics(ref_domain)
            metrics_dict = ref_domain_year.mean()
            print(test_domain_year.getTime().asComponentTime())
            print(test.getTime().asComponentTime())

            result = RefsTestMetrics(test=test_domain_year, refs=refs, metrics=metrics_dict)
            regions_to_data[region] = result
            
        #print(regions_to_data.values())
        # plot(parameter.current_set, data, parameter)
        
        plot(regions_to_data, parameter)
        # TODO: How will this work when there are a bunch of plots for each image?
        # Yes, these files should be saved.
        # utils.general.save_ncfiles(parameter.current_set,
        #                     mv1_domain, mv2_domain, diff, parameter)
    return parameter

def plot(regions_to_data, parameter):
    # Data is a list based on the len of the regions parameter.
    # Each element is a tuple with ref data,
    # test data, and metrics for that region.
    msg = 'We are plotting {} plots in this image,'.format(len(regions_to_data))
    msg += ' because regions = {}'.format(parameter.regions)
    print(msg)

    # plot time series
    plotTitle = {'fontsize': 8.5}
    plotSideTitle = {'fontsize': 6.5}
    
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # The dimensions [left, bottom, width, height] of the new axes. All quantities are in fractions of figure width and height.
    
    panel = [(0.1500, 0.5500, 0.7500, 0.3000),
             (0.1500, 0.1300, 0.7500, 0.3000),
             ]
    
    panel = [(0.1, 0.68, 0.25, 0.25),
             (0.4, 0.68, 0.25, 0.25),
             (0.7, 0.68, 0.25, 0.25),
             (0.1, 0.38, 0.25, 0.25),
             (0.4, 0.38, 0.25, 0.25),
             (0.7, 0.38, 0.25, 0.25),
             (0.1, 0.08, 0.25, 0.25),
             (0.4, 0.08, 0.25, 0.25),
             (0.7, 0.08, 0.25, 0.25),
             ]
    
    
    # Create figure
    figsize = [11.0,8.5]
    dpi = 150
    fig = plt.figure(figsize=figsize, dpi=dpi)
    #fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

    for i_region, data_set_for_region in enumerate(regions_to_data.values()):
        print('*****',i_region)
        refs = data_set_for_region.refs
        test = data_set_for_region.test
        ax1 = fig.add_axes(panel[i_region])
        print('refs',refs)
        print('test',test)
        ax1.plot(test.asma()[:-1], 'k', linewidth=2,label = parameter.test_name +'{0:.1f}'.format(np.mean(test.asma()[:-1])))
        for ref in refs:
            ax1.plot(ref.asma(), 'b', linewidth=2,label = 'label' +'{0:.1f}'.format(np.mean(ref.asma())))
        if i_region % 3 == 0 :
            ax1.set_ylabel('Surface air temperature (K)')
        ax1.legend(loc=1, prop={'size': 6})
        fig.text(panel[i_region][0]+0.12, panel[i_region][1]+panel[i_region][3]-0.015, parameter.regions[i_region],ha='center', color='black')
    plt.savefig('/global/project/projectdirs/acme/www/zhang40/figs/test_pr.png')


#        # You can have multiple reference data, so we
#        # make a list of all the data to plot.s
#        data_to_plot = refs + [test]
#        print('In this plot, we have {} data sets'.format(len(data_to_plot)))
#        #metrics = data_set_for_region.metrics
#    
#        for data in data_to_plot:
#            # Plot each of these data sets.
#            print('data',data)    
#            pass
