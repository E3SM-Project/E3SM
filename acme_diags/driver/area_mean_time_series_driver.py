import os
import collections
import cdms2
import cdutil
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import acme_diags
from acme_diags.driver import utils
from acme_diags.metrics import mean
from acme_diags.driver.utils.general import get_output_dir


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
    area_mean_ref_names = getattr(parameter, 'area_mean_time_series_ref_names', [])

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
            print('test shape',test.shape, test.units)

            # TODO: Remove this if using the A-Prime viewer.
            parameter.viewer_descr[var] = test.long_name if hasattr(
                test, 'long_name') else 'No long_name attr in test data.'
            # Get the name of the data, appended with the years averaged.
            parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)

            refs = []

            for ref_name in area_mean_ref_names:    
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
                ref_domain_year.ref_name = ref_name

                refs.append(ref_domain_year)

            #metrics_dict = create_metrics(ref_domain)
            metrics_dict = ref_domain_year.mean()
            print(test_domain_year.getTime().asComponentTime())
            print(test.getTime().asComponentTime())

            result = RefsTestMetrics(test=test_domain_year, refs=refs, metrics=metrics_dict)
            regions_to_data[region] = result
        #print('test_domain_year',test_domain_year)
            
            
        #print(regions_to_data.values())
        # plot(parameter.current_set, data, parameter)
        
        plot(var, regions_to_data, parameter)
        # TODO: How will this work when there are a bunch of plots for each image?
        # Yes, these files should be saved.
        # utils.general.save_ncfiles(parameter.current_set,
        #                     mv1_domain, mv2_domain, diff, parameter)
    return parameter

def plot(var, regions_to_data, parameter):
    # Data is a list based on the len of the regions parameter.
    # Each element is a tuple with ref data,
    # test data, and metrics for that region.

    # plot time series
    plotTitle = {'fontsize': 8.5}
    plotSideTitle = {'fontsize': 6.5}
    
    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # The dimensions [left, bottom, width, height] of the new axes. All quantities are in fractions of figure width and height.
    line_color = ['r', 'b', 'g', 'm']
    
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
    
    # Create the figure.
    figsize = [11.0, 8.5]
    dpi = 150
    fig = plt.figure(figsize=figsize, dpi=dpi)
    num_year = int(parameter.test_end_yr) - int(parameter.test_start_yr) +1

    for i_region, data_set_for_region in enumerate(regions_to_data.values()):
        refs = data_set_for_region.refs
        test = data_set_for_region.test
        ax1 = fig.add_axes(panel[i_region])
        ax1.plot(test.asma()[:-1], 'k', linewidth=2,label = 'model' +' ({0:.1f})'.format(np.mean(test.asma()[:-1])))
        for i_ref, ref in enumerate(refs):
            ax1.plot(ref.asma(), line_color[i_ref], linewidth=2,label = ref.ref_name +' ({0:.1f})'.format(np.mean(ref.asma())))

        x = np.arange(num_year)
        ax1.set_xticks(x)
        x_ticks_labels = [str(x) for x in range(int(parameter.test_start_yr),int(parameter.test_end_yr)+1)]
        ax1.set_xticklabels(x_ticks_labels, rotation='45', fontsize=6)

        if i_region % 3 == 0 :
            ax1.set_ylabel('variable name (units)')
            #ax1.set_ylabel(test.long_name + ' (' + test.units + ')')
        ax1.legend(loc=1, prop={'size': 6})
        fig.text(panel[i_region][0]+0.12, panel[i_region][1]+panel[i_region][3]-0.015, parameter.regions[i_region],ha='center', color='black')
    # Figure title.
    fig.suptitle('Annual mean ' + var + ' over regions ' + parameter.test_name_yrs, x=0.5, y=0.97, fontsize=15)

    # Save the figure.
    output_file_name = var
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        fnm = os.path.join(get_output_dir(parameter.current_set,
            parameter), output_file_name + '.' + f)
        plt.savefig(fnm)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(get_output_dir(parameter.current_set, parameter,
            ignore_container=True), output_file_name + '.' + f)
        print('Plot saved in: ' + fnm)

    plt.close()
