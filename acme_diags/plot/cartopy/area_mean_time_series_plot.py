import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from acme_diags.driver.utils.general import get_output_dir


def plot(var, regions_to_data, parameter):
    # Data is a list based on the len of the regions parameter.
    # Each element is a tuple with ref data,
    # test data, and metrics for that region.

    # plot time series
    plotTitle = {'fontsize': 8.5}
    plotSideTitle = {'fontsize': 6.5}

    # Position and sizes of subplot axes in page coordinates (0 to 1)
    # The dimensions [left, bottom, width, height] of the new axes. All quantities are in fractions of figure width and height.
    line_color = ['r', 'b', 'g', 'm', 'c', 'y']

    panel = [(0.1500, 0.5500, 0.7500, 0.3000),
             (0.1500, 0.1300, 0.7500, 0.3000),
             ]

    panel = [(0.1, 0.70, 0.25, 0.225),
             (0.4, 0.70, 0.25, 0.225),
             (0.7, 0.70, 0.25, 0.225),
             (0.1, 0.38, 0.25, 0.225),
             (0.4, 0.38, 0.25, 0.225),
             (0.7, 0.38, 0.25, 0.225),
             (0.1, 0.06, 0.25, 0.225),
             (0.4, 0.06, 0.25, 0.225),
             (0.7, 0.06, 0.25, 0.225),
             ]

    # Create the figure.
    figsize = [17.0, 10]
    dpi = 150
    fig = plt.figure(figsize=figsize, dpi=dpi)
    start_time = int(parameter.start_yr)
    end_time = int(parameter.end_yr)
    num_year = end_time - start_time+1

    s = parameter.test_name_yrs
    years = s[s.find("(")+1:s.find(")")]
    test_name = s.split('(')[0].replace(' ', '')
    if test_name is '':
        test_name = 'test data'

    for i_region, data_set_for_region in enumerate(regions_to_data.values()):
        refs = data_set_for_region.refs
        test = data_set_for_region.test
        ax1 = fig.add_axes(panel[i_region])
        ax1.plot(test.asma(), 'k', linewidth=2,label = test_name +'(mean: {0:.2f}, std: {1:.3f})'.format(np.mean(test.asma()), np.std(test.asma())))
        for i_ref, ref in enumerate(refs):
            ax1.plot(ref.asma(), line_color[i_ref], linewidth=2,label = ref.ref_name +'(mean: {0:.2f}, std: {1:.3f})'.format(np.mean(ref.asma()), np.std(ref.asma())))

        x = np.arange(num_year)
        #do Truncation Division to accommodating long time records
        if num_year > 19:
            stepsize = num_year //10
        else:
            stepsize = 1
        ax1.set_xticks(x[::stepsize])
        x_ticks_labels = np.arange(start_time, end_time + 1, stepsize)
        ax1.set_xticklabels(x_ticks_labels, rotation='45', fontsize=8)

        #if i_region % 3 == 0 :
        #    ax1.set_ylabel(var + ' (' + test.units + ')')
        ax1.set_ylabel(var + ' (' + test.units + ')')
        ax1.set_xlabel('Year')
        ax1.legend(loc='best', prop={'size': 7})
        #fig.text(panel[i_region][0]+0.12, panel[i_region][1]+panel[i_region][3]-0.015, parameter.regions[i_region],ha='center', color='black')
    # Figure title.
        ax1.set_title(parameter.regions[i_region],fontsize = 10)
    fig.suptitle('Annual mean ' + var + ' time series over regions ' + parameter.test_name_yrs, x=0.3, y=0.98, fontsize=15)

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

