from __future__ import print_function

import json
import numpy as np
import os
import scipy.fftpack
import cdutil
from acme_diags.derivations import acme, default_regions
from acme_diags.driver import utils
from acme_diags.plot.cartopy.qbo_plot import plot


def process_u_for_time_height(data_region):
    # Average over longitude (i.e., each latitude's average in data_region)
    data_lon_average = cdutil.averager(data_region, axis='x')
    # Average over latitude (i.e., average for entire data_region)
    data_lon_lat_average = cdutil.averager(data_lon_average, axis='y')
    # Get data by vertical level
    level_data = data_lon_lat_average.getAxis(1)
    return data_lon_lat_average, level_data


def deseason(xraw):
    # Calculates the deseasonalized data
    months_per_year = 12
    # Create array to hold climatological values and deseasonalized data
    # Create months_per_year x 1 array of zeros
    xclim = np.zeros((months_per_year, 1))
    # Create array with same shape as xraw
    x_deseasoned = np.zeros(xraw.shape)
    # Iterate through all 12 months.
    for month in np.arange(months_per_year):
        # `xraw[month::12]` will return the data for this month every year (12 months)
        # (i.e., from month until the end of xraw, get every 12th month)
        # Get the mean of this month, using data from every year, ignoring NaNs
        xclim[month] = np.nanmean(xraw[month::months_per_year])
    num_years = int(np.floor(len(x_deseasoned) / months_per_year))
    # Iterate through all years in x_deseasoned (same number as in xraw)
    for year in np.arange(num_years):
        year_index = year * months_per_year
        # Iterate through all months of the year
        for month in np.arange(months_per_year):
            month_index = year_index + month
            # Subtract the month's mean over num_years from xraw's data for this month in this year
            # i.e., get the difference between this month's value and it's "usual" value
            x_deseasoned[month_index] = xraw[month_index] - xclim[month]
    return x_deseasoned


def get_20to40month_fft_amplitude(qboN, levelN):
    # Calculates the amplitude of wind variations in the 20 - 40 month period
    psd_sumN = np.zeros(levelN.shape)
    amplitudeN = np.zeros(levelN.shape)

    for ilev in np.arange(len(levelN)):
        # `qboN[:, ilev]` returns the entire 0th dimension for ilev in the 1st dimension of the array.
        y_input = deseason(np.squeeze(qboN[:, ilev]))
        y = scipy.fftpack.fft(y_input)
        n = len(y)
        frequency = np.arange(n / 2) / n
        period = 1 / frequency
        values = y[0:int(np.floor(n / 2))]
        fyy = values * np.conj(values)
        # Choose the range 20 - 40 months that captures most QBOs (in nature)
        psd_sumN[ilev] = 2 * np.nansum(fyy[(period <= 40) & (period >= 20)])
        amplitudeN[ilev] = np.sqrt(2 * psd_sumN[ilev]) * (frequency[1] - frequency[0])
    return psd_sumN, amplitudeN


def process_u_for_power_spectral_density(data_region):
    # Average over vertical levels and horizontal area
    level_bottom = 22
    level_top = 18
    # Average over lat and lon
    data_lat_lon_average = cdutil.averager(data_region, axis='xy')
    # Average over vertical
    try:
        average = data_lat_lon_average(level=(level_top, level_bottom))
    except:
        raise Exception('No levels found between {}hPa and {}hPa'.format(level_top, level_bottom))
    x0 = np.nanmean(np.array(average), axis=1)
    # x0 should now be 1D
    return x0


def ceil_log2(x):
    """
    Given a number, calculate the exponent for the next power of 2.

    Example:
        ceil_log2(16) = 4
        ceil_log2(17) = 5
    """
    return np.ceil(np.log2(x)).astype('int')


def get_psd_from_deseason(xraw, period_new):
    x_deseasoned = deseason(xraw)

    # Sampling frequency: assumes frequency of sampling = 1 month
    sampling_frequency = 1
    # Calculate the period as a function of frequency
    period0 = 1 / sampling_frequency
    L0 = len(xraw)
    t0 = np.arange(0, L0) * period0
    NFFT0 = 2 ** ceil_log2(L0)

    # Apply fft on x_deseasoned with n = NFFT
    x0 = scipy.fftpack.fft(x_deseasoned, n=NFFT0) / L0
    # Frequency (cycles/month). Frequency will be increasing.
    frequency0 = sampling_frequency * np.arange(0, (NFFT0 / 2 + 1)) / NFFT0
    # Period (months/cycle). Calculate as a function of frequency. Thus, period will be decreasing.
    period0 = 1 / frequency0

    # Calculate amplitude as a function of frequency
    amplitude0 = 2 * abs(x0[0:int(NFFT0 / 2 + 1)])
    # Calculate power spectral density as a function of frequency
    psd_x0 = amplitude0 ** 2 / L0
    period_end0 = period0 * L0
    # Total spectral power
    Pxf0 = period_end0 * np.sum(psd_x0)
    # In the next code block, we will perform an interpolation using the period
    # (interpolating values of amplitude0_flipped and psd_x0_flipped from period0_flipped to period_new).
    # For that interpolation, we want the period to be increasing.
    # Therefore, we will flip the following values:
    period0_flipped = period0[::-1]
    amplitude0_flipped = amplitude0[::-1]
    psd_x0_flipped = psd_x0[::-1]

    amplitude_new0 = np.interp(period_new, period0_flipped[:-1], amplitude0_flipped[:-1])
    psd_x_new0 = np.interp(period_new, period0_flipped[:-1], psd_x0_flipped[:-1])
    return psd_x_new0, amplitude_new0


def run_diag(parameter):
    variables = parameter.variables
    # The region will always be 5S5N
    region = '5S5N'
    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)
    for variable in variables:
        if parameter.print_statements:
            print('Variable={}'.format(variable))
        test_var = test_data.get_timeseries_variable(variable)
        ref_var = ref_data.get_timeseries_variable(variable)
        qbo_region = default_regions.regions_specs[region]['domain']

        test_region = test_var(qbo_region)
        ref_region = ref_var(qbo_region)

        test = {}
        ref = {}

        # Diagnostic 1: average over longitude & latitude to produce time-height array of u field:
        # Richter, J. H., Chen, C. C., Tang, Q., Xie, S., & Rasch, P. J. (2019). Improved Simulation of the QBO in E3SMv1. Journal of Advances in Modeling Earth Systems, 11(11), 3403-3418.
        # https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2019MS001763
        # U = "Monthly mean zonal mean zonal wind averaged between 5S and 5N as a function of pressure and time" (p. 3406)
        test['qbo'], test['level'] = process_u_for_time_height(test_region)
        ref['qbo'], ref['level'] = process_u_for_time_height(ref_region)

        # Diagnostic 2: calculate and plot the amplitude of wind variations with a 20-40 month period
        test['psd_sum'], test['amplitude'] = get_20to40month_fft_amplitude(np.squeeze(np.array(test['qbo'])), test['level'])
        ref['psd_sum'], ref['amplitude'] = get_20to40month_fft_amplitude(np.squeeze(np.array(ref['qbo'])), ref['level'])

        # Diagnostic 3: calculate the Power Spectral Density
        # Pre-process data to average over lat,lon,height
        x_test = process_u_for_power_spectral_density(test_region)
        x_ref = process_u_for_power_spectral_density(ref_region)
        # Calculate the PSD and interpolate to period_new. Specify periods to plot
        period_new = np.concatenate((np.arange(2, 33), np.arange(34, 100, 2)), axis=0)
        test['psd_x_new'], test['amplitude_new'] = get_psd_from_deseason(x_test, period_new)
        ref['psd_x_new'], ref['amplitude_new'] = get_psd_from_deseason(x_ref, period_new)

        parameter.var_id = variable
        parameter.main_title = 'QBO index, amplitude, and power spectral density for {}'.format(variable)
        parameter.viewer_descr[variable] = parameter.main_title

        test['name'] = parameter.test_name
        ref['name'] = parameter.ref_name

        test_nc = {}
        ref_nc = {}
        test_json = {}
        ref_json = {}
        for key in test.keys():
            if key in ['qbo', 'level']:
                test_nc[key] = test[key]
                ref_nc[key] = ref[key]
            else:
                if key == 'name':
                    test_json[key] = test[key]
                    ref_json[key] = ref[key]
                else:
                    test_json[key] = list(test[key])
                    ref_json[key] = list(ref[key])

        # TODO: Check the below works properly by using ncdump on Cori
        utils.general.save_transient_variables_to_netcdf(parameter.current_set, test_nc, 'test', parameter)
        utils.general.save_transient_variables_to_netcdf(parameter.current_set, ref_nc, 'ref', parameter)

        # Saving the other data as json.
        for dict_type in ['test', 'ref']:
            json_output_file_name = os.path.join(
                utils.general.get_output_dir(parameter.current_set, parameter),
                parameter.output_file + '_{}.json'.format(dict_type))
            with open(json_output_file_name, 'w') as outfile:
                if dict_type == 'test':
                    json_dict = test_json
                else:
                    json_dict = ref_json
                json.dump(json_dict,  outfile)
            # Get the file name that the user has passed in and display that.
            # When running in a container, the paths are modified.
            json_output_file_name = os.path.join(
                utils.general.get_output_dir(parameter.current_set,
                                             parameter, ignore_container=True),
                parameter.output_file + '_{}.json'.format(dict_type))
            print('Metrics saved in: {}'.format(json_output_file_name))

        plot(period_new, parameter, test, ref)
    return parameter
