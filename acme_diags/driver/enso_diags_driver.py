from __future__ import print_function

import json
import math
import numpy
import os
import scipy.stats
import cdms2
import cdutil
import acme_diags
from acme_diags.derivations import acme, default_regions
from acme_diags.driver import utils
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean, std
from acme_diags.plot.cartopy.enso_diags_plot import plot_map, plot_scatter


def calculate_nino_index(nino_region_str, parameter):
    """
    Use the built-in HadISST nino index time series from http://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/
    for observation datasets and when model output is not available.

    Relevant files are acme_diags/driver/default_diags/enso_NINO{3, 34, 4}.long.data
    """
    data_file = ''.join(['enso_', nino_region_str, '.long.data'])
    nino_index_path = os.path.join(acme_diags.INSTALL_PATH, 'enso_diags', data_file)
    sst_orig = numpy.loadtxt(nino_index_path, skiprows=1, max_rows=149)  # load data up to year 2018 from 1870
    start = int(parameter.start_yr)
    end = int(parameter.end_yr)
    sst_years = sst_orig[:,0].astype(int)
    try:
        start_ind = numpy.where(sst_years == start)[0][0]
        end_ind = numpy.where(sst_years == end)[0][0]
    except:
        msg = 'Requested years are outside of available sst obs records.'
        raise RuntimeError(msg)

    sst = sst_orig[start_ind:end_ind+1,1:]
    # Get anomaly from annual cycle climatology
    annual_cycle = numpy.mean(sst, axis=0)
    sst_anomaly = numpy.empty_like(sst)
    num_years = end - start + 1
    for iyr in range(num_years):
        sst_anomaly[iyr,:] = sst[iyr,:] - annual_cycle
    nino_index = numpy.reshape(sst_anomaly, (num_years*12))

    if parameter.print_statements:
        print('nino_index_obs', nino_index)
    return nino_index


def calculate_nino_index_model(data, nino_region_str, parameter):
    """
    Calculate nino index based on model output TS. If this is not available use the built-in HadISST nino index
    """

    try:
        sst = data.get_timeseries_variable('SST')
        nino_region = default_regions.regions_specs[nino_region_str]['domain']
        sst_nino = sst(nino_region)
        # Domain average
        sst_avg = cdutil.averager(sst_nino, axis='xy')
        # Get anomaly from annual cycle climatology
        sst_avg_anomaly = cdutil.ANNUALCYCLE.departures(sst_avg)
        nino_index = sst_avg_anomaly
    except:
        nino_index = calculate_nino_index(nino_region_str, parameter)
        print("Simulated surface temperature not found, using built-in HadISST nino index time series.")

    if parameter.print_statements:
        print('nino_index_model', nino_index)

    return nino_index


def perform_regression(data, parameter, var, region, land_frac, ocean_frac, nino_index):
    ts_var = data.get_timeseries_variable(var)
    domain = utils.general.select_region(region, ts_var, land_frac, ocean_frac, parameter)
    # Average over selected region, and average
    # over months to get the yearly mean.
    cdutil.setTimeBoundsMonthly(domain)
    # Get anomaly from annual cycle climatology
    if parameter.print_statements:
        print('domain.shape: {}'.format(domain.shape))
    anomaly = cdutil.ANNUALCYCLE.departures(domain)
    nlat = len(anomaly.getLatitude())
    nlon = len(anomaly.getLongitude())
    reg_coe = anomaly[0, :, :](squeeze=1)
    confidence_levels = cdutil.ANNUALCYCLE.departures(domain)[0, :, :](squeeze=1)
    # Neither of the following methods work, so we just set values in confidence_levels
    # to be explicitly 0 or 1.
    #confidence_levels = anomaly[0, :, :](squeeze=1).fill(0)
    #confidence_levels = numpy.zeros_like(reg_coe)
    for ilat in range(nlat):
        if parameter.print_statements:
            print('ilat: {}'.format(ilat))
        for ilon in range(nlon):
            dependent_var = anomaly[:, ilat, ilon]
            independent_var = nino_index
            # Uncomment the following line to use CDAT/genutil instead
            # (You'll also need to set pvalue)
            #slope, intercept = genutil.statistics.linearregression(dependent_var, x=independent_var)
            slope, _, _, pvalue, _ = scipy.stats.linregress(independent_var, dependent_var)
            reg_coe[ilat, ilon] = slope
            # Set confidence level to 1 if significant and 0 if not
            if pvalue < 0.05:
                # p-value < 5%
                # This implies significance at 95% confidence level
                confidence_levels[ilat, ilon] = 1
            else:
                confidence_levels[ilat, ilon] = 0
    if parameter.print_statements:
        print("confidence in fn:", confidence_levels.shape)
    sst_units = 'degC'
    reg_coe.units = '{}/{}'.format(ts_var.units, sst_units)
    if parameter.print_statements:
        print('reg_coe.shape: {}'.format(reg_coe.shape))
    return domain, reg_coe, confidence_levels


def create_single_metrics_dict(values):
    d = {
        'min': float(min_cdms(values)),
        'max': float(max_cdms(values)),
        'mean': float(mean(values)),
        'std': float(std(values))
    }
    return d


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the relevant metrics in a dictionary"""
    metrics_dict = {}
    metrics_dict['ref'] = create_single_metrics_dict(ref)
    metrics_dict['ref_regrid'] = create_single_metrics_dict(ref_regrid)
    metrics_dict['test'] = create_single_metrics_dict(test)
    metrics_dict['test_regrid'] = create_single_metrics_dict(test_regrid)
    metrics_dict['diff'] = create_single_metrics_dict(diff)
    d = metrics_dict['diff']
    d['rmse'] = float(rmse(test_regrid, ref_regrid))
    d['corr'] = float(corr(test_regrid, ref_regrid))
    return metrics_dict


def run_diag_map(parameter):
    variables = parameter.variables
    seasons = parameter.seasons
    regions = parameter.regions
    nino_region_str = parameter.nino_region
    run_type = parameter.run_type

    if parameter.print_statements:
        print('run_type: {}'.format(run_type))
    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)
    if run_type == 'model_vs_model':
        test_nino_index = calculate_nino_index_model(test_data, nino_region_str, parameter)
        ref_nino_index = calculate_nino_index_model(ref_data, nino_region_str, parameter)
    elif run_type == 'model_vs_obs':
        test_nino_index = calculate_nino_index_model(test_data, nino_region_str, parameter)
        ref_nino_index = calculate_nino_index(nino_region_str, parameter)
    else:
        raise Exception('Invalid run_type={}'.format(run_type))
    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)

    for season in seasons:
        if parameter.print_statements:
            print('Season: {}'.format(season))
        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data, season)
        parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data, season)

        # Get land/ocean fraction for masking.
        try:
            land_frac = test_data.get_climo_variable('LANDFRAC', season)
            ocean_frac = test_data.get_climo_variable('OCNFRAC', season)
        except:
            mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
            with cdms2.open(mask_path) as f:
                land_frac = f('LANDFRAC')
                ocean_frac = f('OCNFRAC')

        for var in variables:
            if parameter.print_statements:
                print('Variable: {}'.format(var))
            parameter.var_id = '{}-regression-over-nino'.format(var)

            for region in regions:
                if parameter.print_statements:
                    print("Selected region: {}".format(region))

                # This will be the title of the plot.
                parameter.main_title = 'Regression coefficient, {} anomaly to {}'.format(var, nino_region_str)
                parameter.viewer_descr[var] = ', '.join(
                    [parameter.main_title, region])

                # Test
                test_domain, test_reg_coe, test_confidence_levels = perform_regression(test_data, parameter, var,
                                                                                       region, land_frac,
                                                                                       ocean_frac, test_nino_index)

                # Reference
                ref_domain, ref_reg_coe, ref_confidence_levels = perform_regression(ref_data, parameter, var,
                                                                                    region, land_frac, ocean_frac,
                                                                                    ref_nino_index)

                # Difference
                # Regrid towards the lower resolution of the two variables for calculating the difference.
                test_reg_coe_regrid, ref_reg_coe_regrid = utils.general.regrid_to_lower_res(test_reg_coe,
                                                                                            ref_reg_coe,
                                                                                            parameter.regrid_tool,
                                                                                            parameter.regrid_method)
                diff = test_reg_coe_regrid - ref_reg_coe_regrid

                # Metrics
                metrics_dict = create_metrics(ref_reg_coe, test_reg_coe, ref_reg_coe_regrid, test_reg_coe_regrid,
                                              diff)
                # If not defined, determine contour_levels
                if not parameter.contour_levels or not parameter.diff_levels:
                    # We want contour levels for the plot,
                    # which uses original (non-regridded) test and ref,
                    # so we use those min and max values.
                    min_contour_level = math.floor(min(
                        metrics_dict['ref']['min'],
                        metrics_dict['test']['min']
                    ))
                    max_contour_level = math.ceil(max(
                        metrics_dict['ref']['max'],
                        metrics_dict['test']['max']
                    ))
                    CENTER_ON_ZERO = True
                    if CENTER_ON_ZERO:
                        bound = max(abs(max_contour_level), abs(min_contour_level))
                        lower_bound = -bound
                        upper_bound = bound
                        contour_level_range = 2 * bound
                    else:
                        lower_bound = min_contour_level
                        upper_bound = max_contour_level
                        contour_level_range = max_contour_level - min_contour_level
                    step_size = (contour_level_range / 10)
                    if step_size > 1:
                        step_size = int(step_size)
                    elif step_size == 0:
                        step_size = 1 / 10
                    contour_levels = list(numpy.arange(lower_bound, upper_bound + 1, step_size))
                    if not parameter.contour_levels:
                        parameter.contour_levels = contour_levels
                    if not parameter.diff_levels:
                        parameter.diff_levels = contour_levels
                parameter.output_file = 'regression-coefficient-{}-over-{}'.format(var.lower(),
                                                                                   nino_region_str.lower())

                # Saving the metrics as a json.
                metrics_dict['unit'] = test_reg_coe_regrid.units
                metrics_output_file_name = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + '.json')
                with open(metrics_output_file_name, 'w') as outfile:
                    json.dump(metrics_dict, outfile)
                # Get the file name that the user has passed in and display that.
                # When running in a container, the paths are modified.
                metrics_output_file_name = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter, ignore_container=True),
                    parameter.output_file + '.json')
                print('Metrics saved in: {}'.format(metrics_output_file_name))

                # Plot
                parameter.var_region = region
                # Plot original ref and test, not regridded versions.
                plot_map(ref_reg_coe, test_reg_coe, diff,
                         metrics_dict, ref_confidence_levels, test_confidence_levels,
                         parameter)
                utils.general.save_ncfiles(parameter.current_set,
                                           test_reg_coe, ref_reg_coe, diff, parameter)
    return parameter


def run_diag_scatter(parameter):
    variables = parameter.variables
    run_type = parameter.run_type
    # We will always use the same regions, so we don't do the following:
    # x['region'] = parameter.nino_region
    # y['region'] = parameter.regions[0]
    x = {'var': 'TS', 'units': 'degC', 'region': 'NINO3'}
    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)
    if parameter.print_statements:
        print('run_type: {}'.format(run_type))
    if run_type == 'model_vs_model':
        x['test'] = calculate_nino_index_model(test_data, x['region'], parameter)
        x['ref'] = calculate_nino_index_model(ref_data, x['region'], parameter)
    elif run_type == 'model_vs_obs':
        x['test'] = calculate_nino_index_model(test_data, x['region'], parameter)
        x['ref'] = calculate_nino_index(x['region'], parameter)
    else:
        raise Exception('Invalid run_type={}'.format(run_type))
    parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)
    parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data)
    for y_var in variables:
        if y_var == 'TAUX':
            regions = ['NINO4']
        else:
            regions = ['NINO3']
        for region in regions:
            y = {'var': y_var, 'region': region}
            test_data_ts = test_data.get_timeseries_variable(y_var)
            ref_data_ts = ref_data.get_timeseries_variable(y_var)
            y_region = default_regions.regions_specs[region]['domain']
            test_data_ts_regional = test_data_ts(y_region)
            ref_data_ts_regional = ref_data_ts(y_region)
            # Domain average
            test_avg = cdutil.averager(test_data_ts_regional, axis='xy')
            ref_avg = cdutil.averager(ref_data_ts_regional, axis='xy')
            # Get anomaly from annual cycle climatology
            y['test'] = cdutil.ANNUALCYCLE.departures(test_avg)
            y['ref'] = cdutil.ANNUALCYCLE.departures(ref_avg)
            y['units'] = test_avg.units
            if y_var == 'TAUX':
                y['test'] *= 1000
                y['ref'] *= 1000
                y['units'] = '10^3 {}'.format(y['units'])
            parameter.var_id = '{}-feedback'.format(y['var'])
            title_tuple = (y['var'], y['region'], x['var'], x['region'])
            parameter.main_title = '{} anomaly ({}) vs. {} anomaly ({})'.format(*title_tuple)
            parameter.viewer_descr[y['var']] = parameter.main_title
            parameter.output_file = 'feedback-{}-{}-{}-{}'.format(*title_tuple)
            plot_scatter(x, y, parameter)
    return parameter


def run_diag(parameter):
    if parameter.plot_type == 'map':
        return run_diag_map(parameter)
    elif parameter.plot_type == 'scatter':
        return run_diag_scatter(parameter)
    else:
        raise Exception('Invalid plot_type={}'.format(parameter.plot_type))
