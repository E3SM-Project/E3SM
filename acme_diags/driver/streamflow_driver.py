from __future__ import print_function

import csv
import numpy
import scipy.io
import time
import cdutil
from acme_diags.driver import utils
from acme_diags.plot.cartopy.streamflow_plot import plot_seasonality


def get_drainage_area_error(radius, resolution, lon_ref, lat_ref, area_upstream, area_ref):
    k_bound = len(range(-radius, radius + 1))
    k_bound *= k_bound
    area_test = numpy.zeros((k_bound, 1))
    error_test = numpy.zeros((k_bound, 1))
    lat_lon_test = numpy.zeros((k_bound, 2))
    k = 0
    for i in range(-radius, radius + 1):
        for j in range(-radius, radius + 1):
            x = ((lon_ref + j * resolution) - (-180 + resolution / 2)) / resolution
            y = ((lat_ref + i * resolution) - (-90 + resolution / 2)) / resolution
            area_test[k] = area_upstream[int(x), int(y)] / 1000000
            error_test[k] = abs(area_test[k] - area_ref) / area_ref
            lat_lon_test[k, 0] = lat_ref + i * resolution
            lat_lon_test[k, 1] = lon_ref + j * resolution
            k += 1
    # The id of the center grid in the searching area
    center_id = (max(error_test.shape) - 1) / 2

    lat_lon_ref = [lat_ref, lon_ref]
    drainage_area_error = error_test[int(center_id)]
    return drainage_area_error, lat_lon_ref


def get_seasonality(monthly):
    # See https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2018MS001603 Equations 1 and 2
    num_years = monthly.shape[1]
    p_k = numpy.zeros((12, 1))
    # The total streamflow for each year (sum of Q_ij in the denominator of Equation 1, for all j)
    # 1 x num_years
    total_streamflow = numpy.sum(monthly, axis=0)
    for month in range(12):
        # The streamflow for this month in each year (Q_ij in the numerator of Equation 1, for all j)
        # 1 x num_years
        streamflow_month_all_years = monthly[month, :]
        # Proportion that this month contributes to streamflow that year.
        # 1 x num_years
        streamflow_proportion = streamflow_month_all_years / total_streamflow
        # The sum is the sum over j in Equation 1.
        # Dividing the sum of proportions by num_years gives the *average* proportion of annual streamflow during
        # this month.
        # Multiplying by 12 makes it so that Pk_i (`p_k[month]`) will be 1 if all months have equal streamflow and
        # 12 if all streamflow occurs in one month.
        # These steps produce the 12/n factor in Equation 1.
        p_k[month] = numpy.nansum(streamflow_proportion) * 12 / num_years
    # From Equation 2
    seasonality_index = numpy.max(p_k)
    # `p_k == numpy.max(p_k)` produces a Boolean matrix, True if the value (i.e., streamflow) is the max value.
    # `np.where(p_k == numpy.max(p_k))` produces the indices (i.e., months) where the max value is reached.
    peak_month = numpy.where(p_k == numpy.max(p_k))[0]
    # If more than one month has peak streamflow, simply define the peak month as the first one of the peak months.
    peak_month = peak_month[0]
    return seasonality_index, peak_month


def run_diag_seasonality(parameter):
    # Set path to the gauge metadata
    with open(parameter.gauges_path) as gauges_path:
        gauges_list = list(csv.reader(gauges_path))
    # Remove headers
    gauges_list.pop(0)
    gauges = numpy.array(gauges_list)
    if parameter.print_statements:
        print('gauges.shape={}'.format(gauges.shape))

    variables = parameter.variables
    for var in variables:

        if not parameter.ref_mat_file:
            ref_data = utils.dataset.Dataset(parameter, ref=True)
            ref_data_ts = ref_data.get_timeseries_variable(var)
            var_array = ref_data_ts(cdutil.region.domain(latitude=(-90., 90, 'ccb')))
            if parameter.print_statements:
                print('ref var original dimensions={}'.format(var_array.shape))
            var_transposed = numpy.transpose(var_array, (2, 1, 0))
            if parameter.print_statements:
                print('ref var transposed dimensions={}'.format(var_transposed.shape))
            ref_array = var_transposed
        else:
            # Load the observed streamflow dataset (GSIM)
            # the data has been reorganized to a 1380 * 30961 matrix. 1380 is the month
            # number from 1901.1 to 2015.12. 30961 include two columns for year and month plus
            # streamflow at 30959 gauge locations reported by GSIM
            ref_mat = scipy.io.loadmat(parameter.ref_mat_file)
            ref_array = ref_mat['GSIM']
        if parameter.print_statements:
            # GSIM: 1380 x 30961
            # wrmflow: 720 x 360 x 360
            print('ref_array.shape={}'.format(ref_array.shape))

        # Load E3SM simulated streamflow dataset
        if not parameter.test_mat_file:
            # `Dataset` will take the time slice from test_start_yr to test_end_yr
            test_data = utils.dataset.Dataset(parameter, test=True)
            test_data_ts = test_data.get_timeseries_variable(var)
            var_array = test_data_ts(cdutil.region.domain(latitude=(-90., 90, 'ccb')))
            if parameter.print_statements:
                print('test var original dimensions={}'.format(var_array.shape))
            var_transposed = numpy.transpose(var_array, (2, 1, 0))
            if parameter.print_statements:
                print('test var transposed dimensions={}'.format(var_transposed.shape))
            test_array = var_transposed
            areatotal2 = test_data.get_static_variable('areatotal2', var)
            area_upstream = numpy.transpose(areatotal2, (1, 0))
            if parameter.print_statements:
                print('area_upstream dimensions={}'.format(area_upstream.shape))
        else:
            data_mat = scipy.io.loadmat(parameter.test_mat_file)
            if 'E3SMflow' in data_mat.keys():
                # `edison` file uses this block
                e3sm_flow = data_mat['E3SMflow']
                if parameter.print_statements:
                    print('e3sm_flow = data_mat["E3SMflow"]')
            else:
                # `test` file uses this block
                e3sm_flow = data_mat
                if parameter.print_statements:
                    print('e3sm_flow = data_mat')
            try:
                if e3sm_flow['uparea'].shape == (1, 1):
                    # `edison` file uses this block
                    area_upstream = e3sm_flow['uparea'][0][0]
                    if parameter.print_statements:
                        print('e3sm_flow["uparea"] was indexed into for later use')
                else:
                    area_upstream = e3sm_flow['uparea']
                    if parameter.print_statements:
                        print('e3sm_flow["uparea"] will be used')
            except KeyError:
                # `test` file uses this block
                area_upstream = None
                if parameter.print_statements:
                    print('WARNING: uparea not found and will thus not be used')
            if e3sm_flow['wrmflow'].shape == (1, 1):
                # `edison` file uses this block
                test_array = e3sm_flow['wrmflow'][0][0]
                if parameter.print_statements:
                    print('e3sm_flow["wrmflow"] was indexed into for later use')
            else:
                # `test` file uses this block
                test_array = e3sm_flow['wrmflow']
                if parameter.print_statements:
                    print('e3sm_flow["wrmflow"] will be used')
        if parameter.print_statements:
            print('test_array.shape={}'.format(test_array.shape))

        # Resolution of MOSART output
        resolution = 0.5
        # Search radius (number of grids around the center point)
        radius = 1
        bins = numpy.floor(gauges[:, 7:9].astype(numpy.float64) / resolution)
        # Move the ref lat lon to grid center
        lat_lon = (bins + 0.5) * resolution
        if parameter.print_statements:
            print('lat_lon.shape={}'.format(lat_lon.shape))

        # Define the export matrix
        # Annual mean of test, annual mean of ref, error for area, lat, lon
        export = numpy.zeros((lat_lon.shape[0], 9))
        if parameter.print_statements:
            print('export.shape={}'.format(export.shape))
        t0 = time.time()
        for i in range(lat_lon.shape[0]):
            if parameter.print_statements and (i % 1000 == 0):
                print('On gauge #{}'.format(i))
            if parameter.max_num_gauges and i > parameter.max_num_gauges:
                break
            lat_ref = lat_lon[i, 1]
            lon_ref = lat_lon[i, 0]
            # Estimated drainage area (km^2) from ref
            area_ref = gauges[i, 13].astype(numpy.float64)

            if area_upstream is not None:
                drainage_area_error, lat_lon_ref = get_drainage_area_error(
                    radius, resolution, lon_ref, lat_ref, area_upstream, area_ref)
            else:
                # Use the center location
                lat_lon_ref = [lat_ref, lon_ref]
            if parameter.ref_mat_file:
                origin_id = gauges[i, 1].astype(numpy.int64)
                # Column 0 -- year
                # Column 1 -- month
                # Column origin_id + 1 -- the ref streamflow from gauge with the corresponding origin_id
                extracted = ref_array[:, [0, 1, origin_id + 1]]
                monthly_mean = numpy.zeros((12,1)) + numpy.nan
                for month in range(12):
                    # Add 1 to month to account for the months being 1-indexed
                    if sum(extracted[:, 1] == month + 1) > 0:
                        # `extracted[:,1]`: for all x, examine `extracted[x,1]`
                        # `extracted[:,1] == m`: Boolean array where 0 means the item in position [x,1] is NOT m,
                        # and 1 means it is m
                        # Example:
                        # a = [[1,2,3,4],
                        #      [5,6,7,8],
                        #      [1,2,3,4]]
                        # a[:,1]: [[2],
                        #          [6],
                        #          [2]]
                        # a[:,1] == 2: [[1], # False is 0, True is 1
                        #               [0],
                        #               [1]]
                        # a[a[:,1] == 2, 2]: [[3],
                        #                     [3]]
                        monthly_mean[month] = numpy.nanmean(extracted[extracted[:, 1] == month + 1, 2])
                # This is ref annual mean streamflow
                annual_mean_ref = numpy.mean(monthly_mean)
            if parameter.ref_mat_file and numpy.isnan(annual_mean_ref):
                # All elements of row i will be nan
                export[i,:] = numpy.nan
            else:
                if parameter.ref_mat_file:
                    # Reshape extracted[:,2] into a 12 x ? matrix; -1 means to calculate the size of the missing dimension.
                    mmat = numpy.reshape(extracted[:, 2], (12, -1))
                    mmat_id = numpy.sum(mmat, axis=0).transpose()
                    if numpy.sum(~numpy.isnan(mmat_id), axis=0) > 0:
                        # There's at least one year of record
                        monthly = mmat[:, ~numpy.isnan(mmat_id)]
                    else:
                        monthly = monthly_mean
                    seasonality_index_ref, peak_month_ref = get_seasonality(monthly)
                else:
                    ref_lon = int(1 + (lat_lon_ref[1] - (-180 + resolution / 2)) / resolution)
                    ref_lat = int(1 + (lat_lon_ref[0] - (-90 + resolution / 2)) / resolution)
                    ref = numpy.squeeze(ref_array[ref_lon, ref_lat, :])
                    mmat = numpy.reshape(ref, (12, -1))
                    monthly_mean_ref = numpy.nanmean(mmat, axis=1)
                    annual_mean_ref = numpy.mean(monthly_mean_ref)
                    if numpy.isnan(annual_mean_ref) == 1:
                        # The identified grid is in the ocean
                        monthly = numpy.ones((12, 1))
                    else:
                        monthly = mmat
                    seasonality_index_ref, peak_month_ref = get_seasonality(monthly)

                test_lon = int(1 + (lat_lon_ref[1] - (-180 + resolution / 2)) / resolution)
                test_lat = int(1 + (lat_lon_ref[0] - (-90 + resolution / 2)) / resolution)
                test = numpy.squeeze(test_array[test_lon, test_lat, :])
                mmat = numpy.reshape(test, (12, -1))
                monthly_mean_test = numpy.nanmean(mmat, axis=1)
                annual_mean_test = numpy.mean(monthly_mean_test)
                if numpy.isnan(annual_mean_test) == 1:
                    # The identified grid is in the ocean
                    monthly = numpy.ones((12, 1))
                else:
                    monthly = mmat
                seasonality_index_test, peak_month_test = get_seasonality(monthly)

                export[i, 0] = annual_mean_ref
                export[i, 1] = annual_mean_test
                if area_upstream is not None:
                    export[i, 2] = drainage_area_error * 100  # From fraction to percentage of the drainage area bias
                export[i, 3] = seasonality_index_ref  # Seasonality index of ref
                export[i, 4] = peak_month_ref  # Max flow month of ref
                export[i, 5] = seasonality_index_test  # Seasonality index of test
                export[i, 6] = peak_month_test  # Max flow month of test
                export[i, 7:9] = lat_lon_ref  # latlon of ref
        t1 = time.time()
        if parameter.print_statements:
            print('Loop time={}s={}m'.format(t1 - t0, (t1 - t0) / 60))
        # Remove the gauges with nan flow
        # `export[:,0]` => get first column of export
        # `numpy.isnan(export[:,0])` => Boolean column, True if value in export[x,0] is nan
        # `export[numpy.isnan(export[:,0]),:]` => rows of `export` where the Boolean column was True
        # Gauges will thus only be plotted if they have a non-nan value for both test and ref.
        if parameter.print_statements:
            print('export.shape before removing ref nan means={}'.format(export.shape))
        export = export[~numpy.isnan(export[:, 0]), :]
        if parameter.print_statements:
            print('export.shape before removing test nan means={}'.format(export.shape))
        export = export[~numpy.isnan(export[:, 1]), :]
        if parameter.print_statements:
            print('export.shape after both nan removals={}'.format(export.shape))

        if area_upstream is not None:
            # Set the max area error (percent) for all plots
            max_area_error = 20
            # `export[:,2]` gives the third column of `export`
            # `export[:,2]<=max_area_error` gives a Boolean array,
            # `True` if the value in the third column of `export` is `<= max_area_error`
            # `export[export[:,2]<=max_area_error,:]` is `export` with only the rows where the above is `True`.
            export = export[export[:, 2] <= max_area_error, :]
            if parameter.print_statements:
                print('export.shape after max_area_error cut={}'.format(export.shape))

        if parameter.print_statements:
            print('Variable: {}'.format(var))
        parameter.var_id = '{}-seasonality'.format(var)

        # This will be the title of the plot.
        parameter.main_title = 'Seasonality'
        parameter.viewer_descr[var] = parameter.main_title
        parameter.output_file = 'seasonality'

        # Plot
        # Plot original ref and test, not regridded versions.
        plot_seasonality(export, parameter)
    return parameter


def run_diag(parameter):
    # TODO: Add bias plot and scatterplot in future pull requests.
    if parameter.plot_type == 'seasonality':
        return run_diag_seasonality(parameter)
    else:
        raise Exception('Invalid plot_type={}'.format(parameter.plot_type))
