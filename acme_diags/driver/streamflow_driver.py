from __future__ import print_function

import csv
import os

import cdms2
import cdutil
import numpy
import scipy.io

from acme_diags.driver import utils
from acme_diags.plot.cartopy.streamflow_plot import (
    plot_annual_map,
    plot_annual_scatter,
    plot_seasonality_map,
)


def get_drainage_area_error(
    radius, resolution, lon_ref, lat_ref, area_upstream, area_ref
):
    k_bound = len(range(-radius, radius + 1))
    k_bound *= k_bound
    area_test = numpy.zeros((k_bound, 1))
    error_test = numpy.zeros((k_bound, 1))
    lat_lon_test = numpy.zeros((k_bound, 2))
    k = 0
    for i in range(-radius, radius + 1):
        for j in range(-radius, radius + 1):
            x = int(
                1
                + ((lon_ref + j * resolution) - (-180 + resolution / 2))
                / resolution
            )
            y = int(
                1
                + ((lat_ref + i * resolution) - (-90 + resolution / 2))
                / resolution
            )
            area_test[k] = area_upstream[x - 1, y - 1] / 1000000
            error_test[k] = numpy.abs(area_test[k] - area_ref) / area_ref
            lat_lon_test[k, 0] = lat_ref + i * resolution
            lat_lon_test[k, 1] = lon_ref + j * resolution
            k += 1
    # The id of the center grid in the searching area
    center_id = (k_bound - 1) / 2

    lat_lon_ref = [lat_ref, lon_ref]
    drainage_area_error = error_test[int(center_id)]
    return drainage_area_error, lat_lon_ref


def get_seasonality(monthly):
    monthly = monthly.astype(numpy.float64)
    # See https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2018MS001603 Equations 1 and 2
    if monthly.shape[0] != 12:
        raise Exception(
            "monthly.shape={} does not include 12 months".format(monthly.shape)
        )
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
        # For all i, divide streamflow_month_all_years[i] by total_streamflow[i]
        streamflow_proportion = numpy.divide(
            streamflow_month_all_years, total_streamflow
        )
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
    # Month 0 is January, Month 1 is February, and so on.
    peak_month = peak_month[0]
    return seasonality_index, peak_month


# FIXME:C901 'run_diag' is too complex (57)
def run_diag(parameter):  # noqa
    # Assume `model` will always be a `nc` file.
    # Assume `obs` will always be a `mat` file.
    using_test_mat_file = False
    if parameter.run_type == "model_vs_model":
        using_ref_mat_file = False
        if parameter.gauges_path is None:
            raise Exception(
                "To use a non-GSIM reference, please specify streamflow_param.gauges_path. This might be {}/{}".format(
                    parameter.reference_data_path.rstrip("/"),
                    "GSIM/GSIM_catchment_characteristics_all_1km2.csv",
                )
            )
        else:
            gauges_path = parameter.gauges_path
    elif parameter.run_type == "model_vs_obs":
        using_ref_mat_file = True
        # The metadata file of GSIM that has observed gauge lat lon and drainage area
        # This file includes 25765 gauges, which is a subset of the entire
        # dataset (30959 gauges). The removed gauges are associated with very
        # small drainage area (<1km2), which is not meaningful to be included.
        gauges_path = (
            "{}/GSIM/GSIM_catchment_characteristics_all_1km2.csv".format(
                parameter.reference_data_path.rstrip("/")
            )
        )
    else:
        raise Exception(
            "parameter.run_type={} not supported".format(parameter.run_type)
        )

    # Set path to the gauge metadata
    with open(gauges_path) as gauges_file:
        gauges_list = list(csv.reader(gauges_file))
    # Remove headers
    gauges_list.pop(0)
    gauges = numpy.array(gauges_list)
    if parameter.print_statements:
        print("gauges.shape={}".format(gauges.shape))

    variables = parameter.variables
    for var in variables:

        if not using_ref_mat_file:
            ref_data = utils.dataset.Dataset(parameter, ref=True)
            parameter.ref_name_yrs = utils.general.get_name_and_yrs(
                parameter, ref_data
            )
            ref_data_ts = ref_data.get_timeseries_variable(var)
            var_array = ref_data_ts(
                cdutil.region.domain(latitude=(-90.0, 90, "ccb"))
            )
            if parameter.print_statements:
                print("ref var original dimensions={}".format(var_array.shape))
            var_transposed = numpy.transpose(var_array, (2, 1, 0))
            if parameter.print_statements:
                print(
                    "ref var transposed dimensions={}".format(
                        var_transposed.shape
                    )
                )
            ref_array = var_transposed.astype(numpy.float64)
        else:
            # Load the observed streamflow dataset (GSIM)
            # the data has been reorganized to a 1380 * 30961 matrix. 1380 is the month
            # number from 1901.1 to 2015.12. 30961 include two columns for year and month plus
            # streamflow at 30959 gauge locations reported by GSIM
            ref_mat_file = "{}/GSIM/GSIM_198601_199512.mat".format(
                parameter.reference_data_path.rstrip("/")
            )
            if parameter.short_ref_name != "":
                ref_name = parameter.short_ref_name
            elif parameter.reference_name != "":
                # parameter.ref_name is used to search though the reference data directories.
                # parameter.reference_name is printed above ref plots.
                ref_name = parameter.reference_name
            else:
                ref_name = "GSIM"
            parameter.ref_name_yrs = "{} ({}-{})".format(
                ref_name, parameter.ref_start_yr, parameter.ref_end_yr
            )
            ref_mat = scipy.io.loadmat(ref_mat_file)
            ref_array = ref_mat["GSIM"].astype(numpy.float64)
        if parameter.print_statements:
            # GSIM: 1380 x 30961
            # wrmflow: 720 x 360 x 360
            print("ref_array.shape={}".format(ref_array.shape))

        # Load E3SM simulated streamflow dataset
        if not using_test_mat_file:
            # `Dataset` will take the time slice from test_start_yr to test_end_yr
            test_data = utils.dataset.Dataset(parameter, test=True)
            parameter.test_name_yrs = utils.general.get_name_and_yrs(
                parameter, test_data
            )
            test_data_ts = test_data.get_timeseries_variable(var)
            var_array = test_data_ts(
                cdutil.region.domain(latitude=(-90.0, 90, "ccb"))
            )
            if parameter.print_statements:
                print(
                    "test var original dimensions={}".format(var_array.shape)
                )
            var_transposed = numpy.transpose(var_array, (2, 1, 0))
            if parameter.print_statements:
                print(
                    "test var transposed dimensions={}".format(
                        var_transposed.shape
                    )
                )
            test_array = var_transposed.astype(numpy.float64)
            areatotal2 = test_data.get_static_variable("areatotal2", var)
            area_upstream = numpy.transpose(areatotal2, (1, 0)).astype(
                numpy.float64
            )
            if parameter.print_statements:
                print(
                    "area_upstream dimensions={}".format(area_upstream.shape)
                )
        else:
            # This block is only for debugging -- i.e., when testing with a `mat` file.
            files_in_test_data_path = os.listdir(parameter.test_data_path)
            mat_files = list(
                filter(
                    lambda file_name: file_name.endswith(".mat"),
                    files_in_test_data_path,
                )
            )
            if len(mat_files) == 1:
                mat_file = mat_files[0]
            elif len(mat_files) > 1:
                raise Exception(
                    "More than one .mat file in parameter.test_data_path={}".format(
                        parameter.test_data_path
                    )
                )
            else:
                raise Exception(
                    "No .mat file in parameter.test_data_path={}".format(
                        parameter.test_data_path
                    )
                )
            test_mat_file = "{}/{}".format(
                parameter.test_data_path.rstrip("/"), mat_file
            )
            # FIXME: F524 '...'.format(...) is missing argument(s) for placeholder(s): 2
            parameter.test_name_yrs = "{} ({}-{})".format(  # noqa
                parameter.test_start_yr, parameter.test_end_yr
            )
            data_mat = scipy.io.loadmat(test_mat_file)
            if "E3SMflow" in data_mat.keys():
                # `edison` file uses this block
                e3sm_flow = data_mat["E3SMflow"]
                if parameter.print_statements:
                    print('e3sm_flow = data_mat["E3SMflow"]')
            else:
                # `test` file uses this block
                e3sm_flow = data_mat
                if parameter.print_statements:
                    print("e3sm_flow = data_mat")
            try:
                if e3sm_flow["uparea"].shape == (1, 1):
                    # `edison` file uses this block
                    area_upstream = e3sm_flow["uparea"][0][0].astype(
                        numpy.float64
                    )
                    if parameter.print_statements:
                        print(
                            'e3sm_flow["uparea"] was indexed into for later use'
                        )
                else:
                    area_upstream = e3sm_flow["uparea"].astype(numpy.float64)
                    if parameter.print_statements:
                        print('e3sm_flow["uparea"] will be used')
            except KeyError:
                # `test` file uses this block
                area_upstream = None
                if parameter.print_statements:
                    print(
                        "WARNING: uparea not found and will thus not be used"
                    )
            if e3sm_flow["wrmflow"].shape == (1, 1):
                # `edison` file uses this block
                test_array = e3sm_flow["wrmflow"][0][0].astype(numpy.float64)
                if parameter.print_statements:
                    print(
                        'e3sm_flow["wrmflow"] was indexed into for later use'
                    )
            else:
                # `test` file uses this block
                test_array = e3sm_flow["wrmflow"].astype(numpy.float64)
                if parameter.print_statements:
                    print('e3sm_flow["wrmflow"] will be used')
        if parameter.print_statements:
            # For edison: 720x360x600
            print("test_array.shape={}".format(test_array.shape))
        if type(area_upstream) == cdms2.tvariable.TransientVariable:
            area_upstream = area_upstream.getValue()

        # Resolution of MOSART output
        resolution = 0.5
        # Search radius (number of grids around the center point)
        radius = 1
        bins = numpy.floor(gauges[:, 7:9].astype(numpy.float64) / resolution)
        # Move the ref lat lon to grid center
        lat_lon = (bins + 0.5) * resolution
        if parameter.print_statements:
            print("lat_lon.shape={}".format(lat_lon.shape))

        # Define the export matrix
        # Annual mean of test, annual mean of ref, error for area, lat, lon
        export = numpy.zeros((lat_lon.shape[0], 9))
        if parameter.print_statements:
            print("export.shape={}".format(export.shape))
        for i in range(lat_lon.shape[0]):
            if parameter.print_statements and (i % 1000 == 0):
                print("On gauge #{}".format(i))
            if parameter.max_num_gauges and i >= parameter.max_num_gauges:
                break
            lat_ref = lat_lon[i, 1]
            lon_ref = lat_lon[i, 0]
            # Estimated drainage area (km^2) from ref
            area_ref = gauges[i, 13].astype(numpy.float64)

            if area_upstream is not None:
                drainage_area_error, lat_lon_ref = get_drainage_area_error(
                    radius,
                    resolution,
                    lon_ref,
                    lat_ref,
                    area_upstream,
                    area_ref,
                )
            else:
                # Use the center location
                lat_lon_ref = [lat_ref, lon_ref]
            if using_ref_mat_file:
                origin_id = gauges[i, 1].astype(numpy.int64)
                # Column 0 -- year
                # Column 1 -- month
                # Column origin_id + 1 -- the ref streamflow from gauge with the corresponding origin_id
                extracted = ref_array[:, [0, 1, origin_id + 1]]
                monthly_mean = numpy.zeros((12, 1)) + numpy.nan
                # For GSIM, shape is (1380,)
                month_array = extracted[:, 1]
                for month in range(12):
                    # Add 1 to month to account for the months being 1-indexed
                    month_array_boolean = month_array == month + 1
                    s = numpy.sum(month_array_boolean)
                    if s > 0:
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
                        monthly_mean[month] = numpy.nanmean(
                            extracted[month_array_boolean, 2]
                        )
                # This is ref annual mean streamflow
                annual_mean_ref = numpy.mean(monthly_mean)
            if using_ref_mat_file and numpy.isnan(annual_mean_ref):
                # All elements of row i will be nan
                export[i, :] = numpy.nan
            else:
                if using_ref_mat_file:
                    # Reshape extracted[:,2] into a 12 x ? matrix; -1 means to
                    # calculate the size of the missing dimension.
                    # Note that `numpy.reshape(extracted[:, 2], (12,-1))` will not work.
                    # We do need to go from (12n x 1) to (12 x n).
                    # `reshape` alone would make the first row [January of year 1, February of year 1,...]
                    # (i.e., 12 sequential rows with n entries)
                    # We actually want the first row to be [January of year 1, January of year 2,...]
                    # (i.e., n sequential columns with 12 entries)
                    # So, we use `reshape` to slice into n segments of length 12 and then we `transpose`.
                    mmat = numpy.transpose(
                        numpy.reshape(extracted[:, 2], (-1, 12))
                    )
                    mmat_id = numpy.sum(mmat, axis=0).transpose()
                    if numpy.sum(~numpy.isnan(mmat_id), axis=0) > 0:
                        # There's at least one year of record
                        monthly = mmat[:, ~numpy.isnan(mmat_id)]
                    else:
                        monthly = monthly_mean
                    seasonality_index_ref, peak_month_ref = get_seasonality(
                        monthly
                    )
                else:
                    ref_lon = int(
                        1
                        + (lat_lon_ref[1] - (-180 + resolution / 2))
                        / resolution
                    )
                    ref_lat = int(
                        1
                        + (lat_lon_ref[0] - (-90 + resolution / 2))
                        / resolution
                    )
                    ref = numpy.squeeze(ref_array[ref_lon - 1, ref_lat - 1, :])
                    # Note that `numpy.reshape(ref, (12,-1))` will not work.
                    # We do need to go from (12n x 1) to (12 x n).
                    # `reshape` alone would make the first row [January of year 1, February of year 1,...]
                    # (i.e., 12 sequential rows with n entries)
                    # We actually want the first row to be [January of year 1, January of year 2,...]
                    # (i.e., n sequential columns with 12 entries)
                    # So, we use `reshape` to slice into n segments of length 12 and then we `transpose`.
                    mmat = numpy.transpose(numpy.reshape(ref, (-1, 12)))
                    monthly_mean_ref = numpy.nanmean(mmat, axis=1)
                    annual_mean_ref = numpy.mean(monthly_mean_ref)
                    if numpy.isnan(annual_mean_ref) == 1:
                        # The identified grid is in the ocean
                        monthly = numpy.ones((12, 1))
                    else:
                        monthly = mmat

                    if type(monthly) == cdms2.tvariable.TransientVariable:
                        monthly = monthly.getValue()

                    seasonality_index_ref, peak_month_ref = get_seasonality(
                        monthly
                    )

                test_lon = int(
                    1 + (lat_lon_ref[1] - (-180 + resolution / 2)) / resolution
                )
                test_lat = int(
                    1 + (lat_lon_ref[0] - (-90 + resolution / 2)) / resolution
                )
                # For edison: 600x1
                test = numpy.squeeze(test_array[test_lon - 1, test_lat - 1, :])
                # For edison: 12x50
                # Note that `numpy.reshape(test, (12,-1))` will not work.
                # We do need to go from (12n x 1) to (12 x n).
                # `reshape` alone would make the first row [January of year 1, February of year 1,...]
                # (i.e., 12 sequential rows with n entries)
                # We actually want the first row to be [January of year 1, January of year 2,...]
                # (i.e., n sequential columns with 12 entries)
                # So, we use `reshape` to slice into n segments of length 12 and then we `transpose`.
                mmat = numpy.transpose(numpy.reshape(test, (-1, 12)))
                monthly_mean_test = numpy.nanmean(mmat, axis=1)
                annual_mean_test = numpy.mean(monthly_mean_test)
                if numpy.isnan(annual_mean_test) == 1:
                    # The identified grid is in the ocean
                    monthly = numpy.ones((12, 1))
                else:
                    monthly = mmat

                if type(monthly) == cdms2.tvariable.TransientVariable:
                    monthly = monthly.getValue()

                seasonality_index_test, peak_month_test = get_seasonality(
                    monthly
                )

                export[i, 0] = annual_mean_ref
                export[i, 1] = annual_mean_test
                if area_upstream is not None:
                    export[i, 2] = (
                        drainage_area_error * 100
                    )  # From fraction to percentage of the drainage area bias
                export[
                    i, 3
                ] = seasonality_index_ref  # Seasonality index of ref
                export[i, 4] = peak_month_ref  # Max flow month of ref
                export[
                    i, 5
                ] = seasonality_index_test  # Seasonality index of test
                export[i, 6] = peak_month_test  # Max flow month of test
                export[i, 7:9] = lat_lon_ref  # latlon of ref

        # Remove the gauges with nan flow
        # `export[:,0]` => get first column of export
        # `numpy.isnan(export[:,0])` => Boolean column, True if value in export[x,0] is nan
        # `export[numpy.isnan(export[:,0]),:]` => rows of `export` where the Boolean column was True
        # Gauges will thus only be plotted if they have a non-nan value for both test and ref.
        if parameter.print_statements:
            print(
                "export.shape before removing ref nan means={}".format(
                    export.shape
                )
            )
        export = export[~numpy.isnan(export[:, 0]), :]
        if parameter.print_statements:
            print(
                "export.shape before removing test nan means={}".format(
                    export.shape
                )
            )
        export = export[~numpy.isnan(export[:, 1]), :]
        if parameter.print_statements:
            print(
                "export.shape after both nan removals={}".format(export.shape)
            )

        if area_upstream is not None:
            # Set the max area error (percent) for all plots
            max_area_error = 20
            # `export[:,2]` gives the third column of `export`
            # `export[:,2]<=max_area_error` gives a Boolean array,
            # `True` if the value in the third column of `export` is `<= max_area_error`
            # `export[export[:,2]<=max_area_error,:]` is `export` with only the rows where the above is `True`.
            export = export[export[:, 2] <= max_area_error, :]
            if parameter.print_statements:
                print(
                    "export.shape after max_area_error cut={}".format(
                        export.shape
                    )
                )

        if parameter.print_statements:
            print("Variable: {}".format(var))

        if parameter.test_title == "":
            parameter.test_title = parameter.test_name_yrs
        if parameter.reference_title == "":
            parameter.reference_title = parameter.ref_name_yrs

        # Seasonality
        # Plot original ref and test, not regridded versions.
        plot_seasonality_map(export, parameter)

        # Bias between test and ref as a percentage
        # (Relative error as a percentage)
        # 100*((annual_mean_test - annual_mean_ref) / annual_mean_ref)
        bias = 100 * ((export[:, 1] - export[:, 0]) / export[:, 0])
        plot_annual_map(export, bias, parameter)

        # Scatterplot
        # These arrays will have fewer entries than the original `export` matrix
        # because of the nan removal steps.
        xs = export[:, 0]
        ys = export[:, 1]
        zs = export[:, 2]
        plot_annual_scatter(xs, ys, zs, parameter)

    return parameter
