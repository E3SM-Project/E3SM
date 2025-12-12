# create_e3sm_gcam_land_scalar_baseline_local.r

# This function reads monthly elm h1 files to obtain pft-level values for baseline of:
#   veg_cf%npp, col_cf%hr, veg_pp%wtgcell for the veg landunit only, and total cell area

# the period averaged data in E3SM-GCAM include zero values (and not uninitialized or large missing data values)
#    because only active pfts are updated to be set to non-zero values and sent through the coupler
# the coupler does not distinguish which values to include in the total average
# the passed pft weights are also zero when the pft is not active, 
#    which reduces the influence of averages that include a pft that is active for only part of a period
# so need to include missing data in the period average as zero values

# It is necessary to use the total period average because the max monthly values do not line up across variables,
#    and when matched to npp max the hr is far from a max

# the current E3SM-GCAM period is 5-years, so use this for the baseline calculation
# even if the GCAM period decrdeases to 1, we may still want a 5 year average for a baseline
# the calculation period is (year_end - year_start + 1)

# Note that the output wtgcell here is the weight of the full grid cell,
#    while veg_pp%wtgcell in the h1 file is the weight of pft on land that may not cover the entire cell
#    because elm assumes that land covers the entire cell

# The monthly values for npp (gC/m^2/s) and hr (gC/m^2/s) and pft weight (fraction) are averages of 30-minute elm calcs

# The default outputs of this script are the average values over the period calculated as:
#    the month-weighted average (calculated via days) of the monthly averages across years
#       this just splits the period average into two steps, reducing the number of weights
#       and possibly reducing the precision error as the weights for 12 months are larger than the weights for all months
#    this matches how it is calculated in E3SM, wich is a total average over the period, based on 30-minute increments
#    note that any missing monthly input data is included in the calculations as zero values


# Don't use these, they have been deprecated, but the code is still here for diagnostic purposes:
# Alternatively, the output can be values corresponding to the month of max monthly npp
# Alternatively, the max monthly value of each can be output for comparison

# The full grid cell area is also output (km^2), only for land cells;
#    the elm h1 file has missing values for non-land cells

# if using a sequence of multiple years, the monthly average across the years will be calculated first,
#   then the annual average will be calculated or the maximum monthly value will be selected

# arguments

# indir:				path to h1 files (requires final "/")
# case_name:			simulation case name that is the main part of the history file name including ".elm.h#."; "yyyy-mm.nc" will be appended as appropriate
# year_start:		start year of h1 files to read; will also be included in output file names
# year_end:			end year of h1 files to read; will also be included in output file names
# outdir:			path to output files (requires final "/")
# out_base_name:    the beginning of the output file names; needs to include E3SM resolution and active component cofig info
# avg_period:       Default=5; number of years in averaging period

# outputs

# <out_base_name>_<avg_tag>_<year_start>-<year end>_hr.csv
# <out_base_name>_<avg_tag>_<year_start>-<year end>_npp.csv
# <out_base_name>_<avg_tag>_<year_start>-<year end>_pft_wt.csv
# <out_base_name>_cell_area.csv

# output file format

# csv with one header file, separator = ","
# four columns for pft data (npp, weight, hr), in column order: pft, lon, lat, value
# three columns for area data, in column order: lon, lat, value
# data are in numerical order sorted according to:
#   lon varies fastest
#   lat varies second fastest
#   pft varies third fastest (if present)
# all pftXlonXlat combinations are needed in the output files

# NOTES

# about 1.5 hours to finish on desktop

# resolution will be determined from the h1 file
# recall that longitude varies fastest
# this sript is agnostic regarding resolution; indices are determined based on input file order

# for f09 and f19:
# elm grid cell edges start at lon=-<half-lon-res> (zero center) and lat=-90; pole cells are half-lat-res (pole full-cell center)

# for r05 and r025 and r0125:
# elm grid cell edges start at lon = -180 lat = -90, with equal spacing between cell centers

# lon-lat values are indices starting at one, increasing positively from the origin defined above
# not all pftXlonXlat combinations are included in the h1 1d arrays

# there is only one time record in each file
# all twelve months must be present for each year

# HR is for the whole column, so all pfts in the veg land unit column have the same HR
#   vegetated land unit is in the first index, with id = 1

# NPP can be negative! or missing value/NA or zero
# HR and pft weight are all positive or zero or missing value/NA (pft area only has no missing values)
# cell area is always positive and non-zero, na/missing value, in the h1 files

# Currently this assumes that there is only one topounit per grid cell and that only veg landunit pfts are mapped
#   this mapping in E3SM does not include topounits or non-veg land units

# this updated calculation, which better matches what is now calculated in E3SM-GCAM gives slightly different results from the previous calculation
# for one example month in one year:
#    npp values have 238 instances of absolute differences greater than 10% (out of 81159 valid values)
#    pft values have max perecent absolute difference of 0.48% (out of 81159 valid values)
#    hr values have 272 instances of absolute differences greater than 1.8% (out of 15349 valid values)

# this code operates only on the veg landunit

# there should be 17 pfts for each grid cellXtopounit
#   pfts in order are:
#      0 - bare
#      1 - needle leaf evergreen temperate tree
#      2 - needle leaf evergreen boreal tree
#      3 - needle leaf deciduous boreal tree
#      4 - broad leaf evergreen tropical tree
#      5 - broad leaf evergreen temperate tree
#      6 - broad leaf deciduaous tropical tree
#      7 - broad leaf deciduaous temperate tree
#      8 - broad leaf deciduaous boreal tree
#      9 - broad leaf evergreen temperate shrub
#      10 - broad leaf deciduaous temperate shrub
#      11 - broad leaf deciduaous boreal shrub
#      12 - C3 arctic grass
#      13 - C3 non-arctic grass
#      14 - C4 grass
#      15 - crop
#      16 - NA


library(ncdf4)

create_e3sm_gcam_land_scalar_baseline_local <- function( indir = "./",
											  case_name,
											  year_start = 2010,
											  year_end = 2014,
											  outdir = "./",
											  out_base_name,
											  avg_period = 5) {
	
	cat("Start create_baseline_for_iesm_scalars.r at", date(), "\n")
	
	veg_lunit_id = 1
	num_years = year_end - year_start + 1

    # stop if start/end years do not match the averaging period
    if(num_years != avg_period) {
    	cat("Number of input years (num_years = year_end - year_start + 1) do not match the averaging period (avg_period): ", num_years, "!=", avg_period, "\n")
		stop()
    }
	
	avg_tag = "_PerAvg"
	
	# error tolerance for comparing weighted avering methods for precision
	# this tolerance shows now differences in all three variables for the test case
	# note that precision calc errors are larger for pft weight than npp and hr
	err_tol = 1e-14
	
	# no leap years in the model
	days_in_months = c(31,28,31,30,31,30,31,31,30,31,30,31)
	days_in_year = sum(days_in_months)
	all_weights = array(dim=c(avg_period, 12))
	
    # open one file to get the constant info
    
    if(year_start < 10){
		ystr = paste0("000", year_start)
	} else if(year_start < 100) {
	   ystr = paste0("00", year_start)
	} else if(year_start < 1000) {
	   ystr = paste0("0", year_start)
	} else {
	   ystr = paste0(year_start)
	}
    
    nid = nc_open(paste0(indir, case_name, ystr, "-01.nc"))
	num_lon = nid$dim$lon$len
	num_lat = nid$dim$lat$len
	num_cells = num_lon * num_lat
	num_col_ind = nid$dim$column$len
	num_pft_ind = nid$dim$pft$len
	num_pft = nid$dim$natpft$len
	num_lunit = nid$dim$ltype$len
	num_gridcell = nid$dim$gridcell$len
	num_topounit = nid$dim$topounit$len
	lon = ncvar_get(nid,varid="lon",start=c(1), count=c(num_lon))
	lat = ncvar_get(nid,varid="lat",start=c(1), count=c(num_lat))
	# these are lon, lat variables
	area = ncvar_get(nid,varid="area")
	landfrac = ncvar_get(nid,varid="landfrac") # this should match landmask
	pftmask = ncvar_get(nid,varid="pftmask") # not all landfrac > 0 have valid pfts
	# these are column index variables
	col_lunit = ncvar_get(nid,varid="cols1d_itype_lunit")
	col_lon_ind = ncvar_get(nid,varid="cols1d_ixy")
	col_lat_ind = ncvar_get(nid,varid="cols1d_jxy")
	col_topounit = ncvar_get(nid,varid="cols1d_topounit") # this is topounit index
	# these are pft index variables
	pft_pft = ncvar_get(nid,varid="pfts1d_itype_veg")
	pft_lon_ind = ncvar_get(nid,varid="pfts1d_ixy")
	pft_lat_ind = ncvar_get(nid,varid="pfts1d_jxy")
	pft_topounit_ind = ncvar_get(nid,varid="pfts1d_topounit") # this is topounit index
	pft_lunit = ncvar_get(nid,varid="pfts1d_itype_lunit")
	pft_active = ncvar_get(nid,varid="pfts1d_active") # this includes pfts on other land units and pft wts > 0 on veg land unit
	nc_close(nid)    
	
	# check to make sure that there is only one topounit per grid cell
	if (num_gridcell != num_topounit) {
		stop("Error: This code does not currently support mapping for multiple topounits within a grid cell\n")
	}
	
	# set up some arrays
	avg_npp = array(dim=c(num_pft_ind))
	avg_hr = array(dim=c(num_col_ind))
	avg_pft_weight = array(dim=c(num_pft_ind))
	pft_cell_weight = array(dim=c(num_pft_ind))
	hr_monthly_avg = array(dim=c(12, num_col_ind))
	npp_monthly_avg = array(dim=c(12, num_pft_ind))
	pft_weight_monthly_avg = array(dim=c(12, num_pft_ind))
	hr_monthly_cnt = array(dim=c(12, num_col_ind))
	npp_monthly_cnt = array(dim=c(12, num_pft_ind))
	pft_weight_monthly_cnt = array(dim=c(12, num_pft_ind))
	
	# these are for period average
	npp_month_values = array(dim=c(num_years, 12, num_pft_ind))
	pft_weight_month_values = array(dim=c(num_years, 12, num_pft_ind))
	hr_month_values = array(dim=c(num_years, 12, num_col_ind))
	npp_month_values[,,] = 0.0
	pft_weight_month_values[,,] = 0.0
	hr_month_values[,,] = 0.0
	
	avg_npp[] = 0
	avg_hr[] = 0
	avg_pft_weight[] = 0
	
	hr_monthly_avg[,] = 0
	npp_monthly_avg[,] = 0
	pft_weight_monthly_avg[,] = 0
	hr_monthly_cnt[,] = 0
	npp_monthly_cnt[,] = 0
	pft_weight_monthly_cnt[,] = 0
	
	# loop over 12 months
	for (m in 1:12) {
		if (m < 10) {mtag = paste0(0,m)
		} else { mtag = paste0(m)}
		
		# loop over years to sum the monthly values
		# assume that all values are valid here, unless pft weight is zero
		for (y in year_start:year_end) {
			cat("Processing month", m, "year", y, "\n")
			yind = y - year_start + 1
			
			if(y < 10){
				ystr = paste0("000", y)
			} else if(y < 100) {
			   ystr = paste0("00", y)
			} else if(y < 1000) {
			   ystr = paste0("0", y)
			} else {
			   ystr = paste0(y)
			}
			
			fname = paste0(indir, case_name, ystr, "-", mtag, ".nc")
			nid = nc_open(fname)
			
			# NPP and pft weight are by pft index
			pft_land_weight = ncvar_get(nid,varid="pfts1d_wtgcell")
			npp = ncvar_get(nid,varid="NPP")
			
			# HR is by column index (and with multiple landunits)
			hr = ncvar_get(nid,varid="HR")
			
			# convert h1 pft frac of land to frac of grid cell
			pft_cell_weight[] = 0
			for (p in 1: num_pft_ind) {
				pft_cell_weight[p] = pft_land_weight[p] * landfrac[pft_lon_ind[p], pft_lat_ind[p]] * pftmask[pft_lon_ind[p], pft_lat_ind[p]]
			}
			
			# since zero-weight pfts are used to determine whether a cell has valid values for processing,
			#   zeroing the pft weights for undesired records effectively filters them out
			
			# zero out the non-veg landunit pft weights and npp values to be able to match with hr uniquely
			nvlu_inds = which(pft_lunit != veg_lunit_id)
			pft_cell_weight[nvlu_inds] = 0
			npp[nvlu_inds] = 0
			
			# if pft weight is not zero, add it into the sum and count it
			# need to figure out the hr index - but at the column level; so include if any pft is non-zero in a cell
			# this also removes any NA/missing pft weight values (of which there are none)
			pft_nonzero_inds = which(pft_cell_weight > 0)
			
			# store the sum for period avg
			pft_weight_monthly_avg[m, pft_nonzero_inds] = pft_weight_monthly_avg[m, pft_nonzero_inds] + pft_cell_weight[pft_nonzero_inds]
			pft_weight_monthly_cnt[m, pft_nonzero_inds] = pft_weight_monthly_cnt[m, pft_nonzero_inds] + 1
			
			# do not count missing values
			na_inds = which(is.na(npp))
			avg_inds = setdiff(pft_nonzero_inds, na_inds)
			
			# store the sum for period avg
			npp_monthly_avg[m, avg_inds] = npp_monthly_avg[m, avg_inds] + npp[avg_inds]
			npp_monthly_cnt[m, avg_inds] = npp_monthly_cnt[m, avg_inds] + 1

			# find the hr columns that are in cells with at least one non-zero pft
			# only process veg land unit
			# use data frames to only do a cell once
			
			# get unique valid cells/topounits - unique essentially selects one pft record for each cell/topounit
			valid_pft_cell = data.frame(valid_pft_lon = pft_lon_ind[pft_nonzero_inds],
										valid_pft_lat = pft_lat_ind[pft_nonzero_inds],
										valid_pft_topounit = pft_topounit_ind[pft_nonzero_inds],
										valid_pft_lunit = pft_lunit[pft_nonzero_inds])
			valid_pft_cell = unique(valid_pft_cell)
			# get unique hr cells for veg land unit
			hr_veglu = data.frame(hr_lon_ind = col_lon_ind, hr_lat_ind = col_lat_ind, hr_topounit = col_topounit,
									hr_lu = col_lunit, hr_value = hr, hr_col_ind = c(1:num_col_ind))
			hr_veglu  = hr_veglu[hr_veglu$hr_lu == veg_lunit_id,]
			# this unique call shouldn't be necessary, but do it anyway
			hr_veglu = hr_veglu[!duplicated(hr_veglu[,c(1:4)]),]
			# merge these on the cell indices, topounit, and landunit
			valid_pft_hr = merge(valid_pft_cell, hr_veglu, by.x = c("valid_pft_lon", "valid_pft_lat", "valid_pft_topounit", "valid_pft_lunit"),
															by.y = c("hr_lon_ind", "hr_lat_ind", "hr_topounit", "hr_lu"), all.x = TRUE)
			# drop any NA/missing values
			valid_pft_hr = valid_pft_hr[!is.na(valid_pft_hr$hr_value),]
			
			# store the sum for period avg
			hr_monthly_avg[m, valid_pft_hr$hr_col_ind] = hr_monthly_avg[m, valid_pft_hr$hr_col_ind] + valid_pft_hr$hr_value
			hr_monthly_cnt[m, valid_pft_hr$hr_col_ind] = hr_monthly_cnt[m, valid_pft_hr$hr_col_ind] + 1
			
			nc_close(nid)
			
			# store all the valid values for calculating average differently
			pft_weight_month_values[yind, m, pft_nonzero_inds] = pft_cell_weight[pft_nonzero_inds]
			npp_month_values[yind, m, avg_inds] = npp[avg_inds]
			hr_month_values[yind, m, valid_pft_hr$hr_col_ind] = valid_pft_hr$hr_value
			all_weights[yind, m] = days_in_months[m] / (avg_period * sum(days_in_months))
			
		} # end y loop over year
		
		# now calc the average monthly values; equal weights cuz month length is the same across years
		# doing this first makes the weighted calc easier
		# use the avg_period instead of the count so that the zero values are included in the averages
		pft_weight_monthly_avg[m,] = pft_weight_monthly_avg[m,] / avg_period
		npp_monthly_avg[m,] = npp_monthly_avg[m,] / avg_period
		hr_monthly_avg[m,] = hr_monthly_avg[m,] / avg_period
		
	} # end m loop over month
	
	# calc period average
	# in E3SM coupler values are summed and averaged without error checking, and zeros are passed for missing/non-active pft data
	# in E3SM outlier npp/hr values are removed before scalar calculation (outliers are determined without zero and nan values)
	# use weights based on days in month
	
	# notify of missing monthly averages
	# but just sum them all up because this includes all pfts in all cells, many of which do not exist
	pft_monthly_zero_inds = NULL
	npp_monthly_zero_inds = NULL
	hr_monthly_zero_inds = NULL
	
	# npp and pft weight
	for(v in 1: num_pft_ind) {
		# pft weight
		# log the non-zero and zero monthly values
		cnt_mask_pft = pft_weight_monthly_cnt[,v]
		pft_monthly_zero_inds = c(pft_monthly_zero_inds, which(cnt_mask_pft == 0))
		# use all weights for all months to include the zero values in the averages
		day_weights_pft = (days_in_months) / sum(days_in_months)
		avg_pft_weight[v] = sum(pft_weight_monthly_avg[,v] * day_weights_pft)
	
		# npp
		# log the non-zero and zero monthly values
		cnt_mask_npp = npp_monthly_cnt[,v]
		npp_monthly_zero_inds = c(npp_monthly_zero_inds, which(cnt_mask_npp == 0))
		# use all weights for all months to include the zero values in the averages
		day_weights_npp = (days_in_months) / sum(days_in_months)
		avg_npp[v] = sum(npp_monthly_avg[,v] * day_weights_npp)
	}
	# hr
	for(v in 1: num_col_ind) {
		# log the non-zero and zero monthly values
		cnt_mask_hr = hr_monthly_cnt[,v]
		hr_monthly_zero_inds = c(hr_monthly_zero_inds, which(cnt_mask_hr == 0))
		# use all weights for all months to include the zero values in the averages
		day_weights_hr = (days_in_months) / sum(days_in_months)
		avg_hr[v] = sum(hr_monthly_avg[,v] * day_weights_hr)
	}
	
	cat("These missing indices include where the pfts do not exist\n")
	if(length(pft_monthly_zero_inds) > 0){
		cat("pft index has", length(pft_monthly_zero_inds), "missing monthly averages\n")
	}
	if(length(npp_monthly_zero_inds) > 0){
		cat("npp index has", length(npp_monthly_zero_inds), "missing monthly averages\n")
	}
	if(length(hr_monthly_zero_inds) > 0){
		cat("hr index has", length(hr_monthly_zero_inds), "missing monthly averages\n")
	}
	
	# output the npp, hr, and pft weight data
	# use a data frame to organize the data, as they are all output by pft, lon, lat
	out_pft_cell = data.frame(pft_id = pft_pft, lon_ind = pft_lon_ind, lat_ind = pft_lat_ind, topounit_ind = pft_topounit_ind, landunit_id = pft_lunit,
								npp_gC_per_m2_per_s = avg_npp, pft_wt_cell_frac = avg_pft_weight)
	# select only the veg landunit
	out_pft_cell = out_pft_cell[out_pft_cell$landunit_id == veg_lunit_id,]
	# get unique hr cells for veg land unit
	out_hr_veglu = data.frame(lon_ind = col_lon_ind, lat_ind = col_lat_ind, topounit_ind = col_topounit, landunit_id = col_lunit, hr_gC_per_m2_per_s = avg_hr)
	out_hr_veglu  = out_hr_veglu[out_hr_veglu$landunit_id == veg_lunit_id,]
	# again, this shouldn't be necessary
	out_hr_veglu = out_hr_veglu[!duplicated(out_hr_veglu[,c(1:4)]),]
	# merge these on the cell indices
	out_pft_hr = merge(out_pft_cell, out_hr_veglu, by = c("lon_ind", "lat_ind", "topounit_ind", "landunit_id"), all.x = TRUE)

	# check for bad values and change them to zero
	# there may not be any bad values cuz the na/NaN values are not included in the weighted average
	#    any missing data would have zero values
	
	# npp
	npp_inf_inds = which((out_pft_hr$npp_gC_per_m2_per_s == Inf) == TRUE)
	npp_neg_inf_inds = which((out_pft_hr$npp_gC_per_m2_per_s == -Inf) == TRUE)
	npp_na_inds = which(is.na(out_pft_hr$npp_gC_per_m2_per_s) == TRUE) # this includes NaN
	npp_bad_inds = c(npp_inf_inds, npp_neg_inf_inds, npp_na_inds)
	if (length(npp_bad_inds) > 0) { 
		cat("Warning: some bad values have been set to zero\n")
		cat("NPP # of bad inds is", length(npp_bad_inds), "\n")
		out_pft_hr$npp_gC_per_m2_per_s[npp_bad_inds] = 0
	}
	
	# also check for negative npp
	# negative values are passed cuz they can contribute to positive averages later when aggregated to gcam regions/types
	# the final calc scalar code filters out negative npp base values after aggregation because GCAM has only positive yields/accumulation
	npp_neg_inds = which((out_pft_hr$npp_gC_per_m2_per_s < 0) == TRUE)
	if (length(npp_neg_inds) > 0) { 
		cat("NPP # of negative inds is", length(npp_neg_inds), "\n")
	}
	
	# pft weight
	pftwt_inf_inds = which((out_pft_hr$pft_wt_cell_frac == Inf) == TRUE)
	pftwt_neg_inf_inds = which((out_pft_hr$pft_wt_cell_frac == -Inf) == TRUE)
	pftwt_na_inds = which(is.na(out_pft_hr$pft_wt_cell_frac) == TRUE) # this includes NaN
	pftwt_bad_inds = c(pftwt_inf_inds, pftwt_neg_inf_inds, pftwt_na_inds)
	if (length(pftwt_bad_inds) > 0) { 
		cat("Warning: some bad values have been set to zero\n")
		cat("pft wt # of bad inds is", length(pftwt_bad_inds), "\n")
		out_pft_hr$pft_wt_cell_frac[pftwt_bad_inds] = 0
	}
	
	# hr
	hr_inf_inds = which((out_pft_hr$hr_gC_per_m2_per_s == Inf) == TRUE)
	hr_neg_inf_inds = which((out_pft_hr$hr_gC_per_m2_per_s == -Inf) == TRUE)
	hr_na_inds = which(is.na(out_pft_hr$hr_gC_per_m2_per_s) == TRUE) # this includes NaN
	hr_bad_inds = c(hr_inf_inds, hr_neg_inf_inds, hr_na_inds)
	if (length(hr_bad_inds) > 0) { 
		cat("Warning: some bad values have been set to zero\n")
		cat("HR # of bad inds is", length(hr_bad_inds), "\n")
		out_pft_hr$hr_gC_per_m2_per_s[hr_bad_inds] = 0
	}

	out_pft_hr = out_pft_hr[order(out_pft_hr$pft_id, out_pft_hr$lat_ind, out_pft_hr$lon_ind, out_pft_hr$topounit_ind, out_pft_hr$landunit_id),]
	
	# make complete pftXlonXlat df for writing output files
	
	# this sets up the full record list
	pft_id_out = NULL
	lon_ind_out = NULL
	lat_ind_out = NULL
	area_out = NULL
	for (p in 0:(num_pft-1)) {
		for (t in 1:num_lat) {
			
			pft_id_out = c(pft_id_out, rep(p,num_lon))
			lon_ind_out = c(lon_ind_out, 1:num_lon)
			lat_ind_out = c(lat_ind_out, rep(t,num_lon))
			# just fill each pft grid here for completeness
			area_out = c(area_out, area[,t])
			
		} # end for t loop over lat
	} # end for loop over pft
	
	# now merge the data into the complete df and set all NA values to zero
	out_df = data.frame(pft_id = pft_id_out, lon_ind = lon_ind_out, lat_ind = lat_ind_out, cell_area_km2 = area_out)
	out_df = merge(out_df, out_pft_hr, by = c("pft_id", "lon_ind","lat_ind"), all.x = TRUE)
	out_df[is.na(out_df)] = 0
	
	# now write the files
	
	# the output file col names should be, in order, "pft_id", "lon_ind", "lat_ind", "value"
	# with lon varying fastest, and pft varying slowest
	
	npp_out_name = paste0(outdir, out_base_name, avg_tag, "_", year_start, "-", year_end, "_npp.csv")
	pft_wt_out_name = paste0(outdir, out_base_name, avg_tag, "_", year_start, "-", year_end, "_pft_wt.csv")
	hr_out_name = paste0(outdir, out_base_name, avg_tag, "_", year_start, "-", year_end, "_hr.csv")
	area_out_name = paste0(outdir, out_base_name, "_cell_area.csv")
	
	# npp
	write_df = out_df[,c("pft_id", "lon_ind", "lat_ind", "npp_gC_per_m2_per_s")]
	write_df = write_df[order(write_df$pft_id, write_df$lat_ind, write_df$lon_ind),]
	write.csv(write_df, file = npp_out_name, row.names=FALSE)
	
	# pft weight
	write_df = out_df[,c("pft_id", "lon_ind", "lat_ind", "pft_wt_cell_frac")]
	write_df = write_df[order(write_df$pft_id, write_df$lat_ind, write_df$lon_ind),]
	write.csv(write_df, file = pft_wt_out_name, row.names=FALSE)
	
	# hr
	write_df = out_df[,c("pft_id", "lon_ind", "lat_ind", "hr_gC_per_m2_per_s")]
	write_df = write_df[order(write_df$pft_id, write_df$lat_ind, write_df$lon_ind),]
	write.csv(write_df, file = hr_out_name, row.names=FALSE)
	
	# cell area
	write_df = out_df[out_df$pft_id == 0,c("lon_ind", "lat_ind", "cell_area_km2")]
	write_df = write_df[order(write_df$lat_ind, write_df$lon_ind),]
	write.csv(write_df, file = area_out_name, row.names=FALSE)
	
	cat("Finish create_baseline_for_iesm_scalars.r at", date(), "\n")
	
}