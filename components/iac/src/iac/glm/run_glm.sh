#! /bin/csh

# Things to change before running GLM:
#    Change simulation options (nodata, future_rate, resolution, logging, priority, etc
#    Set future_run =0 or 1 (in addition to future_rate_option)
#    Change output directory
#    Set casename
#    In GLM.c change the bestcase flags
#    In GLM.c change the NX and NY values
#    In GLM.FUTURE.c change the bestcase flags

set fmain_dir = "/data/"   #root directory

# set future_run to 0 to run the historical time period only and set to 1 to run the future time period only
set future_run = 1

# 1=no urban, 2=urban
foreach urban_option (1)

# set codes that will be used based on settings above
if($urban_option == 1) then
   set progmain = glm.c
   set progfile1 = land_t1-t5.c
   cp $progmain $progfile1
   set progfile2 = land_t6.c
   set progmain2 = glm.future.c
   cp $progmain2 $progfile2
endif

if($urban_option == 2) then
   set progmain = glm_urban.c
   set progfile1 = land_t1-t5.c
   cp $progmain $progfile1
   set progfile2 = land_t6.c
   set progmain2 = glm_urban.future.c
   cp $progmain2 $progfile2
endif

# only half-degree resolution is supported in this version
foreach res_option (2) # 1=1deg resolution, 2=0.5degree resolution

# choose future scenario and input data to use
foreach future_rate_option (3)
# 0 = no future
# 3 = GCAM
# 4 = AIM
# 5 = IMAGE
# 6 = MESSAGE

# choose historical land-use dataset. 5=HYDE 3, 6="No-Data"
foreach nodata_option (5)

# choose priority for clearing and wood harvest. 1=primary, 2=secondary
foreach smart_flow_option (1)
set harvest_option = $smart_flow_option

# choose agricultural residence option. 1=minimum flows only, 5=shifting cultivation within the locations defined by our SC map
foreach adjust_smart_flow_option (1)

# choose whether clearing for agriculture is counted towards meting wood harvest demand (2) or not (1)
foreach converted_forest_land_option (1)

# set option for algorithm for spatial allocation of wood harvest
# only option 1 is supported in this version
foreach zdis_option (1)

# set wood harvest dataset to use
foreach logging_option (1) 
#0 1 4)
# 0=wh=zero,1=standard wh data, 4="nodata"

# set historical time period to simulate (or to build future simulation off)
foreach hist_option (2)
# 1=1700-2000, 2=1500-2000, 3=1850-2005, 4=no historical simulation

set foutput_dir = "data/glm/output/iESM/Expt1/"

if($future_rate_option == 0) set scenario_id = 0
if($future_rate_option == 3) then
   set scenario_id = GCAM
   set num_reg = 14
   set gridded_wh = 2
endif
if($future_rate_option == 4) then
   set scenario_id = AIM
   set num_reg = 24
   set gridded_wh = 0
endif
if($future_rate_option == 5) then
   set scenario_id = IMAGE
   set num_reg = 24
   set gridded_wh = 0
endif
if($future_rate_option == 6) then
   set scenario_id = MESSAGE
   set num_reg = 24
   set gridded_wh = 0
endif

## map glm settings to "model factor" ids
set luh = $nodata_option
set lut = $adjust_smart_flow_option  ## 1=min flows,2=min flows everywhere except 23.5N to 23.5S,
                                     ## 3=min flows everywhere except region boundaries
set lup = $harvest_option  ## 1=v priority, 2=s priority, 3=s for transitions, v for harvest
set luw = $converted_forest_land_option
set luz = $zdis_option
set luf = $future_rate_option
set lul = $logging_option
set lur = $res_option
set lud = $hist_option
set luu = $urban_option

#set casename = h${luh}t${lut}p${lup}w${luw}z${luz}f${luf}l${lul}r${lur}d${lud}u${luu}_${scenario_id}
#set casename = h${luh}t${lut}p${lup}w${luw}z${luz}f${luf}l${lul}r${lur}d${lud}u${luu}
set casename = GCAM_AEZ__noSC
set hist_data = data/glm/output/LUHa.rc2/LUHa.rc2/
#set hist_data = data/glm/output/test_public_glm/h${luh}t${lut}p${lup}w${luw}z${luz}f0l${lul}r${lur}d${lud}u${luu}/
if($hist_option == 4) set hist_data = dummy
echo $hist_data
echo $casename

rm -f *.test
rm -f wh.*
rm -r -f ${fmain_dir}${foutput_dir}${casename}
rm -r -f ${fmain_dir}${foutput_dir}updated_states
rm -r -f ${fmain_dir}${foutput_dir}lu
rm -r -f ${fmain_dir}${foutput_dir}tester
mkdir ${fmain_dir}${foutput_dir}${casename}
mkdir ${fmain_dir}${foutput_dir}updated_states
mkdir ${fmain_dir}${foutput_dir}lu
mkdir ${fmain_dir}${foutput_dir}tester

if($future_run == 0) then
   gfortran -o country.exe country.option.f90
   ./country.exe $smart_flow_option $converted_forest_land_option $zdis_option $adjust_smart_flow_option $harvest_option $nodata_option

   # past runs
   foreach trun ("tone" "ttwo" "tthree" "tfour" "tfive")
      gcc -o test.exe $progfile1 -lm -fpic
      ./test.exe $trun $future_rate_option $nodata_option $logging_option $foutput_dir $res_option $hist_option
   end

endif

# future runs
if($future_run == 1) then
   gfortran -o regional.exe regional.option.f90
   ./regional.exe $smart_flow_option $converted_forest_land_option $zdis_option $adjust_smart_flow_option $harvest_option $num_reg

   foreach trun ("tsix" "tseven")
     if($trun == "tsix") set time_period = 46 
     if($trun == "tseven") set time_period = 51      
     gcc -o test2.exe $progfile2 -lm -fpic
     ./test2.exe $trun $future_rate_option $nodata_option $logging_option $foutput_dir $res_option $time_period $hist_option $gridded_wh $hist_data
   end 

endif

########################

mv ${fmain_dir}${foutput_dir}updated_states ${fmain_dir}${foutput_dir}${casename}
mv ${fmain_dir}${foutput_dir}lu ${fmain_dir}${foutput_dir}${casename}
mv ${fmain_dir}${foutput_dir}tester ${fmain_dir}${foutput_dir}${casename}


# global
mv global.test ${fmain_dir}${foutput_dir}${casename}/tester/global.${casename}.txt
mv global.future.test ${fmain_dir}${foutput_dir}${casename}/tester/global.future.${casename}.txt

foreach country (us br cn ak aus)
  mv $country.country.test ${fmain_dir}${foutput_dir}${casename}/tester/${country}.country.${casename}.txt
  mv $country.future.test ${fmain_dir}${foutput_dir}${casename}/tester/${country}.future.${casename}.txt
end 

mv global.primeflow.txt ${fmain_dir}${foutput_dir}${casename}/tester/global.primeflow.${casename}.txt
mv country.primeflow.txt ${fmain_dir}${foutput_dir}${casename}/tester/country.primeflow.${casename}.txt

foreach continent (na sa eu asia africa aus)
  mv $continent.continent.test ${fmain_dir}${foutput_dir}${casename}/tester/${continent}.continent.${casename}.txt
  mv $continent.primeflow.txt ${fmain_dir}${foutput_dir}${casename}/tester/${continent}.primeflow.${casename}.txt
end


foreach regional (canada USA_NO_Alaska Alaska Japan former_USSR central_america south_america north_africa west_africa east_africa south_africa OECD_Europe central_eastern_europe middle_east south_asia east_asia southeast_asia oceania)
  mv $regional.regional.test ${fmain_dir}${foutput_dir}${casename}/tester/${regional}.regional.${casename}.txt
end

foreach ipcc_region (oecd90 ref alm asia)
  mv ${ipcc_region}.future.test ${fmain_dir}${foutput_dir}${casename}/tester/${ipcc_region}.future.${casename}.txt
end

mv wh.* ${fmain_dir}${foutput_dir}${casename}/tester/.
mv country.final.stats.txt ${fmain_dir}${foutput_dir}${casename}/tester/.

endif


end
end
end
end
end
end
end 

rm -f $progfile1 $progfile2

exit
