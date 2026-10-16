include(ScreamUtils)

if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  CreateUnitTest(mam_photo_table_test
    SOURCES mam_photo_table_test.cpp
    LIBS mam eamxx_scorpio_interface eamxx_io ekat::YamlParser
    LABELS "mam;physics"
    MPI_RANKS 1
    THREADS 1
  )
  target_compile_definitions(mam_photo_table_test PRIVATE
    SCREAM_DATA_DIR="${SCREAM_DATA_DIR}"
  )
  # Ensure test input files are present in the data dir
  GetInputFile(scream/mam4xx/photolysis/table_photo_input_ts_355.yaml)
  GetInputFile(scream/mam4xx/photolysis/jlong_input_ts_355.yaml)

  # The controlled surface-emission fixture has 72 levels, like the official
  # ne2np4 standalone case. Temporary NetCDF copies live under /tmp.
  if (SCREAM_NUM_VERTICAL_LEV EQUAL 72)
    find_library(MAM_TEST_NETCDF_C netcdf HINTS ${NetCDF_C_PATH}
      PATH_SUFFIXES lib lib64)
    if (NOT MAM_TEST_NETCDF_C)
      message(FATAL_ERROR "mam_surface_emissions_test requires the NetCDF C library")
    endif()
    CreateUnitTest(mam_surface_emissions_test
      SOURCES mam_surface_emissions_test.cpp
      LIBS mam eamxx_scorpio_interface eamxx_io ${MAM_TEST_NETCDF_C}
      LABELS "mam;physics;emissions"
      MPI_RANKS 1
      THREADS 1
    )
    target_include_directories(mam_surface_emissions_test PRIVATE ${NetCDF_C_PATH}/include)
    target_compile_definitions(mam_surface_emissions_test PRIVATE
      SCREAM_DATA_DIR="${SCREAM_DATA_DIR}"
    )
    set(emission_input_dir scream/mam4xx/emissions/ne2np4)
    set(surface_input_dir ${emission_input_dir}/surface)
    GetInputFile(${emission_input_dir}/dst_ne2np4_c20241028.nc)
    set(marine_organics_file
      monthly_macromolecules_0.1deg_bilinear_year01_merge_ne2np4_c20260807.nc)
    GetInputFile(${emission_input_dir}/${marine_organics_file})
    set(dms_file
      DMSflux.2010.ne2np4_conserv.POPmonthlyClimFromACES4BGC_c20260730.nc)
    GetInputFile(${surface_input_dir}/${dms_file})
    GetInputFile(${surface_input_dir}/cmip6_mam4_so2_surf_ne2np4_2010_clim_c20260730.nc)
    foreach(species bc_a4 num_a1 num_a2 num_a4 pom_a4 so4_a1 so4_a2)
      set(species_file
        cmip6_mam4_${species}_surf_ne2np4_2010_clim_c20260730.nc)
      GetInputFile(${surface_input_dir}/${species_file})
    endforeach()
  endif()
endif()
