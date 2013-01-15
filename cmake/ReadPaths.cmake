

FILE(STRINGS ${PROJECT_SOURCE_DIR}/systemPaths/${Homme_Hostname}.${CMAKE_Fortran_COMPILER_ID}
     Homme_Raw_Paths
     LIMIT_COUNT 4)

# Split each 
FOREACH(LINE ${Homme_Raw_Paths})

  # Turn the line into a list
  STRING(REPLACE "=" ";" LINE_LIST ${LINE})

  # Get the 0th and 1th items of the list
  LIST(GET LINE_LIST 0 Homme_Dep)
  LIST(GET LINE_LIST 1 Homme_Dep_Path)

  # Remove whitespace
  STRING(STRIP ${Homme_Dep} Homme_Dep)
  STRING(STRIP ${Homme_Dep_Path} Homme_Dep_Path)

  # Create a new variable called Homme_PackageName_Path
  SET(Homme_${Homme_Dep}_DIR "${Homme_Dep_Path}")

endforeach()

MESSAGE(STATUS "Homme_NETCDF_DIR = ${Homme_NETCDF_DIR}")
MESSAGE(STATUS "Homme_PNETCDF_DIR = ${Homme_PNETCDF_DIR}")
