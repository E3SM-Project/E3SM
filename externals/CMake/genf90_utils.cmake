# Utility for invoking genf90 on a template file.
#
# If ENABLE_GENF90 is set to a true value, the functions here will behave
# as described below. In this case, the variable GENF90 must be defined and
# contain the genf90.pl command.
#
# If ENABLE_GENF90 is not true, no source code generation or other side
# effects will occur, but output variables will be set as if the generation
# had occurred.
#
#==========================================================================
#
# process_genf90_source_list
#
# Arguments:
#    genf90_file_list - A list of template files to process.
#    output_directory - Directory where generated sources will be placed.
#    fortran_list_name - The name of a list used as output.
#
# Produces generated sources for each of the input templates. Then
# this function *appends* the location of each generated file to the output
# list.
#
# As a side effect, this function will add a target for each generated
# file. For a generated file named "foo.F90", the target will be named
# "generate_foo".
#
# Limitations:
#    This function adds targets to work around a deficiency in CMake (see
#    "declare_generated_dependencies" in Sourcelist_utils). Unfortunately,
#    this means that you cannot use this function to generate two files
#    with the same name in a single project.
#
#==========================================================================

#==========================================================================
# Copyright (c) 2013-2014, University Corporation for Atmospheric Research
#
# This software is distributed under a two-clause BSD license, with no
# warranties, express or implied. See the accompanying LICENSE file for
# details.
#==========================================================================

if(ENABLE_GENF90)

  # Notify CMake that a Fortran file can be generated from a genf90
  # template.
  function(preprocess_genf90_template genf90_file fortran_file)

    add_custom_command(OUTPUT ${fortran_file}
      COMMAND ${GENF90} ${genf90_file} >${fortran_file}
      MAIN_DEPENDENCY ${genf90_file})

    get_filename_component(stripped_name ${fortran_file} NAME_WE)

    add_custom_target(generate_${stripped_name} DEPENDS ${fortran_file})

  endfunction(preprocess_genf90_template)

else()

  # Stub if genf90 is off.
  function(preprocess_genf90_template)
  endfunction()

endif()

# Auto-generate source names.
function(process_genf90_source_list genf90_file_list output_directory
    fortran_list_name)

  foreach(genf90_file IN LISTS genf90_file_list)

    # If a file is a relative path, expand it (relative to current source
    # directory.
    get_filename_component(genf90_file "${genf90_file}" ABSOLUTE)

    # Get extensionless base name from input.
    get_filename_component(genf90_file_stripped "${genf90_file}" NAME_WE)

    # Add generated file to the test list.
    set(fortran_file ${output_directory}/${genf90_file_stripped}.F90)
    preprocess_genf90_template(${genf90_file} ${fortran_file})
    list(APPEND ${fortran_list_name} ${fortran_file})
  endforeach()

  # Export ${fortran_list_name} to the caller.
  set(${fortran_list_name} "${${fortran_list_name}}" PARENT_SCOPE)

endfunction(process_genf90_source_list)
