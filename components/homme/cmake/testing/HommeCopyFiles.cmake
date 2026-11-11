# Copy some files from a dir to another. Usage:
#   cmake -DSRC_DIR=<dir> -DDST_DIR=<dir> -DFILES=<f1,f2,...> -P HommeCopyFiles.cmake
# or include this file and call homme_copy_files(<src-dir> <dst-dir> <f1,f2,...>).
# Used to copy test outputs into the baselines dir. It is an error if a file is missing.

function (homme_copy_files src_dir dst_dir files)
  string (REPLACE "," ";" files "${files}")
  file (MAKE_DIRECTORY ${dst_dir})
  foreach (f IN LISTS files)
    if (NOT EXISTS ${src_dir}/${f})
      message (FATAL_ERROR "Cannot copy ${src_dir}/${f}: file not found")
    endif()
    file (COPY ${src_dir}/${f} DESTINATION ${dst_dir})
    message ("Copied ${src_dir}/${f} to ${dst_dir}")
  endforeach()
endfunction()

if (CMAKE_SCRIPT_MODE_FILE STREQUAL CMAKE_CURRENT_LIST_FILE)
  homme_copy_files ("${SRC_DIR}" "${DST_DIR}" "${FILES}")
endif()
