# Add all subdirectories to cmake build

file(GLOB FILES "*")                              # get list of all files
foreach(file ${FILES})                            # for each file
  if(IS_DIRECTORY ${file})                        # if the file is a dir
    add_subdirectory(${file})                     # add subdir to cmake
  endif()
endforeach(file)