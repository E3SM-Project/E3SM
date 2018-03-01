set (PFUNIT_PREPROCESSOR python $ENV{PFUNIT}/bin/pFUnitParser.py)
set (PFUNIT_SUFFIX "\\.pfunit")
# - Pass a list of files through the pFUnit macro processor
#
# ADD_PFUNIT_SOURCES( OUTVAR source1 ... sourceN )
#
#    OUTVAR    A list containing all the output file names, suitable
#                    to be passed to add_executable or add_library.
#
# If the source files have a .m4 suffix it is stripped from the output
# file name. The output files are placed in the same relative location
# to CMAKE_CURRENT_BINARY_DIR as they are to CMAKE_CURRENT_SOURCE_DIR.
#
# Example:
#    add_pfunit_sources( SRCS src/test1.cxx.pfunit src/test2.cxx.pfunit )
#    add_executable( test ${SRCS} )
function( ADD_PFUNIT_SOURCES OUTVAR )
     set( outfiles )
     foreach( f ${ARGN} )
         # first we might need to make the input file absolute
         get_filename_component( f "${f}" ABSOLUTE )
         # get the relative path of the file to the current source dir
         file( RELATIVE_PATH baseFile "${CMAKE_CURRENT_SOURCE_DIR}" "${f}" )
         # strip the .pfunit off the end if present
         string( REGEX REPLACE "${PFUNIT_SUFFIX}" ".F90" outFile "${CMAKE_CURRENT_BINARY_DIR}/${baseFile}" )
         # append the output file to the list of outputs
         list( APPEND outfiles "${outFile}" )
         # create the output directory if it doesn't exist
         get_filename_component( dir "${outFile}" PATH )
         if( NOT IS_DIRECTORY "${dir}" )
             file( MAKE_DIRECTORY "${dir}" )
         endif( NOT IS_DIRECTORY "${dir}" )
         # now add the custom command to generate the output file
         add_custom_command( OUTPUT "${outFile}"
             COMMAND ${PFUNIT_PREPROCESSOR} "${f}" "${outFile}"
             DEPENDS "${f}"
             )
     endforeach( f )
     # set the output list in the calling scope
     set( ${OUTVAR} ${outfiles} PARENT_SCOPE )
endfunction( ADD_PFUNIT_SOURCES )


