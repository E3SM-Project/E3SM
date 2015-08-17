include (CMakeParseArguments)

#==============================================================================
#
#  FUNCTIONS TO HELP WITH Check* MODULES
#
#==============================================================================

#______________________________________________________________________________
# - Basic function to check a property of a package using a try_compile step
#
# SYNTAX:  check_property (<return_variable>
#                          NAME <filename>
#                          HINTS <path> <path> ...
#                          DEFINITIONS <definition1> <definition> ...
#                          COMMENT <string_comment>)
#                         
function (check_property VARIABLE)

    # Parse the input arguments
    set (options)
    set (oneValueArgs COMMENT NAME)
    set (multiValueArgs HINTS DEFINITIONS)
    cmake_parse_arguments (${VARIABLE} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})    
    
    # If the return variable is defined, already, don't continue
    if (NOT DEFINED ${VARIABLE})
        
        message (STATUS "Checking ${${VARIABLE}_COMMENT}")
        find_file (${VARIABLE}_TRY_FILE
                   NAMES ${${VARIABLE}_NAME}
                   HINTS ${${VARIABLE}_HINTS})
        if (${VARIABLE}_TRY_FILE)
            try_compile (COMPILE_RESULT
                         ${CMAKE_CURRENT_BINARY_DIR}/try${VARIABLE}
                         SOURCES ${${VARIABLE}_TRY_FILE}
                         COMPILE_DEFINITIONS ${${VARIABLE}_DEFINITIONS}
                         OUTPUT_VARIABLE TryOUT)
            if (COMPILE_RESULT)
                message (STATUS "Checking ${${VARIABLE}_COMMENT} - yes")
            else ()
                message (STATUS "Checking ${${VARIABLE}_COMMENT} - no")
            endif ()
            
            set (${VARIABLE} ${COMPILE_RESULT}
                 CACHE BOOL "${${VARIABLE}_COMMENT}")
                 
        else ()
            message (STATUS "Checking ${${VARIABLE}_COMMENT} - failed")
        endif ()
     
        unset (${VARIABLE}_TRY_FILE CACHE)
    endif ()
            
endfunction ()

#______________________________________________________________________________
# - Basic function to check the version of a package using a try_run step
#
# SYNTAX:  check_version (<pkg>
#                         NAME <try_version_file>
#                         HINTS <path> <path> ...
#                         DEFINITIONS <definition1> <definition> ...)
#  
function (check_version PKG)

    # Parse the input arguments
    set (options)
    set (oneValueArgs NAME)
    set (multiValueArgs HINTS DEFINITIONS)
    cmake_parse_arguments (${PKG} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})    

    # If the return variable is defined, already, don't continue
    if (NOT DEFINED ${PKG}_VERSION)
        
        message (STATUS "Checking ${PKG} version")
        find_file (${PKG}_VERSION_TRY_FILE
                   NAMES ${${PKG}_NAME}
                   HINTS ${${PKG}_HINTS})
        if (${PKG}_VERSION_TRY_FILE)
            try_run (RUN_RESULT COMPILE_RESULT
                     ${CMAKE_CURRENT_BINARY_DIR}/try${PKG}_VERSION
                     ${${PKG}_VERSION_TRY_FILE}
                     COMPILE_DEFINITIONS ${${PKG}_DEFINITIONS}
                     COMPILE_OUTPUT_VARIABLE TryCompileOUT
                     RUN_OUTPUT_VARIABLE TryRunOUT)
            if (COMPILE_RESULT AND NOT RUN_RESULT STREQUAL FAILED_TO_RUN)
                message (STATUS "Checking ${PKG} version - ${TryRunOUT}")
                set (${PKG}_VERSION ${TryRunOUT}
                     CACHE STRING "${PKG} version string")
                if (${PKG}_VERSION VERSION_LESS ${PKG}_FIND_VERSION})
                    message (FATAL_ERROR "${PKG} version insufficient")
                endif ()
            else ()
                message (STATUS "Checking ${PKG} version - failed")
            endif ()
            
            unset (COMPILE_RESULT CACHE)
            unset (RUN_RESULT CACHE)
            
        else ()
            message (STATUS "Checking ${PKG} version - failed")
        endif ()
        
        unset (${PKG}_VERSION_TRY_FILE CACHE)
     
    endif ()

endfunction ()