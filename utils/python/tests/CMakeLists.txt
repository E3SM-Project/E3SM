cmake_minimum_required(VERSION 2.8)

    MACRO(STRING_UNQUOTE var str)

	# ';' and '\' are tricky, need to be encoded.
	# '\' => '#B'
	# '#' => '#H'
	# ';' => '#S'
	STRING(REGEX REPLACE "#" "#H" _ret "${str}")
	STRING(REGEX REPLACE "\\\\" "#B" _ret "${_ret}")
	STRING(REGEX REPLACE ";" "#S" _ret "${_ret}")

	IF(_ret MATCHES "^[ \t\r\n]+")
	    STRING(REGEX REPLACE "^[ \t\r\n]+" "" _ret "${_ret}")
	ENDIF(_ret MATCHES "^[ \t\r\n]+")
	IF(_ret MATCHES "^\"")
	    # Double quote
	    STRING(REGEX REPLACE "\"\(.*\)\"[ \t\r\n]*$" "\\1" _ret "${_ret}")
	ELSEIF(_ret MATCHES "^'")
	    # Single quote
	    STRING(REGEX REPLACE "'\(.*\)'[ \t\r\n]*$" "\\1" _ret "${_ret}")
	ELSE(_ret MATCHES "^\"")
	    SET(_ret "")
	ENDIF(_ret MATCHES "^\"")

	# Unencoding
	STRING(REGEX REPLACE "#B" "\\\\" _ret "${_ret}")
	STRING(REGEX REPLACE "#H" "#" _ret "${_ret}")
	STRING(REGEX REPLACE "#S" "\\\\;" ${var} "${_ret}")
    ENDMACRO(STRING_UNQUOTE var str)

include(CTest)

execute_process(COMMAND "python" "list_tests" OUTPUT_VARIABLE STR_TESTS
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_STRIP_TRAILING_WHITESPACE)
separate_arguments(TEST_LIST UNIX_COMMAND ${STR_TESTS} )

foreach(ATEST ${TEST_LIST})
    # this assignment prevents quotes being added to testname in add_test
    set(fulltest "${ATEST}")
    add_test(NAME ${ATEST} COMMAND ./scripts_regression_tests.py ${fulltest}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach(ATEST)

