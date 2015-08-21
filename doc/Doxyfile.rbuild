# This file contains the customizations required to override
# settings in Doxyfile to enable building documentation from
# the build directory

# The INPUT tag is used to specify the the location of the source
INPUT               = @CMAKE_SOURCE_DIR@

# The EXAMPLE_PATH tag can be used to specify one or more files or directories
# that contain example code fragments that are included (see the \include
# command).

EXAMPLE_PATH           = @CMAKE_SOURCE_DIR@/doc/example \
                         @CMAKE_SOURCE_DIR@/testpio

