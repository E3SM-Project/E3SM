#
# Makefile to build clm testing documentation
#

# Get list of tests_ files
SOURCES    = $(wildcard tests_*)

all: test_table.html

test_table.html: $(SOURCES)
	gen_test_table.sh

