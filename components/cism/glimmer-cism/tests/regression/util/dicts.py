"""A set of dictionaries for the regression tests.

These dictionaries will determine what tests are run and what options to pass
to each test. See descriptions of each dictionary.
"""

# SPECIFIC TESTS DICTIONARIES
# ===========================
# These dictionaries describe the NON DEFAULT options to pass to the test run scripts found within
# the CISM tests directory. These modified test cases are only called if the performance option is
# specified.
# 
# Dictionaries for use in the MAIN TEST DICTIONARY (at bottom). They should
# consist of key-value pairs like:
#   key: anything descriptive -- used only internally
#   value: a string containing the options to pass to the script. 
# See the dictionaries below for an example.
# NOTE: these options should not include the -e/--executable, --hpc, -o/--output-dir, -m/--modifier, 
#       or -s/--setup-only options because the build_and_test script will include those automatically.

# The higher-order/dome test
# --------------------------
# for tests with -n N for N < 16
dome_perf_small = { 
        's1': '--scale 1 -n 1',
        's2': '--scale 2 -n 1',
        'p1': '-n 2',
        'p2': '-n 4',
        'b1': '--scale 1 -n 4',
        'p3': '-n 8',
        'b2': '--scale 2 -n 16',
        }

# for tests with -n N for N > 16
dome_perf_large = {
        's3': '--scale 2 -n 256',
        's4': '--scale 3 -n 256',
        'b3': '--scale 3 -n 64',
        'b4': '--scale 4 -n 256',
        }

# The higher-order/shelf tests
# ----------------------------
# NOTE: empty dict because no performance testing for confined shelf. Leaving here for possible
#       future expansion. 
shelfConfined_perf_small = {}

# NOTE: empty dict because no performance testing for circular shelf. Leaving here for possible
#       future expansion. 
shelfCircular_perf_small = {}

# NOTE: empty dict because no performance testing for ISMIP-HOM. Leaving here for possible
#       future expansion. 
ismip_perf_small = {}

# NOTE: empty dict because no performance testing for stream. Leaving here for possible
#       future expansion. 
stream_perf_small = {}

# MAIN TEST DICTIONARY
# ====================
# This is the main dictionary that describes what tests to run.
# Each dictionary item should consist of key-value pairs like:
#   key: path to test from $CISM/tests
#       example: 'higher-order/dome'
#       NOTE: key can be a space separated list with the first entry the path to the test directory
#             and the rest of the list used to define uniqueness. This is useful for tests that have
#             multiple run scripts like shelf.
#           example: 'higher-order/shelf Confined'
#       
#   value: tuple of (run_script, options_dict) where run_script is the test run
#          script that can be found within the directory specified by the key
#          and options_dict is a dictionary containing the options to pass to the
#          run_script (as described above)
#        example: ('runDome.py', dome_perf_small)
test_dict = {
        'higher-order/dome': ('runDome.py', dome_perf_small),
        'higher-order/shelf Confined': ('runShelfConfined.py', shelfConfined_perf_small),
        'higher-order/shelf Circular': ('runShelfCircular.py', shelfConfined_perf_small),
        'higher-order/ismip-hom': ('runISMIP_HOM.py -r a c f', ismip_perf_small),
        'higher-order/stream': ('runStream.py', stream_perf_small),
        }

perf_dict = {
        'higher-order/dome': ('runDome.py', dome_perf_large),
        }

# HPC PLATFORM DICTIONARIES
# =========================
# There should be a dictionary for each supported HPC platform which specifies 
# the default batch scheduler options to use. 
hopper_dict = {
        'PBS_A': 'm1795',
        'PBS_q': 'regular',
        'PBS_N': 'reg_test_all',
        'PBS_RES': 'mppwidth',
        'RES_NUM': '24',
        'PBS_walltime': '01:00:00',
        }


titan_dict = {
        'PBS_A': 'cli106',
        'PBS_q': 'batch',
        'PBS_N': 'reg_test_all',
        'PBS_RES': 'nodes',
        'RES_NUM': '1',
        'PBS_walltime': '01:00:00',
        }


# MAIN HPC DICTIONARY
# ===================
# Collection of all the HPC platform dictionaries. 
hpc_dict = {
        'titan': titan_dict,
        'hopper': hopper_dict,
        }
