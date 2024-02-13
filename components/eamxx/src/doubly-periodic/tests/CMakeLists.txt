INCLUDE (ScreamUtils)

set(DP_TESTS_SRCS
  dp_unit_tests.cpp
    dp_advance_iop_forcing_tests.cpp
    dp_advance_iop_nudging_tests.cpp
    dp_advance_iop_subsidence_tests.cpp
    dp_iop_setinitial_tests.cpp
    dp_iop_broadcast_tests.cpp
    dp_apply_iop_forcing_tests.cpp
    dp_iop_domain_relaxation_tests.cpp
    dp_crm_resolved_turb_tests.cpp
    dp_iop_default_opts_tests.cpp
    dp_iop_setopts_tests.cpp
    dp_setiopupdate_init_tests.cpp
    dp_setiopupdate_tests.cpp
    dp_readiopdata_tests.cpp
    dp_iop_intht_tests.cpp
    ) # DP_TESTS_SRCS

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  CreateUnitTest(dp_tests "${DP_TESTS_SRCS}"
    LIBS dp
    THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
  )
endif()
