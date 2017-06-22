cmake_minimum_required(VERSION 2.8)

#
# Top level CMakeLists.txt file for ACME
# All logic is currently in included file in CIME
# so that the same build system files can be developed
# and tested in CIME for data/stub models, and 
# then tested here subsequently. (Subject to change.)
#
SET(ACME_BUILD ON)
SET(CIME_REL_DIR cime)
SET(ACME_COMPONENTS_DIR components)

INCLUDE(cime/scripts/cmake/includedCMakeLists.txt)

