03/2013 CGB and KJE and JER

HOMME now has a CMake build option for sweqx, swim, and preqx.
It is in the BETA test mode, alert Chris Baker or Kate Evans or Jennifer Ribbeckof problems/unclear instructions
(swim is the SW version of HOMME that uses trilinos in the implicit solve option)

1. mkdir BUILD_DIR somewhere, usually the main trunk directory
2. cp /bld/cmake-script/$APPROPRIATE_BUILD_SCRIPT into $BUILD_DIR 
3. export HOMME_ROOT="location_of trunk"
3. JAGUAR ONLY: export XTPE_LINK_TYPE='dynamic' needed right now TODO: put into script build process
3. Linux box ONLY: export Z_DIR='/usr/lib64' needed right now TODO: put into script build process
4. Modify script as appropriate (DEBUG or not etc)
 a. PLEV=# vertical levels
 b. NUM_POINTS=np
 c. -D CMAKE_INSTALL_PREFIX=where /bin/$EXE will sit 
5. ./$APPROPRIATE_BUILD_SCRIPT
5. make -j4
6. make install (where you told it to go)

most build failures are do to residual build info. In the build directory:
rm -rf CMakeCache.txt CMakeFiles src utils cmake_install.cmake Makefile bin install_manifest.txt *.h
to get a fresh build starting point.

