9/2012 CGB and KJE 

HOMME now has a CMake build option for both sweqx and swim.
It is in the BETA test mode
(swim is the SW version of HOMME that uses trilinos in the implicit solve option)

To use this:
* right now, either build PIO and timing separately, then comment the build lines out of the build script OR
* grab the changes needed to build those external libraries using cmake from the swim infrastructure developers

1. mkdir BUILD_DIR somewhere
2. cp /bld/cmake-script/$APPROPRIATE_BUILD_SCRIPT into $BUILD_DIR 
(remember, right now comment out the PIO timing cmake build lines)
3. ./$APPROPRIATE_BUILD_SCRIPT
4. Modify as appropriate
5. make -j4
6. make install (if you want the executable installed in /opt/homme/trunk-$DATE)

