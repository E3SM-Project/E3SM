Use of Zoltan2 within HOMME.

- Compile/install Trilinos with Zoltan2 enabled. 
  Make sure CMAKE config files are installed.
  Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=ON

- When configuring HOMME, provide HOMME_USE_TRILINOS=TRUE, 
and Trilinos_DIR=<path_to_trilinos_installation>

- Zoltan2 can perform partitioning and/or task mapping.
  Partitioning: is the assigment of the tasks to the processors.
	     	This is important for load balancing and minimization 
		of interprocessor communication.
  Task Mapping: is architecture aware task placement. The algorithm
		is aware of the allocated nodes within the network,
		and it will place the tasks based on this information.
		This aims to reduce the network link contention and 
		the distance messages travel. This reduces the interprocessor 
		communication cost, by utilizing network links in better ways.

  It can be considered as a 2 step method. 
	1) partitioning is done either using space filling curves that is 
	   implemented in HOMME, or any Zoltan2 methods.
  
	2) an optional architecture aware task mapping can be done on the 
	   result of partitioning regardless of the method used. 

- Parameters controlling zoltan usage.

  partmethod: value-4 was for previous space filling curves (SFC).
	      values from 5 to 22 are the Zoltan2 partitioning methods. 
	      It provides different algorithms for testing purposes but,
	      currently best working ones are 5,6,7,8 (geometric methods.)

	Suggested Parameter: 5 for quality. 5 is okay probably upto 1M tasks (columns).
	   	  	     6 for scalability. 
			     4 for when Zoltan2 is not enabled. Zoltan methods will throw a run time error if it is not enabled..

				   
  z2_map_method: Task Mapping method that will be used by zoltan. 
		 1 - No task mapping, in this case no architecture aware task placement is done. 
			Only partitioning is performed. 
		 2 - Task mapping is performed.
		 3 - Optimized task mapping is performed. 
	Suggested Parameter: 
	    	 3 - when Zoltan is enabled.
		 1 - when Zoltan is not enabled. [2-3] will throw run time error if zoltan is not enabled.
	
  coord_transform_method: Coordinate transformation method. Zoltan will use
                1 - Sphere coordinates
                2 - Projected cube coordinates
                3 - 2D face coordinates that is obtained by opening a box corresponding to cube coordinates.
        Suggested Parameter:
                - Does not matter, value is discarded, if Zoltan is not enabled.
                - 3 if zoltan methods are used.
		- 2 if SFC is used for partitioning, and Zoltan2 is used for mapping.

  OVERAL SUGGESTED PARAMETERS: partmethod=5 coord_transform_method=3 z2_map_method=3 WITH ZOLTAN
			       partmethod=4 z2_map_method=1 without zoltan

- TRILINOS INSTALLATION FOR TASK MAPPING.
  1) If you want to use task mapping for blue-gene/Q machines (mira, vulcan):
	1) Download and install https://github.com/bhatele/topomgr
	2) Download an Trilinos. When configuring trilinos
	    *) set TPL_ENABLE_TopoManager:BOOL=ON 
	   If cmake should find the topomanager folders. If it does not you might need to provide 
	    *)  TPL_TopoManager_INCLUDE_DIRS:PATH="<path_to_topomanager_include_path>" \
		TPL_TopoManager_LIBRARIES:FILEPATH="<path_to_topomanagerlibfolder/libtopomgr.a>
	3) Then configure HOMME with trilinos as before.
  2) For Titan, Edison or Bluewaters, trilinos use rca library to capture the underlying network.
  This is upcoming, the setup of this will require a configure of trilinos using rca as above. 
 
  3) If you do not provide underlying network library to trilinos, 
     it will do the mapping based on mpi ranks, assuming that consecutive ranks are closer to each other.

If you have problems, contact
Mehmet Deveci
mndevec@sandia.gov


