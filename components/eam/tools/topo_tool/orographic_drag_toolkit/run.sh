tool_root=/lcrc/group/e3sm/ac.xie7/gwd/E3SMv2/code/20240712/components/eam/tools/topo_tool/gwd/
	${tool_root}/cube_to_target \
	  --target-grid  ${tool_root}/homme_production/topo2/ne30pg2_scrip.nc \
	  --input-topography /lcrc/group/e3sm/data/inputdata/atm/cam/hrtopo/USGS-topo-cube3000.nc \
	  --output-topography  ${tool_root}/homme_production/topo2/ne30pg4_c2t_topo.nc

