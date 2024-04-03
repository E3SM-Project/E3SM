

source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
tempest_root=~/.conda/envs/jinbo

cd /lcrc/group/e3sm/ac.xie7/gwd/E3SMv2_maint2.0_20230513/components/eam/tools/topo_tool/cube_to_target
# Generate the element mesh.
${tempest_root}/bin/GenerateCSMesh --alt --res 4 --file topo2/ne4/ne4.g
#exit
# Generate the target physgrid mesh.
${tempest_root}/bin/GenerateVolumetricMesh --in topo2/ne4/ne4.g --out topo2/ne4/ne4pg2.g --np 2 --uniform
# Generate a high-res target physgrid mesh for cube_to_target.
${tempest_root}/bin/GenerateVolumetricMesh --in topo2/ne4/ne4.g --out topo2/ne4/ne4pg4.g --np 4 --uniform
# Generate SCRIP files for cube_to_target.
${tempest_root}/bin/ConvertMeshToSCRIP --in topo2/ne4/ne4pg4.g --out topo2/ne4/ne4pg4_scrip.nc
${tempest_root}/bin/ConvertMeshToSCRIP --in topo2/ne4/ne4pg2.g --out topo2/ne4/ne4pg2_scrip.nc
