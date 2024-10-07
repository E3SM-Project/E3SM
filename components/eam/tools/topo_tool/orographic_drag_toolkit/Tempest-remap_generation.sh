

source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
tempest_root=~/.conda/envs/jinbo
# Generate the element mesh.
${tempest_root}/bin/GenerateCSMesh --alt --res 30 --file topo2/ne30.g
# Generate the target physgrid mesh.
${tempest_root}/bin/GenerateVolumetricMesh --in topo2/ne30.g --out topo2/ne30pg2.g --np 2 --uniform
# Generate a high-res target physgrid mesh for cube_to_target.
${tempest_root}/bin/GenerateVolumetricMesh --in topo2/ne30.g --out topo2/ne30pg4.g --np 4 --uniform
# Generate SCRIP files for cube_to_target.
${tempest_root}/bin/ConvertMeshToSCRIP --in topo2/ne30pg4.g --out topo2/ne30pg4_scrip.nc
${tempest_root}/bin/ConvertMeshToSCRIP --in topo2/ne30pg2.g --out topo2/ne30pg2_scrip.nc
