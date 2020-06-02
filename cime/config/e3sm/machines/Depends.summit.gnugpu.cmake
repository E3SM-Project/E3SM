
set(FILES_DECLARE_CUDA
cam/src/physics/crm/samxx/abcoefs.cpp
cam/src/physics/crm/samxx/accelerate_crm.cpp
cam/src/physics/crm/samxx/adams.cpp
cam/src/physics/crm/samxx/advect2_mom_xy.cpp
cam/src/physics/crm/samxx/advect2_mom_z.cpp
cam/src/physics/crm/samxx/advect_all_scalars.cpp
cam/src/physics/crm/samxx/advect_mom.cpp
cam/src/physics/crm/samxx/advect_scalar.cpp
cam/src/physics/crm/samxx/advect_scalar2D.cpp
cam/src/physics/crm/samxx/advect_scalar3D.cpp
cam/src/physics/crm/samxx/bound_duvdt.cpp
cam/src/physics/crm/samxx/bound_exchange.cpp
cam/src/physics/crm/samxx/boundaries.cpp
cam/src/physics/crm/samxx/buoyancy.cpp
cam/src/physics/crm/samxx/cloud.cpp
cam/src/physics/crm/samxx/coriolis.cpp
cam/src/physics/crm/samxx/crm.cpp
cam/src/physics/crm/samxx/crmsurface.cpp
cam/src/physics/crm/samxx/damping.cpp
cam/src/physics/crm/samxx/diagnose.cpp
cam/src/physics/crm/samxx/diffuse_mom.cpp
cam/src/physics/crm/samxx/diffuse_mom2D.cpp
cam/src/physics/crm/samxx/diffuse_mom3D.cpp
cam/src/physics/crm/samxx/diffuse_scalar.cpp
cam/src/physics/crm/samxx/diffuse_scalar2D.cpp
cam/src/physics/crm/samxx/diffuse_scalar3D.cpp
cam/src/physics/crm/samxx/forcing.cpp
cam/src/physics/crm/samxx/ice_fall.cpp
cam/src/physics/crm/samxx/kurant.cpp
cam/src/physics/crm/samxx/microphysics.cpp
cam/src/physics/crm/samxx/periodic.cpp
cam/src/physics/crm/samxx/post_icycle.cpp
cam/src/physics/crm/samxx/post_timeloop.cpp
cam/src/physics/crm/samxx/pre_timeloop.cpp
cam/src/physics/crm/samxx/precip_init.cpp
cam/src/physics/crm/samxx/precip_proc.cpp
cam/src/physics/crm/samxx/press_grad.cpp
cam/src/physics/crm/samxx/press_rhs.cpp
cam/src/physics/crm/samxx/pressure.cpp
cam/src/physics/crm/samxx/setparm.cpp
cam/src/physics/crm/samxx/setperturb.cpp
cam/src/physics/crm/samxx/sgs.cpp
cam/src/physics/crm/samxx/shear_prod2D.cpp
cam/src/physics/crm/samxx/shear_prod3D.cpp
cam/src/physics/crm/samxx/task_init.cpp
cam/src/physics/crm/samxx/timeloop.cpp
cam/src/physics/crm/samxx/tke_full.cpp
cam/src/physics/crm/samxx/uvw.cpp
cam/src/physics/crm/samxx/vars.cpp
cam/src/physics/crm/samxx/zero.cpp
../externals/YAKL/YAKL.cpp )


foreach(ITEM IN LISTS FILES_DECLARE_CUDA)
  e3sm_add_flags("${ITEM}" "--expt-extended-lambda --expt-relaxed-constexpr -arch sm_70 -O3 --use_fast_math -std=c++14 -D__USE_CUDA__ -I${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/cub")
  e3sm_declare_cuda("${ITEM}")
endforeach()

