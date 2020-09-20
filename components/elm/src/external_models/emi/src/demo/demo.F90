program demo

  use ExternalModelInterfaceMod
  use clm_varctl                , only : iulog
  use decompMod                 , only : bounds_type, get_proc_bounds, get_proc_clumps, get_clump_bounds
  use elm_instMod               , only : elm_inst_biogeophys
  use clm_varpar                , only : clm_varpar_init
  use clm_varcon                , only : clm_varcon_init
  use clm_varpar                , only : nlevdecomp_full, ndecomp_pools
  use spmdMod                   , only : spmd_init
  use ExternalModelConstants    , only : EM_ID_STUB, EM_STUB_SOIL_HYDRO_STAGE, EM_STUB_SOIL_THERMAL_STAGE
  use elm_instMod               , only : soilstate_vars, waterstate_vars, waterflux_vars
  use elm_instMod               , only : energyflux_vars, temperature_vars, carbonstate_vars
  use shr_kind_mod              , only : r8 => shr_kind_r8, SHR_KIND_CL

  implicit none

  type(bounds_type) :: bounds_proc, bounds_clump
  integer           :: clump_rank, nclumps
  integer           :: c, num_hydrologyc, num_nolakec_and_nourbanc
  integer, pointer  :: filter_hydrologyc(:), filter_nolakec_and_nourbanc(:)

  nclumps = 1
  write(iulog,*)''
  write(iulog,*)'This is a demo for the External Model Interface (EMI)'
  write(iulog,*)''

  call spmd_init()
  call set_namelist_variables()
  call clm_varpar_init()
  call clm_varcon_init()
  call decompInit()

  call get_proc_bounds(bounds_proc)

  call elm_inst_biogeophys(bounds_proc)
  call initialize_clm_data_structures(bounds_proc)

  call EMI_Determine_Active_EMs()

  write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(iulog,*)'1. Lets initialize the Stub EM'
  write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call EMI_Init_EM(EM_ID_STUB)

  num_hydrologyc = 1
  num_nolakec_and_nourbanc = 1
  allocate(filter_hydrologyc(num_hydrologyc))
  allocate(filter_nolakec_and_nourbanc(num_nolakec_and_nourbanc))
  do c = 1, num_hydrologyc
     filter_hydrologyc(c) = c
     filter_nolakec_and_nourbanc(c) = c
  end do

  write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(iulog,*)'2. Lets now timestep the Stub EM'
  write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  write(*,*)'num_hydrologyc : ',num_hydrologyc
  write(*,*)'nlevdecomp_full: ',nlevdecomp_full
  write(*,*)'ndecomp_pools  : ',ndecomp_pools

  !$OMP PARALLEL DO PRIVATE (clump_rank, bounds_clump)
  do clump_rank = 1, nclumps

     call get_clump_bounds(clump_rank, bounds_clump)

     write(iulog,*)'  Stub EM will solve hydrologic processes'
     call EMI_Driver(                                                 &
          em_id             = EM_ID_STUB                            , &
          em_stage          = EM_STUB_SOIL_HYDRO_STAGE              , &
          dt                = 1800._r8                              , &
          clump_rank        = bounds_clump%clump_index              , &
          num_hydrologyc    = num_hydrologyc                        , &
          filter_hydrologyc = filter_hydrologyc                     , &
          num_nolakec_and_nourbanc = num_nolakec_and_nourbanc       , &
          filter_nolakec_and_nourbanc = filter_nolakec_and_nourbanc , &
          waterstate_vars   = waterstate_vars                       , &
          waterflux_vars    = waterflux_vars                        , &
          energyflux_vars   = energyflux_vars                       , &
          soilstate_vars    = soilstate_vars                        , &
          temperature_vars  = temperature_vars)

     write(iulog,*)''
     write(iulog,*)''
     write(iulog,*)'  Stub EM will solve thermal processes'
     call EMI_Driver(                                                 &
          em_id             = EM_ID_STUB                            , &
          em_stage          = EM_STUB_SOIL_THERMAL_STAGE            , &
          dt                = 1800._r8                              , &
          clump_rank        = bounds_clump%clump_index              , &
          num_hydrologyc    = num_hydrologyc                        , &
          filter_hydrologyc = filter_hydrologyc                     , &
          num_nolakec_and_nourbanc = num_nolakec_and_nourbanc       , &
          filter_nolakec_and_nourbanc = filter_nolakec_and_nourbanc , &
          waterstate_vars   = waterstate_vars                       , &
          waterflux_vars    = waterflux_vars                        , &
          energyflux_vars   = energyflux_vars                       , &
          soilstate_vars    = soilstate_vars                        , &
          temperature_vars  = temperature_vars                      , &
          carbonstate_vars  = carbonstate_vars)
  enddo
  !$OMP END PARALLEL DO

  
end program demo

!-----------------------------------------------------------------------
subroutine set_namelist_variables()

  use clm_varctl, only : use_em_stub

  implicit none

  use_em_stub = .true.

end subroutine set_namelist_variables
!-----------------------------------------------------------------------
subroutine decompInit ()
  !
  use decompMod
  use clm_varctl, only : iulog
  use abortutils      , only : endrun
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use spmdMod , only : iam
  !
  implicit none
  integer :: ier                    ! error code
  integer :: ncells, ntopounits, nlunits, ncols, npfts, nCohorts

  clump_pproc = 1
  nclumps     = 1

  ncells      = 1
  ntopounits  = 1
  nlunits     = 1
  ncols       = 1
  npfts       = 16
  nCohorts    = 1

  allocate(procinfo%cid(clump_pproc), stat=ier)
  if (ier /= 0) then
     write(iulog,*) 'decompInit_lnd(): allocation error for procinfo%cid'
     call endrun(msg=errMsg(__FILE__, __LINE__))
  endif
  procinfo%nclumps    = clump_pproc
  procinfo%cid(:)     = 1
  procinfo%ncells     = ncells
  procinfo%ntopounits = ntopounits
  procinfo%nlunits    = nlunits
  procinfo%ncols      = ncols

  procinfo%npfts      = npfts
  procinfo%nCohorts   = nCohorts
  procinfo%begg       = 1
  procinfo%begt       = 1
  procinfo%begl       = 1
  procinfo%begc       = 1
  procinfo%begp       = 1
  procinfo%begCohort  = 1
  procinfo%endg       = ncells
  procinfo%endt       = ntopounits
  procinfo%endl       = nlunits
  procinfo%endc       = ncols
  procinfo%endp       = npfts
  procinfo%endCohort  = nCohorts

  allocate(clumps(nclumps), stat=ier)
  if (ier /= 0) then
     write(iulog,*) 'decompInit_lnd(): allocation error for clumps'
     call endrun(msg=errMsg(__FILE__, __LINE__))
  end if
  clumps(:)%owner      = iam
  clumps(:)%ncells     = ncells
  clumps(:)%ntopounits = ntopounits
  clumps(:)%nlunits    = nlunits
  clumps(:)%ncols      = ncols
  clumps(:)%npfts      = npfts
  clumps(:)%nCohorts   = nCohorts
  clumps(:)%begg       = 1
  clumps(:)%begt       = 1
  clumps(:)%begl       = 1
  clumps(:)%begc       = 1
  clumps(:)%begp       = 1
  clumps(:)%begCohort  = 1
  clumps(:)%endg       = ncells
  clumps(:)%endt       = ntopounits
  clumps(:)%endl       = nlunits
  clumps(:)%endc       = ncols
  clumps(:)%endp       = npfts
  clumps(:)%endCohort  = nCohorts

end subroutine decompInit

!-----------------------------------------------------------------------
subroutine initialize_clm_data_structures(bounds_proc)

  use decompMod     , only : bounds_type
  use elm_instMod   , only : soilstate_vars!, waterstate_vars, waterflux_vars
  use elm_instMod   , only : carbonstate_vars
  !use elm_instMod   , only : energyflux_vars
  use ColumnDataType, only : col_ws, col_wf, col_ef
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use clm_varpar    , only : nlevgrnd
  use clm_varpar    , only : nlevdecomp_full, ndecomp_pools
  use shr_const_mod , only : SHR_CONST_PI
  !
  implicit none
  !
  type(bounds_type) :: bounds_proc
  !
  integer :: begc, endc, ncol
  integer :: c,j,k
  real(r8) :: counter

  begc = bounds_proc%begc; endc = bounds_proc%endc;
  ncol = endc-begc+1;
  counter = 0.d0

  do c = begc, endc
     do j = 1, nlevgrnd
        soilstate_vars%hksat_col(c,j)   = 10._r8*(c**2._r8) + j
        soilstate_vars%bsw_col(c,j)     = 0.5_r8*(c**2._r8) + j
        col_ws%h2osoi_liq(c,j) = (0.3_r8*(c**2._r8) + j)/100._r8
        col_ws%h2osoi_ice(c,j) = 1._r8 - (0.3_r8*(c**2._r8) + j)/100._r8
        col_wf%mflx_et(c,j) = 20._r8*(c**2._r8) + j
     enddo
     col_ef%eflx_hs_soil(c) = 30._r8*(c**2._r8)
     do j = 1, nlevdecomp_full
        do k = 1, ndecomp_pools
           carbonstate_vars%decomp_cpools_vr_col(c,j,k) = counter;
           counter = counter + 1.d0
        end do
     end do
     
  enddo

end subroutine initialize_clm_data_structures
