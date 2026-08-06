!================================================================================
! manages mapping of CLM generated wild fire emissions to chemical constituents
!================================================================================
module fire_emissions

  use shr_kind_mod,      only : r8 => shr_kind_r8, shr_kind_cl
  use shr_fire_emis_mod, only : shr_fire_emis_mechcomps, shr_fire_emis_mechcomps_n, shr_fire_emis_elevated
  use shr_const_mod,     only : pi => SHR_CONST_PI
  use shr_const_mod,     only : avogad => SHR_CONST_AVOGAD ! Avogadro's number ~ molecules/kmole
  use cam_abortutils,    only : endrun
  use cam_cpl_indices,   only : index_x2a_Fall_flxfire
  use cam_history,       only : addfld, horiz_only, outfld, fieldname_len, add_default
  use cam_logfile,       only : iulog
  use ppgrid,            only : pver, pverp
  use constituents,      only : cnst_get_ind
  use rad_constituents,  only : rad_cnst_get_aer_props, rad_cnst_num_name
  use mo_chem_utls,      only : get_spc_ndx, get_extfrc_ndx
  use chem_mods,         only : adv_mass ! g/mole
  use infnan,            only : nan, assignment(=)
!LXu@05/20
  use spmd_utils,        only : masterproc

  implicit none
  private
  save

  public :: fire_emissions_init
  public :: fire_emissions_srf
  public :: fire_emissions_vrt

  ! for surface emissions
  integer, allocatable :: fire_emis_indices_map(:) 
  integer, allocatable :: fire_emis_num_map(:)

  ! for vertically distributed forcings
  integer,  allocatable :: frc_spc_map(:)
  integer,  allocatable :: chm_spc_map(:)
  integer,  allocatable :: frc_num_map(:)
  real(r8), allocatable :: spc_mass_factor(:)
  real(r8), allocatable :: num_mass_factor(:)
  character(len=fieldname_len), allocatable :: fire_frc_name(:)
  character(len=fieldname_len), allocatable :: fire_numfrc_name(:)
  character(len=fieldname_len), allocatable :: fire_sflx_name(:)
  character(len=fieldname_len), allocatable :: fire_vflx_name(:)

!================================================================================
contains
!================================================================================

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine fire_emissions_init()


    ! local vars
    integer :: n, ii 
    
    integer :: frc_ndx, spc_ndx, ndx
    integer :: mode, spec
    character(len=16) :: name
    character(len=32) :: spc_name
    character(len=32) :: num_name

    real(r8), parameter :: demis_acc = 0.134e-6_r8 ! meters 
    ! volume-mean emissions diameter of primary BC/OM aerosols, see :
    ! Liu et al, Toward a minimal representation of aerosols in climate models: 
    ! Description and evaluation in the Community Atmosphere Model CAM5. 
    ! Geosci. Model Dev., 5, 709–739, doi:10.5194/gmd-5-709-2012
    ! and Table S1 in Supplement: http://www.geosci-model-dev.net/5/709/2012/gmd-5-709-2012-supplement.pdf

    real(r8), parameter :: x_numfact = 1.e-6_r8 * avogad * 6.0_r8 / (pi*(demis_acc**3))   ! 1.e-6 converts m-3 to cm-3. 
    real(r8) :: specdens  ! kg/m3
    logical :: found

    if (shr_fire_emis_mechcomps_n<1) return

    if (shr_fire_emis_elevated) then ! initialize elevated forcings

       allocate( frc_spc_map(shr_fire_emis_mechcomps_n) )
       allocate( chm_spc_map(shr_fire_emis_mechcomps_n) )
       allocate( frc_num_map(shr_fire_emis_mechcomps_n) )
       allocate( spc_mass_factor(shr_fire_emis_mechcomps_n) )
       allocate( num_mass_factor(shr_fire_emis_mechcomps_n) )
       allocate( fire_frc_name(shr_fire_emis_mechcomps_n) )
       allocate( fire_numfrc_name(shr_fire_emis_mechcomps_n) )
       allocate( fire_sflx_name(shr_fire_emis_mechcomps_n) )
       allocate( fire_vflx_name(shr_fire_emis_mechcomps_n) )

       frc_spc_map(:) = -1
       frc_num_map(:) = -1
       spc_mass_factor(:) = nan
       num_mass_factor(:) = nan

!       call addfld ('Fire_ZTOP', horiz_only, 'A', 'm', 'top of vertical fire emissions' )
!       call add_default ('Fire_ZTOP', 1, ' ')

       do n=1,shr_fire_emis_mechcomps_n

          name = shr_fire_emis_mechcomps(n)%name
          fire_frc_name(n)  = 'FireFrc_'//trim(name)
          fire_sflx_name(n) = 'FireSFLX_'//trim(name)
          fire_vflx_name(n) = 'FireVFLX_'//trim(name)

!          call addfld (fire_frc_name(n), (/'lev'/), 'A','molecules/cm^3/s', 'vertical fire emissions for '//trim(name))
!          call addfld (fire_sflx_name(n),horiz_only,'A','kg/m^2/s', 'surface fire emissions for '//trim(name))
!          call addfld (fire_vflx_name(n),horiz_only,'A','kg/m^2/s', 'vertically integrated fire emissions for '//trim(name))
!          call add_default (fire_frc_name(n), 1, ' ')
!          call add_default (fire_sflx_name(n), 1, ' ')
!          call add_default (fire_vflx_name(n), 1, ' ')

          frc_ndx = get_extfrc_ndx( name )
          spc_ndx = get_spc_ndx( name )

          if (frc_ndx>0 .and. spc_ndx>0) then
             frc_spc_map(n) = frc_ndx
             chm_spc_map(n) = spc_ndx
          else
             write(iulog,*) 'fire_emissions_init: not able to map '//trim(name)//' to chem species/forcing '
             write(iulog,*) 'fire_emissions_init:  ... frc_ndx = ',frc_ndx
             write(iulog,*) 'fire_emissions_init:  ... spc_ndx = ',spc_ndx
             call endrun('fire_emissions_init: not able to map '//trim(name)//' to chem species/forcing ')
          endif

          spc_mass_factor(n) = 1.e-6_r8 * avogad / adv_mass(spc_ndx) ! 1.e-6 converts m-3 to cm-3. 
          ! (molecules/kmole) / (g/mole) --> molecules/kg

          ! for MAM need to include cooresponding forcings of number densities 

          found = rad_cnst_num_name(0, name, num_name, mode_out=mode, spec_out=spec )

          if ( found ) then

             frc_ndx = get_extfrc_ndx( num_name )

             call rad_cnst_get_aer_props(0, mode, spec, density_aer=specdens)
             frc_num_map(n) = frc_ndx 
             num_mass_factor(n) = x_numfact / specdens

             fire_numfrc_name(n) = 'FireFrc_'//trim(name)//'_'//trim(num_name)
!             call addfld (fire_numfrc_name(n),(/'lev'/), 'A', 'molecules/cm^3/s', &
!                  'vertical fire emissions for '//trim(num_name)//' due to component '//trim(name))
!             call add_default (fire_numfrc_name(n), 1, ' ')
        endif

       enddo

    else ! initialize surface emissions

       allocate( fire_emis_indices_map(shr_fire_emis_mechcomps_n) )
       allocate( fire_emis_num_map(shr_fire_emis_mechcomps_n) )
       fire_emis_indices_map(:) = -1
       fire_emis_num_map(:) = -1

       do n=1,shr_fire_emis_mechcomps_n
          call cnst_get_ind (shr_fire_emis_mechcomps(n)%name,  fire_emis_indices_map(n), abrtf=.false. )

          ii = get_spc_ndx(shr_fire_emis_mechcomps(n)%name)
          if (ii<1) then
             call endrun('gas_phase_chemdr_inti: Fire emissions compound not in chemistry mechanism : '&
                  //trim(shr_fire_emis_mechcomps(n)%name))
          endif

          ! Fire emis history fields
          call addfld( 'FireSF_'//trim(shr_fire_emis_mechcomps(n)%name),horiz_only,'A','kg/m2/sec',&
               trim(shr_fire_emis_mechcomps(n)%name)//' Fire emissions flux')
          call add_default ('FireSF_'//trim(shr_fire_emis_mechcomps(n)%name), 1, ' ')

!added by LXu@05/26/20
          name = shr_fire_emis_mechcomps(n)%name
          found = rad_cnst_num_name(0, name, num_name, mode_out=mode, spec_out=spec )
!          if (masterproc) write(iulog,*) 'find fire aerosol number indices: ', name, num_name, found

          if ( found ) then

!             frc_ndx = get_extfrc_ndx( num_name )
             
!	     fire_emis_num_map(n) = frc_ndx 
	     
	     call cnst_get_ind (num_name,  fire_emis_num_map(n), abrtf=.false. )

          endif

       end do

!added by LXu@05/21+++++      
!      call cnst_get_ind ('num_a4',  fire_emis_indices_map(shr_fire_emis_mechcomps_n+1), abort=.false. )
!      ii = get_spc_ndx('num_a4')
!      if (ii<1) then
!         call endrun('gas_phase_chemdr_inti: Fire emissions compound not in chemistry mechanism : num_a4')
!      endif
!      if (masterproc) write(iulog,*) 'initializing fire num_a4', fire_emis_indices_map(shr_fire_emis_mechcomps_n+1), ii

!      call cnst_get_ind ('num_a1',  fire_emis_indices_map(shr_fire_emis_mechcomps_n+2), abort=.false. )
!      ii = get_spc_ndx('num_a1')
!      if (ii<1) then
!         call endrun('gas_phase_chemdr_inti: Fire emissions compound not in chemistry mechanism : num_a1')
!      endif

    endif

  end subroutine fire_emissions_init

  !------------------------------------------------------------------------------
  ! sets surface emissions 
  !------------------------------------------------------------------------------
  subroutine fire_emissions_srf( lchnk, ncol, fireflx, sflx )

    ! dummy args
    integer,  intent(in) :: lchnk, ncol
    real(r8), pointer, intent(in) :: fireflx(:,:)
    real(r8), intent(inout) :: sflx(:,:)

    ! local vars
    integer :: i, n
!added by LXu@05/20
    real(r8) :: x_mass_factor, x_mton
!    real(r8), parameter :: spec_density(6) = (/ 1.7e3_r8,  1.0e3_r8,     1.77e3_r8, 1.77e3_r8,  1.0e3_r8,  1.0e3_r8/) !kg/m3, BC/OC/SO4/SO2/SOAG/CO
!    real(r8), parameter :: spec_mweight(6) = (/12.011_r8, 12.011_r8, 115.107340_r8, 64.066_r8, 12.011_r8,   28.0_r8/)
!    real(r8), parameter :: spec_volavg_d(6) = (/0.134e-6_r8, 0.134e-6_r8, 0.134e-6_r8, 0.134e-6_r8, 0.134e-6_r8, 0.134e-6_r8/) !micronmeter
    real(r8), parameter :: spec_density(8) = (/ 1.7e3_r8,  1.0e3_r8,     1.77e3_r8, 1.77e3_r8,  1.0e3_r8,  1.0e3_r8, &
                                                1.823_r8,  1.823_r8/) !kg/m3, BC/OC/SO4/SO2/SOAG/CO/po4_a4/po4_a3
    real(r8), parameter :: spec_mweight(8) = (/12.011_r8, 12.011_r8, 115.107340_r8, 64.066_r8, 12.011_r8,   28.0_r8, &
                                                 95.0_r8,   95.0_r8/)
    real(r8), parameter :: spec_volavg_d(8) = (/0.134e-6_r8, 0.134e-6_r8, 0.134e-6_r8, 0.134e-6_r8, 0.134e-6_r8, 0.134e-6_r8, &
                                                0.134e-6_r8, 0.134e-6_r8/) !micronmeter
    real(r8), parameter :: demis_acc_2 = 0.134e-6_r8 ! meters 
!    real(r8), parameter :: x_numfact_2 = 1.e-4_r8 * avogad * 6.0_r8 / (pi*(demis_acc_2**3))   ! 1.e-4 converts m-2 to cm-2. 
    real(r8), parameter :: x_numfact_2 = 6.0_r8 / (pi*(demis_acc_2**3))   ! m-3

    
    if (masterproc) write(iulog,*) 'Fire Aerosol indices and number: ',index_x2a_Fall_flxfire, shr_fire_emis_mechcomps_n
    ! fire surface emissions if not elevated forcing
    if ((.not.shr_fire_emis_elevated) .and. index_x2a_Fall_flxfire>0 .and. shr_fire_emis_mechcomps_n>0 ) then

       ! set Fire Emis fluxes ( add to other emis )
       do i =1,ncol
          do n = 1,shr_fire_emis_mechcomps_n
             sflx(i,fire_emis_indices_map(n)) &
                  = sflx(i,fire_emis_indices_map(n)) + fireflx(i,n) 

!added by LXu@03/20+++++
!adjust number emission from bc and pom, so4
             ! fire species number emission (molecules/cm2/s)
!       sflx(:ncol) = sflx(:ncol) * 1.e4_r8 * adv_mass(chm_spc_map(n))/avogad ! molecules/cm2/sec --> kg/m2/sec
                                                                              ! / avogad     --> kmoles/cm2/sec
                                                                              ! * adv_mass   -->     kg/cm2/sec
                                                                              ! * 1.e4       -->     kg/ m2/sec
!             num_mass_factor(n) = x_numfact / specdens
!             fire_frc(:ncol,k) = fire_sflx(:ncol,n) * vertical_fire(:ncol,k) * num_mass_factor(n)  ! molecules/cm3/s 

!	     if ( n < 4 ) then
!	          if ( n < 3 ) then
!		  x_mton = 1.e-4_r8 * avogad * 6.0_r8 / (spec_density(n)*pi*(spec_volavg_d(n)**3))
!                  sflx(i,fire_emis_indices_map(shr_fire_emis_mechcomps_n+1)) &
!                       = sflx(i,fire_emis_indices_map(shr_fire_emis_mechcomps_n+1)) + fireflx(i,n) * x_numfact_2 / spec_density(n)
!		  else
!		  x_mton = 1.e-4_r8 * avogad * 6.0_r8 / (spec_density(n)*pi*(spec_volavg_d(n)**3))
!                  sflx(i,fire_emis_indices_map(shr_fire_emis_mechcomps_n+2)) &
!                       = sflx(i,fire_emis_indices_map(shr_fire_emis_mechcomps_n+2)) + fireflx(i,n) * x_numfact_2 / spec_density(n)
!		  end if
!	     end if	  
!added by LXu@03/20-----
!added by LXu@05/26/20+++++
!          frcing(:ncol,:,frc_num_map(n)) = frcing(:ncol,:,frc_num_map(n)) + fire_frc(:ncol,:)
             if (masterproc) write(iulog,*) 'Aerosol spc and number indices: ',fire_emis_indices_map(n), fire_emis_num_map(n)
	     if (fire_emis_num_map(n)>0) then
                  sflx(i,fire_emis_num_map(n)) &
                       = sflx(i,fire_emis_num_map(n)) + fireflx(i,n) * x_numfact_2 / spec_density(n)
	     end if
!added by LXu@05/26/20-----
          end do
       end do

       ! output fire emis fluxes to history
       do n = 1,shr_fire_emis_mechcomps_n
          if (masterproc) write(iulog,*) 'Aerosol spc and number indices: ',fire_emis_indices_map(n), fire_emis_num_map(n)
          call outfld('FireSF_'//trim(shr_fire_emis_mechcomps(n)%name), fireflx(:ncol,n), ncol, lchnk)
       enddo

    endif

  end subroutine fire_emissions_srf

  !------------------------------------------------------------------------------
  ! sets vertical emissions (forcings)
  ! vertically distributes wild fire emissions 
  !------------------------------------------------------------------------------
  subroutine fire_emissions_vrt( ncol, lchnk, zint, fire_sflx, fire_ztop, frcing )

   ! args
    integer,          intent(in) :: ncol,lchnk
    real(r8),         intent(in) :: zint(:,:)      ! interface geopot above surface (km)
    real(r8),pointer, intent(in) :: fire_sflx(:,:) ! fire surface emissions (kg/m2/sec)
    real(r8),pointer, intent(in) :: fire_ztop(:)   ! top of vert distribution of fire surface emissions (m)
    real(r8),      intent(inout) :: frcing(:,:,:)  ! insitu forcings (molecules/cm3/sec)

    ! local vars
    real(r8) :: vertical_fire(ncol,pver), ztop
    integer  :: n, i,k
    real(r8) :: fire_frc(ncol,pver)
    real(r8) :: sflx(ncol)

    if (.not.shr_fire_emis_elevated) return
    if (shr_fire_emis_mechcomps_n<1) return

    !   define vertical_fire from Dentener  units /m
    do k=1,pver
       do i=1,ncol
          ztop = fire_ztop(i)*1.e-3_r8 ! convert m to km
          if(zint(i,k)<ztop)then
             vertical_fire(i,k)=1.e-3_r8/(ztop-zint(i,pverp))
          elseif(zint(i,k)>ztop.and.zint(i,k+1)<ztop)then
             vertical_fire(i,k)=1.e-3_r8*(ztop-zint(i,k+1))/(zint(i,k)-zint(i,k+1))/(ztop-zint(i,pverp))
          else
             vertical_fire(i,k)=0._r8
          endif
       end do
    end do

    ! extend the surface emissions in the vertical ....

    do n=1,shr_fire_emis_mechcomps_n

       fire_frc(:,:) = 0._r8

       do k = 1,pver !        (kg/m2/sec)        (/m)                    (molecules/kg)*1.e-6   --> molecules/cm3/s
          fire_frc(:ncol,k) = fire_sflx(:ncol,n) * vertical_fire(:ncol,k) * spc_mass_factor(n)  ! molecules/cm3/s
       enddo
!       call outfld( fire_frc_name(n), fire_frc, ncol, lchnk )
       frcing(:ncol,:,frc_spc_map(n)) = frcing(:ncol,:,frc_spc_map(n)) + fire_frc(:ncol,:)

       ! for debugging ...
       ! vertical intergration of the forcing should get back the surface flux
       sflx(:) = 0._r8
       do k = 1,pver
          sflx(:ncol) = sflx(:ncol) + 1.e5_r8*(zint(:ncol,k)-zint(:ncol,k+1))*fire_frc(:ncol,k) ! molecules/cm3/s --> molecules/cm2/sec 
       enddo
       sflx(:ncol) = sflx(:ncol) * 1.e4_r8 * adv_mass(chm_spc_map(n))/avogad ! molecules/cm2/sec --> kg/m2/sec
                                                                             ! / avogad     --> kmoles/cm2/sec
                                                                             ! * adv_mass   -->     kg/cm2/sec
                                                                             ! * 1.e4       -->     kg/ m2/sec
!       call outfld( fire_vflx_name(n),      sflx(:ncol  ), ncol, lchnk )
!       call outfld( fire_sflx_name(n), fire_sflx(:ncol,n), ncol, lchnk )

       ! for MAM need to include corresponding forcings of number densities 
       if (frc_num_map(n)>0) then
          do k = 1,pver
             fire_frc(:ncol,k) = fire_sflx(:ncol,n) * vertical_fire(:ncol,k) * num_mass_factor(n)  ! molecules/cm3/s 
          enddo
!          call outfld( fire_numfrc_name(n), fire_frc, ncol, lchnk )
          frcing(:ncol,:,frc_num_map(n)) = frcing(:ncol,:,frc_num_map(n)) + fire_frc(:ncol,:)
       endif
       
    enddo

!    call outfld( 'Fire_ZTOP', fire_ztop(:ncol), ncol, lchnk )

  end subroutine fire_emissions_vrt

end module fire_emissions
