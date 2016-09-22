! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
#define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_singles(dminfo,size(A),A)

!=================================================================================================================
 module mpas_atmphys_init_microphysics
 use mpas_dmpar
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_utilities
 use module_mp_thompson, only: is_aerosol_aware,naCCN0,naCCN1,naIN0,naIN1,ntb_arc,ntb_arw,ntb_art,ntb_arr, &
                               ntb_ark,tnccn_act

 implicit none
 private
 public:: init_thompson_clouddroplets_forMPAS,    &
          init_thompson_aerosols_forMPAS

! MPAS main initialization of the Thompson parameterization of cloud microphysics with nucleation of cloud
! droplets based on distributions of CCNs and INs (aerosol-aware parameterization).
! Laura D. Fowler (send comments to laura@ucar.edu).
! 2016-03-28.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added "use mpas_dmpar" at the top of the module.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-04.
! * modified the initialization of nifa and nwfa.If nifa and nwfa are already available in the initial conditions
!   using the climatological GOCART data,do not recalculate nifa and nwfa using an exponential profile of CCN and
!   IN as a function of height.
!   Laura D. Fowler (laura@ucar.edu) / 2016-05-27.


 contains


!=================================================================================================================
 subroutine init_thompson_clouddroplets_forMPAS(mesh,sfc_input,diag_physics)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input

!inout variables:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer:: nt_c,mu_c,nwfa2d

 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine init_thompson_clouddroplets_forMPAS:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(sfc_input,'landmask',landmask)
 call mpas_pool_get_array(diag_physics,'nt_c',nt_c)
 call mpas_pool_get_array(diag_physics,'mu_c',mu_c)

!... initialize the prescribed number of cloud droplets, and mu_c (parameter in the exponential of the generalized
!gamma distribution) as a function of the land-cean mask. as set in the thompson cloud microphysics scheme, nt_c
!is set to 100 per cc (100.E6 m^-3) for maritime cases and 300 per cc (300.E6 m^-3) for continental cases.
 do iCell = 1, nCellsSolve
    if(landmask(iCell) .eq. 1) then
       nt_c(iCell) = 300.e6 
    elseif(landmask(iCell) .eq. 0) then
       nt_c(iCell) = 100.e6
    endif
    mu_c(iCell) = MIN(15., (1000.e6/nt_c(iCell) + 2.))
 enddo

!write(0,*) '--- end subroutine init_thompson_clouddroplets_forMPAS.'

 end subroutine init_thompson_clouddroplets_forMPAS 

!=================================================================================================================
 subroutine init_thompson_aerosols_forMPAS(dminfo,mesh,state,time_lev,diag_physics)
!=================================================================================================================

!input variables:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 integer,intent(in):: time_lev

!inout variables:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: state

!local variables and pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 integer,pointer:: index_nifa,index_nwfa

 real(kind=RKIND),dimension(:),pointer    :: areaCell
 real(kind=RKIND),dimension(:),pointer    :: nwfa2d
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid
 real(kind=RKIND),dimension(:,:),pointer  :: nifa,nwfa
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 character(len=StrKIND):: mess

 integer:: iCell, k

 real(kind=RKIND):: max_test
 real(kind=RKIND):: h_01
 real(kind=RKIND):: niIN3,niCCN3
 real(kind=RKIND):: nifa_max,nifa_min,global_nifa_max,global_nifa_min
 real(kind=RKIND):: nwfa_max,nwfa_min,global_nwfa_max,global_nwfa_min
 real(kind=RKIND),dimension(:,:),allocatable:: hgt

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine init_thompson_aerosols_forMPAS:'

 is_aerosol_aware = .true.

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(mesh,'areaCell',areaCell)
 call mpas_pool_get_array(mesh,'zgrid'   ,zgrid   )

 call mpas_pool_get_array(diag_physics,'nwfa2d',nwfa2d)

 call mpas_pool_get_dimension(state,'index_nifa'  ,index_nifa  )
 call mpas_pool_get_dimension(state,'index_nwfa'  ,index_nwfa  )

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 nifa   => scalars(index_nifa,:,:)
 nwfa   => scalars(index_nwfa,:,:)

 if(.not.allocated(hgt)) allocate(hgt(1:nVertLevels,1:nCellsSolve))
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       hgt(k,iCell) = 0.5_RKIND * (zgrid(k,iCell)+zgrid(k+1,iCell))
    enddo
 enddo

!... initialize the distribution of hygroscopic ("water friendly") aerosols if not already initialized using
!    GOCART data:
 global_nwfa_min = 0._RKIND
 global_nwfa_max = 0._RKIND
 nwfa_min = minval(nwfa(:,1:nCellsSolve))
 nwfa_max = maxval(nwfa(:,1:nCellsSolve))
 call mpas_dmpar_min_real(dminfo,nwfa_min,global_nwfa_min)
 call mpas_dmpar_max_real(dminfo,nwfa_max,global_nwfa_max)
 write(0,*) '--- global_nwfa_min = ', global_nwfa_min
 write(0,*) '--- global_nwfa_max = ', global_nwfa_max

 if(global_nwfa_min == 0._RKIND .and. global_nwfa_max == 0._RKIND) then
    write(mess,*) '--- initialize nwfa using an exponential distribution of CCN as a function of height.'
    call physics_message(mess)
    do iCell = 1, nCellsSolve
       if(hgt(1,iCell).le.1000.0) then
          h_01 = 0.8
       elseif(hgt(1,iCell).ge.2500.0) then
          h_01 = 0.01
       else
          h_01 = 0.8*cos(hgt(1,iCell)*0.001 - 1.0)
       endif
       niCCN3 = -1.0*ALOG(naCCN1/naCCN0)/h_01
       nwfa(1,iCell) = naCCN1+naCCN0*exp(-((hgt(2,iCell)-hgt(1,iCell))/1000.)*niCCN3)
       do k = 2, nVertLevels
          nwfa(k,iCell) = naCCN1+naCCN0*exp(-((hgt(k,iCell)-hgt(1,iCell))/1000.)*niCCN3)
       enddo
    enddo
 else
    write(mess,*) '--- initialize nwfa using the climatological GOCART data.'
    call physics_message(mess)
 endif

!... initialize the distribution of nonhygroscopic ("ice friendly") aerosols if not already initialized using
!    GOCART data:
 global_nifa_min = 0._RKIND
 global_nifa_max = 0._RKIND
 nifa_min = minval(nifa(:,1:nCellsSolve))
 nifa_max = maxval(nifa(:,1:nCellsSolve))
 call mpas_dmpar_min_real(dminfo,nifa_min,global_nifa_min)
 call mpas_dmpar_max_real(dminfo,nifa_max,global_nifa_max)
 write(0,*) '--- global_nifa_min = ', global_nifa_min
 write(0,*) '--- global_nifa_max = ', global_nifa_max

 if(global_nwfa_min == 0._RKIND .and. global_nwfa_max == 0._RKIND) then
    write(mess,*) '--- initialize nifa using an exponential distribution of IN as a function of height.'
    do iCell = 1, nCellsSolve
       if(hgt(1,iCell).le.1000.0) then
          h_01 = 0.8
       elseif(hgt(1,iCell).ge.2500.0) then
          h_01 = 0.01
       else
          h_01 = 0.8*cos(hgt(1,iCell)*0.001 - 1.0)
       endif
       niIN3 = -1.0*ALOG(naIN1/naIN0)/h_01
       nifa(1,iCell) = naIN1+naIN0*exp(-((hgt(2,iCell)-hgt(1,iCell))/1000.)*niIN3)
       do k = 2, nVertLevels
          nifa(k,iCell) = naIN1+naIN0*exp(-((hgt(k,iCell)-hgt(1,iCell))/1000.)*niIN3)
       enddo
    enddo
 else
    write(mess,*) '--- initialize nifa using the climatological GOCART data.'
    call physics_message(mess)
 endif

!do iCell = 1, nCellsSolve
!   do k = 1, nVertLevels
!      write(0,202) iCell,k,hgt(k,iCell),nwfa(k,iCell),nifa(k,iCell)
!   enddo
!   write(0,*)
!enddo
!202 format(2i8,10(1x,e15.8))
!stop

!... capture initial state lowest level CCN aerosol data in 2D array:
!do iCell = 1,nCellsSolve
!   nwfa2d(iCell) = nwfa(1,iCell)
!enddo

!... scale the lowest level aerosol data into an emissions rate.  This is very far from ideal, but need higher
!... emissions where larger amount of existing and lesser emissions where not already lots of aerosols for first
!... -order simplistic approach. Later, proper connection to emission inventory would be better, but, for now,
!... scale like this:
!... where: Nwfa=50 per cc, emit 0.875E4 aerosols per kg per second
!...        Nwfa=500 per cc, emit 0.875E5 aerosols per kg per second
!...        Nwfa=5000 per cc, emit 0.875E6 aerosols per kg per second
!... for a grid with 20km spacing and scale accordingly for other spacings.

 do iCell = 1, nCellsSolve
    if(sqrt(areaCell(iCell))/20000.0 .ge. 1.0) then
       h_01 = 0.875
    else
       h_01 = (0.875 + 0.125*((20000.-sqrt(areaCell(iCell)))/16000.)) * sqrt(areaCell(iCell))/20000.
    endif
    nwfa2d(iCell) = 10.0**(LOG10(nwfa(1,iCell)*1.E-6)-3.69897)
    nwfa2d(iCell) = nwfa2d(iCell)*h_01*1.E6
 enddo

!... read a static file containing CCN activation of aerosols. The data were created from a parcel model by
!... Feingold & Heymsfield with further changes by Eidhammer and Kriedenweis.
 call table_ccnAct(dminfo)
 write(0,*) '--- end read table_ccnAct:'

!... deallocate local arrays:
 if(allocated(hgt)) deallocate(hgt)

 write(0,*) '--- end subroutine init_thompson_aerosols_forMPAS.'

 end subroutine init_thompson_aerosols_forMPAS

!=================================================================================================================
 subroutine table_ccnAct(dminfo)
!=================================================================================================================

!input variables:
 type(dm_info),intent(in):: dminfo

!local variables:
 logical:: opened
 integer:: ccn_unit,i,istat
 character(len=StrKIND):: errmess
!-----------------------------------------------------------------------------------------------------------------

 if(.not.allocated(tnccn_act)) allocate(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark))

!get a unit to open binary file:
 istat = -999
 if(dminfo % my_proc_id == IO_NODE) then
    do i = 10,99
       inquire(i,opened = opened,iostat=istat)
       if(.not. opened ) then
          ccn_unit = i
          exit
       endif
    enddo
    if(istat /= 0) &
       call physics_error_fatal('mpas_atmphys_init_microphysics table_ccnAct: Can not '// &
                                'find unused fortran unit to read in lookup table.' )
 endif

!distribute unit to other processors:
 call mpas_dmpar_bcast_int(dminfo,ccn_unit)

!open binary file:
 istat = -999
 if(dminfo % my_proc_id == IO_NODE) then
    open(ccn_unit,file='CCN_ACTIVATE_DATA',form='UNFORMATTED',status='OLD',iostat=istat)
    if(istat /= 0) then
       write(errmess,'(A,I4)') 'mpas_atmphys_init_microphysics table_ccnAct:: '// &
                                'error opening CCN_ACTIVATE_DATA on unit', ccn_unit
       call physics_error_fatal(errmess)
    endif
 endif

!read and broadcast data to all nodes:
 istat = -999
 if(dminfo % my_proc_id == IO_NODE) then
    read(ccn_unit,iostat=istat) tnccn_act
    if(istat /= 0) then
       write(errmess,'(A,I4)') 'mpas_atmphys_init_microphysics table_ccnAct:: '// &
                                'error reading tnccn_act on unit', ccn_unit
       call physics_error_fatal(errmess)
    endif
 endif
 
 DM_BCAST_MACRO(tnccn_act)

 end subroutine table_ccnAct

!=================================================================================================================
 end module mpas_atmphys_init_microphysics
!=================================================================================================================

 
 
