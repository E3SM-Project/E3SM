module initVerticalMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initialize vertical components of column datatype
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_sys_mod    , only : shr_sys_abort
  use decompMod      , only : bounds_type
  use spmdMod        , only : masterproc
  use clm_varpar     , only : more_vertlayers, nlevsno, nlevgrnd, nlevlak
  use clm_varpar     , only : toplev_equalspace, nlev_equalspace
  use clm_varpar     , only : nlevsoi, nlevsoifl, nlevurb, nlevslp 
  use clm_varctl     , only : fsurdat, iulog, use_var_soil_thick
  use clm_varctl     , only : use_vancouver, use_mexicocity, use_vertsoilc, use_extralakelayers
  use clm_varctl     , only : use_erosion
  use elm_varcon     , only : zlak, dzlak, zsoi, dzsoi, zisoi, dzsoi_decomp, spval, grlnd 
  use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv
  use landunit_varcon, only : istdlak, istice_mec
  use fileutils      , only : getfil
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  use ncdio_pio
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initVertical
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine initVertical(bounds, snow_depth, thick_wall, thick_roof)
    !
    ! !ARGUMENTS:
    type(bounds_type)   , intent(in)    :: bounds
    real(r8)            , intent(in)    :: snow_depth(bounds%begc:)
    real(r8)            , intent(in)    :: thick_wall(bounds%begl:)
    real(r8)            , intent(in)    :: thick_roof(bounds%begl:)
    !
    ! LOCAL VARAIBLES:
    integer               :: c,l,g,i,j,lev     ! indices 
    type(file_desc_t)     :: ncid              ! netcdf id
    logical               :: readvar 
    integer               :: dimid             ! dimension id
    character(len=256)    :: locfn             ! local filename
    real(r8) ,pointer     :: std (:)           ! read in - topo_std 
    real(r8) ,pointer     :: tslope (:)        ! read in - topo_slope 
    real(r8) ,pointer     :: hslp_p10 (:,:)    ! read in - hillslope slope percentiles
    real(r8) ,pointer     :: dtb (:)           ! read in - DTB
    real(r8)              :: beddep            ! temporary
    integer               :: nlevbed           ! temporary
    real(r8)              :: zimid             ! temporary
    real(r8)              :: slope0            ! temporary
    real(r8)              :: slopebeta         ! temporary
    real(r8)              :: slopemax          ! temporary
    integer               :: ier               ! error status
    real(r8)              :: scalez = 0.025_r8 ! Soil layer thickness discretization (m)
    real(r8)              :: thick_equal = 0.2
    real(r8) ,pointer     :: lakedepth_in(:)   ! read in - lakedepth 
    real(r8), allocatable :: zurb_wall(:,:)    ! wall (layer node depth)
    real(r8), allocatable :: zurb_roof(:,:)    ! roof (layer node depth)
    real(r8), allocatable :: dzurb_wall(:,:)   ! wall (layer thickness)
    real(r8), allocatable :: dzurb_roof(:,:)   ! roof (layer thickness)
    real(r8), allocatable :: ziurb_wall(:,:)   ! wall (layer interface)
    real(r8), allocatable :: ziurb_roof(:,:)   ! roof (layer interface)
    real(r8)              :: depthratio        ! ratio of lake depth to standard deep lake depth 
    integer               :: begc, endc
    integer               :: begl, endl
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl

    SHR_ASSERT_ALL((ubound(snow_depth)  == (/endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(thick_wall)  == (/endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(thick_roof)  == (/endl/)), errMsg(__FILE__, __LINE__))

    ! Open surface dataset to read in data below 

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_inqdlen(ncid, dimid, nlevsoifl, name='nlevsoi')
    if ( .not. more_vertlayers )then
       if ( nlevsoifl /= nlevsoi )then
          call shr_sys_abort(' ERROR: Number of soil layers on file does NOT match the number being used'//&
               errMsg(__FILE__, __LINE__))
       end if
    else
       ! read in layers, interpolate to high resolution grid later
    end if

    ! --------------------------------------------------------------------
    ! Define layer structure for soil, lakes, urban walls and roof 
    ! Vertical profile of snow is not initialized here - but below
    ! --------------------------------------------------------------------
    
    ! Soil layers and interfaces (assumed same for all non-lake patches)
    ! "0" refers to soil surface and "nlevsoi" refers to the bottom of model soil
    
    if ( more_vertlayers )then
       ! replace standard exponential grid with a grid that starts out exponential, 
       ! then has several evenly spaced layers, then finishes off exponential. 
       ! this allows the upper soil to behave as standard, but then continues 
       ! with higher resolution to a deeper depth, so that, for example, permafrost
       ! dynamics are not lost due to an inability to resolve temperature, moisture, 
       ! and biogeochemical dynamics at the base of the active layer
       do j = 1, toplev_equalspace
          zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
       enddo

       do j = toplev_equalspace+1,toplev_equalspace + nlev_equalspace
          zsoi(j) = zsoi(j-1) + thick_equal
       enddo

       do j = toplev_equalspace + nlev_equalspace +1, nlevgrnd
          zsoi(j) = scalez*(exp(0.5_r8*((j - nlev_equalspace)-0.5_r8))-1._r8) + nlev_equalspace * thick_equal
       enddo
    else

       do j = 1, nlevgrnd
          zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
       enddo
    end if

    dzsoi(1) = 0.5_r8*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
    do j = 2,nlevgrnd-1
       dzsoi(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
    enddo
    dzsoi(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)

    zisoi(0) = 0._r8
    do j = 1, nlevgrnd-1
       zisoi(j) = 0.5_r8*(zsoi(j)+zsoi(j+1))         !interface depths
    enddo
    zisoi(nlevgrnd) = zsoi(nlevgrnd) + 0.5_r8*dzsoi(nlevgrnd)

    if (masterproc) then
       write(iulog, *) 'zsoi', zsoi(:) 
       write(iulog, *) 'zisoi: ', zisoi(:)
       write(iulog, *) 'dzsoi: ', dzsoi(:)
    end if

    ! define a vertical grid spacing such that it is the normal dzsoi if nlevdecomp =nlevgrnd, or else 1 meter
    if (use_vertsoilc) then
       dzsoi_decomp(1) = 0.5_r8*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
       do j = 2,nlevgrnd-1
          dzsoi_decomp(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
       enddo
       dzsoi_decomp(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)
    else
       dzsoi_decomp(1) = 1.
    end if
    if (masterproc) then
       write(iulog, *) 'dzsoi_decomp', dzsoi_decomp(:) 
    end if

    if (nlevurb > 0) then
       allocate(zurb_wall(bounds%begl:bounds%endl,nlevurb),    &
            zurb_roof(bounds%begl:bounds%endl,nlevurb),    &
            dzurb_wall(bounds%begl:bounds%endl,nlevurb),   &
            dzurb_roof(bounds%begl:bounds%endl,nlevurb),   &
            ziurb_wall(bounds%begl:bounds%endl,0:nlevurb), &
            ziurb_roof(bounds%begl:bounds%endl,0:nlevurb), &
            stat=ier)
       if (ier /= 0) then
          call shr_sys_abort(' ERROR allocation error for '//&
               'zurb_wall,zurb_roof,dzurb_wall,dzurb_roof,ziurb_wall,ziurb_roof'//&
               errMsg(__FILE__, __LINE__))
       end if
    end if


!    ! (FATES-INTERF)  Move this call to inside alm_fates%init()
!    if(use_fates) then
!       call ed_hist_scpfmaps()
!    end if

    ! Column level initialization for urban wall and roof layers and interfaces
    do l = bounds%begl,bounds%endl

       ! "0" refers to urban wall/roof surface and "nlevsoi" refers to urban wall/roof bottom
       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then       
             zurb_wall(l,1) = 0.010_r8/2._r8
             zurb_wall(l,2) = zurb_wall(l,1) + 0.010_r8/2._r8 + 0.020_r8/2._r8
             zurb_wall(l,3) = zurb_wall(l,2) + 0.020_r8/2._r8 + 0.070_r8/2._r8
             zurb_wall(l,4) = zurb_wall(l,3) + 0.070_r8/2._r8 + 0.070_r8/2._r8
             zurb_wall(l,5) = zurb_wall(l,4) + 0.070_r8/2._r8 + 0.030_r8/2._r8

             zurb_roof(l,1) = 0.010_r8/2._r8
             zurb_roof(l,2) = zurb_roof(l,1) + 0.010_r8/2._r8 + 0.010_r8/2._r8
             zurb_roof(l,3) = zurb_roof(l,2) + 0.010_r8/2._r8 + 0.010_r8/2._r8
             zurb_roof(l,4) = zurb_roof(l,3) + 0.010_r8/2._r8 + 0.010_r8/2._r8
             zurb_roof(l,5) = zurb_roof(l,4) + 0.010_r8/2._r8 + 0.030_r8/2._r8

             dzurb_wall(l,1) = 0.010_r8
             dzurb_wall(l,2) = 0.020_r8
             dzurb_wall(l,3) = 0.070_r8
             dzurb_wall(l,4) = 0.070_r8
             dzurb_wall(l,5) = 0.030_r8
             write(iulog,*)'Total thickness of wall: ',sum(dzurb_wall(l,:))
             write(iulog,*)'Wall layer thicknesses: ',dzurb_wall(l,:)

             dzurb_roof(l,1) = 0.010_r8
             dzurb_roof(l,2) = 0.010_r8
             dzurb_roof(l,3) = 0.010_r8
             dzurb_roof(l,4) = 0.010_r8
             dzurb_roof(l,5) = 0.030_r8
             write(iulog,*)'Total thickness of roof: ',sum(dzurb_roof(l,:))
             write(iulog,*)'Roof layer thicknesses: ',dzurb_roof(l,:)

             ziurb_wall(l,0) = 0.
             ziurb_wall(l,1) = dzurb_wall(l,1)
             do j = 2,nlevurb
                ziurb_wall(l,j) = sum(dzurb_wall(l,1:j))
             end do
             write(iulog,*)'Wall layer interface depths: ',ziurb_wall(l,:)

             ziurb_roof(l,0) = 0.
             ziurb_roof(l,1) = dzurb_roof(l,1)
             do j = 2,nlevurb
                ziurb_roof(l,j) = sum(dzurb_roof(l,1:j))
             end do
             write(iulog,*)'Roof layer interface depths: ',ziurb_roof(l,:)
          else if (use_mexicocity) then
             zurb_wall(l,1) = 0.015_r8/2._r8
             zurb_wall(l,2) = zurb_wall(l,1) + 0.015_r8/2._r8 + 0.120_r8/2._r8
             zurb_wall(l,3) = zurb_wall(l,2) + 0.120_r8/2._r8 + 0.150_r8/2._r8
             zurb_wall(l,4) = zurb_wall(l,3) + 0.150_r8/2._r8 + 0.150_r8/2._r8
             zurb_wall(l,5) = zurb_wall(l,4) + 0.150_r8/2._r8 + 0.015_r8/2._r8

             zurb_roof(l,1) = 0.010_r8/2._r8
             zurb_roof(l,2) = zurb_roof(l,1) + 0.010_r8/2._r8 + 0.050_r8/2._r8
             zurb_roof(l,3) = zurb_roof(l,2) + 0.050_r8/2._r8 + 0.050_r8/2._r8
             zurb_roof(l,4) = zurb_roof(l,3) + 0.050_r8/2._r8 + 0.050_r8/2._r8
             zurb_roof(l,5) = zurb_roof(l,4) + 0.050_r8/2._r8 + 0.025_r8/2._r8

             dzurb_wall(l,1) = 0.015_r8
             dzurb_wall(l,2) = 0.120_r8
             dzurb_wall(l,3) = 0.150_r8
             dzurb_wall(l,4) = 0.150_r8
             dzurb_wall(l,5) = 0.015_r8
             write(iulog,*)'Total thickness of wall: ',sum(dzurb_wall(l,:))
             write(iulog,*)'Wall layer thicknesses: ',dzurb_wall(l,:)

             dzurb_roof(l,1) = 0.010_r8
             dzurb_roof(l,2) = 0.050_r8
             dzurb_roof(l,3) = 0.050_r8
             dzurb_roof(l,4) = 0.050_r8
             dzurb_roof(l,5) = 0.025_r8
             write(iulog,*)'Total thickness of roof: ',sum(dzurb_roof(l,:))
             write(iulog,*)'Roof layer thicknesses: ',dzurb_roof(l,:)

             ziurb_wall(l,0) = 0.
             ziurb_wall(l,1) = dzurb_wall(l,1)
             do j = 2,nlevurb
                ziurb_wall(l,j) = sum(dzurb_wall(l,1:j))
             end do
             write(iulog,*)'Wall layer interface depths: ',ziurb_wall(l,:)

             ziurb_roof(l,0) = 0.
             ziurb_roof(l,1) = dzurb_roof(l,1)
             do j = 2,nlevurb
                ziurb_roof(l,j) = sum(dzurb_roof(l,1:j))
             end do
             write(iulog,*)'Roof layer interface depths: ',ziurb_roof(l,:)
          else
             do j = 1, nlevurb
                zurb_wall(l,j) = (j-0.5)*(thick_wall(l)/float(nlevurb))  !node depths
             end do
             do j = 1, nlevurb
                zurb_roof(l,j) = (j-0.5)*(thick_roof(l)/float(nlevurb))  !node depths
             end do

             dzurb_roof(l,1) = 0.5*(zurb_roof(l,1)+zurb_roof(l,2))    !thickness b/n two interfaces
             do j = 2,nlevurb-1
                dzurb_roof(l,j)= 0.5*(zurb_roof(l,j+1)-zurb_roof(l,j-1)) 
             enddo
             dzurb_roof(l,nlevurb) = zurb_roof(l,nlevurb)-zurb_roof(l,nlevurb-1)

             dzurb_wall(l,1) = 0.5*(zurb_wall(l,1)+zurb_wall(l,2))    !thickness b/n two interfaces
             do j = 2,nlevurb-1
                dzurb_wall(l,j)= 0.5*(zurb_wall(l,j+1)-zurb_wall(l,j-1)) 
             enddo
             dzurb_wall(l,nlevurb) = zurb_wall(l,nlevurb)-zurb_wall(l,nlevurb-1)

             ziurb_wall(l,0) = 0.
             do j = 1, nlevurb-1
                ziurb_wall(l,j) = 0.5*(zurb_wall(l,j)+zurb_wall(l,j+1))          !interface depths
             enddo
             ziurb_wall(l,nlevurb) = zurb_wall(l,nlevurb) + 0.5*dzurb_wall(l,nlevurb)

             ziurb_roof(l,0) = 0.
             do j = 1, nlevurb-1
                ziurb_roof(l,j) = 0.5*(zurb_roof(l,j)+zurb_roof(l,j+1))          !interface depths
             enddo
             ziurb_roof(l,nlevurb) = zurb_roof(l,nlevurb) + 0.5*dzurb_roof(l,nlevurb)
          end if
       end if
    end do

    do c = bounds%begc,bounds%endc
       l = col_pp%landunit(c)

       if (lun_pp%urbpoi(l)) then
          if (col_pp%itype(c)==icol_sunwall .or. col_pp%itype(c)==icol_shadewall) then
             col_pp%z(c,1:nlevurb)  = zurb_wall(l,1:nlevurb)
             col_pp%zi(c,0:nlevurb) = ziurb_wall(l,0:nlevurb)
             col_pp%dz(c,1:nlevurb) = dzurb_wall(l,1:nlevurb)
             if (nlevurb < nlevgrnd) then
                col_pp%z(c,nlevurb+1:nlevgrnd)  = spval
                col_pp%zi(c,nlevurb+1:nlevgrnd) = spval
                col_pp%dz(c,nlevurb+1:nlevgrnd) = spval
             end if
          else if (col_pp%itype(c)==icol_roof) then
             col_pp%z(c,1:nlevurb)  = zurb_roof(l,1:nlevurb)
             col_pp%zi(c,0:nlevurb) = ziurb_roof(l,0:nlevurb)
             col_pp%dz(c,1:nlevurb) = dzurb_roof(l,1:nlevurb)
             if (nlevurb < nlevgrnd) then
                col_pp%z(c,nlevurb+1:nlevgrnd)  = spval
                col_pp%zi(c,nlevurb+1:nlevgrnd) = spval
                col_pp%dz(c,nlevurb+1:nlevgrnd) = spval
             end if
          else
             col_pp%z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
             col_pp%zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
             col_pp%dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
          end if
       else if (lun_pp%itype(l) /= istdlak) then
          col_pp%z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
          col_pp%zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
          col_pp%dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
       end if
    end do

    if (nlevurb > 0) then
       deallocate(zurb_wall, zurb_roof, dzurb_wall, dzurb_roof, ziurb_wall, ziurb_roof)
    end if

    !-----------------------------------------------
    ! Set lake levels and layers (no interfaces)
    !-----------------------------------------------

    allocate(lakedepth_in(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='LAKEDEPTH', flag='read', data=lakedepth_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          write(iulog,*) 'WARNING:: LAKEDEPTH not found on surface data set. All lake columns will have lake depth', &
               ' set equal to default value.'
       end if
       lakedepth_in(:) = spval
    end if
    do c = begc, endc
       g = col_pp%gridcell(c)
       col_pp%lakedepth(c) = lakedepth_in(g)
    end do
    deallocate(lakedepth_in)

    ! Lake layers
    if (.not. use_extralakelayers) then
       dzlak(1) = 0.1_r8
       dzlak(2) = 1._r8
       dzlak(3) = 2._r8
       dzlak(4) = 3._r8
       dzlak(5) = 4._r8
       dzlak(6) = 5._r8
       dzlak(7) = 7._r8
       dzlak(8) = 7._r8
       dzlak(9) = 10.45_r8
       dzlak(10)= 10.45_r8

       zlak(1) =  0.05_r8
       zlak(2) =  0.6_r8
       zlak(3) =  2.1_r8
       zlak(4) =  4.6_r8
       zlak(5) =  8.1_r8
       zlak(6) = 12.6_r8
       zlak(7) = 18.6_r8
       zlak(8) = 25.6_r8
       zlak(9) = 34.325_r8
       zlak(10)= 44.775_r8
    else
       dzlak(1) =0.1_r8
       dzlak(2) =0.25_r8
       dzlak(3) =0.25_r8
       dzlak(4) =0.25_r8
       dzlak(5) =0.25_r8
       dzlak(6) =0.5_r8
       dzlak(7) =0.5_r8
       dzlak(8) =0.5_r8
       dzlak(9) =0.5_r8
       dzlak(10) =0.75_r8
       dzlak(11) =0.75_r8
       dzlak(12) =0.75_r8
       dzlak(13) =0.75_r8
       dzlak(14) =2_r8
       dzlak(15) =2_r8
       dzlak(16) =2.5_r8
       dzlak(17) =2.5_r8
       dzlak(18) =3.5_r8
       dzlak(19) =3.5_r8
       dzlak(20) =3.5_r8
       dzlak(21) =3.5_r8
       dzlak(22) =5.225_r8
       dzlak(23) =5.225_r8
       dzlak(24) =5.225_r8
       dzlak(25) =5.225_r8

       zlak(1) = dzlak(1)/2._r8
       do i=2,nlevlak
          zlak(i) = zlak(i-1) + (dzlak(i-1)+dzlak(i))/2._r8
       end do
    end if

    do c = bounds%begc,bounds%endc
       l = col_pp%landunit(c)

       if (lun_pp%itype(l) == istdlak) then

          if (col_pp%lakedepth(c) == spval) then
             col_pp%lakedepth(c)         = zlak(nlevlak) + 0.5_r8*dzlak(nlevlak)
             col_pp%z_lake(c,1:nlevlak)  = zlak(1:nlevlak)
             col_pp%dz_lake(c,1:nlevlak) = dzlak(1:nlevlak)

          else if (col_pp%lakedepth(c) > 1._r8 .and. col_pp%lakedepth(c) < 5000._r8) then

             depthratio                 = col_pp%lakedepth(c) / (zlak(nlevlak) + 0.5_r8*dzlak(nlevlak)) 
             col_pp%z_lake(c,1)            = zlak(1)
             col_pp%dz_lake(c,1)           = dzlak(1)
             col_pp%dz_lake(c,2:nlevlak-1) = dzlak(2:nlevlak-1)*depthratio
             col_pp%dz_lake(c,nlevlak)     = dzlak(nlevlak)*depthratio - (col_pp%dz_lake(c,1) - dzlak(1)*depthratio)
             do lev=2,nlevlak
                col_pp%z_lake(c,lev) = col_pp%z_lake(c,lev-1) + (col_pp%dz_lake(c,lev-1)+col_pp%dz_lake(c,lev))/2._r8
             end do

          else if (col_pp%lakedepth(c) > 0._r8 .and. col_pp%lakedepth(c) <= 1._r8) then

             col_pp%dz_lake(c,:) = col_pp%lakedepth(c) / nlevlak;
             col_pp%z_lake(c,1)  = col_pp%dz_lake(c,1) / 2._r8;
             do lev=2,nlevlak
                col_pp%z_lake(c,lev) = col_pp%z_lake(c,lev-1) + (col_pp%dz_lake(c,lev-1)+col_pp%dz_lake(c,lev))/2._r8
             end do

          else

             write(iulog,*)'Bad lake depth: lakedepth: ', col_pp%lakedepth(c)
             call shr_sys_abort(errmsg(__FILE__, __LINE__))

          end if

          col_pp%z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
          col_pp%zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
          col_pp%dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
       end if
    end do

    !-----------------------------------------------
    ! Set cold-start values for snow levels, snow layers and snow interfaces 
    !-----------------------------------------------

    associate(snl => col_pp%snl) ! Output: [integer (:)    ]  number of snow layers   

      do c = bounds%begc,bounds%endc
         l = col_pp%landunit(c)

         col_pp%dz(c,-nlevsno+1: 0) = spval
         col_pp%z (c,-nlevsno+1: 0) = spval
         col_pp%zi(c,-nlevsno  :-1) = spval

         if (.not. lun_pp%lakpoi(l)) then
            if (snow_depth(c) < 0.01_r8) then
               snl(c)             = 0
               col_pp%dz(c,-nlevsno+1:0) = 0._r8
               col_pp%z (c,-nlevsno+1:0) = 0._r8
               col_pp%zi(c,-nlevsno+0:0) = 0._r8
            else
               if ((snow_depth(c) >= 0.01_r8) .and. (snow_depth(c) <= 0.03_r8)) then
                  snl(c)  = -1
                  col_pp%dz(c,0) = snow_depth(c)
               else if ((snow_depth(c) > 0.03_r8) .and. (snow_depth(c) <= 0.04_r8)) then
                  snl(c)   = -2
                  col_pp%dz(c,-1) = snow_depth(c)/2._r8
                  col_pp%dz(c, 0) = col_pp%dz(c,-1)
               else if ((snow_depth(c) > 0.04_r8) .and. (snow_depth(c) <= 0.07_r8)) then
                  snl(c)   = -2
                  col_pp%dz(c,-1) = 0.02_r8
                  col_pp%dz(c, 0) = snow_depth(c) - col_pp%dz(c,-1)
               else if ((snow_depth(c) > 0.07_r8) .and. (snow_depth(c) <= 0.12_r8)) then
                  snl(c)   = -3
                  col_pp%dz(c,-2) = 0.02_r8
                  col_pp%dz(c,-1) = (snow_depth(c) - 0.02_r8)/2._r8
                  col_pp%dz(c, 0) = col_pp%dz(c,-1)
               else if ((snow_depth(c) > 0.12_r8) .and. (snow_depth(c) <= 0.18_r8)) then
                  snl(c)   = -3
                  col_pp%dz(c,-2) = 0.02_r8
                  col_pp%dz(c,-1) = 0.05_r8
                  col_pp%dz(c, 0) = snow_depth(c) - col_pp%dz(c,-2) - col_pp%dz(c,-1)
               else if ((snow_depth(c) > 0.18_r8) .and. (snow_depth(c) <= 0.29_r8)) then
                  snl(c)   = -4
                  col_pp%dz(c,-3) = 0.02_r8
                  col_pp%dz(c,-2) = 0.05_r8
                  col_pp%dz(c,-1) = (snow_depth(c) - col_pp%dz(c,-3) - col_pp%dz(c,-2))/2._r8
                  col_pp%dz(c, 0) = col_pp%dz(c,-1)
               else if ((snow_depth(c) > 0.29_r8) .and. (snow_depth(c) <= 0.41_r8)) then
                  snl(c)   = -4
                  col_pp%dz(c,-3) = 0.02_r8
                  col_pp%dz(c,-2) = 0.05_r8
                  col_pp%dz(c,-1) = 0.11_r8
                  col_pp%dz(c, 0) = snow_depth(c) - col_pp%dz(c,-3) - col_pp%dz(c,-2) - col_pp%dz(c,-1)
               else if ((snow_depth(c) > 0.41_r8) .and. (snow_depth(c) <= 0.64_r8)) then
                  snl(c)   = -5
                  col_pp%dz(c,-4) = 0.02_r8
                  col_pp%dz(c,-3) = 0.05_r8
                  col_pp%dz(c,-2) = 0.11_r8
                  col_pp%dz(c,-1) = (snow_depth(c) - col_pp%dz(c,-4) - col_pp%dz(c,-3) - col_pp%dz(c,-2))/2._r8
                  col_pp%dz(c, 0) = col_pp%dz(c,-1)
               else if (snow_depth(c) > 0.64_r8) then
                  snl(c)   = -5
                  col_pp%dz(c,-4) = 0.02_r8
                  col_pp%dz(c,-3) = 0.05_r8
                  col_pp%dz(c,-2) = 0.11_r8
                  col_pp%dz(c,-1) = 0.23_r8
                  col_pp%dz(c, 0) = snow_depth(c)-col_pp%dz(c,-4)-col_pp%dz(c,-3)-col_pp%dz(c,-2)-col_pp%dz(c,-1)
               endif
            end if
            do j = 0, snl(c)+1, -1
               col_pp%z(c,j)    = col_pp%zi(c,j) - 0.5_r8*col_pp%dz(c,j)
               col_pp%zi(c,j-1) = col_pp%zi(c,j) - col_pp%dz(c,j)
            end do
         else !lake
            snl(c)             = 0
            col_pp%dz(c,-nlevsno+1:0) = 0._r8
            col_pp%z (c,-nlevsno+1:0) = 0._r8
            col_pp%zi(c,-nlevsno+0:0) = 0._r8
         end if
      end do

      !-----------------------------------------------
      ! Read in topographic index and slope
      !-----------------------------------------------

      allocate(tslope(bounds%begg:bounds%endg))
      call ncd_io(ncid=ncid, varname='SLOPE', flag='read', data=tslope, dim1name=grlnd, readvar=readvar)
      if (.not. readvar) then
         call shr_sys_abort(' ERROR: TOPOGRAPHIC SLOPE NOT on surfdata file'//&
              errMsg(__FILE__, __LINE__)) 
      end if
      do c = begc,endc
         g = col_pp%gridcell(c)
         ! check for near zero slopes, set minimum value
         col_pp%topo_slope(c) = max(tslope(g), 0.2_r8)
      end do
      deallocate(tslope)

      allocate(std(bounds%begg:bounds%endg))
      call ncd_io(ncid=ncid, varname='STD_ELEV', flag='read', data=std, dim1name=grlnd, readvar=readvar)
      if (.not. readvar) then
         call shr_sys_abort(' ERROR: TOPOGRAPHIC STDdev (STD_ELEV) NOT on surfdata file'//&
              errMsg(__FILE__, __LINE__)) 
      end if
      do c = begc,endc
         g = col_pp%gridcell(c)
         ! Topographic variables
         col_pp%topo_std(c) = std(g)
      end do
      deallocate(std)

      if (use_erosion) then
         allocate(hslp_p10(bounds%begg:bounds%endg,nlevslp))
         call ncd_io(ncid=ncid, varname='SLP_P10', flag='read', data=hslp_p10, dim1name=grlnd, readvar=readvar)
         if (.not. readvar) then
            call shr_sys_abort(' ERROR: hillslope slope percentiles NOT on surfdata file'//&
                 errMsg(__FILE__, __LINE__))
         end if
         do c = begc,endc
            g = col_pp%gridcell(c)
            col_pp%hslp_p10(c,:) = hslp_p10(g,:)
         end do
         deallocate(hslp_p10)
      else
         do c = begc,endc
            g = col_pp%gridcell(c)
            col_pp%hslp_p10(c,:) = 0._r8
         end do
      end if

      !-----------------------------------------------
      ! Read in depth to bedrock
      !-----------------------------------------------

      if (use_var_soil_thick) then
         allocate(dtb(bounds%begg:bounds%endg))
         call ncd_io(ncid=ncid, varname='aveDTB', flag='read', data=dtb, dim1name=grlnd, readvar=readvar)
         if (.not. readvar) then
            write(iulog,*) 'aveDTB not in surfdata: reverting to default 10 layers.'
            do c = begc,endc
               col_pp%nlevbed(c) = nlevsoi
	       col_pp%zibed(c) = zisoi(nlevsoi)
	    end do
         else
	    do c = begc,endc
               g = col_pp%gridcell(c)
               l = col_pp%landunit(c)
               if (lun_pp%urbpoi(l) .and. col_pp%itype(c) /= icol_road_imperv .and. col_pp%itype(c) /= icol_road_perv) then
               	  col_pp%nlevbed(c) = nlevurb
               else if (lun_pp%itype(l) == istdlak) then
               	  col_pp%nlevbed(c) = nlevlak
               else if (lun_pp%itype(l) == istice_mec) then
               	  col_pp%nlevbed(c) = 5
               else
                  ! check for near zero DTBs, set minimum value
	          beddep = max(dtb(g), 0.2_r8)
	          j = 0
	          zimid = 0._r8
                  do while (zimid < beddep .and. j < nlevgrnd)
	             zimid = 0.5_r8*(zisoi(j)+zisoi(j+1))
	             if (beddep > zimid) then
	                nlevbed = j + 1
	             else
	                nlevbed = j
                     end if
	             j = j + 1
                  enddo
	          nlevbed = max(nlevbed, 5)
	          nlevbed = min(nlevbed, nlevgrnd)
                  col_pp%nlevbed(c) = nlevbed
	          col_pp%zibed(c) = zisoi(nlevbed)
               end if
            end do
	 end if
         deallocate(dtb)
      else
         do c = begc,endc
            col_pp%nlevbed(c) = nlevsoi
	    col_pp%zibed(c) = zisoi(nlevsoi)
	 end do
      end if

      !-----------------------------------------------
      ! SCA shape function defined
      !-----------------------------------------------

      do c = begc,endc
         l = col_pp%landunit(c)

         if (lun_pp%itype(l)==istice_mec) then
            ! ice_mec columns already account for subgrid topographic variability through
            ! their use of multiple elevation classes; thus, to avoid double-accounting for
            ! topographic variability in these columns, we ignore topo_std and use a value
            ! of n_melt that assumes little topographic variability within the column
            col_pp%n_melt(c) = 10._r8
         else
            col_pp%n_melt(c) = 200.0/max(10.0_r8, col_pp%topo_std(c))
         end if

         ! microtopographic parameter, units are meters (try smooth function of slope)

         slopebeta = 3._r8
         slopemax = 0.4_r8
         slope0 = slopemax**(-1._r8/slopebeta)
         col_pp%micro_sigma(c) = (col_pp%topo_slope(c) + slope0)**(-slopebeta)
      end do

    end associate

    call ncd_pio_closefile(ncid)

  end subroutine initVertical

end module initVerticalMod
