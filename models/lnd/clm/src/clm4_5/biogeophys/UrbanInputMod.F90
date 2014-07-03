module UrbanInputMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: UrbanInputMod
! 
! !DESCRIPTION: 
! Read in input urban data - fill in data structure urbinp
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun  
  use shr_sys_mod , only : shr_sys_flush 
!
! !PUBLIC TYPES:
  implicit none
  save

  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanInput         ! Read in urban input data

  type urbinp_t
     real(r8), pointer :: canyon_hwr(:,:)  
     real(r8), pointer :: wtlunit_roof(:,:)  
     real(r8), pointer :: wtroad_perv(:,:)  
     real(r8), pointer :: em_roof(:,:)   
     real(r8), pointer :: em_improad(:,:)  
     real(r8), pointer :: em_perroad(:,:)  
     real(r8), pointer :: em_wall(:,:)  
     real(r8), pointer :: alb_roof_dir(:,:,:)  
     real(r8), pointer :: alb_roof_dif(:,:,:)  
     real(r8), pointer :: alb_improad_dir(:,:,:)  
     real(r8), pointer :: alb_improad_dif(:,:,:)  
     real(r8), pointer :: alb_perroad_dir(:,:,:)  
     real(r8), pointer :: alb_perroad_dif(:,:,:)  
     real(r8), pointer :: alb_wall_dir(:,:,:)  
     real(r8), pointer :: alb_wall_dif(:,:,:)  
     real(r8), pointer :: ht_roof(:,:)
     real(r8), pointer :: wind_hgt_canyon(:,:)
     real(r8), pointer :: tk_wall(:,:,:)
     real(r8), pointer :: tk_roof(:,:,:)
     real(r8), pointer :: tk_improad(:,:,:)
     real(r8), pointer :: cv_wall(:,:,:)
     real(r8), pointer :: cv_roof(:,:,:)
     real(r8), pointer :: cv_improad(:,:,:)
     real(r8), pointer :: thick_wall(:,:)
     real(r8), pointer :: thick_roof(:,:)
     integer,  pointer :: nlev_improad(:,:)
     real(r8), pointer :: t_building_min(:,:)
     real(r8), pointer :: t_building_max(:,:)
  end type urbinp_t
  public urbinp_t

  type (urbinp_t)   , public :: urbinp        ! urban input derived type
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanInput
!
! !INTERFACE:
  subroutine UrbanInput(mode)
!
! !DESCRIPTION: 
! Allocate memory and read in urban input data
!
! !USES:
    use clm_varpar, only : numrad, nlevurb, numurbl
    use clm_varctl, only : iulog, fsurdat, single_column
    use fileutils , only : getavu, relavu, getfil, opnfil
    use spmdMod   , only : masterproc
    use clmtype   , only : grlnd
    use decompMod , only : get_proc_bounds
    use domainMod , only : ldomain
    use ncdio_pio 
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: mode
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein July 2004
! Revised by Keith Oleson for netcdf input Jan 2008
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn      ! local file name
    type(file_desc_t)  :: ncid       ! netcdf id
    integer :: dimid                 ! netCDF id
    integer :: begg,endg             ! start/stop gridcells
    integer :: nw,n,k,i,j,ni,nj,ns   ! indices
    integer :: nlevurb_i             ! input grid: number of urban vertical levels
    integer :: numrad_i              ! input grid: number of solar bands (VIS/NIR)
    integer :: numurbl_i             ! input grid: number of urban landunits
    integer :: ier,ret               ! error status
    logical :: isgrid2d              ! true => file is 2d 
    logical :: readvar               ! true => variable is on dataset
    logical :: has_numurbl           ! true => numurbl dimension is on dataset
    character(len=32) :: subname = 'UrbanInput' ! subroutine name
!-----------------------------------------------------------------------

    if ( nlevurb == 0 ) return

    call get_proc_bounds(begg,endg)

    if (mode == 'initialize') then

       ! Read urban data
       
       if (masterproc) then
          write(iulog,*)' Reading in urban input data from fsurdat file ...'
       end if
       
       call getfil (fsurdat, locfn, 0)
       call ncd_pio_openfile (ncid, locfn, 0)

       if (masterproc) then
          write(iulog,*) subname,trim(fsurdat)
       end if

       ! Check whether this file has new-format urban data
       call ncd_inqdid(ncid, 'numurbl', dimid, dimexist=has_numurbl)

       ! If file doesn't have numurbl, then it is old-format urban;
       ! in this case, set nlevurb to zero
       if (.not. has_numurbl) then
         nlevurb = 0
         write(iulog,*)'PCT_URBAN is not multi-density, nlevurb set to 0'
       end if

       if ( nlevurb == 0 ) return

       ! Allocate dynamic memory
       allocate(urbinp%canyon_hwr(begg:endg,numurbl), &  
                urbinp%wtlunit_roof(begg:endg,numurbl), &  
                urbinp%wtroad_perv(begg:endg,numurbl), &
                urbinp%em_roof(begg:endg,numurbl), &     
                urbinp%em_improad(begg:endg,numurbl), &    
                urbinp%em_perroad(begg:endg,numurbl), &    
                urbinp%em_wall(begg:endg,numurbl), &    
                urbinp%alb_roof_dir(begg:endg,numurbl,numrad), &    
                urbinp%alb_roof_dif(begg:endg,numurbl,numrad), &    
                urbinp%alb_improad_dir(begg:endg,numurbl,numrad), &    
                urbinp%alb_perroad_dir(begg:endg,numurbl,numrad), &    
                urbinp%alb_improad_dif(begg:endg,numurbl,numrad), &    
                urbinp%alb_perroad_dif(begg:endg,numurbl,numrad), &    
                urbinp%alb_wall_dir(begg:endg,numurbl,numrad), &    
                urbinp%alb_wall_dif(begg:endg,numurbl,numrad), &
                urbinp%ht_roof(begg:endg,numurbl), &
                urbinp%wind_hgt_canyon(begg:endg,numurbl), &
                urbinp%tk_wall(begg:endg,numurbl,nlevurb), &
                urbinp%tk_roof(begg:endg,numurbl,nlevurb), &
                urbinp%tk_improad(begg:endg,numurbl,nlevurb), &
                urbinp%cv_wall(begg:endg,numurbl,nlevurb), &
                urbinp%cv_roof(begg:endg,numurbl,nlevurb), &
                urbinp%cv_improad(begg:endg,numurbl,nlevurb), &
                urbinp%thick_wall(begg:endg,numurbl), &
                urbinp%thick_roof(begg:endg,numurbl), &
                urbinp%nlev_improad(begg:endg,numurbl), &
                urbinp%t_building_min(begg:endg,numurbl), &
                urbinp%t_building_max(begg:endg,numurbl), &
                stat=ier)
       if (ier /= 0) then
          write(iulog,*)'initUrbanInput: allocation error '; call endrun()
       endif

       call ncd_inqfdims (ncid, isgrid2d, ni, nj, ns)
       if (ldomain%ns /= ns .or. ldomain%ni /= ni .or. ldomain%nj /= nj) then
          write(iulog,*)trim(subname), 'ldomain and input file do not match dims '
          write(iulog,*)trim(subname), 'ldomain%ni,ni,= ',ldomain%ni,ni
          write(iulog,*)trim(subname), 'ldomain%nj,nj,= ',ldomain%nj,nj
          write(iulog,*)trim(subname), 'ldomain%ns,ns,= ',ldomain%ns,ns
          call endrun()
       end if

       call ncd_inqdid(ncid, 'nlevurb', dimid)
       call ncd_inqdlen(ncid, dimid, nlevurb_i)
       if (nlevurb_i /= nlevurb) then
          write(iulog,*)trim(subname)// ': parameter nlevurb= ',nlevurb, &
               'does not equal input dataset nlevurb= ',nlevurb_i
          call endrun
       endif

       call ncd_inqdid(ncid, 'numrad', dimid)
       call ncd_inqdlen(ncid, dimid, numrad_i)
       if (numrad_i /= numrad) then
          write(iulog,*)trim(subname)// ': parameter numrad= ',numrad, &
               'does not equal input dataset numrad= ',numrad_i
          call endrun
       endif
       call ncd_inqdid(ncid, 'numurbl', dimid)
       call ncd_inqdlen(ncid, dimid, numurbl_i)
       if (numurbl_i /= numurbl) then
          write(iulog,*)trim(subname)// ': parameter numurbl= ',numurbl, &
               'does not equal input dataset numurbl= ',numurbl_i
          call endrun
       endif
       call ncd_io(ncid=ncid, varname='CANYON_HWR', flag='read', data=urbinp%canyon_hwr,&
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: CANYON_HWR NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='WTLUNIT_ROOF', flag='read', data=urbinp%wtlunit_roof, &
            dim1name=grlnd,  readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: WTLUNIT_ROOF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='WTROAD_PERV', flag='read', data=urbinp%wtroad_perv, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: WTROAD_PERV NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='EM_ROOF', flag='read', data=urbinp%em_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: EM_ROOF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='EM_IMPROAD', flag='read', data=urbinp%em_improad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: EM_IMPROAD NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='EM_PERROAD', flag='read', data=urbinp%em_perroad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: EM_PERROAD NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='EM_WALL', flag='read', data=urbinp%em_wall, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: EM_WALL NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='HT_ROOF', flag='read', data=urbinp%ht_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: HT_ROOF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='WIND_HGT_CANYON', flag='read', data=urbinp%wind_hgt_canyon, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: WIND_HGT_CANYON NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='THICK_WALL', flag='read', data=urbinp%thick_wall, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: THICK_WALL NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='THICK_ROOF', flag='read', data=urbinp%thick_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: THICK_ROOF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='NLEV_IMPROAD', flag='read', data=urbinp%nlev_improad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: NLEV_IMPROAD NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='T_BUILDING_MIN', flag='read', data=urbinp%t_building_min, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: T_BUILDING_MIN NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='T_BUILDING_MAX', flag='read', data=urbinp%t_building_max, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: T_BUILDING_MAX NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='ALB_IMPROAD_DIR', flag='read', data=urbinp%alb_improad_dir, &
            dim1name=grlnd, readvar=readvar)
       if (.not.readvar) call endrun( trim(subname)//' ERROR: ALB_IMPROAD_DIR NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='ALB_IMPROAD_DIF', flag='read', data=urbinp%alb_improad_dif, &
            dim1name=grlnd, readvar=readvar)
       if (.not.readvar) call endrun( trim(subname)//' ERROR: ALB_IMPROAD_DIF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='ALB_PERROAD_DIR', flag='read',data=urbinp%alb_perroad_dir, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: ALB_PERROAD_DIR NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='ALB_PERROAD_DIF', flag='read',data=urbinp%alb_perroad_dif, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: ALB_PERROAD_DIF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='ALB_ROOF_DIR', flag='read', data=urbinp%alb_roof_dir,  &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: ALB_ROOF_DIR NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='ALB_ROOF_DIF', flag='read', data=urbinp%alb_roof_dif,  &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: ALB_ROOF_DIF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='ALB_WALL_DIR', flag='read', data=urbinp%alb_wall_dir, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: ALB_WALL_DIR NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='ALB_WALL_DIF', flag='read', data=urbinp%alb_wall_dif, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: ALB_WALL_DIF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='TK_IMPROAD', flag='read', data=urbinp%tk_improad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TK_IMPROAD NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='TK_ROOF', flag='read', data=urbinp%tk_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TK_ROOF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='TK_WALL', flag='read', data=urbinp%tk_wall, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TK_WALL NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='CV_IMPROAD', flag='read', data=urbinp%cv_improad, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: CV_IMPROAD NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='CV_ROOF', flag='read', data=urbinp%cv_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: CV_ROOF NOT on fsurdat file' )

       call ncd_io(ncid=ncid, varname='CV_WALL', flag='read', data=urbinp%cv_wall, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: CV_WALL NOT on fsurdat file' )

       call ncd_pio_closefile(ncid)
       if (masterproc) then
          write(iulog,*)' Sucessfully read urban input data' 
          write(iulog,*)
       end if

    else if (mode == 'finalize') then

       if ( nlevurb == 0 ) return

       deallocate(urbinp%canyon_hwr, &
                  urbinp%wtlunit_roof, &
                  urbinp%wtroad_perv, &
                  urbinp%em_roof, &
                  urbinp%em_improad, &
                  urbinp%em_perroad, &
                  urbinp%em_wall, &
                  urbinp%alb_roof_dir, &
                  urbinp%alb_roof_dif, &
                  urbinp%alb_improad_dir, &
                  urbinp%alb_perroad_dir, &
                  urbinp%alb_improad_dif, &
                  urbinp%alb_perroad_dif, &
                  urbinp%alb_wall_dir, &
                  urbinp%alb_wall_dif, &
                  urbinp%ht_roof, &
                  urbinp%wind_hgt_canyon, &
                  urbinp%tk_wall, &
                  urbinp%tk_roof, &
                  urbinp%tk_improad, &
                  urbinp%cv_wall, &
                  urbinp%cv_roof, &
                  urbinp%cv_improad, &
                  urbinp%thick_wall, &
                  urbinp%thick_roof, &
                  urbinp%nlev_improad, &
                  urbinp%t_building_min, &
                  urbinp%t_building_max, &
                  stat=ier)
       if (ier /= 0) then
          write(iulog,*)'initUrbanInput: deallocation error '; call endrun()
       endif

    else
       write(iulog,*)'initUrbanInput error: mode ',trim(mode),' not supported '
       call endrun()
    end if

  end subroutine UrbanInput

end module UrbanInputMod

