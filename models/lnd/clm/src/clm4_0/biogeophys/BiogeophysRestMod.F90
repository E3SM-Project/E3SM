module BiogeophysRestMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: BiogeophysRestMod
!
! !DESCRIPTION:
! Reads from or biogeophysics restart/initial data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use spmdMod     , only : masterproc
!
! !PUBLIC TYPES:
  implicit none

  private
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: BiogeophysRest
!
! !REVISION HISTORY:
! 2005-06-12: Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

  private :: weights_exactly_the_same
  private :: weights_within_roundoff_different
  private :: weights_tooDifferent

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BiogeophysRest
!
! !INTERFACE:
  subroutine BiogeophysRest( ncid, flag )
!
! !DESCRIPTION:
! Read/Write biogeophysics information to/from restart file.
!
! !USES:
    use ncdio_pio
    use clmtype
    use decompMod       , only : get_proc_bounds
    use clm_varpar      , only : nlevgrnd, nlevsno, nlevlak, nlevurb
    use clm_varcon      , only : istcrop
    use clm_varcon      , only : denice, denh2o, istdlak, istslak, isturb, &
                                 istsoil, pondmx, watmin, spval
    use clm_varctl      , only : allocate_all_vegpfts, nsrest, fpftdyn,    &
                                 iulog, nsrContinue, nsrStartup, nsrBranch, &
                                 use_cndv, use_cn, use_snicar_frc
    use initSurfAlbMod  , only : do_initsurfalb
    use clm_time_manager, only : is_first_step
    use SNICARMod       , only : snw_rds_min
    use shr_infnan_mod  , only : shr_infnan_isnan
    use clm_time_manager, only : is_restart
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in)  :: flag  ! 'read' or 'write'
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 12/11/2003, Peter Thornton: Added cps%coszen, pps%gdir, and pps%omega
!   for new sunlit/shaded canopy algorithm (in SUNSHA ifdef block)
! 4/25/2005, Peter Thornton: Removed the SUNSHA ifdefs, since this is now the
!   default code behavior.
! 6/12/2005, Moved to netcdf format and renamed file
!
!
! !LOCAL VARIABLES:
!EOP
!
! local pointers to implicit in arguments
!
    real(r8) :: maxwatsat                 !maximum porosity    
    real(r8) :: excess                    !excess volumetric soil water
    real(r8) :: totwat                    !total soil water (mm)
    real(r8) :: maxdiff                   !maximum difference in PFT weights
    real(r8), pointer :: wtgcell(:)       ! Grid cell weights for PFT
    real(r8), pointer :: wtlunit(:)       ! Land-unit weights for PFT
    real(r8), pointer :: wtcol(:)         ! Column weights for PFT
    integer :: p,c,l,g,j    ! indices
    integer :: nlevs        ! number of layers
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    logical :: readvar      ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    integer , pointer :: clandunit(:)     ! landunit of corresponding column
    integer , pointer :: ltype(:)         ! landunit type
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
    real(r8), pointer   :: temp2d(:,:)    ! temporary for zisno
    real(r8), parameter :: adiff = 5.e-04_r8   ! tolerance of acceptible difference
    character(len=7)  :: filetypes(0:3)
    character(len=32) :: fileusing
    character(len=*), parameter :: sub="BiogeophysRest"
!-----------------------------------------------------------------------
    filetypes(:)           = "missing"
    filetypes(nsrStartup)  = "finidat"
    filetypes(nsrContinue) = "restart"
    filetypes(nsrBranch)   = "nrevsn"

    ! Set pointers into derived type

    gptr       => grc
    lptr       => lun
    cptr       => col
    pptr       => pft
    ltype      => lptr%itype
    clandunit  => cptr%landunit
    clandunit  => cptr%landunit

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    !
    ! Read in weights if allocating all vegetation types
    !

    if (allocate_all_vegpfts) then

       ! pft weight wrt gridcell 

       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='PFT_WTGCELL', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='pft weight relative to corresponding gridcell', units='')
       else if (flag == 'read' .or. flag == 'write') then
          ! Copy weights calculated from fsurdat/fpftdyn to temp array for comparision
          ! Don't read directly into temp array -- so that answers are identical with clm3.6.58. EBK 1/9/2010
          if (flag == 'read' )then
             allocate( wtgcell(begp:endp) )
             wtgcell(:) = pptr%wtgcell(:)
          end if
          call ncd_io(varname='PFT_WTGCELL', data=pptr%wtgcell, &
               dim1name=namep, &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun()
          end if
       end if

       ! pft weight wrt landunit

       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='PFT_WTLUNIT', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='pft weight relative to corresponding landunit', units='')
       else if (flag == 'read' .or. flag == 'write') then
          ! Copy weights calculated from fsurdat/fpftdyn to temp array for comparision
          ! Don't read directly into temp array -- so that answers are identical with clm3.6.58. EBK 1/9/2010
          if (flag == 'read' )then
             allocate( wtlunit(begp:endp) )
             wtlunit(:) = pptr%wtlunit(:)
          end if
          call ncd_io(varname='PFT_WTLUNIT', data=pptr%wtlunit, &
               dim1name=namep, &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun()
          end if
       end if

       ! pft weight wrt column

       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='PFT_WTCOL', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='pft weight relative to corresponding column', units='')
       else if (flag == 'read' .or. flag == 'write') then
          ! Copy weights calculated from fsurdat/fpftdyn to temp array for comparision
          ! Don't read directly into temp array -- so that answers are identical with clm3.6.58. EBK 1/9/2010
          if (flag == 'read' )then
             allocate( wtcol(begp:endp)   )
             wtcol(:) = pptr%wtcol(:)
          end if
          call ncd_io(varname='PFT_WTCOL', data=pptr%wtcol, &
               dim1name=namep, &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun()
          end if
       end if

       if (flag == 'read' )then

          if ( fpftdyn /= ' ' )then
             fileusing = "fsurdat/fpftdyn"
          else
             fileusing = "fsurdat"
          end if
          !
          ! Note: Do not compare weights if restart or if dynamic-pft branch
          !
          if ( nsrest == nsrContinue .or. fpftdyn /= ' ' )then
             ! Do NOT do any testing for restart or a pftdyn case
          !
          ! Otherwise test and make sure weights agree to reasonable tolerence
          !
          else if ( .not.weights_exactly_the_same( pptr, wtgcell, wtlunit, wtcol ) )then
             if (.not. use_cndv) then
                if (      weights_within_roundoff_different( pptr, wtgcell, wtlunit, wtcol ) )then
                   write(iulog,*) sub//"::NOTE, PFT weights from ", filetypes(nsrest),      &
                        " file and ", trim(fileusing), " file(s) are different to roundoff -- using ", &
                        trim(fileusing), " values."
                else if ( weights_tooDifferent( begp, endp, pptr, wtgcell, adiff, maxdiff ) )then
                   write(iulog,*) "ERROR:: PFT weights are SIGNIFICANTLY different from the input ", &
                        filetypes(nsrest), " file and ", trim(fileusing), " file(s)."
                   write(iulog,*) "ERROR:: maximum difference is ", maxdiff, " max allowed = ", adiff
                   write(iulog,*) "ERROR:: Run interpinic on your initial condition file to interpolate to the new surface dataset"
                   call endrun( sub//"::ERROR:: Weights between initial condition file and surface dataset are too different" )
                else
                   write(iulog,*) sub//"::NOTE, PFT weights from ", filetypes(nsrest),      &
                        " file and ", trim(fileusing), " file(s) are different to < ", &
                        adiff, " -- using ", trim(fileusing), " values."
                end if
                write(iulog,*) sub//"::WARNING, weights different between ", filetypes(nsrest), &
                     " file and ", trim(fileusing), " file(s), but close enough -- using ",    &
                     trim(fileusing), " values."
                ! Copy weights from fsurdat file back in -- they are only off by roundoff to 1% or so...
                pptr%wtgcell(:) = wtgcell(:)
                pptr%wtlunit(:) = wtlunit(:)
                pptr%wtcol(:)   = wtcol(:)
             end if
          end if
 
          deallocate( wtgcell )
          deallocate( wtlunit )
          deallocate( wtcol   )

       end if

    end if

    ! Note - for the snow interfaces, are only examing the snow interfaces
    ! above zi=0 which is why zisno and zsno have the same level dimension below
    ! (Note - for zisno, zi(0) is set to 0 in routine iniTimeConst)
    
    ! pft energy flux - eflx_lwrad_out

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='EFLX_LWRAD_OUT', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='emitted infrared (longwave) radiation', units='watt/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='EFLX_LWRAD_OUT', data=pef%eflx_lwrad_out, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - snow levels

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='SNLSNO', xtype=ncd_int,  &
            dim1name='column', &
            long_name='number of snow layers', units='unitless')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='SNLSNO', data=cps%snl, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - snowdp

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='SNOWDP', xtype=ncd_double, &
            dim1name='column', &
            long_name='snow depth', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='SNOWDP', data=cps%snowdp, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - wa

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='WA', xtype=ncd_double,  &
            dim1name='column', &
            long_name='water in the unconfined aquifer', units='mm')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='WA', data=cws%wa, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - wt

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='WT', xtype=ncd_double,  &
            dim1name='column', &
            long_name='total water storage', units='mm')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='WT', data=cws%wt, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - zwt

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ZWT', xtype=ncd_double,  &
            dim1name='column', &
            long_name='water table depth', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='ZWT', data=cws%zwt, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - frac_sno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frac_sno', xtype=ncd_double,  &
            dim1name='column',&
            long_name='fraction of ground covered by snow (0 to 1)',units='unitless')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frac_sno', data=cps%frac_sno, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - dzsno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='DZSNO', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer thickness', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       allocate(temp2d(begc:endc,-nlevsno+1:0))
       if (flag == 'write') then 
          temp2d(begc:endc,-nlevsno+1:0) = cps%dz(begc:endc,-nlevsno+1:0)
       end if
       call ncd_io(varname='DZSNO', data=temp2d, &
            dim1name=namec, switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
       if (flag == 'read') then
          cps%dz(begc:endc,-nlevsno+1:0) = temp2d(begc:endc,-nlevsno+1:0) 
       end if
       deallocate(temp2d)
    end if

    ! column type physical state variable - zsno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ZSNO', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer depth', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       allocate(temp2d(begc:endc,-nlevsno+1:0))
       if (flag == 'write') then 
          temp2d(begc:endc,-nlevsno+1:0) = cps%z(begc:endc,-nlevsno+1:0)
       end if
       call ncd_io(varname='ZSNO', data=temp2d, &
            dim1name=namec, switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
       if (flag == 'read') then
          cps%z(begc:endc,-nlevsno+1:0) = temp2d(begc:endc,-nlevsno+1:0) 
       end if
       deallocate(temp2d)
    end if

    ! column type physical state variable - zisno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ZISNO', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow interface depth', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       allocate(temp2d(begc:endc,-nlevsno:-1))
       if (flag == 'write') then 
          temp2d(begc:endc,-nlevsno:-1) = cps%zi(begc:endc,-nlevsno:-1)
       end if
       call ncd_io(varname='ZISNO', data=temp2d, &
            dim1name=namec, switchdim=.true., &
            lowerb2=-nlevsno, upperb2=-1, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
       if (flag == 'read') then 
          cps%zi(begc:endc,-nlevsno:-1) = temp2d(begc:endc,-nlevsno:-1)
       end if
       deallocate(temp2d)	
    end if

    ! column type physical state variable - coszen

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='coszen', xtype=ncd_double,  &
            dim1name='column', &
            long_name='cosine of solar zenith angle', units='unitless')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='coszen', data=cps%coszen, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - gdir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gdir', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='leaf projection in solar direction (0 to 1)', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='gdir', data=pps%gdir, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - omega

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='omega', xtype=ncd_double,  &
            dim1name='pft', dim2name='numrad', switchdim=.true., &
            long_name='fraction of intercepted radiation that is scattered (0 to 1)', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='omega', data=pps%omega, &
            dim1name=namep, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_roof_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_roof_dir', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='direct solar absorbed by roof per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_roof_dir', data=lps%sabs_roof_dir, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_roof_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_roof_dif', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='diffuse solar absorbed by roof per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_roof_dif', data=lps%sabs_roof_dif, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_sunwall_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_sunwall_dir', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='direct solar absorbed by sunwall per unit wall area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_sunwall_dir', data=lps%sabs_sunwall_dir, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_sunwall_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_sunwall_dif', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='diffuse solar absorbed by sunwall per unit wall area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_sunwall_dif', data=lps%sabs_sunwall_dif, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_shadewall_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_shadewall_dir', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='direct solar absorbed by shadewall per unit wall area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_shadewall_dir', data=lps%sabs_shadewall_dir, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_shadewall_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_shadewall_dif', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='diffuse solar absorbed by shadewall per unit wall area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_shadewall_dif', data=lps%sabs_shadewall_dif, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_improad_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_improad_dir', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='direct solar absorbed by impervious road per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_improad_dir', data=lps%sabs_improad_dir, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_improad_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_improad_dif', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='diffuse solar absorbed by impervious road per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_improad_dif', data=lps%sabs_improad_dif, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_perroad_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_perroad_dir', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='direct solar absorbed by pervious road per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_perroad_dir', data=lps%sabs_perroad_dir, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_perroad_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_perroad_dif', xtype=ncd_double,  &
            dim1name='landunit', dim2name='numrad', switchdim=.true., &
            long_name='diffuse solar absorbed by pervious road per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sabs_perroad_dif', data=lps%sabs_perroad_dif, &
            dim1name=namel, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_sr

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_sr', xtype=ncd_double,  &
            dim1name='landunit', &
            long_name='view factor of sky for road',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='vf_sr', data=lps%vf_sr, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_wr

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_wr', xtype=ncd_double,  &
            dim1name='landunit', &
            long_name='view factor of one wall for road',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='vf_wr', data=lps%vf_wr, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_sw

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_sw', xtype=ncd_double,  &
            dim1name='landunit', &
            long_name='view factor of sky for one wall',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='vf_sw', data=lps%vf_sw, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_rw

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_rw', xtype=ncd_double,  &
            dim1name='landunit', &
            long_name='view factor of road for one wall',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='vf_rw', data=lps%vf_rw, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_ww

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_ww', xtype=ncd_double,  &
            dim1name='landunit', &
            long_name='view factor of opposing wall for one wall',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='vf_ww', data=lps%vf_ww, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - taf

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='taf', xtype=ncd_double,  &
            dim1name='landunit', &
            long_name='urban canopy air temperature',units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='taf', data=lps%taf, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - qaf

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='qaf', xtype=ncd_double,  &
            dim1name='landunit', &
            long_name='urban canopy specific humidity',units='kg/kg')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='qaf', data=lps%qaf, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - albd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albd', xtype=ncd_double,  &
            dim1name='pft', dim2name='numrad', switchdim=.true., &
            long_name='surface albedo (direct) (0 to 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='albd', data=pps%albd, &
            dim1name=namep, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             if (nsrest == nsrStartup) do_initsurfalb = .true.
          end if
       end if
    end if

    ! pft type physical state variable - albi

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albi', xtype=ncd_double,  &
            dim1name='pft', dim2name='numrad', switchdim=.true., &
            long_name='surface albedo (diffuse) (0 to 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='albi', data=pps%albi, &
            dim1name=namep, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             if (nsrest == nsrStartup) do_initsurfalb = .true.
          end if
       end if
    end if

    ! column type physical state variable - albgrd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albgrd', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo (direct) (0 to 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='albgrd', data=cps%albgrd, &
            dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - albgri

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albgri', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo (indirect) (0 to 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='albgri', data=cps%albgri, &
            dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (use_snicar_frc) then
       ! column type physical state variable - albgrd_bc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='albgrd_bc', xtype=ncd_double,  &
               dim1name='column', dim2name='numrad', switchdim=.true., &
               long_name='ground albedo without BC (direct) (0 to 1)',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='albgrd_bc', data=cps%albgrd_bc, &
               dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_bc in restart (or initial) file..."
             if (masterproc) write(iulog,*) "Initialize albgrd_bc to albgrd"
             do c=begc,endc
                cps%albgrd_bc(c,:) = cps%albgrd(c,:)
             enddo
          end if
       end if
       ! column type physical state variable - albgri_bc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='albgri_bc', xtype=ncd_double,  &
               dim1name='column', dim2name='numrad', switchdim=.true., &
               long_name='ground albedo without BC (diffuse) (0 to 1)',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='albgri_bc', data=cps%albgri_bc, &
               dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (masterproc) write(iulog,*) "SNICAR: can't find albgri_bc in restart (or initial) file..."
             if (masterproc) write(iulog,*) "Initialize albgri_bc to albgri"
             do c=begc,endc
                cps%albgri_bc(c,:) = cps%albgri(c,:)
             enddo
          end if
       end if
       ! column type physical state variable - albgrd_pur
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='albgrd_pur', xtype=ncd_double,  &
               dim1name='column', dim2name='numrad', switchdim=.true., &
               long_name='pure snow ground albedo (direct) (0 to 1)',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='albgrd_pur', data=cps%albgrd_pur, &
               dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_pur in restart (or initial) file..."
             if (masterproc) write(iulog,*) "Initialize albgrd_pur to albgrd"
             do c=begc,endc
                cps%albgrd_pur(c,:) = cps%albgrd(c,:)
             enddo
          end if
       end if
       ! column type physical state variable - albgri_pur
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='albgri_pur', xtype=ncd_double,  &
               dim1name='column', dim2name='numrad', switchdim=.true., &
               long_name='pure snow ground albedo (diffuse) (0 to 1)',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='albgri_pur', data=cps%albgri_pur, &
               dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (masterproc) write(iulog,*) "SNICAR: can't find albgri_pur in restart (or initial) file..."
             if (masterproc) write(iulog,*) "Initialize albgri_pur to albgri"
             do c=begc,endc
                cps%albgri_pur(c,:) = cps%albgri(c,:)
             enddo
          end if
       end if
       ! column type physical state variable - albgrd_oc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='albgrd_oc', xtype=ncd_double,  &
               dim1name='column', dim2name='numrad', switchdim=.true., &
               long_name='ground albedo without OC (direct) (0 to 1)',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='albgrd_oc', data=cps%albgrd_oc, &
               dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_oc in restart (or initial) file..."
             if (masterproc) write(iulog,*) "Initialize albgrd_oc to albgrd"
             do c=begc,endc
                cps%albgrd_oc(c,:) = cps%albgrd(c,:)
             enddo
          end if
       end if
       ! column type physical state variable - albgri_oc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='albgri_oc', xtype=ncd_double,  &
               dim1name='column', dim2name='numrad', switchdim=.true., &
               long_name='ground albedo without OC (diffuse) (0 to 1)',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='albgri_oc', data=cps%albgri_oc, &
               dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (masterproc) write(iulog,*) "SNICAR: can't find albgri_oc in restart (or initial) file..."
             if (masterproc) write(iulog,*) "Initialize albgri_oc to albgri"
             do c=begc,endc
                cps%albgri_oc(c,:) = cps%albgri(c,:)
             enddo
          end if
       end if
       ! column type physical state variable - albgrd_dst
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='albgrd_dst', xtype=ncd_double,  &
               dim1name='column', dim2name='numrad', switchdim=.true., &
               long_name='ground albedo without dust (direct) (0 to 1)',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='albgrd_dst', data=cps%albgrd_dst, &
               dim1name=namec, switchdim=.true.,  ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_dst in restart (or initial) file..."
             if (masterproc) write(iulog,*) "Initialize albgrd_dst to albgrd"
             do c=begc,endc
                cps%albgrd_dst(c,:) = cps%albgrd(c,:)
             enddo
          end if
       end if
       ! column type physical state variable - albgri_dst
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='albgri_dst', xtype=ncd_double,  &
               dim1name='column', dim2name='numrad', switchdim=.true., &
               long_name='ground albedo without dust (diffuse) (0 to 1)',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='albgri_dst', data=cps%albgri_dst, &
               dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (masterproc) write(iulog,*) "SNICAR: can't find albgri_dst in restart (or initial) file..."
             if (masterproc) write(iulog,*) "Initialize albgri_dst to albgri"
             do c=begc,endc
                cps%albgri_dst(c,:) = cps%albgri(c,:)
             enddo
          end if
       end if
    end if  !end of if use_snicar_frc

   ! column water state variable - h2osno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='H2OSNO', xtype=ncd_double,  &
            dim1name='column', &
            long_name='snow water', units='kg/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='H2OSNO', data=cws%h2osno, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - h2osoi_liq

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='H2OSOI_LIQ', xtype=ncd_double,  &
            dim1name='column', dim2name='levtot', switchdim=.true., &
            long_name='liquid water', units='kg/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='H2OSOI_LIQ', data=cws%h2osoi_liq, &
            dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - h2osoi_ice

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='H2OSOI_ICE', xtype=ncd_double,   &
            dim1name='column', dim2name='levtot', switchdim=.true., &
            long_name='ice lens', units='kg/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='H2OSOI_ICE', data=cws%h2osoi_ice, &
            dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column energy state variable - t_grnd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_GRND', xtype=ncd_double,  &
            dim1name='column', &
            long_name='ground temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_GRND', data=ces%t_grnd, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column urban energy state variable - eflx_urban_ac

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='URBAN_AC', xtype=ncd_double,  &
            dim1name='column', &
            long_name='urban air conditioning flux', units='watt/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='URBAN_AC', data=cef%eflx_urban_ac, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column urban energy state variable - eflx_urban_heat

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='URBAN_HEAT', xtype=ncd_double,  &
            dim1name='column', &
            long_name='urban heating flux', units='watt/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='URBAN_HEAT', data=cef%eflx_urban_heat, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! pft energy state variable - t_ref2m_min

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='daily minimum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MIN', data=pes%t_ref2m_min, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='daily maximum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MAX', data=pes%t_ref2m_max, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_min_inst

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_INST', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='instantaneous daily min of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MIN_INST', data=pes%t_ref2m_min_inst, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_inst

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_INST', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='instantaneous daily max of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MAX_INST', data=pes%t_ref2m_max_inst, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

    ! pft energy state variable - t_ref2m_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname="T_REF2M_U", xtype=ncd_double,  &
            dim1name='pft', &
            long_name='Urban 2m height surface air temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname="T_REF2M_U", data=pes%t_ref2m_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column energy state variable - t_grnd_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_GRND_U', xtype=ncd_double,  &
            dim1name='column', &
            long_name='urban ground temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_GRND_U', data=ces%t_grnd_u, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! pft energy state variable - t_ref2m_min_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_U', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='urban daily minimum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MIN_U', data=pes%t_ref2m_min_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_U', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='urban daily maximum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MAX_U', data=pes%t_ref2m_max_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_min_inst_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_INST_U', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='urban instantaneous daily min of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MIN_INST_U', data=pes%t_ref2m_min_inst_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_inst_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_INST_U', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='urban instantaneous daily max of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MAX_INST_U', data=pes%t_ref2m_max_inst_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

    ! pft energy state variable - t_ref2m_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname="T_REF2M_R", xtype=ncd_double,  &
            dim1name='pft', &
            long_name='Rural 2m height surface air temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname="T_REF2M_R", data=pes%t_ref2m_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column energy state variable - t_grnd_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_GRND_R', xtype=ncd_double,  &
            dim1name='column', &
            long_name='rural ground temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_GRND_R', data=ces%t_grnd_r, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! pft energy state variable - t_ref2m_min_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_R', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='rural daily minimum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MIN_R', data=pes%t_ref2m_min_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_R', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='rural daily maximum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MAX_R', data=pes%t_ref2m_max_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_min_inst_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_INST_R', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='rural instantaneous daily min of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MIN_INST_R', data=pes%t_ref2m_min_inst_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_inst_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_INST_R', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='rural instantaneous daily max of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_REF2M_MAX_INST_R', data=pes%t_ref2m_max_inst_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

    ! column energy state variable - t_soisno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_SOISNO', xtype=ncd_double,   &
            dim1name='column', dim2name='levtot', switchdim=.true., &
            long_name='soil-snow temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_SOISNO', data=ces%t_soisno, &
            dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type energy state variable - t_lake

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_LAKE', xtype=ncd_double,  &
            dim1name='column', dim2name='levlak', switchdim=.true., &
            long_name='lake temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_LAKE', data=ces%t_lake, &
            dim1name=namec, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft physical state variable - frac_veg_nosno_alb

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FRAC_VEG_NOSNO_ALB', xtype=ncd_int,  &
            dim1name='pft',&
            long_name='fraction of vegetation not covered by snow (0 or 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='FRAC_VEG_NOSNO_ALB', data=pps%frac_veg_nosno_alb, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - fwet

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FWET', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='fraction of canopy that is wet (0 to 1)', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='FWET', data=pps%fwet, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - tlai

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tlai', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='one-sided leaf area index, no burying by snow', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tlai', data=pps%tlai, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - tsai

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tsai', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='one-sided stem area index, no burying by snow', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tsai', data=pps%tsai, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mlaidiff', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='difference between lai month one and month two',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mlaidiff', data=pps%mlaidiff, &
            dim1name='pft', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - elai

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='elai', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='one-sided leaf area index, with burying by snow', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='elai', data=pps%elai, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - esai

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='esai', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='one-sided stem area index, with burying by snow', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='esai', data=pps%esai, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - fsun

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fsun', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='sunlit fraction of canopy', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fsun', data=pps%fsun, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' )then
          if ( .not. readvar) then
             if (is_restart()) call endrun()
          else
             do p = begp, endp
                if ( shr_infnan_isnan( pps%fsun(p) ) )then
                   pps%fsun(p) = spval
                end if
             end do
          end if
       end if
    end if

    ! pft type physical state variable - htop

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='htop', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='canopy top', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='htop', data=pps%htop, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - hbot

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='hbot', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='canopy botton', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='hbot', data=pps%hbot, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - fabd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fabd', xtype=ncd_double,  &
            dim1name='pft', dim2name='numrad', switchdim=.true., &
            long_name='flux absorbed by veg per unit direct flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fabd', data=pps%fabd, &
            dim1name=namep, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - fabi

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fabi', xtype=ncd_double,  &
            dim1name='pft', dim2name='numrad', switchdim=.true., &
            long_name='flux absorbed by veg per unit diffuse flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fabi', data=pps%fabi, &
            dim1name=namep, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - ftdd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ftdd', xtype=ncd_double,  &
            dim1name='pft', dim2name='numrad', switchdim=.true., &
            long_name='down direct flux below veg per unit direct flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='ftdd', data=pps%ftdd, &
            dim1name=namep, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - ftid

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ftid', xtype=ncd_double,  &
            dim1name='pft', dim2name='numrad', switchdim=.true., &
            long_name='down diffuse flux below veg per unit direct flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='ftid', data=pps%ftid, &
            dim1name=namep, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - ftii

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ftii', xtype=ncd_double,  &
            dim1name='pft', dim2name='numrad', switchdim=.true., &
            long_name='down diffuse flux below veg per unit diffuse flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='ftii', data=pps%ftii, &
            dim1name=namep, switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft energy state variable - t_veg

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_VEG', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='vegetation temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_VEG', data=pes%t_veg, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft energy state variable - t_ref2m

    varname = 'T_REF2M'
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_double,  &
            dim1name='pft', &
            long_name='2m height surface air temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname=varname, data=pes%t_ref2m, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (allocate_all_vegpfts) then
             call endrun()
          else
             if (is_restart()) call endrun()
          end if
       end if
    end if

    ! pft type water state variable - h2ocan

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='H2OCAN', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='canopy water', units='kg/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='H2OCAN', data=pws%h2ocan, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column irrigation variable - n_irrig_steps_left

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='n_irrig_steps_left', xtype=ncd_int,  &
            dim1name='column', &
            long_name='number of irrigation time steps left', units='#')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='n_irrig_steps_left', data=cps%n_irrig_steps_left, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
          cps%n_irrig_steps_left = 0
       end if
    end if

   ! column irrigation variable - irrig_rate

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='irrig_rate', xtype=ncd_double,  &
            dim1name='column', &
            long_name='irrigation rate', units='mm/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='irrig_rate', data=cps%irrig_rate, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
          cps%irrig_rate = 0.0_r8
       end if
    end if

    ! ------------------------------------------------------------
    ! Determine volumetric soil water (for read only)
    ! ------------------------------------------------------------

    if (flag == 'read' ) then
       do c = begc,endc
          l = clandunit(c)
          if ( ltype(l) == istdlak .or. ltype(l) == istslak )then
             nlevs = nlevlak
          else if ( ltype(l) == isturb )then
             nlevs = nlevurb
          else
             nlevs = nlevgrnd
          end if
          ! NOTE: THIS IS A MEMORY INEFFICIENT COPY
          do j = 1,nlevs
             cws%h2osoi_vol(c,j) = cws%h2osoi_liq(c,j)/(cps%dz(c,j)*denh2o) &
                                      + cws%h2osoi_ice(c,j)/(cps%dz(c,j)*denice)
          end do
       end do


       ! ------------------------------------------------------------
       ! If initial run -- ensure that water is properly bounded
       ! ------------------------------------------------------------

       if ( is_first_step() )then
          do c = begc,endc
             l = clandunit(c)
             if (      ltype(l) == istdlak .or. ltype(l) == istslak )then
                nlevs = nlevlak
             else if ( ltype(l) == isturb )then
                nlevs = nlevurb
             else
                nlevs = nlevgrnd
             end if
             do j = 1,nlevs
                l = clandunit(c)
                if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
                   cws%h2osoi_liq(c,j) = max(0._r8,cws%h2osoi_liq(c,j))
                   cws%h2osoi_ice(c,j) = max(0._r8,cws%h2osoi_ice(c,j))
                   cws%h2osoi_vol(c,j) = cws%h2osoi_liq(c,j)/(cps%dz(c,j)*denh2o) &
                                              + cws%h2osoi_ice(c,j)/(cps%dz(c,j)*denice)
                   if (j == 1) then
                      maxwatsat = (cps%watsat(c,j)*cps%dz(c,j)*1000.0_r8 + pondmx) / &
                                  (cps%dz(c,j)*1000.0_r8)
                   else
                      maxwatsat = cps%watsat(c,j)
                   end if
                   if (cws%h2osoi_vol(c,j) > maxwatsat) then 
                      excess = (cws%h2osoi_vol(c,j) - maxwatsat)*cps%dz(c,j)*1000.0_r8
                      totwat = cws%h2osoi_liq(c,j) + cws%h2osoi_ice(c,j)
                      cws%h2osoi_liq(c,j) = cws%h2osoi_liq(c,j) - &
                        (cws%h2osoi_liq(c,j)/totwat) * excess
                      cws%h2osoi_ice(c,j) = cws%h2osoi_ice(c,j) - &
                        (cws%h2osoi_ice(c,j)/totwat) * excess
                   end if
                   cws%h2osoi_liq(c,j) = max(watmin,cws%h2osoi_liq(c,j))
                   cws%h2osoi_ice(c,j) = max(watmin,cws%h2osoi_ice(c,j))
                   cws%h2osoi_vol(c,j) = cws%h2osoi_liq(c,j)/(cps%dz(c,j)*denh2o) &
                                              + cws%h2osoi_ice(c,j)/(cps%dz(c,j)*denice)
                end if
             end do
          end do
       end if
    endif
    !
    ! variables needed for SNICAR
    !
    ! column type physical state variable - snw_rds
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='snw_rds', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer effective radius', units='um')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='snw_rds', data=cps%snw_rds, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find snw_rds in restart (or initial) file..."
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize snw_rds
             if (masterproc) then
                write(iulog,*) "SNICAR: This is an initial run (not a restart), and grain size/aerosol " // &
                               "mass data are not defined in initial condition file. Initialize snow " // &
                               "effective radius to fresh snow value, and snow/aerosol masses to zero."
             endif
             do c=begc,endc
                if (cps%snl(c) < 0) then
                   cps%snw_rds(c,cps%snl(c)+1:0) = snw_rds_min
                   cps%snw_rds(c,-nlevsno+1:cps%snl(c)) = 0._r8
                   cps%snw_rds_top(c) = snw_rds_min
                   cps%sno_liq_top(c) = cws%h2osoi_liq(c,cps%snl(c)+1) / &
                        (cws%h2osoi_liq(c,cps%snl(c)+1)+cws%h2osoi_ice(c,cps%snl(c)+1))
                elseif (cws%h2osno(c) > 0._r8) then
                   cps%snw_rds(c,0) = snw_rds_min
                   cps%snw_rds(c,-nlevsno+1:-1) = 0._r8
                   cps%snw_rds_top(c) = spval
                   cps%sno_liq_top(c) = spval
                else
                   cps%snw_rds(c,:) = 0._r8
                   cps%snw_rds_top(c) = spval
                   cps%sno_liq_top(c) = spval
                endif
             enddo
          endif
       end if
    end if

    ! column type physical state variable - mss_bcpho
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mss_bcpho', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer hydrophobic black carbon mass', units='kg m-2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mss_bcpho', data=cps%mss_bcpho, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize mss_bcpho to zero
             do c=begc,endc
                cps%mss_bcpho(c,-nlevsno+1:0) = 0._r8
             enddo
          endif
       end if
    end if

    ! column type physical state variable - mss_bcphi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mss_bcphi', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer hydrophilic black carbon mass', units='kg m-2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mss_bcphi', data=cps%mss_bcphi, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize mss_bcphi to zero
             do c=begc,endc
                cps%mss_bcphi(c,-nlevsno+1:0) = 0._r8
             enddo
          endif
       end if
    end if

    ! column type physical state variable - mss_ocpho
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mss_ocpho', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer hydrophobic organic carbon mass', units='kg m-2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mss_ocpho', data=cps%mss_ocpho, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize mss_ocpho to zero
             do c=begc,endc
                cps%mss_ocpho(c,-nlevsno+1:0) = 0._r8
             enddo
          endif
       end if
    end if

    ! column type physical state variable - mss_ocphi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mss_ocphi', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer hydrophilic organic carbon mass', units='kg m-2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mss_ocphi', data=cps%mss_ocphi, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize mss_ocphi to zero
             do c=begc,endc
                cps%mss_ocphi(c,-nlevsno+1:0) = 0._r8
             enddo
          endif
       end if
    end if

    ! column type physical state variable - mss_dst1
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mss_dst1', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer dust species 1 mass', units='kg m-2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mss_dst1', data=cps%mss_dst1, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize mss_dst1 to zero
             do c=begc,endc
                cps%mss_dst1(c,-nlevsno+1:0) = 0._r8
             enddo
          endif
       end if
    end if

    ! column type physical state variable - mss_dst2
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mss_dst2', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer dust species 2 mass', units='kg m-2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mss_dst2', data=cps%mss_dst2, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize mss_dst2 to zero
             do c=begc,endc
                cps%mss_dst2(c,-nlevsno+1:0) = 0._r8
             enddo
          endif
       end if
    end if

    ! column type physical state variable - mss_dst3
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mss_dst3', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer dust species 3 mass', units='kg m-2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mss_dst3', data=cps%mss_dst3, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize mss_dst3 to zero
             do c=begc,endc
                cps%mss_dst3(c,-nlevsno+1:0) = 0._r8
             enddo
          endif
       end if
    end if

    ! column type physical state variable - mss_dst4
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mss_dst4', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer dust species 4 mass', units='kg m-2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mss_dst4', data=cps%mss_dst4, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize mss_dst4 to zero
             do c=begc,endc
                cps%mss_dst4(c,-nlevsno+1:0) = 0._r8
             enddo
          endif
       end if
    end if

    ! column type physical state variable - flx_absdv
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='flx_absdv', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno1', switchdim=.true., &
            long_name='snow layer flux absorption factors (direct, VIS)', units='fraction')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='flx_absdv', data=cps%flx_absdv, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=1, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
          ! SNICAR, via SurfaceAlbedo, will define the needed flux absorption factors
          if (nsrest == nsrStartup) do_initsurfalb = .true.
       end if
    end if

    ! column type physical state variable - flx_absdn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='flx_absdn', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno1', switchdim=.true., &
            long_name='snow layer flux absorption factors (direct, NIR)', units='fraction')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='flx_absdn', data=cps%flx_absdn, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=1, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
          ! SNICAR, via SurfaceAlbedo, will define the needed flux absorption factors
          if (nsrest == nsrStartup) do_initsurfalb = .true.
       end if
    end if

    ! column type physical state variable - flx_absiv
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='flx_absiv', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno1', switchdim=.true., &
            long_name='snow layer flux absorption factors (diffuse, VIS)', units='fraction')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='flx_absiv', data=cps%flx_absiv, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=1, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
          ! SNICAR, via SurfaceAlbedo, will define the needed flux absorption factors
          if (nsrest == nsrStartup) do_initsurfalb = .true.
       end if
    end if

    ! column type physical state variable - flx_absin
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='flx_absin', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno1', switchdim=.true., &
            long_name='snow layer flux absorption factors (diffuse, NIR)', units='fraction')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='flx_absin', data=cps%flx_absin, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=1, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
          ! SNICAR, via SurfaceAlbedo, will define the needed flux absorption factors
          if (nsrest == nsrStartup) do_initsurfalb = .true.
       end if
    end if

    ! column type physical state variable - albsnd_hst
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albsnd_hst', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='snow albedo (direct) (0 to 1)',units='proportion')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='albsnd_hst', data=cps%albsnd_hst, &
            dim1name='column', switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - albsni_hst
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albsni_hst', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='snow albedo (diffuse) (0 to 1)',units='proportion')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='albsni_hst', data=cps%albsni_hst, &
            dim1name='column', switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type water flux variable - qflx_snofrz_lyr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='qflx_snofrz_lyr', xtype=ncd_double,  &
            dim1name='column', dim2name='levsno', switchdim=.true., &
            long_name='snow layer ice freezing rate', units='kg m-2 s-1')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='qflx_snofrz_lyr', data=cwf%qflx_snofrz_lyr, &
            dim1name='column', switchdim=.true., &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize qflx_snofrz_lyr to zero
             do c=begc,endc
                cwf%qflx_snofrz_lyr(c,-nlevsno+1:0) = 0._r8
             enddo
          endif
       end if
    end if

    ! column type water flux variable - qflx_snow_melt
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='qflx_snow_melt', xtype=ncd_double,  &
            dim1name='column', long_name='net snow melt', units='mm/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='qflx_snow_melt', data=cwf%qflx_snow_melt, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             ! initial run, not restart: initialize qflx_snow_melt to zero
             cwf%qflx_snow_melt = 0._r8
          endif
       end if
    end if

    ! initialize other variables that are derived from those
    ! stored in the restart buffer. (there may be a more appropriate
    ! place to do this, but functionally this works)
   if (flag == 'read' ) then
       do j = -nlevsno+1,0
          do c = begc,endc
             ! mass concentrations of aerosols in snow
             if (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j) > 0._r8) then
                cps%mss_cnc_bcpho(c,j) = cps%mss_bcpho(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_bcphi(c,j) = cps%mss_bcphi(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_ocpho(c,j) = cps%mss_ocpho(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_ocphi(c,j) = cps%mss_ocphi(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
             
                cps%mss_cnc_dst1(c,j) = cps%mss_dst1(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_dst2(c,j) = cps%mss_dst2(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_dst3(c,j) = cps%mss_dst3(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_dst4(c,j) = cps%mss_dst4(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
             else
                cps%mss_cnc_bcpho(c,j) = 0._r8
                cps%mss_cnc_bcphi(c,j) = 0._r8
                cps%mss_cnc_ocpho(c,j) = 0._r8
                cps%mss_cnc_ocphi(c,j) = 0._r8
             
                cps%mss_cnc_dst1(c,j) = 0._r8
                cps%mss_cnc_dst2(c,j) = 0._r8
                cps%mss_cnc_dst3(c,j) = 0._r8
                cps%mss_cnc_dst4(c,j) = 0._r8
             endif
          enddo
       enddo
    endif
    !-- SNICAR variables


  end subroutine BiogeophysRest

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: weights_exactly_the_same
!
! !INTERFACE:
  logical function weights_exactly_the_same( pptr, wtgcell, wtlunit, wtcol )
!
! !DESCRIPTION:
! Determine if the weights read in are exactly the same as those from surface dataset
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    type(pft_type), pointer :: pptr       ! pointer to pft derived subtype
    real(r8), intent(IN)    :: wtgcell(:) ! grid cell weights for each PFT
    real(r8), intent(IN)    :: wtlunit(:) ! land-unit weights for each PFT
    real(r8), intent(IN)    :: wtcol(:)   ! column weights for each PFT
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!-----------------------------------------------------------------------

     ! Check that weights are identical for all PFT's and all weight types
     if (  all( pptr%wtgcell(:) == wtgcell ) .and. all( pptr%wtlunit(:) == wtlunit ) &
     .and. all( pptr%wtcol(:) == wtcol ) )then
        weights_exactly_the_same = .true.
     else
        weights_exactly_the_same = .false.
     end if

  end function weights_exactly_the_same

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: weights_within_roundoff_different
!
! !INTERFACE:
  logical function weights_within_roundoff_different( pptr, wtgcell, wtlunit, wtcol )
!
! !DESCRIPTION:
! Determine if the weights are within roundoff different from each other
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    type(pft_type), pointer :: pptr       ! pointer to pft derived subtype
    real(r8), intent(IN)    :: wtgcell(:) ! grid cell weights for each PFT
    real(r8), intent(IN)    :: wtlunit(:) ! land-unit weights for each PFT
    real(r8), intent(IN)    :: wtcol(:)   ! column weights for each PFT
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!-----------------------------------------------------------------------
     real(r8), parameter :: rndVal  = 1.e-13_r8

     ! If differences between all weights for each PFT and each weight type is
     ! less than or equal to double precision roundoff level -- weights are close
     if (  all( abs(pptr%wtgcell(:) - wtgcell) <= rndVal ) &
     .and. all( abs(pptr%wtlunit(:) - wtlunit) <= rndVal ) &
     .and. all( abs(pptr%wtcol(:)   - wtcol  ) <= rndVal ) )then
        weights_within_roundoff_different = .true.
     else
        weights_within_roundoff_different = .false.
     end if

  end function weights_within_roundoff_different

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: weights_tooDifferent
!
! !INTERFACE:
  logical function weights_tooDifferent( begp, endp, pptr, wtgcell, adiff, maxdiff )
!
! !DESCRIPTION:
! Determine if the weights read in are too different and should flag an error
!
! !USES:
    use clmtype
    implicit none
!
! !ARGUMENTS:
    integer, intent(IN)     :: begp, endp         ! per-proc beginning and ending pft indices
    type(pft_type), pointer :: pptr               ! pointer to pft derived subtype
    real(r8), intent(IN)    :: wtgcell(begp:endp) ! grid cell weights for each PFT
    real(r8), intent(IN)    :: adiff              ! tolerance of acceptible difference
    real(r8), intent(OUT)   :: maxdiff            ! maximum difference found
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!-----------------------------------------------------------------------
     integer  :: p        ! PFT index
     real(r8) :: diff     ! difference in weights

     ! Assume weights are NOT different and only change if find weights too different
     weights_tooDifferent = .false.
     maxdiff = 0.0_r8
     do p = begp, endp

        diff = abs(pptr%wtgcell(p) - wtgcell(p))
        if ( diff > maxdiff ) maxdiff = diff
        if ( diff > adiff   ) weights_tooDifferent = .true.
     end do

  end function weights_tooDifferent


end module BiogeophysRestMod
