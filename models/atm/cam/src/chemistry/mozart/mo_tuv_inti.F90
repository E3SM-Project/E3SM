
module mo_tuv_inti

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  private
  public :: nj
  public :: nlng, nzen, ncof
  public :: tuv_inti

  save

  integer :: nj
  integer :: nlng
  integer :: nzen, ncof

contains

  subroutine tuv_inti( nz, tuv_xsect_file, lng_indexer )
    !-----------------------------------------------------------------------------
    !   purpose:
    !   read 17 bins data outputed from tuv
    !-----------------------------------------------------------------------------
    !   parameters:
    !   nw     - integer, number of specified intervals + 1 in working
    !            wavelength grid
    !   wl     - real(r8), vector of lower limits of wavelength intervals in
    !            working wavelength grid
    !   wc     - real(r8), vector of center  of wavelength intervals in
    !            working wavelength grid
    !   wu     - real(r8), vector of upper limits of wavelength intervals in
    !            working wavelength grid
    !   f      - real(r8), spectral irradiance at the top of the atmosphere at
    !            each specified wavelength
    !-----------------------------------------------------------------------------
    !   edit history:        
    !   10/2000  similified by xuexi
    !-----------------------------------------------------------------------------

    use spmd_utils,    only : masterproc
    use cam_logfile,   only : iulog
    use abortutils,    only : endrun
    use mo_params,     only : kj, kw, smallest, largest
    use mo_waveall,    only : r01g1, r01g2, r01g3, r01g4, &
         r04g, r08g, r06g1, r06g2, &
         r10g1, r10g2, r10g3, r10g4, r10g5, &
         r11g, r11g1, r11g2, r11g3, r11g4, &
         r14g, r14g1, r14g2, &
         r15g, r15g1, r15g2, r15g3, &
         r17g, r17g1, &
         r18g, r18g2
    use mo_wavelab,    only : sj
    use mo_wavelen,    only : nw, deltaw, delw_bin, sflx, wc, wl, wu
    use mo_waveo3,     only : xso3, s226, s263, s298
    use mo_zadj,       only : adj_coeffs
    use mo_schu,       only : schu_inti
    use mo_xsections,  only : r44_inti, r08_inti
    use chem_mods,     only : phtcnt, pht_alias_lst, rxt_tag_lst
    use ioFileMod,     only : getfil
    use cam_pio_utils, only : cam_pio_openfile
    use pio,           only : file_desc_t, pio_nowrite, pio_closefile, &
         pio_inq_dimid, pio_inq_varid, pio_inq_dimlen, pio_get_var, &
         pio_seterrorhandling, pio_bcast_error, pio_internal_error, pio_noerr
    implicit none

    !-----------------------------------------------------------------------------
    !	... dummy arguments
    !-----------------------------------------------------------------------------
    integer, intent(in)    :: nz
    integer, intent(inout) :: lng_indexer(phtcnt)
    character(len=*), intent(in) :: tuv_xsect_file

    !-----------------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------------
    type(file_desc_t) :: ncid
    integer :: ndx
    integer :: dimid, vid
    integer :: iw, ios, iret
    integer :: k, m
    integer :: ind_wrk(4)
    integer :: wrk_ndx(phtcnt)
    real(r8), allocatable  :: coeff_adj(:,:)
    character(len=256) :: filespec
    character(len=256) :: locfn
    character(len=20)  :: coeff_tag

    !------------------------------------------------------------------------
    !     for wl(iw) .lt. 150.01                                susim_hi.flx
    !     for wl(iw) .ge. 150.01 and wl(iw) .le. 400            atlas3.flx 
    !     for wl(iw) .gt. 400                                   neckel & labs 
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    !     input data files including:
    !           (0) wavelength nw,wl,wc,wu
    !           (1) solar flux
    !           (2) o3 cross sections
    !           (3) other cross 
    !           (4) t dependence parameter of cross section
    !------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !       ... open netcdf file
    !---------------------------------------------------------------------------
    filespec = trim( tuv_xsect_file )
    call getfil( filespec, locfn, 0 )
    call cam_pio_openfile( ncid, trim(locfn), PIO_NOWRITE)
    !---------------------------------------------------------------------------
    ! 	... get the dimensions
    !---------------------------------------------------------------------------
    iret = pio_inq_dimid( ncid, 'nw', dimid )
    iret = pio_inq_dimlen( ncid, dimid, nw )
    iret = pio_inq_dimid( ncid, 'nzen', dimid )
    iret = pio_inq_dimlen( ncid, dimid, nzen )
    iret = pio_inq_dimid( ncid, 'ncof', dimid )
    iret = pio_inq_dimlen( ncid, dimid, ncof )
    !---------------------------------------------------------------------------
    ! 	... read the wave bin coordinates
    !---------------------------------------------------------------------------
    iret = pio_inq_varid( ncid, 'wl', vid )
    iret = pio_get_var( ncid, vid, wl(1:nw) )
    iret = pio_inq_varid( ncid, 'wc', vid )
    iret = pio_get_var( ncid, vid, wc(1:nw) )
    iret = pio_inq_varid( ncid, 'wu', vid )
    iret = pio_get_var( ncid, vid, wu(1:nw) )
    wl(nw+1) = wu(nw)
    write(iulog,*) ' '
    write(iulog,*) 'tuv_inti: wl(nw+1) = ',wl(nw+1)
    !------------------------------------------------------------------------
    !  	... solar flux
    !------------------------------------------------------------------------
    iret = pio_inq_varid( ncid, 'sflx', vid )
    iret = pio_get_var( ncid, vid, sflx(1:nw) )
    !------------------------------------------------------------------------
    !    	... o3 cross (t dependence)
    !------------------------------------------------------------------------
    iret = pio_inq_varid( ncid, 'xso3', vid )
    iret = pio_get_var( ncid, vid, xso3(1:nw) )
    iret = pio_inq_varid( ncid, 's226', vid )
    iret = pio_get_var( ncid, vid, s226(1:nw) )
    iret = pio_inq_varid( ncid, 's263', vid )
    iret = pio_get_var( ncid, vid, s263(1:nw) )
    iret = pio_inq_varid( ncid, 's298', vid )
    iret = pio_get_var( ncid, vid, s298(1:nw) )
    !---------------------------------------------------------------------------
    !   	... temperature dependent cross section parameters
    !---------------------------------------------------------------------------
    iret = pio_inq_varid( ncid, 'r01g1', vid )
    iret = pio_get_var( ncid, vid, r01g1(1:nw) )
    iret = pio_inq_varid( ncid, 'r01g2', vid )
    iret = pio_get_var( ncid, vid, r01g2(1:nw) )
    iret = pio_inq_varid( ncid, 'r01g3', vid )
    iret = pio_get_var( ncid, vid, r01g3(1:nw) )
    iret = pio_inq_varid( ncid, 'r01g4', vid )
    iret = pio_get_var( ncid, vid, r01g4(1:nw) )

    iret = pio_inq_varid( ncid, 'r04g', vid )
    iret = pio_get_var( ncid, vid, r04g(1:nw) )

    iret = pio_inq_varid( ncid, 'r08g', vid )
    iret = pio_get_var( ncid, vid, r08g(1:nw) )

    iret = pio_inq_varid( ncid, 'r06g1', vid )
    iret = pio_get_var( ncid, vid, r06g1(1:nw) )
    iret = pio_inq_varid( ncid, 'r06g2', vid )
    iret = pio_get_var( ncid, vid, r06g2(1:nw) )

    iret = pio_inq_varid( ncid, 'r10g1', vid )
    iret = pio_get_var( ncid, vid, r10g1(1:nw) )
    iret = pio_inq_varid( ncid, 'r10g2', vid )
    iret = pio_get_var( ncid, vid, r10g2(1:nw) )
    iret = pio_inq_varid( ncid, 'r10g3', vid )
    iret = pio_get_var( ncid, vid, r10g3(1:nw) )
    iret = pio_inq_varid( ncid, 'r10g4', vid )
    iret = pio_get_var( ncid, vid, r10g4(1:nw) )
    iret = pio_inq_varid( ncid, 'r10g5', vid )
    iret = pio_get_var( ncid, vid, r10g5(1:nw) )

    iret = pio_inq_varid( ncid, 'r11g', vid )
    iret = pio_get_var( ncid, vid, r11g(1:nw) )
    iret = pio_inq_varid( ncid, 'r11g1', vid )
    iret = pio_get_var( ncid, vid, r11g1(1:nw) )
    iret = pio_inq_varid( ncid, 'r11g2', vid )
    iret = pio_get_var( ncid, vid, r11g2(1:nw) )
    iret = pio_inq_varid( ncid, 'r11g3', vid )
    iret = pio_get_var( ncid, vid, r11g3(1:nw) )
    iret = pio_inq_varid( ncid, 'r11g4', vid )
    iret = pio_get_var( ncid, vid, r11g4(1:nw) )

    iret = pio_inq_varid( ncid, 'r14g', vid )
    iret = pio_get_var( ncid, vid, r14g(1:nw) )
    iret = pio_inq_varid( ncid, 'r14g1', vid )
    iret = pio_get_var( ncid, vid, r14g1(1:nw) )
    iret = pio_inq_varid( ncid, 'r14g2', vid )
    iret = pio_get_var( ncid, vid, r14g2(1:nw) )

    iret = pio_inq_varid( ncid, 'r15g', vid )
    iret = pio_get_var( ncid, vid, r15g(1:nw) )
    iret = pio_inq_varid( ncid, 'r15g1', vid )
    iret = pio_get_var( ncid, vid, r15g1(1:nw) )
    iret = pio_inq_varid( ncid, 'r15g2', vid )
    iret = pio_get_var( ncid, vid, r15g2(1:nw) )
    iret = pio_inq_varid( ncid, 'r15g3', vid )
    iret = pio_get_var( ncid, vid, r15g3(1:nw) )

    iret = pio_inq_varid( ncid, 'r17g', vid )
    iret = pio_get_var( ncid, vid, r17g(1:nw) )
    iret = pio_inq_varid( ncid, 'r17g1', vid )
    iret = pio_get_var( ncid, vid, r17g1(1:nw) )

    iret = pio_inq_varid( ncid, 'r18g', vid )
    iret = pio_get_var( ncid, vid, r18g(1:nw) )
    iret = pio_inq_varid( ncid, 'r18g2', vid )
    iret = pio_get_var( ncid, vid, r18g2(1:nw) )
    !------------------------------------------------------------------------------
    !       ... check for cross section in dataset
    !------------------------------------------------------------------------------
    call pio_seterrorhandling(ncid, pio_bcast_error)
    do m = 1,phtcnt
       if( pht_alias_lst(m,2) == ' ' ) then
          iret = pio_inq_varid( ncid, rxt_tag_lst(m), vid )
          if( iret == pio_noerr ) then 
             lng_indexer(m) = vid
          end if
       else if( pht_alias_lst(m,2) == 'userdefined' ) then
          lng_indexer(m) = -1
       else
          iret = pio_inq_varid( ncid, trim(pht_alias_lst(m,2)), vid )
          if( iret == pio_noerr ) then 
             lng_indexer(m) = vid
          else
             write(iulog,*) 'tuv_inti : ',rxt_tag_lst(m)(:len_trim(rxt_tag_lst(m))),' alias ', &
                  pht_alias_lst(m,2)(:len_trim(pht_alias_lst(m,2))),' not in dataset'            
             call endrun
          end if
       end if
    end do
    call pio_seterrorhandling(ncid, pio_internal_error)
    nlng = 0
    do m = 1,phtcnt
       if( lng_indexer(m) > 0 ) then
          if( any( lng_indexer(:m-1) == lng_indexer(m) ) ) then
             cycle
          end if
          nlng = nlng + 1
       end if
    end do
    !---------------------------------------------------------------------------
    !	... allocate the cross section array
    !---------------------------------------------------------------------------
    allocate( sj(nw,nz,nlng), adj_coeffs(ncof,nlng,nzen), coeff_adj(ncof,nzen), stat=ios )
    if( ios /= 0 ) then
       write(iulog,*) 'tuv_inti: failed to allocate sj ... coeff_adj; error = ',ios
       call endrun
    end if
    sj(:,:,:)         = 0._r8
    adj_coeffs(:,:,:) = 0._r8
    write(iulog,*) 'tuv_inti: nlng = ',nlng
    write(iulog,*) 'tuv_inti: lng_indexer'
    write(iulog,'(10i5)') lng_indexer(:)
    if( nlng > 0 ) then
       write(iulog,*) ' '
       write(iulog,*) 'tuv_inti: photo xsect analysis'
       do m = 1,phtcnt
          if( lng_indexer(m) > 0 ) then
             write(iulog,*) trim(rxt_tag_lst(m)),lng_indexer(m)
          end if
       end do
    end if
    ndx = 0
    do m = 1,phtcnt
       if( lng_indexer(m) > 0 ) then
          if( any( lng_indexer(:m-1) == lng_indexer(m) ) ) then
             cycle
          end if
          ndx = ndx + 1
          iret = pio_get_var( ncid, lng_indexer(m), sj(1:nw,1,ndx) )

          do k = 2,nz
             sj(:,k,ndx) = sj(:,1,ndx)
          end do
          coeff_tag = trim(rxt_tag_lst(m)) // '_adj'
          iret = pio_inq_varid( ncid, trim(coeff_tag), vid )
          iret = pio_get_var( ncid, vid, coeff_adj )

          adj_coeffs(:,ndx,1:nzen) = coeff_adj(:,1:nzen)
       end if
    end do
    if( ndx /= nlng ) then
       write(iulog,*) 'tuv_inti : ndx count /= cross section count'
       call endrun
    end if
    !------------------------------------------------------------------------------
    !       ... setup final lng_indexer
    !------------------------------------------------------------------------------
    ndx = 0
    wrk_ndx(:) = lng_indexer(:)
    do m = 1,phtcnt
       if( wrk_ndx(m) > 0 ) then
          ndx = ndx + 1
          k = wrk_ndx(m)
          where( wrk_ndx(:) == k )
             lng_indexer(:) = ndx
             wrk_ndx(:)     = -100000
          end where
       end if
    end do
    if( nlng > 0 ) then
       write(iulog,*) ' '
       write(iulog,*) 'tuv_inti: photo xsect analysis'
       do m = 1,phtcnt
          if( lng_indexer(m) > 0 ) then
             write(iulog,*) trim(rxt_tag_lst(m)),lng_indexer(m)
          end if
       end do
    end if
    !---------------------------------------------------------------------------
    ! 	... close netcdf file
    !---------------------------------------------------------------------------
    call pio_closefile( ncid )
    deallocate( coeff_adj )


    delw_bin(:nw) = wu(:nw) - wl(:nw)
    deltaw(:nw)   = delw_bin(:nw) * wc(:nw) * 5.039e11_r8
    delw_bin(:nw) = 1._r8/delw_bin(:nw)
    largest  = huge( largest )
    smallest = tiny( largest )

    write(iulog,'(''tuv_inti: smallest,largest = '',1p,2e21.13)') smallest,largest

    call schu_inti
    call r44_inti( nw, wc )
    call r08_inti( nw, wl, wc )

  end subroutine tuv_inti

end module mo_tuv_inti
