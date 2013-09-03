#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module aquaplanet_io_mod
  use kinds, only : real_kind
  use common_io_mod, only : nf_selectedvar, nf_handle, nfsizekind
#ifdef PIO_INTERP
  use pio_io_mod, only :  nf_output_register_variables, &
       nf_put_var, nf_get_frame
#else
  use netcdf_io_mod, only :  nf_output_register_variables, nf_put_var, nf_get_frame !_EXTERNAL
#endif
  use dof_mod, only : UniquePoints

  implicit none
  private

  integer, parameter :: varcnt=8
#ifdef PIO_INTERP
  integer, parameter :: maxdims=4
#else
  integer, parameter :: maxdims=3
#endif
 public :: aq_movie_init
 public :: aq_movie_output
 public :: aq_set_varnames
! public :: aq_movie_finish
 character*(*), parameter :: varnames(varcnt)=(/'udrag', &
                                                'vdrag', &
                                                'tsflx', &
                                                'qsflx', &
                                                'usfrc', &
                                                'vsfrc', &
                                                'tsfrc', &
                                                'qsfrc'/)
contains
  subroutine aq_movie_init(ncdf)
    type(nf_handle), intent(inout) :: ncdf(:)
#ifdef PIO_INTERP
    ! Refer to interp_movie_mod for registered dims
    integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,2,5,0,  &
                                                                 1,2,5,0,  &
                                                                 1,2,5,0,  &
                                                                 1,2,5,0,  &
                                                                 1,2,3,5,  &
                                                                 1,2,3,5,  &
                                                                 1,2,3,5,  &
                                                                 1,2,3,5/),  &
                                                                 shape=(/maxdims,varcnt/))
#else
    ! Refer to prim_movie_mod for registered dims
    integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,4,0,  &
                                                                 1,4,0,  &
                                                                 1,4,0,  &
                                                                 1,4,0,  &
                                                                 1,2,4,  &
                                                                 1,2,4,  &
                                                                 1,2,4,  &
                                                                 1,2,4/),  &
                                                                 shape=(/maxdims,varcnt/))
#endif
  call nf_output_register_variables(ncdf,varcnt,varnames,vardims)

  end subroutine aq_movie_init

  subroutine aq_set_varnames(lvarcnt,nlvarnames)
    character*(*), intent(out) :: nlvarnames(:)
    integer, intent(inout) :: lvarcnt
    lvarcnt=lvarcnt+varcnt
    nlvarnames=varnames
  end subroutine aq_set_varnames

#ifdef PIO_INTERP
  subroutine aq_movie_output(ncdf, elem, interpdata, output_varnames, nc, nlev)
    use aquaplanet, only : udrag, vdrag, qsflx, tsflx, usf, vsf, tsf, qsf
    use dimensions_mod, only : np, nelemd
    use element_mod, only : element_t
    use interpolate_mod, only: interpdata_t, interpolate_scalar
    type(nf_handle), intent(inout) :: ncdf
    type(element_t), intent(in) :: elem(:)
    type(interpdata_t), intent(in) :: interpdata(:)
    integer, intent(in) :: nc , nlev
    character*(*), intent(in) :: output_varnames(:)
    real(kind=real_kind) :: var2d(nc),var3d(nc,nlev)
    integer(kind=nfsizekind) :: start(4), count(4)
    integer :: ie, st, en

! start and count are only used for time dimension in pio
    start(3)=nf_get_frame(ncdf)
    start(4)=nf_get_frame(ncdf)
    if(nf_selectedvar('udrag', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),udrag(:,:,ie), &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp
       enddo
       call nf_put_var(ncdf, var2d,start(1:3), count(1:3), name='udrag')
    endif
    if(nf_selectedvar('vdrag', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),vdrag(:,:,ie), &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp
       enddo
       call nf_put_var(ncdf, var2d,start(1:3), count(1:3), name='vdrag')
    endif
    if(nf_selectedvar('tsflx', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),tsflx(:,:,ie), &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp
       enddo
       call nf_put_var(ncdf, var2d,start(1:3), count(1:3), name='tsflx')
    endif
    if(nf_selectedvar('qsflx', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),qsflx(:,:,ie), &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp
       enddo
       call nf_put_var(ncdf, var2d,start(1:3), count(1:3), name='qsflx')
    endif


    if(nf_selectedvar('usfrc', output_varnames)) then
       st=1	
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),usf(:,:,:,ie), &
               np, nlev, var3d(st:en,:))
          st=st+interpdata(ie)%n_interp
       end do
       call nf_put_var(ncdf, var3d,start, count, name='usfrc')
    end if

    if(nf_selectedvar('vsfrc', output_varnames)) then
       st=1	
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),vsf(:,:,:,ie), &
               np, nlev, var3d(st:en,:))
          st=st+interpdata(ie)%n_interp
       end do
       call nf_put_var(ncdf, var3d,start, count, name='vsfrc')
    end if

    if(nf_selectedvar('tsfrc', output_varnames)) then
       st=1	
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),tsf(:,:,:,ie), &
               np, nlev, var3d(st:en,:))
          st=st+interpdata(ie)%n_interp
       end do
       call nf_put_var(ncdf, var3d,start, count, name='tsfrc')
    end if



    if(nf_selectedvar('qsfrc', output_varnames)) then
       st=1	
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),qsf(:,:,:,ie), &
               np, nlev, var3d(st:en,:))
          st=st+interpdata(ie)%n_interp
       end do
       call nf_put_var(ncdf, var3d,start, count, name='qsfrc')
    end if

  end subroutine aq_movie_output

#else

  subroutine aq_movie_output(ncdf, elem, output_varnames, nc, nlev)
    use aquaplanet, only : udrag, vdrag, qsflx, tsflx, usf, vsf, tsf, qsf
    use dimensions_mod, only : nelemd, nelemdmax
    use element_mod, only : element_t
    type(nf_handle), intent(inout) :: ncdf
    type(element_t), intent(in) :: elem(:)
    integer, intent(in) :: nc , nlev
    character*(*), intent(in) :: output_varnames(:)
    real(kind=real_kind) :: var2d(nc),var3d(nc,nlev)
    integer(kind=nfsizekind) :: start(3), count(3)
    integer :: ie, st, en

! start and count are only used for time dimension in pio
    start(2)=nf_get_frame(ncdf)
    start(3)=nf_get_frame(ncdf)
    if(nf_selectedvar('udrag', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,udrag(:,:,ie),var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d,start(1:2), count(1:2), name='udrag')
    endif
    if(nf_selectedvar('vdrag', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,vdrag(:,:,ie),var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d,start(1:2), count(1:2), name='vdrag')
    endif
    if(nf_selectedvar('tsflx', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,tsflx(:,:,ie),var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d,start(1:2), count(1:2), name='tsflx')
    endif
    if(nf_selectedvar('qsflx', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,qsflx(:,:,ie),var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d,start(1:2), count(1:2), name='qsflx')
    endif


    if(nf_selectedvar('usfrc', output_varnames)) then
       st=1	
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,nlev,usf(:,:,:,ie),var3d(st:en,:))
          st=en+1
       end do
       call nf_put_var(ncdf, var3d,start, count, name='usfrc')
    end if

    if(nf_selectedvar('vsfrc', output_varnames)) then
       st=1	
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,nlev,vsf(:,:,:,ie),var3d(st:en,:))
          st=en+1
       end do
       call nf_put_var(ncdf, var3d,start, count, name='vsfrc')
    end if

    if(nf_selectedvar('tsfrc', output_varnames)) then
       st=1	
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,nlev,tsf(:,:,:,ie),var3d(st:en,:))
          st=en+1
       end do
       call nf_put_var(ncdf, var3d,start, count, name='tsfrc')
    end if



    if(nf_selectedvar('qsfrc', output_varnames)) then
       st=1	
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,nlev,qsf(:,:,:,ie),var3d(st:en,:))
          st=en+1
       end do
       call nf_put_var(ncdf, var3d,start, count, name='qsfrc')
    end if

  end subroutine aq_movie_output
#endif

! obsolete.  HOMME now only supports PIO/native grid or PIO_INTERP
#if defined(NETCDFXXX) 
  subroutine aq_movie_output(ncdf, elem, output_varnames, nc, nlev)
    use aquaplanet, only : udrag, vdrag, qsflx, tsflx, usf, vsf, tsf, qsf
    use dimensions_mod, only : nelemd, nelemdmax
    use element_mod, only : element_t
    type(nf_handle), intent(inout) :: ncdf
    type(element_t), intent(in) :: elem(:)
    integer, intent(in) :: nc , nlev
    character*(*), intent(in) :: output_varnames(:)
    real(kind=real_kind) :: var2d(nc),var3d(nc,nlev)
    integer(kind=nfsizekind) :: start(3), count(3)
    integer :: ie, ncnt

    if(nf_selectedvar('udrag', output_varnames)) then
       start(1)=1
       start(2)=nf_get_frame(ncdf)
       count(1)=0
       count(2)=1
       do ie=1,nelemdmax          
	  if(ie <= nelemd) then
             start(1)=elem(ie)%idxP%UniquePtOffset
             count(1)=elem(ie)%idxP%NumUniquePts
             call UniquePoints(elem(ie)%idxP,udrag(:,:,ie),var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d,start, count, name='udrag')
       enddo
    endif
    if(nf_selectedvar('vdrag', output_varnames)) then
       start(1)=1
       start(2)=nf_get_frame(ncdf)
       count(1)=0
       count(2)=1

       do ie=1,nelemdmax
	  if(ie <= nelemd) then
             start(1)=elem(ie)%idxP%UniquePtOffset
             count(1)=elem(ie)%idxP%NumUniquePts
             call UniquePoints(elem(ie)%idxP,vdrag(:,:,ie), var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d,start, count, name='vdrag')
       enddo
    endif
    if(nf_selectedvar('tsflx', output_varnames)) then
       start(1)=1
       start(2)=nf_get_frame(ncdf)
       count(1)=0
       count(2)=1

       do ie=1,nelemdmax
	  if(ie <= nelemd) then
             start(1)=elem(ie)%idxP%UniquePtOffset
             count(1)=elem(ie)%idxP%NumUniquePts
             call UniquePoints(elem(ie)%idxP,tsflx(:,:,ie),var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d,start, count, name='tsflx')
       enddo
    endif
    if(nf_selectedvar('qsflx', output_varnames)) then
       start(1)=1
       start(2)=nf_get_frame(ncdf)
       count(1)=0
       count(2)=1

       do ie=1,nelemdmax
	  if(ie <= nelemd) then
             start(1)=elem(ie)%idxP%UniquePtOffset
             count(1)=elem(ie)%idxP%NumUniquePts
             call UniquePoints(elem(ie)%idxP,qsflx(:,:,ie),var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d,start, count, name='qsflx')
       enddo
    endif

    if(nf_selectedvar('usfrc', output_varnames)) then
       start(1)=1
       start(2)=1
       start(3)=nf_get_frame(ncdf)
       count(1)=0
       count(2)=nlev
       count(3)=1

       do ie=1,nelemdmax
	  if(ie <= nelemd) then
             start(1)=elem(ie)%idxP%UniquePtOffset
             count(1)=elem(ie)%idxP%NumUniquePts
             ncnt=count(1)
             call UniquePoints(elem(ie)%idxP,nlev,usf(:,:,:,ie),var3d)
          else
             ncnt=1
             count=0
          end if
          call nf_put_var(ncdf, var3d(1:ncnt,:),start, count, name='usfrc')
       enddo
    end if
    if(nf_selectedvar('vsfrc', output_varnames)) then
       start(1)=1
       start(2)=1
       start(3)=nf_get_frame(ncdf)
       count(1)=0
       count(2)=nlev
       count(3)=1

       do ie=1,nelemdmax
	  if(ie <= nelemd) then
             start(1)=elem(ie)%idxP%UniquePtOffset
             count(1)=elem(ie)%idxP%NumUniquePts
             ncnt=count(1)
             call UniquePoints(elem(ie)%idxP,nlev,vsf(:,:,:,ie),var3D)
          else
                   ncnt=1
             count=0
          end if
          call nf_put_var(ncdf, var3d(1:ncnt,:),start, count, name='vsfrc')
       enddo
    end if
    if(nf_selectedvar('tsfrc', output_varnames)) then
       start(1)=1
       start(2)=1
       start(3)=nf_get_frame(ncdf)
       count(1)=0
       count(2)=nlev
       count(3)=1

       do ie=1,nelemdmax
	  if(ie <= nelemd) then
             start(1)=elem(ie)%idxP%UniquePtOffset
             count(1)=elem(ie)%idxP%NumUniquePts
             ncnt=count(1)
             call UniquePoints(elem(ie)%idxP,nlev,tsf(:,:,:,ie),var3D)
          else
                   ncnt=1
             count=0
          end if
          call nf_put_var(ncdf, var3d(1:ncnt,:),start, count, name='tsfrc')
       enddo
    end if
    if(nf_selectedvar('qsfrc', output_varnames)) then
       start(1)=1
       start(2)=1
       start(3)=nf_get_frame(ncdf)
       count(1)=0
       count(2)=nlev
       count(3)=1

       do ie=1,nelemdmax
	  if(ie <= nelemd) then
             start(1)=elem(ie)%idxP%UniquePtOffset
             count(1)=elem(ie)%idxP%NumUniquePts
             ncnt=count(1)
             call UniquePoints(elem(ie)%idxP,nlev,qsf(:,:,:,ie),var3D)
          else
             ncnt=1
             count=0
          end if
          call nf_put_var(ncdf, var3d(1:ncnt,:),start, count, name='qsfrc')
       enddo
    end if
  end subroutine aq_movie_output
#endif
end module aquaplanet_io_mod
