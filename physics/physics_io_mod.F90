#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module physics_io_mod
  use kinds, only : real_kind
  use common_io_mod, only : nf_selectedvar, nf_handle, nfsizekind
#ifdef PIO_INTERP
  use pio_io_mod, only :  nf_output_register_variables,  &
       nf_put_var, nf_get_frame
#else
  use netcdf_io_mod, only :  nf_output_register_variables, nf_put_var, nf_get_frame 
#endif
  implicit none
  private

  public :: physics_movie_init
  public :: physics_movie_output
  public :: physics_set_varnames
  integer, parameter :: varcnt=6
#ifdef PIO_INTERP
  integer, parameter :: maxdims=3
#else
  integer, parameter :: maxdims=2
#endif
  character*(*), parameter :: varnames(varcnt)=(/'precip ', &
                                                 'tprecip', &
                                                 'CBMF   ', &
                                                 'Wd     ', &
                                                 'Tprime ', &
                                                 'Qprime '/)
contains
  subroutine physics_movie_init(ncdf)
    type(nf_handle), intent(inout) :: ncdf(:)

#ifdef PIO_INTERP
    integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,2,5,  &
                                                                 1,2,5,  &
                                                                 1,2,5,  &
                                                                 1,2,5,  &
                                                                 1,2,5,  &
                                                                 1,2,5/),  &
                                                               shape=(/maxdims,varcnt/))
#else
    integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,4,  &
                                                                 1,4,  &
                                                                 1,4,  &
                                                                 1,4,  &
                                                                 1,4,  &
                                                                 1,4/),  &
                                                               shape=(/maxdims,varcnt/))
#endif 

    call nf_output_register_variables(ncdf,varcnt,varnames,vardims)

  end subroutine physics_movie_init

  subroutine physics_set_varnames(lvarcnt,nlvarnames)
    character*(*), intent(out) :: nlvarnames(:)
    integer, intent(inout) :: lvarcnt
    lvarcnt=lvarcnt+varcnt
    nlvarnames=varnames
  end subroutine physics_set_varnames

#ifdef PIO_INTERP
  subroutine physics_movie_output(ncdf, elem, interpdata, output_varnames, nxyp)
    use dimensions_mod, only : nelemd, nelemdmax, nlev, np
    use dof_mod, only : Uniquepoints
    use physics_types_mod, only : pelem
    use element_mod, only : element_t
    use interpolate_mod, only: interpdata_t, interpolate_scalar

    type(nf_handle), intent(inout) :: ncdf
    type(element_t), intent(in) :: elem(:)
    type(interpdata_t), intent(in) :: interpdata(:)
    character*(*), intent(in) :: output_varnames(:)
    integer, intent(in) :: nxyp
    real(kind=real_kind) :: var2d(nxyp)
    integer :: ie, st, en
    integer(kind=nfsizekind) :: start(3), count(3)

    start(3)=nf_get_frame(ncdf)

    if(nf_selectedvar('precip', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie), pelem(ie)%surfc%precip, &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='precip')
    endif

    if(nf_selectedvar('tprecip', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),pelem(ie)%accum%precip, &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='tprecip')
    endif

    if(nf_selectedvar('CBMF', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),pelem(ie)%state%CBMF, &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp

       enddo
       call nf_put_var(ncdf, var2d, start, count, name='CBMF')
    endif

    if(nf_selectedvar('Wd', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),pelem(ie)%surfc%Wd, &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='Wd')
    endif
    if(nf_selectedvar('Tprime', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie), pelem(ie)%surfc%Tprime, &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='Tprime')
    endif

    if(nf_selectedvar('Qprime', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+interpdata(ie)%n_interp
          call interpolate_scalar(interpdata(ie),pelem(ie)%surfc%Qprime, &
               np, var2d(st:en))
          st=st+interpdata(ie)%n_interp
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='Qprime')
    endif

  end subroutine physics_movie_output

#else

  subroutine physics_movie_output(ncdf, elem, output_varnames, nxyp)
    use dimensions_mod, only : nelemd, nelemdmax, nlev
    use dof_mod, only : Uniquepoints
    use physics_types_mod, only : pelem
    use element_mod, only : element_t
    type(nf_handle), intent(inout) :: ncdf
    type(element_t), intent(in) :: elem(:)
    character*(*), intent(in) :: output_varnames(:)
    integer, intent(in) :: nxyp
    real(kind=real_kind) :: var2d(nxyp)
    integer :: ie, st, en
    integer(kind=nfsizekind) :: start(2), count(2)

    start(2)=nf_get_frame(ncdf)

    if(nf_selectedvar('precip', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,pelem(ie)%surfc%precip,var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='precip')
    endif

    if(nf_selectedvar('tprecip', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,pelem(ie)%accum%precip,var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='tprecip')
    endif

    if(nf_selectedvar('CBMF', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,pelem(ie)%state%CBMF,var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='CBMF')
    endif

    if(nf_selectedvar('Wd', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,pelem(ie)%surfc%Wd,var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='Wd')
    endif
    if(nf_selectedvar('Tprime', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,pelem(ie)%surfc%Tprime,var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='Tprime')
    endif

    if(nf_selectedvar('Qprime', output_varnames)) then
       st=1
       do ie=1,nelemd
          en=st+elem(ie)%idxp%NumUniquePts-1
          call UniquePoints(elem(ie)%idxP,pelem(ie)%surfc%Qprime,var2D(st:en))
          st=en+1
       enddo
       call nf_put_var(ncdf, var2d, start, count, name='Qprime')
    endif

  end subroutine physics_movie_output
#endif

#ifdef XXXXXX
! obsolete. HOMME only supports PIO or PIO_INTERP   defined(NETCDF) || defined(PNETCDF)
  subroutine physics_movie_output(ncdf, elem, output_varnames, nxyp)
    use dimensions_mod, only : nelemd, nelemdmax, nlev
    use dof_mod, only : Uniquepoints
    use physics_types_mod, only : pelem
    use element_mod, only : element_t
    type(nf_handle), intent(inout) :: ncdf
    type(element_t), intent(in) :: elem(:)
    character*(*), intent(in) :: output_varnames(:)
    integer, intent(in) :: nxyp
    real(kind=real_kind) :: var2d(nxyp)
    integer :: ie
    integer(kind=nfsizekind) :: start(2), count(2)

    start(2)=nf_get_frame(ncdf)

    if(nf_selectedvar('precip', output_varnames)) then
       start(1)=1
       count(1)=0
       count(2)=1
       do ie=1,nelemdmax
          if(ie<=nelemd) then
             start(1)=elem(ie)%idxp%UniquePtOffset
             count(1)=elem(ie)%idxp%numUniquePts
             call UniquePoints(elem(ie)%idxP,pelem(ie)%surfc%precip,var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d, start, count, name='precip')
       enddo
    endif
    if(nf_selectedvar('tprecip', output_varnames)) then
       start(1)=1
       count(1)=0
       count(2)=1
       do ie=1,nelemdmax
          if(ie<=nelemd) then
             start(1)=elem(ie)%idxp%UniquePtOffset
             count(1)=elem(ie)%idxp%numUniquePts
             call UniquePoints(elem(ie)%idxP,pelem(ie)%accum%precip,var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d,start, count,name='tprecip')
       enddo
    endif
    if(nf_selectedvar('CBMF', output_varnames)) then
       start(1)=1
       count(1)=0
       count(2)=1
       do ie=1,nelemdmax
          if(ie<=nelemd) then
             start(1)=elem(ie)%idxp%UniquePtOffset
             count(1)=elem(ie)%idxp%numUniquePts
             call UniquePoints(elem(ie)%idxP,pelem(ie)%state%CBMF,var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d,start, count,name='CBMF')
       enddo
    endif
    if(nf_selectedvar('Wd', output_varnames)) then
       start(1)=1
       count(1)=0
       count(2)=1
       do ie=1,nelemdmax
          if(ie<=nelemd) then
             start(1)=elem(ie)%idxp%UniquePtOffset
             count(1)=elem(ie)%idxp%numUniquePts
             call UniquePoints(elem(ie)%idxP,pelem(ie)%surfc%Wd,var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d,start, count, name='Wd')
       enddo
    endif
    if(nf_selectedvar('Tprime', output_varnames)) then
       start(1)=1
       count(1)=0
       count(2)=1
       do ie=1,nelemdmax
          if(ie<=nelemd) then
             start(1)=elem(ie)%idxp%UniquePtOffset
             count(1)=elem(ie)%idxp%numUniquePts
             call UniquePoints(elem(ie)%idxP,pelem(ie)%surfc%Tprime,var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d,start, count, name='Tprime')
       enddo
    endif
    if(nf_selectedvar('Qprime', output_varnames)) then
       start(1)=1
       count(1)=0
       count(2)=1
       do ie=1,nelemdmax
          if(ie<=nelemd) then
             start(1)=elem(ie)%idxp%UniquePtOffset
             count(1)=elem(ie)%idxp%numUniquePts
             call UniquePoints(elem(ie)%idxP,pelem(ie)%surfc%Qprime,var2D)
          else
             count=0
          end if
          call nf_put_var(ncdf, var2d,start, count, name='Qprime')
       enddo

    endif
  end subroutine physics_movie_output
#endif
end module physics_io_mod
