!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  diag_dynvar_ic --- record state variables to IC file
!
! !INTERFACE:
  subroutine diag_dynvar_ic(grid, phis, ps, t3, u3s, v3s, tracer)

! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    use cam_history  , only: outfld, write_inithist
    use constituents , only: cnst_name
    use dynamics_vars, only: T_FVDYCORE_GRID

    implicit none
!
!-----------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    type (T_FVDYCORE_GRID), intent(in) :: grid
! surface geopotential (grav*zs)
    real(r8), intent(in) :: phis(grid%ifirstxy:grid%ilastxy,            &
                                 grid%jfirstxy:grid%jlastxy)
! Surface pressure (pa)   
    real(r8), intent(in) :: ps  (grid%ifirstxy:grid%ilastxy,            &
                                 grid%jfirstxy:grid%jlastxy)
! Temperature (K)
    real(r8), intent(in) :: t3  (grid%ifirstxy:grid%ilastxy,            &
                                 grid%jfirstxy:grid%jlastxy, grid%km )  
! u wind velocities, staggered grid
    real(r8), intent(in) :: u3s (grid%ifirstxy:grid%ilastxy,            &
                                 grid%jfirstxy:grid%jlastxy, grid%km)
! v wind velocities, staggered grid
    real(r8), intent(in) :: v3s (grid%ifirstxy:grid%ilastxy,            &
                                 grid%jfirstxy:grid%jlastxy, grid%km)
! Tracers
    real(r8), intent(in) :: tracer (grid%ifirstxy:grid%ilastxy,          &
                                    grid%jfirstxy:grid%jlastxy,          &
                                    grid%km,grid%ntotq)

! !HISTORY:
!   01.01.01   XXXXXX  Delivery
!   05.07.12   Sawyer  Simplified interface with grid argument, ProTeX
!   06.03.22   Sawyer  Rewritten for XY decomposition
!   06.06.28   Sawyer  T3 changed from IKJ to IJK indexing
!   06.07.01   Sawyer  Transitioned tracers q3 to T_TRACERS
!
!EOP
!-----------------------------------------------------------------------
!BOC
!---------------------------Local workspace-----------------------------
    integer :: i, j, k, m   ! indices
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy, km, ntotq, idim
    real(r8):: tmp(grid%ifirstxy:grid%ilastxy,grid%km)
!
!-----------------------------------------------------------------------
!

    ifirstxy     = grid%ifirstxy
    ilastxy      = grid%ilastxy
    jfirstxy     = grid%jfirstxy
    jlastxy      = grid%jlastxy
    km           = grid%km
    ntotq        = grid%ntotq
    idim         = ilastxy - ifirstxy + 1

    if( write_inithist() ) then

!$OMP PARALLEL DO PRIVATE (I, J, K, M, TMP)
       do j = jfirstxy, jlastxy

          call outfld ('PS&IC      ', ps  (:,j) , idim, j)

          do k = 1, km
             do i = ifirstxy, ilastxy
                tmp(i,k) = t3(i,j,k)
             enddo
          enddo
          call outfld ('T&IC       ', tmp       , idim, j) 

          do k = 1, km
             do i = ifirstxy, ilastxy
                tmp(i,k) = u3s(i,j,k)
             enddo
          enddo
          call outfld ('US&IC      ', tmp       , idim, j)

          do k = 1, km
             do i = ifirstxy, ilastxy
                tmp(i,k) = v3s(i,j,k)
             enddo
          enddo
          call outfld ('VS&IC      ', tmp       , idim, j)

          do m = 1, ntotq
             do k = 1, km
                do i = ifirstxy, ilastxy
                   tmp(i,k) = tracer(i,j,k,m)
                enddo
             enddo
             call outfld(trim(cnst_name(m))//'&IC' , tmp  , idim, j)
          end do

       enddo

    end if

    return
!EOC
  end subroutine diag_dynvar_ic
