
module damping_mod
  use params, only: asyncid
  use task_util_mod
  implicit none

contains

  subroutine damping(ncrms)
    !  "Spange"-layer damping at the domain top region
    use vars
    use microphysics, only: micro_field, index_water_vapor
    use params, only: crm_rknd
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) tau_min    ! minimum damping time-scale (at the top)
    real(crm_rknd) tau_max    ! maxim damping time-scale (base of damping layer)
    real(crm_rknd) damp_depth ! damping depth as a fraction of the domain height
    parameter(tau_min=60.D0, tau_max=450.D0, damp_depth=0.4D0)
    real(crm_rknd) tau(ncrms,nzm), tmp
    integer, allocatable :: n_damp(:)
    integer :: i, j, k, icrm
    integer :: numgangs  !For working around PGI OpenACC bug where it didn't create enough gangs
    ! crjones tests: make changes to u0, v0, t0 local instead of shared with vars
    real(crm_rknd), allocatable :: t0loc(:,:)
    real(crm_rknd), allocatable :: u0loc(:,:)
    real(crm_rknd), allocatable :: v0loc(:,:)

    allocate( n_damp(ncrms) )
    allocate( t0loc(ncrms,nzm) )
    allocate( u0loc(ncrms,nzm) )
    allocate( v0loc(ncrms,nzm) )
    call prefetch( n_damp)
    call prefetch( t0loc )
    call prefetch( u0loc )
    call prefetch( v0loc )
   
    if(tau_min.lt.2*dt) then
      print*,'Error: in damping() tau_min is too small!'
      call task_abort()
    end if

    !$acc parallel loop async(asyncid)
    do icrm = 1 , ncrms
      do k=nzm,1,-1
        if(z(icrm,nzm)-z(icrm,k).lt.damp_depth*z(icrm,nzm)) then
          n_damp(icrm)=nzm-k+1
        endif
      end do
    end do

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nzm
      do icrm = 1 , ncrms
        if ( k <= nzm .and. k >= nzm-n_damp(icrm) ) then
          tau(icrm,k) = tau_min *(tau_max/tau_min)**((z(icrm,nzm)-z(icrm,k))/(z(icrm,nzm)-z(icrm,nzm-n_damp(icrm))))
          tau(icrm,k)=1.D0/tau(icrm,k)
        endif
      end do
    end do

    ! recalculate grid-mean u0, v0, t0 first,
    ! as t has been updated. No need for qv0, as
    ! qv has not been updated yet the calculation of qv0.
    !$acc parallel loop collapse(2) async(asyncid)
    do k=1, nzm
      do icrm = 1 , ncrms
        u0loc(icrm,k)=0.0
        v0loc(icrm,k)=0.0
        t0loc(icrm,k)=0.0
      end do
    end do
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1, nzm
      do j=1, ny
        do i=1, nx
          do icrm = 1 , ncrms
            tmp = u(icrm,i,j,k)/(nx*ny)
            !$acc atomic update
            u0loc(icrm,k) = u0loc(icrm,k) + tmp
            tmp = v(icrm,i,j,k)/(nx*ny)
            !$acc atomic update
            v0loc(icrm,k) = v0loc(icrm,k) + tmp
            tmp = t(icrm,i,j,k)/(nx*ny)
            !$acc atomic update
            t0loc(icrm,k) = t0loc(icrm,k) + tmp
          end do
        end do
      end do
    end do

   !For working around PGI OpenACC bug where it didn't create enough gangs 
    numgangs = ceiling(ncrms*ny*nx/128.D0)
    !$acc parallel loop collapse(4) vector_length(128) num_gangs(numgangs) async(asyncid)
    do k = 1 , nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            if ( k <= nzm .and. k >= nzm-n_damp(icrm) ) then
              dudt(icrm,i,j,k,na)= dudt(icrm,i,j,k,na)-(u(icrm,i,j,k)-u0loc(icrm,k)) * tau(icrm,k)
              dvdt(icrm,i,j,k,na)= dvdt(icrm,i,j,k,na)-(v(icrm,i,j,k)-v0loc(icrm,k)) * tau(icrm,k)
              dwdt(icrm,i,j,k,na)= dwdt(icrm,i,j,k,na)-w(icrm,i,j,k) * tau(icrm,k)
              t(icrm,i,j,k)= t(icrm,i,j,k)-dtn*(t(icrm,i,j,k)-t0loc(icrm,k)) * tau(icrm,k)
              micro_field(icrm,i,j,k,index_water_vapor)= micro_field(icrm,i,j,k,index_water_vapor)-dtn*(qv(icrm,i,j,k)-qv0(icrm,k)) * tau(icrm,k)
            endif
          end do! i
        end do! j
      end do ! k
    end do

    deallocate( n_damp )
    deallocate( t0loc )
    deallocate( u0loc )
    deallocate( v0loc )

  end subroutine damping

end module damping_mod
