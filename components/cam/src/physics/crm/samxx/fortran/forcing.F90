
module forcing_mod
  implicit none

contains

  subroutine forcing(ncrms)
    use vars
    use params
    use microphysics, only: micro_field, index_water_vapor
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd), allocatable :: qneg(:,:)
    real(crm_rknd), allocatable :: qpoz(:,:)
    integer       , allocatable :: nneg(:,:)
    real(crm_rknd) :: coef, factor
    integer        :: i,j,k,icrm

    allocate( qneg(ncrms,nzm) )
    allocate( qpoz(ncrms,nzm) )
    allocate( nneg(ncrms,nzm) )
    call prefetch( qneg )
    call prefetch( qpoz )
    call prefetch( nneg )

    coef = 1./3600.

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nzm
      do icrm = 1 , ncrms
        qpoz(icrm,k) = 0.
        qneg(icrm,k) = 0.
        nneg(icrm,k) = 0
      enddo
    enddo

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            t(icrm,i,j,k)=t(icrm,i,j,k) + ttend(icrm,k) * dtn
            micro_field(icrm,i,j,k,index_water_vapor)=micro_field(icrm,i,j,k,index_water_vapor) + qtend(icrm,k) * dtn
            if(micro_field(icrm,i,j,k,index_water_vapor).lt.0.) then
              !$acc atomic update
              nneg(icrm,k) = nneg(icrm,k) + 1
              !$acc atomic update
              qneg(icrm,k) = qneg(icrm,k) + micro_field(icrm,i,j,k,index_water_vapor)
            else
              !$acc atomic update
              qpoz(icrm,k) = qpoz(icrm,k) + micro_field(icrm,i,j,k,index_water_vapor)
            end if
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na) + utend(icrm,k)
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na) + vtend(icrm,k)
          end do
        end do
      end do
    end do

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            if(nneg(icrm,k).gt.0.and.qpoz(icrm,k)+qneg(icrm,k).gt.0.) then
              factor = 1. + qneg(icrm,k)/qpoz(icrm,k)
              micro_field(icrm,i,j,k,index_water_vapor) = max(real(0.,crm_rknd),micro_field(icrm,i,j,k,index_water_vapor)*factor)
            end if
          end do
        end do
      end do
    end do

    deallocate( qneg )
    deallocate( qpoz )
    deallocate( nneg )

  end subroutine forcing

end module forcing_mod
