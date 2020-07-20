
module forcing_mod
  implicit none

contains

  subroutine forcing(ncrms)
    use vars
    use params
    use microphysics, only: micro_field, index_water_vapor
#if defined(_OPENACC)
    use openacc_utils
#endif
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
#if defined(_OPENACC)
    call prefetch( qneg )
    call prefetch( qpoz )
    call prefetch( nneg )
#elif defined(_OPENMP)
    !$omp target enter data map(alloc: qneg)
    !$omp target enter data map(alloc: qpoz)
    !$omp target enter data map(alloc: nneg)
#endif

    coef = 1./3600.
#if defined(_OPENACC)
    !$acc parallel loop collapse(2) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(2)
#endif
    do k=1,nzm
      do icrm = 1 , ncrms
        qpoz(icrm,k) = 0.
        qneg(icrm,k) = 0.
        nneg(icrm,k) = 0
      enddo
    enddo
#if defined(_OPENACC)
    !$acc parallel loop collapse(4) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(4)
#endif
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif
            t(icrm,i,j,k)=t(icrm,i,j,k) + ttend(icrm,k) * dtn
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif
            micro_field(icrm,i,j,k,index_water_vapor)=micro_field(icrm,i,j,k,index_water_vapor) + qtend(icrm,k) * dtn
            if(micro_field(icrm,i,j,k,index_water_vapor).lt.0.) then
#if defined(_OPENACC)
              !$acc atomic update
#elif defined(_OPENMP)
              !$omp atomic update
#endif
              nneg(icrm,k) = nneg(icrm,k) + 1
#if defined(_OPENACC)
              !$acc atomic update
#elif defined(_OPENMP)
              !$omp atomic update
#endif
              qneg(icrm,k) = qneg(icrm,k) + micro_field(icrm,i,j,k,index_water_vapor)
            else
#if defined(_OPENACC)
              !$acc atomic update
#elif defined(_OPENMP)
              !$omp atomic update
#endif
              qpoz(icrm,k) = qpoz(icrm,k) + micro_field(icrm,i,j,k,index_water_vapor)
            end if
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na) + utend(icrm,k)
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na) + vtend(icrm,k)
          end do
        end do
      end do
    end do
#if defined(_OPENACC)
    !$acc parallel loop collapse(4) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(4)
#endif
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
#if defined(_OPENMP)
    !$omp target exit data map(delete: qneg)
    !$omp target exit data map(delete: qpoz)
    !$omp target exit data map(delete: nneg)
#endif
    deallocate( qneg )
    deallocate( qpoz )
    deallocate( nneg )

  end subroutine forcing

end module forcing_mod
