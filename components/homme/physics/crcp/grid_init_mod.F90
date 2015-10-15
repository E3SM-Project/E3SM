module grid_init_mod
  real :: dt,dti,dx,dxi,dz,dzi
  real :: time
  ! exponent of stretching function (set in namelist):   
  real :: ex1

  real, allocatable :: height(:), gac(:), xx(:), zz(:)
contains
  subroutine grid_init(dtin, dxin, dzin, nx,nz)
    implicit none
    real, intent(in) :: dtin, dxin, dzin
    integer, intent(in) :: nx,nz
  
    integer :: i, k


    dt = dtin
    dti = 1./dt
    dx = dxin
    dxi= 1./dx
    dz = dzin
    dzi= 1./dz

    if(.not.allocated(height)) then
       allocate(height(nz), gac(nz), xx(nx), zz(nz))
    end if
    time = 0.0

    do i=1,nx
       xx(i)=float(i-1)*dx
    enddo
    ! zz is regular (ustretched) dzeta grid
    do k=1,nz
       zz(k)=float(k-1)*dz
    enddo

    call zstrtch(zz,nz,dz)

  end subroutine grid_init

  SUBROUTINE zstrtch(zz, nz, dz) 
    implicit none
    !c define vertical grid for stretched coordinate, derive Jacobian       
    integer, intent(in) :: nz
    real, intent(in) :: zz(nz), dz 

    real :: zb(100) 
    real :: sum1 , coe, top, aa
                                  
    integer :: i, k

    IF(nz.gt.100) stop 'grid formulation' 

    top = float(nz - 1) * dz 


    !    ex1 = 4.1 / 2.9 
    !      ex1=1.0                                                          

    aa = top / top**ex1 

    !c vertical coordinate:                                                 
    DO k = 1, nz 
       zb(k) = aa * zz(k) **ex1 
       height(k) = zb(k) 
    enddo

    !c jacobian:  
    gac(1) =(zb(2) - zb(1) ) *dzi   
    gac(nz) =(zb(nz) - zb(nz - 1) ) *dzi 
    DO k = 2, nz-1 
!       IF(k.eq.1) gac(k) =(zb(2) - zb(1) ) / dz 
!       IF(k.eq.nz) gac(k) =(zb(nz) - zb(nz - 1) ) / dz 
!       IF(k.ne.1.and.k.ne.nz) gac(k) =(zb(k + 1) - zb(k - 1) )      &
!            /(2. * dz)                                                       
       gac(k) =(zb(k + 1) - zb(k - 1) )*dzi*0.5
    end do
    do k=1,nz
       PRINT *, k, zb(k), gac(k) 
    enddo

    !c check consistency(int gac ddzeta = H)                               
    sum1 = .5 *(gac(1) * dz + gac(nz) * dz) 
    DO i = 2, nz - 1 
       sum1 = sum1 + gac(i) * dz 
    enddo
    PRINT * , 'int Jacobian before adj: ', sum1 

    !c adjust numerical jacobian:                                           
    coe = float(nz - 1) * dz / sum1 
    DO i = 1, nz 
       gac(i) = gac(i) * coe 
    enddo

    !c check:                                                               
    sum1 = .5 *(gac(1) * dz + gac(nz) * dz) 
    DO i = 2, nz - 1 
       sum1 = sum1 + gac(i) * dz 
    enddo
    PRINT * , 'int Jacobian after adj: ', sum1

    RETURN 
  END SUBROUTINE zstrtch
end module grid_init_mod
