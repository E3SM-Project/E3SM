module Gauss_module

  use Grid_Unstructured_Cell_module
  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"

  private
  
  PetscInt, parameter, public :: LINE_TYPE          = 7
  
  type, public :: gauss_type
    PetscInt :: dim                 ! dimension
    PetscInt :: EleType             ! Element type
    PetscInt :: NGPTS               ! Number of gauss points
    PetscReal, pointer :: r(:,:)    ! location of points
    PetscReal, pointer :: w(:)      ! weights
  end type gauss_type
    
  public :: GaussCalculatePoints, GaussDestroy, GaussInitialize
  
  contains

! ************************************************************************** !

subroutine GaussInitialize(gauss)
  ! 
  ! Initializes Gauss type
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 6/19/2013
  ! 

  type(gauss_type) :: gauss
  
  gauss%dim = 0 
  gauss%EleType = 0
  gauss%NGPTS = 0
  nullify(gauss%r)
  nullify(gauss%w)

end subroutine GaussInitialize   

! ************************************************************************** !

subroutine GaussCalculatePoints(gauss)
  ! 
  ! Calculates Gauss points
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  use Utility_module, only: DeallocateArray

  type(gauss_type) :: gauss
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  
  select case(gauss%dim)
    case(ONE_DIM_GRID)
      call Gauss1D(gauss%EleType,gauss%NGPTS,r,w)
    case(TWO_DIM_GRID)
      call Gauss2D(gauss%EleType,gauss%NGPTS,r,w)
    case(THREE_DIM_GRID)
      call Gauss3D(gauss%EleType,gauss%NGPTS,r,w)
    case default
      print *, 'Error: Invalid dimension for Gauss point calculation'
      stop
  end select  
  
  allocate(gauss%r(size(r,1),size(r,2)))
  allocate(gauss%w(size(w)))
  
  gauss%r = r
  gauss%w = w
  
  deallocate(r)
  deallocate(w)
    
end subroutine GaussCalculatePoints  

! ************************************************************************** !

subroutine Gauss1D(EleType,NGPTS,r,w)
  ! 
  ! Calculates Gauss points for 1D elements
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  PetscInt :: EleType
  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  
  allocate(r(NGPTS,1))
  allocate(w(NGPTS))
  
  if (EleType /= LINE_TYPE) then
    print *, 'Error: in Element type. Only L2 ' // &
             '(line type) can be used for 1D Gauss quadrature.'
  endif
    
  select case(NGPTS)
    case(1)

    !------------------------------
    ! No of Gauss Points = 1
    !------------------------------
        
      r(1,1) = 0.0
      !
      w(1) = 2.0

    !------------------------------
    ! No of Gauss Points = 2
    !------------------------------

    case(2)
        
      r(1,1) = -1.0/sqrt(3.0)
      r(2,1) = -r(1,1) 
      !
      w(1) = 1.0
      w(2) = 1.0

    !-------------------------------
    ! No of Gauss Points = 3
    !-------------------------------

    case(3)
        
      r(1,1) = -sqrt(0.6)
      r(2,1) = 0.0
      r(3,1) = -r(1,1)
      !
      w(1) = 5.0/9.0
      w(2) = 8.0/9.0
      w(3) = w(1)

    !--------------------------------
    ! No of Gauss Points = 4
    !--------------------------------

    case(4)
        
      r(1,1) = -0.861136311594053
      r(2,1) = -0.339981043584856
      r(3,1) =  0.339981043584856
      r(4,1) =  0.861136311594053
      ! 
      w(1) = 0.347854845137454
      w(2) = 0.652145154862546
      w(3) = 0.652145154862546
      w(4) = 0.347854845137454

    !----------------------------------
    ! No of Gauss Points = 5
    !----------------------------------

    case(5)
        
      r(1,1) = -0.906179845938664
      r(2,1) = -0.538469310105683
      r(3,1) =  0.000000000000000
      r(4,1) =  0.538469310105683
      r(5,1) =  0.906179845938664
      !
      w(1) =  0.236926885056189
      w(2) =  0.478628670499366 
      w(3) =  0.568888888888889
      w(4) =  0.478628670499366
      w(5) =  0.236926885056189

    !----------------------------------
    ! No of Gauss Points = 6
    !----------------------------------
   
    case(6)
        
      r(1,1) = -0.932469514203152
      r(2,1) = -0.661209386466265
      r(3,1) = -0.238619186083197
      r(4,1) =  0.238619186083197
      r(5,1) =  0.661209386466265
      r(6,1) =  0.932469514203152
      !
      w(1) =  0.171324492379170
      w(2) =  0.360761573048139
      w(3) =  0.467913934572691
      w(4) =  0.467913934572691
      w(5) =  0.360761573048139
      w(6) =  0.171324492379170

    !------------------------------------
    ! No of Gauss Points = 7
    !------------------------------------

    case(7)
        
      r(1,1) = -0.949107912342759
      r(2,1) = -0.741531185599394
      r(3,1) = -0.405845151377397
      r(4,1) =  0.000000000000000
      r(5,1) =  0.405845151377397
      r(6,1) =  0.741531185599394
      r(7,1) =  0.949107912342759
      !
      w(1) =  0.129484966168870
      w(2) =  0.279705391489277
      w(3) =  0.381830050505119
      w(4) =  0.417959183673469
      w(5) =  0.381830050505119
      w(6) =  0.279705391489277
      w(7) =  0.129484966168870

    !------------------------------------
    ! No of Gauss Points = 8
    !------------------------------------
    
    case(8)
        
      r(1,1) = -0.960289856497536
      r(2,1) = -0.796666477413627
      r(3,1) = -0.525532409916329
      r(4,1) = -0.183434642495650
      r(5,1) =  0.183434642495650
      r(6,1) =  0.525532409916329
      r(7,1) =  0.796666477413627
      r(8,1) =  0.960289856497536
      !
      w(1) =  0.101228536290376
      w(2) =  0.222381034453374
      w(3) =  0.313706645877887
      w(4) =  0.362683783378362
      w(5) =  0.362683783378362
      w(6) =  0.313706645877887
      w(7) =  0.222381034453374
      w(8) =  0.101228536290376

    !------------------------------------
    ! No of Gauss Points = 9
    !------------------------------------
 
    case(9)
        
      r(1,1) = -0.968160239507626
      r(2,1) = -0.836031170326636
      r(3,1) = -0.613371432700590
      r(4,1) = -0.324253423403809
      r(5,1) =  0.000000000000000
      r(6,1) =  0.324253423403809
      r(7,1) =  0.613371432700590
      r(8,1) =  0.836031107326636
      r(9,1) =  0.968160239507626

      w(1) =  0.081274388361574
      w(2) =  0.180648160694857
      w(3) =  0.260610696402935
      w(4) =  0.312347077040003
      w(5) =  0.330239355001260
      w(6) =  0.312347077040003
      w(7) =  0.260610696402935
      w(8) =  0.180648160694857
      w(9) =  0.081274388361574

   case default
     print *, 'Error in NGPTS for 1D Gauss quadrature'
     stop
   end select

end subroutine Gauss1D  

! ************************************************************************** !

subroutine Gauss2D(EleType,NGPTS,r,w)
  ! 
  ! Calculates Gauss points for 2D elements
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  PetscInt :: EleType
  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
    
  select case(EleType)
    case(QUAD_TYPE)
      call GaussSquare(NGPTS,r,w)
    case(TRI_TYPE)
      call GaussTriangle(NGPTS,r,w)
    case default
      print *, 'Error: Only T3 and Q4 elements available for 2D.'
      stop
  end select

end subroutine Gauss2D

! ************************************************************************** !

subroutine GaussSquare(NGPTS,r,w)
  ! 
  ! Calculates Gauss points for Q4 element
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal, pointer :: l(:,:)
  PetscReal, pointer :: m(:)
  PetscInt :: counter,i,j 
  
  allocate(r(NGPTS*NGPTS,2))
  allocate(w(NGPTS*NGPTS))
  allocate(l(NGPTS,1))
  allocate(m(NGPTS))

  ! 1D Gauss points are stored in l vector and weights are stored in m vector

  call Gauss1D(LINE_TYPE,NGPTS,l,m)
  
  ! Generate the Q4 Gauss points and weights using for loops
  counter = 1
  do i = 1, NGPTS
    do j = 1, NGPTS
      r(counter,1) = l(i,1)
      r(counter,2) = l(j,1)
      w(counter) = m(i)*m(j)
      counter = counter + 1
    enddo
  enddo

  deallocate(l)
  deallocate(m)

end subroutine GaussSquare

! ************************************************************************** !

subroutine GaussTriangle(NGPTS,r,w)
  ! 
  ! Calculates Gauss points for T3 elements
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  PetscInt :: EleType
  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  
  allocate(r(NGPTS,2))
  allocate(w(NGPTS))
    
  select case(NGPTS)
    case(1)

    !------------------------------
    ! No of Gauss Points = 1
    !------------------------------
      r(1,:) = 1.0/3.0*(/1.0,1.0/)
      !
      w(1) = 0.5

    !-------------------------------
    ! No of Gauss Points = 3
    !-------------------------------

    case(3)
        
      r(1,:) = 0.5*(/1.0,1.0/)
      r(2,:) = 0.5*(/1.0,0.0/)
      r(3,:) = 0.5*(/0.0,1.0/) 
      !
      w = 1.0/6.0*(/1.0,1.0,1.0/)

    !--------------------------------
    ! No of Gauss Points = 4
    !--------------------------------

    case(4)
    
      r(1,:) = (/0.333333333333333,0.333333333333333/)
      r(2,:) = (/0.600000000000000,0.200000000000000/)
      r(3,:) = (/0.200000000000000,0.600000000000000/)
      r(4,:) = (/0.200000000000000,0.200000000000000/)
      ! 
      w = 0.5*(/-0.562500000000000, &
                 0.520833333333333, &
                 0.520833333333333, &
                 0.520833333333333/)

    !------------------------------------
    ! No of Gauss Points = 7
    !------------------------------------

    case(7)
    
      r(1,:) = (/0.333333333333333,0.333333333333333/)
      r(2,:) = (/0.797426985353087,0.101286507323456/)
      r(3,:) = (/0.101286507323456,0.797426985353087/)
      r(4,:) = (/0.101286507323456,0.101286507323456/)
      r(5,:) = (/0.470142064105115,0.059715871789770/)
      r(6,:) = (/0.059715871789770,0.470142064105115/)
      r(7,:) = (/0.470142064105115,0.470142064105115/)
      !
      w = 0.5*(/0.225000000000000, &
                0.125939180544827, &
                0.125939180544827, &
                0.125939180544827, &
                0.132394152788506, &
                0.132394152788506, &
                0.132394152788506/)

    !------------------------------------
    ! No of Gauss Points = 9
    !------------------------------------
 
    case(9)
       
      r(1,:) = (/0.124949503233232,0.437525248383384/)
      r(2,:) = (/0.437525248383384,0.124949503233232/)
      r(3,:) = (/0.437525248383384,0.437525248383384/)
      r(4,:) = (/0.797112651860071,0.165409927389841/)
      r(5,:) = (/0.797112651860071,0.037477420750088/)
      r(6,:) = (/0.165409927389841,0.797112651860071/)
      r(7,:) = (/0.165409927389841,0.037477420750088/)
      r(8,:) = (/0.037477420750088,0.797112651860071/)
      r(9,:) = (/0.037477420750088,0.165409927389841/)

      w = 0.5*(/0.205950504760887, &
                0.205950504760887, &
                0.205950504760887, &
                0.063691414286223, &
                0.063691414286223, &
                0.063691414286223, &
                0.063691414286223, &
                0.063691414286223, &
                0.063691414286223/)
                
    !------------------------------------
    ! No of Gauss Points = 12
    !------------------------------------
 
    case(12)
       
      r(1,:) = (/0.873821971016996,0.063089014491502/)
      r(2,:) = (/0.063089014491502,0.873821971016996/) 
      r(3,:) = (/0.063089014491502,0.063089014491502/)
      r(4,:) = (/0.501426509658179,0.249286745170910/)
      r(5,:) = (/0.249286745170910,0.501426509658179/) 
      r(6,:) = (/0.249286745170910,0.249286745170910/)
      r(7,:) = (/0.636502499121399,0.310352451033785/)
      r(8,:) = (/0.636502499121399,0.053145049844816/)
      r(9,:) = (/0.310352451033785,0.636502499121399/)
      r(10,:) = (/0.310352451033785,0.053145049844816/)
      r(11,:) = (/0.053145049844816,0.636502499121399/) 
      r(12,:) = (/0.053145049844816,0.310352451033785/)

      w = 0.5*(/0.050844906370207, &
                0.050844906370207, &
                0.050844906370207, &
                0.116786275726379, &
                0.116786275726379, &
                0.116786275726379, &
                0.082851075618374, &
                0.082851075618374, &
                0.082851075618374, &
                0.082851075618374, &
                0.082851075618374, &
                0.082851075618374/)
                
    !------------------------------------
    ! No of Gauss Points = 13
    !------------------------------------
 
    case(13)
       
      r(1,:) = (/0.333333333333333,0.333333333333333/)
      r(2,:) = (/0.479308067841923,0.260345966079038/)
      r(3,:) = (/0.260345966079038,0.479308067841923/)
      r(4,:) = (/0.260345966079038,0.260345966079038/)
      r(5,:) = (/0.869739794195568,0.065130102902216/)
      r(6,:) = (/0.065130102902216,0.869739794195568/)
      r(7,:) = (/0.065130102902216,0.065130102902216/)
      r(8,:) = (/0.638444188569809,0.312865496004875/)
      r(9,:) = (/0.638444188569809,0.086903154253160/)
      r(10,:) = (/0.312865496004875,0.638444188569809/)
      r(11,:) = (/0.312865496004875,0.086903154253160/)
      r(12,:) = (/0.086903154253160,0.638444188569809/)
      r(13,:) = (/0.086903154253160,0.312865496004875/) 
      
      w = 0.5*(/-0.149570044467670, &
                +0.175615257433204, &
                +0.175615257433204, &
                +0.175615257433204, &
                +0.053347235608839, &
                +0.053347235608839, &
                +0.053347235608839, &
                +0.077113760890257, &
                +0.077113760890257, &
                +0.077113760890257, &
                +0.077113760890257, &
                +0.077113760890257, &
                +0.077113760890257/)
              
   case default
     print *, 'Invalid NGPTS for T3 Gauss quadrature'
     stop
   end select
  
  
end subroutine GaussTriangle

! ************************************************************************** !

subroutine GaussTetrahedra(NGPTS,r,w)
  ! 
  ! Calculates Gauss points for tetrahedra elements
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 7/11/2013
  ! 

  PetscInt :: EleType
  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal :: p1,p2,p3,p4,p5,p6,p7
  PetscReal :: q1,q2,q3,q4
  
  allocate(r(NGPTS,3))
  allocate(w(NGPTS))
    
  select case(NGPTS)
    case(1)

    !------------------------------
    ! No of Gauss Points = 1
    !------------------------------
      r(1,:) = 1.0/4.0*(/1.0,1.0,1.0/)
      !
      w(1) = 1.0/6.0

    !--------------------------------
    ! No of Gauss Points = 4
    !--------------------------------

    case(4)
    
      p1 = 0.5854101966249638
      p2 = 0.1381966011250105
      r(1,:) = (/p1,p2,p2/)
      r(2,:) = (/p2,p1,p2/)
      r(3,:) = (/p2,p2,p1/)
      r(4,:) = (/p1,p1,p1/)
      ! 
      w = 1.0/(6.0*4.0)*(/1.0,1.0,1.0,1.0/)
      
    !------------------------------------
    ! No of Gauss Points = 5
    !------------------------------------     
    
    case(5)
    
      r(1,:) = 1.0/4.0*(/1.0,1.0,1.0/)
      r(2,:) = (/1.0/2.0,1.0/6.0,1.0/6.0/)
      r(3,:) = (/1.0/6.0,1.0/2.0,1.0/6.0/)
      r(4,:) = (/1.0/6.0,1.0/6.0,1.0/2.0/)
      r(5,:) = (/1.0/6.0,1.0/6.0,1.0/6.0/)
      !
      w = 1.0/6.0*(/-4.0/5.0,9.0/20.0,9.0/20.0,9.0/20.0,9.0/20.0/)

    !------------------------------------
    ! No of Gauss Points = 11
    !------------------------------------

    case(11)
    
      p1 = 0.250000000000000  
      p2 = 0.785714285714286
      p3 = 0.071428571428571
      p4 = 0.399403576166799
      p5 = 0.100596423833201

      r(1,:) = (/p1,p1,p1/)
      r(2,:) = (/p2,p3,p3/)
      r(3,:) = (/p3,p2,p3/)
      r(4,:) = (/p3,p3,p2/)
      r(5,:) = (/p3,p3,p3/)
      r(6,:)  = (/p4,p5,p5/)
      r(7,:)  = (/p5,p4,p5/)
      r(8,:)  = (/p5,p5,p4/)
      r(9,:)  = (/p5,p4,p4/)
      r(10,:) = (/p4,p5,p4/)
      r(11,:) = (/p4,p4,p5/)
      !
      q1 = -0.013155555555556
      q2 =  0.007622222222222
      q3 =  0.024888888888889

      w = (/q1,q2,q2,q2,q2,q3,q3,q3,q3,q3,q3/)
      

    !------------------------------------
    ! No of Gauss Points = 15
    !------------------------------------
 
    case(15)
 
      p1 = 0.250000000000000
      p2 = 0.000000000000000
      p3 = 0.333333333333333
      p4 = 0.727272727272727
      p5 = 0.090909090909091
      p6 = 0.066550153573664
      p7 = 0.433449846426336

      r(1,:) = (/p1,p1,p1/)
      r(2,:) = (/p2,p3,p3/)
      r(3,:) = (/p3,p2,p3/)
      r(4,:) = (/p3,p3,p2/)
      r(5,:) = (/p3,p3,p3/)
      r(6,:) = (/p4,p5,p5/)
      r(7,:) = (/p5,p4,p5/)
      r(8,:) = (/p5,p5,p4/)
      r(9,:) = (/p5,p5,p5/)
      r(10,:) = (/p6,p7,p7/)
      r(11,:) = (/p7,p6,p7/)
      r(12,:) = (/p7,p7,p6/)
      r(13,:) = (/p7,p6,p6/)
      r(14,:) = (/p6,p7,p6/)
      r(15,:) = (/p6,p6,p7/)
      !
      q1 = 0.030283678097089
      q2 = 0.006026785714286
      q3 = 0.011645249086029
      q4 = 0.010949141561386

      w = (/q1,q2,q2,q2, q2,q3,q3,q3,q3,q4,q4,q4,q4,q4,q4/)

    case default
     print *, 'Invalid NGPTS for Tetrahedra Gauss quadrature'
     stop
   end select
  
  
end subroutine GaussTetrahedra

! ************************************************************************** !

subroutine GaussPyramid(NGPTS,r,w)
  ! 
  ! Calculates Gauss points for tetrahedra elements
  ! Reference:
  ! http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_pyramid/quadrature_rules_pyramid.html
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 7/11/2013
  ! 

  PetscInt :: EleType
  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  
  allocate(r(NGPTS,3))
  allocate(w(NGPTS))
    
  select case(NGPTS)
    case(1)

    !------------------------------
    ! No of Gauss Points = 1
    !------------------------------
      r(1,:) = (/0.0,0.0,0.25/)
      !
      w(1) = 1.0

    !--------------------------------
    ! No of Gauss Points = 5
    !--------------------------------

    case(5)
    
    
      r(1,:) = (/-0.48686449556014765641,-0.48686449556014765641,0.16666666666666666667/) 
      r(2,:) = (/ 0.48686449556014765641,-0.48686449556014765641,0.16666666666666666667/) 
      r(3,:) = (/ 0.48686449556014765641, 0.48686449556014765641,0.16666666666666666667/) 
      r(4,:) = (/-0.48686449556014765641, 0.48686449556014765641,0.16666666666666666667/)
      r(5,:) = (/ 0.00000000000000000000, 0.00000000000000000000,0.70000000000000000000/)
      !
      w = (/0.2109375, &
            0.2109375, &
            0.2109375, &
            0.2109375, &
            0.15625/)
      
    !------------------------------------
    ! No of Gauss Points = 6
    !------------------------------------     
    
    case(6)
    
      r(1,:) = (/-0.48795003647426658968,-0.48795003647426658968,0.16666666666666666667/)
      r(2,:) = (/ 0.48795003647426658968,-0.48795003647426658968,0.16666666666666666667/)
      r(3,:) = (/ 0.48795003647426658968, 0.48795003647426658968,0.16666666666666666667/)
      r(4,:) = (/-0.48795003647426658968, 0.48795003647426658968,0.16666666666666666667/)
      r(5,:) = (/ 0.00000000000000000000, 0.00000000000000000000,0.58333333333333333333/)
      r(6,:) = (/ 0.00000000000000000000, 0.00000000000000000000,0.75000000000000000000/)
      !
      w = (/0.21000000000000000000, &
            0.21000000000000000000, &
            0.21000000000000000000, &
            0.21000000000000000000, &
            0.06000000000000000000, &
            0.10000000000000000000/)

     case default
     print *, 'Invalid NGPTS for Tetrahedra Gauss quadrature'
     stop
   end select
  
  
end subroutine GaussPyramid

! ************************************************************************** !

subroutine Gauss3D(EleType,NGPTS,r,w)
  ! 
  ! Calculates Gauss points for 3D element
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  PetscInt :: EleType
  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  
  select case(EleType)
    case(HEX_TYPE)
      call GaussBrick(NGPTS,r,w)
    case(WEDGE_TYPE)
      call GaussWedge(NGPTS,r,w)
    case(TET_TYPE)
      call GaussTetrahedra(NGPTS,r,w)
    case(PYR_TYPE)
      call GaussPyramid(NGPTS,r,w)
    case default
      print *, 'Error: Only B8, W6, P5 and TET4 elements available for 3D.'
      stop
  end select
  
end subroutine Gauss3D

! ************************************************************************** !

subroutine GaussBrick(NGPTS,r,w)
  ! 
  ! Calculates Gauss points for B8 element
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal, pointer :: l(:,:)
  PetscReal, pointer :: m(:)
  PetscInt :: counter, i, j, k
  
  allocate(r(NGPTS*NGPTS*NGPTS,3))
  allocate(w(NGPTS*NGPTS*NGPTS))
  
  call Gauss1D(LINE_TYPE,NGPTS,l,m)
  
  ! Generate the B8 Gauss points and weights using for loops

  counter = 1
  do i = 1, NGPTS
    do j = 1, NGPTS
      do k = 1, NGPTS
        r(counter,1) = l(i,1)
        r(counter,2) = l(j,1)
        r(counter,3) = l(k,1)
        w(counter) = m(i)*m(j)*m(k)
        counter = counter + 1
      enddo
    enddo
  enddo
  
  deallocate(l)
  deallocate(m)
  
end subroutine GaussBrick

! ************************************************************************** !

subroutine GaussWedge(NGPTS,r,w)
  ! 
  ! Calculates Gauss points for wedge element
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 7/10//2013
  ! 

  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal, pointer :: rT3(:,:),rL2(:,:)  
  PetscReal, pointer :: wT3(:),wL2(:)
  PetscInt :: counter, i, j, k
  
  allocate(r(NGPTS*NGPTS,3))
  allocate(w(NGPTS*NGPTS))
  
  call Gauss1D(LINE_TYPE,NGPTS,rL2,wL2)
  call Gauss2D(TRI_TYPE,NGPTS,rT3,wT3)
  
  ! Generate the wedge Gauss points and weights using for loops
  do i = 1, NGPTS
    do j = 1, NGPTS
      r((i-1)*NGPTS+j,1) = rT3(i,1)
      r((i-1)*NGPTS+j,2) = rT3(i,2)
      r((i-1)*NGPTS+j,3) = rL2(j,1)
      w((i-1)*NGPTS+j) = wT3(i)*wL2(j)
    enddo
  enddo
  
  deallocate(rL2)
  deallocate(rT3)
  deallocate(wL2)
  deallocate(wT3)
  
end subroutine GaussWedge

! ************************************************************************** !

subroutine GaussDestroy(gauss)
  ! 
  ! Deallocate gauss type
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  ! 

  type(gauss_type) :: gauss
  
  deallocate(gauss%r)
  nullify(gauss%r)
  deallocate(gauss%w)
  nullify(gauss%w)

end subroutine GaussDestroy
     
end module Gauss_module
