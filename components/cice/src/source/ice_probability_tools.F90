module ice_probability_tools

  use ice_kinds_mod
  use ice_blocks, only: nx_block, ny_block, nblocks_tot

  implicit none

    integer, parameter :: numCoeff = 10

    real (dbl_kind) :: perfmodel(numCoeff)
    character (len=char_len) :: perfmodel_name
    !=========================================================
    ! This performance model derived timing data and least-sq 
    ! approximation.
    !=========================================================
    ! minimizes: dynamics subcycling 
!    data perfmodel / 6.111e-05_dbl_kind, 3.0994e-2_dbl_kind, &
!                    -2.1127e-2_dbl_kind, 2.5825e-4_dbl_kind, &
!                     5.3287e-4_dbl_kind, 0.0000e-1_dbl_kind, &
!                     0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind, &
!                     0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind /
    ! minimizes: Column 
!    data perfmodel / 1.2253e-03_dbl_kind, 9.1685e-2_dbl_kind, &
!                    -7.4930e-2_dbl_kind,  5.7960e-2_dbl_kind, &
!                     5.0205e-3_dbl_kind,  2.0189e-3_dbl_kind /

    !---------------------------------------------------------
    ! minimizes: Column v2
    !---------------------------------------------------------
!    data perfmodel / 1.2253e-03_dbl_kind, 9.1685e-2_dbl_kind, &
!                    -7.4930e-2_dbl_kind,  1.0000e-3_dbl_kind, &
!                     5.0205e-3_dbl_kind,  2.0189e-3_dbl_kind, &
!		     0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind, &
!		     0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind /
!     data perfmodel_name /'Column v2'/
   
    
    !---------------------------------------------------------
    ! name: subycle v2
    ! minimizes: dynamics subcycle
    !---------------------------------------------------------
!    data perfmodel / 0.0000e-01_dbl_kind, 0.0000e-1_dbl_kind, &
!                     5.3586e-2_dbl_kind,  0.0000e-1_dbl_kind, &
!                     4.6250e-4_dbl_kind,  0.0000e-1_dbl_kind, &
!                    0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind, &
!                    0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind /
!     data perfmodel_name /'subcycle v2'/

    !---------------------------------------------------------
    ! name: subycle v3
    ! minimizes: dynamics subcycle
    !---------------------------------------------------------
    data perfmodel / 9.9777e-06_dbl_kind, 0.0000e-1_dbl_kind, &
                     1.9347e-2_dbl_kind,  0.0000e-1_dbl_kind, &
                     0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind, &
                     0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind, &
                     0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind /
     data perfmodel_name /'subcycle v3'/


    !---------------------------------------------------------
    ! name: combo v7
    ! minimizes: dynamics subcycle + Advection + Column
    !            Thermo + Shortwave + Ridging + Bound  
    !---------------------------------------------------------
!    data perfmodel / 9.9215e-04_dbl_kind, 0.0000e-1_dbl_kind, &
!                     1.1514e-1_dbl_kind,  0.0000e-1_dbl_kind, &
!                     1.6484e-2_dbl_kind,  8.7421e-4_dbl_kind, &
!		     0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind, &
!		     0.0000e-1_dbl_kind,  0.0000e-1_dbl_kind /
    
    ! minimizes: Thermo 
!    data perfmodel / 9.4472e-04_dbl_kind, 8.6489e-2_dbl_kind, &
!                    -6.6769-2_dbl_kind,   2.5593e-2_dbl_kind, &
!                     3.4758e-3_dbl_kind,  1.8534e-3_dbl_kind /

    ! minimizes: dynamics subcycling  + Column + Thermo + Ridging
!    data perfmodel / 6.1443e-04_dbl_kind, 1.0804e-1_dbl_kind, &
!                    -7.5396e-2_dbl_kind,  1.8444e-3_dbl_kind, &
!                     0.0000e-1_dbl_kind, -6.1136e-4_dbl_kind /
    ! minimizes: dynamics subcycling  + Column + Thermo + Ridging
!    data perfmodel / 6.2215e-03_dbl_kind, 6.3478e-2_dbl_kind, &
!                     0.0000e-1_dbl_kind,  1.2852e-1_dbl_kind, &
!                     -8.4564e-3_dbl_kind, 1.2788e-2_dbl_kind /

    ! manual selection 
!    data perfmodel /0.000e-1_dbl_kind, 0.0e-1_dbl_kind, &
!                    0.000e-1_dbl_kind, 0.0e-1_dbl_kind, &
!                    1.000e-2_dbl_kind, 3.0e-1_dbl_kind /    

    ! minimizes: dynamics
    !data perfmodel / 1.4455e-6_dbl_kind, 2.7671e-05_dbl_kind, &
    !                -4.7234e-05_dbl_kind,2.0304e-2_dbl_kind, &
    !                1.6523e-5_dbl_kind /
    ! minimizes: column 
    ! data perfmodel / 2.3000e-04_dbl_kind, 1.6847e-2_dbl_kind, &
    !                -1.2453e-2_dbl_kind, 0.0000e-0_dbl_kind, &
    !                 2.7042e-2_dbl_kind /

   integer (int_kind), allocatable, dimension(:) :: &
      nocn,                &! number of ocean points per block
      nice005,             &! number of ice points P > 0.005
      nice010,             &! number of ice points P > 0.01
      nice050,             &! number of ice points P > 0.05
      nice100,             &! number of ice points P > 0.10     
      nice250,		   &! number of ice points P > 0.25
      nice500              ! number of ice points P > 0.50
 
   public BuildProbabilityStats2

contains 

 subroutine BuildProbabilityStats2(blockLocation,coeffMatrix)

     integer (int_kind) :: blockLocation(:)
     real (dbl_kind)    :: coeffMatrix(:,:)

     integer (int_kind) :: maxNice,n,ip

       coeffMatrix=0.0
       maxNice = MAXVAL(nice005)
       do n=1,nblocks_tot
          ip = blockLocation(n)
          if(ip > 0) then
             coeffMatrix(1,ip) = coeffMatrix(1,ip) + real(nocn(n),dbl_kind)
             coeffMatrix(2,ip) = coeffMatrix(2,ip) + real(nice005(n),dbl_kind)
             coeffMatrix(3,ip) = coeffMatrix(3,ip) + real(nice010(n),dbl_kind)
             coeffMatrix(4,ip) = maxNice
             coeffMatrix(5,ip) = coeffMatrix(5,ip) + 2*(nx_block + ny_block)
             if( nice005(n) > 0) then
                 coeffMatrix(6,ip) = coeffMatrix(6,ip) + real(nx_block*ny_block,kind=dbl_kind)
             endif
             coeffMatrix(7,ip)  = coeffMatrix(7,ip) + real(nice050(n),dbl_kind)
             coeffMatrix(8,ip)  = coeffMatrix(8,ip) + real(nice100(n),dbl_kind)
             coeffMatrix(9,ip)  = coeffMatrix(9,ip) + real(nice250(n),dbl_kind)
             coeffMatrix(10,ip) = coeffMatrix(10,ip) + real(nice500(n),dbl_kind)
     
          endif
       enddo

 end subroutine BuildProbabilityStats2

   
end module ice_probability_tools

