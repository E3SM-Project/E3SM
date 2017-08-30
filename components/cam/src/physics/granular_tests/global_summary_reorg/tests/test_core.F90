module test_core

  use physics_types,  only: physics_state, physics_tend
  use glb_verif_smry, only: tp_stat_smry

  implicit none

  type(physics_state), pointer :: phys_state(:) => null()     ! shape: (begchunk:endchunk)
  type(physics_tend ), pointer :: phys_tend(:)  => null()     ! shape: (begchunk:endchunk)

  type(tp_stat_smry), pointer :: chunk_smry(:,:) => null()   ! shape: (begchunk:endchunk,nfld)
  type(tp_stat_smry), pointer :: domain_smry(:)  => null()   ! shape: (nfld)

  integer :: ncol = PCOLS-1  ! deliberately making ncol and pcols different to make sure 
                             ! the code works in such situations

contains

!@before
  !================================================================
  ! Initialization subroutine
  !  - allocation memory
  !  - read initial conditions
  !  - register tracers
  !  - register fields for which global summaries will be obtained
  !================================================================
  subroutine initialize

    use shr_kind_mod,  only: r8=>shr_kind_r8
    use cam_logfile,   only: iulog
    use ppgrid,        only: begchunk, endchunk, pver

    use constituents,  only: cnst_add
    use physconst,     only: mwdry, cpair, mwh2o, cpwv

    use physpkg,       only: phys_init

    use glb_verif_smry,only: SMALLER_THAN, GREATER_EQ, ABS_SMALLER_THAN, ABS_GREATER_EQ, CLIPPING, &
                             add_smry_field
    use glb_verif_smry,only: glb_verif_smry_frq, glb_verif_smry_level, l_print_smry_for_all_fields

    integer :: nstep = STEP
    integer :: idummy, icol, ichnk

    begchunk = BCHNK
    endchunk = ECHNK

    !-------------------------------
    ! open an ASCII file for output
    !-------------------------------
   !open(unit=iulog,file='unitTest_1',status='unknown')

    write(iulog,*)
    write(iulog,*) 'Active domain size: ',pver,' levels, ',ncol,' columns, ',(endchunk-begchunk+1),' chunks.'
    write(iulog,*)
    write(iulog,*) '# of cells per chunk: ',pver*ncol
    write(iulog,*) '# of cells in domain: ',pver*ncol*(endchunk-begchunk+1)
    write(iulog,*)

    !-------------------------------------
    ! Use the smry module in verbose mode
    !-------------------------------------
    glb_verif_smry_level        = 2      ! also print chunk smry
    l_print_smry_for_all_fields = .true.

    !-------------------------------
    ! Initialize tracer indices
    !-------------------------------
    call cnst_add('Q',      mwh2o, cpwv,  1.E-12_r8, idummy, longname='Specific humidity')
    call cnst_add('CLDLIQ', mwdry, cpair, 0._r8,     idummy, longname='Grid box averaged cloud liquid amount')

    call cnst_add('COLIDX',    mwdry, cpair, 0._r8,     idummy, longname='Column index as real value, for unit test')
    call cnst_add('NEGCOLIDX', mwdry, cpair, 0._r8,     idummy, longname='Negative column index as real value, for unit test')

    call cnst_add('NEGCOLIDX_FIX', mwdry, cpair, 0._r8, idummy, longname='Negative column index as real value, for unit test')
    !-------------------------------------------------------------------------
    ! Register fields for calculating global statistics summary.
    ! This has to be done before 'call phys_init' which allocates memory for 
    ! chunk_smry and domain_smry.
    !-------------------------------------------------------------------------
    call add_smry_field('Q @test_part_1','kg/kg',GREATER_EQ,  1.E-4_r8)
    call add_smry_field('Q @test_part_2','kg/kg',SMALLER_THAN,1.E-4_r8)

    call add_smry_field('CLDLIQ @test_part_1','kg/kg',GREATER_EQ,  1.E-9_r8)
    call add_smry_field('CLDLIQ @test_part_2','kg/kg',SMALLER_THAN,1.E-9_r8)

    call add_smry_field('COLIDX @test_part_1','-',SMALLER_THAN,5._r8)
    call add_smry_field('COLIDX @test_part_2','-',GREATER_EQ,  5._r8)

    call add_smry_field('NEGCOLIDX @test_part_1','-',ABS_SMALLER_THAN,5._r8)
    call add_smry_field('NEGCOLIDX @test_part_2','-',ABS_GREATER_EQ,  5._r8)

    call add_smry_field('NEGCOLIDX_FIX @test_part_1','-',SMALLER_THAN,-4.5_r8,fixer=CLIPPING)
    call add_smry_field('NEGCOLIDX_FIX @test_part_2','-',GREATER_EQ,  -4.0_r8,fixer=CLIPPING)

    !-------------------------------------------------------------------------------
    ! Allocate memory for state, tend, and stat vectors; read in initial conditions.
    !-------------------------------------------------------------------------------
    call phys_init(phys_state, phys_tend, chunk_smry, domain_smry, nstep)

    ! Re-assign values to CLDICE to facilitate verification
    do ichnk=begchunk,endchunk
    do icol = 1,ncol
       phys_state(ichnk)%q(icol,:,3) =  icol*1._r8
       phys_state(ichnk)%q(icol,:,4) = -icol*1._r8
       phys_state(ichnk)%q(icol,:,5) = -icol*1._r8
    end do
    end do

  end subroutine initialize

!@test
  !================================================================
  ! Main body of the test
  !================================================================
  subroutine test_glb_verif_smry

    use ppgrid,         only: begchunk, endchunk, pver
    use constituents, only: cnst_add, cnst_name
    use glb_verif_smry, only: get_smry_field_idx, get_chunk_smry, get_global_smry
    use cam_logfile,   only: iulog
    use cam_abortutils,only: endrun

    integer :: nstep = STEP
    integer :: itr, istat, istat1, istat2
    integer :: icnst, ichnk, nchnk, k
    integer :: n_tot_cnt_in_chunk (begchunk:endchunk,PCNST)
    integer :: n_tot_cnt_in_domain(PCNST)


    ! Test functionalities for getting global statistics summary

    ! The following chunk loop mimics the corresponding loop in phys_run1 (in which tphysbc is called)
    do ichnk=begchunk,endchunk

       write(iulog,*) '-----------------------------'
       write(iulog,*) '  chunk ',ichnk
       write(iulog,*) '-----------------------------'

       do icnst = 1,PCNST

         n_tot_cnt_in_chunk(ichnk,icnst) = 0

         itr = icnst
         call get_chunk_smry( trim(cnst_name(icnst))//' @test_part_1',       &! intent:in
                              ncol, pver,                                   &! intent: in
                              phys_state(ichnk)%q(:ncol,:,itr),             &! intent:inout, in
                              phys_state(ichnk)%lat, phys_state(ichnk)%lon, &! intent:in
                              chunk_smry(ichnk,:), istat )                   ! intent:inout, out

         n_tot_cnt_in_chunk(ichnk,icnst) = n_tot_cnt_in_chunk(ichnk,icnst) + chunk_smry(ichnk,istat)%count

         call get_chunk_smry( trim(cnst_name(icnst))//' @test_part_2',       &! intent:in
                              ncol, pver,                                   &! intent:in
                              phys_state(ichnk)%q(:ncol,:,itr),             &! intent:inout, in
                              phys_state(ichnk)%lat, phys_state(ichnk)%lon, &! intent:in
                              chunk_smry(ichnk,:), istat )                   ! intent:inout, out

         n_tot_cnt_in_chunk(ichnk,icnst) = n_tot_cnt_in_chunk(ichnk,icnst) + chunk_smry(ichnk,istat)%count

       end do

       !@assert
       ! check clipping.

       itr = PCNST
       istat = PCNST*2-1

       if (minval(phys_state(ichnk)%q(:ncol,:,itr)) .ne. chunk_smry(ichnk,istat)%threshold ) then
          write(iulog,*) "Expected min. = ",chunk_smry(ichnk,istat)%threshold
          write(iulog,*) "Actual   min. = ",minval(phys_state(ichnk)%q(:ncol,:,itr))
          call endrun('Found clipping error! -- (1)')
       end if

       istat = PCNST*2

       if (maxval(phys_state(ichnk)%q(:ncol,:,itr)) .ne. chunk_smry(ichnk,istat)%threshold ) then
          write(iulog,*) "Expected max. = ",chunk_smry(ichnk,istat)%threshold
          write(iulog,*) "Actual   max. = ",maxval(phys_state(ichnk)%q(:ncol,:,itr))
          call endrun('Found clipping error! -- (2)')
       end if

      !do k=1,PLEV
      !   write(iulog,'(a,i4,20f12.3)') "lev = ",k,phys_state(ichnk)%q(:ncol,k,itr) 
      !end do

    end do  ! ichnk=begchunk,endchunk

    !@assert

    if (any(n_tot_cnt_in_chunk/=ncol*pver)) then
       write(iulog,*) "Expected value of n_tot_cnt_in_chunk (all cols): ",ncol*pver
       write(iulog,*) "Actual values  of n_tot_cnt_in_chunk (all cols): ",n_tot_cnt_in_chunk
       call endrun('Test error in chunk_smry.')
    end if

    ! After all chunks have been processed, get the domain statistics summary

    nchnk = endchunk-begchunk+1
    write(iulog,*) '-----------------------------'
    write(iulog,*) '  entire domain on this CPU'
    write(iulog,*) '-----------------------------'

    call get_global_smry( chunk_smry, domain_smry, nstep)

    !----------------------------------
    ! Check if the results are correct
    !----------------------------------
    ! For each constituent, n_tot_cnt_in_domain should match the total number of cells in the domain

    do icnst = 1,PCNST

       itr = icnst
       call get_smry_field_idx(trim(cnst_name(icnst))//' @test_part_1',istat1)
       call get_smry_field_idx(trim(cnst_name(icnst))//' @test_part_2',istat2)

       n_tot_cnt_in_domain(icnst) = domain_smry(istat1)%count &
                                  + domain_smry(istat2)%count
    end do

    !@assert

    if (any(n_tot_cnt_in_domain/=ncol*nchnk*pver)) then

       write(iulog,*)
       write(iulog,*) '----'
       write(iulog,*) 'cnst idx, n_tot_cnt_in_domain, ncol*nchnk*pver: do they match?'
       write(iulog,*)

       do icnst = 1,PCNST
          write(iulog,*) icnst, n_tot_cnt_in_domain(icnst), ncol*nchnk*pver
       end do

       call endrun('Test error in domain_smry.')
    end if

    !@assert
    ! For each consitituent and stat_type, the total count for the domain shound match the sum of 
    ! the counts from all chunks.

    if (any(domain_smry(:)%count/=sum(chunk_smry(:,:)%count,dim=1))) then
       write(iulog,*) domain_smry(:)%count
       write(iulog,*) sum( chunk_smry(:,:)%count, dim=1 )
       call endrun('Test error. chunk_smry and domain_smry do not match.')
    end if

  end subroutine test_glb_verif_smry

!@after
  !================================================================
  ! Clean-up subroutine
  !================================================================
  subroutine finalize

   !use physpkg,      only: phys_final
    use cam_logfile,   only: iulog

    write(iulog,*) '==========================================='
    write(iulog,*) ' all results conformed with expected values'
    write(iulog,*) '==========================================='

   !close(unit=iulog)

  end subroutine finalize

!-------------------
end module test_core

