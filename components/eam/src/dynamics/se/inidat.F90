module inidat
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: Read initial dataset and spectrally truncate as appropriate.
  !
  ! Method: Initialize one or a few fields at a time, to minimize the 
  !         memory  requirements
  ! 
  ! Author: 
  ! Modified: P. Worley, to implement initialization of subsets
  !           of fields. (8/03)
  !
  !           A. Gettelman and C. Craig (Nov 2010) - put micro/macro physics 
  !           into separate routines
  ! 
  !-----------------------------------------------------------------------
  use cam_logfile, only : iulog
  use element_mod, only : element_t
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: iam, masterproc
  use cam_control_mod, only : ideal_phys, aqua_planet, pertlim, seed_custom, seed_clock, new_random
  use random_xgc, only: init_ranx, ranx
  use scamMod, only: single_column, precip_off, scmlat, scmlon, scm_multcols, dp_crm, iop_perturb_high
  use perf_mod, only: t_startf, t_stopf

  implicit none
  private
  public read_inidat

contains

  subroutine read_inidat( ncid_ini, ncid_topo, dyn_in)
  
    use dyn_comp,                only: dyn_import_t, hvcoord, dom_mt
    use parallel_mod,            only: par
    use bndry_mod,               only: bndry_exchangev
    use constituents,            only: cnst_name, cnst_read_iv, qmin
    use dimensions_mod,          only: nelemd, nlev, np, npsq
    use dof_mod,                 only: putUniquePoints
    use edge_mod,                only : edgevpack_nlyr, edgevunpack_nlyr, edge_g
    use ncdio_atm,               only: infld
    use shr_vmath_mod,           only: shr_vmath_log
    use hycoef,                  only: ps0, hyam, hybm
    use cam_abortutils,          only: endrun
    use pio,                     only: file_desc_t, io_desc_t, pio_double, &
                                       pio_get_local_array_size, pio_freedecomp
    use dyn_grid,                only: get_horiz_grid_dim_d, dyn_decomp, fv_nphys
    use chemistry,               only: chem_implements_cnst, chem_init_cnst
    use tracers,                 only: tracers_implements_cnst, tracers_init_cnst
    use aoa_tracers,             only: aoa_tracers_implements_cnst, aoa_tracers_init_cnst
    use clubb_intr,              only: clubb_implements_cnst, clubb_init_cnst
    use stratiform,              only: stratiform_implements_cnst, stratiform_init_cnst
    use microp_driver,           only: microp_driver_implements_cnst, microp_driver_init_cnst
    use phys_control,            only: phys_getopts
    use co2_cycle,               only: co2_implements_cnst, co2_init_cnst
    use cam_history_support,     only: max_fieldname_len
    use cam_grid_support,        only: cam_grid_get_local_size, cam_grid_get_gcid
    use cam_map_utils,           only: iMap
    use shr_const_mod,           only: SHR_CONST_PI
    use scamMod,                 only: setiopupdate, readiopdata
    use se_single_column_mod,    only: scm_setinitial, scm_broadcast
    use element_ops,             only: set_thermostate
    use gllfvremap_mod,          only: gfr_fv_phys_to_dyn_topo

    implicit none
    type(file_desc_t),intent(inout) :: ncid_ini, ncid_topo
    type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

    real(r8), parameter :: rad2deg = 180.0 / SHR_CONST_PI
    type(element_t), pointer :: elem(:)
    real(r8), allocatable :: tmp(:,:,:)    ! (npsq,nlev,nelemd)
    real(r8), allocatable :: tmp_point(:,:)! (npsq,nlev)
    real(r8), allocatable :: qtmp(:,:)     ! (npsq*nelemd,nlev)
    real(r8) :: ps(np,np)     
    logical,  allocatable :: tmpmask(:,:)  ! (npsq,nlev,nelemd) unique grid val
    real(r8), allocatable :: phis_tmp(:,:) ! (nphys_sq,nelemd)
    integer :: nphys_sq                    ! # of fv physics columns per element
    integer :: ie, k, t
    integer :: indx_scm, ie_scm, i_scm, j_scm
    character(len=max_fieldname_len) :: fieldname
    logical :: found
    integer :: kptr, m_cnst
    integer :: lsize
    real(r8) :: p_ref(nlev)

    integer,parameter :: pcnst = PCNST
    integer(iMap), pointer :: ldof(:) => NULL() ! Basic (2D) grid dof
    integer,       pointer :: gcid(:) => NULL() ! ID based on ldof with no holes

    character(len=max_fieldname_len) :: ncol_name
    character(len=max_fieldname_len) :: grid_name
    logical :: read_pg_grid
    integer :: rndm_seed_sz
    integer, allocatable :: rndm_seed(:)
    real(r8) :: pertval
    integer :: sysclk
    integer :: i, j, indx, tl
    real(r8), parameter :: D0_0 = 0.0_r8
    real(r8), parameter :: D0_5 = 0.5_r8
    real(r8), parameter :: D1_0 = 1.0_r8
    real(r8), parameter :: D2_0 = 2.0_r8
    real(r8) :: scmposlon, minpoint, testlat, testlon, testval 
    character*16 :: subname='READ_INIDAT'
    integer :: nlev_tot

    logical :: iop_update_surface

    tl = 1

    if(par%dynproc) then
       elem=> dyn_in%elem
    else
       nullify(elem)
    end if

    lsize = cam_grid_get_local_size(dyn_decomp)	

    if (lsize /= (np*np*nelemd)) then
      call endrun(trim(subname)//': mismatch in local input array size')
    end if
    allocate(tmp(npsq,nlev,nelemd))
    allocate(tmp_point(1,nlev)) ! To find input at a single location
    allocate(qtmp(npsq*nelemd,nlev))

    if (fv_nphys>0) then
      nphys_sq = fv_nphys*fv_nphys
      allocate(phis_tmp(nphys_sq,nelemd))
    end if

    if (par%dynproc) then
      if(elem(1)%idxP%NumUniquePts <=0 .or. elem(1)%idxP%NumUniquePts > np*np) then
         write(iulog,*)  elem(1)%idxP%NumUniquePts
         call endrun(trim(subname)//': invalid idxP%NumUniquePts')
      end if
    end if

!   Determine column closest to SCM point
    if (single_column .and. .not. scm_multcols .and. par%dynproc) then
      if (scmlon .lt. 0._r8) then
        scmposlon=scmlon+360._r8
      else
        scmposlon=scmlon
      endif 
      minpoint=10000.0_r8
      ie_scm=0
      i_scm=0
      j_scm=0
      indx_scm=0
      do ie=1, nelemd
        indx=1
        do j=1, np
          do i=1, np
            testlat=elem(ie)%spherep(i,j)%lat * rad2deg
            testlon=elem(ie)%spherep(i,j)%lon * rad2deg
            if (testlon .lt. 0._r8) testlon=testlon+360._r8
            testval=abs(scmlat-testlat)+abs(scmposlon-testlon)
            if (testval .lt. minpoint) then
              ie_scm=ie
              indx_scm=indx
              i_scm=i
              j_scm=j
              minpoint=testval
            endif 
            indx=indx+1                   
          enddo
        enddo
      enddo

      if (ie_scm == 0 .or. i_scm == 0 .or. j_scm == 0 .or. indx_scm == 0) then
        call endrun('Could not find closest SCM point on input datafile')
      endif

    endif ! single_column

    if (scm_multcols) then
      indx_scm = 1
    endif

    grid_name = 'GLL'
    if (fv_nphys > 0) then
      ncol_name = 'ncol_d'
    else
      ncol_name = 'ncol'
    endif
    
    fieldname = 'U'
    tmp = 0.0_r8

    call t_startf('read_inidat_infld')
    if (.not. scm_multcols) then
      ! Standard Call
      call infld(fieldname, ncid_ini, ncol_name, 'lev', 1, npsq,          &
           1, nlev, 1, nelemd, tmp, found, gridname=grid_name)
    else
      ! Else find input just for the location of interest
      ! This logic follows for the rest of the input fields
      tmp_point = 0.0_r8
      call infld(fieldname, ncid_ini, ncol_name, 'lev', 1, 1,          &
           1, nlev, tmp_point, found, gridname=grid_name)
    endif
    call t_stopf('read_inidat_infld')

    if(.not. found) then
       call endrun('Could not find U field on input datafile')
    end if
    
    do ie=1,nelemd
       elem(ie)%state%v=0.0_r8
       indx = 1
       do j = 1, np
          do i = 1, np
             elem(ie)%state%v(i,j,1,:,tl) = tmp(indx,:,ie)
             ! If in SCM mode, put data of our column of interest
             ! in all dynamics columncs
             if (single_column .and. .not. scm_multcols) elem(ie)%state%v(i,j,1,:,tl)=tmp(indx_scm,:,ie_scm)
             if (scm_multcols) elem(ie)%state%v(i,j,1,:,tl)=tmp_point(indx_scm,:)
             indx = indx + 1
          end do
       end do
    end do

    fieldname = 'V'
    tmp = 0.0_r8

    call t_startf('read_inidat_infld')
    if (.not. scm_multcols) then
      call infld(fieldname, ncid_ini, ncol_name, 'lev', 1, npsq,          &
           1, nlev, 1, nelemd, tmp, found, gridname=grid_name)
    else
      tmp_point = 0.0_r8
      call infld(fieldname, ncid_ini, ncol_name, 'lev', 1, 1,          &
           1, nlev, tmp_point, found, gridname=grid_name)
    endif
    call t_stopf('read_inidat_infld')

    if(.not. found) then
       call endrun('Could not find V field on input datafile')
    end if

    do ie=1,nelemd
       indx = 1
       do j = 1, np
          do i = 1, np
             elem(ie)%state%v(i,j,2,:,tl) = tmp(indx,:,ie)
             if (single_column .and. .not. scm_multcols) elem(ie)%state%v(i,j,2,:,tl) = tmp(indx_scm,:,ie_scm)
             if (scm_multcols) elem(ie)%state%v(i,j,2,:,tl)=tmp_point(indx_scm,:)
             indx = indx + 1
          end do
       end do
    end do

    fieldname = 'T'
    tmp = 0.0_r8

    call t_startf('read_inidat_infld')
    if (.not. scm_multcols) then
      call infld(fieldname, ncid_ini, ncol_name, 'lev', 1, npsq,          &
           1, nlev, 1, nelemd, tmp, found, gridname=grid_name)
    else
      tmp_point = 0.0_r8
      call infld(fieldname, ncid_ini, ncol_name, 'lev', 1, 1,          &
           1, nlev, tmp_point, found, gridname=grid_name)
    endif
    call t_stopf('read_inidat_infld')

    if(.not. found) then
       call endrun('Could not find T field on input datafile')
    end if

    do ie=1,nelemd
#ifdef MODEL_THETA_L
       elem(ie)%derived%FT=0.0_r8
#else
       elem(ie)%state%T=0.0_r8
#endif
       indx = 1
       do j = 1, np
          do i = 1, np
#ifdef MODEL_THETA_L
             elem(ie)%derived%FT(i,j,:) = tmp(indx,:,ie)

             if (scm_multcols) elem(ie)%derived%FT(i,j,:) = tmp_point(indx_scm,:)
             if (single_column .and. .not. scm_multcols) elem(ie)%derived%FT(i,j,:) = tmp(indx_scm,:,ie_scm)
#else
             elem(ie)%state%T(i,j,:,tl) = tmp(indx,:,ie)

             if (single_column .and. .not. scm_multcols) elem(ie)%state%T(i,j,:,tl) = tmp(indx_scm,:,ie_scm)
             if (scm_multcols) elem(ie)%state%T(i,j,:,tl) = tmp_point(indx_scm,:)
#endif
             indx = indx + 1
          end do
       end do
    end do

    if (associated(ldof)) then
       call endrun(trim(subname)//': ldof should not be associated')
    end if
    call cam_grid_get_gcid(dyn_decomp, ldof)
    if (associated(gcid)) then
       call endrun(trim(subname)//': gcid should not be associated')
    end if
    allocate(gcid(size(ldof)))
    where (ldof == 0)
       gcid = 1
    elsewhere
       gcid = ldof
    end where

    ! qmin = 1e-12,0,0

    do m_cnst=1,pcnst
       found = .false.
       if(cnst_read_iv(m_cnst)) then

        ! If precip processes are turned off, do not initialize the field	
          if (precip_off .and. (cnst_name(m_cnst) .eq. 'RAINQM' .or. cnst_name(m_cnst) .eq. 'SNOWQM' &
            .or. cnst_name(m_cnst) .eq. 'NUMRAI' .or. cnst_name(m_cnst) .eq. 'NUMSNO')) then	    
            found = .false.
    
          else

            tmp = 0.0_r8

            call t_startf('read_inidat_infld')
            if (.not. scm_multcols) then
              call infld(cnst_name(m_cnst), ncid_ini, ncol_name, 'lev',      &
                   1, npsq, 1, nlev, 1, nelemd, tmp, found, gridname=grid_name)
            else
              tmp_point = 0.0_r8
              call infld(cnst_name(m_cnst), ncid_ini, ncol_name, 'lev', 1, 1,          &
                            1, nlev, tmp_point, found, gridname=grid_name)
            endif
            call t_stopf('read_inidat_infld')
    
          endif
       end if
       if(.not. found) then

          if(par%masterproc  ) write(iulog,*) 'Field ',cnst_name(m_cnst),' not found on initial dataset'

          if (microp_driver_implements_cnst(cnst_name(m_cnst))) then
             call microp_driver_init_cnst(cnst_name(m_cnst),qtmp , gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "microp_driver_init_cnst"'
          else if (clubb_implements_cnst(cnst_name(m_cnst))) then
             call clubb_init_cnst(cnst_name(m_cnst), qtmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "clubb_init_cnst"'
          else if (stratiform_implements_cnst(cnst_name(m_cnst))) then
             call stratiform_init_cnst(cnst_name(m_cnst), qtmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "stratiform_init_cnst"'
          else if (chem_implements_cnst(cnst_name(m_cnst))) then
             call chem_init_cnst(cnst_name(m_cnst), qtmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "chem_init_cnst"'
          else if (tracers_implements_cnst(cnst_name(m_cnst))) then
             call tracers_init_cnst(cnst_name(m_cnst), qtmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "tracers_init_cnst"'
          else if (aoa_tracers_implements_cnst(cnst_name(m_cnst))) then
             call aoa_tracers_init_cnst(cnst_name(m_cnst), qtmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "aoa_tracers_init_cnst"'
          else if (co2_implements_cnst(cnst_name(m_cnst))) then
             call co2_init_cnst(cnst_name(m_cnst), qtmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "co2_init_cnst"'
          else
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' set to 0.'
              qtmp = 0.0_r8
          end if
          ! Since the rest of processing uses tmp, copy qtmp into tmp
          do ie = 1, nelemd
            do k=1,nlev
              do i = 1, npsq
                ! Implicit reshape (qtmp is (np*np*nelemd, nlev)
                tmp(i,k,ie) = qtmp(i+((ie-1)*npsq),k)
              end do
            end do
          end do
       end if
       indx = 0
       do ie = 1, nelemd
          do k=1,nlev
             do i = 1, npsq
                ! Zero out the tmp values which might have been set
                ! erroneously by <param>_init_const
                if (ldof(indx + i) /= 0) then
                   ! Implicit reshape (qtmp is (np*np*nelemd, nlev)
                   tmp(i,k,ie)=max(qmin(m_cnst),tmp(i,k,ie))
                else
                   tmp(i,k,ie) = 0._r8
                end if
             end do
          end do
          indx = indx + npsq
       end do
       
       do ie=1,nelemd
          elem(ie)%state%Q(:,:,:,m_cnst)=0.0_r8
          indx = 1
          do j = 1, np
             do i = 1, np
                elem(ie)%state%Q(i,j,:,m_cnst) = tmp(indx,:,ie)
                if (single_column .and. .not. scm_multcols) elem(ie)%state%Q(i,j,:,m_cnst) = tmp(indx_scm,:,ie_scm)
                if (scm_multcols) elem(ie)%state%Q(i,j,:,m_cnst) = tmp_point(indx_scm,:)
                indx = indx + 1
             end do
          end do
       end do
    end do
    ! Cleanup
    if (associated(gcid)) then
      deallocate(gcid)
      nullify(gcid)
    end if

    fieldname = 'PS'
    tmp(:,1,:) = 0.0_r8
    call t_startf('read_inidat_infld')
    if (.not. scm_multcols) then
      call infld(fieldname, ncid_ini, ncol_name,      &
           1, npsq, 1, nelemd, tmp(:,1,:), found, gridname=grid_name)
    else
      call infld(fieldname, ncid_ini, ncol_name,      &
           1, 1, 1, 1, tmp(:,1,:), found, gridname=grid_name)
    endif
    call t_stopf('read_inidat_infld')
    if(.not. found) then
       call endrun('Could not find PS field on input datafile')
    end if

    ! Check read-in data to make sure it is in the appropriate units
    allocate(tmpmask(npsq,nelemd))
    tmpmask = (reshape(ldof, (/npsq,nelemd/)) /= 0)

    if(minval(tmp(:,1,:), mask=tmpmask) < 10000._r8 .and. .not. scm_multcols) then
       call endrun('Problem reading ps field')
    end if

    if (scm_multcols) then
      if (tmp(1,1,1) < 10000._r8) then
        call endrun('Problem reading ps field')
      endif
    endif

    deallocate(tmpmask)

    do ie=1,nelemd
       elem(ie)%state%ps_v=0.0_r8
          indx = 1
          do j = 1, np
             do i = 1, np
                elem(ie)%state%ps_v(i,j,tl) = tmp(indx,1,ie)
                if (single_column .and. .not. scm_multcols) elem(ie)%state%ps_v(i,j,tl) = tmp(indx_scm,1,ie_scm)
                if (scm_multcols) elem(ie)%state%ps_v(i,j,tl) = tmp(1,1,1)
                indx = indx + 1
             end do
          end do
    end do

    read_pg_grid = .false.
    if ( (ideal_phys .or. aqua_planet)) then
       tmp(:,1,:) = 0._r8
       if (fv_nphys > 0) phis_tmp(:,:) = 0._r8
    else    
      fieldname = 'PHIS'
      tmp(:,1,:) = 0.0_r8
      if (fv_nphys == 0) then
         call t_startf('read_inidat_infld')
         if (.not. scm_multcols) then
           call infld(fieldname, ncid_topo, ncol_name,      &
              1, npsq, 1, nelemd, tmp(:,1,:), found, gridname=grid_name)
         else
           call infld(fieldname, ncid_topo, ncol_name,      &
              1, 1, 1, 1, tmp(:,1,:), found, gridname=grid_name)
         endif
         call t_stopf('read_inidat_infld')

      else
         ! Attempt to read a mixed GLL-FV topo file, which contains PHIS_d in
         ! addition to PHIS.

         call t_startf('read_inidat_infld')
         call infld(trim(fieldname) // '_d', ncid_topo, ncol_name, &
              1, npsq, 1, nelemd, tmp(:,1,:), found, gridname=grid_name)
         call t_stopf('read_inidat_infld')

         if (found) then
            if (masterproc) then
               write(iulog,*) 'reading GLL ', trim(fieldname) // '_d', &
                    ' on gridname ', trim(grid_name)
            end if
         else
            ! Pure-FV topo file, so read FV PHIS and map it to GLL.
            if (masterproc) then
               write(iulog,*) 'reading FV ', trim(fieldname), &
                    ' on gridname physgrid_d'
            end if
            read_pg_grid = .true.

            call t_startf('read_inidat_infld')
            call infld(fieldname, ncid_topo, 'ncol', 1, nphys_sq, &
                 1, nelemd, phis_tmp, found, gridname='physgrid_d')
            call t_stopf('read_inidat_infld')

            call gfr_fv_phys_to_dyn_topo(par, dom_mt, elem, phis_tmp)
         end if
      end if
      if(.not. found) then
         call endrun('Could not find PHIS field on input datafile')
      end if
    end if

    if (.not. read_pg_grid) then
      do ie=1,nelemd
         elem(ie)%state%phis=0.0_r8
         indx = 1
         do j = 1, np
            do i = 1, np
               elem(ie)%state%phis(i,j) = tmp(indx,1,ie)
               if (single_column .and. .not. scm_multcols) elem(ie)%state%phis(i,j) = tmp(indx_scm,1,ie_scm)
               if (scm_multcols) elem(ie)%state%phis(i,j) = tmp(1,1,1)
               indx = indx + 1
            end do
         end do
      end do
    end if ! not read_pg_grid
    
    if (single_column) then
      iop_update_surface = .false.
      if (masterproc) call setiopupdate()
      if (masterproc) call readiopdata(iop_update_surface,hyam,hybm)
      if (scm_multcols) call scm_broadcast()
      call scm_setinitial(elem)
    endif

    if (pertlim .ne. D0_0) then
      if(masterproc) then
        write(iulog,*) trim(subname), ': Adding random perturbation bounded', &
                       'by +/- ', pertlim, ' to initial temperature field'
      end if

      if (dp_crm) then
        ! Define reference pressure, to potentially restrict initial perturbations
        !  to a certain height if requested
        do k=1,nlev
          p_ref(k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*hvcoord%ps0
        enddo
      endif

      if (new_random) then
        rndm_seed_sz = 1
      else
        call random_seed(size=rndm_seed_sz)
      endif
      allocate(rndm_seed(rndm_seed_sz))

      do ie=1,nelemd
        ! seed random number generator based on element ID
        ! (possibly include a flag to allow clock-based random seeding)
        rndm_seed(:) = elem(ie)%GlobalId
        if (seed_custom > 0) rndm_seed(:) = ieor( rndm_seed(1) , int(seed_custom,kind(rndm_seed(1))) )
        if (seed_clock) then
          call system_clock(sysclk)
          rndm_seed(:) = ieor( sysclk , int(rndm_seed(1),kind(sysclk)) )
        endif
        if (new_random) then
          call init_ranx(rndm_seed(1))
        else
          call random_seed(put=rndm_seed)
        endif
        do i=1,np
          do j=1,np
            do k=1,nlev

              if (new_random) then
                pertval = ranx()
              else
                call random_number(pertval)
              endif
              pertval = D2_0*pertlim*(D0_5 - pertval)

              ! If DP-CRM mode potentially only perturb a portion of the profile
              if (.not. dp_crm .or. p_ref(k) .gt. iop_perturb_high*100._r8) then

#ifdef MODEL_THETA_L
                elem(ie)%derived%FT(i,j,k) = elem(ie)%derived%FT(i,j,k)*(D1_0 + pertval)
#else
                elem(ie)%state%T(i,j,k,tl) = elem(ie)%state%T(i,j,k,tl)*(D1_0 + pertval)
#endif
              endif

            end do
          end do
        end do
      end do

      deallocate(rndm_seed)
    end if

    if (.not. single_column) then

      ! once we've read all the fields we do a boundary exchange to 
      ! update the redundent columns in the dynamics
      nlev_tot=(3+pcnst)*nlev+2

#ifdef MODEL_THETA_L
      do ie=1,nelemd
        kptr=0
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%ps_v(:,:,tl),1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%v(:,:,:,:,tl),2*nlev,kptr,nlev_tot)
        kptr=kptr+2*nlev
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FT(:,:,:),nlev,kptr,nlev_tot)
        kptr=kptr+nlev
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
      end do
#else
      do ie=1,nelemd
        kptr=0
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%ps_v(:,:,tl),1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%v(:,:,:,:,tl),2*nlev,kptr,nlev_tot)
        kptr=kptr+2*nlev
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%T(:,:,:,tl),nlev,kptr,nlev_tot)
        kptr=kptr+nlev
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
      end do
#endif
      if(par%dynproc) then
        call bndry_exchangeV(par,edge_g)
      end if
#ifdef MODEL_THETA_L
      do ie=1,nelemd
        kptr=0
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%ps_v(:,:,tl),1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%v(:,:,:,:,tl),2*nlev,kptr,nlev_tot)
        kptr=kptr+2*nlev
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FT(:,:,:),nlev,kptr,nlev_tot)
        kptr=kptr+nlev
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
      end do
#else
      do ie=1,nelemd
        kptr=0
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%ps_v(:,:,tl),1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%v(:,:,:,:,tl),2*nlev,kptr,nlev_tot)
        kptr=kptr+2*nlev
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%T(:,:,:,tl),nlev,kptr,nlev_tot)
        kptr=kptr+nlev
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
      end do
#endif    
    endif

!$omp parallel do private(ie, ps, t, m_cnst)
    do ie=1,nelemd
       ps=elem(ie)%state%ps_v(:,:,tl)
#ifdef MODEL_THETA_L
       elem(ie)%state%w_i = 0.0
       call set_thermostate(elem(ie),ps,elem(ie)%derived%FT,hvcoord,elem(ie)%state%Q(:,:,:,1))
       !FT used as tmp array - reset
       elem(ie)%derived%FT = 0.0
#else
       call set_thermostate(elem(ie),ps,elem(ie)%state%T(:,:,:,tl),hvcoord)
#endif
    end do

    deallocate(tmp)
    deallocate(tmp_point)
    deallocate(qtmp)
    if (fv_nphys>0) then
      deallocate(phis_tmp)
    end if

  end subroutine read_inidat

end module inidat
