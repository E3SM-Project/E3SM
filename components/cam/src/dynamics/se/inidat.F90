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
  use scamMod, only: single_column, precip_off, scmlat, scmlon
  implicit none
  private
  public read_inidat

contains



  subroutine read_inidat( ncid_ini, ncid_topo, dyn_in)
    use dyn_comp,      only: dyn_import_t
    use parallel_mod,     only: par
    use bndry_mod,     only: bndry_exchangev
    use constituents, only: cnst_name, cnst_read_iv, qmin
    use dimensions_mod,     only: nelemd, nlev, np, npsq
    use dof_mod, only           : putUniquePoints
    use edge_mod, only : edgevpack, edgevunpack, InitEdgeBuffer, FreeEdgeBuffer
    use edgetype_mod, only : EdgeBuffer_t
    use ncdio_atm, only : infld
    use shr_vmath_mod, only: shr_vmath_log
    use hycoef,           only: ps0, hyam, hybm
    use cam_abortutils,     only: endrun
    use pio, only : file_desc_t, io_desc_t, pio_double, pio_get_local_array_size, pio_freedecomp
    use dyn_grid, only : get_horiz_grid_dim_d, dyn_decomp
    use chemistry   , only: chem_implements_cnst, chem_init_cnst
    use carma_intr,   only: carma_implements_cnst, carma_init_cnst
    use tracers     , only: tracers_implements_cnst, tracers_init_cnst
    use aoa_tracers , only: aoa_tracers_implements_cnst, aoa_tracers_init_cnst
    use clubb_intr,         only: clubb_implements_cnst, clubb_init_cnst
    use stratiform,   only: stratiform_implements_cnst, stratiform_init_cnst
    use microp_driver, only: microp_driver_implements_cnst, microp_driver_init_cnst
    use phys_control,  only: phys_getopts
    use co2_cycle   , only: co2_implements_cnst, co2_init_cnst
    use unicon_cam,          only: unicon_implements_cnst, unicon_init_cnst
    use nctopo_util_mod, only: nctopo_util_inidat
    use cam_history_support, only: max_fieldname_len
    use cam_grid_support,    only: cam_grid_get_local_size, cam_grid_get_gcid
    use cam_map_utils,       only: iMap
    use shr_const_mod,       only: SHR_CONST_PI
    use scamMod,             only: setiopupdate, readiopdata
    use se_single_column_mod, only: scm_setinitial
    implicit none
    type(file_desc_t),intent(inout) :: ncid_ini, ncid_topo
    type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

    real(r8), parameter :: rad2deg = 180.0 / SHR_CONST_PI
    type(element_t), pointer :: elem(:)
    real(r8), allocatable :: tmp(:,:,:)    ! (npsp,nlev,nelemd)
    real(r8), allocatable :: qtmp(:,:)     ! (npsp*nelemd,nlev)
    logical,  allocatable :: tmpmask(:,:)  ! (npsp,nlev,nelemd) unique grid val
    integer :: ie, k, t
    integer :: indx_scm, ie_scm, i_scm, j_scm
    character(len=max_fieldname_len) :: fieldname
    logical :: found
    integer :: kptr, m_cnst
    type(EdgeBuffer_t) :: edge
    integer :: lsize

    integer,parameter :: pcnst = PCNST
    integer(iMap), pointer :: ldof(:) => NULL() ! Basic (2D) grid dof
    integer,       pointer :: gcid(:) => NULL() ! ID based on ldof with no holes

    integer :: rndm_seed_sz
    integer, allocatable :: rndm_seed(:)
    real(r8) :: pertval
    integer :: sysclk
    integer :: i, j, indx
    real(r8), parameter :: D0_0 = 0.0_r8
    real(r8), parameter :: D0_5 = 0.5_r8
    real(r8), parameter :: D1_0 = 1.0_r8
    real(r8), parameter :: D2_0 = 2.0_r8
    real(r8) :: scmposlon, minpoint, testlat, testlon, testval 
    character*16 :: subname='READ_INIDAT'

    logical :: iop_update_surface

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
    tmp = 0.0_r8
    allocate(qtmp(npsq*nelemd,nlev))

    if (par%dynproc) then
      if(elem(1)%idxP%NumUniquePts <=0 .or. elem(1)%idxP%NumUniquePts > np*np) then
         write(iulog,*)  elem(1)%idxP%NumUniquePts
         call endrun(trim(subname)//': invalid idxP%NumUniquePts')
      end if
    end if

!   Determine column closest to SCM point
    if (single_column) then
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

    endif

    fieldname = 'U'
    tmp = 0.0_r8
    call infld(fieldname, ncid_ini, 'ncol', 'lev', 1, npsq,          &
         1, nlev, 1, nelemd, tmp, found, gridname='GLL')
    if(.not. found) then
       call endrun('Could not find U field on input datafile')
    end if
    
    do ie=1,nelemd
       elem(ie)%state%v=0.0_r8
       indx = 1
       do j = 1, np
          do i = 1, np
             elem(ie)%state%v(i,j,1,:,1) = tmp(indx,:,ie)
             if (single_column) elem(ie)%state%v(i,j,1,:,1)=tmp(indx_scm,:,ie_scm)
             indx = indx + 1
          end do
       end do
    end do

    fieldname = 'V'
    tmp = 0.0_r8
    call infld(fieldname, ncid_ini, 'ncol', 'lev', 1, npsq,          &
         1, nlev, 1, nelemd, tmp, found, gridname='GLL')
    if(.not. found) then
       call endrun('Could not find V field on input datafile')
    end if

    do ie=1,nelemd
       indx = 1
       do j = 1, np
          do i = 1, np
             elem(ie)%state%v(i,j,2,:,1) = tmp(indx,:,ie)
             if (single_column) elem(ie)%state%v(i,j,2,:,1) = tmp(indx_scm,:,ie_scm)
             indx = indx + 1
          end do
       end do
    end do

    fieldname = 'T'
    tmp = 0.0_r8
    call infld(fieldname, ncid_ini, 'ncol', 'lev', 1, npsq,          &
         1, nlev, 1, nelemd, tmp, found, gridname='GLL')
    if(.not. found) then
       call endrun('Could not find T field on input datafile')
    end if

    do ie=1,nelemd
       elem(ie)%state%T=0.0_r8
       indx = 1
       do j = 1, np
          do i = 1, np
             elem(ie)%state%T(i,j,:,1) = tmp(indx,:,ie)
             if (single_column) elem(ie)%state%T(i,j,:,1) = tmp(indx_scm,:,ie_scm)
             indx = indx + 1
          end do
       end do
    end do

    if (pertlim .ne. D0_0) then
      if(masterproc) then
        write(iulog,*) trim(subname), ': Adding random perturbation bounded', &
                       'by +/- ', pertlim, ' to initial temperature field'
      end if

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
              elem(ie)%state%T(i,j,k,1) = elem(ie)%state%T(i,j,k,1)*(D1_0 + pertval)
            end do
          end do
        end do
      end do

      deallocate(rndm_seed)
    end if

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
            call infld(cnst_name(m_cnst), ncid_ini, 'ncol', 'lev',      &
                 1, npsq, 1, nlev, 1, nelemd, tmp, found, gridname='GLL')
	    
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
          else if (carma_implements_cnst(cnst_name(m_cnst))) then
             call carma_init_cnst(cnst_name(m_cnst), qtmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "carma_init_cnst"'
          else if (co2_implements_cnst(cnst_name(m_cnst))) then
             call co2_init_cnst(cnst_name(m_cnst), qtmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "co2_init_cnst"'
          else if (unicon_implements_cnst(cnst_name(m_cnst))) then
             call unicon_init_cnst(cnst_name(m_cnst), qtmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "unicon_init_cnst"'
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
                if (single_column) elem(ie)%state%Q(i,j,:,m_cnst) = tmp(indx_scm,:,ie_scm)
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
    call infld(fieldname, ncid_ini, 'ncol',      &
         1, npsq, 1, nelemd, tmp(:,1,:), found, gridname='GLL')
    if(.not. found) then
       call endrun('Could not find PS field on input datafile')
    end if

    ! Check read-in data to make sure it is in the appropriate units
    allocate(tmpmask(npsq,nelemd))
    tmpmask = (reshape(ldof, (/npsq,nelemd/)) /= 0)

    if(minval(tmp(:,1,:), mask=tmpmask) < 10000._r8) then
       call endrun('Problem reading ps field')
    end if
    deallocate(tmpmask)

    do ie=1,nelemd
       elem(ie)%state%ps_v=0.0_r8
          indx = 1
          do j = 1, np
             do i = 1, np
                elem(ie)%state%ps_v(i,j,1) = tmp(indx,1,ie)
                if (single_column) elem(ie)%state%ps_v(i,j,1) = tmp(indx_scm,1,ie_scm)
                indx = indx + 1
             end do
          end do
    end do

    if ( (ideal_phys .or. aqua_planet)) then
       tmp(:,1,:) = 0._r8
    else    
       fieldname = 'PHIS'
       tmp(:,1,:) = 0.0_r8
       call infld(fieldname, ncid_topo, 'ncol',      &
            1, npsq, 1, nelemd, tmp(:,1,:), found, gridname='GLL')
       if(.not. found) then
          call endrun('Could not find PHIS field on input datafile')
       end if
    end if

    do ie=1,nelemd
       elem(ie)%state%phis=0.0_r8
       indx = 1
       do j = 1, np
          do i = 1, np
             elem(ie)%state%phis(i,j) = tmp(indx,1,ie)
             if (single_column) elem(ie)%state%phis(i,j) = tmp(indx_scm,1,ie_scm)
             indx = indx + 1
          end do
       end do
    end do
    
    if (single_column) then
      iop_update_surface = .false.
      call setiopupdate()
      call readiopdata(iop_update_surface,hyam,hybm)
      call scm_setinitial(elem)
    endif

    if (.not. single_column) then    

      ! once we've read all the fields we do a boundary exchange to 
      ! update the redundent columns in the dynamics
      if(par%dynproc) then
        call initEdgeBuffer(par, edge, elem, (3+pcnst)*nlev+2)
      end if
      do ie=1,nelemd
        kptr=0
        call edgeVpack(edge, elem(ie)%state%ps_v(:,:,1),1,kptr,ie)
        kptr=kptr+1
        call edgeVpack(edge, elem(ie)%state%phis,1,kptr,ie)
        kptr=kptr+1
        call edgeVpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,ie)
        kptr=kptr+2*nlev
        call edgeVpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,ie)
        kptr=kptr+nlev
        call edgeVpack(edge, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,ie)
      end do
      if(par%dynproc) then
        call bndry_exchangeV(par,edge)
      end if
      do ie=1,nelemd
        kptr=0
        call edgeVunpack(edge, elem(ie)%state%ps_v(:,:,1),1,kptr,ie)
        kptr=kptr+1
        call edgeVunpack(edge, elem(ie)%state%phis,1,kptr,ie)
        kptr=kptr+1
        call edgeVunpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,ie)
        kptr=kptr+2*nlev
        call edgeVunpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,ie)
        kptr=kptr+nlev
        call edgeVunpack(edge, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,ie)
      end do
    
    endif

!$omp parallel do private(ie, t, m_cnst)
    do ie=1,nelemd
       do t=2,3
          elem(ie)%state%ps_v(:,:,t)=elem(ie)%state%ps_v(:,:,1)
          elem(ie)%state%v(:,:,:,:,t)=elem(ie)%state%v(:,:,:,:,1)
          elem(ie)%state%T(:,:,:,t)=elem(ie)%state%T(:,:,:,1)
       end do
    end do

    if (.not. single_column) then
      if(par%dynproc) then
        call FreeEdgeBuffer(edge)
      end if
    endif

    !
    ! This subroutine is used to create nc_topo files, if requested
    ! 

    call nctopo_util_inidat(ncid_topo,elem)

    deallocate(tmp)

  end subroutine read_inidat

end module inidat
