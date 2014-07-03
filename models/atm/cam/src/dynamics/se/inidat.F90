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
  use cam_control_mod, only : ideal_phys, aqua_planet, pertlim
  implicit none
  private
  public read_inidat

contains



  subroutine read_inidat( ncid_ini, ncid_topo, dyn_in)
    use dyn_comp,      only: dyn_import_t
    use parallel_mod,     only: par
    use bndry_mod,     only: bndry_exchangev
    use constituents, only: cnst_name, cnst_read_iv, qmin
    use dimensions_mod,     only: nelemd, nlev, np
    use dof_mod, only           : putUniquePoints
    use edge_mod, only : edgevpack, edgevunpack, InitEdgeBuffer, FreeEdgeBuffer, EdgeBuffer_t
    use ncdio_atm, only : infld
    use shr_vmath_mod, only: shr_vmath_log
    use hycoef,           only: ps0
    use abortutils,     only: endrun
    use pio, only : file_desc_t, io_desc_t, pio_double, pio_get_local_array_size, pio_freedecomp
    use dyn_grid, only : get_horiz_grid_dim_d
    use chemistry   , only: chem_implements_cnst, chem_init_cnst
    use carma_intr,   only: carma_implements_cnst, carma_init_cnst
    use tracers     , only: tracers_implements_cnst, tracers_init_cnst
    use aoa_tracers , only: aoa_tracers_implements_cnst, aoa_tracers_init_cnst
    use stratiform,   only: stratiform_implements_cnst, stratiform_init_cnst
    use microp_driver, only: microp_driver_implements_cnst, microp_driver_init_cnst
    use phys_control,  only: phys_getopts
    use co2_cycle   , only: co2_implements_cnst, co2_init_cnst
    use nctopo_util_mod, only: nctopo_util_inidat
    implicit none
    type(file_desc_t),intent(inout) :: ncid_ini, ncid_topo
    type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

    type(element_t), pointer :: elem(:)
    real(r8), allocatable :: tmp(:,:)
    integer, allocatable :: gcols(:)
    integer :: tlncols, ig, ie, start, j, t, k
    character(len=40) :: fieldname
    logical :: found
    integer :: kptr, m_cnst
    type(EdgeBuffer_t) :: edge
    type(io_desc_t) :: iodesc
    integer :: lsize

    integer,parameter :: pcnst = PCNST
    integer, pointer :: gcid(:)

    integer :: rndm_seed_sz
    integer, allocatable :: rndm_seed(:)
    real(r8) :: pertval
    integer :: i
    real(r8), parameter :: D0_0 = 0.0_r8
    real(r8), parameter :: D0_5 = 0.5_r8
    real(r8), parameter :: D1_0 = 1.0_r8
    real(r8), parameter :: D2_0 = 2.0_r8
    character*16 :: subname='READ_INIDAT'

    if(iam < par%nprocs) then
       elem=> dyn_in%elem
    else
       nullify(elem)
    end if

    call get_dyn_decomp(elem, nlev, pio_double, iodesc)

    lsize = pio_get_local_array_size(iodesc)	

    tlncols = lsize/nlev

    allocate(tmp(tlncols,nlev))    

    if (iam < par%nprocs) then
      if(elem(1)%idxP%NumUniquePts <=0 .or. elem(1)%idxP%NumUniquePts > np*np) then
         write(iulog,*)  elem(1)%idxP%NumUniquePts
         call endrun('inidat')
      end if
    end if

    fieldname = 'U'
    call infld(fieldname, ncid_ini, iodesc, tlncols,'lev',  tmp, found)
    if(.not. found) then
       call endrun('Could not find U field on input datafile')
    end if
    
    start=1
    do ie=1,nelemd
       elem(ie)%state%v=0.0_r8
       call putUniquePoints(elem(ie)%idxP, nlev, tmp(start:,:), &
            elem(ie)%state%v(:,:,1,:,1))
       start=start+elem(ie)%idxP%numUniquePts
    end do

    fieldname = 'V'
    call infld(fieldname, ncid_ini, iodesc, tlncols,'lev',  tmp, found)
    if(.not. found) then
       call endrun('Could not find V field on input datafile')
    end if
    start=1
    do ie=1,nelemd
       call putUniquePoints(elem(ie)%idxP, nlev, tmp(start:,:), &
            elem(ie)%state%v(:,:,2,:,1))
       start=start+elem(ie)%idxP%numUniquePts
    end do

    fieldname = 'T'
    call infld(fieldname, ncid_ini, iodesc, tlncols ,'lev',  tmp, found)
    if(.not. found) then
       call endrun('Could not find T field on input datafile')
    end if
    start=1
    do ie=1,nelemd
       elem(ie)%state%T=0.0_r8
       call putUniquePoints(elem(ie)%idxP, nlev, tmp(start:,:), &
            elem(ie)%state%T(:,:,:,1))
       start=start+elem(ie)%idxP%numUniquePts
    end do

    if (pertlim .ne. D0_0) then
      if(masterproc) then
        write(iulog,*) trim(subname), ': Adding random perturbation bounded', &
                       'by +/- ', pertlim, ' to initial temperature field'
      end if

      call random_seed(size=rndm_seed_sz)
      allocate(rndm_seed(rndm_seed_sz))

      do ie=1,nelemd
        ! seed random number generator based on element ID
        ! (possibly include a flag to allow clock-based random seeding)
        rndm_seed = elem(ie)%GlobalId
        call random_seed(put=rndm_seed)
        do i=1,np
          do j=1,np
            do k=1,nlev
              call random_number(pertval)
              pertval = D2_0*pertlim*(D0_5 - pertval)
              elem(ie)%state%T(i,j,k,1) = elem(ie)%state%T(i,j,k,1)*(D1_0 + pertval)
            end do
          end do
        end do
      end do

      deallocate(rndm_seed)
    end if

    gcid => get_ldof(elem, 1)

    ! qmin = 1e-12,0,0

    do m_cnst=1,pcnst
       found = .false.
       if(cnst_read_iv(m_cnst)) then
          call infld(cnst_name(m_cnst), ncid_ini, iodesc, tlncols, 'lev', tmp, found)
       end if
       if(.not. found) then

          if(par%masterproc  ) write(iulog,*) 'Field ',cnst_name(m_cnst),' not found on initial dataset'

          if (microp_driver_implements_cnst(cnst_name(m_cnst))) then
             call microp_driver_init_cnst(cnst_name(m_cnst),tmp , gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "microp_driver_init_cnst"'
          else if (stratiform_implements_cnst(cnst_name(m_cnst))) then
             call stratiform_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "stratiform_init_cnst"'
          else if (chem_implements_cnst(cnst_name(m_cnst))) then
             call chem_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "chem_init_cnst"'
          else if (tracers_implements_cnst(cnst_name(m_cnst))) then
             call tracers_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "tracers_init_cnst"'
          else if (aoa_tracers_implements_cnst(cnst_name(m_cnst))) then
             call aoa_tracers_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "aoa_tracers_init_cnst"'
          else if (carma_implements_cnst(cnst_name(m_cnst))) then
             call carma_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "carma_init_cnst"'
          else if (co2_implements_cnst(cnst_name(m_cnst))) then
             call co2_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "co2_init_cnst"'
          else
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' set to 0.'
          end if
       end if
       do k=1,nlev
          do ig=1,tlncols
             tmp(ig,k)=max(qmin(m_cnst),tmp(ig,k))
          end do
       end do
       
       start=1
       do ie=1,nelemd
          elem(ie)%state%Q(:,:,:,m_cnst)=0.0_r8
          call putUniquePoints(elem(ie)%idxP, nlev, tmp(start:,:), &
               elem(ie)%state%Q(:,:,:,m_cnst))
          start=start+elem(ie)%idxP%numUniquePts
       end do
    end do
    deallocate(gcid)

    call pio_freedecomp(ncid_ini, iodesc)

    call get_dyn_decomp(elem, 1, pio_double, iodesc)

    fieldname = 'PS'
    call infld(fieldname, ncid_ini, iodesc, tmp(:,1), found)
    if(.not. found) then
       call endrun('Could not find PS field on input datafile')
    end if
    start=1

    if(minval(tmp(:,1)) < 10000) then
       call endrun('Problem reading ps field')
    end if

    do ie=1,nelemd
       elem(ie)%state%ps_v=0.0_r8
       call putUniquePoints(elem(ie)%idxP, tmp(start:start+elem(ie)%idxP%numUniquePts-1,1), &
            elem(ie)%state%ps_v(:,:,1))
       start=start+elem(ie)%idxP%numUniquePts
    end do

    if ( (ideal_phys .or. aqua_planet)) then
       tmp(:,:) = 0._r8
    else    
       fieldname = 'PHIS'
       call infld(fieldname, ncid_topo, iodesc, tmp(:,1), found)
       if(.not. found) then
          call endrun('Could not find PHIS field on input datafile')
       end if
    end if
    start=1
    do ie=1,nelemd
       elem(ie)%state%phis=0.0_r8
       call putUniquePoints(elem(ie)%idxP, tmp(start:,1), &
            elem(ie)%state%phis(:,:))
       start=start+elem(ie)%idxP%numUniquePts
    end do
    
    ! once we've read all the fields we do a boundary exchange to 
    ! update the redundent columns in the dynamics
    if(iam < par%nprocs) then
       call initEdgeBuffer(edge, (3+pcnst)*nlev+2)
    end if
    do ie=1,nelemd
       kptr=0
       call edgeVpack(edge, elem(ie)%state%ps_v(:,:,1),1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVpack(edge, elem(ie)%state%phis,1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,elem(ie)%desc)
       kptr=kptr+2*nlev
       call edgeVpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,elem(ie)%desc)
       kptr=kptr+nlev
       call edgeVpack(edge, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,elem(ie)%desc)
    end do
    if(iam < par%nprocs) then
       call bndry_exchangeV(par,edge)
    end if
    do ie=1,nelemd
       kptr=0
       call edgeVunpack(edge, elem(ie)%state%ps_v(:,:,1),1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVunpack(edge, elem(ie)%state%phis,1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVunpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,elem(ie)%desc)
       kptr=kptr+2*nlev
       call edgeVunpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,elem(ie)%desc)
       kptr=kptr+nlev
       call edgeVunpack(edge, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,elem(ie)%desc)
    end do

!$omp parallel do private(ie, t, m_cnst)
    do ie=1,nelemd
       do t=2,3
          elem(ie)%state%ps_v(:,:,t)=elem(ie)%state%ps_v(:,:,1)
          elem(ie)%state%v(:,:,:,:,t)=elem(ie)%state%v(:,:,:,:,1)
          elem(ie)%state%T(:,:,:,t)=elem(ie)%state%T(:,:,:,1)
       end do
       call shr_vmath_log(elem(ie)%state%ps_v,elem(ie)%state%lnps,size(elem(ie)%state%lnps))
    end do

    if(iam < par%nprocs) then
       call FreeEdgeBuffer(edge)
    end if

    !
    ! This subroutine is used to create nc_topo files, if requested
    ! 

    call nctopo_util_inidat(ncid_topo,iodesc,elem)

    call pio_freedecomp(ncid_topo, iodesc)
    deallocate(tmp)

  end subroutine read_inidat



  subroutine get_dyn_decomp(elem, nlev, datatype, iodesc)
    use pio, only : io_desc_t, pio_initdecomp
    use cam_pio_utils, only : pio_subsystem
    use dyn_grid, only : get_horiz_grid_dim_d

    type(element_t), pointer :: elem(:)
    integer, intent(in) :: nlev, datatype
    type(io_desc_t), intent(out) :: iodesc
    integer, pointer :: ldof(:)
    integer :: dimlens(2), dimcnt

    dimcnt=1
    call get_horiz_grid_dim_d(dimlens(1)) 
    if(nlev>1) then
       dimlens(2) = nlev
       dimcnt=dimcnt+1
    end if

    ldof => get_ldof(elem, nlev)
    call pio_initdecomp(pio_subsystem, datatype, dimlens(1:dimcnt), ldof, iodesc)

    deallocate(ldof)

  end subroutine get_dyn_decomp


  function get_ldof(elem, nlev) result(ldof)
    use dimensions_mod,     only: nelemd
    use dyn_grid, only : get_horiz_grid_dim_d
    use abortutils,     only: endrun

    type(element_t), pointer :: elem(:)
    integer, intent(in) :: nlev
    integer, pointer :: ldof(:)

    integer :: lcnt, ie, j, k, ig, numpts, offset, hdim


    call get_horiz_grid_dim_d(hdim)

    lcnt = 0
    do ie=1,nelemd
       lcnt = lcnt+nlev*elem(ie)%idxP%NumUniquePts
    end do
    allocate(ldof(lcnt))
    ig=1
    ldof(:) = 0
    do k=1,nlev
       do ie=1,nelemd
          numpts = elem(ie)%idxP%NumUniquePts
          offset = elem(ie)%idxP%UniquePtOffset
          do j=1,numpts
             ldof(ig)=offset+(j-1)+(k-1)*hdim
             ig=ig+1
          end do
       end do
    end do


  end function get_ldof





end module inidat
