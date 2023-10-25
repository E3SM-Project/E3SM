module physconst

   ! Physical constants.  Use CCSM shared values whenever available.

   use shr_kind_mod,  only: r8 => shr_kind_r8
   use shr_const_mod, only: shr_const_g,      shr_const_stebol, shr_const_tkfrz,  &
                            shr_const_mwdair, shr_const_rdair,  shr_const_mwwv,   &
                            shr_const_latice, shr_const_latvap, shr_const_cpdair, &
                            shr_const_rhofw,  shr_const_cpwv,   shr_const_rgas,   &
                            shr_const_karman, shr_const_pstd,   shr_const_rhodair,&
                            shr_const_avogad, shr_const_boltz,  shr_const_cpfw,   &
                            shr_const_rwv,    shr_const_zvir,   shr_const_pi,     &
                            shr_const_rearth, shr_const_sday,   shr_const_cday,   &
                            shr_const_spval,  shr_const_omega,  shr_const_cpvir,  &
                            shr_const_tktrip
   use ppgrid,        only: pcols, pver, pverp, begchunk, endchunk   ! Dimensions and chunk bounds

   implicit none

   private
   public  :: physconst_init
   public  :: physconst_readnl
   public  :: physconst_update
   save
   
   ! Constantants for MAM spciesi classes
   integer, public, parameter :: spec_class_undefined  = 0
   integer, public, parameter :: spec_class_cldphysics = 1
   integer, public, parameter :: spec_class_aerosol    = 2
   integer, public, parameter :: spec_class_gas        = 3
   integer, public, parameter :: spec_class_other      = 4

   ! Constants based off share code or defined in physconst

   real(r8), public, parameter :: avogad      = shr_const_avogad     ! Avogadro's number (molecules/kmole)
   real(r8), public, parameter :: boltz       = shr_const_boltz      ! Boltzman's constant (J/K/molecule)
   real(r8), public, parameter :: cday        = shr_const_cday       ! sec in calendar day ~ sec
   real(r8), public, parameter :: cpair       = shr_const_cpdair     ! specific heat of dry air (J/K/kg)
   real(r8), public, parameter :: cpliq       = shr_const_cpfw       ! specific heat of fresh h2o (J/K/kg)
   real(r8), public, parameter :: karman      = shr_const_karman     ! Von Karman constant
   real(r8), public, parameter :: latice      = shr_const_latice     ! Latent heat of fusion (J/kg)
   real(r8), public, parameter :: latvap      = shr_const_latvap     ! Latent heat of vaporization (J/kg)
   real(r8), public, parameter :: pi          = shr_const_pi         ! 3.14...
   real(r8), public, parameter :: pstd        = shr_const_pstd       ! Standard pressure (Pascals)
   real(r8), public, parameter :: r_universal = shr_const_rgas       ! Universal gas constant (J/K/kmol)
   real(r8), public, parameter :: rhoh2o      = shr_const_rhofw      ! Density of liquid water (STP)
   real(r8), public, parameter :: spval       = shr_const_spval      !special value 
   real(r8), public, parameter :: stebol      = shr_const_stebol     ! Stefan-Boltzmann's constant (W/m^2/K^4)
   real(r8), public, parameter :: h2otrip     = shr_const_tktrip     ! Triple point temperature of water (K)

   real(r8), public, parameter :: c0          = 2.99792458e8_r8      ! Speed of light in a vacuum (m/s)
   real(r8), public, parameter :: planck      = 6.6260755e-34_r8     ! Planck's constant (J.s)

   ! Molecular weights
   real(r8), public, parameter :: mwco2       =  44._r8             ! molecular weight co2
   real(r8), public, parameter :: mwn2o       =  44._r8             ! molecular weight n2o
   real(r8), public, parameter :: mwch4       =  16._r8             ! molecular weight ch4
   real(r8), public, parameter :: mwf11       = 136._r8             ! molecular weight cfc11
   real(r8), public, parameter :: mwf12       = 120._r8             ! molecular weight cfc12
   real(r8), public, parameter :: mwo3        =  48._r8             ! molecular weight O3
   real(r8), public, parameter :: mwso2       =  64._r8
   real(r8), public, parameter :: mwso4       =  96._r8
   real(r8), public, parameter :: mwh2o2      =  34._r8
   real(r8), public, parameter :: mwdms       =  62._r8
   real(r8), public, parameter :: mwnh4       =  18._r8


   ! modifiable physical constants for aquaplanet or doubly periodic mode

   real(r8), public           :: gravit       = shr_const_g     ! gravitational acceleration (m/s**2)
   real(r8), public           :: sday         = shr_const_sday  ! sec in siderial day ~ sec
   real(r8), public           :: mwh2o        = shr_const_mwwv  ! molecular weight h2o
   real(r8), public           :: cpwv         = shr_const_cpwv  ! specific heat of water vapor (J/K/kg)
   real(r8), public           :: mwdry        = shr_const_mwdair! molecular weight dry air
   real(r8), public           :: rearth       = shr_const_rearth! radius of earth (m)
   real(r8), public           :: tmelt        = shr_const_tkfrz ! Freezing point of water (K)
   real(r8), public           :: omega        = shr_const_omega ! earth rot ~ rad/sec   

!---------------  Variables below here are derived from those above -----------------------

   real(r8), public           :: rga          = 1._r8/shr_const_g                 ! reciprocal of gravit
   real(r8), public           :: ra           = 1._r8/shr_const_rearth            ! reciprocal of earth radius
   real(r8), public           :: rh2o         = shr_const_rwv                     ! Water vapor gas constant ~ J/K/kg
   real(r8), public           :: rair         = shr_const_rdair   ! Dry air gas constant     ~ J/K/kg
   real(r8), public           :: epsilo       = shr_const_mwwv/shr_const_mwdair   ! ratio of h2o to dry air molecular weights 
   real(r8), public           :: zvir         = shr_const_zvir                    ! (rh2o/rair) - 1
   real(r8), public           :: cpvir        = shr_const_cpvir                   ! CPWV/CPDAIR - 1.0
   real(r8), public           :: rhodair      = shr_const_rhodair                 ! density of dry air at STP  ~ kg/m^3
   real(r8), public           :: cappa        = (shr_const_rgas/shr_const_mwdair)/shr_const_cpdair  ! R/Cp
   real(r8), public           :: ez           ! Coriolis expansion coeff -> omega/sqrt(0.375)   
   real(r8), public           :: Cpd_on_Cpv   = shr_const_cpdair/shr_const_cpwv

!---------------  Variables below here are for WACCM-X -----------------------
   real(r8), public, dimension(:,:,:), pointer :: cpairv ! composition dependent specific heat at constant pressure
   real(r8), public, dimension(:,:,:), pointer :: rairv  ! composition dependent gas "constant"
   real(r8), public, dimension(:,:,:), pointer :: cappav ! rairv/cpairv
   real(r8), public, dimension(:,:,:), pointer :: mbarv  ! composition dependent atmosphere mean mass
   real(r8), public, dimension(:,:,:), pointer :: kmvis  ! molecular viscosity      kg/m/s
   real(r8), public, dimension(:,:,:), pointer :: kmcnd  ! molecular conductivity   J/m/s/K
                         
!---------------  Variables below here are for turbulent mountain stress -----------------------

   real(r8), public           :: tms_orocnst
   real(r8), public           :: tms_z0fac

!================================================================================================
contains
!================================================================================================

  subroutine physconst_init()

    use cam_abortutils,   only: endrun

    implicit none
    
    integer :: ierr

!-------------------------------------------------------------------------------
!  Allocate constituent dependent properties 
!-------------------------------------------------------------------------------
    allocate( cpairv(pcols,pver,begchunk:endchunk), &
              rairv(pcols,pver,begchunk:endchunk),  &
              cappav(pcols,pver,begchunk:endchunk), &
              mbarv(pcols,pver,begchunk:endchunk),  &
              kmvis(pcols,pverp,begchunk:endchunk), &
              kmcnd(pcols,pverp,begchunk:endchunk), stat=ierr )
    if ( ierr /= 0 ) call endrun('physconst: allocate failed in physconst_init')
    

!-------------------------------------------------------------------------------
!  Initialize constituent dependent properties 
!-------------------------------------------------------------------------------    
    cpairv(:pcols,:pver,begchunk:endchunk) = cpair
    rairv(:pcols,:pver,begchunk:endchunk) = rair
    cappav(:pcols,:pver,begchunk:endchunk) = rair/cpair
    mbarv(:pcols,:pver,begchunk:endchunk) = mwdry

    return

  end subroutine physconst_init

!==============================================================================     

   ! Read namelist variables.
   subroutine physconst_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand
      use spmd_utils,      only: masterproc
      use cam_abortutils,      only: endrun
      use cam_logfile,     only: iulog

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'physconst_readnl'
      logical       newg, newsday, newmwh2o, newcpwv, newmwdry, newrearth, newtmelt, newomega

      ! Physical constants needing to be reset (ie. for aqua planet experiments)
      namelist /physconst_nl/  cpwv, gravit, mwdry, mwh2o, rearth, sday, tmelt, tms_orocnst, tms_z0fac, omega

      !-----------------------------------------------------------------------------

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'physconst_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, physconst_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(cpwv,      1,                   mpir8,   0, mpicom)
      call mpibcast(gravit,    1,                   mpir8,   0, mpicom)
      call mpibcast(mwdry,     1,                   mpir8,   0, mpicom)
      call mpibcast(mwh2o,     1,                   mpir8,   0, mpicom)
      call mpibcast(rearth,    1,                   mpir8,   0, mpicom)
      call mpibcast(sday,      1,                   mpir8,   0, mpicom)
      call mpibcast(tmelt,     1,                   mpir8,   0, mpicom)
      call mpibcast(tms_orocnst, 1,                 mpir8,   0, mpicom)
      call mpibcast(tms_z0fac, 1,                   mpir8,   0, mpicom)
      call mpibcast(omega,     1,                   mpir8,   0, mpicom)
#endif


      
      newg     =  gravit .ne. shr_const_g 
      newsday  =  sday   .ne. shr_const_sday
      newmwh2o =  mwh2o  .ne. shr_const_mwwv
      newcpwv  =  cpwv   .ne. shr_const_cpwv
      newmwdry =  mwdry  .ne. shr_const_mwdair
      newrearth=  rearth .ne. shr_const_rearth
      newtmelt =  tmelt  .ne. shr_const_tkfrz
      newomega =  omega  .ne. shr_const_omega
      
      if (newg .or. newsday .or. newmwh2o .or. newcpwv .or. newmwdry .or. newrearth .or. newtmelt .or. newomega) then
         if (masterproc) then
            write(iulog,*)'****************************************************************************'
            write(iulog,*)'***    New Physical Constant Values set via namelist                     ***'
            write(iulog,*)'***                                                                      ***'
            write(iulog,*)'***    Physical Constant    Old Value                  New Value         ***'
            if (newg)       write(iulog,*)'***       GRAVITY   ',shr_const_g,gravit,'***'
            if (newsday)    write(iulog,*)'***       SDAY      ',shr_const_sday,sday,'***'
            if (newmwh2o)   write(iulog,*)'***       MWH20     ',shr_const_mwwv,mwh2o,'***'
            if (newcpwv)    write(iulog,*)'***       CPWV      ',shr_const_cpwv,cpwv,'***'
            if (newmwdry)   write(iulog,*)'***       MWDRY     ',shr_const_mwdair,mwdry,'***'
            if (newrearth)  write(iulog,*)'***       REARTH    ',shr_const_rearth,rearth,'***'
            if (newtmelt)   write(iulog,*)'***       TMELT     ',shr_const_tkfrz,tmelt,'***'
            if (newomega)   write(iulog,*)'***       OMEGA     ',shr_const_omega,omega,'***'
            write(iulog,*)'****************************************************************************'
         end if
         rga         = 1._r8/gravit 
         ra          = 1._r8/rearth
         cpvir       = cpwv/cpair - 1._r8
         epsilo      = mwh2o/mwdry      

         if (.not.newomega) omega = 2.0_r8*pi/sday
         
         !  rair and rh2o have to be defined before any of the variables that use them
         
         rair        = r_universal/mwdry
         rh2o        = r_universal/mwh2o  
         
         cappa       = rair/cpair       
         rhodair     = pstd/(rair*tmelt)
         zvir        =  (rh2o/rair)-1.0_R8
         ez          = omega / sqrt(0.375_r8)
         Cpd_on_Cpv  = cpair/cpwv
         
      else
         ez          = omega / sqrt(0.375_r8)
      end if
      
    end subroutine physconst_readnl
    
!===============================================================================
  
  subroutine physconst_update(mmr, t, cnst_mw_o, cnst_mw_o2, cnst_mw_h, cnst_mw_n, ixo, ixo2, ixh, ncnst, lchnk, ncol)
  
!-----------------------------------------------------------------------
! Update the physics "constants" that vary
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------------------------------------

    integer, intent(in)  :: ncnst           ! number of constituents in mmr array
    real(r8), intent(in) :: mmr(pcols,pver,ncnst) ! constituents q array from state structure
    real(r8), intent(in) :: t(pcols,pver)   ! temperature t array from state structure

    real(r8), intent(in) :: cnst_mw_o, cnst_mw_o2, cnst_mw_h, cnst_mw_n  ! Molecular weight of O, O2, H, and N

    integer, intent(in)  :: ixo, ixo2, ixh  ! indices for O, O2, and H
    integer, intent(in)  :: lchnk           ! Chunk number
    integer, intent(in)  :: ncol            ! number of columns
!
!---------------------------Local storage-------------------------------------------------------------
    integer :: i,k                                 ! column,level,constituent indices

    real(r8):: mmro, mmro2, mmrh, mmrn2            ! Mass mixing ratios of O, O2, H, and N
    real(r8):: mbarvi, tint                        ! Mean mass, temperature, and specific heat on interface levels
    real(r8):: dof1, dof2                          ! Degress of freedom for cpairv calculation
    real(r8):: kv1, kv2, kv3, kv4                  ! Coefficients for kmvis calculation
    real(r8):: kc1, kc2, kc3, kc4                  ! Coefficients for kmcnd calculation
    !--------------------------------------------
    ! Set constants needed for updates
    !--------------------------------------------
    dof1 = 5._r8 
    dof2 = 7._r8
    kv1  = 4.03_r8
    kv2  = 3.42_r8
    kv3  = 3.9_r8
    kv4  = 0.69_r8
    kc1  = 56._r8
    kc2  = 56._r8
    kc3  = 75.9_r8
    kc4  = 0.69_r8

     !--------------------------------------------
     ! update cpairv, rairv, mbarv, and cappav
     !--------------------------------------------
     do k=1,pver
        do i=1,ncol
           mmro = mmr(i,k, ixo)
           mmro2 = mmr(i,k, ixo2)
           mmrh = mmr(i,k, ixh)
           mmrn2 = 1._r8-mmro-mmro2-mmrh
           mbarv(i,k,lchnk) = 1._r8/(mmro/cnst_mw_o +      &
                                     mmro2/cnst_mw_o2 +    &
                                     mmrn2/cnst_mw_n/2._r8 + &
                                     mmrh/cnst_mw_h)
           rairv(i,k,lchnk) = shr_const_rgas / mbarv(i,k,lchnk)
           cpairv(i,k,lchnk) = 0.5_r8*shr_const_rgas*    &
                                 (dof1*mmro/cnst_mw_o+        &
                                  dof2*mmro2/cnst_mw_o2+     &
                                  dof2*mmrn2/cnst_mw_n/2._r8+ &
                                     dof1*mmrh/cnst_mw_h)
           cappav(i,k,lchnk) = rairv(i,k,lchnk)/cpairv(i,k,lchnk)
        enddo
     enddo

     do k=2,pver
        do i=1,ncol
           mmro = .5_r8*(mmr(i,k-1, ixo)+mmr(i,k,ixo))
           mmro2 = .5_r8*(mmr(i,k-1, ixo2)+mmr(i,k,ixo2))
           mmrn2 = 1._r8-mmro-mmro2
           mbarvi = .5_r8*(mbarv(i,k-1,lchnk)+mbarv(i,k,lchnk))
           tint = .5_r8*(t(i,k-1)+t(i,k))
 
           kmvis(i,k,lchnk) = (kv1*mmro2/cnst_mw_o2+             &
                               kv2*mmrn2/cnst_mw_n/2._r8+        &
                               kv3*mmro/cnst_mw_o)*mbarvi*       &
                               tint**kv4 * 1.e-7_r8
           kmcnd(i,k,lchnk) = (kc1*mmro2/cnst_mw_o2+             &
                               kc2*mmrn2/cnst_mw_n/2._r8+        &
                               kc3*mmro/cnst_mw_o)*mbarvi*   &
                               tint**kc4 * 1.e-5_r8
        enddo
     enddo
     do i=1,ncol
        kmvis(i,1,lchnk) = 1.5_r8*kmvis(i,2,lchnk)-.5_r8*kmvis(i,3,lchnk)
        kmcnd(i,1,lchnk) = 1.5_r8*kmcnd(i,2,lchnk)-.5_r8*kmcnd(i,3,lchnk)
        kmvis(i,pverp,lchnk) = kmvis(i,pver,lchnk)
        kmcnd(i,pverp,lchnk) = kmcnd(i,pver,lchnk)
     enddo

   end subroutine physconst_update

end module physconst












