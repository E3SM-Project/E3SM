       module mo_jeuv

       use shr_kind_mod, only : r8 => shr_kind_r8
       use abortutils,   only : endrun
       use cam_logfile,  only : iulog

       implicit none

       private
       public :: jeuv_init, jeuv, heuv, neuv
       public :: nIonRates

       save

!------------------------------------------------------------------------------
!      define EUV photolysis cross sections, branching ratios,
!      wavelength parameters,etc
!------------------------------------------------------------------------------
       integer, parameter :: neuv = 26               ! number of photolysis/ionization reactions
       integer, parameter :: nIonRates = 11          ! number of photo-ionizations rates needed for waccmx
       integer, parameter :: nmaj = 3                ! number of major neutral species (O,O2,N2)
       integer, parameter :: nspecies = 5            ! number of neutral species(O,N,O2,N2)
       integer, parameter :: nstat = 6               ! maximum number of ionization/dissociation
                                                     ! /excitation states for each speies
       integer, parameter :: lmax = 23               ! number of wavelength bins in EUV 
       real(r8), parameter :: heat_eff_fac = .05_r8  ! heating efficiency factor
       real(r8), parameter :: hc = 6.62608e-34_r8 * 2.9979e8_r8 / 1.e-9_r8

       real(r8) :: d2r
       real(r8) :: sigabs(lmax,nspecies)               ! absorption cross sections of major species
       real(r8) :: branch_p(lmax,nstat,nmaj) = 0._r8   ! branching ratios for photoionization/dissociation
       real(r8) :: branch_e(lmax,nstat,nmaj) = 0._r8   ! branching ratios for photoelectron ionization/dissociation/excitation
       real(r8) :: energy(lmax)                        ! solar energy

       contains

       subroutine jeuv_init (euvacdat_file, photon_file, electron_file, indexer)
!==============================================================================
!   Purpose:
!      read tabulated data:
!           (1) thermosphere neutral species' absorption cross sections, 
!               photoionization/dissociation branching ratios
!           (2) read photoelectron enhancement factor, photoelectron ionization/
!               dissociation/excitation branching ratios
!           (3) read solar flux
!==============================================================================

       use physconst,     only : pi
       use units,         only : getunit, freeunit
       use ppgrid,        only : pver
       use ioFileMod,     only : getfil
       use mo_chem_utls,  only : get_rxt_ndx

       implicit none

       character(len=*), intent(in) :: euvacdat_file
       character(len=*), intent(in) :: photon_file
       character(len=*), intent(in) :: electron_file
       integer, optional,intent(inout) :: indexer(:)

!------------------------------------------------------------------------------
!       local variables
!------------------------------------------------------------------------------
        integer  :: i, j, m                      ! loop indicies
        integer  :: unit                         ! fortran i/o unit number
        integer  :: istat                        ! file op status
        real(r8) :: waves(lmax)                  ! wavelength array for short bound of bins (A)
        real(r8) :: wavel(lmax)                  ! wavelength array for long bound of bins (A)
        real(r8) :: wc(lmax)                     ! wavelength bin center (nm)
        real(r8) :: sflux(lmax)                  ! solar flux (photon cm-2 s-1)
        real(r8) :: dummy                        ! temp variable
        character(len=200) :: str,fmt            ! string for comments in data file
        character(len=256) :: locfn

        integer :: jeuv_1_ndx

        jeuv_1_ndx = get_rxt_ndx( 'jeuv_1' )
        
        d2r = pi/180._r8
!------------------------------------------------------------------------------
!       read from data file the absorption cross sections for neutral species,
!       braching ratios for photoionization/dissociation, and braching ratios 
!       for photoelectron ionization/dissociation/excitation
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! read neutral species' absorption cross section and 
! photoionization/dissociation branching ratio
!------------------------------------------------------------------------------
	unit     = getunit()
        call getfil( photon_file, locfn, 0 )
        open( unit, file = trim(locfn), status='UNKNOWN', iostat=istat )
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to open ',trim(locfn),'; error = ',istat
           call endrun
        end if
!------------------------------------------------------------------------------
! read O
!------------------------------------------------------------------------------
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
        fmt = "(f7.2,2x,f7.2,2x,f9.3,4(2x,f7.3))"
	do i = 1,lmax
	   read(unit,fmt,iostat=istat) waves(i), wavel(i), sigabs(i,1), (branch_p(i,j,1),j=1,4)
           if( istat /= 0 ) then
              write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
              call endrun
           end if
        end do
!------------------------------------------------------------------------------
! read O2
!------------------------------------------------------------------------------
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	fmt = "(f7.2,2x,f7.2,2x,f9.3,5(2x,f7.3))"
	do i = 1,lmax
	   read(unit,fmt,iostat=istat) waves(i), wavel(i), sigabs(i,2), (branch_p(i,j,2),j=1,5)
           if( istat /= 0 ) then
              write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
              call endrun
           end if
        end do
!------------------------------------------------------------------------------
! read N2
!------------------------------------------------------------------------------
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	do i = 1,lmax
	   read(unit,fmt,iostat=istat) waves(i), wavel(i), sigabs(i,3), (branch_p(i,j,3),j=1,5)
           if( istat /= 0 ) then
              write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
              call endrun
           end if
        end do
!------------------------------------------------------------------------------
! read N
!------------------------------------------------------------------------------
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	fmt = "(f7.2,2x,f7.2,2x,f9.3)"
	do i = 1,lmax
	   read(unit,fmt,iostat=istat) waves(i), wavel(i), sigabs(i,4)
           if( istat /= 0 ) then
              write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
              call endrun
           end if
        end do

!------------------------------------------------------------------------------
! read CO2
!------------------------------------------------------------------------------
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	fmt = "(f7.2,2x,f7.2,2x,f9.3)"
	do i = 1,lmax
	   read(unit,fmt,iostat=istat) waves(i), wavel(i), sigabs(i,5)
           if( istat /= 0 ) then
              write(iulog,*) 'jeuv_init: failed to read CO2 data ',trim(locfn),'; error = ',istat
              call endrun
           end if
        end do

	close( unit )

!------------------------------------------------------------------------------
! unit transfer for absorption cross sections
! from Megabarns to cm^2
!------------------------------------------------------------------------------
	sigabs(:,:) = sigabs(:,:)*1.e-18_r8

!------------------------------------------------------------------------------
! read photoelectron ionization/dissociation/excitation branching ratio
!------------------------------------------------------------------------------
        call getfil( electron_file, locfn, 0 )
        open( unit, file = trim(locfn), status='UNKNOWN', iostat=istat )
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to open ',trim(locfn),'; error = ',istat
           call endrun
        end if
!------------------------------------------------------------------------------
! read O
!------------------------------------------------------------------------------
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	fmt="(f7.2,2x,f7.2,5(2x,f8.3))"
	do i = 1,lmax
	   read(unit,fmt,iostat=istat) waves(i), wavel(i), (branch_e(i,j,1),j=1,5)
           if( istat /= 0 ) then
              write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
              call endrun
           end if
        end do
!------------------------------------------------------------------------------
! read O2
!------------------------------------------------------------------------------
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	fmt = "(f7.2,2x,f7.2,6(2x,f8.3))"
	do i = 1,lmax
	   read(unit,fmt,iostat=istat) waves(i), wavel(i), (branch_e(i,j,2),j=1,6)
           if( istat /= 0 ) then
              write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
              call endrun
           end if
        end do
!------------------------------------------------------------------------------
! read N2
!------------------------------------------------------------------------------
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	do i = 1,lmax
	   read(unit,fmt,iostat=istat) waves(i), wavel(i), (branch_e(i,j,3),j=1,6)
           if( istat /= 0 ) then
              write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
              call endrun
           end if
        end do

	close( unit )

!------------------------------------------------------------------------------
! get solar flux
!------------------------------------------------------------------------------
        call getfil( euvacdat_file, locfn, 0 )
        open( unit, file = trim(locfn), status='UNKNOWN', iostat=istat )
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to open ',trim(locfn),'; error = ',istat
           call endrun
        end if

	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	read(unit,*,iostat=istat) str 
        if( istat /= 0 ) then
           write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
           call endrun
        end if
	do i = 1,lmax
	   read(unit,*,iostat=istat) waves(i), wavel(i), sflux(i), dummy
           if( istat /= 0 ) then
              write(iulog,*) 'jeuv_init: failed to read ',trim(locfn),'; error = ',istat
              call endrun
           end if
        end do

	close( unit )
	call freeunit( unit )

	wc(:)     = .1_r8*(waves(:) + wavel(:))          ! A to nm
	energy(:) = heat_eff_fac*hc/wc(:)

        do m = 1,neuv
           indexer(jeuv_1_ndx+m-1) = m
        enddo

	end subroutine jeuv_init

        subroutine jeuv( nlev, zen, occ, o2cc, n2cc, &
                         zkm, euv_prates)
!==============================================================================
!   Purpose:
!      Calculate euv photolysis/ionization rates (in s-1)
!==============================================================================
!   Arguments:
!       nlev: number of model vertical levels
!       zen:  solar zenith angle in degree
!       occ:  atmoic oxygen number density (#/cm3)
!       o2cc: molecular oxygen number density (#/cm3)
!       n2cc: molecular nitrogen number density (#/cm3)
!       zkm:  altitude of model levels in KM
!       euv_prates: array for EUV photolysis/ionization rates
!==============================================================================
!   Approach:
!       call sphers 
!            input: zenith angle
!            output: dsdh and nid used in slant column routine
!       call slant_col
!            input: dsdh and nid, etc
!            output: slant column density
!       calculate photon production rates and photoelectron production rates
!==============================================================================
!  Photolysis/ionization considered in EUV calculation
!
!  O + hv --> O+ (4S) + e*                        ! J1
!  O + hv --> O+ (2D) + e*                        ! J2
!  O + hv --> O+ (2P) + e*                        ! J3
!  N (4S) + hv --> N+ + e*                        ! J4
!  O2 + hv --> O2+ + e*                           ! J5
!  N2 + hv --> N2+ + e*                           ! J6
!  O2 + hv --> O + O+(4S) + e*                    ! J7
!  O2 + hv --> O + O+(2D) + e*                    ! J8
!  O2 + hv --> O + O+(2P) + e*                    ! J9
!  N2 + hv --> N (4S) + N+ + e*                   ! J10
!  N2 + hv --> N (2D) + N+ + e*                   ! J11
!  O2 + hv --> O (3P) + O (3P)                    ! J12
!  N2 + hv --> N (4S) + N (2D)                    ! J13
!
!  O + e* --> O+ (4S) + e* + e                    ! J14
!  O + e* --> O+ (2D) + e*+ e                     ! J15
!  O + e* --> O+ (2P) + e*+ e                     ! J16
!  O2 + e* --> O2+ + e*+ e                        ! J17
!  N2 + e*--> N2+ + e*+ e                         ! J18
!  O2 + e* --> O + O+(4S) + e*+ e                 ! J19
!  O2 + e*--> O + O+(2D) + e*+ e                  ! J20
!  O2 + e* --> O + O+(2P) + e*+ e                 ! J21
!  N2 + e* --> N (4S) + N+ + e*+ e                ! J22
!  N2 + e* --> N (2D) + N+ + e*+ e                ! J23
!  O2 + e* --> O (3P) + O (3P) + e*               ! J24
!  N2 + e* --> N (4S) + N (2D) + e*               ! J25
!==============================================================================

        use euvac,       only : euvac_etf
        use mo_jshort,   only : sphers, slant_col
        use cam_history, only : outfld

	implicit none

!------------------------------------------------------------------------------
!       dummy arguments
!------------------------------------------------------------------------------
	integer, intent(in)     :: nlev                        ! model vertical levels
	real(r8), intent(in)    :: zen	                       ! Zenith  angle in degree
        real(r8), intent(in)    :: occ(nlev)                   ! atmic oxygen number density (#/cm3) 
	real(r8), intent(in)    :: o2cc(nlev)		       ! Molecular oxygen number density (#/cm3) 
	real(r8), intent(in)    :: n2cc(nlev)		       ! molecular nitrogen number density(#/cm3) 
	real(r8), intent(in)    :: zkm(nlev)		       ! Altitude, km,from top to bottom
	real(r8), intent(out)   :: euv_prates(:,:)             ! EUV photolysis/ionization rates (s-1)   

!------------------------------------------------------------------------------
!       local variables
!------------------------------------------------------------------------------
        real(r8), parameter :: km2cm = 1.e5_r8
        integer  :: l, k, m, n              ! loop indecies
        real(r8) :: tau(lmax)               ! wavelength dependant optical depth
        real(r8) :: delz(nlev)              ! layer thickness (cm)
	real(r8) :: scol(nlev,nmaj)         ! major species' (O,O2,N2) Slant Column Density
	real(r8) :: nflux(nlev,lmax)
	real(r8) :: wrk(nmaj)               ! temporary array for photoabsorption rate
	real(r8) :: absorp(nlev,lmax)       ! temporary array for photoabsorption rate
	real(r8) :: ioniz(nlev,lmax)        ! temporary array for photoionization rate
	real(r8) :: dsdh(0:nlev,nlev)       ! Slant path of direct beam through each layer 
	                                    ! crossed  when travelling from the top of the atmosphere 
				            ! to layer i; dsdh(i,j), i = 0..NZ-1, j = 1..NZ-1
        integer :: nid(0:nlev)              ! Number of layers crossed by the direct
	                                    ! beam when travelling from the top of the 
                                            ! atmosphere to layer i; NID(i), i = 0..NZ-1
        real(r8) :: p_photon(nlev,nstat,nspecies)   !  photoionization/dissociation rates(s-1) (O,O2,N2,N)
        real(r8) :: p_electron(nlev,nstat,nmaj) !  photoelectron ionization/dissociation rates(s-1) (O,O2,N2)

        real(r8) :: prates(nlev,neuv)
!------------------------------------------------------------------------------
! zero arrays
!------------------------------------------------------------------------------
        p_photon(:,:,:)   = 0._r8
        p_electron(:,:,:) = 0._r8

        call sphers( nlev, zkm, zen, dsdh, nid )
	delz(1:nlev-1) = km2cm*(zkm(1:nlev-1) - zkm(2:nlev))
	call slant_col( nlev, delz, dsdh, nid, occ, scol )
	call slant_col( nlev, delz, dsdh, nid, o2cc, scol(1,2) )
	call slant_col( nlev, delz, dsdh, nid, n2cc, scol(1,3) )

!------------------------------------------------------------------------------
! The calculation is in the order from model bottom to model top
! because scol array is in this order.
!------------------------------------------------------------------------------
        do k = 1,nlev
           wrk(:) = scol(k,:)
           tau(:) = matmul( sigabs(:,:nmaj),wrk )
           where( tau(:) < 20._r8 )
	      nflux(k,:) = euvac_etf(:) * exp( -tau(:) )
           elsewhere
	      nflux(k,:) = 0._r8
           endwhere
        end do

!------------------------------------------------------------------------------
! remember occ,o2cc and n2cc is from top to bottom
!------------------------------------------------------------------------------
        do m = 1,nspecies
           do l = 1,lmax
              absorp(:,l) = sigabs(l,m) * nflux(:,l)
           end do
	   if( m <= nmaj ) then
              do l = 1,lmax
                 ioniz(:,l) = absorp(:,l) * branch_p(l,1,m)
              end do
	      do n = 1,nstat
	         p_photon(:,n,m)   = matmul( absorp,branch_p(:,n,m) )
	         p_electron(:,n,m) = matmul( ioniz,branch_e(:,n,m) )
              end do
           else
              p_photon(:,1,m) = matmul( nflux,sigabs(:,m) )
           end if
        end do

!------------------------------------------------------------------------------
! set photolysis/ionization rate for each reaction
!------------------------------------------------------------------------------
       prates(:,1)  = p_photon(:,2,1)  
       prates(:,2)  = p_photon(:,3,1)
       prates(:,3)  = p_photon(:,4,1)
       prates(:,4)  = p_photon(:,1,4)
       prates(:,5)  = p_photon(:,2,2) + p_photon(:,3,2)
       prates(:,6)  = p_photon(:,2,3) + p_photon(:,3,3)
       prates(:,7)  = .54_r8*p_photon(:,4,2)
       prates(:,8)  = .24_r8*p_photon(:,4,2)
       prates(:,9)  = .22_r8*p_photon(:,4,2)
       prates(:,10) = .2_r8*p_photon(:,4,3)
       prates(:,11) = .8_r8*p_photon(:,4,3)
       prates(:,12) = p_photon(:,5,2)
       prates(:,13) = p_photon(:,5,3)
       prates(:,14) = p_electron(:,2,1)
       prates(:,15) = p_electron(:,3,1)
       prates(:,16) = p_electron(:,4,1)
       prates(:,17) = p_electron(:,2,2) + p_electron(:,3,2)
       prates(:,18) = p_electron(:,2,3) + p_electron(:,3,3)
       prates(:,19) = .54_r8*p_electron(:,4,2)
       prates(:,20) = .24_r8*p_electron(:,4,2)
       prates(:,21) = .22_r8*p_electron(:,4,2)
       prates(:,22) = .2_r8*p_electron(:,4,3)
       prates(:,23) = .8_r8*p_electron(:,4,3)
       prates(:,24) = p_electron(:,5,2)
       prates(:,25) = p_electron(:,5,3)
       prates(:,26) = p_photon(:,1,5) 
 
      do m = 1,neuv
          euv_prates(:,m) = prates(nlev:1:-1,m)
       enddo

       end subroutine jeuv

       subroutine heuv( nlev, zen, occ, o2cc, n2cc, &
                        o_vmr, o2_vmr, n2_vmr, cparg, mw, &
			zkm, euv_hrates, kbot )
!==============================================================================
!   Purpose:
!      Calculate euv photolysis heating rates
!==============================================================================
!   Arguments:
!       nlev: number of model vertical levels
!       zen:  solar zenith angle in degree
!       occ:  atmoic oxygen number density (#/cm3)
!       o2cc: molecular oxygen number density (#/cm3)
!       n2cc: molecular nitrogen number density (#/cm3)
!       zkm:  altitude of model levels in KM
!       euv_prates: array for EUV photolysis/ionization rates
!==============================================================================
!   Approach:
!       call sphers 
!            input: zenith angle
!            output: dsdh and nid used in slant column routine
!       call slant_col
!            input: dsdh and nid, etc
!            output: slant column density
!       calculate photon production rates and photoelectron production rates
!==============================================================================
!  Photolysis/ionization considered in EUV heating rate calculation
!  O + hv --> O+ (4S) + e*                        ! J1
!  O + hv --> O+ (2D) + e*                        ! J2
!  O + hv --> O+ (2P) + e*                        ! J3
!  N (4S) + hv --> N+ + e*                        ! J4
!  O2 + hv --> O2+ + e*                           ! J5
!  N2 + hv --> N2+ + e*                           ! J6
!  O2 + hv --> O + O+(4S) + e*                    ! J7
!  O2 + hv --> O + O+(2D) + e*                    ! J8
!  O2 + hv --> O + O+(2P) + e*                    ! J9
!  N2 + hv --> N (4S) + N+ + e*                   ! J10
!  N2 + hv --> N (2D) + N+ + e*                   ! J11
!  O2 + hv --> O (3P) + O (3P)                    ! J12
!  N2 + hv --> N (4S) + N (2D)                    ! J13
!==============================================================================

        use euvac,         only : euvac_etf
        use mo_jshort,     only : sphers, slant_col
        use physconst,     only : avogad

	implicit none

!------------------------------------------------------------------------------
!       dummy arguments
!------------------------------------------------------------------------------
	integer, intent(in)     :: nlev                        ! model vertical levels
	integer, intent(in)     :: kbot                        ! heating vertical levels
	real(r8), intent(in)    :: zen	                       ! zenith  angle (degrees)
        real(r8), intent(in)    :: occ(nlev)                   ! atomic oxygen number density (molecules/cm^3)
	real(r8), intent(in)    :: o2cc(nlev)		       ! molecular oxygen concentration (molecules/cm^3)
	real(r8), intent(in)    :: n2cc(nlev)		       ! molecular nitrogen concentration (molecules/cm^3)
        real(r8), intent(in)    :: o_vmr(nlev)                 ! atomic oxygen concentration (mol/mol)
	real(r8), intent(in)    :: o2_vmr(nlev)		       ! molecular oxygen concentration (mol/mol)
	real(r8), intent(in)    :: n2_vmr(nlev)		       ! molecular nitrogen concentration (mol/mol)
	real(r8), intent(in)    :: zkm(nlev)		       ! midpoint geopotential (km)
	real(r8), intent(in)    :: cparg(nlev)		       ! specific heat capacity
	real(r8), intent(in)    :: mw(nlev)		       ! atm mean mass
	real(r8), intent(out)   :: euv_hrates(:)               ! euv heating rates

!------------------------------------------------------------------------------
!       local variables
!------------------------------------------------------------------------------
        real(r8), parameter :: km2cm = 1.e5_r8
        integer  :: l, k, k1, m, n          ! indicies
        real(r8) :: tau(lmax)               ! wavelength dependant optical depth
        real(r8) :: delz(nlev)              ! layer thickness (cm)
        real(r8) :: hfactor(kbot)
	real(r8) :: scol(nlev,nmaj)         ! major species' (O,O2,N2) Slant Column Density
	real(r8) :: nflux(kbot,lmax)
	real(r8) :: prates(kbot,13)         ! working photorates array
	real(r8) :: wrk(nmaj)               ! temporary array for photoabsorption rate
	real(r8) :: absorp(kbot,lmax)       ! temporary array for photoabsorption rate
	real(r8) :: ioniz(kbot,lmax)        ! temporary array for photoionization rate
	real(r8) :: dsdh(0:nlev,nlev)       ! Slant path of direct beam through each layer 
	                                    ! crossed  when travelling from the top of the atmosphere 
				            ! to layer i; dsdh(i,j), i = 0..NZ-1, j = 1..NZ-1
        integer :: nid(0:nlev)              ! Number of layers crossed by the direct
	                                    ! beam when travelling from the top of the 
                                            ! atmosphere to layer i; NID(i), i = 0..NZ-1
        real(r8) :: p_photon(kbot,nstat,nspecies)   !  photoionization/dissociation rates(s-1) (O,O2,N2,N)

!------------------------------------------------------------------------------
! zero arrays
!------------------------------------------------------------------------------
        p_photon(:,:,:)   = 0._r8

        call sphers( nlev, zkm, zen, dsdh, nid )
	delz(1:nlev-1) = km2cm*(zkm(1:nlev-1) - zkm(2:nlev))
	call slant_col( nlev, delz, dsdh, nid, occ, scol )
	call slant_col( nlev, delz, dsdh, nid, o2cc, scol(1,2) )
	call slant_col( nlev, delz, dsdh, nid, n2cc, scol(1,3) )

!------------------------------------------------------------------------------
! The calculation is in the order from model bottom to model top
! because scol array is in this order.
!------------------------------------------------------------------------------
        do k = 1,kbot
           k1 = nlev - kbot + k
           wrk(:) = scol(k1,:)
           tau(:) = matmul( sigabs(:,:nmaj),wrk )
           where( tau(:) < 20._r8 )
	      nflux(k,:) = energy(:) * euvac_etf(:) * exp( -tau(:) )
           elsewhere
	      nflux(k,:) = 0._r8
           endwhere
        end do

!------------------------------------------------------------------------------
! remember occ,o2cc and n2cc is from top to bottom
!------------------------------------------------------------------------------
        do m = 1,nspecies
           do l = 1,lmax
              absorp(:,l) = sigabs(l,m) * nflux(:,l)
           end do
	   if( m <= nmaj ) then
	      do n = 1,nstat
	         p_photon(:,n,m) = matmul( absorp,branch_p(:,n,m) )
              end do
           else
              p_photon(:,1,m) = matmul( nflux,sigabs(:,m) )
           end if
        end do

!------------------------------------------------------------------------------
! set photolysis rate
!------------------------------------------------------------------------------
       prates(:,1)  = p_photon(:,2,1)
       prates(:,2)  = p_photon(:,3,1)
       prates(:,3)  = p_photon(:,4,1)
       prates(:,5)  = p_photon(:,2,2) + p_photon(:,3,2)
       prates(:,6)  = p_photon(:,2,3) + p_photon(:,3,3)
       prates(:,7)  = .54_r8*p_photon(:,4,2)
       prates(:,8)  = .24_r8*p_photon(:,4,2)
       prates(:,9)  = .22_r8*p_photon(:,4,2)
       prates(:,10) = .2_r8*p_photon(:,4,3)
       prates(:,11) = .8_r8*p_photon(:,4,3)
       prates(:,12) = p_photon(:,5,2)
       prates(:,13) = p_photon(:,5,3)
       hfactor(:)   = avogad/(cparg(:kbot)*mw(:kbot))
       euv_hrates(kbot+1:nlev) = 0._r8
       euv_hrates(:kbot) = ((prates(kbot:1:-1,1) + prates(kbot:1:-1,2) + prates(kbot:1:-1,3))*o_vmr(:kbot) &
		     + (prates(kbot:1:-1,5) + prates(kbot:1:-1,7) + prates(kbot:1:-1,8) &
			+ prates(kbot:1:-1,9) + prates(kbot:1:-1,12))*o2_vmr(:kbot) &
		     + (prates(kbot:1:-1,6) + prates(kbot:1:-1,10) &
			+ prates(kbot:1:-1,11) + prates(kbot:1:-1,13))*n2_vmr(:kbot))*hfactor(:)

       end subroutine heuv

       end module mo_jeuv
