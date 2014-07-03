module cosp_share

!------------------------------------------------------------------------------------------------------------------

! Purpose: CAM interface to provide vars to both CAM_history.F90 and cospsimulator_intr.F90 
! This module is needed to break circular dependencies.
! Programming notes: 
! 1) cam_history.F90 does not have an init procedure.  
! This is why variables here are public instead of using a get method.
! 2) Public variables defined here: nprs_cosp,ntau_cosp,ndbze_cosp,nsr_cosp,nhtmisr_cosp,nsza_cosp,
! prlim_cosp,taulim_cosp,dbzelim_cosp,srlim_cosp,htmisrlim_cosp,sza_cosp
! 3) Variables set here by setcospvalues using namelist input: nht_cosp, htlim_cosp, nscol_cosp,docosp_camhist
! 4) Also defined here are mid-points and mixed dimension variables for cam_history.F90
!
! Author:  J. Kay (jenkay@ucar.edu)
! Created: August 2009
! Last modified: August 13, 2010
!------------------------------------------------------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: pver

   implicit none
   private
   save					
	
! Public subroutines
   public :: setcospvalues

! Public variables to be used by both cam_history.F90 and cospsimulator_intr.F90

   ! whether or not to do CAM history variable/dimension creation for COSP
   logical, public :: docosp_camhist

   ! frequency at which cosp is called, every cosp_nradsteps radiation timestep
   integer, public :: cosp_nradsteps

   ! number of dimensions
   integer, public, parameter :: &
      	nprs_cosp = 7,    	&! number of pressure ranges
      	ntau_cosp  = 7,    	&! number of optical depth ranges
      	ntau_cosp_modis  = 6,   &! number of optical depth ranges MODIS
      	ndbze_cosp = 15,   	&! number of dBZe ranges for COSP radar simulator
      	nsr_cosp = 15,      	&! number of scattering ranges for COSP lidar simulator
      	nhtmisr_cosp = 16,  	&! number of heights for misr output (per Marchand)
      	nbnds_cosp = 2,  	&! number of bounds of cosp output (for cam_history.F90)
      	nsza_cosp = 5     	! number of solar zenith angle for COSP parasol output
   integer,public,parameter :: nhtml_cosp = pver  ! number of model levels is pver
   integer,public ::  nscol_cosp   		! number of subcolumns for COSP outputs.  
						! use namelist input Ncolumns to set
   integer,public ::  nht_cosp   		! number of height for COSP radar and lidar simulator outputs.  
						! set to 40 if csat_vgrid=.true., else set to Nlr

   ! limits of dimensions, used here to find mid-points, sza_cosp passed to cam_history.F90
   real(r8), public, parameter :: &
      	prslim_cosp_1d(nprs_cosp+1) = (/1000._r8, 800._r8, 680._r8, 560._r8, 440._r8,310._r8,180._r8,0._r8/),  &
      	taulim_cosp_1d(ntau_cosp+1) = (/0._r8, 0.3_r8, 1.3_r8, 3.6_r8, 9.4_r8, 23._r8, 60._r8, 379._r8/), &
      	taulim_cosp_modis_1d(ntau_cosp_modis+1) = (/0.3_r8, 1.3_r8, 3.6_r8, 9.4_r8, 23._r8, 60._r8, 100000._r8/), &
      	dbzelim_cosp_1d(ndbze_cosp+1) = (/-50._r8, -45._r8, -40._r8, -35._r8, -30._r8, -25._r8, -20._r8, -15._r8, &
		-10.0_r8, -5.0_r8, 0.0_r8, 5.0_r8, 10._r8, 15._r8, 20._r8, 25._r8/), &
      	srlim_cosp_1d(nsr_cosp+1) = (/0.01_r8, 1.2_r8, 3.0_r8, 5.0_r8, 7.0_r8, 10.0_r8, 15.0_r8, 20.0_r8, 25.0_r8, 30.0_r8, &
		40.0_r8, 50.0_r8, 60.0_r8, 80.0_r8, 999.0_r8, 1009.0_r8/), &
      	htmisrlim_cosp_1d(nhtmisr_cosp+1) = (/-99.0_r8, 0.0_r8, 0.5_r8, 1.0_r8, 1.5_r8, 2.0_r8, 2.5_r8, & 
		3.0_r8, 4.0_r8, 5.0_r8, 7.0_r8, 9.0_r8, 11.0_r8, 13.0_r8, 15.0_r8, 17.0_r8, 99.0_r8/), &
     	sza_cosp(nsza_cosp) = (/0.0_r8, 15.0_r8, 30.0_r8, 45.0_r8, 60.0_r8/)

   ! limits of dimensions, used by cam_history.F90, 2,nX_cosp  where 1,nX_cosp = min, 2,nX_cosp = max
   real(r8), public :: prslim_cosp(2,nprs_cosp)
   real(r8), public :: taulim_cosp(2,ntau_cosp)
   real(r8), public :: taulim_cosp_modis(2,ntau_cosp_modis)
   real(r8), public :: dbzelim_cosp(2,ndbze_cosp)
   real(r8), public :: srlim_cosp(2,nsr_cosp)
   real(r8), public :: htmisrlim_cosp(2,nhtmisr_cosp)

   ! variable declarations - known sizes

   ! mid points and dimension indices of dimensions
   real(r8),public :: prsmid_cosp(nprs_cosp)      		! pressure midpoints of COSP ISCCP output
   real(r8),public :: taumid_cosp(ntau_cosp)      		! optical depth midpoints of COSP ISCCP output
   real(r8),public :: taumid_cosp_modis(ntau_cosp_modis)      	! optical depth midpoints of COSP MODIS output
   real(r8),public :: dbzemid_cosp(ndbze_cosp)      		! dbze midpoints of COSP radar output
   real(r8),public :: srmid_cosp(nsr_cosp)      		! sr midpoints of COSP lidar output					
   real(r8),public :: htmisrmid_cosp(nhtmisr_cosp)      	! htmisr midpoints of COSP misr simulator output
   real(r8),public :: htmlmid_cosp(nhtml_cosp)      		! model level height midpoints for output

   integer,public :: prstau_cosp(nprs_cosp*ntau_cosp)		! ISCCP mixed output dimension index
   integer,public :: prstau_cosp_modis(nprs_cosp*ntau_cosp_modis)	! MODIS mixed output dimension index
   integer,public :: htmisrtau_cosp(nhtmisr_cosp*ntau_cosp)	! misr mixed output dimension index

   ! real values associated with the collapsed mixed dimensions
   real(r8),public :: prstau_prsmid_cosp(nprs_cosp*ntau_cosp)
   real(r8),public :: prstau_taumid_cosp(nprs_cosp*ntau_cosp)
   real(r8),public :: prstau_prsmid_cosp_modis(nprs_cosp*ntau_cosp_modis)
   real(r8),public :: prstau_taumid_cosp_modis(nprs_cosp*ntau_cosp_modis)
   real(r8),public :: htmisrtau_htmisrmid_cosp(nhtmisr_cosp*ntau_cosp)
   real(r8),public :: htmisrtau_taumid_cosp(nhtmisr_cosp*ntau_cosp)

   ! variable declarations - allocatable sizes
   real(r8),public,allocatable :: htlim_cosp(:,:)		! height limits for COSP outputs (nht_cosp+1)
   real(r8),public,allocatable :: htlim_cosp_1d(:)		! height limits for COSP outputs (nht_cosp+1)
   real(r8),public,allocatable :: htmid_cosp(:)		      	! height midpoints of COSP radar/lidar output (nht_cosp)
   integer,public,allocatable :: scol_cosp(:)      		! sub-column number (nscol_cosp)
   integer,public,allocatable :: htdbze_cosp(:)			! radar CFAD mixed output dimension index (nht_cosp*ndbze_cosp)
   integer,public,allocatable :: htsr_cosp(:)			! lidar CFAD mixed output dimension index (nht_cosp*nsr_cosp)
   integer,public,allocatable :: htmlscol_cosp(:)		! html-subcolumn mixed output dimension index (nhtml_cosp*nscol_cosp)
   real(r8),public,allocatable :: htdbze_htmid_cosp(:)		! (nht_cosp*ndbze_cosp)
   real(r8),public,allocatable :: htdbze_dbzemid_cosp(:)	! (nht_cosp*ndbze_cosp)
   real(r8),public,allocatable :: htsr_htmid_cosp(:)		! (nht_cosp*nsr_cosp)
   real(r8),public,allocatable :: htsr_srmid_cosp(:)		! (nht_cosp*nsr_cosp)
   real(r8),public,allocatable:: htmlscol_htmlmid_cosp(:)	! (nhtml_cosp*nscol_cosp)
   real(r8),public,allocatable :: htmlscol_scol_cosp(:)		! (nhtml_cosp*nscol_cosp)

   ! Local variables
   integer :: i,k		! indices
   real(r8) :: zstep

contains

!===============================================================================

subroutine setcospvalues(Nlr_in,use_vgrid_in,csat_vgrid_in,Ncolumns_in,docosp_in,cosp_nradsteps_in)

   ! input arguments
   integer, intent(in) :: Nlr_in
   logical, intent(in) :: use_vgrid_in
   logical, intent(in) :: csat_vgrid_in
   integer, intent(in) :: Ncolumns_in
   logical, intent(in) :: docosp_in
   integer, intent(in) :: cosp_nradsteps_in

   ! Local variables
   integer :: i,k		! indices
   real(r8) :: zstep

   ! set vertical grid, reference code from cosp_types.F90, line 549
   ! used to set vgrid_bounds in cosp_io.f90, line 844
   if (use_vgrid_in) then		!! using fixed vertical grid
   	if (csat_vgrid_in) then
     	   nht_cosp = 40
     	   zstep = 480.0_r8
   	else
     	   nht_cosp = Nlr_in
     	   zstep = 20000.0_r8/Nlr_in  ! constant vertical spacing, top at 20 km
   	end if
   end if

!  if (use_vgrid_in=.false.) then    !using the model vertical height grid
!	nht_cosp = pver
!	htlim_cosp = (/0._r8/) ##2check##
!  end if

   ! set number of sub-columns using namelist input
   nscol_cosp=Ncolumns_in

  ! set whether or not to do variable/dimension creation in CAM history
  if (docosp_in) then
	docosp_camhist = .true.
  else
	docosp_camhist = .false.
  end if

  cosp_nradsteps = cosp_nradsteps_in
  
  ! need to allocate memory for these variables
   allocate(htlim_cosp(2,nht_cosp),htlim_cosp_1d(nht_cosp+1),htmid_cosp(nht_cosp),scol_cosp(nscol_cosp),&
	htdbze_cosp(nht_cosp*ndbze_cosp),htsr_cosp(nht_cosp*nsr_cosp),htmlscol_cosp(nhtml_cosp*nscol_cosp),&
	htdbze_htmid_cosp(nht_cosp*ndbze_cosp),htdbze_dbzemid_cosp(nht_cosp*ndbze_cosp),&
	htsr_htmid_cosp(nht_cosp*nsr_cosp),htsr_srmid_cosp(nht_cosp*nsr_cosp),&
	htmlscol_htmlmid_cosp(nhtml_cosp*nscol_cosp),htmlscol_scol_cosp(nhtml_cosp*nscol_cosp))

   if (use_vgrid_in) then		!! using fixed vertical grid
      	htlim_cosp_1d(1)= 0.0_r8
   	do i=2,nht_cosp+1
      	   htlim_cosp_1d(i)=(i-1)*zstep      !! based on cosp_types.F90 line 556
   	enddo
   end if

   ! calculate mid-points and bounds for cam_history.F90

   do k=1,nprs_cosp
	prsmid_cosp(k) = 0.5_r8*(prslim_cosp_1d(k) + prslim_cosp_1d(k+1))
        prslim_cosp(1,k) = prslim_cosp_1d(k)
        prslim_cosp(2,k) = prslim_cosp_1d(k+1)
   end do

   do k=1,ntau_cosp
	taumid_cosp(k) = 0.5_r8*(taulim_cosp_1d(k) + taulim_cosp_1d(k+1))
        taulim_cosp(1,k) = taulim_cosp_1d(k)
        taulim_cosp(2,k) = taulim_cosp_1d(k+1)
   end do

   do k=1,ntau_cosp_modis
	taumid_cosp_modis(k) = 0.5_r8*(taulim_cosp_modis_1d(k) + taulim_cosp_modis_1d(k+1))
        taulim_cosp_modis(1,k) = taulim_cosp_modis_1d(k)
        taulim_cosp_modis(2,k) = taulim_cosp_modis_1d(k+1)
   end do

   do k=1,ndbze_cosp
	dbzemid_cosp(k) = 0.5_r8*(dbzelim_cosp_1d(k) + dbzelim_cosp_1d(k+1))
        dbzelim_cosp(1,k) = dbzelim_cosp_1d(k)
        dbzelim_cosp(2,k) = dbzelim_cosp_1d(k+1)
   end do

   do k=1,nsr_cosp
	srmid_cosp(k) = 0.5_r8*(srlim_cosp_1d(k) + srlim_cosp_1d(k+1))
        srlim_cosp(1,k) = srlim_cosp_1d(k)
        srlim_cosp(2,k) = srlim_cosp_1d(k+1)
   end do

   htmisrmid_cosp(1) = -99.0_r8
   htmisrlim_cosp(1,1) = htmisrlim_cosp_1d(1)
   htmisrlim_cosp(2,1) = htmisrlim_cosp_1d(2)
   do k=2,nhtmisr_cosp
	htmisrmid_cosp(k) = 0.5_r8*(htmisrlim_cosp_1d(k) + htmisrlim_cosp_1d(k+1))
        htmisrlim_cosp(1,k) = htmisrlim_cosp_1d(k)
        htmisrlim_cosp(2,k) = htmisrlim_cosp_1d(k+1)
   end do

   do k=1,nht_cosp
	htmid_cosp(k) = 0.5_r8*(htlim_cosp_1d(k) + htlim_cosp_1d(k+1))
        htlim_cosp(1,k) = htlim_cosp_1d(k)
        htlim_cosp(2,k) = htlim_cosp_1d(k+1)
   end do

   do k=1,nscol_cosp
	scol_cosp(k) = k
   end do

   !  Just using an index here, model height is a prognostic variable
   do k=1,nhtml_cosp
	htmlmid_cosp(k) = k
   end do

   ! assign mixed dimensions an integer index for cam_history.F90
   do k=1,nprs_cosp*ntau_cosp
	 prstau_cosp(k) = k
   end do
   do k=1,nprs_cosp*ntau_cosp_modis
	 prstau_cosp_modis(k) = k
   end do
   do k=1,nht_cosp*ndbze_cosp
	htdbze_cosp(k) = k
   end do
   do k=1,nht_cosp*nsr_cosp
	htsr_cosp(k) = k
   end do
   do k=1,nhtml_cosp*nscol_cosp
	htmlscol_cosp(k) = k
   end do
   do k=1,nhtmisr_cosp*ntau_cosp
	htmisrtau_cosp(k) = k
   end do

   ! next, assign collapsed reference vectors for cam_history.F90
   ! convention for saving output = prs1,tau1 ... prs1,tau7 ; prs2,tau1 ... prs2,tau7 etc.
   ! actual output is specified in cospsimulator_intr.F90
   do k=1,nprs_cosp
	prstau_taumid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=taumid_cosp(1:ntau_cosp)
	prstau_prsmid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=prsmid_cosp(k)
	prstau_taumid_cosp_modis(ntau_cosp_modis*(k-1)+1:k*ntau_cosp_modis)=taumid_cosp_modis(1:ntau_cosp_modis)
	prstau_prsmid_cosp_modis(ntau_cosp_modis*(k-1)+1:k*ntau_cosp_modis)=prsmid_cosp(k)
   enddo

   do k=1,nht_cosp
	htdbze_dbzemid_cosp(ndbze_cosp*(k-1)+1:k*ndbze_cosp)=dbzemid_cosp(1:ndbze_cosp)
	htdbze_htmid_cosp(ndbze_cosp*(k-1)+1:k*ndbze_cosp)=htmid_cosp(k)
   enddo

   do k=1,nht_cosp
	htsr_srmid_cosp(nsr_cosp*(k-1)+1:k*nsr_cosp)=srmid_cosp(1:nsr_cosp)
	htsr_htmid_cosp(nsr_cosp*(k-1)+1:k*nsr_cosp)=htmid_cosp(k)
   enddo


   do k=1,nhtml_cosp
	htmlscol_scol_cosp(nscol_cosp*(k-1)+1:k*nscol_cosp)=scol_cosp(1:nscol_cosp)
	htmlscol_htmlmid_cosp(nscol_cosp*(k-1)+1:k*nscol_cosp)=htmlmid_cosp(k)
   enddo

   do k=1,nhtmisr_cosp
	htmisrtau_taumid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=taumid_cosp(1:ntau_cosp)
	htmisrtau_htmisrmid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=htmisrmid_cosp(k)
   enddo

end subroutine setcospvalues

!===============================================================================


end module cosp_share
