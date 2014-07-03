module ch4varcon

  !-----------------------------------------------------------------------
  ! Module containing CH4 parameters and logical switches and routine to read constants from CLM namelist.
  !
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
  use clm_varctl  , only : NLFileName_in
  implicit none
  save
  !
  ! Methane Model Parameters
  !

  logical :: use_aereoxid_prog = .true. ! if false then aereoxid is read off of
  ! the parameter file and may be modifed by the user (default aereoxid on the
  ! file is 0.0).

  logical :: transpirationloss = .true. ! switch for activating CH4 loss from transpiration
                                      ! Transpiration loss assumes that the methane concentration in dissolved soil
                                      ! water remains constant through the plant and is released when the water evaporates
                                      ! from the stomata.
                                      ! Currently hard-wired to true; impact is < 1 Tg CH4/yr

  logical :: allowlakeprod = .false. ! Switch to allow production under lakes based on soil carbon dataset
                                     ! (Methane can be produced, and CO2 produced from methane oxidation,
                                     ! which will slowly reduce the available carbon stock, if ! replenishlakec, but no other biogeochem is done.)
                                     ! Note: switching this off turns off ALL lake methane biogeochem. However, 0 values
                                     ! will still be averaged into the concentration _sat history fields.

  logical :: usephfact = .false. ! Switch to use pH factor in methane production

  logical :: replenishlakec = .true. ! Switch for keeping carbon storage under lakes constant
                                      ! so that lakes do not affect the carbon balance
                                      ! Good for long term rather than transient warming experiments
               ! NOTE SWITCHING THIS OFF ASSUMES TRANSIENT CARBON SUPPLY FROM LAKES; COUPLED MODEL WILL NOT CONSERVE CARBON
               ! IN THIS MODE.

  ! New namelists added 6/12/11

  logical :: fin_use_fsat = .false. ! Use fsat rather than the inversion to Prigent satellite inundation obs. (applied to
                                    ! CLM water table depth and surface runoff) to calculated finundated which is
                                    ! used in methane code and potentially soil code
                                    !!!! Attn EK: Set this to true when Sean Swenson's prognostic, tested
                                       ! fsat is integrated. (CLM4 fsat is bad for these purposes.)

  logical :: usefrootc = .false.    ! Use CLMCN fine root C rather than ann NPP & LAI based parameterization to
                                    ! calculate tiller C for aerenchyma area calculation.
                                    ! The NPP & LAI param. was based on Wania for Arctic sedges and may not be
                                    ! appropriate for woody PFTs, although nongrassporosratio above partly adjusts
                                    ! for this.  However, using fine root C reduces the aerenchyma area by a large
                                    ! factor.

  logical :: ch4offline = .true.    ! true --> Methane is not passed between the land & atmosphere.
                                    ! NEM is not added to NEE flux to atm. to correct for methane production,
                                    ! and ambient CH4 is set to constant 2009 value.

  logical :: ch4rmcnlim = .false.   ! Remove the N and low moisture limitations on SOM HR when calculating
                                    ! methanogenesis.
                                    ! Note: this option has not been extensively tested.
                                    ! Currently hardwired off.

  logical :: anoxicmicrosites = .false. ! Use Arah & Stephen 1998 expression to allow production above the water table
                                        ! Currently hardwired off; expression is crude.

  logical :: ch4frzout = .false.    ! Exclude CH4 from frozen fraction of soil pore H2O, to simulate "freeze-out" pulse
                                    ! as in Mastepanov 2008.
                                    ! Causes slight increase in emissions in the fall and decrease in the spring.
                                    ! Currently hardwired off; small impact.

  public :: ch4conrd ! Read and initialize CH4 constants
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ch4conrd ()
    !
    ! !DESCRIPTION:
    ! Read and initialize CH4 constants
    !
    ! !USES:
    use fileutils   , only : relavu, getavu
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8, MPI_LOGICAL
    use shr_nl_mod  , only : shr_nl_find_group_name
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    implicit none
    !
    integer :: i,j,n                ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'ch4conrd'  ! subroutine name
    !-----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! Driver
    namelist /ch4par_in/ &
         ch4offline, fin_use_fsat, replenishlakec, allowlakeprod

    ! Production
    namelist /ch4par_in/ &
         usephfact 

    ! Methane
    namelist /ch4par_in/ &
         use_aereoxid_prog, usefrootc

       ! ----------------------------------------------------------------------
       ! Read namelist from standard input.
       ! ----------------------------------------------------------------------

    if (masterproc) then

       write(iulog,*) 'Attempting to read CH4 parameters .....'
       unitn = getavu()
       write(iulog,*) 'Read in ch4par_in namelist from: ', trim(NLFilename_in)
       open( unitn, file=trim(NLFilename_in), status='old' )
       call shr_nl_find_group_name(unitn, 'ch4par_in', status=ierr)
       if (ierr == 0) then
          read(unitn, ch4par_in, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg='error in reading in ch4par_in namelist'//&
                  errMsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if ! masterproc


    call mpi_bcast ( use_aereoxid_prog, 1 , MPI_LOGICAL, 0, mpicom, ierr )          
    call mpi_bcast (allowlakeprod, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (usephfact, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (replenishlakec, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (fin_use_fsat, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (usefrootc, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (ch4offline, 1 , MPI_LOGICAL, 0, mpicom, ierr)            

    if (masterproc) then
       write(iulog,*) 'Successfully read CH4 namelist'
       write(iulog,*)' '
       write(iulog,*)'allowlakeprod = ', allowlakeprod
       write(iulog,*)'usephfact = ', usephfact
       write(iulog,*)'replenishlakec = ', replenishlakec
       write(iulog,*)'fin_use_fsat = ', fin_use_fsat
       write(iulog,*)'usefrootc = ', usefrootc
       write(iulog,*)'ch4offline = ', ch4offline

       if (ch4offline) write(iulog,*)'CH4 Model will be running offline and not affect fluxes to atmosphere'

       write(iulog,*)'use_aereoxid_prog = ', use_aereoxid_prog 
       if ( .not. use_aereoxid_prog ) then
          write(iulog,*) 'Aerenchyma oxidation (aereoxid) value is being read from '//&
            ' the parameters file'
       endif

       if (.not. allowlakeprod) write(iulog,*) 'Lake production has been disabled. '// &
          '  Lakes will not factor into CH4 BGC.  "Sat" history fields will not average'// &
          '  over lakes except for concentrations, which will average zero from lakes.'
       if (.not. replenishlakec .and. .not. ch4offline) write(iulog,*)'LAKE SOIL CARBON WILL NOT BE REPLENISHED BUT INSTEAD ',&
             'WILL BE TRANSIENTLY RELEASED: COUPLED MODEL WILL NOT CONSERVE CARBON IN THIS MODE!'
       write(iulog,*)'Successfully initialized CH4 parameters from namelist.'
       write(iulog,*)

    end if

  end subroutine ch4conrd

end module ch4varcon

