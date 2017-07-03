module Well_TOilIms_class
#ifdef WELL_CLASS

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class
  use Well_Base_class
  use Well_Flow_class
  use Well_FlowEnergy_class
  use Well_WaterInjector_class
  use Well_OilProducer_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_water_injector_type) :: well_toil_ims_wat_inj_type
    ! ......................
  contains
    procedure, public :: PrintMsg => PrintTOilImsWatInj
    procedure, public :: ExplRes => TOilImsWatInjExplRes
    procedure, public  :: PrintOutputHeader => TOilImsWatInjOutputHeader
    procedure, public :: Output => TOilImsWatInjOutput
  end type well_toil_ims_wat_inj_type

  type, public, extends(well_oil_producer_type) :: well_toil_ims_oil_prod_type
    ! ......................
  contains
    procedure, public :: ExplRes => TOilImsOilProdExplRes
    procedure, public  :: PrintOutputHeader => TOilImsOilProdOutputHeader
    procedure, public :: Output => TOilImsOilProdOutput
  end type well_toil_ims_oil_prod_type

  !Below another data structure considered for the well
  !discarded in the attempt to build a water_injector, oil_producer, 
  ! etc reusable by other modules  
  !type, public, extends(well_flow_energy_type) :: well_toil_ims_type
  !  !class(auxvar_toil_ims_type), pointer :: toil_ims_auxvars(:,:)
  !  ! .................
  !contains  ! add here type-bound procedure 
  !  procedure, public  :: PrintOutputHeader => PrintOutputHeaderWellTOilIms
  !end type  well_toil_ims_type
  !type, public, extends(well_toil_ims_type) :: well_toil_ims_oil_prod_type
  ! ......................
  !end type

  public :: CreateTOilImsWell

contains

! ************************************************************************** !

subroutine PrintTOilImsWatInj(this)

  implicit none

  class(well_toil_ims_wat_inj_type) :: this

  write(*,*) "Well PrintTOilImsWatInj Printing message"

end subroutine PrintTOilImsWatInj

! ************************************************************************** !

subroutine TOilImsWatInjOutputHeader(this,output_option,file_unit)
  ! 
  ! Write header for well_TOilIms output file
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 
  use Output_Aux_module

  implicit none

  !class(well_toil_ims_type) :: this
  class(well_toil_ims_wat_inj_type) :: this
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: file_unit

  character(len=MAXSTRINGLENGTH) :: tunit

  tunit = trim(output_option%tunit)

  !TODO: can do something more clever than this: 
  !      e.g. small loop to add well vars
  write(IUNIT_TEMP,*) " VARIABLES = " // &
      '"Time [' // trim(tunit) // ']", ' // &
                '""Pw[Pa]"", ' // &
        '"Tw[C]", "dh2o[kg/m3]", ' // &
        '"Qwat[m3/' // trim(tunit) // ']", ' // &
        '"Mwat[kg/' // trim(tunit) // ']" ' 

!below to print all vars - 
!  write(IUNIT_TEMP,*) " VARIABLES = " // &
!      '"Time [' // trim(tunit) // ']", ' // &
!                '""Pw[Pa]"", ' // &
!        '"Tw[C]", "dh2o[kg/m3]", ' // &
!        '"doil[kg/' // trim(tunit) // ']", ' // &
!        '"Qwat[m3/' // trim(tunit) // ']", ' // &
!        '"Qoil[m3/' // trim(tunit) //  ']", ' // &
!        '"Mwat[kg/' // trim(tunit) // ']", ' // &
!        '"Moil[kg/' // trim(tunit) // ']"' 



end subroutine TOilImsWatInjOutputHeader

! ************************************************************************** !

subroutine TOilImsOilProdOutputHeader(this,output_option,file_unit)
  ! 
  ! Write header for TOilIms_oil_producer output file
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 06/12/16
  ! 
  use Output_Aux_module

  implicit none

  class(well_toil_ims_oil_prod_type) :: this
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: file_unit

  character(len=MAXSTRINGLENGTH) :: tunit

  tunit = trim(output_option%tunit)

  !TODO: can do something more clever than this: 
  !      e.g. small loop to add well vars
  write(IUNIT_TEMP,*) " VARIABLES = " // &
      '"Time [' // trim(tunit) // ']", ' // &
                '"Pw[Pa]", ' // &
        '"Tw[C]", "den_well_fluid[kg/m3]", ' // &
        '"Qwat[m3/' // trim(tunit) // ']", ' // &
        '"Mwat[kg/' // trim(tunit) // ']" '  // &
        '"Qoil[m3/' // trim(tunit) // ']", ' // &
        '"Moil[kg/' // trim(tunit) // ']" '  // &
        '"vol_WOR[-]" '

end subroutine TOilImsOilProdOutputHeader

! ************************************************************************** !

subroutine TOilImsWatInjOutput(this,output_file_unit,output_option,option)
  ! 
  ! Write output file for TOilImsWatInj
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 
  use Option_module
  use Output_Aux_module

  implicit none

  !class(well_toil_ims_type) :: this
  class(well_toil_ims_wat_inj_type) :: this
  PetscInt, intent(in) :: output_file_unit
  type(output_option_type), intent(in) :: output_option
  type(option_type) :: option

  !character(len=MAXWORDLENGTH) :: src_name
  !type(output_option_type), intent(in) :: output_option

  !character(len=MAXWORDLENGTH) :: wfile_name
  !PetscMPIInt :: cur_w_myrank
  !PetscInt :: ios, ierr

  !if( this%connection_set%num_connections > 0 ) then
  !  call MPI_Comm_rank(this%comm, cur_w_myrank, ierr )  
  !  if(this%cntr_rank == cur_w_myrank ) then
  !    wfile_name = trim(option%global_prefix) // "_" // &
  !                      trim(src_name) // ".tec" 
  !    open(unit=IUNIT_TEMP,file=wfile_name,action="write", &
  !         position="append",status="old",iostat=ios)

  write(output_file_unit,"(9(E10.4,1x))") option%time/output_option%tconv , &
                                          this%pw_ref, &
                                          this%tw_ref, &
                                          this%dw_kg_ref(LIQUID_PHASE), &
                                          this%q_fld(LIQUID_PHASE) * &
                                          output_option%tconv, &
                                          this%mr_fld(LIQUID_PHASE) * &
                                          output_option%tconv

  !  end if
  !end if 

end subroutine TOilImsWatInjOutput

! ************************************************************************** !

subroutine TOilImsOilProdOutput(this,output_file_unit,output_option,option)
  ! 
  ! Write output file for TOilImsWatInj
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 
  use Option_module
  use Output_Aux_module

  implicit none

  class(well_toil_ims_oil_prod_type) :: this
  PetscInt, intent(in) :: output_file_unit
  type(output_option_type), intent(in) :: output_option
  type(option_type) :: option

  PetscReal :: vol_WOR

  vol_WOR = this%q_fld(option%liquid_phase) / &
            ( this%q_fld(option%liquid_phase) + this%q_fld(option%oil_phase) )

  write(output_file_unit,"(9(E10.4,1x))") option%time/output_option%tconv , &
                                          this%pw_ref, &
                                          this%tw_ref, &
                                          this%dw_kg_ref(option%liquid_phase),&
                                          this%q_fld(option%liquid_phase) * &
                                          output_option%tconv, &
                                          this%mr_fld(option%liquid_phase) * &
                                          output_option%tconv, &
                                          this%q_fld(option%oil_phase) * &
                                          output_option%tconv, &
                                          this%mr_fld(option%oil_phase) * &
                                          output_option%tconv, &
                                          vol_WOR 


end subroutine TOilImsOilProdOutput

! ************************************************************************** !

function CreateTOilImsWell(well_spec,option)
  ! 
  ! Create a toil ims well object based on type specified in the well_spec
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  use Option_module
  use WellSpec_Base_class

  implicit none

  class(well_spec_base_type), pointer :: well_spec
  type(option_type) :: option

  class(well_flow_energy_type), pointer :: CreateTOilImsWell

  class(well_toil_ims_wat_inj_type), pointer :: well_toil_ims_wat_inj
  class(well_toil_ims_oil_prod_type), pointer :: well_toil_ims_oil_prod


  select case(well_spec%itype)
    case( WATER_INJ_WELL_TYPE )
      allocate(well_toil_ims_wat_inj);
      CreateTOilImsWell => well_toil_ims_wat_inj;
    case( OIL_PROD_WELL_TYPE )
      allocate(well_toil_ims_oil_prod);
      CreateTOilImsWell => well_toil_ims_oil_prod

    ! need to add water producer and oil injector
    case default
      option%io_buffer = 'Well type not recognize in CreateTOilImsWell'
      call printErrMsg(option)
  end select

  !initialise different well components 
  call WellBaseInit(CreateTOilImsWell,well_spec,option);
  call WellFlowInit(CreateTOilImsWell,option);
  call WellFlowEnergyInit(CreateTOilImsWell,option);

  !anything to initialise at injector/producer level?

end function CreateTOilImsWell

! ************************************************************************** !
subroutine TOilImsWatInjExplRes(this,iconn,ss_flow_vol_flux,isothermal, &
                                ghosted_id, dof,option,res)
  ! 
  ! Compute residual term for a TOilIms Water injector
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 06/06/16
  ! 
  use PM_TOilIms_Aux_module
  use EOS_Water_module
  use Option_module
  
  implicit none

  class(well_toil_ims_wat_inj_type) :: this
  PetscInt :: iconn
  PetscBool :: isothermal
  PetscInt :: ghosted_id, dof
  type(option_type) :: option
  !PetscReal :: Res(1:option%nflowdof)
  !PetscReal :: ss_flow_vol_flux(1:option%nphase)
  PetscReal :: Res(:)
  PetscReal :: ss_flow_vol_flux(:)


  !why not using a pointer to avoid the copy?
  PetscReal :: dw_kg, dw_h2o_mol
  PetscReal :: dphi, vol_flux, cfact, mob, hc 
  PetscReal :: enth_src_h2o
  PetscInt :: ierr

  Res = 0.0d0
  vol_flux = 0.0d0

  hc = this%conn_h(iconn)
  cfact = this%conn_factors(iconn)

  mob = this%ConnMob(this%flow_auxvars(dof,ghosted_id)%mobility, &
                                       option%liquid_phase)
  !mob = 1754.0d0
  dphi = this%pw_ref + hc - & 
            this%flow_auxvars(dof,ghosted_id)%pres(option%liquid_phase)

  ! it is assumed that the temperature is uniform throughout the well
  call EOSWaterDensity(this%tw_ref,this%pw_ref+hc, &
                       dw_kg,dw_h2o_mol,ierr) 

  if (.not.isothermal) then  
    call EOSWaterEnthalpy(this%tw_ref,this%pw_ref+hc,enth_src_h2o,ierr)     
    !enth_src_h2o = enth_src_h2o * option%scale
    enth_src_h2o = enth_src_h2o * 1.d-6 ! J/kmol -> MJ/kmol 
   end if

  if(cfact * mob > wfloweps) then

    if ( dphi < 0.0d0 .and. dof==ZERO_INTEGER ) &
      write(*,"('TOilImsWatInj reverse flow at gh = ',I5,' dp = ',e10.4)") &
            ghosted_id, dphi

    !         m^3 * 1/(Pa.s) * Pa = m^3/s
    vol_flux = cfact * mob * dphi
    !vol_flux = 0.00015
    ss_flow_vol_flux(option%liquid_phase) = vol_flux
    !no cross-flow allowed with this model
    ss_flow_vol_flux(option%oil_phase) = 0.0d0
    ! H2O equation       !m^3/s * kmol/m^3 = Kmol/sec 
    Res(option%water_id) = vol_flux* dw_h2o_mol 
    ! energy equation             !m^3/s * kmol/m^3 * MJ/Kmol = MJ/s
    if (.not.isothermal) Res(3) = vol_flux*dw_h2o_mol * enth_src_h2o
   end if

#ifdef WELL_DEBUG
  if ( dof==ZERO_INTEGER ) then
    write(*,*) 'ExplRes dof = ', dof
    write(*,*) 'ExplRes gh = ', ghosted_id
    write(*,"('ExplRes gh press = ',e46.40)") &
        this%flow_auxvars(dof,ghosted_id)%pres(option%liquid_phase)
    write(*,"('ExplRes sat_wat = ',e46.40)") &
       this%flow_auxvars(dof,ghosted_id)%sat(option%liquid_phase)
    write(*,"('ExplRes sat_oil = ',e46.40)") &
       this%flow_auxvars(dof,ghosted_id)%sat(option%oil_phase)
    write(*,"('ExplRes mob_wat = ',e46.40)") &
       this%flow_auxvars(dof,ghosted_id)%mobility(option%liquid_phase)
    write(*,"('ExplRes mob_oil = ',e46.40)") &
       this%flow_auxvars(dof,ghosted_id)%mobility(option%oil_phase)
    write(*,"('ExplRes mob = ',e16.10)") mob
    write(*,"('ExplRes dphi = ',e16.10)") dphi
    write(*,"('ExplRes hc = ',e46.40)") hc
    write(*,"('ExplRes conn_den_kg = ',e46.40)") this%conn_den_kg(iconn)
    write(*,"('ExplRes pw_ref = ',e46.40)") this%pw_ref 
    write(*,"('ExplRes vol_flux = ',e10.4)") vol_flux
    write(*,"('ExplRes dw_h2o_mol = ',e10.4)") dw_h2o_mol
    write(*,"('ExplRes Res(water_id) = ',e10.4)") Res(option%water_id)
  end if
#endif

end subroutine TOilImsWatInjExplRes

! ************************************************************************** !

subroutine TOilImsOilProdExplRes(this,iconn,ss_flow_vol_flux,isothermal, &
                                ghosted_id, dof,option,res)
  ! 
  ! Compute residual term for a TOilIms Oil Producer
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 06/12/16
  ! 
  use PM_TOilIms_Aux_module
  use Option_module
  
  implicit none

  class(well_toil_ims_oil_prod_type) :: this
  PetscInt :: iconn
  PetscBool :: isothermal
  PetscInt :: ghosted_id, dof
  type(option_type) :: option
  !PetscReal :: Res(1:option%nflowdof)
  !PetscReal :: ss_flow_vol_flux(1:option%nphase)
  PetscReal :: Res(:)
  PetscReal :: ss_flow_vol_flux(:)


  call TOilImsProducerExplRes(this,iconn,ss_flow_vol_flux,isothermal, &
                               ghosted_id, dof,option,res)

end subroutine TOilImsOilProdExplRes

! ************************************************************************** !

subroutine TOilImsProducerExplRes(this,iconn,ss_flow_vol_flux,isothermal, &
                                  ghosted_id,dof,option,res)

  ! 
  ! Compute residual term for a TOilIms Producers
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 06/12/16
  ! 

  use PM_TOilIms_Aux_module
  use Option_module
  use EOS_Oil_module
  use EOS_Water_module

  implicit none

  !using well_flow_energy_type to pass whaterver type of producer 
  class(well_flow_energy_type) :: this
  PetscInt :: iconn
  PetscBool :: isothermal
  PetscInt :: ghosted_id, dof, ierr
  type(option_type) :: option
  !PetscReal :: Res(1:option%nflowdof)
  !PetscReal :: ss_flow_vol_flux(1:option%nphase)
  PetscReal :: Res(:)
  PetscReal :: ss_flow_vol_flux(:)



  PetscInt :: i_ph
  PetscReal :: dphi, vol_flux, cfact, mob, hc, temp
  PetscReal :: well_oil_mol_den, mol_den_av, dw_kg
  PetscReal :: phase_mol_den(1:option%nphase)
  PetscReal :: phase_ent(1:option%nphase)

  hc = this%conn_h(iconn)
  temp = this%conn_temp(iconn)
  cfact = this%conn_factors(iconn)
  
  Res = 0.d0
  vol_flux = 0.0d0 
  ss_flow_vol_flux = 0.0d0
  phase_mol_den = 0.0d0 
  phase_ent = 0.0d0 
  dw_kg = 0.d0

  do i_ph = 1, option%nphase
    !pressure gradient positive for flow entering the well

    dphi = this%flow_auxvars(dof,ghosted_id)%pres(i_ph) - &
           this%pw_ref - hc 

    !upwind for den_mol 
    if (dphi >= 0.0d0) then
      phase_mol_den(i_ph) = this%flow_auxvars(dof,ghosted_id)%den(i_ph)
    else if (dphi < 0.0d0 ) then
      if (i_ph == option%liquid_phase) then
         call EOSWaterDensity(temp,this%pw_ref+hc, &
                              dw_kg,phase_mol_den(i_ph),ierr) 
      else if (i_ph == option%oil_phase) then 
        call EOSOilDensity(temp,this%pw_ref+hc,phase_mol_den(i_ph),ierr)
      end if
    end if

    !upwind for enthalpy
    if (.not.isothermal) then
      if (dphi >= 0.0d0) then
        phase_ent(i_ph) = this%flow_energy_auxvars(dof,ghosted_id)%H(i_ph)
      else if (dphi < 0.0d0) then 
        if (i_ph == option%liquid_phase) then  
          call EOSWaterEnthalpy(temp,this%pw_ref+hc,phase_ent(i_ph),ierr)
        else if (i_ph == option%oil_phase) then
          call EOSOilEnthalpy(temp,this%pw_ref+hc,phase_ent(i_ph),ierr)
        end if
        phase_ent = phase_ent * 1.d-6 ! J/kmol -> whatever units
      end if 
    end if
    
    mob = this%ConnMob(this%flow_auxvars(dof,ghosted_id)%mobility,i_ph)

    !if ( dabs(dphi) < 1.d-5 ) dphi = 0.0d0 !cut off noise (Pa)
    !if ( dphi < 0.0d0 ) &
    !  write(*,"('TOilImsWellProd reverse flow at gh = ',I5,' dp = ',e10.4)") &
    !  ghosted_id, dphi
     
    !if(cfact * mob > wfloweps) then
    if( mob > wfloweps) then
      !!         m^3 * 1/(Pa.s) * Pa = m^3/s 
      vol_flux = cfact * mob * dphi

      if ( vol_flux < 0.0d0 .and. dof==ZERO_INTEGER ) &
       write(*,"('TOilImsWellProd reverse flow at gh = ',I5,' dp = ',e10.4)") &
        ghosted_id, dphi
   
      !stopping reversing flows to occur - they cannot be handled with this model
      if ( vol_flux < wfloweps ) cycle
      !if( dabs(vol_flux) > 1.d-10 ) then !try to cut som noise
      !if( dabs(dphi/this%pw_ref) > 1.d-7 ) then !try to cut som noise
        ! the minus sign indicate component fluxes out the reservoir
        ss_flow_vol_flux(i_ph) = -1.d0 * vol_flux

        !call EOSOilDensity(temp,this%pw_ref+hc,well_oil_mol_den,ierr)  
        !mol_den_av = ( this%flow_auxvars(dof,ghosted_id)%den(i_ph) + &
        !               well_oil_mol_den ) * 0.5d0 
        !obs: can use Res(i_ph) here because i_ph conicide with equation indices
        !the minus sign indicate component fluxes out the reservoir
        Res(i_ph) = - vol_flux * phase_mol_den(i_ph)
        !Res(i_ph) = 0.d0
        if (.not.isothermal) then
          Res(TOIL_IMS_ENERGY_EQUATION_INDEX) = &
                Res(TOIL_IMS_ENERGY_EQUATION_INDEX) - &
                vol_flux * phase_mol_den(i_ph) * phase_ent(i_ph)   
                !vol_flux * this%flow_auxvars(dof,ghosted_id)%den(i_ph) * &
                !this%flow_energy_auxvars(dof,ghosted_id)%H(i_ph)          
        end if
      !end if
    end if  

#ifdef WELL_DEBUG
  if ( dof==ZERO_INTEGER ) then
    write(*,*) 'ExplRes dof = ', dof
    write(*,*) 'ExplRes gh = ', ghosted_id
    write(*,*) 'ExplRes i_ph = ', i_ph
    write(*,"('ExplRes cell press = ',e26.20)") &
        this%flow_auxvars(dof,ghosted_id)%pres(i_ph)
    write(*,"('ExplRes mob = ',e46.40)") mob
    write(*,"('ExplRes sat_i_ph = ',e46.40)") &
       this%flow_auxvars(dof,ghosted_id)%sat(i_ph)
    write(*,"('ExplRes dphi = ',e46.40)") dphi
    write(*,"('ExplRes dp/Pw = ',e46.40)") dphi/this%pw_ref
    write(*,"('ExplRes hc = ',e26.20)") hc
    write(*,"('ExplRes pw_ref = ',e16.10)") this%pw_ref 
    write(*,"('ExplRes vol_flux = ',e26.20)") vol_flux
    write(*,"('ExplRes Res(i_ph) = ',e46.40)") Res(i_ph)
  end if
#endif

  end do


end subroutine TOilImsProducerExplRes

! ************************************************************************** !

#endif  
end module Well_TOilIms_class
!end of WELL_CLASS

