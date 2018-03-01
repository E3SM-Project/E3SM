module BgcCentCnpIndexType

  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_varcon    , only : spinup_state => bspinup_state
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  integer, parameter :: loc_name_len=64

  type, private :: list_t
    character(len=loc_name_len) :: name
    type(list_t), pointer :: next => null()
  end type list_t

  type, public :: centurybgc_index_type
     integer           :: nom_pools                              !not include coarse wood debris

     integer           :: nom_tot_elms
     integer           :: lit1, lit1_dek_reac
     integer           :: lit2, lit2_dek_reac
     integer           :: lit3, lit3_dek_reac
     integer           :: som1, som1_dek_reac
     integer           :: som2, som2_dek_reac
     integer           :: som3, som3_dek_reac
     integer           :: cwd,  cwd_dek_reac
     integer           :: lwd,  lwd_dek_reac
     integer           :: fwd,  fwd_dek_reac
     integer           :: litr_beg, litr_end  !litr group
     integer           :: wood_beg, wood_end  !wood group
     integer           :: som_beg,  som_end   !som group
     integer           :: dom_beg,  dom_end   !dom group
     integer           :: Bm_beg,  Bm_end   !dom group

     integer           :: c_loc
     integer           :: n_loc
     integer           :: p_loc
     integer           :: c13_loc = 0
     integer           :: c14_loc = 0
     integer           :: nelms                                  !number of chemical elements in an om pool
                                                                 !reactive primary variables
     integer           :: lid_nh4, lid_nh4_nit_reac              !local position of nh4 in the state variable vector
     integer           :: lid_no3, lid_no3_den_reac              !local position of no3 in the state variable vector
     integer           :: lid_plant_minn_nh4, lid_plant_minn_nh4_up_reac !local position of plant uptake of mineral nitrogen NH4 in the state variable vector
     integer           :: lid_plant_minn_no3, lid_plant_minn_no3_up_reac !
     integer           :: lid_plant_minp, lid_plant_minp_up_reac !local position of plant uptake of mineral P in the state variable vector
     integer           :: lid_minp_soluble, lid_minp_soluble_to_secp_reac    !conversation of adsorbed into secondary phase
     integer           :: lid_minp_secondary,lid_minp_secondary_to_sol_occ_reac   !local position of secondary P in the state variable vector

     integer           :: lid_minp_occlude      !local position of occluded P in the state variable vector

     integer           :: lid_autr_rt, lid_autr_rt_reac             !root autotrophic respiration

                                                                 !non reactive primary variables
     integer           :: lid_ar, lid_ar_aren_reac               !local position of ar in the state variable vector
     integer           :: lid_ch4, lid_ch4_aren_reac             !nonreactive primary variables

                                                                 !secondary variables
     integer           :: lid_o2,  lid_o2_aren_reac              !local position of o2 in the state variable vector
     integer           :: lid_co2, lid_co2_aren_reac             !local position of co2 in the state variable vector
     integer           :: lid_n2,  lid_n2_aren_reac
     integer           :: lid_n2o, lid_n2o_aren_reac
                                                                 !diagnostic variables
     integer           :: lid_n2o_nit                            !n2o production from nitrification, used to for mass balance book keeping
     integer           :: lid_co2_hr                             !co2 production from heterotrophic respiration
     integer           :: lid_c13_co2, lid_c13_co2_aren_reac
     integer           :: lid_c14_co2, lid_c14_co2_aren_reac
     integer           :: lid_no3_den                            !no3 consumption due to denitrification
     integer           :: lid_minn_nh4_immob                     !net mineral NH4 immobilization for decomposition
     integer           :: lid_minn_no3_immob                     !net mineral NO3 immobilization for decomposition
     integer           :: lid_nh4_nit
     integer           :: lid_minp_secondary_trc
     integer           :: lid_minp_occlude_trc
                                                                 !aerechyma transport, diagnostic efflux
     integer           :: lid_minp_immob                         !net P immobilization by aerobic decomposer

     integer           :: lid_ar_paere
     integer           :: lid_n2_paere
     integer           :: lid_o2_paere
     integer           :: lid_co2_paere
     integer           :: lid_c13_co2_paere
     integer           :: lid_c14_co2_paere
     integer           :: lid_ch4_paere
     integer           :: lid_n2o_paere

     integer           :: lid_nh4_supp
     integer           :: nprimvars                              !total number of primary variables
     integer           :: nstvars                                !number of equations for the state variabile vector
     integer           :: nreactions                             !seven decomposition pathways plus nitrification, denitrification and plant immobilization

     integer , pointer :: primvarid(:)   => null()
     logical , pointer :: is_aerobic_reac(:)=> null()

     integer, pointer :: lid_plant_minn_no3_pft(:)=> null()
     integer, pointer :: lid_plant_minn_nh4_pft(:)=> null()
     integer, pointer :: lid_plant_minp_pft(:)=> null()
     logical :: debug
     character(len=loc_name_len), allocatable :: varnames(:)
   contains
     procedure, public  :: Init
     procedure, private :: InitPars
     procedure, private :: InitAllocate
     procedure, private :: set_primvar_reac_ids
  end type centurybgc_index_type

  contains
  !-----------------------------------------------------------------------
  subroutine list_init(self, name)
  implicit none
  type(list_t), pointer :: self
  character(len=*), intent(in) :: name

  allocate(self)
  nullify(self%next)
  write(self%name,'(A)')trim(name)

  end subroutine list_init
  !-----------------------------------------------------------------------
  subroutine list_insert(self, name)

  implicit none
  type(list_t), pointer :: self
  character(len=*), intent(in) :: name
  type(list_t), pointer :: next

  allocate(next)
  write(next%name,'(A)')trim(name)
  next%next=> self
  self => next

  end subroutine list_insert
  !-----------------------------------------------------------------------
  subroutine list_free(self)
  implicit none
  type(list_t), pointer :: self
  type(list_t), pointer :: current
  type(list_t), pointer :: elem

  elem => self
  do while(associated(elem))
    current => elem
    elem => current%next
    deallocate(current)
  enddo
  end subroutine list_free
  !-----------------------------------------------------------------------
  function list_next(self)result(next)

  implicit none
  type(list_t), pointer :: self
  type(list_t), pointer :: next

  next => self%next
  end function list_next

  !-------------------------------------------------------------------------------
  subroutine copy_name(num_names, list_name, outnames)

  implicit none
  integer, intent(in) :: num_names
  type(list_t), pointer :: list_name
  character(len=loc_name_len), intent(out) :: outnames(num_names)

  type(list_t), pointer :: next
  integer :: jj
  next => list_name
  do jj = num_names, 1, -1
    write(outnames(jj),'(A)')trim(next%name)
    next=>list_next(next)
  enddo
  end subroutine copy_name
  !-------------------------------------------------------------------------------
  subroutine add_ompool_name(list_name, prefix, use_c13, use_c14, do_init)
  implicit none
  type(list_t), pointer :: list_name
  character(len=*), intent(in) :: prefix
  logical, intent(in) :: use_c13, use_c14
  logical, intent(in) :: do_init

  if(do_init)then
    call list_init(list_name, trim(prefix)//'_c')
  else
    call list_insert(list_name, trim(prefix)//'_c')
  endif
  call list_insert(list_name, trim(prefix)//'_n')
  call list_insert(list_name, trim(prefix)//'_p')
  if(use_c13)call list_insert(list_name, trim(prefix)//'_c13')
  if(use_c14)call list_insert(list_name, trim(prefix)//'_c14')

  end subroutine add_ompool_name

  !-------------------------------------------------------------------------------
  subroutine Init(this, use_c13, use_c14, maxpft)
    !
    ! DESCRIPTION:
    ! Initialize centurybgc type
    ! !USES:
  implicit none
  ! !ARGUMENTS:
  class(centurybgc_index_type), intent(inout) :: this
  logical, intent(in) :: use_c13
  logical, intent(in) :: use_c14
  integer, optional, intent(in) :: maxpft

  ! !LOCAL VARIABLES:
  integer :: maxpft_loc

  maxpft_loc = 0
  if(present(maxpft))maxpft_loc=maxpft
  call this%InitPars(maxpft_loc, use_c14, use_c13)

  call this%InitAllocate()

  call this%set_primvar_reac_ids()

  this%debug = .false.
  end subroutine Init
  !-------------------------------------------------------------------------------

  subroutine InitPars(this, maxpft, use_c14, use_c13)
    !
    ! !DESCRIPTION:
    !  describe the layout of the stoichiometric matrix for the reactions
    !           r{1} r{2} r{3} r{4} ... r{n}
    ! s{1}
    ! s{2}
    ! s{3}
    ! s{4}
    ! ...
    ! s{n}
    ! s{n+1}  nonreactive primary variables
    ! s{n+2}
    ! ...
    ! s{m}
    ! s{m+1} diagnostic variables
    ! s{p}
    ! each reaction is associated with a primary species, the secondary species follows after primary species
    ! for the century model, the primary species are seven om pools and nh4, no3 and plant nitrogen
    !
    ! !USES:
    use MathfuncMod   , only : addone, countelm
    use betr_utils    , only : num2str
    use betr_constants, only : betr_string_length_long
    implicit none
    class(centurybgc_index_type) :: this
    integer, intent(in) :: maxpft
    logical, intent(in) :: use_c13
    logical, intent(in) :: use_c14
    ! !LOCAL VARIABLES:
    integer :: itemp
    integer :: ireac   !counter of reactions
    integer :: itemp1
    integer :: ielem
    integer :: jj
    type(list_t), pointer :: list_name => null()
    character(len=loc_name_len) :: postfix

    itemp = 0
    ireac = 0
    ielem= 0

    this%c_loc = addone(ielem)
    this%n_loc = addone(ielem)
    this%p_loc = addone(ielem)
    if(use_c13)then
      this%c13_loc = addone(ielem)
    endif
    if(use_c14)then
      this%c14_loc = addone(ielem)
    endif
    this%nelms = ielem   !carbon and nitrogen

    !litter group
    this%litr_beg=1
    this%lit1 = addone(itemp); this%lit1_dek_reac = addone(ireac)
    call add_ompool_name(list_name, 'lit1', use_c13, use_c14, do_init=.true.)
    this%lit2 = addone(itemp); this%lit2_dek_reac = addone(ireac)
    call add_ompool_name(list_name, 'lit2', use_c13, use_c14, do_init=.false.)
    this%lit3 = addone(itemp); this%lit3_dek_reac = addone(ireac)
    call add_ompool_name(list_name, 'lit3', use_c13, use_c14, do_init=.false.)
    this%litr_end=this%litr_beg-1+3*this%nelms

    !woody group
    this%wood_beg=this%litr_end+1
    this%cwd  = addone(itemp); this%cwd_dek_reac  = addone(ireac)
    call add_ompool_name(list_name, 'cwd', use_c13, use_c14, do_init=.false.)
    this%lwd  = addone(itemp); this%lwd_dek_reac  = addone(ireac)
    call add_ompool_name(list_name, 'lwd', use_c13, use_c14, do_init=.false.)
    this%fwd  = addone(itemp); this%fwd_dek_reac  = addone(ireac)
    call add_ompool_name(list_name, 'fwd', use_c13, use_c14, do_init=.false.)
    this%wood_end=this%wood_beg-1+3*this%nelms

    !som group
    this%Bm_beg=this%wood_end+1
    this%som1 = addone(itemp); this%som1_dek_reac = addone(ireac)
    call add_ompool_name(list_name, 'Bm', use_c13, use_c14, do_init=.false.)
    this%Bm_end=this%Bm_beg-1+this%nelms

    this%som_beg=this%Bm_end+1
    this%som3 = addone(itemp); this%som3_dek_reac = addone(ireac)
    call add_ompool_name(list_name, 'som3', use_c13, use_c14, do_init=.false.)
    this%som_end=this%som_beg-1+this%nelms

    !dom group
    this%dom_beg=this%som_end+1
    this%som2 = addone(itemp); this%som2_dek_reac = addone(ireac)  !put som2 at the end because it is defined as dom
    call add_ompool_name(list_name, 'som2', use_c13, use_c14, do_init=.false.)
    this%dom_end=this%dom_beg-1+this%nelms

    this%nom_pools = (countelm(this%litr_beg, this%litr_end)+&
       countelm(this%wood_beg,this%wood_end) + &
       countelm(this%som_beg,this%som_end) + &
       countelm(this%Bm_beg,this%Bm_end) + &
       countelm(this%dom_beg,this%dom_end))/this%nelms   !include coarse wood debris

    itemp               = this%nom_pools*this%nelms

    this%nom_tot_elms    = itemp
    this%lid_nh4        = addone(itemp); this%lid_nh4_nit_reac = addone(ireac)       !this is also used to indicate the nitrification reaction
    call list_insert(list_name, 'nh4')
    this%lid_no3        = addone(itemp); this%lid_no3_den_reac = addone(ireac)       !this is also used to indicate the denitrification reaction
    call list_insert(list_name, 'no3')

    this%lid_minp_soluble=addone(itemp); this%lid_minp_soluble_to_secp_reac = addone(ireac)
    call list_insert(list_name, 'minp_soluble')

    this%lid_minp_secondary = addone(itemp); this%lid_minp_secondary_to_sol_occ_reac=addone(ireac)
    call list_insert(list_name, 'minp_secondary')

    this%lid_minp_occlude = addone(itemp);
    call list_insert(list_name, 'minp_occlude')

    this%lid_plant_minn_nh4_up_reac = addone(ireac)
    this%lid_plant_minn_no3_up_reac = addone(ireac)
    this%lid_plant_minp_up_reac = addone(ireac)
    this%lid_autr_rt_reac = addone(ireac)

    !non-reactive primary variables
    this%lid_ch4        = addone(itemp);call list_insert(list_name, 'ch4')
    this%lid_ar         = addone(itemp);call list_insert(list_name, 'ar')

    !second primary variables
    this%lid_o2         = addone(itemp);call list_insert(list_name, 'o2')

    this%lid_co2        = addone(itemp);call list_insert(list_name, 'co2')

    if(use_c13)then
      this%lid_c13_co2        = addone(itemp);call list_insert(list_name, 'co2_c13')
    endif
    if(use_c14)then
      this%lid_c14_co2        = addone(itemp);call list_insert(list_name, 'co2_c14')
    endif

    this%lid_n2o        = addone(itemp);call list_insert(list_name, 'n2o')
    this%lid_n2         = addone(itemp);call list_insert(list_name, 'n2')
    this%nprimvars      = itemp     !primary state variables 14 + 6

    !diagnostic variables
    this%lid_plant_minn_nh4 = addone(itemp)  !this is used to indicate plant mineral nitrogen uptake
    call list_insert(list_name, 'plant_minn_nh4')

    this%lid_plant_minn_no3 = addone(itemp)  !this is used to indicate plant mineral nitrogen uptake
    call list_insert(list_name, 'plant_minn_no3')

    this%lid_plant_minp = addone(itemp)
    call list_insert(list_name, 'plant_minp')

    this%lid_autr_rt      = addone(itemp)           !this is used to indicate plant autotrophic root respiration
    call list_insert(list_name, 'autr_rt')

    this%lid_o2_aren_reac  = addone(ireac)

    this%lid_n2o_nit        = addone(itemp); call list_insert(list_name, 'n2o_nit')

    this%lid_co2_hr         = addone(itemp); call list_insert(list_name, 'co2_hr')

    this%lid_no3_den        = addone(itemp); call list_insert(list_name, 'no3_den')

    this%lid_minn_nh4_immob = addone(itemp); call list_insert(list_name, 'minn_nh4_immob')

    this%lid_minn_no3_immob = addone(itemp); call list_insert(list_name, 'minn_no3_immob')

    this%lid_minp_immob = addone(itemp); call list_insert(list_name, 'minp_immob')

    this%lid_nh4_nit        = addone(itemp); call list_insert(list_name, 'nh4_nit')

    !aerechyma transport
    this%lid_o2_paere   = addone(itemp); call list_insert(list_name, 'o2_paere')
    if ( spinup_state /= 1 ) then
       this%lid_ar_paere   = addone(itemp);  this%lid_ar_aren_reac  = addone(ireac)   !
       call list_insert(list_name, 'ar_paere')

       this%lid_n2_paere   = addone(itemp);  this%lid_n2_aren_reac  = addone(ireac)   !
       call list_insert(list_name, 'n2_paere')

       this%lid_co2_paere  = addone(itemp);  this%lid_co2_aren_reac = addone(ireac)   !
       call list_insert(list_name, 'co2_paere')

       if(use_c13)then
         this%lid_c13_co2_paere  = addone(itemp);  this%lid_c13_co2_aren_reac = addone(ireac)   !
         call list_insert(list_name, 'c13_co2_paere')
       endif

       if(use_c14)then
         this%lid_c14_co2_paere  = addone(itemp);  this%lid_c14_co2_aren_reac = addone(ireac)   !
         call list_insert(list_name, 'c14_co2_paere')
       endif

       this%lid_ch4_paere  = addone(itemp);  this%lid_ch4_aren_reac = addone(ireac)   !
       call list_insert(list_name, 'ch4_paere')

       this%lid_n2o_paere  = addone(itemp);  this%lid_n2o_aren_reac = addone(ireac)   !
       call list_insert(list_name, 'n2o_paere')
    endif

    if(maxpft>0)then
      allocate(this%lid_plant_minn_no3_pft(maxpft));
      allocate(this%lid_plant_minn_nh4_pft(maxpft));
      allocate(this%lid_plant_minp_pft(maxpft));
      do jj = 1, maxpft
        this%lid_plant_minn_no3_pft(jj) = addone(itemp)
        postfix = num2str(jj,'(I2.2)')
        call list_insert(list_name, 'plant_minn_no3_pft_'//trim(postfix))
        this%lid_plant_minn_nh4_pft(jj) = addone(itemp)
        call list_insert(list_name, 'plant_minn_nh4_pft_'//trim(postfix))
        this%lid_plant_minp_pft(jj)     = addone(itemp)
        call list_insert(list_name, 'plant_minp_pft_'//trim(postfix))
      enddo
    endif

    this%nstvars          = itemp          !totally 14+32 state variables
    this%nreactions = ireac            !seven decomposition pathways plus root auto respiration, nitrification, denitrification and plant immobilization

    allocate(this%primvarid(ireac)); this%primvarid(:) = -1
    allocate(this%is_aerobic_reac(ireac)); this%is_aerobic_reac(:)=.false.

    allocate(this%varnames(this%nstvars))
    call copy_name(this%nstvars, list_name, this%varnames(1:this%nstvars))
    call list_free(list_name)

  end subroutine InitPars
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this)
    !
    ! !DESCRIPTION:
    ! memory allocation for the data type specified by this
    !
  implicit none
    ! !ARGUMENTS:
  class(centurybgc_index_type), intent(inout) :: this



  end subroutine InitAllocate
  !-------------------------------------------------------------------------------
  subroutine set_primvar_reac_ids(this)

  implicit none
  class(centurybgc_index_type), intent(inout)  :: this

  integer :: reac

  associate(                                      &
    lit1  => this%lit1                       ,    &
    lit2  => this%lit2                       ,    &
    lit3  => this%lit3                       ,    &
    cwd  => this%cwd                         ,    &
    som1  => this%som1                       ,    &
    som2  => this%som2                       ,    &
    som3  => this%som3                       ,    &
    nelms => this%nelms                      ,    &
    c_loc => this%c_loc                      ,    &
    n_loc => this%n_loc                      ,    &
    p_loc => this%p_loc                      ,    &
    lid_nh4=> this%lid_nh4                   ,    &
    lid_no3 => this%lid_no3                  ,    &
    lid_o2 => this%lid_o2                    ,    &
    lid_minp_soluble=> this%lid_minp_soluble ,    &
    lid_minp_secondary => this%lid_minp_secondary &
  )
  !reaction1, lit1 -> s1
  reac=this%lit1_dek_reac;     this%primvarid(reac) = (lit1-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.

  !reaction 2, lit2 -> s1
  reac =this%lit2_dek_reac;   this%primvarid(reac) = (lit2-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.

  !reaction 3, lit3->s2
  reac =this%lit3_dek_reac; this%primvarid(reac) = (lit3-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.

  !reaction 4, SOM1 -> f1*SOM2 + f2*SOm3
  reac =this%som1_dek_reac; this%primvarid(reac) = (som1-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.

  !reaction 5, som2->som1, som3
  reac =this%som2_dek_reac;  this%primvarid(reac) = (som2-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.

  !reaction 6, s3-> s1
  reac = this%som3_dek_reac;  this%primvarid(reac) = (som3-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.

  !reaction 7, cwd -> lit1 + lit2
  reac = this%cwd_dek_reac; this%primvarid(reac) = (cwd-1)*nelms+c_loc
  !x is_aerobic_reac(reac) = .true.

  !reaction 8, nitrification
  reac = this%lid_nh4_nit_reac; this%primvarid(reac) = this%lid_nh4
  !x is_aerobic_reac(reac) = .true.

  !reaction 9, denitrification
  reac = this%lid_no3_den_reac; this%primvarid(reac) = lid_no3

  !reaction 10, inorganic P non-equilibrium adsorption
  !P_solution -> p_secondary
  reac = this%lid_minp_soluble_to_secp_reac; this%primvarid(reac) = lid_minp_soluble

  !reaction 11, inorganic P non-equilibrium desorption
  ! p_secondary -> P_solution
  reac = this%lid_minp_secondary_to_sol_occ_reac; this%primvarid(reac) = lid_minp_secondary

  !reaction 12, plant mineral nitrogen nh4 uptake
  reac = this%lid_plant_minn_nh4_up_reac; this%primvarid(reac) = lid_nh4

  !reaction 13, plant mineral nitrogen no3 uptake
  reac = this%lid_plant_minn_no3_up_reac; this%primvarid(reac) = lid_no3

  !reaction 14, plant mineral phosphorus uptake
  reac = this%lid_plant_minp_up_reac; this%primvarid(reac) = lid_minp_soluble

  !reaction 15, autotrophic respiration, ar + o2 -> co2
  !x reac = lid_ar_rt_reac; is_aerobic_reac(reac) = .true.

  !reaction 15, o2 transport through arenchyma
  reac = this%lid_o2_aren_reac; this%primvarid(reac) = lid_o2
  end associate
  end subroutine set_primvar_reac_ids
end module BgcCentCnpIndexType
