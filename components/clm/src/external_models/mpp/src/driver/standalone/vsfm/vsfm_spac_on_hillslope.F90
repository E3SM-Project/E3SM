!------------------------------------------------------------------------
module spac_component
  
#include <petsc/finclude/petsc.h>
  
  implicit none

  type, public :: spac_component_mesh_type
     PetscInt            :: ncells
     PetscInt  , pointer :: id(:,:,:)
     PetscReal , pointer :: xc(:)                          !
     PetscReal , pointer :: yc(:)                          !
     PetscReal , pointer :: zc(:)                          !
     PetscReal , pointer :: dx(:)                          !
     PetscReal , pointer :: dy(:)                          !
     PetscReal , pointer :: dz(:)                          !
     PetscReal , pointer :: area(:)                        !
     PetscReal , pointer :: vol(:)                         !
     PetscInt  , pointer :: filter(:)                      !
   contains
     procedure, public :: Init     => MeshInit
     procedure, public :: Copy     => MeshCopy
     procedure, public :: AddToMPP => MeshAddToMPP
  end type spac_component_mesh_type

  type, public :: spac_component_pp_type
     PetscInt            :: ncells
     PetscReal , pointer :: por(:)
     PetscReal , pointer :: perm(:)
     PetscReal , pointer :: lambda(:)
     PetscReal , pointer :: alpha(:)
     PetscReal , pointer :: eff_porosity(:)
     PetscReal , pointer :: residual_sat(:)
     PetscInt  , pointer :: satfunc_type(:)
     PetscInt  , pointer :: relperm_type(:)
     PetscReal , pointer :: relperm_param_1(:)
     PetscReal , pointer :: relperm_param_2(:)
   contains
     procedure, public :: Init => PPInit
     procedure, public :: Copy => PPCopy
  end type spac_component_pp_type

  type, public :: spac_component_conn_type

     PetscInt            :: nconn

     PetscInt            :: up_eqn_id
     PetscInt            :: dn_eqn_id
     character(len=256)  :: up_eqn_name
     character(len=256)  :: dn_eqn_name

     PetscInt  , pointer :: id_up(:)   !
     PetscInt  , pointer :: id_dn(:)   !
     PetscReal , pointer :: dist_up(:) !
     PetscReal , pointer :: dist_dn(:) !
     PetscReal , pointer :: area(:)    !
     PetscInt  , pointer :: itype(:)    !

     PetscReal , pointer :: unit_vec(:,:)

     PetscInt  , pointer :: flux_type(:)    !
     PetscInt  , pointer :: cond_type(:)    !
     PetscReal , pointer :: cond     (:)    !
     PetscReal , pointer :: cond_up  (:)    !
     PetscReal , pointer :: cond_dn  (:)    !

     PetscInt  , pointer :: satparam_up_itype(:)
     PetscReal , pointer :: satparam_up_param_1(:)
     PetscReal , pointer :: satparam_up_param_2(:)
     PetscReal , pointer :: satparam_up_param_3(:)
     PetscInt  , pointer :: relperm_up_itype(:)
     PetscReal , pointer :: relperm_up_param_1(:)
     PetscReal , pointer :: relperm_up_param_2(:)
     PetscReal , pointer :: relperm_up_param_3(:)

     PetscInt  , pointer :: satparam_dn_itype(:)
     PetscReal , pointer :: satparam_dn_param_1(:)
     PetscReal , pointer :: satparam_dn_param_2(:)
     PetscReal , pointer :: satparam_dn_param_3(:)
     PetscInt  , pointer :: relperm_dn_itype(:)
     PetscReal , pointer :: relperm_dn_param_1(:)
     PetscReal , pointer :: relperm_dn_param_2(:)
     PetscReal , pointer :: relperm_dn_param_3(:)

   contains

     procedure, public :: Init => ConnInit
     procedure, public :: Copy => ConnCopy
     procedure, public :: Destroy => ConnDestroy
     procedure, public :: AddToMPPAsInternalConn
     procedure, public :: AddToMPPCouplingBC

  end type spac_component_conn_type


contains

  !------------------------------------------------------------------------
  subroutine MeshInit(this, ncells)
    !
    implicit none
    !
    class (spac_component_mesh_type) :: this
    PetscInt                         :: ncells

    this%ncells = ncells

    allocate (this%xc     (ncells)); this%xc     (:) = 0.d0
    allocate (this%yc     (ncells)); this%yc     (:) = 0.d0
    allocate (this%zc     (ncells)); this%zc     (:) = 0.d0
    allocate (this%dx     (ncells)); this%dx     (:) = 0.d0
    allocate (this%dy     (ncells)); this%dy     (:) = 0.d0
    allocate (this%dz     (ncells)); this%dz     (:) = 0.d0
    allocate (this%area   (ncells)); this%area   (:) = 0.d0
    allocate (this%filter (ncells)); this%filter (:) = 0
    allocate (this%vol    (ncells)); this%vol    (:) = 0.d0

  end subroutine MeshInit

  !------------------------------------------------------------------------
  subroutine MeshCopy(this, idx_beg, idx_end, mesh)
    !
    implicit none
    !
    class (spac_component_mesh_type) :: this
    class (spac_component_mesh_type) :: mesh
    PetscInt                         :: idx_beg
    PetscInt                         :: idx_end

    this%xc     (idx_beg:idx_end) = mesh%xc     (:)
    this%yc     (idx_beg:idx_end) = mesh%yc     (:)
    this%zc     (idx_beg:idx_end) = mesh%zc     (:)
    this%dx     (idx_beg:idx_end) = mesh%dx     (:)
    this%dy     (idx_beg:idx_end) = mesh%dy     (:)
    this%dz     (idx_beg:idx_end) = mesh%dz     (:)
    this%area   (idx_beg:idx_end) = mesh%area   (:)
    this%vol    (idx_beg:idx_end) = mesh%vol    (:)
    this%filter (idx_beg:idx_end) = mesh%filter (:)

  end subroutine MeshCopy

  !------------------------------------------------------------------------
  subroutine MeshAddToMPP(this, vsfm_mpp, imesh, mesh_name, nz)
    !
    use MultiPhysicsProbVSFM      , only : mpp_vsfm_type
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : VAR_VOLUME
    !
    implicit none
    !
    class (spac_component_mesh_type) :: this
    class (mpp_vsfm_type)            :: vsfm_mpp
    PetscInt                         :: imesh
    character(len =*)                :: mesh_name
    PetscInt                         :: ncells_ghost
    PetscInt                         :: nz
    
    ncells_ghost = 0

    write(*,*)'Adding following mesh to MPP: ' // trim(mesh_name)
    call vsfm_mpp%MeshSetName                (imesh, trim(mesh_name)          )
    call vsfm_mpp%MeshSetOrientation         (imesh, MESH_ALONG_GRAVITY       )
    call vsfm_mpp%MeshSetID                  (imesh, MESH_CLM_SOIL_COL        )
    call vsfm_mpp%MeshSetDimensions          (imesh, this%ncells, ncells_ghost, nz )

    call vsfm_mpp%MeshSetGridCellFilter      (imesh, this%filter              )

    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC     , this%xc     )
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC     , this%yc     )
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC     , this%zc     )
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX     , this%dx     )
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY     , this%dy     )
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ     , this%dz     )
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA   , this%area   )
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME , this%vol    )

  end subroutine MeshAddToMPP

  !------------------------------------------------------------------------
  subroutine PPInit(this, ncells)
    !
    implicit none
    !
    class (spac_component_pp_type) :: this
    PetscInt                       :: ncells

    this%ncells = ncells

    allocate(this%por              (ncells)); this%por             (:) = 0.d0
    allocate(this%perm             (ncells)); this%perm            (:) = 0.d0
    allocate(this%alpha            (ncells)); this%alpha           (:) = 0.d0
    allocate(this%lambda           (ncells)); this%lambda          (:) = 0.d0
    allocate(this%eff_porosity     (ncells)); this%eff_porosity    (:) = 0.d0
    allocate(this%residual_sat     (ncells)); this%residual_sat    (:) = 0.d0
    allocate(this%satfunc_type     (ncells)); this%satfunc_type    (:) = 0
    allocate(this%relperm_type     (ncells)); this%relperm_type    (:) = 0
    allocate(this%relperm_param_1  (ncells)); this%relperm_param_1 (:) = 0.d0
    allocate(this%relperm_param_2  (ncells)); this%relperm_param_2 (:) = 0.d0

  end subroutine PPInit

  !------------------------------------------------------------------------
  subroutine PPCopy(this, idx_beg, idx_end, pp)
    !
    implicit none
    !
    class (spac_component_pp_type) :: this
    PetscInt                       :: idx_beg
    PetscInt                       :: idx_end
    class (spac_component_pp_type) :: pp

    this%por             (idx_beg:idx_end) = pp%por             (:)
    this%perm            (idx_beg:idx_end) = pp%perm            (:)
    this%alpha           (idx_beg:idx_end) = pp%alpha           (:)
    this%lambda          (idx_beg:idx_end) = pp%lambda          (:)
    this%relperm_type    (idx_beg:idx_end) = pp%relperm_type    (:)
    this%relperm_param_1 (idx_beg:idx_end) = pp%relperm_param_1 (:)
    this%relperm_param_2 (idx_beg:idx_end) = pp%relperm_param_2 (:)
    this%residual_sat    (idx_beg:idx_end) = pp%residual_sat    (:)
    this%satfunc_type    (idx_beg:idx_end) = pp%satfunc_type    (:)

  end subroutine PPCopy

  !------------------------------------------------------------------------
  subroutine ConnInit(this, nconn)
    !
    implicit none
    !
    class (spac_component_conn_type) :: this
    PetscInt                         :: nconn

    this%nconn = nconn

    allocate (this%id_up               (nconn)); this%id_up               (:) = 0
    allocate (this%id_dn               (nconn)); this%id_dn               (:) = 0
    allocate (this%dist_up             (nconn)); this%dist_up             (:) = 0.d0
    allocate (this%dist_dn             (nconn)); this%dist_dn             (:) = 0.d0
    allocate (this%area                (nconn)); this%area                (:) = 0.d0
    allocate (this%itype               (nconn)); this%itype               (:) = 0

    allocate (this%unit_vec            (nconn,3)); this%unit_vec          (:,:) = 0.d0

    allocate (this%flux_type           (nconn)); this%flux_type           (:) = 0
    allocate (this%cond_type           (nconn)); this%cond_type           (:) = 0
    allocate (this%cond                (nconn)); this%cond                (:) = 0.d0
    allocate (this%cond_up             (nconn)); this%cond_up             (:) = 0.d0
    allocate (this%cond_dn             (nconn)); this%cond_dn             (:) = 0.d0

    allocate (this%satparam_up_itype   (nconn)); this%satparam_up_itype   (:) = 0
    allocate (this%satparam_up_param_1 (nconn)); this%satparam_up_param_1 (:) = 0.d0
    allocate (this%satparam_up_param_2 (nconn)); this%satparam_up_param_2 (:) = 0.d0
    allocate (this%satparam_up_param_3 (nconn)); this%satparam_up_param_3 (:) = 0.d0
    allocate (this%relperm_up_itype    (nconn)); this%relperm_up_itype    (:) = 0
    allocate (this%relperm_up_param_1  (nconn)); this%relperm_up_param_1  (:) = 0.d0
    allocate (this%relperm_up_param_2  (nconn)); this%relperm_up_param_2  (:) = 0.d0
    allocate (this%relperm_up_param_3  (nconn)); this%relperm_up_param_3  (:) = 0.d0

    allocate (this%satparam_dn_itype   (nconn)); this%satparam_dn_itype   (:) = 0
    allocate (this%satparam_dn_param_1 (nconn)); this%satparam_dn_param_1 (:) = 0.d0
    allocate (this%satparam_dn_param_2 (nconn)); this%satparam_dn_param_2 (:) = 0.d0
    allocate (this%satparam_dn_param_3 (nconn)); this%satparam_dn_param_3 (:) = 0.d0
    allocate (this%relperm_dn_itype    (nconn)); this%relperm_dn_itype    (:) = 0
    allocate (this%relperm_dn_param_1  (nconn)); this%relperm_dn_param_1  (:) = 0.d0
    allocate (this%relperm_dn_param_2  (nconn)); this%relperm_dn_param_2  (:) = 0.d0
    allocate (this%relperm_dn_param_3  (nconn)); this%relperm_dn_param_3  (:) = 0.d0

  end subroutine ConnInit

  !------------------------------------------------------------------------
  subroutine ConnDestroy(this)
    !
    implicit none
    !
    class (spac_component_conn_type) :: this

    deallocate (this%id_up               )
    deallocate (this%id_dn               )
    deallocate (this%dist_up             )
    deallocate (this%dist_dn             )
    deallocate (this%area                )
    deallocate (this%itype               )

    deallocate (this%unit_vec            )

    deallocate (this%flux_type           )
    deallocate (this%cond_type           )
    deallocate (this%cond                )
    deallocate (this%cond_up             )
    deallocate (this%cond_dn             )

    deallocate (this%satparam_up_itype   )
    deallocate (this%satparam_up_param_1 )
    deallocate (this%satparam_up_param_2 )
    deallocate (this%satparam_up_param_3 )
    deallocate (this%relperm_up_itype    )
    deallocate (this%relperm_up_param_1  )
    deallocate (this%relperm_up_param_2  )
    deallocate (this%relperm_up_param_3  )

    deallocate (this%satparam_dn_itype   )
    deallocate (this%satparam_dn_param_1 )
    deallocate (this%satparam_dn_param_2 )
    deallocate (this%satparam_dn_param_3 )
    deallocate (this%relperm_dn_itype    )
    deallocate (this%relperm_dn_param_1 )
    deallocate (this%relperm_dn_param_2 )
    deallocate (this%relperm_dn_param_3 )

  end subroutine ConnDestroy

  !------------------------------------------------------------------------
  subroutine ConnCopy(this, idx_beg, idx_end, conn, reverse_updn)
    !
    use petscsys
    !    
    implicit none
    !
    class (spac_component_conn_type) :: this
    class (spac_component_conn_type) :: conn
    PetscInt                         :: idx_beg
    PetscInt                         :: idx_end
    PetscBool , optional             :: reverse_updn 
    !
    PetscBool                        :: reverse

    idx_end = idx_beg + conn%nconn - 1
    reverse = PETSC_FALSE

    if (present(reverse_updn)) reverse = reverse_updn
    
    if (.not. reverse) then
       this%id_up               (idx_beg:idx_end) = conn%id_up               (:)
       this%id_dn               (idx_beg:idx_end) = conn%id_dn               (:)
       this%dist_up             (idx_beg:idx_end) = conn%dist_up             (:)
       this%dist_dn             (idx_beg:idx_end) = conn%dist_dn             (:)
       this%area                (idx_beg:idx_end) = conn%area                (:)
       this%itype               (idx_beg:idx_end) = conn%itype               (:)

       this%unit_vec            (idx_beg:idx_end,:)= conn%unit_vec           (:,:)

       this%flux_type           (idx_beg:idx_end) = conn%flux_type           (:)
       this%cond_type           (idx_beg:idx_end) = conn%cond_type           (:)
       this%cond                (idx_beg:idx_end) = conn%cond                (:)
       this%cond_up             (idx_beg:idx_end) = conn%cond_up             (:)
       this%cond_dn             (idx_beg:idx_end) = conn%cond_dn             (:)

       this%satparam_up_itype   (idx_beg:idx_end) = conn%satparam_up_itype   (:)
       this%satparam_up_param_1 (idx_beg:idx_end) = conn%satparam_up_param_1 (:)
       this%satparam_up_param_2 (idx_beg:idx_end) = conn%satparam_up_param_2 (:)
       this%satparam_up_param_3 (idx_beg:idx_end) = conn%satparam_up_param_3 (:)
       this%relperm_up_itype    (idx_beg:idx_end) = conn%relperm_up_itype    (:)
       this%relperm_up_param_1  (idx_beg:idx_end) = conn%relperm_up_param_1  (:)
       this%relperm_up_param_2  (idx_beg:idx_end) = conn%relperm_up_param_2  (:)

       this%satparam_dn_itype   (idx_beg:idx_end) = conn%satparam_dn_itype   (:)
       this%satparam_dn_param_1 (idx_beg:idx_end) = conn%satparam_dn_param_1 (:)
       this%satparam_dn_param_2 (idx_beg:idx_end) = conn%satparam_dn_param_2 (:)
       this%satparam_dn_param_3 (idx_beg:idx_end) = conn%satparam_dn_param_3 (:)
       this%relperm_dn_itype    (idx_beg:idx_end) = conn%relperm_dn_itype    (:)
       this%relperm_dn_param_1  (idx_beg:idx_end) = conn%relperm_dn_param_1  (:)
       this%relperm_dn_param_2  (idx_beg:idx_end) = conn%relperm_dn_param_2  (:)
       this%relperm_dn_param_3  (idx_beg:idx_end) = conn%relperm_dn_param_3  (:)
    else
       this%id_dn               (idx_beg:idx_end) = conn%id_up               (:)
       this%id_up               (idx_beg:idx_end) = conn%id_dn               (:)
       this%dist_dn             (idx_beg:idx_end) = conn%dist_up             (:)
       this%dist_up             (idx_beg:idx_end) = conn%dist_dn             (:)
       this%area                (idx_beg:idx_end) = conn%area                (:)
       this%itype               (idx_beg:idx_end) = conn%itype               (:)

       this%unit_vec            (idx_beg:idx_end,:)= -conn%unit_vec           (:,:)

       this%flux_type           (idx_beg:idx_end) = conn%flux_type           (:)
       this%cond_type           (idx_beg:idx_end) = conn%cond_type           (:)
       this%cond                (idx_beg:idx_end) = conn%cond                (:)
       this%cond_dn             (idx_beg:idx_end) = conn%cond_up             (:)
       this%cond_up             (idx_beg:idx_end) = conn%cond_dn             (:)

       this%satparam_dn_itype   (idx_beg:idx_end) = conn%satparam_up_itype   (:)
       this%satparam_dn_param_1 (idx_beg:idx_end) = conn%satparam_up_param_1 (:)
       this%satparam_dn_param_2 (idx_beg:idx_end) = conn%satparam_up_param_2 (:)
       this%satparam_dn_param_3 (idx_beg:idx_end) = conn%satparam_up_param_3 (:)
       this%relperm_dn_itype    (idx_beg:idx_end) = conn%relperm_up_itype    (:)
       this%relperm_dn_param_1  (idx_beg:idx_end) = conn%relperm_up_param_1  (:)
       this%relperm_dn_param_2  (idx_beg:idx_end) = conn%relperm_up_param_2  (:)

       this%satparam_up_itype   (idx_beg:idx_end) = conn%satparam_dn_itype   (:)
       this%satparam_up_param_1 (idx_beg:idx_end) = conn%satparam_dn_param_1 (:)
       this%satparam_up_param_2 (idx_beg:idx_end) = conn%satparam_dn_param_2 (:)
       this%satparam_up_param_3 (idx_beg:idx_end) = conn%satparam_dn_param_3 (:)
       this%relperm_up_itype    (idx_beg:idx_end) = conn%relperm_dn_itype    (:)
       this%relperm_up_param_1  (idx_beg:idx_end) = conn%relperm_dn_param_1  (:)
       this%relperm_up_param_2  (idx_beg:idx_end) = conn%relperm_dn_param_2  (:)
       this%relperm_up_param_3  (idx_beg:idx_end) = conn%relperm_dn_param_3  (:)
    endif

  end subroutine ConnCopy

  !------------------------------------------------------------------------
  subroutine AddToMPPAsInternalConn(this, vsfm_mpp, imesh)
    !
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbVSFM      , only : mpp_vsfm_type
    !
    implicit none
    !
    class (spac_component_conn_type) :: this
    class (mpp_vsfm_type)            :: vsfm_mpp
    PetscInt                         :: imesh

    call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
         this%nconn, &
         this%id_up, &
         this%id_dn,       &
         this%dist_up, &
         this%dist_dn,  &
         this%area,  &
         this%itype)

  end subroutine AddToMPPAsInternalConn

  !------------------------------------------------------------------------
  subroutine AddToMPPCouplingBC(this, vsfm_mpp, reverse_up2dn)
    !
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbVSFM      , only : mpp_vsfm_type
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetDestroy
    use MeshType                  , only : MeshCreateConnectionSet
    use petscsys
    !
    implicit none
    !
    class (spac_component_conn_type) :: this
    class (mpp_vsfm_type)            :: vsfm_mpp
    PetscBool                        :: reverse_up2dn
    !
    type(connection_set_type) , pointer :: conn_set
    PetscInt                  , pointer :: id(:)
    PetscInt                            :: iconn
    PetscInt                            :: eqn_id
    PetscInt                  , pointer :: ieqn_other(:)
    PetscReal                 , pointer :: unit_vec(:,:)
    character(len=256)                  :: name


    allocate(ieqn_other(1));

    allocate(id(this%nconn)); id(:) = 0

    if (.not.reverse_up2dn) then
       ieqn_other(1) = this%up_eqn_id
       eqn_id        = this%dn_eqn_id
       name          = trim(this%up_eqn_name) // ' BC in ' // trim(this%dn_eqn_name) // ' equation'

       call MeshCreateConnectionSet(  &
            vsfm_mpp%meshes(eqn_id) , &
            this%nconn              , &
            this%id_up              , &
            this%id_dn              , &
            this%dist_up            , &
            this%dist_dn            , &
            this%area               , &
            this%itype              , &
            this%unit_vec           , &
            skip_id_check=PETSC_TRUE, &
            conn_set=conn_set )
    else
       ieqn_other(1) = this%dn_eqn_id
       eqn_id        = this%up_eqn_id
       name          = trim(this%dn_eqn_name) // ' BC in ' // trim(this%up_eqn_name) // ' equation'

       allocate(unit_vec(this%nconn, 3))
       do iconn = 1, this%nconn
          unit_vec(iconn,1) = -this%unit_vec(iconn,1)
          unit_vec(iconn,2) = -this%unit_vec(iconn,2)
          unit_vec(iconn,3) = -this%unit_vec(iconn,3)
       end do

       call MeshCreateConnectionSet(  &
            vsfm_mpp%meshes(eqn_id) , &
            this%nconn              , &
            this%id_dn              , &
            this%id_up              , &
            this%dist_dn            , &
            this%dist_up            , &
            this%area               , &
            this%itype              , &
            unit_vec                , &
            skip_id_check=PETSC_TRUE, &
            conn_set=conn_set )
       deallocate(unit_vec)
    endif

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn( &
         eqn_id                             , &
         name               = name          , &
         unit               = 'Pa'          , &
         region_type        = 0             , &
         num_other_goveqs   = 1             , &
         id_of_other_goveqs = ieqn_other    , &
         conn_set           = conn_set)

    call ConnectionSetDestroy(conn_set)

    deallocate(id)
    deallocate(ieqn_other)

  end subroutine AddToMPPCouplingBC

end module spac_component

!------------------------------------------------------------------------
module soil_parameters

#include <petsc/finclude/petsc.h>
  
  use spac_component

  implicit none

  PetscInt  , parameter :: soil_nx                          = 2          !
  PetscInt  , parameter :: soil_ny                          = 1          !
  PetscInt  , parameter :: soil_nz                          = 20         !
  PetscReal , parameter :: soil_dx                          = 10.d0      ! [m]
  PetscReal , parameter :: soil_dy                          = 10.d0      ! [m]
  PetscReal , parameter :: soil_dz                          = 0.25d0     ! [m]
  PetscReal             :: slope                                         ! [m m^{-1}]

  PetscReal, parameter  :: perm_xy_top                       = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: perm_z_top                        = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: sat_res_top                       = 0.06d0    ! [-]
  PetscReal, parameter  :: alpha_top                         = 0.00005d0 ! [Pa^{-1}]
  PetscReal, parameter  :: vg_m_top                          = 0.33d0    ! [-]
  PetscReal, parameter  :: por_top                           = 0.5d0     ! [-]

  PetscReal, parameter  :: perm_xy_mid                       = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: perm_z_mid                        = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: sat_res_mid                       = 0.06d0    ! [-]
  PetscReal, parameter  :: alpha_mid                         = 0.00005d0 ! [Pa^{-1}]
  PetscReal, parameter  :: vg_m_mid                          = 0.33d0    ! [-]
  PetscReal, parameter  :: por_mid                           = 0.5d0     ! [-]

  PetscReal, parameter  :: perm_xy_bot                       = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: perm_z_bot                        = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: sat_res_bot                       = 0.06d0    ! [-]
  PetscReal, parameter  :: alpha_bot                         = 0.00005d0 ! [Pa^{-1}]
  PetscReal, parameter  :: vg_m_bot                          = 0.33d0    ! [-]
  PetscReal, parameter  :: por_bot                           = 0.5d0     ! [-]
  
  PetscInt              :: soil_ncells                                   ! Total number of active grid cells
  PetscInt              :: soil_mesh_nconn                               ! Total number of all (vertical + horizontal) connections

  PetscBool             :: is_soil_horizontally_disconnected             ! Are soil grid cells horizontally disconnected?

  type (spac_component_mesh_type) :: soil_mesh
  type (spac_component_pp_type)   :: soil_pp
  type (spac_component_conn_type) :: s2s_conn

  PetscInt  , pointer   :: soil_id(:,:,:)                                ! [-] Id > 0 indicates active grid cell, otherwise inactive
  PetscReal , pointer   :: soil_xc3d(:,:,:)                              ! [m]
  PetscReal , pointer   :: soil_yc3d(:,:,:)                              ! [m]
  PetscReal , pointer   :: soil_zc3d(:,:,:)                              ! [m]

  PetscInt  , pointer   :: top_active_layer_kk_index(:,:)                ! [-] Index of the first active grid cell in z-dimension
  PetscReal , pointer   :: elevation(:,:)                                ! [m]

end module soil_parameters

!------------------------------------------------------------------------
module overstory_parameters
  
#include <petsc/finclude/petsc.h>

  use spac_component
  
  implicit none

  PetscReal , parameter :: overstory_root_depth             = 2.0d0     ! [m]
  PetscReal , parameter :: overstory_canopy_height          = 17.d0     ! [m]
  PetscReal , parameter :: overstory_area_canopy            = 12.25d0   ! [m^2]
  PetscReal , parameter :: overstory_area_sapwood           = 0.013d0   ! [m^2]
  PetscReal , parameter :: overstory_LAI                    = 5.52d0    ! [m^2 m^{-2}]
  PetscReal , parameter :: overstory_area_tappering_coeff   = 0.75d0    ! [-]
  PetscReal , parameter :: overstory_branch_length          = 1.75d0    ! [m]
  PetscReal , parameter :: overstory_branch_area_ratio      = 0.15d0    ! [m^2 branch m^{-2} sapwood]
  PetscReal , parameter :: overstory_root_radius            = 2.9d-4    ! [m]
  PetscReal , parameter :: overstory_root_conductance       = 3.d-11    ! [s^{-1}]
  PetscReal , parameter :: overstory_xylem_Kmax             = 2.5d-5    ! [m s^{-1}]
  PetscReal , parameter :: overstory_xylem_vulnerability_c  = 3.5d0     ! [-]
  PetscReal , parameter :: overstory_xylem_vulnerability_d  = 480.d0    ! [m]
  PetscReal , parameter :: overstory_xylem_phi0             = -2.87d6   ! [Pa]
  PetscReal , parameter :: overstory_xylem_p                = 100.0d0   ! [-]
  PetscReal , parameter :: overstory_xylem_porosity         = 0.57d0    ! [-]
  PetscReal , parameter :: overstory_lad_profile(68)        = (/ &
       0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,  &
       0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,  &
       0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,  &
       0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,  &
       0.00d0, 0.00d0, 0.01d0, 0.03d0, 0.04d0,  &
       0.05d0, 0.06d0, 0.07d0, 0.08d0, 0.09d0,  &
       0.13d0, 0.21d0, 0.30d0, 0.38d0, 0.46d0,  &
       0.52d0, 0.59d0, 0.65d0, 0.71d0, 0.74d0,  &
       0.78d0, 0.81d0, 0.84d0, 0.85d0, 0.85d0,  &
       0.84d0, 0.84d0, 0.83d0, 0.81d0, 0.79d0,  &
       0.77d0, 0.74d0, 0.72d0, 0.69d0, 0.67d0,  &
       0.64d0, 0.61d0, 0.57d0, 0.54d0, 0.51d0,  &
       0.46d0, 0.42d0, 0.37d0, 0.32d0, 0.27d0,  &
       0.20d0, 0.13d0, 0.05d0  &
       /)
  PetscReal , parameter :: overstory_B_profile(8) = (/ &
       45.73d0, 42.82d0, 43.02d0, 39.23d0, 36.80d0,  &
       36.52d0, 21.94d0, 22.83d0  &
       /)

  PetscInt                        :: overstory_xylem_nz                 !
  PetscInt                        :: overstory_leaf_nz                  !
  PetscInt                        :: overstory_root_nz                  !
  PetscInt                        :: overstory_nconn

  PetscReal , pointer             :: overstory_xylem_area_profile(:)    !
  PetscReal , pointer             :: overstory_branch_length_profile(:) !
  PetscBool , pointer             :: overstory_xylem_has_branches(:)    !
  PetscReal , pointer             :: overstory_root_length_profile(:)
  PetscReal , pointer             :: overstory_root_area_profile(:)
  PetscReal , pointer             :: overstory_root_vol_profile(:)
  PetscInt  , pointer             :: overstory_branch_2_xylem_index(:)

  type (spac_component_mesh_type) :: o_root_mesh
  type (spac_component_mesh_type) :: o_xylem_mesh
  type (spac_component_mesh_type) :: o_leaf_mesh

  type (spac_component_pp_type)   :: o_root_pp
  type (spac_component_pp_type)   :: o_xylem_pp
  type (spac_component_pp_type)   :: o_leaf_pp

  type (spac_component_conn_type) :: o_r2s_conn
  type (spac_component_conn_type) :: o_x2r_conn
  type (spac_component_conn_type) :: o_x2x_conn
  type (spac_component_conn_type) :: o_x2l_conn

end module overstory_parameters

!------------------------------------------------------------------------
module understory_parameters
  !
#include <petsc/finclude/petsc.h>
  !
  use spac_component
  !
  implicit none
  !
  PetscReal , parameter :: understory_root_depth            = 0.5d0     ! [m]
  PetscReal , parameter :: understory_canopy_height         = 1.d0      ! [m]
  PetscReal , parameter :: understory_area_canopy           = 9.0d0     ! [m^2]
  PetscReal , parameter :: understory_area_sapwood          = 0.010d0   ! [m^2]
  PetscReal , parameter :: understory_LAI                   = 3.19d0    ! [m^2 m^{-2}]
  PetscReal , parameter :: understory_area_tappering_coeff  = 0.75d0    ! [-]
  PetscReal , parameter :: understory_branch_length         = 1.50d0    ! [m]
  PetscReal , parameter :: understory_branch_area_ratio     = 0.15d0    ! [m^2 branch m^{-2} sapwood]
  PetscReal , parameter :: understory_root_radius           = 2.9d-4    ! [m]
  PetscReal , parameter :: understory_root_conductance      = 3.d-11    ! [s^{-1}]
  PetscReal , parameter :: understory_xylem_Kmax            = 2.5d-5    ! [m s^{-1}]
  PetscReal , parameter :: understory_xylem_vulnerability_c = 3.5d0     ! [-]
  PetscReal , parameter :: understory_xylem_vulnerability_d = 480.d0    ! [m]
  PetscReal , parameter :: understory_xylem_phi0            = -2.87d6   ! [Pa]
  PetscReal , parameter :: understory_xylem_p               = 100.0d0   ! [-]
  PetscReal , parameter :: understory_xylem_porosity        = 0.57d0    ! [-]
  PetscReal , parameter :: understory_lad_profile(24) = (/ &
       0.00d0, 0.07d0, 0.21d0, 0.35d0, 0.49d0,  &
       0.54d0, 0.57d0, 0.61d0, 0.64d0, 0.66d0,  &
       0.67d0, 0.69d0, 0.70d0, 0.70d0, 0.69d0,  &
       0.68d0, 0.66d0, 0.65d0, 0.61d0, 0.58d0,  &
       0.54d0, 0.50d0, 0.39d0, 0.28d0  &
       /)
  PetscReal , parameter :: understory_B_profile(2) = (/ &
       0.76d0, 0.16d0  &
       /)

  PetscInt                        :: understory_xylem_nz                 !
  PetscInt                        :: understory_leaf_nz                  !
  PetscInt                        :: understory_root_nz                  !
  PetscInt                        :: understory_nconn
  
  PetscReal , pointer             :: understory_xylem_area_profile(:)    !
  PetscReal , pointer             :: understory_branch_length_profile(:) !
  PetscBool , pointer             :: understory_xylem_has_branches(:)    !
  PetscReal , pointer             :: understory_root_length_profile(:)   !
  PetscReal , pointer             :: understory_root_area_profile(:)     !
  PetscReal , pointer             :: understory_root_vol_profile(:)      !
  PetscInt  , pointer             :: understory_branch_2_xylem_index(:)

  type (spac_component_mesh_type) :: u_root_mesh
  type (spac_component_mesh_type) :: u_xylem_mesh
  type (spac_component_mesh_type) :: u_leaf_mesh

  type (spac_component_pp_type)   :: u_root_pp
  type (spac_component_pp_type)   :: u_xylem_pp
  type (spac_component_pp_type)   :: u_leaf_pp

  type (spac_component_conn_type) :: u_r2s_conn
  type (spac_component_conn_type) :: u_x2r_conn
  type (spac_component_conn_type) :: u_x2x_conn
  type (spac_component_conn_type) :: u_x2l_conn

end module understory_parameters

!------------------------------------------------------------------------
module problem_parameters

  use spac_component
  use soil_parameters
  use overstory_parameters
  use understory_parameters

  implicit none

  PetscReal , parameter :: PI       = 4 * atan (1.0_8) ! [-]
  PetscReal , parameter :: vish2o   = 0.001002d0       ! [N s/m^2] @ 20 degC

  PetscReal , parameter :: init_wtd = 3.d0             ! Initial water table depth below surface [m]
  
  PetscInt              :: ncells

  PetscBool             :: multi_goveqns_formulation

  type(spac_component_mesh_type) :: combined_mesh
  type(spac_component_pp_type)   :: combined_pp
  type(spac_component_conn_type) :: combined_conn

  PetscInt, parameter :: soil_eqn_id    = 1
  PetscInt, parameter :: o_root_eqn_id  = 2
  PetscInt, parameter :: o_xylem_eqn_id = 3
  PetscInt, parameter :: o_leaf_eqn_id  = 4
  PetscInt, parameter :: u_root_eqn_id  = 5
  PetscInt, parameter :: u_xylem_eqn_id = 6
  PetscInt, parameter :: u_leaf_eqn_id  = 7

  PetscInt, parameter :: soil_mesh_id    = 1
  PetscInt, parameter :: o_root_mesh_id  = 2
  PetscInt, parameter :: o_xylem_mesh_id = 3
  PetscInt, parameter :: o_leaf_mesh_id  = 4
  PetscInt, parameter :: u_root_mesh_id  = 5
  PetscInt, parameter :: u_xylem_mesh_id = 6
  PetscInt, parameter :: u_leaf_mesh_id  = 7

end module problem_parameters

!------------------------------------------------------------------------
program vsfm_spac_on_hillslope

#include <petsc/finclude/petsc.h>

  use petscsys
  !
  implicit none
  !
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call run_vsfm_spac_on_hillslope()

  call PetscFinalize(ierr)

end program vsfm_spac_on_hillslope

!------------------------------------------------------------------------
subroutine run_vsfm_spac_on_hillslope()
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
  use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
  use problem_parameters
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda
  !
  implicit none
  !
  PetscBool      :: converged
  PetscInt       :: converged_reason
  PetscErrorCode :: ierr
  PetscReal      :: time
  PetscReal      :: dtime
  PetscInt       :: nstep
  PetscInt       :: istep
  PetscViewer    :: viewer
  Vec            :: sat
  PetscReal, pointer :: sat_p(:)
  PetscReal, pointer :: data(:)
  PetscBool          :: flg

  ! Initialize data
  time  = 0.d0

  ! Set default settings
  dtime                             = 180.d0
  nstep                             = 1
  slope                             = 0.05d0
  slope                             = 0.0d0
  is_soil_horizontally_disconnected = PETSC_FALSE
  multi_goveqns_formulation         = PETSC_FALSE
  
  call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-slope',slope,flg,ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-soil_horizontally_disconnected',is_soil_horizontally_disconnected,flg,ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-multi_goveqns_formulation',multi_goveqns_formulation,flg,ierr)

  call initialize_problem()

  call vsfm_mpp%soe%SetDtime(1.d0)
  call vsfm_mpp%soe%PreSolve()
  call vsfm_mpp%soe%PostSolve()

  allocate(data(soil_ncells))
  call VecDuplicate(vsfm_mpp%soe%solver%soln, sat, ierr); CHKERRQ(ierr)

  call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL, VAR_LIQ_SAT, -1, data)

  call VecGetArrayF90(sat, sat_p, ierr); CHKERRQ(ierr)
  sat_p(:) = data(:)
  call VecRestoreArrayF90(sat, sat_p, ierr); CHKERRQ(ierr)

  call PetscViewerBinaryOpen(PETSC_COMM_SELF, 'sat_init.bin', FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
  call VecView(sat, viewer, ierr); CHKERRQ(ierr)
  call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  do istep = 1, nstep

     !call set_bondary_conditions(istep)
     time = time + dtime

     ! Run the model
     call vsfm_mpp%soe%StepDT(dtime, istep, &
          converged, converged_reason, ierr); CHKERRQ(ierr)

  end do

  call PetscViewerBinaryOpen(PETSC_COMM_SELF, 'pressure_final.bin', FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
  call VecView(vsfm_mpp%soe%solver%soln, viewer, ierr); CHKERRQ(ierr)
  call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL, VAR_LIQ_SAT, -1, data)

  call VecGetArrayF90(sat, sat_p, ierr); CHKERRQ(ierr)
  sat_p(:) = data(:)
  call VecRestoreArrayF90(sat, sat_p, ierr); CHKERRQ(ierr)

  call PetscViewerBinaryOpen(PETSC_COMM_SELF, 'sat_final.bin', FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
  call VecView(sat, viewer, ierr); CHKERRQ(ierr)
  call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  deallocate(data)
  call VecDestroy(sat, ierr)

end subroutine run_vsfm_spac_on_hillslope

!------------------------------------------------------------------------
subroutine initialize_problem()
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use problem_parameters
  !
  implicit none

  ! 1. Initialize the multi-physics-problem (MPP)
  call initialize_mpp()
  
  ! 2. Setup meshes
  call setup_soil_mesh()
  call setup_overstory_mesh()
  call setup_understory_mesh()

  if (.not. multi_goveqns_formulation) then
     ! 2.1 Add the combined mesh
     call add_single_mesh()

     ! 3. Add all governing equations
     call add_single_goveqn()

  else
     ! 2.1 Add all meshes
     call add_multiple_meshes()

     ! 3. Add all governing equations
     call add_multiple_goveqns()
  end if


  ! 4. Add boundary and source-sink conditions to all governing equations
  call add_conditions_to_goveqns()

  ! 5. Allocate memory to hold auxvars
  call allocate_auxvars()

  ! 6. Setup the MPP
  call vsfm_mpp%SetupProblem()

  ! 7. Add material properities associated with all governing equations
  if (.not. multi_goveqns_formulation) then
     call set_material_properties()
  else
     call set_material_properties_multi_goveqns()
  endif

  ! 8. Set initial conditions
  call set_initial_conditions()

end subroutine initialize_problem

!------------------------------------------------------------------------
subroutine initialize_mpp()
  !
  ! !DESCRIPTION:
  ! Initialization VSFM
  !
  !
#include <petsc/finclude/petsc.h>
  ! !USES:
  use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use petscsys
  !
  ! !ARGUMENTS
  implicit none
  !
  PetscInt       :: iam
  PetscErrorCode :: ierr

  call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

  !
  ! Set up the multi-physics problem
  !
  call vsfm_mpp%Init       ()
  call vsfm_mpp%SetName    ('Tree hydrosoil_dynamics on a hillslope')
  call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
  call vsfm_mpp%SetMPIRank (iam)

end subroutine initialize_mpp

!------------------------------------------------------------------------
subroutine add_single_mesh()
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
  use problem_parameters
  !
  implicit none
  !
#include <petsc/finclude/petsc.h>
  !
  PetscInt            :: imesh
  PetscInt            :: soil_ncells_ghost
  PetscInt            :: nconn
  PetscInt            :: idx_beg, idx_end  
  PetscErrorCode      :: ierr

  !
  ! Meshes
  !
  ncells = soil_mesh%ncells             + &
           o_root_mesh%ncells   + &
           o_xylem_mesh%ncells  + &
           o_leaf_mesh%ncells   + &
           u_root_mesh%ncells  + &
           u_xylem_mesh%ncells + &
           u_leaf_mesh%ncells 

  call combined_mesh%Init(ncells)
  call combined_pp%Init(  ncells)

  ! Add soil grid cells
  idx_beg = 1;  idx_end = idx_beg + soil_mesh%ncells - 1
  call combined_mesh%Copy(idx_beg, idx_end, soil_mesh)
  call combined_pp%Copy(  idx_beg, idx_end, soil_pp  )

  ! Overstory: Add root grid cells
  idx_beg = idx_end + 1; idx_end = idx_beg + o_root_mesh%ncells - 1
  call combined_mesh%Copy(idx_beg, idx_end, o_root_mesh)
  call combined_pp%Copy(  idx_beg, idx_end, o_root_pp  )

  ! Overstory: Add xylem grid cells
  idx_beg = idx_end + 1; idx_end = idx_beg + o_xylem_mesh%ncells - 1
  call combined_mesh%Copy(idx_beg, idx_end, o_xylem_mesh)
  call combined_pp%Copy(  idx_beg, idx_end, o_xylem_pp  )

  ! Overstory: Add leaf grid cells
  idx_beg = idx_end + 1; idx_end = idx_beg + o_leaf_mesh%ncells - 1
  call combined_mesh%Copy(idx_beg, idx_end, o_leaf_mesh)
  call combined_pp%Copy(  idx_beg, idx_end, o_leaf_pp  )

  ! Understory: Add root grid cells
  idx_beg = idx_end + 1; idx_end = idx_beg + u_root_mesh%ncells - 1
  call combined_mesh%Copy(idx_beg, idx_end, u_root_mesh)
  call combined_pp%Copy(  idx_beg, idx_end, u_root_pp  )

  ! Understory: Add xylem grid cells
  idx_beg = idx_end + 1; idx_end = idx_beg + u_xylem_mesh%ncells - 1
  call combined_mesh%Copy(idx_beg, idx_end, u_xylem_mesh)
  call combined_pp%Copy(  idx_beg, idx_end, u_xylem_pp  )

  ! Understory: Add leaf grid cells
  idx_beg = idx_end + 1; idx_end = idx_beg + u_leaf_mesh%ncells - 1
  call combined_mesh%Copy(idx_beg, idx_end, u_leaf_mesh)
  call combined_pp%Copy(  idx_beg, idx_end, u_leaf_pp  )

  !
  ! Set up the meshes
  !    
  call vsfm_mpp%SetNumMeshes(1)

  imesh             = 1
  soil_ncells_ghost = 0

  call combined_mesh%AddToMPP(vsfm_mpp, imesh, &
       'Combined mesh for Soil-Overstory-Understory', soil_nz)

  !
  ! Add connections
  !
  nconn = s2s_conn%nconn   + &
          o_r2s_conn%nconn + &
          o_x2r_conn%nconn + &
          o_x2x_conn%nconn + &
          o_x2l_conn%nconn + &
          u_r2s_conn%nconn + &
          u_x2r_conn%nconn + &
          u_x2x_conn%nconn + &
          u_x2l_conn%nconn

  call combined_conn%Init(nconn)

  ! Add soil-2-soil connections
  idx_beg = 1; call combined_conn%Copy(idx_beg, idx_end, s2s_conn)

  ! Add overstory root-2-soil connections
  idx_beg = idx_end + 1; call combined_conn%Copy(idx_beg, idx_end, o_r2s_conn)

  ! Add overstory xylem-2-root connections
  idx_beg = idx_end + 1; call combined_conn%Copy(idx_beg, idx_end, o_x2r_conn)

  ! Add overstory xylem-2-xylem connections
  idx_beg = idx_end + 1; call combined_conn%Copy(idx_beg, idx_end, o_x2x_conn)

  ! Add overstory xylem-2-leaf connections
  idx_beg = idx_end + 1; call combined_conn%Copy(idx_beg, idx_end, o_x2l_conn)

  ! Add understory root-2-soil connections
  idx_beg = idx_end + 1; call combined_conn%Copy(idx_beg, idx_end, u_r2s_conn)

  ! Add understory xylem-2-root connections
  idx_beg = idx_end + 1; call combined_conn%Copy(idx_beg, idx_end, u_x2r_conn)

  ! Add understory xylem-2-xylem connections
  idx_beg = idx_end + 1; call combined_conn%Copy(idx_beg, idx_end, u_x2x_conn)

  ! Add understory xylem-2-leaf connections
  idx_beg = idx_end + 1; call combined_conn%Copy(idx_beg, idx_end, u_x2l_conn)

  imesh = 1
  call combined_conn%AddToMPPAsInternalConn(vsfm_mpp, imesh)

end subroutine add_single_mesh

!------------------------------------------------------------------------
subroutine add_multiple_meshes()
  !
#include <petsc/finclude/petsc.h>
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
  use problem_parameters
  !
  implicit none
  !
  PetscInt :: num_cells
  PetscInt :: ncells_ghost
  PetscInt :: nz
  PetscInt :: imesh

  ! Create memory for 7 meshes
  call vsfm_mpp%SetNumMeshes(7)

  ! 1. Add mesh for soil
  imesh = 1
  call soil_mesh%AddToMPP(vsfm_mpp, imesh, 'Soil mesh', soil_nz)

  ! 2. Add mesh for overstory root
  imesh     = imesh + 1
  nz        = overstory_root_nz
  call o_root_mesh%AddToMPP(vsfm_mpp, imesh, 'Overstory root mesh', nz)

  ! 3. Add mesh for overstory xylem
  imesh     = imesh + 1
  nz        = overstory_xylem_nz
  call o_xylem_mesh%AddToMPP(vsfm_mpp, imesh, 'Overstory xylem mesh', nz)

  ! 4. Add mesh for overstory leaf
  imesh     = imesh + 1
  nz        = overstory_leaf_nz
  call o_leaf_mesh%AddToMPP(vsfm_mpp, imesh, 'Overstory leaf mesh', nz)

  ! 5. Add mesh for understory root
  imesh     = imesh + 1
  nz        = understory_root_nz
  call u_root_mesh%AddToMPP(vsfm_mpp, imesh, 'Understory root mesh', nz)

  ! 6. Add mesh for understory xylem
  imesh     = imesh + 1
  nz        = understory_xylem_nz
  call u_xylem_mesh%AddToMPP(vsfm_mpp, imesh, 'Understory xylem mesh', nz)

  ! 7. Add mesh for understory leaf
  imesh     = imesh + 1
  nz        = understory_leaf_nz
  call u_leaf_mesh%AddToMPP(vsfm_mpp, imesh, 'Understory leaf mesh', nz)

  !
  ! Add connections
  !
  imesh = 1
  call s2s_conn%AddToMPPAsInternalConn(vsfm_mpp, imesh)

  imesh = 3
  call o_x2x_conn%AddToMPPAsInternalConn(vsfm_mpp, imesh)
  
  imesh = 6
  call u_x2x_conn%AddToMPPAsInternalConn(vsfm_mpp, imesh)
  
end subroutine add_multiple_meshes

!------------------------------------------------------------------------
subroutine setup_soil_mesh()
  !
  use SaturationFunction , only : RELPERM_FUNC_MUALEM
  use SaturationFunction , only : SAT_FUNC_VAN_GENUCHTEN
  use MultiPhysicsProbConstants , only : CONN_VERTICAL
  use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
  use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
  use problem_parameters
  !
  implicit none
  !
#include <petsc/finclude/petsc.h>
  !
  PetscInt :: count
  PetscInt :: ii
  PetscInt :: jj
  PetscInt :: kk
  PetscInt :: iconn
  
  allocate(soil_id                   (soil_nx,soil_ny,soil_nz ))
  allocate(soil_xc3d                 (soil_nx,soil_ny,soil_nz ))
  allocate(soil_yc3d                 (soil_nx,soil_ny,soil_nz ))
  allocate(soil_zc3d                 (soil_nx,soil_ny,soil_nz ))
  allocate(top_active_layer_kk_index (soil_nx,soil_ny         )); top_active_layer_kk_index(:,:) = 0
  allocate(elevation                 (soil_nx,soil_ny         ))

  soil_id(:,:,:) = 0

  count = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, soil_nz
           soil_xc3d(ii,jj,kk) =  soil_dx/2 + soil_dx*(ii-1)
           soil_yc3d(ii,jj,kk) =  soil_dy/2 + soil_dy*(jj-1)
           soil_zc3d(ii,jj,kk) = -soil_dz/2 - soil_dz*(kk-1)
           if (soil_zc3d(ii,jj,kk) <= soil_zc3d(1,1,1) -slope * soil_dx*(ii-1) ) then
              count = count + 1
              soil_id(ii,jj,kk) = count
           endif
        end do
     end do
  end do

  soil_ncells       = count

  call soil_mesh%Init(soil_ncells)
  call soil_pp%Init(  soil_ncells)

  count           = 0
  soil_mesh_nconn = 0

  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, soil_nz

           ! is the grid cell active?
           if (soil_id(ii,jj,kk) > 0) then

              count = count + 1

              soil_mesh%xc(count)      = soil_xc3d(ii,jj,kk)
              soil_mesh%yc(count)      = soil_yc3d(ii,jj,kk)
              soil_mesh%zc(count)      = soil_zc3d(ii,jj,kk)
              soil_mesh%dx(count)      = soil_dx
              soil_mesh%dy(count)      = soil_dy
              soil_mesh%dz(count)      = soil_dz
              soil_mesh%area(count)    = soil_dx*soil_dy
              soil_mesh%vol(count)     = soil_dx*soil_dy*soil_dz
              soil_mesh%filter(count)  = 1

              ! add vertical connection except if the grid cell is last in the z-dir
              if (kk < soil_nz) soil_mesh_nconn = soil_mesh_nconn + 1

              ! should horizontal connections be added?
              if (.not.is_soil_horizontally_disconnected) then
                 ! add horizontal connection except if the grid cell is last in x-dir
                 if (ii < soil_nx .and. soil_id(ii+1,jj,kk) > 0) soil_mesh_nconn = soil_mesh_nconn + 1
              end if

              ! is the index of first active grid cell in z-dir already set?
              if (top_active_layer_kk_index(ii,jj) == 0) then
                 top_active_layer_kk_index(ii,jj) = kk
                 elevation                (ii,jj) = soil_zc3d(ii,jj,kk) + soil_dz/2.d0
              end if

           end if
        end do
     end do
  end do

  soil_pp%por             (1:soil_ncells) = por_top
  soil_pp%perm            (1:soil_ncells) = perm_z_top
  soil_pp%alpha           (1:soil_ncells) = alpha_top
  soil_pp%lambda          (1:soil_ncells) = vg_m_top
  soil_pp%relperm_type    (1:soil_ncells) = RELPERM_FUNC_MUALEM
  soil_pp%relperm_param_1 (1:soil_ncells) = 0.d0
  soil_pp%relperm_param_2 (1:soil_ncells) = 0.d0
  soil_pp%residual_sat    (1:soil_ncells) = sat_res_top
  soil_pp%satfunc_type    (1:soil_ncells) = SAT_FUNC_VAN_GENUCHTEN

  call s2s_conn%Init(soil_mesh_nconn)

  !  
  ! Soil-to-Soil Connections: Add vertical connections
  !
  iconn = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, soil_nz-1
           if (soil_id(ii,jj,kk) > 0 .and. soil_id(ii,jj,kk+1) > 0) then
              iconn                                  = iconn + 1
              s2s_conn%id_up(iconn)     = soil_id(ii,jj,kk  )
              s2s_conn%id_dn(iconn)     = soil_id(ii,jj,kk+1)
              s2s_conn%dist_up(iconn)   = 0.5d0*soil_dz
              s2s_conn%dist_dn(iconn)   = 0.5d0*soil_dz
              s2s_conn%area(iconn)      = soil_dx*soil_dy
              s2s_conn%itype(iconn)      = CONN_VERTICAL
              s2s_conn%flux_type(iconn) = DARCY_FLUX_TYPE
           end if
        end do
     end do
  end do
  
  !  
  ! Soil-to-Soil Connections: Add horizontal connections in x-direction
  !
  if (.not.is_soil_horizontally_disconnected) then
     do ii = 1, soil_nx-1
        do jj = 1, soil_ny
           do kk = 1, soil_nz
              if (soil_id(ii,jj,kk) > 0 .and. soil_id(ii+1,jj,kk) > 0) then
                 iconn                                  = iconn + 1
                 s2s_conn%id_up(iconn)     = soil_id(ii  ,jj,kk)
                 s2s_conn%id_dn(iconn)     = soil_id(ii+1,jj,kk)
                 s2s_conn%dist_up(iconn)   = 0.5d0*soil_dx
                 s2s_conn%dist_dn(iconn)   = 0.5d0*soil_dx
                 s2s_conn%area(iconn)      = soil_dy*soil_dz
                 s2s_conn%itype(iconn)      = CONN_HORIZONTAL
                 s2s_conn%flux_type(iconn) = DARCY_FLUX_TYPE
              end if
           end do
        end do
     end do

     if (soil_ny > 1) then
        write(*,*)'Add code to for horizontal connections in y-direciton'
        stop
     endif
  endif

end subroutine setup_soil_mesh

!------------------------------------------------------------------------
subroutine setup_overstory_mesh()
  !
#include <petsc/finclude/petsc.h>
  !
  use mpp_varcon         , only : grav
  use mpp_varcon         , only : denh2o
  use SaturationFunction , only : RELPERM_FUNC_WEIBULL
  use SaturationFunction , only : SAT_FUNC_CHUANG
  use MultiPhysicsProbConstants , only : CONN_VERTICAL
  use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
  use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
  use MultiPhysicsProbConstants , only : CONDUCTANCE_CAMPBELL_TYPE
  use MultiPhysicsProbConstants , only : CONDUCTANCE_MANOLI_TYPE
  use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
  use petscsys
  use problem_parameters
  !
  implicit none
  !
  PetscInt  :: ii,jj,kk
  PetscInt  :: idx
  PetscInt  :: id_value
  PetscInt  :: count
  PetscInt  :: iconn
  PetscReal :: zz
  PetscReal :: soil_volume
  PetscReal :: root_length
  PetscInt  :: root_id_offset
  PetscInt  :: xylem_id_offset
  PetscInt  :: leaf_id_offset
  PetscReal :: dist_x, dist_y, dist_z, dist

  overstory_leaf_nz  = 0
  overstory_xylem_nz = int(overstory_canopy_height/soil_dz)
  overstory_root_nz  = int(overstory_root_depth   /soil_dz)

  allocate (overstory_xylem_area_profile    (overstory_xylem_nz))
  allocate (overstory_xylem_has_branches    (overstory_xylem_nz))
  allocate (overstory_branch_length_profile (overstory_xylem_nz))
  allocate (overstory_branch_2_xylem_index  (overstory_xylem_nz))
  allocate (overstory_root_length_profile   (overstory_root_nz ))
  allocate (overstory_root_area_profile     (overstory_root_nz ))
  allocate (overstory_root_vol_profile      (overstory_root_nz ))
  
  overstory_branch_length_profile(:) = 0.d0
  count = 0
  
  do kk = 1, overstory_xylem_nz
     zz = (kk-1)*soil_dz + soil_dz/2.d0
     overstory_xylem_area_profile(kk) = overstory_area_sapwood * &
          ( 1.d0 - overstory_area_tappering_coeff * zz/overstory_canopy_height) ** 2.d0

     if (overstory_lad_profile(kk) > 0.d0) then
        count                               = count + 1
        overstory_leaf_nz                   = overstory_leaf_nz + 1;
        overstory_xylem_has_branches(kk)    = PETSC_TRUE
        overstory_branch_length_profile(kk) = overstory_xylem_area_profile(kk) * &
                                              overstory_branch_area_ratio
        overstory_branch_2_xylem_index(count) = kk
     else
        overstory_xylem_has_branches(kk) = PETSC_FALSE
     end if
  end do

  soil_volume = soil_dx * soil_dy * soil_dz
  do kk = 1, overstory_root_nz
     root_length                       = overstory_B_profile(kk) * soil_volume
     overstory_root_length_profile(kk) = root_length
     overstory_root_area_profile(kk)   = 2*PI*overstory_root_radius*root_length
     overstory_root_vol_profile(kk)    = PI*(overstory_root_radius**2.d0)*root_length
  end do


  call o_root_mesh%Init(soil_nx*soil_ny*overstory_root_nz)
  call o_root_pp%Init(  soil_nx*soil_ny*overstory_root_nz)

  call o_xylem_mesh%Init(soil_nx*soil_ny*overstory_xylem_nz)
  call o_xylem_pp%Init(  soil_nx*soil_ny*overstory_xylem_nz)

  call o_leaf_mesh%Init(soil_nx*soil_ny*overstory_leaf_nz)
  call o_leaf_pp%Init(  soil_nx*soil_ny*overstory_leaf_nz)

  allocate(o_root_mesh%id               (soil_nx, soil_ny, overstory_root_nz ))
  allocate(o_xylem_mesh%id              (soil_nx, soil_ny, overstory_xylem_nz))
  allocate(o_leaf_mesh%id               (soil_nx, soil_ny, overstory_leaf_nz ))

  if (multi_goveqns_formulation) then
     id_value = 0
  else
     id_value = soil_ncells
  end if
  root_id_offset = id_value
  
  ! Add root grid cells
  count = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, overstory_root_nz

           id_value                                  = id_value + 1
           o_root_mesh%id(ii,jj,kk)          = id_value

           count                                     = count + 1

           o_root_mesh%xc            (count) = soil_xc3d(ii,jj,1) - 0.1d0
           o_root_mesh%yc            (count) = soil_yc3d(ii,jj,1)
           o_root_mesh%zc            (count) = elevation(ii,jj) - soil_dz/2.d0 - soil_dz*(kk-1)
           o_root_mesh%dx            (count) = soil_dx
           o_root_mesh%dy            (count) = soil_dy
           o_root_mesh%dz            (count) = soil_dz
           o_root_mesh%area          (count) = overstory_root_area_profile(kk)
           o_root_mesh%vol           (count) = overstory_root_vol_profile(kk)
           o_root_mesh%filter        (count) = 1

           o_root_pp%por             (count) = 0.d0
           o_root_pp%perm            (count) = 0.d0
           o_root_pp%alpha           (count) = overstory_xylem_phi0
           o_root_pp%lambda          (count) = overstory_xylem_p
           o_root_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           o_root_pp%relperm_param_1 (count) = overstory_xylem_vulnerability_d * grav * denh2o
           o_root_pp%relperm_param_2 (count) = overstory_xylem_vulnerability_c
           o_root_pp%residual_sat    (count) = 0.d0
           o_root_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  ! Add xylem grid cells
  count = 0
  if (multi_goveqns_formulation) id_value = 0
  xylem_id_offset = id_value

  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, overstory_xylem_nz
           id_value                                   = id_value + 1
           o_xylem_mesh%id        (ii,jj,kk)  = id_value

           count                                      = count + 1
           o_xylem_mesh%xc            (count) = soil_xc3d(ii,jj,1) - 0.1d0
           o_xylem_mesh%yc            (count) = soil_yc3d(ii,jj,1)
           o_xylem_mesh%zc            (count) = elevation(ii,jj) + soil_dz/2.d0 + (kk-1)*soil_dz
           o_xylem_mesh%dx            (count) = soil_dx
           o_xylem_mesh%dy            (count) = soil_dy
           o_xylem_mesh%dz            (count) = soil_dz
           o_xylem_mesh%area          (count) = overstory_xylem_area_profile(kk)
           o_xylem_mesh%vol           (count) = overstory_xylem_area_profile(kk) * soil_dz
           o_xylem_mesh%filter        (count) = 1

           o_xylem_pp%por             (count) = overstory_xylem_porosity
           o_xylem_pp%perm            (count) = overstory_xylem_Kmax * vish2o / (denh2o * grav)
           o_xylem_pp%alpha           (count) = overstory_xylem_phi0
           o_xylem_pp%lambda          (count) = overstory_xylem_p
           o_xylem_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           o_xylem_pp%relperm_param_1 (count) = overstory_xylem_vulnerability_d * grav * denh2o
           o_xylem_pp%relperm_param_2 (count) = overstory_xylem_vulnerability_c
           o_xylem_pp%residual_sat    (count) = 0.d0
           o_xylem_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  ! Add leaf grid cells
  count = 0
  if (multi_goveqns_formulation) id_value = 0
  leaf_id_offset = id_value

  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, overstory_leaf_nz
           id_value                                  = id_value + 1
           o_leaf_mesh%id        (ii,jj,kk)  = id_value

           count                                     = count + 1
           idx                                       = overstory_branch_2_xylem_index(kk)

           o_leaf_mesh%xc            (count) = soil_xc3d(ii,jj,1)-overstory_branch_length_profile(idx) - 0.1d0
           o_leaf_mesh%yc            (count) = soil_yc3d(ii,jj,1)
           o_leaf_mesh%zc            (count) = elevation(ii,jj) + soil_dz/2.d0 + (kk-1)*soil_dz + (overstory_xylem_nz-overstory_leaf_nz)*soil_dz
           o_leaf_mesh%dx            (count) = soil_dx
           o_leaf_mesh%dy            (count) = soil_dy
           o_leaf_mesh%dz            (count) = soil_dz
           o_leaf_mesh%area          (count) = overstory_xylem_area_profile(kk) * overstory_branch_area_ratio

           o_leaf_mesh%vol           (count) = o_leaf_mesh%area(count) * overstory_branch_length_profile(idx)
           o_leaf_mesh%filter        (count) = 1

           o_leaf_pp%por             (count) = 0.d0
           o_leaf_pp%perm            (count) = overstory_xylem_Kmax * vish2o / (denh2o * grav)
           o_leaf_pp%alpha           (count) = overstory_xylem_phi0
           o_leaf_pp%lambda          (count) = overstory_xylem_p
           o_leaf_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           o_leaf_pp%relperm_param_1 (count) = overstory_xylem_vulnerability_d * grav * denh2o
           o_leaf_pp%relperm_param_2 (count) = overstory_xylem_vulnerability_c
           o_leaf_pp%residual_sat    (count) = 0.d0
           o_leaf_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  call o_r2s_conn%Init (overstory_root_nz*        soil_nx * soil_ny)
  call o_x2r_conn%Init (overstory_root_nz*        soil_nx * soil_ny)
  call o_x2x_conn%Init ((overstory_xylem_nz - 1)* soil_nx * soil_ny)
  call o_x2l_conn%Init (overstory_leaf_nz*        soil_nx * soil_ny)

  ! Root-to-Soil: Condunctance flux type
  o_r2s_conn%up_eqn_id   = o_root_eqn_id
  o_r2s_conn%dn_eqn_id   = soil_eqn_id
  o_r2s_conn%up_eqn_name = 'Overstory root'
  o_r2s_conn%dn_eqn_name = 'Soil'

  iconn = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, overstory_root_nz
           iconn                                  = iconn + 1
           o_r2s_conn%id_up               (iconn) = o_root_mesh%id(ii,jj,kk)
           o_r2s_conn%id_dn               (iconn) = soil_id(ii,jj,kk-1+top_active_layer_kk_index(ii,jj))

           o_r2s_conn%dist_up             (iconn) = 0.d0
           o_r2s_conn%dist_dn             (iconn) = overstory_root_length_profile(kk)
           o_r2s_conn%area                (iconn) = overstory_root_area_profile  (kk)
           o_r2s_conn%itype               (iconn) = CONN_HORIZONTAL
           o_r2s_conn%flux_type           (iconn) = CONDUCTANCE_FLUX_TYPE
           o_r2s_conn%cond_type           (iconn) = CONDUCTANCE_MANOLI_TYPE
           o_r2s_conn%cond_up             (iconn) = overstory_root_conductance
           o_r2s_conn%cond_dn             (iconn) = perm_z_top /vish2o * (denh2o * grav) / & ! [m/s]
                                                              overstory_root_length_profile(kk)        ! [m]

           ! Set unit vector
           dist_x = soil_mesh%xc(o_r2s_conn%id_dn(iconn)) - o_root_mesh%xc(o_r2s_conn%id_up(iconn) - root_id_offset)
           dist_y = soil_mesh%yc(o_r2s_conn%id_dn(iconn)) - o_root_mesh%yc(o_r2s_conn%id_up(iconn) - root_id_offset)
           dist_z = soil_mesh%zc(o_r2s_conn%id_dn(iconn)) - o_root_mesh%zc(o_r2s_conn%id_up(iconn) - root_id_offset)

           dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0
           o_r2s_conn%unit_vec(iconn,1) = dist_x/dist
           o_r2s_conn%unit_vec(iconn,2) = dist_y/dist
           o_r2s_conn%unit_vec(iconn,3) = dist_z/dist
           
           ! Saturation up: CHUANG
           o_r2s_conn%satparam_up_itype   (iconn) = o_root_pp%satfunc_type    (o_r2s_conn%id_up(iconn) - root_id_offset )
           o_r2s_conn%satparam_up_param_1 (iconn) = o_root_pp%alpha           (o_r2s_conn%id_up(iconn) - root_id_offset )
           o_r2s_conn%satparam_up_param_2 (iconn) = o_root_pp%lambda          (o_r2s_conn%id_up(iconn) - root_id_offset )
           o_r2s_conn%satparam_up_param_3 (iconn) = o_root_pp%residual_sat    (o_r2s_conn%id_up(iconn) - root_id_offset )
           ! Relative Perm. up: WEIBULL
           o_r2s_conn%relperm_up_itype    (iconn) = o_root_pp%relperm_type    (o_r2s_conn%id_up(iconn) - root_id_offset )
           o_r2s_conn%relperm_up_param_1  (iconn) = o_root_pp%relperm_param_1 (o_r2s_conn%id_up(iconn) - root_id_offset )
           o_r2s_conn%relperm_up_param_2  (iconn) = o_root_pp%relperm_param_2 (o_r2s_conn%id_up(iconn) - root_id_offset )

           ! Saturation dn: CHUANG
           o_r2s_conn%satparam_dn_itype   (iconn) = soil_pp%satfunc_type    (o_r2s_conn%id_dn(iconn) )
           o_r2s_conn%satparam_dn_param_1 (iconn) = soil_pp%alpha           (o_r2s_conn%id_dn(iconn) )
           o_r2s_conn%satparam_dn_param_2 (iconn) = soil_pp%lambda          (o_r2s_conn%id_dn(iconn) )
           o_r2s_conn%satparam_dn_param_3 (iconn) = soil_pp%residual_sat    (o_r2s_conn%id_dn(iconn) )
           ! Relative Perm. dn: MUALEM
           o_r2s_conn%relperm_dn_itype    (iconn) = soil_pp%relperm_type    (o_r2s_conn%id_dn(iconn) )
           o_r2s_conn%relperm_dn_param_1  (iconn) = soil_pp%alpha           (o_r2s_conn%id_dn(iconn) )
           o_r2s_conn%relperm_dn_param_2  (iconn) = soil_pp%lambda          (o_r2s_conn%id_dn(iconn) )
           o_r2s_conn%relperm_dn_param_3  (iconn) = soil_pp%residual_sat    (o_r2s_conn%id_dn(iconn) )
        end do
     end do
  end do

  ! Xylem-to-Root: Conductance flux type
  o_x2r_conn%up_eqn_id   = o_xylem_eqn_id
  o_x2r_conn%dn_eqn_id   = o_root_eqn_id
  o_x2r_conn%up_eqn_name = 'Overstory xylem'
  o_x2r_conn%dn_eqn_name = 'Overstory root'
  iconn = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, overstory_root_nz
           iconn                 = iconn + 1
           o_x2r_conn%id_up(iconn)     = o_xylem_mesh%id(ii,jj,1)
           o_x2r_conn%id_dn(iconn)     = o_root_mesh%id (ii,jj,kk)
           o_x2r_conn%dist_up(iconn)   = 0.1d0
           o_x2r_conn%dist_dn(iconn)   = 0.1d0
           o_x2r_conn%area(iconn)      = overstory_area_sapwood
           o_x2r_conn%itype(iconn)      = CONN_VERTICAL
           o_x2r_conn%flux_type(iconn) = CONDUCTANCE_FLUX_TYPE
           o_x2r_conn%cond_type(iconn) = CONDUCTANCE_CAMPBELL_TYPE
           o_x2r_conn%cond(iconn)      = overstory_root_conductance
           o_x2r_conn%cond_up(iconn)   = overstory_root_conductance
           o_x2r_conn%cond_dn(iconn)   = overstory_root_conductance

           ! Set unit vector
           dist_x = o_root_mesh%xc( o_x2r_conn%id_dn(iconn) - root_id_offset ) - &
                    o_xylem_mesh%xc(o_x2r_conn%id_up(iconn) - xylem_id_offset)
           dist_y = o_root_mesh%yc( o_x2r_conn%id_dn(iconn) - root_id_offset ) - &
                    o_xylem_mesh%yc(o_x2r_conn%id_up(iconn) - xylem_id_offset)
           dist_z = o_root_mesh%zc( o_x2r_conn%id_dn(iconn) - root_id_offset ) - &
                    o_xylem_mesh%zc(o_x2r_conn%id_up(iconn) - xylem_id_offset)
           dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0
           o_x2r_conn%unit_vec(iconn,1) = dist_x/dist
           o_x2r_conn%unit_vec(iconn,2) = dist_y/dist
           o_x2r_conn%unit_vec(iconn,3) = dist_z/dist
           
           ! Saturation up: CHUANG
           o_x2r_conn%satparam_up_itype   (iconn) = o_xylem_pp%satfunc_type    (o_x2r_conn%id_up(iconn) - xylem_id_offset )
           o_x2r_conn%satparam_up_param_1 (iconn) = o_xylem_pp%alpha           (o_x2r_conn%id_up(iconn) - xylem_id_offset )
           o_x2r_conn%satparam_up_param_2 (iconn) = o_xylem_pp%lambda          (o_x2r_conn%id_up(iconn) - xylem_id_offset )
           o_x2r_conn%satparam_up_param_3 (iconn) = o_xylem_pp%residual_sat    (o_x2r_conn%id_up(iconn) - xylem_id_offset )
           ! Relative Perm. up: WEIBULL
           o_x2r_conn%relperm_up_itype    (iconn) = o_xylem_pp%relperm_type    (o_x2r_conn%id_up(iconn) - xylem_id_offset )
           o_x2r_conn%relperm_up_param_1  (iconn) = o_xylem_pp%relperm_param_1 (o_x2r_conn%id_up(iconn) - xylem_id_offset )
           o_x2r_conn%relperm_up_param_2  (iconn) = o_xylem_pp%relperm_param_2 (o_x2r_conn%id_up(iconn) - xylem_id_offset )

           ! Saturation dn: CHUANG
           o_x2r_conn%satparam_dn_itype   (iconn) = o_root_pp%satfunc_type    (o_x2r_conn%id_dn(iconn) - root_id_offset )
           o_x2r_conn%satparam_dn_param_1 (iconn) = o_root_pp%alpha           (o_x2r_conn%id_dn(iconn) - root_id_offset )
           o_x2r_conn%satparam_dn_param_2 (iconn) = o_root_pp%lambda          (o_x2r_conn%id_dn(iconn) - root_id_offset )
           o_x2r_conn%satparam_dn_param_3 (iconn) = o_root_pp%residual_sat    (o_x2r_conn%id_dn(iconn) - root_id_offset )
           ! Relative Perm. dn: WEIBULL
           o_x2r_conn%relperm_dn_itype    (iconn) = o_root_pp%relperm_type    (o_x2r_conn%id_dn(iconn) - root_id_offset )
           o_x2r_conn%relperm_dn_param_1  (iconn) = o_root_pp%relperm_param_1 (o_x2r_conn%id_dn(iconn) - root_id_offset )
           o_x2r_conn%relperm_dn_param_2  (iconn) = o_root_pp%relperm_param_2 (o_x2r_conn%id_dn(iconn) - root_id_offset )
        end do
     end do
  end do

  ! Xylem-to-Xylem: Darcy flux
  iconn = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny
        
        do kk = 1, overstory_xylem_nz-1
           iconn                 = iconn + 1
           o_x2x_conn%id_up(iconn)     = o_xylem_mesh%id(ii,jj,kk  )
           o_x2x_conn%id_dn(iconn)     = o_xylem_mesh%id(ii,jj,kk+1)
           o_x2x_conn%dist_up(iconn)   = 0.5d0*soil_dz
           o_x2x_conn%dist_dn(iconn)   = 0.5d0*soil_dz
           o_x2x_conn%area(iconn)      = overstory_area_sapwood
           o_x2x_conn%itype(iconn)      = CONN_VERTICAL
           o_x2x_conn%flux_type(iconn) = DARCY_FLUX_TYPE
        end do
     end do
  end do

  ! Xylem-to-Leaf
  o_x2l_conn%up_eqn_id   = o_xylem_eqn_id
  o_x2l_conn%dn_eqn_id   = o_leaf_eqn_id
  o_x2l_conn%up_eqn_name = 'Overstory xylem'
  o_x2l_conn%dn_eqn_name = 'Overstory leaf'
  iconn = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, overstory_leaf_nz
           idx                   = overstory_branch_2_xylem_index(kk)
           
           iconn                       = iconn + 1
           o_x2l_conn%id_up(iconn)     = o_xylem_mesh%id(ii,jj,idx)
           o_x2l_conn%id_dn(iconn)     = o_leaf_mesh%id (ii,jj,kk )
           o_x2l_conn%dist_up(iconn)   = 0.5d0*overstory_branch_length_profile(idx)
           o_x2l_conn%dist_dn(iconn)   = 0.5d0*overstory_branch_length_profile(idx)
           o_x2l_conn%area(iconn)      = overstory_xylem_area_profile(idx)*overstory_branch_area_ratio
           o_x2l_conn%itype(iconn)     = CONN_HORIZONTAL
           o_x2l_conn%flux_type(iconn) = DARCY_FLUX_TYPE

           ! Set unit vector
           dist_x = o_leaf_mesh%xc( o_x2l_conn%id_dn(iconn) - leaf_id_offset ) - &
                    o_xylem_mesh%xc(o_x2l_conn%id_up(iconn) - xylem_id_offset)
           dist_y = o_leaf_mesh%yc( o_x2l_conn%id_dn(iconn) - leaf_id_offset ) - &
                    o_xylem_mesh%yc(o_x2l_conn%id_up(iconn) - xylem_id_offset)
           dist_z = o_leaf_mesh%zc( o_x2l_conn%id_dn(iconn) - leaf_id_offset ) - &
                    o_xylem_mesh%zc(o_x2l_conn%id_up(iconn) - xylem_id_offset)
           dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0
           o_x2l_conn%unit_vec(iconn,1) = dist_x/dist
           o_x2l_conn%unit_vec(iconn,2) = dist_y/dist
           o_x2l_conn%unit_vec(iconn,3) = dist_z/dist

        end do

     end do
  end do

end subroutine setup_overstory_mesh

!------------------------------------------------------------------------
subroutine setup_understory_mesh()
  !
#include <petsc/finclude/petsc.h>
  !
  use mpp_varcon                , only : grav
  use mpp_varcon                , only : denh2o
  use SaturationFunction        , only : RELPERM_FUNC_WEIBULL
  use SaturationFunction        , only : SAT_FUNC_CHUANG
  use MultiPhysicsProbConstants , only : CONN_VERTICAL
  use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
  use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
  use MultiPhysicsProbConstants , only : CONDUCTANCE_CAMPBELL_TYPE
  use MultiPhysicsProbConstants , only : CONDUCTANCE_MANOLI_TYPE
  use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
  use petscsys
  use problem_parameters
  !
  implicit none
  !
  PetscInt  :: ii,jj,kk
  PetscInt  :: idx
  PetscInt  :: id_value
  PetscInt  :: count
  PetscInt  :: iconn
  PetscReal :: zz
  PetscReal :: soil_volume
  PetscReal :: root_length
  PetscInt  :: root_id_offset
  PetscInt  :: xylem_id_offset
  PetscInt  :: leaf_id_offset
  PetscReal :: dist_x, dist_y, dist_z, dist

  understory_leaf_nz  = 0
  understory_xylem_nz = int(understory_canopy_height/soil_dz)
  understory_root_nz  = int(understory_root_depth   /soil_dz)

  allocate (understory_xylem_area_profile    (understory_xylem_nz))
  allocate (understory_xylem_has_branches    (understory_xylem_nz))
  allocate (understory_branch_length_profile (understory_xylem_nz))
  allocate (understory_branch_2_xylem_index  (understory_xylem_nz))
  allocate (understory_root_length_profile   (understory_root_nz ))
  allocate (understory_root_area_profile     (understory_root_nz ))
  allocate (understory_root_vol_profile      (understory_root_nz ))

  understory_branch_length_profile(:) = 0.d0
  count = 0

  do kk = 1, understory_xylem_nz
     zz = (kk-1)*soil_dz + soil_dz/2.d0
     understory_xylem_area_profile(kk) = understory_area_sapwood * &
          ( 1.d0 - understory_area_tappering_coeff * zz/understory_canopy_height) ** 2.d0

     if (understory_lad_profile(kk) > 0.d0) then
        count                                = count + 1
        understory_leaf_nz                   = understory_leaf_nz + 1;
        understory_xylem_has_branches(kk)    = PETSC_TRUE
        understory_branch_length_profile(kk) = understory_xylem_area_profile(kk) * &
                                               understory_branch_area_ratio
        understory_branch_2_xylem_index(count) = kk
     else
        understory_xylem_has_branches(kk) = PETSC_FALSE
     end if
     
  end do

  soil_volume = soil_dx * soil_dy * soil_dz
  do kk = 1, understory_root_nz
     root_length                        = understory_B_profile(kk) * soil_volume
     understory_root_length_profile(kk) = root_length
     understory_root_area_profile(kk)   = 2*PI*understory_root_radius*root_length
     understory_root_vol_profile(kk)    = PI*(overstory_root_radius**2.d0)*root_length
  end do

  call u_root_mesh%Init(soil_nx*soil_ny*understory_root_nz)
  call u_root_pp%Init(  soil_nx*soil_ny*understory_root_nz)

  call u_xylem_mesh%Init(soil_nx*soil_ny*understory_xylem_nz)
  call u_xylem_pp%Init(  soil_nx*soil_ny*understory_xylem_nz)

  call u_leaf_mesh%Init(soil_nx*soil_ny*understory_leaf_nz)
  call u_leaf_pp%Init(  soil_nx*soil_ny*understory_leaf_nz)

  allocate(u_root_mesh%id  (soil_nx, soil_ny, understory_root_nz  ))
  allocate(u_xylem_mesh%id (soil_nx, soil_ny, understory_xylem_nz ))
  allocate(u_leaf_mesh%id  (soil_nx, soil_ny, understory_leaf_nz  ))
  
  if (multi_goveqns_formulation) then
     id_value = 0
  else
     id_value = soil_ncells + soil_nx * soil_ny *( overstory_root_nz  + overstory_xylem_nz  + overstory_leaf_nz)
  end if
  root_id_offset = id_value

  ! Add root grid cells
  count = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, understory_root_nz
           id_value                          = id_value + 1
           u_root_mesh%id        (ii,jj,kk)  = id_value

           count                             = count + 1
           u_root_mesh%xc            (count) = soil_xc3d(ii,jj,1) + 0.1d0
           u_root_mesh%yc            (count) = soil_yc3d(ii,jj,1)
           u_root_mesh%zc            (count) = elevation(ii,jj) - soil_dz/2.d0 - soil_dz*(kk-1)
           u_root_mesh%dx            (count) = soil_dx
           u_root_mesh%dy            (count) = soil_dy
           u_root_mesh%dz            (count) = soil_dz
           u_root_mesh%area          (count) = understory_root_area_profile(kk)
           u_root_mesh%vol           (count) = understory_root_vol_profile(kk)
           u_root_mesh%filter        (count) = 1

           u_root_pp%por             (count) = 0.d0
           u_root_pp%perm            (count) = 0.d0
           u_root_pp%alpha           (count) = understory_xylem_phi0
           u_root_pp%lambda          (count) = understory_xylem_p
           u_root_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           u_root_pp%relperm_param_1 (count) = understory_xylem_vulnerability_d * grav * denh2o
           u_root_pp%relperm_param_2 (count) = understory_xylem_vulnerability_c
           u_root_pp%residual_sat    (count) = 0.d0
           u_root_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  ! Add xylem grid cells
  count = 0
  if (multi_goveqns_formulation) id_value = 0
  xylem_id_offset = id_value

  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, understory_xylem_nz
           id_value                           = id_value + 1
           u_xylem_mesh%id         (ii,jj,kk) = id_value

           count                              = count + 1
           u_xylem_mesh%xc            (count) = soil_xc3d(ii,jj,1) + 0.1d0
           u_xylem_mesh%yc            (count) = soil_yc3d(ii,jj,1)
           u_xylem_mesh%zc            (count) = elevation(ii,jj) + soil_dz/2.d0 + (kk-1)*soil_dz
           u_xylem_mesh%dx            (count) = soil_dx
           u_xylem_mesh%dy            (count) = soil_dy
           u_xylem_mesh%dz            (count) = soil_dz
           u_xylem_mesh%area          (count) = understory_xylem_area_profile(kk)
           u_xylem_mesh%vol           (count) = understory_xylem_area_profile(kk) * soil_dz
           u_xylem_mesh%filter        (count) = 1

           u_xylem_pp%por             (count) = understory_xylem_porosity
           u_xylem_pp%perm            (count) = understory_xylem_Kmax * vish2o / (denh2o * grav)
           u_xylem_pp%alpha           (count) = understory_xylem_phi0
           u_xylem_pp%lambda          (count) = understory_xylem_p
           u_xylem_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           u_xylem_pp%relperm_param_1 (count) = understory_xylem_vulnerability_d * grav * denh2o
           u_xylem_pp%relperm_param_2 (count) = understory_xylem_vulnerability_c
           u_xylem_pp%residual_sat    (count) = 0.d0
           u_xylem_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  ! Add branch grid cells
  count = 0
  if (multi_goveqns_formulation) id_value = 0
  leaf_id_offset = id_value

  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, understory_leaf_nz
           id_value                          = id_value + 1
           u_leaf_mesh%id        (ii,jj,kk)  = id_value

           count                             = count + 1
           idx                               = understory_branch_2_xylem_index(kk)

           u_leaf_mesh%xc            (count) = soil_xc3d(ii,jj,1)-understory_branch_length_profile(idx) + 0.1d0
           u_leaf_mesh%yc            (count) = soil_yc3d(ii,jj,1)
           u_leaf_mesh%zc            (count) = elevation(ii,jj) + soil_dz/2.d0 + (kk-1)*soil_dz + (understory_xylem_nz-understory_leaf_nz)*soil_dz
           u_leaf_mesh%dx            (count) = soil_dx
           u_leaf_mesh%dy            (count) = soil_dy
           u_leaf_mesh%dz            (count) = soil_dz
           u_leaf_mesh%area          (count) = understory_xylem_area_profile(kk) * understory_branch_area_ratio

           u_leaf_mesh%vol           (count) = u_leaf_mesh%area(count) * understory_branch_length_profile(idx)
           u_leaf_mesh%filter        (count) = 1

           u_leaf_pp%por             (count) = 0.d0
           u_leaf_pp%perm            (count) = understory_xylem_Kmax * vish2o / (denh2o * grav)
           u_leaf_pp%alpha           (count) = understory_xylem_phi0
           u_leaf_pp%lambda          (count) = understory_xylem_p
           u_leaf_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           u_leaf_pp%relperm_param_1 (count) = understory_xylem_vulnerability_d * grav * denh2o
           u_leaf_pp%relperm_param_2 (count) = understory_xylem_vulnerability_c
           u_leaf_pp%residual_sat    (count) = 0.d0
           u_leaf_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  call u_r2s_conn%Init   (understory_root_nz*        soil_nx * soil_ny)
  call u_x2r_conn%Init  (understory_root_nz*        soil_nx * soil_ny)
  call u_x2x_conn%Init ((understory_xylem_nz - 1)* soil_nx * soil_ny)
  call u_x2l_conn%Init  (understory_leaf_nz*        soil_nx * soil_ny)

  ! Root-to-Soil: Condunctance flux type
  u_r2s_conn%up_eqn_id   = u_root_eqn_id
  u_r2s_conn%dn_eqn_id   = soil_eqn_id
  u_r2s_conn%up_eqn_name = 'Understory root'
  u_r2s_conn%dn_eqn_name = 'Soil'
  iconn = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, understory_root_nz
           iconn                 = iconn + 1
           u_r2s_conn%id_up(iconn)     = u_root_mesh%id(ii,jj,kk)
           u_r2s_conn%id_dn(iconn)     = soil_id(ii,jj,kk-1+top_active_layer_kk_index(ii,jj))
           u_r2s_conn%dist_up(iconn)   = 0.d0
           u_r2s_conn%dist_dn(iconn)   = understory_root_length_profile(kk)
           u_r2s_conn%area(iconn)      = understory_root_area_profile  (kk)
           u_r2s_conn%itype(iconn)      = CONN_HORIZONTAL
           u_r2s_conn%flux_type(iconn) = CONDUCTANCE_FLUX_TYPE
           u_r2s_conn%cond_type(iconn) = CONDUCTANCE_MANOLI_TYPE
           u_r2s_conn%cond_up(iconn)   = understory_root_conductance
           u_r2s_conn%cond_dn(iconn)   = perm_z_top /vish2o * (denh2o * grav) / & ! [m/s]
                understory_root_length_profile(kk)       ! [m]

           ! Set unit vector
           dist_x = soil_mesh%xc(u_r2s_conn%id_dn(iconn)) - u_root_mesh%xc(u_r2s_conn%id_up(iconn) - root_id_offset)
           dist_y = soil_mesh%yc(u_r2s_conn%id_dn(iconn)) - u_root_mesh%yc(u_r2s_conn%id_up(iconn) - root_id_offset)
           dist_z = soil_mesh%zc(u_r2s_conn%id_dn(iconn)) - u_root_mesh%zc(u_r2s_conn%id_up(iconn) - root_id_offset)

           dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0
           u_r2s_conn%unit_vec(iconn,1) = dist_x/dist
           u_r2s_conn%unit_vec(iconn,2) = dist_y/dist
           u_r2s_conn%unit_vec(iconn,3) = dist_z/dist
           
           ! Saturation up: CHUANG
           u_r2s_conn%satparam_up_itype   (iconn) = u_root_pp%satfunc_type    (u_r2s_conn%id_up(iconn) - root_id_offset )
           u_r2s_conn%satparam_up_param_1 (iconn) = u_root_pp%alpha           (u_r2s_conn%id_up(iconn) - root_id_offset )
           u_r2s_conn%satparam_up_param_2 (iconn) = u_root_pp%lambda          (u_r2s_conn%id_up(iconn) - root_id_offset )
           u_r2s_conn%satparam_up_param_3 (iconn) = u_root_pp%residual_sat    (u_r2s_conn%id_up(iconn) - root_id_offset )
           ! Relative Perm. up: MUALEM
           u_r2s_conn%relperm_up_itype    (iconn) = u_root_pp%relperm_type    (u_r2s_conn%id_up(iconn) - root_id_offset )
           u_r2s_conn%relperm_up_param_1  (iconn) = u_root_pp%relperm_param_1 (u_r2s_conn%id_up(iconn) - root_id_offset )
           u_r2s_conn%relperm_up_param_2  (iconn) = u_root_pp%relperm_param_2 (u_r2s_conn%id_up(iconn) - root_id_offset )

           ! Saturation dn: CHUANG
           u_r2s_conn%satparam_dn_itype   (iconn) = soil_pp%satfunc_type    (u_r2s_conn%id_dn(iconn) )
           u_r2s_conn%satparam_dn_param_1 (iconn) = soil_pp%alpha           (u_r2s_conn%id_dn(iconn) )
           u_r2s_conn%satparam_dn_param_2 (iconn) = soil_pp%lambda          (u_r2s_conn%id_dn(iconn) )
           u_r2s_conn%satparam_dn_param_3 (iconn) = soil_pp%residual_sat    (u_r2s_conn%id_dn(iconn) )
           ! Relative Perm. dn: MUALEM
           u_r2s_conn%relperm_dn_itype    (iconn) = soil_pp%relperm_type    (u_r2s_conn%id_dn(iconn) )
           u_r2s_conn%relperm_dn_param_1  (iconn) = soil_pp%alpha           (u_r2s_conn%id_dn(iconn) )
           u_r2s_conn%relperm_dn_param_2  (iconn) = soil_pp%lambda          (u_r2s_conn%id_dn(iconn) )
           u_r2s_conn%relperm_dn_param_3  (iconn) = soil_pp%residual_sat    (u_r2s_conn%id_dn(iconn) )
        end do
     end do
  end do

  ! Xylem-to-Root: Conductance flux type
  iconn = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

      do kk = 1, understory_root_nz
           iconn                 = iconn + 1
           u_x2r_conn%id_up(iconn)     = u_xylem_mesh%id(ii,jj,1)
           u_x2r_conn%id_dn(iconn)     = u_root_mesh%id (ii,jj,kk)
           u_x2r_conn%dist_up(iconn)   = 0.1d0
           u_x2r_conn%dist_dn(iconn)   = 0.1d0
           u_x2r_conn%area(iconn)      = understory_area_sapwood
           u_x2r_conn%itype(iconn)      = CONN_VERTICAL
           u_x2r_conn%flux_type(iconn) = CONDUCTANCE_FLUX_TYPE
           u_x2r_conn%cond_type(iconn) = CONDUCTANCE_CAMPBELL_TYPE
           u_x2r_conn%cond_up          = understory_root_conductance
           u_x2r_conn%cond_up(iconn)   = understory_root_conductance
           u_x2r_conn%cond_dn(iconn)   = understory_root_conductance

           ! Set unit vector
           dist_x = u_root_mesh%xc( u_x2r_conn%id_dn(iconn) - root_id_offset ) - &
                    u_xylem_mesh%xc(u_x2r_conn%id_up(iconn) - xylem_id_offset)
           dist_y = u_root_mesh%yc( u_x2r_conn%id_dn(iconn) - root_id_offset ) - &
                    u_xylem_mesh%yc(u_x2r_conn%id_up(iconn) - xylem_id_offset)
           dist_z = u_root_mesh%zc( u_x2r_conn%id_dn(iconn) - root_id_offset ) - &
                    u_xylem_mesh%zc(u_x2r_conn%id_up(iconn) - xylem_id_offset)
           dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0
           u_x2r_conn%unit_vec(iconn,1) = dist_x/dist
           u_x2r_conn%unit_vec(iconn,2) = dist_y/dist
           u_x2r_conn%unit_vec(iconn,3) = dist_z/dist

           ! Saturation up: CHUANG
           u_x2r_conn%satparam_up_itype   (iconn) = u_xylem_pp%satfunc_type    (u_x2r_conn%id_up(iconn) - xylem_id_offset )
           u_x2r_conn%satparam_up_param_1 (iconn) = u_xylem_pp%alpha           (u_x2r_conn%id_up(iconn) - xylem_id_offset )
           u_x2r_conn%satparam_up_param_2 (iconn) = u_xylem_pp%lambda          (u_x2r_conn%id_up(iconn) - xylem_id_offset )
           u_x2r_conn%satparam_up_param_3 (iconn) = u_xylem_pp%residual_sat    (u_x2r_conn%id_up(iconn) - xylem_id_offset )
           ! Relative Perm. up: WEIBULL
           u_x2r_conn%relperm_up_itype    (iconn) = u_xylem_pp%relperm_type    (u_x2r_conn%id_up(iconn) - xylem_id_offset )
           u_x2r_conn%relperm_up_param_1  (iconn) = u_xylem_pp%relperm_param_1 (u_x2r_conn%id_up(iconn) - xylem_id_offset )
           u_x2r_conn%relperm_up_param_2  (iconn) = u_xylem_pp%relperm_param_2 (u_x2r_conn%id_up(iconn) - xylem_id_offset )

           ! Saturation dn: CHUANG
           u_x2r_conn%satparam_dn_itype   (iconn) = u_root_pp%satfunc_type    (u_x2r_conn%id_dn(iconn) - root_id_offset )
           u_x2r_conn%satparam_dn_param_1 (iconn) = u_root_pp%alpha           (u_x2r_conn%id_dn(iconn) - root_id_offset )
           u_x2r_conn%satparam_dn_param_2 (iconn) = u_root_pp%lambda          (u_x2r_conn%id_dn(iconn) - root_id_offset )
           u_x2r_conn%satparam_dn_param_3 (iconn) = u_root_pp%residual_sat    (u_x2r_conn%id_dn(iconn) - root_id_offset )
           ! Relative Perm. dn: WEIBULL
           u_x2r_conn%relperm_dn_itype    (iconn) = u_root_pp%relperm_type    (u_x2r_conn%id_dn(iconn) - root_id_offset )
           u_x2r_conn%relperm_dn_param_1  (iconn) = u_root_pp%relperm_param_1 (u_x2r_conn%id_dn(iconn) - root_id_offset )
           u_x2r_conn%relperm_dn_param_2  (iconn) = u_root_pp%relperm_param_2 (u_x2r_conn%id_dn(iconn) - root_id_offset )
        end do
     end do
  end do

  ! Xylem-to-Xylem: Darcy flux
  u_x2r_conn%up_eqn_id   = u_xylem_eqn_id
  u_x2r_conn%dn_eqn_id   = u_root_eqn_id
  u_x2r_conn%up_eqn_name = 'Understory xylem'
  u_x2r_conn%dn_eqn_name = 'Understory root'
  iconn = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk                          = 1, understory_xylem_nz-1
           iconn                       = iconn + 1
           u_x2x_conn%id_up(iconn)     = u_xylem_mesh%id(ii,jj,kk  )
           u_x2x_conn%id_dn(iconn)     = u_xylem_mesh%id(ii,jj,kk+1)
           u_x2x_conn%dist_up(iconn)   = 0.5d0*soil_dz
           u_x2x_conn%dist_dn(iconn)   = 0.5d0*soil_dz
           u_x2x_conn%area(iconn)      = understory_area_sapwood
           u_x2x_conn%itype(iconn)     = CONN_VERTICAL
           u_x2x_conn%flux_type(iconn) = DARCY_FLUX_TYPE
        end do
     end do
  end do

  ! Xylem-to-Leaf
  u_x2l_conn%up_eqn_id   = u_xylem_eqn_id
  u_x2l_conn%dn_eqn_id   = u_leaf_eqn_id
  u_x2l_conn%up_eqn_name = 'Understory xylem'
  u_x2l_conn%dn_eqn_name = 'Understory leaf'
  iconn = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, understory_leaf_nz
           idx                   = understory_branch_2_xylem_index(kk)
           
           iconn                 = iconn + 1
           u_x2l_conn%id_up(iconn)     = u_xylem_mesh%id(ii,jj,idx)
           u_x2l_conn%id_dn(iconn)     = u_leaf_mesh%id (ii,jj,kk )
           u_x2l_conn%dist_up(iconn)   = 0.5d0*understory_branch_length_profile(idx)
           u_x2l_conn%dist_dn(iconn)   = 0.5d0*understory_branch_length_profile(idx)

           u_x2l_conn%area(iconn)      = understory_xylem_area_profile(idx)*understory_branch_area_ratio
           u_x2l_conn%itype(iconn)      = CONN_HORIZONTAL
           u_x2l_conn%flux_type(iconn) = DARCY_FLUX_TYPE


           ! Set unit vector
           dist_x = u_leaf_mesh%xc( u_x2l_conn%id_dn(iconn) - leaf_id_offset ) - &
                    u_xylem_mesh%xc(u_x2l_conn%id_up(iconn) - xylem_id_offset)
           dist_y = u_leaf_mesh%yc( u_x2l_conn%id_dn(iconn) - leaf_id_offset ) - &
                    u_xylem_mesh%yc(u_x2l_conn%id_up(iconn) - xylem_id_offset)
           dist_z = u_leaf_mesh%zc( u_x2l_conn%id_dn(iconn) - leaf_id_offset ) - &
                    u_xylem_mesh%zc(u_x2l_conn%id_up(iconn) - xylem_id_offset)
           dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0
           u_x2l_conn%unit_vec(iconn,1) = dist_x/dist
           u_x2l_conn%unit_vec(iconn,2) = dist_y/dist
           u_x2l_conn%unit_vec(iconn,3) = dist_z/dist

        end do

     end do
  end do

end subroutine setup_understory_mesh

!------------------------------------------------------------------------
subroutine add_single_goveqn()
  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
  use MultiPhysicsProbVSFM , only : vsfm_mpp
  use MultiPhysicsProbConstants , only : GE_RE
  use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
  !
  ! !ARGUMENTS
  implicit none

  call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE', MESH_CLM_SOIL_COL)

  call vsfm_mpp%SetMeshesOfGoveqns()

end subroutine add_single_goveqn

!------------------------------------------------------------------------
subroutine add_multiple_goveqns()
  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
  use MultiPhysicsProbVSFM , only : vsfm_mpp
  use MultiPhysicsProbConstants , only : GE_RE
  !
  ! !ARGUMENTS
  implicit none

  call vsfm_mpp%AddGovEqnWithMeshRank(GE_RE, 'Richards Equation ODE for soil'             , 1)
  call vsfm_mpp%AddGovEqnWithMeshRank(GE_RE, 'Richards Equation ODE for overstory root'   , 2)
  call vsfm_mpp%AddGovEqnWithMeshRank(GE_RE, 'Richards Equation ODE for overstory xylem'  , 3)
  call vsfm_mpp%AddGovEqnWithMeshRank(GE_RE, 'Richards Equation ODE for overstory leaf'   , 4)
  call vsfm_mpp%AddGovEqnWithMeshRank(GE_RE, 'Richards Equation ODE for understory root'  , 5)
  call vsfm_mpp%AddGovEqnWithMeshRank(GE_RE, 'Richards Equation ODE for understory xylem' , 6)
  call vsfm_mpp%AddGovEqnWithMeshRank(GE_RE, 'Richards Equation ODE for understory leaf'  , 7)

  call vsfm_mpp%SetMeshesOfGoveqnsByMeshRank()

end subroutine add_multiple_goveqns

!------------------------------------------------------------------------
subroutine add_conditions_to_goveqns()
  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use MultiPhysicsProbConstants , only : ALL_CELLS
  use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
  use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
  use MultiPhysicsProbConstants , only : COND_BC
  use MultiPhysicsProbConstants , only : COND_SS
  use MultiPhysicsProbConstants , only : COND_DIRICHLET
  use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_FETCH2
  use MultiPhysicsProbConstants , only : CONN_VERTICAL
  use ConnectionSetType         , only : connection_set_type
  use ConnectionSetType         , only : ConnectionSetDestroy
  use MeshType                  , only : MeshCreateConnectionSet
  use petscsys
  use problem_parameters
  !
  ! !ARGUMENTS
  implicit none
  !
  PetscInt                            :: ieqn

  if (multi_goveqns_formulation) call add_coupling_conditions_to_goveqns()
  
end subroutine add_conditions_to_goveqns

!------------------------------------------------------------------------
subroutine add_coupling_conditions_to_goveqns()
  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use MultiPhysicsProbConstants , only : ALL_CELLS
  use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
  use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
  use MultiPhysicsProbConstants , only : COND_BC
  use MultiPhysicsProbConstants , only : COND_SS
  use MultiPhysicsProbConstants , only : COND_DIRICHLET
  use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_FETCH2
  use MultiPhysicsProbConstants , only : CONN_VERTICAL
  use ConnectionSetType         , only : connection_set_type
  use ConnectionSetType         , only : ConnectionSetDestroy
  use MeshType                  , only : MeshCreateConnectionSet
  use problem_parameters
  use petscsys
  !
  ! !ARGUMENTS
  implicit none
  !

  call o_r2s_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_FALSE)
  call o_r2s_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_TRUE)
  
  call u_r2s_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_FALSE)
  call u_r2s_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_TRUE)

  call o_x2r_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_FALSE)
  call o_x2r_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_TRUE)

  call o_x2l_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_FALSE)
  call o_x2l_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_TRUE)

  call u_x2r_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_FALSE)
  call u_x2r_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_TRUE)

  call u_x2l_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_FALSE)
  call u_x2l_conn%AddToMPPCouplingBC(vsfm_mpp, reverse_up2dn = PETSC_TRUE)

end subroutine add_coupling_conditions_to_goveqns

!------------------------------------------------------------------------
subroutine allocate_auxvars()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbVSFM     , only : vsfm_mpp
  use MultiPhysicsProbConstants, only : VAR_PRESSURE
  use problem_parameters
  !
  implicit none
  !
  integer           :: ieqn
  integer           :: nvars_for_coupling
  integer, pointer  :: var_ids_for_coupling(:)
  integer, pointer  :: goveqn_ids_for_coupling(:)
  integer, pointer  :: is_bc(:)

  !
  ! Allocate auxvars
  !
  if (multi_goveqns_formulation) then
     !
     ! Soil <---> Overstory root
     ! Soil <---> Understory root
     !
     ieqn               = soil_eqn_id
     nvars_for_coupling = 2

     allocate (var_ids_for_coupling    (nvars_for_coupling))
     allocate (goveqn_ids_for_coupling (nvars_for_coupling))

     var_ids_for_coupling    (1) = VAR_PRESSURE
     var_ids_for_coupling    (2) = VAR_PRESSURE
     goveqn_ids_for_coupling (1) = o_root_eqn_id
     goveqn_ids_for_coupling (2) = u_root_eqn_id

     call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
          var_ids_for_coupling, goveqn_ids_for_coupling)
       
     deallocate(var_ids_for_coupling   )
     deallocate(goveqn_ids_for_coupling)

     !
     ! Overstory root <---> Soil
     ! Overstory root <---> Overstory xylem
     !
     ieqn               = o_root_eqn_id
     nvars_for_coupling = 2

     allocate (var_ids_for_coupling    (nvars_for_coupling))
     allocate (goveqn_ids_for_coupling (nvars_for_coupling))

     var_ids_for_coupling    (1) = VAR_PRESSURE
     var_ids_for_coupling    (2) = VAR_PRESSURE
     goveqn_ids_for_coupling (1) = soil_eqn_id
     goveqn_ids_for_coupling (2) = o_xylem_eqn_id

     call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
          var_ids_for_coupling, goveqn_ids_for_coupling)
       
     deallocate(var_ids_for_coupling   )
     deallocate(goveqn_ids_for_coupling)

     !
     ! Overstory xylem <---> Overstory root
     ! Overstory xylem <---> Overstory leaf
     !
     ieqn               = o_xylem_eqn_id
     nvars_for_coupling = 2

     allocate (var_ids_for_coupling    (nvars_for_coupling))
     allocate (goveqn_ids_for_coupling (nvars_for_coupling))

     var_ids_for_coupling    (1) = VAR_PRESSURE
     var_ids_for_coupling    (2) = VAR_PRESSURE
     goveqn_ids_for_coupling (1) = o_root_eqn_id
     goveqn_ids_for_coupling (2) = o_leaf_eqn_id

     call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
          var_ids_for_coupling, goveqn_ids_for_coupling)
       
     deallocate(var_ids_for_coupling   )
     deallocate(goveqn_ids_for_coupling)

     !
     ! Overstory leaf <---> Overstory xylem
     !
     ieqn               = o_leaf_eqn_id
     nvars_for_coupling = 1

     allocate (var_ids_for_coupling    (nvars_for_coupling))
     allocate (goveqn_ids_for_coupling (nvars_for_coupling))

     var_ids_for_coupling    (1) = VAR_PRESSURE
     goveqn_ids_for_coupling (1) = o_xylem_eqn_id

     call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
          var_ids_for_coupling, goveqn_ids_for_coupling)
       
     deallocate(var_ids_for_coupling   )
     deallocate(goveqn_ids_for_coupling)

     !
     ! Understory root <---> Soil
     ! Understory root <---> Understory xylem
     !
     ieqn               = u_root_eqn_id
     nvars_for_coupling = 2

     allocate (var_ids_for_coupling    (nvars_for_coupling))
     allocate (goveqn_ids_for_coupling (nvars_for_coupling))

     var_ids_for_coupling    (1) = VAR_PRESSURE
     var_ids_for_coupling    (2) = VAR_PRESSURE
     goveqn_ids_for_coupling (1) = soil_eqn_id
     goveqn_ids_for_coupling (2) = u_xylem_eqn_id

     call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
          var_ids_for_coupling, goveqn_ids_for_coupling)
       
     deallocate(var_ids_for_coupling   )
     deallocate(goveqn_ids_for_coupling)

     !
     ! Understory xylem <---> Understory root
     ! Understory xylem <---> Understory leaf
     !
     ieqn               = u_xylem_eqn_id
     nvars_for_coupling = 2

     allocate (var_ids_for_coupling    (nvars_for_coupling))
     allocate (goveqn_ids_for_coupling (nvars_for_coupling))

     var_ids_for_coupling    (1) = VAR_PRESSURE
     var_ids_for_coupling    (2) = VAR_PRESSURE
     goveqn_ids_for_coupling (1) = u_root_eqn_id
     goveqn_ids_for_coupling (2) = u_leaf_eqn_id

     call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
          var_ids_for_coupling, goveqn_ids_for_coupling)
       
     deallocate(var_ids_for_coupling   )
     deallocate(goveqn_ids_for_coupling)

     !
     ! Understory leaf <---> Understory xylem
     !
     ieqn               = u_leaf_eqn_id
     nvars_for_coupling = 1

     allocate (var_ids_for_coupling    (nvars_for_coupling))
     allocate (goveqn_ids_for_coupling (nvars_for_coupling))

     var_ids_for_coupling    (1) = VAR_PRESSURE
     goveqn_ids_for_coupling (1) = u_xylem_eqn_id

     call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
          var_ids_for_coupling, goveqn_ids_for_coupling)
       
     deallocate(var_ids_for_coupling   )
     deallocate(goveqn_ids_for_coupling)

  end if

  call vsfm_mpp%AllocateAuxVars()

end subroutine allocate_auxvars

!------------------------------------------------------------------------
subroutine set_material_properties()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use EOSWaterMod               , only : DENSITY_IFC67
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetDensityType
  use MultiPhysicsProbConstants , only : AUXVAR_CONN_INTERNAL
  use problem_parameters
  use petscsys
  !
#include <petsc/finclude/petsc.h>
  !
  implicit none
  !
  PetscInt                     :: igoveqn
  PetscBool, pointer           :: set_upwind_auxvar(:)
  !-----------------------------------------------------------------------

  igoveqn = 1

  call VSFMMPPSetDensityType(vsfm_mpp, 1, DENSITY_IFC67)

  call set_physical_property_for_eqn(igoveqn, combined_pp  )
  call set_conn_property_for_eqn(    igoveqn, combined_conn, AUXVAR_CONN_INTERNAL)

end subroutine set_material_properties

!------------------------------------------------------------------------
subroutine set_material_properties_multi_goveqns()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use EOSWaterMod               , only : DENSITY_IFC67
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetDensityType
  use MultiPhysicsProbConstants , only : AUXVAR_CONN_INTERNAL
  use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
  use problem_parameters
  use petscsys
  !
#include <petsc/finclude/petsc.h>
  !
  implicit none
  !
  PetscInt                       :: idx_beg, idx_end
  type(spac_component_conn_type) :: soil_bc_conn
  type(spac_component_conn_type) :: o_root_bc_conn
  type(spac_component_conn_type) :: o_xylem_bc_conn
  type(spac_component_conn_type) :: u_root_bc_conn
  type(spac_component_conn_type) :: u_xylem_bc_conn
  !-----------------------------------------------------------------------

  call VSFMMPPSetDensityType(vsfm_mpp, soil_eqn_id, DENSITY_IFC67)
  call VSFMMPPSetDensityType(vsfm_mpp, o_root_eqn_id, DENSITY_IFC67)
  call VSFMMPPSetDensityType(vsfm_mpp, o_xylem_eqn_id, DENSITY_IFC67)
  call VSFMMPPSetDensityType(vsfm_mpp, o_leaf_eqn_id, DENSITY_IFC67)
  call VSFMMPPSetDensityType(vsfm_mpp, u_root_eqn_id, DENSITY_IFC67)
  call VSFMMPPSetDensityType(vsfm_mpp, u_xylem_eqn_id, DENSITY_IFC67)
  call VSFMMPPSetDensityType(vsfm_mpp, u_leaf_eqn_id, DENSITY_IFC67)

  call set_physical_property_for_eqn(soil_eqn_id    , soil_pp    )
  call set_physical_property_for_eqn(o_root_eqn_id  , o_root_pp  )
  call set_physical_property_for_eqn(o_xylem_eqn_id , o_xylem_pp )
  call set_physical_property_for_eqn(o_leaf_eqn_id  , o_leaf_pp  )
  call set_physical_property_for_eqn(u_root_eqn_id  , u_root_pp  )
  call set_physical_property_for_eqn(u_xylem_eqn_id , u_xylem_pp )
  call set_physical_property_for_eqn(u_leaf_eqn_id  , u_leaf_pp  )

  ! Create soil coupling BC connections
  call soil_bc_conn%Init(o_r2s_conn%nconn + u_r2s_conn%nconn)
  idx_beg = 1          ; call soil_bc_conn%Copy(idx_beg, idx_end, o_r2s_conn)
  idx_beg = idx_end + 1; call soil_bc_conn%Copy(idx_beg, idx_end, u_r2s_conn)

  ! Create overstory root coupling BC connections
  call o_root_bc_conn%Init(o_r2s_conn%nconn + o_x2r_conn%nconn)
  idx_beg = 1          ; call o_root_bc_conn%Copy(idx_beg, idx_end, o_r2s_conn, reverse_updn=PETSC_TRUE)
  idx_beg = idx_end + 1; call o_root_bc_conn%Copy(idx_beg, idx_end, o_x2r_conn)

  ! Create overstory xylem coupling BC connections
  call o_xylem_bc_conn%Init(o_x2r_conn%nconn + o_x2l_conn%nconn)
  idx_beg = 1          ; call o_xylem_bc_conn%Copy(idx_beg, idx_end, o_x2r_conn, reverse_updn=PETSC_TRUE)
  idx_beg = idx_end + 1; call o_xylem_bc_conn%Copy(idx_beg, idx_end, o_x2l_conn, reverse_updn=PETSC_TRUE)

  ! Create understory root coupling BC connections
  call u_root_bc_conn%Init(u_r2s_conn%nconn + u_x2r_conn%nconn)
  idx_beg = 1          ; call u_root_bc_conn%Copy(idx_beg, idx_end, u_r2s_conn, reverse_updn=PETSC_TRUE)
  idx_beg = idx_end + 1; call u_root_bc_conn%Copy(idx_beg, idx_end, u_x2r_conn)

  ! Create understory xylem coupling BC connections
  call u_xylem_bc_conn%Init(u_x2r_conn%nconn + u_x2l_conn%nconn)
  idx_beg = 1          ; call u_xylem_bc_conn%Copy(idx_beg, idx_end, u_x2r_conn, reverse_updn=PETSC_TRUE)
  idx_beg = idx_end + 1; call u_xylem_bc_conn%Copy(idx_beg, idx_end, u_x2l_conn, reverse_updn=PETSC_TRUE)


  ! Soil equation
  call set_conn_property_for_eqn(soil_eqn_id    , s2s_conn        , AUXVAR_CONN_INTERNAL)
  call set_conn_property_for_eqn(soil_eqn_id    , soil_bc_conn    , AUXVAR_CONN_BC      )

  ! Overstory root
  call set_conn_property_for_eqn(o_root_eqn_id  , o_root_bc_conn  , AUXVAR_CONN_BC      )

  ! Overstory xylem
  call set_conn_property_for_eqn(o_xylem_eqn_id , o_x2x_conn      , AUXVAR_CONN_INTERNAL)
  call set_conn_property_for_eqn(o_xylem_eqn_id , o_xylem_bc_conn , AUXVAR_CONN_BC      )

  ! Overstory leaf
  call set_conn_property_for_eqn(o_leaf_eqn_id  , o_x2l_conn      , AUXVAR_CONN_BC      )

  ! Understory root
  call set_conn_property_for_eqn(u_root_eqn_id  , u_root_bc_conn  , AUXVAR_CONN_BC      )

  ! Understory xylem
  call set_conn_property_for_eqn(u_xylem_eqn_id , u_x2x_conn      , AUXVAR_CONN_INTERNAL)
  call set_conn_property_for_eqn(u_xylem_eqn_id , u_xylem_bc_conn , AUXVAR_CONN_BC      )

  ! Understory leaf
  call set_conn_property_for_eqn(u_leaf_eqn_id  , u_x2l_conn      , AUXVAR_CONN_BC      )

  ! Free up the memory
  call soil_bc_conn%Destroy()
  call o_root_bc_conn%Destroy()
  call o_xylem_bc_conn%Destroy()
  call u_root_bc_conn%Destroy()
  call u_xylem_bc_conn%Destroy()

end subroutine set_material_properties_multi_goveqns

!------------------------------------------------------------------------
subroutine set_physical_property_for_eqn(eqn_id, pp)
  !
  use MultiPhysicsProbVSFM , only : VSFMMPPSetSoilPorosity
  use MultiPhysicsProbVSFM , only : VSFMMPPSetSaturationFunction
  use MultiPhysicsProbVSFM , only : VSFMMPPSetSoilPermeability
  use MultiPhysicsProbVSFM , only : VSFMMPPSetRelativePermeability
  use MultiPhysicsProbVSFM , only : vsfm_mpp
  use spac_component
  !
  implicit none
  !
  PetscInt                      :: eqn_id
  type (spac_component_pp_type) :: pp

  call VSFMMPPSetSoilPorosity(vsfm_mpp, eqn_id, pp%por)

  call VSFMMPPSetSaturationFunction(vsfm_mpp, eqn_id, pp%satfunc_type, &
         pp%alpha, pp%lambda, pp%residual_sat)

  call VSFMMPPSetSoilPermeability(vsfm_mpp, eqn_id, pp%perm, pp%perm, pp%perm)

  call VSFMMPPSetRelativePermeability(vsfm_mpp, eqn_id, pp%relperm_type, &
       pp%relperm_param_1, pp%relperm_param_2)

end subroutine set_physical_property_for_eqn

!------------------------------------------------------------------------
subroutine set_conn_property_for_eqn(eqn_id, conn, conn_type)
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunctionAuxVarConn
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnIntValue
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnRealValue
  use MultiPhysicsProbConstants , only : VAR_FLUX_TYPE
  use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE
  use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_TYPE
  use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_UP
  use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_DN
  use spac_component
  use petscsys
  !
  implicit none
  !
  PetscInt                        :: eqn_id
  type (spac_component_conn_type) :: conn
  PetscInt                        :: conn_type
  !
  PetscBool, pointer            :: set_upwind_auxvar(:)

  allocate(set_upwind_auxvar(conn%nconn))

  set_upwind_auxvar(:) = PETSC_TRUE
  ! For upwind cells, set saturation parameters
  call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, eqn_id, conn_type, &
       set_upwind_auxvar, conn%satparam_up_itype, conn%satparam_up_param_1, &
       conn%satparam_up_param_2, conn%satparam_up_param_3)

  ! For upwind cells, set relative permeability parameters
  call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, eqn_id, conn_type, &
       set_upwind_auxvar, conn%relperm_up_itype, conn%relperm_up_param_1, &
       conn%relperm_up_param_2, conn%relperm_up_param_2)

  set_upwind_auxvar(:) = PETSC_FALSE
  ! For downwind cells, set saturation parameters
  call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, eqn_id, conn_type, &
       set_upwind_auxvar, conn%satparam_dn_itype, conn%satparam_dn_param_1, &
       conn%satparam_dn_param_2, conn%satparam_dn_param_3)

  ! For downwind cells, set relative permeability parameters
  call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, eqn_id, conn_type, &
       set_upwind_auxvar, conn%relperm_dn_itype, conn%relperm_dn_param_1, &
       conn%relperm_dn_param_2, conn%relperm_dn_param_3)

  ! Set connection flux type
  call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, conn_type, &
       VAR_FLUX_TYPE, conn%flux_type)

  ! Set conductance type
  call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, conn_type, &
       VAR_CONDUCTANCE_TYPE, conn%cond_type)

  call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, eqn_id, conn_type, &
       VAR_CONDUCTANCE, conn%cond)

  call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, eqn_id, conn_type, &
       VAR_CONDUCTANCE_UP, conn%cond_up)

  call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, eqn_id, conn_type, &
       VAR_CONDUCTANCE_DN, conn%cond_dn)

  deallocate(set_upwind_auxvar)

end subroutine set_conn_property_for_eqn

!------------------------------------------------------------------------
subroutine set_initial_conditions()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbVSFM , only : vsfm_mpp
  use problem_parameters
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda
  !
  implicit none
  !
  PetscInt           :: ii,jj,kk
  PetscInt           :: count
  PetscReal          :: wtd_for_column
  PetscBool          :: wtd_for_column_unset
  PetscReal, pointer :: press_ic(:)
  PetscErrorCode     :: ierr

  call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)

  press_ic(:) = 91325.d0
  count = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny
        wtd_for_column_unset = PETSC_TRUE
        do kk = 1, soil_nz
           if (wtd_for_column_unset .and. soil_id(ii,jj,kk) > 0 ) then
              wtd_for_column = -init_wtd - soil_dz*(kk-1)
              wtd_for_column_unset = PETSC_FALSE
           end if
           if (soil_id(ii,jj,kk) > 0) then
              count = count + 1
              press_ic(count) = 101325 + (wtd_for_column - soil_zc3d(ii,jj,kk))* 1000.d0 * 9.81d0
           endif
        end do
     end do
  end do
  call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
  call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
  call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

end subroutine set_initial_conditions
