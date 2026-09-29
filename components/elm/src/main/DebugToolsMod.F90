module DebugToolsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for mass balance debug.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use elm_varcon          , only : spval, ispval
  use elm_varctl          , only : iulog, use_debug
  use ColumnDataType      , only : column_carbon_state, column_carbon_flux
  use VegetationDataType  , only : vegetation_carbon_state,vegetation_carbon_flux
  use VegetationType      , only : veg_pp
  use ColumnType          , only : col_pp  
  use elm_varpar          , only : max_patch_per_col  
  use decompMod           , only : bounds_type
  use elm_time_manager    , only : get_nstep, get_step_size
  use abortutils          , only : endrun
  implicit none
  save
  private

  public :: EnterCMassDebug
  public :: DebugCMassBal
  public :: Init_Debug_tools
  public :: DebugWritep
  type, public :: column_debug_state
     real(r8), pointer :: totpltc_beg           (:) => null()  !total vegetation C
     real(r8), pointer :: totdsomc_beg          (:) => null()  !total decomposible C
     real(r8), pointer :: totprodc_beg          (:) => null()  !
     real(r8), pointer :: seedc_beg             (:) => null()  !
     real(r8), pointer :: ctrunc_beg            (:) => null()  !
     real(r8), pointer :: totpftc_beg           (:) => null()  !
     real(r8), pointer :: cropseedc_deficit_beg (:) => null()  !
            
  contains
     procedure, public :: Init   => col_dbgs_init
     procedure, public :: Clean  => col_dbgs_clean
  end type column_debug_state

  type, public :: pft_debug_state
    real(r8), pointer :: totpftc_beg (:) => null() !
    real(r8), pointer :: totvegc_beg(:) => null() !
    real(r8), pointer :: xsmrpool_beg(:) => null() 
    real(r8), pointer :: ctrunc_beg(:) => null() 
    real(r8), pointer :: dispvegc_beg(:) => null() !
    real(r8), pointer :: storvegc_beg(:) => null() !
    real(r8), pointer :: leafc_beg(:) => null() !
    real(r8), pointer :: frootc_beg(:) => null() !
    real(r8), pointer :: livestemc_beg(:) => null() !
    real(r8), pointer :: deadstemc_beg(:) => null() !
    real(r8), pointer :: livecrootc_beg(:) => null() !
    real(r8), pointer :: deadcrootc_beg(:) => null() !
    real(r8), pointer :: grainc_beg(:) => null() !
  contains
    procedure, public :: Init => pft_dbgs_Init
    procedure, public :: Clean => pft_dbgs_clean    
  end type pft_debug_state
  type(column_debug_state)  , private, target :: col_dbgs   ! column debug state 
  type(pft_debug_state)     , private, target :: pft_dbgs   ! pft debug state    
contains

  subroutine Init_Debug_tools(bounds)

  type(bounds_type), intent(in)    :: bounds


  call col_dbgs%Init(bounds%begc,bounds%endc)

  call pft_dbgs%Init(bounds%begp, bounds%endp)

  end subroutine Init_Debug_tools

!------------------------------------------------------------------------
  subroutine pft_dbgs_Init(this,begp,endp)

    ! !ARGUMENTS:
    class(pft_debug_state) :: this
    integer, intent(in) :: begp,endp
    integer :: p

  allocate(pft_dbgs%totpftc_beg(begp:endp)); pft_dbgs%totpftc_beg(:) = spval
  allocate(pft_dbgs%totvegc_beg(begp:endp)); pft_dbgs%totvegc_beg(:) = spval
  allocate(pft_dbgs%xsmrpool_beg(begp:endp));pft_dbgs%xsmrpool_beg(:)=spval
  allocate(pft_dbgs%ctrunc_beg(begp:endp)); pft_dbgs%ctrunc_beg(:)=spval
  allocate(pft_dbgs%dispvegc_beg(begp:endp));pft_dbgs%dispvegc_beg(:)=spval
  allocate(pft_dbgs%storvegc_beg(begp:endp));pft_dbgs%storvegc_beg(:)=spval
  allocate(pft_dbgs%leafc_beg(begp:endp));pft_dbgs%leafc_beg(:)=spval
  allocate(pft_dbgs%frootc_beg(begp:endp));pft_dbgs%frootc_beg(:)=spval
  allocate(pft_dbgs%livestemc_beg(begp:endp));pft_dbgs%livestemc_beg(:)=spval
  allocate(pft_dbgs%deadstemc_beg(begp:endp));pft_dbgs%deadstemc_beg(:)=spval
  allocate(pft_dbgs%livecrootc_beg(begp:endp));pft_dbgs%livecrootc_beg(:)=spval
  allocate(pft_dbgs%deadcrootc_beg(begp:endp));pft_dbgs%deadcrootc_beg(:)=spval
  allocate(pft_dbgs%grainc_beg(begp:endp));pft_dbgs%grainc_beg(:)=spval

  do p = begp,endp
    pft_dbgs%totpftc_beg(p) = 0._r8
    pft_dbgs%totvegc_beg(p) = 0._r8
    pft_dbgs%xsmrpool_beg(p)= 0._r8
    pft_dbgs%ctrunc_beg(p) = 0._r8
    pft_dbgs%dispvegc_beg(p) = 0._r8
    pft_dbgs%storvegc_beg(p) = 0._r8
    pft_dbgs%leafc_beg(p) = 0._r8
    pft_dbgs%frootc_beg(p) = 0._r8
    pft_dbgs%livestemc_beg(p) = 0._r8
    pft_dbgs%deadstemc_beg(p) = 0._r8
    pft_dbgs%livecrootc_beg(p) = 0._r8
    pft_dbgs%deadcrootc_beg(p) = 0._r8
    pft_dbgs%grainc_beg(p) = 0._r8
  enddo 
  end subroutine pft_dbgs_Init

!------------------------------------------------------------------------
  subroutine pft_dbgs_clean(this)
    !
    ! !ARGUMENTS:
    class(pft_debug_state) :: this
    
  end subroutine pft_dbgs_clean  
  !------------------------------------------------------------------------
  !subroutine to initialize column debug state data structure
  !------------------------------------------------------------------------ 
  subroutine col_dbgs_init(this,begc,endc)
    !
    !
    ! !ARGUMENTS:
    class(column_debug_state) :: this
    integer, intent(in) :: begc,endc
    integer :: c
    
    allocate(this%totpltc_beg          (begc:endc))     ; this%totpltc_beg          (:)     = spval
    allocate(this%totdsomc_beg          (begc:endc))     ; this%totdsomc_beg          (:)     = spval
    allocate(this%totprodc_beg          (begc:endc)) ; this%totprodc_beg (:)=spval
    allocate(this%seedc_beg             (begc:endc)) ; this%seedc_beg(:) = spval
    allocate(this%ctrunc_beg            (begc:endc)) ; this%ctrunc_beg(:)=spval
    allocate(this%totpftc_beg           (begc:endc)) ; this%totpftc_beg(:)=spval
    allocate(this%cropseedc_deficit_beg (begc:endc)) ; this%cropseedc_deficit_beg(:)=spval

    do c = begc, endc
       this%totpltc_beg(c)=0._r8
       this%totdsomc_beg(c)=0._r8
       this%totprodc_beg (:)=0._r8
       this%seedc_beg(:) = 0._r8
       this%ctrunc_beg(:)=0._r8
       this%totpftc_beg(:)=0._r8
       this%cropseedc_deficit_beg(:)=0._r8       
    enddo

                  
  end subroutine col_dbgs_init
!------------------------------------------------------------------------
  subroutine col_dbgs_clean(this)
    !
    ! !ARGUMENTS:
    class(column_debug_state) :: this
    
  end subroutine col_dbgs_clean
!------------------------------------------------------------------------
  subroutine EnterCMassDebug(bounds, num_soilc, filter_soilc, col_cs, veg_cs)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning carbon balance for mass
    ! conservation debug.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_carbon_state) , intent(in) :: col_cs
    type(vegetation_carbon_state), intent(in) :: veg_cs
    ! !LOCAL VARIABLES:
    integer :: c, p, pi  ! indices
    integer :: fc        ! lake filter indices
    
    !-----------------------------------------------------------------------
    
    associate(                                &
      totabgc      => col_cs%totabgc        , &  !total vegetation C
      cwdc         => col_cs%cwdc           , &  !coarse woody debris C
      totlitc      => col_cs%totlitc        , &  !total non cwd litter C
      totsomc      => col_cs%totsomc        , &  !total soil organic C
      totdsomc_beg => col_dbgs%totdsomc_beg , &
      totpltc_beg  => col_dbgs%totpltc_beg    &
    )
    if(.not.use_debug)return
      do fc = 1,num_soilc
         c = filter_soilc(fc)      
         totdsomc_beg(c) = cwdc(c)     + &
            totlitc(c)  + &
            totsomc(c)             
         totpltc_beg(c) = totabgc(c)

        col_dbgs%totprodc_beg (c) = col_cs%totprodc(c)
        col_dbgs%seedc_beg(c)     = col_cs%seedc(c)
        col_dbgs%ctrunc_beg(c)    = col_cs%ctrunc(c)
        col_dbgs%totpftc_beg(c)   = col_cs%totpftc(c)
        col_dbgs%cropseedc_deficit_beg(c) = col_cs%cropseedc_deficit(c)

        do pi = 1,max_patch_per_col
          if (pi <= col_pp%npfts(c)) then
            p = col_pp%pfti(c) + pi - 1
            if (veg_pp%active(p)) then    
              pft_dbgs%totpftc_beg(p) = veg_cs%totpftc(p)
              pft_dbgs%totvegc_beg(p) = veg_cs%totvegc(p)
              pft_dbgs%xsmrpool_beg(p)= veg_cs%xsmrpool(p)
              pft_dbgs%ctrunc_beg(p) = veg_cs%ctrunc(p)
              pft_dbgs%dispvegc_beg(p) =veg_cs%dispvegc(p)
              pft_dbgs%storvegc_beg(p) = veg_cs%storvegc(p)

              pft_dbgs%leafc_beg(p)   = veg_cs%leafc(p)
              pft_dbgs%frootc_beg(p)   = veg_cs%frootc(p)
              pft_dbgs%livestemc_beg(p)  =veg_cs%livestemc(p)
              pft_dbgs%deadstemc_beg(p)  = veg_cs%deadstemc(p)
              pft_dbgs%livecrootc_beg(p) = veg_cs%livecrootc(p)
              pft_dbgs%deadcrootc_beg(p) = veg_cs%deadcrootc(p)
              pft_dbgs%grainc_beg(p) = veg_cs%grainc(p)
            endif
          endif
        enddo

      enddo
  end associate
end subroutine EnterCMassDebug

    !-----------------------------------------------------------------------
  subroutine DebugCMassBal(c, col_cs, col_cf,veg_cs,veg_cf)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, compare the end-of-step decomposable soil
    ! carbon (cwdc + totlitc + totsomc) and total vegetation carbon against
    ! the beginning-of-step values stored by EnterCMassDebug in col_dbgs,
    ! for the failing column c only, for mass conservation debugging.
    !
    ! !ARGUMENTS:
    integer                   , intent(in) :: c               ! failing soil column index
    type(column_carbon_state) , intent(in) :: col_cs
    type(column_carbon_flux)  , intent(in) :: col_cf          ! available for flux debugging
    type(vegetation_carbon_state), intent(in) :: veg_cs    
    type(vegetation_carbon_flux), intent(in) :: veg_cf
    !
    ! !LOCAL VARIABLES:
    integer  :: nstep, p, pi
    real(r8) :: dsomc_end, dsomc_beg, dsomc_dif
    real(r8) :: totpltc_dif
    real(r8) :: dt,prodc_loss
    !-----------------------------------------------------------------------

    if (.not. use_debug) return

    nstep = get_nstep()
    dt = real( get_step_size(), r8 )
    dsomc_end   = col_cs%cwdc(c) + col_cs%totlitc(c) + col_cs%totsomc(c)
    dsomc_beg   = col_dbgs%totdsomc_beg(c)
    dsomc_dif   = dsomc_end - dsomc_beg
    totpltc_dif = col_cs%totabgc(c) - col_dbgs%totpltc_beg(c)

    prodc_loss= col_cf%prod1c_loss(c)*dt                 &
            + col_cf%hrv_deadstemc_to_prod10c(c)*dt      & ! from harvest
            - col_cf%prod10c_loss(c)*dt                  & ! from decomposition
            + col_cf%hrv_deadstemc_to_prod100c(c)*dt     & ! from harvest
            - col_cf%prod100c_loss(c)*dt                   ! from decomposition            

    write(iulog,*)'DebugCMassBal: c=',c,' step=',nstep
    write(iulog,*)'  dsomc  (end, beg, delta) = ',dsomc_end, dsomc_beg, dsomc_dif
    write(iulog,*)'  totpltc(end, beg, delta) = ',col_cs%totabgc(c), col_dbgs%totpltc_beg(c), totpltc_dif
    write(iulog,*)'somc err =',dsomc_dif+totpltc_dif,dsomc_dif-col_cf%plant_to_litter_cflux(c)*dt+col_cf%hr(c)*dt
    write(iulog,*)'pft err=',totpltc_dif+col_cf%plant_to_litter_cflux(c)*dt+col_cf%ar(c)*dt+col_cf%product_closs(c)*dt
    write(iulog,*)'totc (beg,end)=',dsomc_beg+col_dbgs%totpltc_beg(c),&
      dsomc_end+col_cs%totabgc(c)
    write(iulog,*)'totc11(beg,end)=',dsomc_beg+col_dbgs%totpltc_beg(c)+col_dbgs%cropseedc_deficit_beg(c),&
      dsomc_end+col_cs%totabgc(c)+col_cs%cropseedc_deficit(c)      
    write(iulog,*)'veg2litr  =', col_cf%plant_to_litter_cflux(c)*dt
    write(iulog,*)'veg2cwd   =', col_cf%plant_to_cwd_cflux(c)*dt
    write(iulog,*)'prodc (in, out,netloss)=', col_cf%totprodc_in(c)*dt,col_cf%totprodc_out(c)*dt,(col_cf%totprodc_in(c)-col_cf%totprodc_out(c))*dt
    write(iulog,*)'leach     =', col_cf%som_c_leached(c)*dt
    write(iulog,*)'totprodc (end, beg, delta) =', col_cs%totprodc(c), col_dbgs%totprodc_beg (c),&
     col_cs%totprodc(c)-col_dbgs%totprodc_beg (c)-(col_cf%totprodc_in(c)-col_cf%totprodc_out(c))*dt
    write(iulog,*)'seedc (end, beg) =', col_cs%seedc(c), col_dbgs%seedc_beg(c)
    write(iulog,*)'ctrunc (end,beg) = ',col_cs%ctrunc(c), col_dbgs%ctrunc_beg(c)
    write(iulog,*)'totpftc (end, beg)=',col_cs%totpftc(c),col_dbgs%totpftc_beg(c)
    write(iulog,*)'cropseedc_deficit (end,beg) =', col_cs%cropseedc_deficit(c),col_dbgs%cropseedc_deficit_beg(c)
    write(iulog,*)'crop_seedc_to_leaf=',col_cf%crop_seedc_to_leaf(c)*dt
    write(iulog,*)'totpftc (in,out,net)=',col_cf%totpftc_in(c)*dt,col_cf%totpftc_out(c)*dt,&
      (col_cf%totpftc_in(c)-col_cf%totpftc_out(c))*dt
    write(iulog,*)'dtotpftc =',col_cs%totpftc(c)-col_dbgs%totpftc_beg(c)-(col_cf%totpftc_in(c)-col_cf%totpftc_out(c))*dt &
      + col_cf%crop_seedc_to_leaf(c)*dt 
    do pi = 1,max_patch_per_col
      if (pi <= col_pp%npfts(c)) then
         p = col_pp%pfti(c) + pi - 1
         if (veg_pp%active(p)) then    
            write(iulog,*)'toltpftc (end,beg,delta)=',p,veg_cs%totpftc(p),pft_dbgs%totpftc_beg(p),veg_cs%totpftc(p)-pft_dbgs%totpftc_beg(p)

            write(iulog,*)'totvegc (end,beg,delta)=',p,pft_dbgs%totvegc_beg(p),veg_cs%totvegc(p),pft_dbgs%totvegc_beg(p)-veg_cs%totvegc(p)
            write(iulog,*)'xsmrpool (end,beg,delta)=',p,  pft_dbgs%xsmrpool_beg(p),veg_cs%xsmrpool(p),  pft_dbgs%xsmrpool_beg(p)-veg_cs%xsmrpool(p)
            write(iulog,*)'ctrunc (end,beg,delta)=',p, pft_dbgs%ctrunc_beg(p),veg_cs%ctrunc(p),pft_dbgs%ctrunc_beg(p)-veg_cs%ctrunc(p)
            write(iulog,*)'dispvegc (end,beg,detal)=',p,pft_dbgs%dispvegc_beg(p),veg_cs%dispvegc(p),pft_dbgs%dispvegc_beg(p)-veg_cs%dispvegc(p)
            write(iulog,*)'storvegc (end,beg,delta)=',p,pft_dbgs%storvegc_beg(p),veg_cs%storvegc(p),pft_dbgs%storvegc_beg(p)-veg_cs%storvegc(p)
            write(iulog,*)'leafc (end,beg,delta)=',p,pft_dbgs%leafc_beg(p), veg_cs%leafc(p), pft_dbgs%leafc_beg(p) - veg_cs%leafc(p)
            write(iulog,*)'frootc (end,beg,delta)=',p,pft_dbgs%frootc_beg(p),veg_cs%frootc(p),pft_dbgs%frootc_beg(p)-veg_cs%frootc(p)
            write(iulog,*)'livestemc (end,beg,delta)=',p,pft_dbgs%livestemc_beg(p),veg_cs%livestemc(p),pft_dbgs%livestemc_beg(p)-veg_cs%livestemc(p)
            write(iulog,*)'deadstemc (end,beg,delta)=',p,pft_dbgs%deadstemc_beg(p), veg_cs%deadstemc(p),pft_dbgs%deadstemc_beg(p)- veg_cs%deadstemc(p)
            write(iulog,*)'livecrootc (end,beg,delta)=',p,pft_dbgs%livecrootc_beg(p), veg_cs%livecrootc(p),pft_dbgs%livecrootc_beg(p) - veg_cs%livecrootc(p)
            write(iulog,*)'deadcrootc (end,beg,delta)=',p,pft_dbgs%deadcrootc_beg(p), veg_cs%deadcrootc(p),pft_dbgs%deadcrootc_beg(p) - veg_cs%deadcrootc(p)
            write(iulog,*)'grainc (end,beg,delta)=',p,pft_dbgs%grainc_beg(p), veg_cs%grainc(p),pft_dbgs%grainc_beg(p) - veg_cs%grainc(p)
            write(iulog,*)'frootc_loss =',p,veg_cf%frootc_loss(p)*dt 
            write(iulog,*)'m_frootc_to_litter=',p,veg_cf%m_frootc_to_litter(p)*dt
            write(iulog,*)'m_frootc_to_fire=',p,veg_cf%m_frootc_to_fire(p)*dt
            write(iulog,*)'m_frootc_to_litter_fire=',p,veg_cf%m_frootc_to_litter_fire(p)*dt
            write(iulog,*)'frootc_to_litter=',p,veg_cf%frootc_to_litter(p)*dt
            write(iulog,*)'hrv_frootc_to_litter=',p,veg_cf%hrv_frootc_to_litter(p)*dt      
            write(iulog,*)'dt=',dt
          endif
        endif
     enddo

  end subroutine DebugCMassBal

    !-----------------------------------------------------------------------
  subroutine DebugWritep(p,loc,varname,val,laction)
  implicit none
  integer, intent(in) :: p
  character(len=*), intent(in) :: loc
  character(len=*), intent(in) :: varname
  real(r8), dimension(:), intent(in) :: val
  logical, optional, intent(in) :: laction

  if(.not.use_debug)return
  if(present(laction))then
    if(.not.laction)return
  endif

  
  write(iulog,*)trim(loc),'p '//trim(varname)//'=',p,val(p)

  end subroutine DebugWritep

end module DebugToolsMod
