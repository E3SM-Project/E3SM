
      module mo_imp_sol

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods,    only : clscnt4, gas_pcnst, clsmap
      use ppgrid,       only : pcols

      implicit none

!-----------------------------------------------------------------------      
!   	    newton-Raphson iteration limits
!-----------------------------------------------------------------------      
      integer, parameter  :: itermax   = 11
      integer, parameter  :: cut_limit = 5
      integer, parameter  :: vec_len   = pcols
      real(r8), parameter :: sol_min   = 1.e-20_r8
      real(r8), parameter :: small     = 1.e-40_r8

      save

      real(r8) :: epsilon(clscnt4)
      logical  :: factor(itermax)

      private
      public :: imp_slv_inti
      public :: imp_sol

      contains

      subroutine imp_slv_inti
!-----------------------------------------------------------------------      
!	... initialize the implict solver
!-----------------------------------------------------------------------      

      use m_spc_id

      implicit none

!-----------------------------------------------------------------------      
!	... local variables
!-----------------------------------------------------------------------      
      integer  :: m
      real(r8) :: eps(gas_pcnst)

      factor(:) = .true.
      eps(:)    = .001_r8
      eps((/id_o3,id_no,id_no2,id_no3,id_hno3,id_ho2no2,id_n2o5,id_oh,id_ho2/)) = .0001_r8
      do m = 1,clscnt4
         epsilon(m) = eps(clsmap(m,4))
      end do

      end subroutine imp_slv_inti

      subroutine imp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, &
                          ncol, lchnk, chnkpnts )
!-----------------------------------------------------------------------
!      	... imp_sol advances the volumetric mixing ratio
!           forward one time step via the fully implicit euler scheme.
!           this source is meant for vector architectures such as the
!           nec sx6 and cray x1
!-----------------------------------------------------------------------

      use chem_mods,     only : rxntot, extcnt, nzcnt, diag_map, permute, cls_rxt_cnt
      use mo_tracname,   only : solsym
      use ppgrid,        only : pver
      use pmgrid,        only : iam
      use mo_lin_matrix, only : linmat
      use mo_nln_matrix, only : nlnmat
      use mo_lu_factor,  only : lu_fac
      use mo_lu_solve,   only : lu_slv
      use mo_prod_loss,  only : imp_prod_loss
      use mo_indprd,     only : indprd
      use time_manager,  only : get_nstep

      implicit none

!-----------------------------------------------------------------------
!     	... dummy args
!-----------------------------------------------------------------------
      integer, intent(in)     :: ncol                                        ! columns in chunck
      integer, intent(in)     :: lchnk                                       ! chunk id
      integer, intent(in)     :: chnkpnts                                    ! total spatial points in chunk; ncol*pver
      real(r8), intent(in)    :: delt                                        ! time step (s)
      real(r8), intent(in)    :: reaction_rates(chnkpnts,max(1,rxntot))      ! rxt rates (1/cm^3/s)
      real(r8), intent(in)    :: extfrc(chnkpnts,max(1,extcnt))              ! external in-situ forcing (1/cm^3/s)
      real(r8), intent(in)    :: het_rates(chnkpnts,max(1,gas_pcnst))           ! washout rates (1/s)
      real(r8), intent(inout) :: base_sol(chnkpnts,gas_pcnst)                ! species mixing ratios (vmr)

!-----------------------------------------------------------------------
!     	... local variables
!-----------------------------------------------------------------------
      integer  ::  nr_iter
      integer  ::  ofl
      integer  ::  ofu
      integer  ::  bndx                                                      ! base index
      integer  ::  cndx                                                      ! class index
      integer  ::  pndx                                                      ! permuted class index
      integer  ::  m
      integer  ::  fail_cnt
      integer  ::  cut_cnt
      integer  ::  stp_con_cnt
      integer  ::  nstep
      real(r8) ::  interval_done
      real(r8) ::  dt
      real(r8) ::  dti
      real(r8) ::  max_delta(max(1,clscnt4))
      real(r8) ::  sys_jac(chnkpnts,max(1,nzcnt))
      real(r8) ::  lin_jac(chnkpnts,max(1,nzcnt))
      real(r8) ::  solution(chnkpnts,max(1,clscnt4))
      real(r8) ::  forcing(chnkpnts,max(1,clscnt4))
      real(r8) ::  iter_invariant(chnkpnts,max(1,clscnt4))
      real(r8) ::  prod(chnkpnts,max(1,clscnt4))
      real(r8) ::  loss(chnkpnts,max(1,clscnt4))
      real(r8) ::  ind_prd(chnkpnts,max(1,clscnt4))
      real(r8) ::  sbase_sol(chnkpnts,gas_pcnst)
      real(r8) ::  wrk(chnkpnts)
      logical  ::  convergence
      logical  ::  spc_conv(chnkpnts,max(1,clscnt4))
      logical  ::  cls_conv(chnkpnts)
      logical  ::  converged(max(1,clscnt4))

!-----------------------------------------------------------------------      
!        ... class independent forcing
!-----------------------------------------------------------------------      
      if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
         call indprd( 4, ind_prd, base_sol, extfrc, reaction_rates )
      else
         do m = 1,clscnt4
            ind_prd(:,m) = 0._r8
         end do
      end if

      ofl = 1
chnkpnts_loop : &
      do 
         ofu = min( chnkpnts,ofl + vec_len - 1 )
         do m = 1,gas_pcnst
            sbase_sol(ofl:ofu,m) = base_sol(ofl:ofu,m)
         end do
!-----------------------------------------------------------------------      
!        ... time step loop
!-----------------------------------------------------------------------      
         dt            = delt
         cut_cnt       = 0
         fail_cnt      = 0
         stp_con_cnt   = 0
         interval_done = 0._r8
time_step_loop : &
         do
            dti = 1._r8 / dt
!-----------------------------------------------------------------------      
!        ... transfer from base to class array
!-----------------------------------------------------------------------      
            do cndx = 1,clscnt4
               bndx = clsmap(cndx,4)
               pndx = permute(cndx,4)
               solution(ofl:ofu,pndx) = base_sol(ofl:ofu,bndx)
            end do
!-----------------------------------------------------------------------      
!        ... set the iteration invariant part of the function f(y)
!-----------------------------------------------------------------------      
            if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
               do m = 1,clscnt4
                  iter_invariant(ofl:ofu,m) = dti * solution(ofl:ofu,m) + ind_prd(ofl:ofu,m)
               end do
            else
               do m = 1,clscnt4
                  iter_invariant(ofl:ofu,m) = dti * solution(ofl:ofu,m)
               end do
            end if
!-----------------------------------------------------------------------      
!        ... the linear component
!-----------------------------------------------------------------------      
            if( cls_rxt_cnt(2,4) > 0 ) then
               call linmat( ofl, ofu, lin_jac, base_sol, reaction_rates, het_rates )
            end if
!=======================================================================
!        the newton-raphson iteration for f(y) = 0
!=======================================================================
            cls_conv(ofl:ofu) = .false.
iter_loop : do nr_iter = 1,itermax
!-----------------------------------------------------------------------      
!        ... the non-linear component
!-----------------------------------------------------------------------      
               if( factor(nr_iter) ) then
                  call nlnmat( ofl, ofu, sys_jac, base_sol, reaction_rates, &
                               lin_jac, dti )
!-----------------------------------------------------------------------      
!         ... factor the "system" matrix
!-----------------------------------------------------------------------      
                  call lu_fac( ofl, ofu, sys_jac )
               end if      
!-----------------------------------------------------------------------      
!   	... form f(y)
!-----------------------------------------------------------------------      
               call imp_prod_loss( ofl, ofu, prod, loss, base_sol, &
                                   reaction_rates, het_rates )
               do m = 1,clscnt4
                  forcing(ofl:ofu,m) = solution(ofl:ofu,m)*dti &
                                     - (iter_invariant(ofl:ofu,m) + prod(ofl:ofu,m) - loss(ofl:ofu,m))
               end do
!-----------------------------------------------------------------------      
!         ... solve for the mixing ratio at t(n+1)
!-----------------------------------------------------------------------      
               call lu_slv( ofl, ofu, sys_jac, forcing )
               do m = 1,clscnt4
                  where( .not. cls_conv(ofl:ofu) )
                     solution(ofl:ofu,m) = solution(ofl:ofu,m) + forcing(ofl:ofu,m)
                  elsewhere
                     forcing(ofl:ofu,m)  = 0._r8
                  endwhere
               end do
!-----------------------------------------------------------------------      
!    	... convergence measures and test
!-----------------------------------------------------------------------      
conv_chk :     if( nr_iter > 1 ) then
!-----------------------------------------------------------------------      
!    	... check for convergence
!-----------------------------------------------------------------------      
                  do cndx = 1,clscnt4
                     pndx = permute(cndx,4)
                     bndx = clsmap(cndx,4)
                     where( abs( solution(ofl:ofu,pndx) ) > sol_min )
                        wrk(ofl:ofu) = abs( forcing(ofl:ofu,pndx)/solution(ofl:ofu,pndx) )
                     elsewhere
                        wrk(ofl:ofu) = 0._r8
                     endwhere
                     max_delta(cndx)        = maxval( wrk(ofl:ofu) )
                     solution(ofl:ofu,pndx) = max( 0._r8,solution(ofl:ofu,pndx) )
                     base_sol(ofl:ofu,bndx) = solution(ofl:ofu,pndx)
                     where( abs( forcing(ofl:ofu,pndx) ) > small )
                        spc_conv(ofl:ofu,cndx) = abs(forcing(ofl:ofu,pndx)) <= epsilon(cndx)*abs(solution(ofl:ofu,pndx))
                     elsewhere
                        spc_conv(ofl:ofu,cndx) =  .true.
                     endwhere
                     converged(cndx) = all( spc_conv(ofl:ofu,cndx) )
                  end do
                  convergence = all( converged(:) )
                  if( convergence ) then
                     exit iter_loop
                  end if
                  do m = ofl,ofu
                     if( .not. cls_conv(m) ) then
                        cls_conv(m) = all( spc_conv(m,:) )
                     end if
                  end do
               else conv_chk
!-----------------------------------------------------------------------      
!   	... limit iterate
!-----------------------------------------------------------------------      
                  do m = 1,clscnt4
                     solution(ofl:ofu,m) = max( 0._r8,solution(ofl:ofu,m) )
                  end do
!-----------------------------------------------------------------------      
!   	... transfer latest solution back to base array
!-----------------------------------------------------------------------      
                  do cndx = 1,clscnt4
                     pndx = permute(cndx,4)
                     bndx = clsmap(cndx,4)
                     base_sol(ofl:ofu,bndx) = solution(ofl:ofu,pndx)
                  end do
               end if conv_chk
            end do iter_loop

!-----------------------------------------------------------------------      
!    	... check for newton-raphson convergence
!-----------------------------------------------------------------------      
non_conv :  if( .not. convergence ) then
!-----------------------------------------------------------------------      
!   	... non-convergence
!-----------------------------------------------------------------------      
               fail_cnt = fail_cnt + 1
               nstep    = get_nstep()
               write(*,'('' imp_sol: time step '',1p,g15.7,'' failed to converge @ (lchnk,lev,nstep) = '',3i6)') &
                              dt,lchnk,lev,nstep
               stp_con_cnt = 0
step_reduction : &
               if( cut_cnt < cut_limit ) then
                  cut_cnt = cut_cnt + 1
                  if( cut_cnt < cut_limit ) then
                     dt = .5_r8 * dt
                  else
                     dt = .1_r8 * dt
                  end if
                  do m = 1,gas_pcnst
                     base_sol(ofl:ofu,m) = sbase_sol(ofl:ofu,m)
                  end do
                  cycle time_step_loop
               else step_reduction
                  write(*,'('' imp_sol: step failed to converge @ (lchnk,lev,nstep,dt,time) = '',3i6,1p,2g15.7)') &
                        lchnk,lev,nstep,dt,interval_done+dt
                  do m = 1,clscnt4
                     if( .not. converged(m) ) then
                        write(*,'(1x,a8,1x,1pe10.3)') solsym(clsmap(m,4)), max_delta(m)
                     end if
                  end do
               end if step_reduction
            end if non_conv
!-----------------------------------------------------------------------      
!   	... check for interval done
!-----------------------------------------------------------------------      
            interval_done = interval_done + dt
time_step_done : &
            if( abs( delt - interval_done ) <= .0001_r8 ) then
               if( fail_cnt > 0 ) then
                  write(*,*) 'imp_sol : @ (lchnk,lev) = ',lchnk,lev,' failed ',fail_cnt,' times'
               end if
               exit time_step_loop
            else time_step_done
!-----------------------------------------------------------------------      
!   	... transfer latest solution back to base array
!-----------------------------------------------------------------------      
               if( convergence ) then
                  stp_con_cnt = stp_con_cnt + 1
               end if
               do m = 1,gas_pcnst
                  sbase_sol(ofl:ofu,m) = base_sol(ofl:ofu,m)
               end do
               if( stp_con_cnt >= 2 ) then
                  dt = 2._r8*dt
                  stp_con_cnt = 0
               end if
               dt = min( dt,delt-interval_done )
            end if time_step_done
         end do time_step_loop
         ofl = ofu + 1
         if( ofl > chnkpnts ) then
            exit chnkpnts_loop
         end if
      end do chnkpnts_loop

      end subroutine imp_sol

      end module mo_imp_sol
