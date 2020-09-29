
      module mo_imp_sol

      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods,    only : clscnt4, gas_pcnst, clsmap

      implicit none

!-----------------------------------------------------------------------      
!   	    Newton-Raphson iteration limits
!-----------------------------------------------------------------------      
      integer, parameter :: itermax      = 11
      integer, parameter :: cut_limit    = 5

      save

      real(r8) :: small
      real(r8) :: epsilon(clscnt4)
      logical  :: factor(itermax)

      private
      public :: imp_slv_inti, imp_sol

      contains

      subroutine imp_slv_inti
!-----------------------------------------------------------------------      
!	... Initialize the implict solver
!-----------------------------------------------------------------------      

      use m_spc_id

      implicit none

!-----------------------------------------------------------------------      
!	... Local variables
!-----------------------------------------------------------------------      
      integer  :: m
      real(r8) :: eps(gas_pcnst)

      small     = 1.e6_r8 * tiny( small )
      factor(:) = .true.
      eps(:)    = .001_r8
      eps((/id_o3,id_no,id_no2,id_no3,id_hno3,id_ho2no2,id_n2o5,id_oh,id_ho2/)) = .0001_r8
      do m = 1,clscnt4
         epsilon(m) = eps(clsmap(m,4))
      end do

      end subroutine imp_slv_inti

      subroutine imp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, &
			  ncol, lchnk )
!-----------------------------------------------------------------------
!      	... Imp_sol advances the volumetric mixing ratio
!           forward one time step via the fully implicit
!           Euler scheme
!-----------------------------------------------------------------------

      use chem_mods,     only : rxntot, extcnt, nzcnt, clsze, diag_map, permute, cls_rxt_cnt
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
!     	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)   ::   ncol                                    ! columns in chunck
      integer, intent(in)   ::   lchnk                                   ! chunk id
      real(r8), intent(in)  ::   delt                                    ! time step (s)
      real(r8), intent(in)  ::   reaction_rates(ncol,pver,rxntot), &     ! rxt rates (1/cm^3/s)
                                 extfrc(ncol,pver,extcnt), &             ! external in-situ forcing (1/cm^3/s)
                                 het_rates(ncol,pver,max(1,gas_pcnst))      ! washout rates (1/s)
      real(r8), intent(inout) :: base_sol(ncol,pver,gas_pcnst)           ! species mixing ratios (vmr)

!-----------------------------------------------------------------------
!     	... Local variables
!-----------------------------------------------------------------------
      integer ::   nr_iter, &
                   lev, &
                   ofl, ofu, &
                   i, isec, isecu, &
                   j, &
                   k, l, &
                   m, cols
      integer ::   nstep
      integer ::   stp_con_cnt, cut_cnt, fail_cnt
      real(r8) :: interval_done, dt, dti
      real(r8) :: max_delta(clscnt4)
      real(r8), dimension(clsze,nzcnt) :: &
                   sys_jac, &
                   lin_jac
      real(r8), dimension(clsze,clscnt4) :: &
                   solution, &
                   forcing, &
                   iter_invariant, &
                   prod, &
                   loss
      real(r8) :: lrxt(clsze,rxntot)
      real(r8) :: lhet(clsze,max(1,gas_pcnst))
      real(r8) :: lsol(clsze,gas_pcnst)
      real(r8), dimension(clsze) :: &
                   wrk
      real(r8), dimension(ncol,pver,clscnt4) :: &
                   ind_prd
      logical ::   convergence
      logical ::   iter_conv(clsze)
      logical ::   converged(clscnt4)
      logical ::   do_diag

!-----------------------------------------------------------------------      
!        ... Class independent forcing
!-----------------------------------------------------------------------      
#ifdef DEBUG
      call t_startf('indprd')
#endif
      if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
         call indprd( 4, ind_prd, clscnt4, base_sol, extfrc, &
                      reaction_rates, ncol )
      else
         do m = 1,max(1,clscnt4)
            ind_prd(:,:,m) = 0._r8
         end do
      end if
#ifdef DEBUG
      call t_stopf('indprd')
#endif
      isecu = (ncol - 1)/clsze + 1
Level_loop : &
      do lev = 1,pver
Column_loop : &
	 do isec = 1,isecu
	    ofl = (isec - 1)*clsze + 1
	    ofu = min( ncol,ofl + clsze - 1 )
	    cols = ofu - ofl + 1
!-----------------------------------------------------------------------      
!        ... Transfer from base to local work arrays
!-----------------------------------------------------------------------      
            do m = 1,rxntot
	       lrxt(:cols,m) = reaction_rates(ofl:ofu,lev,m) 
	    end do
            do m = 1,gas_pcnst
	       lhet(:cols,m) = het_rates(ofl:ofu,lev,m) 
	    end do
!           do_diag = lev == 5 .and. isec == 2 .and. lchnk == 290
            do_diag = .false.
            if( do_diag ) then
	      write(*,*) ' '
	      write(*,*) 'imp_sol: lchnk,lev,isec,cols = ',lchnk,lev,isec,cols
	      write(*,*) ' '
	      write(*,*) 'imp_sol: lrxt'
	      write(*,'(1p,4g20.10)') lrxt(1,:)
	      write(*,*) 'imp_sol: lhet'
	      write(*,'(1p,4g20.10)') lhet(1,:)
	    end if
!-----------------------------------------------------------------------      
!        ... Time step loop
!-----------------------------------------------------------------------      
	    dt            = delt
	    cut_cnt       = 0
	    stp_con_cnt   = 0
	    fail_cnt      = 0
	    interval_done = 0._r8
Time_step_loop : &
	    do
	       dti = 1._r8 / dt
!-----------------------------------------------------------------------      
!        ... Transfer from base to local work arrays
!-----------------------------------------------------------------------      
               do m = 1,gas_pcnst
	          lsol(:cols,m) = base_sol(ofl:ofu,lev,m) 
	       end do
!-----------------------------------------------------------------------      
!        ... Transfer from base to class array
!-----------------------------------------------------------------------      
               do k = 1,clscnt4
                  j = clsmap(k,4)
                  m = permute(k,4)
                  solution(:cols,m) = lsol(:cols,j)
               end do
               if( do_diag ) then
	         write(*,*) ' '
	         write(*,*) 'imp_sol: solution'
	         write(*,'(1p,4g20.10)') solution(1,:)
	       end if
!-----------------------------------------------------------------------      
!        ... Set the iteration invariant part of the function F(y)
!        ... If there is "independent" production put it in the forcing
!-----------------------------------------------------------------------      
               if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
                  do m = 1,clscnt4
                     iter_invariant(:cols,m) = dti * solution(:cols,m) + ind_prd(ofl:ofu,lev,m)
                  end do
               else
                  do m = 1,clscnt4
                     iter_invariant(:cols,m) = dti * solution(:cols,m)
                  end do
               end if
!-----------------------------------------------------------------------      
!        ... The linear component
!-----------------------------------------------------------------------      
#ifdef DEBUG
	       call t_startf('lin_mat')
#endif
	       if( cls_rxt_cnt(2,4) > 0 ) then
                  call linmat( lin_jac, lsol, lrxt, lhet, cols )
	       end if
#ifdef DEBUG
	       call t_stopf('lin_mat')
#endif

!=======================================================================
!        The Newton-Raphson iteration for F(y) = 0
!=======================================================================
Iteration_loop : &
               do nr_iter = 1,itermax
!-----------------------------------------------------------------------      
!        ... The non-linear component
!-----------------------------------------------------------------------      
                  if( factor(nr_iter) ) then
#ifdef DEBUG
		     call t_startf('nln_mat')
#endif
                     call nlnmat( sys_jac, lsol, lrxt, lin_jac, dti, cols )
#ifdef DEBUG
		     call t_stopf('nln_mat')
		     call t_startf('lu_fac')
#endif
!-----------------------------------------------------------------------      
!         ... Factor the "system" matrix
!-----------------------------------------------------------------------      
#ifdef DEBUG_LU
		     if( do_diag ) then
			write(*,*) 'imp_sol: before lu_fac - lchnk,lev,isec,iter,dt = ',lchnk,lev,isec,nr_iter,dt
			write(*,*) 'imp_sol: lu(110,20,21,109)'
			write(*,'(1p,4g20.10)') sys_jac(:cols,110)
			write(*,'(1p,4g20.10)') sys_jac(:cols,20)
			write(*,'(1p,4g20.10)') sys_jac(:cols,21)
			write(*,'(1p,4g20.10)') sys_jac(:cols,109)
			write(*,*) 'imp_sol: maxval sys_jac(3,:) = ',maxval( sys_jac(3,:) )
			write(*,*) 'imp_sol: maxval sys_jac(4,:) = ',maxval( sys_jac(4,:) )
		     end if
#endif
	             call lu_fac( sys_jac, cols )
#ifdef DEBUG_LU
		     if( do_diag ) then
			write(*,*) 'imp_sol: after lu_fac - lchnk,lev,isec,iter,dt = ',lchnk,lev,isec,nr_iter,dt
			write(*,*) 'imp_sol: lu(110,20,21,109)'
			write(*,'(1p,4g20.10)') sys_jac(:cols,110)
			write(*,'(1p,4g20.10)') sys_jac(:cols,20)
			write(*,'(1p,4g20.10)') sys_jac(:cols,21)
			write(*,'(1p,4g20.10)') sys_jac(:cols,109)
			write(*,*) 'imp_sol: maxval sys_jac(3,:) = ',maxval( sys_jac(3,:) )
			write(*,*) 'imp_sol: maxval sys_jac(4,:) = ',maxval( sys_jac(4,:) )
		     end if
#endif
#ifdef DEBUG
		     call t_stopf('lu_fac')
#endif
                  end if      
!-----------------------------------------------------------------------      
!   	... Form F(y)
!-----------------------------------------------------------------------      
#ifdef DEBUG
		  call t_startf('frcing')
#endif
                  call imp_prod_loss( prod, loss, lsol, lrxt, lhet, cols )
	          do m = 1,clscnt4
                     forcing(:cols,m) = solution(:cols,m)*dti - (iter_invariant(:cols,m) + prod(:cols,m) - loss(:cols,m))
	          end do

                  if( do_diag ) then
	            write(*,*) ' '
	            write(*,*) 'imp_sol: frcing @ iter,dt = ',nr_iter,dt
	            write(*,'(1p,4g20.10)') forcing(1,:)
	            write(*,*) ' '
	            write(*,*) 'imp_sol: iter_invariant @ iter,dt = ',nr_iter,dt
	            write(*,'(1p,4g20.10)') iter_invariant(1,:)
	            write(*,*) ' '
	            write(*,*) 'imp_sol: prod @ iter,dt = ',nr_iter,dt
	            write(*,'(1p,4g20.10)') prod(1,:)
	            write(*,*) ' '
	            write(*,*) 'imp_sol: loss @ iter,dt = ',nr_iter,dt
	            write(*,'(1p,4g20.10)') loss(1,:)
	          end if
#ifdef DEBUG
		  call t_stopf('frcing')
		  call t_startf('lu_slv')
#endif
!-----------------------------------------------------------------------      
!         ... Solve for the mixing ratio at t(n+1)
!-----------------------------------------------------------------------      
	          call lu_slv( sys_jac, forcing, cols )
                  if( do_diag ) then
	            write(*,*) ' '
	            write(*,*) 'imp_sol: frcing @ iter,dt = ',nr_iter,dt
	            write(*,'(1p,4g20.10)') forcing(1,:)
	          end if
	          do m = 1,clscnt4
                     solution(:cols,m) = solution(:cols,m) + forcing(:cols,m)
	          end do
                  if( do_diag ) then
	            write(*,*) ' '
	            write(*,*) 'imp_sol: solution @ iter,dt = ',nr_iter,dt
	            write(*,'(1p,4g20.10)') solution(1,:)
	          end if
#ifdef DEBUG
		  call t_stopf('lu_slv')
#endif
!-----------------------------------------------------------------------      
!    	... Convergence measures
!-----------------------------------------------------------------------      
#ifdef DEBUG
                  call t_startf('inner1_conv')
#endif
                  if( nr_iter > 1 ) then
	             do k = 1,clscnt4
		        m = permute(k,4)
			do i = 1,cols
		           if( abs(solution(i,m)) > 1.e-40_r8 ) then
		              wrk(i) = abs( forcing(i,m)/solution(i,m) )
			   else
		              wrk(i) = 0._r8
			   end if
			end do
!	        where( abs(solution(:cols,m)) > 1.e-40_r8 )
!	           wrk(:cols) = abs( forcing(:cols,m)/solution(:cols,m) )
!	        elsewhere
!	           wrk(:cols) = 0._r8
!	        endwhere
		        max_delta(k) = maxval( wrk(:cols) )
	             end do
	          end if
#ifdef DEBUG
                  call t_stopf('inner1_conv')
                  call t_startf('inner2_conv')
#endif
!-----------------------------------------------------------------------      
!   	... Limit iterate
!-----------------------------------------------------------------------      
	          do m = 1,clscnt4
!                    where( solution(:cols,m) < 0._r8 )
!                solution(:cols,m) = 0._r8
!             endwhere
                     do i = 1,cols
			if( solution(i,m) < 0. ) then
			   solution(i,m) = 0.
			end if
		     end do
	          end do
!-----------------------------------------------------------------------      
!   	... Transfer latest solution back to work array
!-----------------------------------------------------------------------      
                  do k = 1,clscnt4
                     j = clsmap(k,4)
                     m = permute(k,4)
		     do i = 1,cols
                        lsol(i,j) = solution(i,m)
		     end do
!                    lsol(:cols,j) = solution(:cols,m)
                  end do
#ifdef DEBUG
                  call t_stopf('inner2_conv')
                  call t_startf('inner3_conv')
#endif
!-----------------------------------------------------------------------      
!    	... Check for convergence
!-----------------------------------------------------------------------      
	          if( nr_iter > 1 ) then
	             do k = 1,clscnt4
		        m = permute(k,4)
			do i = 1,cols
			   if( abs( forcing(i,m) ) > small ) then
			      iter_conv(i) =  abs(forcing(i,m)) <= epsilon(k)*abs(solution(i,m))
			   else
			      iter_conv(i) =  .true.
			   end if
			end do
		        converged(k) =  all( iter_conv(:cols) )
	             end do
	             convergence = all( converged(:clscnt4) )
	             if( convergence ) then
#ifdef DEBUG
                        call t_stopf('inner3_conv')
#endif
	                exit Iteration_loop
	             end if
	          end if
#ifdef DEBUG
                  call t_stopf('inner3_conv')
#endif
               end do Iteration_loop

!-----------------------------------------------------------------------      
!    	... Check for Newton-Raphson convergence
!-----------------------------------------------------------------------      
#ifdef DEBUG
               call t_startf('outer_conv')
#endif
               if( .not. convergence ) then
!-----------------------------------------------------------------------      
!   	... Non-convergence
!-----------------------------------------------------------------------      
		  fail_cnt = fail_cnt + 1
		  nstep = get_nstep()
 	          write(*,'('' imp_sol: Time step '',1p,e21.13,'' failed to converge @ (lchnk,lev,isec,nstep) = '',4i6)') dt,lchnk,lev,isec,nstep
		  stp_con_cnt = 0
                  if( cut_cnt < cut_limit ) then
		     cut_cnt = cut_cnt + 1
		     if( cut_cnt < cut_limit ) then
		        dt = .5_r8 * dt
		     else
		        dt = .1_r8 * dt
		     end if
		     cycle Time_step_loop
		  else
	             write(*,'('' imp_sol: Failed to converge @ (lchnk,lev,isec,nstep,dt,time) = '',4i6,1p,2e21.13)') &
				       lchnk,lev,isec,nstep,dt,interval_done+dt
	             do m = 1,clscnt4
	                if( .not. converged(m) ) then
	                   write(*,'(1x,a8,1x,1pe10.3)') solsym(clsmap(m,4)), max_delta(m)
	                end if
	             end do
		  end if
	       end if
!-----------------------------------------------------------------------      
!   	... Check for interval done
!-----------------------------------------------------------------------      
	       interval_done = interval_done + dt
	       if( abs( delt - interval_done ) <= .0001_r8 ) then
		  if( fail_cnt > 0 ) then
	             write(*,*) 'imp_sol : @ (lchnk,lev,isec) = ',lchnk,lev,isec,' failed ',fail_cnt,' times'
		  end if
#ifdef DEBUG
                  call t_stopf('outer_conv')
#endif
		  exit Time_step_loop
	       else
!-----------------------------------------------------------------------      
!   	... Transfer latest solution back to base array
!-----------------------------------------------------------------------      
		  if( convergence ) then
		     stp_con_cnt = stp_con_cnt + 1
		  end if
                  do m = 1,gas_pcnst
                     base_sol(ofl:ofu,lev,m) = lsol(:cols,m)
                  end do
		  if( stp_con_cnt >= 2 ) then
		     dt = 2._r8*dt
		     stp_con_cnt = 0
		  end if
		  dt = min( dt,delt-interval_done )
!	  write(*,'('' imp_sol: New time step '',1p,e21.13)') dt
               end if
#ifdef DEBUG
               call t_stopf('outer_conv')
#endif
            end do Time_step_loop
!-----------------------------------------------------------------------      
!   	... Transfer latest solution back to base array
!-----------------------------------------------------------------------      
            do k = 1,clscnt4
               j = clsmap(k,4)
               m = permute(k,4)
               base_sol(ofl:ofu,lev,j) = solution(:cols,m)
            end do
         end do Column_loop
      end do Level_loop

      end subroutine imp_sol

      end module mo_imp_sol
