
      module charge_neutrality

      use shr_kind_mod,      only : r8 => shr_kind_r8
      use cam_logfile,       only : iulog

      implicit none

      private
      public :: charge_balance
      public :: charge_fix     ! temporary, for fixing charge balance after vertical diffusion
                               ! without converting mass mixing ratios to volume
                               ! mean mass assumed to be mwdry

      contains

      subroutine charge_balance( ncol, conc )
!-----------------------------------------------------------------------      
!        ... force ion/electron balance
!-----------------------------------------------------------------------      

        use ppgrid,       only : pver
        use mo_chem_utls, only : get_spc_ndx
        use chem_mods,    only : gas_pcnst

        implicit none

!-----------------------------------------------------------------------      
!        ... dummy arguments
!-----------------------------------------------------------------------      
      integer,  intent(in)          :: ncol
      real(r8), intent(inout)       :: conc(ncol,pver,gas_pcnst)         ! concentration

!-----------------------------------------------------------------------      
!        ... local variables
!-----------------------------------------------------------------------      
      integer  :: k, n
      integer  :: elec_ndx
      real(r8) :: wrk(ncol,pver)

      elec_ndx = get_spc_ndx('e')
#ifdef CB_DEBUG
      write(iulog,*) ' '
      write(iulog,*) '------------------------------------------------------------------'
      write(iulog,*) 'charge_balance: e ndx,offset = ',elec_ndx,offset
      write(iulog,*) 'charge_balance: size of conc = ',size(conc,dim=1),' x ',size(conc,dim=2),' x ',size(conc,dim=3)
#endif
      if( elec_ndx > 0 ) then
	 wrk(:,:) = 0._r8
         n = get_spc_ndx('Np')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
         n = get_spc_ndx('N2p')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
         n = get_spc_ndx('Op')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
         n = get_spc_ndx('O2p')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
         n = get_spc_ndx('NOp')
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	    end do
         end if
#ifdef CB_DEBUG
         write(iulog,*) 'charge_balance: electron concentration before balance'
         write(iulog,'(1p,5g15.7)') conc(1,:,elec_ndx)
         write(iulog,*) 'charge_balance: electron concentration after  balance'
         write(iulog,'(1p,5g15.7)') wrk(1,:)
         write(iulog,*) '------------------------------------------------------------------'
         write(iulog,*) ' '
#endif
         conc(:ncol,:,elec_ndx) = wrk(:ncol,:)
      end if

      end subroutine charge_balance


      subroutine charge_fix( ncol, q )
!-----------------------------------------------------------------------      
!        ... force ion/electron balance
!-----------------------------------------------------------------------      

      use ppgrid,       only : pver
      use constituents, only : cnst_get_ind, cnst_mw
      use physconst,    only : mwdry                   ! molecular weight of dry air
 
      implicit none

!-----------------------------------------------------------------------      
!        ... dummy arguments
!-----------------------------------------------------------------------      
      integer, intent(in)           :: ncol
      real(r8), intent(inout)       :: q(:,:,:)        ! model mass mixing ratios

!-----------------------------------------------------------------------      
!        ... local variables
!-----------------------------------------------------------------------      
      integer  :: k, n
      integer  :: elec_ndx
      real(r8) :: wrk(ncol,pver)
      real(r8) :: mbar(ncol,pver)  ! mean mass (=mwdry) used to fake out optimizer to get
                                   ! identical answers to old code

!-----------------------------------------------------------------------      

! assume that mbar = mwdry
      mbar(:,:) = mwdry

      call cnst_get_ind( 'e', elec_ndx, abort=.false. )
      if( elec_ndx > 0 ) then
	 wrk(:,:) = 0._r8
         call cnst_get_ind( 'Np', n, abort=.false. )
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'N2p', n, abort=.false. )
         if( n > 0 ) then
	    do k = 1,pver
	      wrk(:,k) = wrk(:,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'Op', n, abort=.false. )
         if( n > 0 ) then
	    do k = 1,pver
              wrk(:,k) = wrk(:,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'O2p', n, abort=.false. )
         if( n > 0 ) then
	    do k = 1,pver
              wrk(:,k) = wrk(:,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'NOp', n, abort=.false. )
         if( n > 0 ) then
	    do k = 1,pver
              wrk(:,k) = wrk(:,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
	 do k = 1,pver
	   q(:ncol,k,elec_ndx) = cnst_mw(elec_ndx) * wrk(:ncol,k) / mbar(:ncol,k)
	 end do
      end if

      end subroutine charge_fix

      end module charge_neutrality
