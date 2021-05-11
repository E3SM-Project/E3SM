
      module mo_chem

      implicit none

      character(len=256) :: buff
      character(len=256) :: buffh

      private
      public :: chem

      contains

      subroutine chem
!-----------------------------------------------------------------------
!	... Scan chemical reactions and produce base chemistry maps
!-----------------------------------------------------------------------

      use io,      only : lin, lout
      use var_mod, only : solsym, fixsym, pcesym, &
			  nq, nfs, spcsym, spccnt, var_lim
      use rxt_mod, only : rxno => rxntot, irc => rxmcnt, &
			  prdcnt, ipcep, ipcel, fixcnt, prdmap, &
			  fixmap, pcep, pcel, rxmap, hetcnt, &
			  hetmap, usrcnt, usrmap, rates => rxparm, &
			  troe_rxparm, troetab, troecnt, troe_sym_rates, &
			  rattab => rxptab, rateno => rxpcnt, &
			  pcoeff_cnt, pcoeff_ind, pcoeff, sym_rates, &
			  phtsym, phtcnt, pht_alias, pht_alias_mult, rxt_lim, rxt_tag, &
                          prd_lim, rxtnt_lim
      use rxt_mod, only : has_cph => cph_flg
      use rxt_mod, only : frc_from_dataset

      implicit none

!-----------------------------------------------------------------------
!	... Local variables
!           nsr = number of solution reactants
!           nsp = number of solution products
!           nf  = number of fixed reactants
!           npr = number of pce reactants
!           npp = number of pce products
!-----------------------------------------------------------------------
      integer, parameter :: photolysis = 1, gas_phase = 2
      integer, parameter :: heterogeneous = 3, extraneous = 4

      character(len=16) :: param
      character(len=16) :: rxparms(prd_lim)
      character(len=16) :: sym_rate(5)
      character(len=16) :: keywords(4) = (/ 'PHOTOLYSIS      ', &
	                                    'REACTIONS       ', &
	                                    'HETEROGENEOUS   ', &
	                                    'EXTFORCING      ' /)

!-----------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------
      integer  ::  nchar, k, nr, np, nsr, nsp, nf, &
                   npr, npp, ic, kc, i, npl, l, j, m, ipp, im, &
                   photo, rxtn, npce = 0, ipl
      integer  ::  il, iu, istat
      integer  ::  beg_alias
      integer  ::  rxttab(5,prd_lim)
      integer  ::  parsw(4), &
                   count1(var_lim), &
                   count2(var_lim)
      integer, allocatable :: toklen(:)
      integer  ::  tokcnt
      
      character(len=128) :: line
      character(len=32) ::  loc_rxt_tag
      character(len=32) ::  loc_pht_alias(2)
      character(len=16) ::  wrk_char
      character(len=16)  ::  rxtsym(rxtnt_lim)
      character(len=16)  ::  prdsym(prd_lim)
      character(len=1)  ::  char
      character(len=24), allocatable ::  ext_tokens(:)
      character(len=16),  allocatable ::  tokens(:)

      real  ::     number
      real  ::     rate(5), pcoeffs(prd_lim)

      logical  ::  cph_flg
      logical  ::  coeff_flg
      logical  ::  found
       
      photo  = 0
      usrcnt = 0
      rxtn   = 1
      rxno   = 0
      ipl    = 0
      rateno = 0
      parsw  = 0
      
keyword_loop : &
      do
         call cardin( lin, buff, nchar )
         buffh = buff
         call upcase( buffh )
         if( buffh == 'ENDCHEMISTRY' ) then
	    exit
         end if
         found = .false.
         do i = 1,4
	    if( buffh == keywords(i) ) then
	       if( parsw(i) /= 0 ) then
                  call errmes ( ' # Keyword already used@', &
                                lout, &
                                keywords(i), &
                                len_trim(keywords(i)), &
                                buff )
	       else if( i == 1 .and. parsw(2) /= 0 ) then
                  call errmes ( 'Must specify Photolysis before Reactions@', &
                                lout, char, 1, buff )
	       end if
	       parsw(i) = 1
	       found = .true.
	       exit
	    end if
         end do
         if( .not. found ) then
            call errmes( 'CHEM: # is not a keyword@', &
                          lout, &
                          buff, &
                          len_trim(buff), &
                          buff )
	 end if
         select case( i )
	    case( photolysis )
!=======================================================================
!     	... The photolysis chemistry processing
!=======================================================================
               write(lout,*) ' '
               write(lout,240)
photolysis_loop : &
	       do
                  call cardin( lin, buff, nchar )
                  buffh = buff
                  call upcase( buffh )
                  if( buffh == 'ENDPHOTOLYSIS' ) then
                     phtcnt = rxno
!-----------------------------------------------------------------------
!        ... Check that all photorates have a reaction tag
!-----------------------------------------------------------------------
                     if( phtcnt > 0 ) then
                        k = count(  rxt_tag(:phtcnt) /= ' ' )
                        if( k /= phtcnt ) then
                           call errmes( 'All photoreactions must have a reaction tag@', lout, char, 1, buff )
                        end if
                        do ic = 1,phtcnt
                           if( any( pht_alias(ic,:) /= ' ' ) ) then
                              line = ' '
                              m = len_trim(rxt_tag(ic))
                              line = rxt_tag(ic)(:m) // ' -> '
                              do k = 1,2
                                 if( pht_alias(ic,k) /= ' ' ) then
                                    if( k == 2 .and. pht_alias(ic,1) /= ' ' ) then
                                       line(len_trim(line)+1:) = ','
                                    end if
                                    m = len_trim(pht_alias_mult(ic,k))
                                    line(len_trim(line)+2:) = pht_alias_mult(ic,k)(:m) // ' * '
                                    m = len_trim(pht_alias(ic,k))
                                    line(len_trim(line)+2:) = pht_alias(ic,k)(:m)
                                 end if
                              end do
                           end if
                        end do
                     end if
                     cycle keyword_loop
                  end if
                  rxtsym(1:3) = ' '
                  prdsym(1:4) = ' '
                  sym_rate(:) = ' '
!-----------------------------------------------------------------------
!        ... Reaction parsing routine
!-----------------------------------------------------------------------
                  call rxtprs( nchar, nr, np, rxtsym, prdsym, &
                               rate, pcoeffs, coeff_flg, rxparms, sym_rate, &
                               loc_rxt_tag, cph_flg, .true. )
!-----------------------------------------------------------------------
!        ... Check for reaction string errors from parsing routine
!-----------------------------------------------------------------------
                  if( nr < 0 ) then
                     call errmes ( 'gross syntax errors in reaction string@', lout, char, 1, buff )
                  end if

!-----------------------------------------------------------------------
!        ... Reaction mapping routine
!-----------------------------------------------------------------------
                  call mapper( nq,       nfs,      npce,     nr,       np, &
                               rxtsym,   prdsym,   solsym,   fixsym,   pcesym, &
                               nsr,      nsp,      nf,       npr,      npp, &
                               rxttab,   photo,    coeff_flg,pcoeffs )

!-----------------------------------------------------------------------
!        ... Check for logic errors in reaction
!-----------------------------------------------------------------------
                  if( (nf + nsr + npr) == 0 ) then
                     call errmes ( 'photo-reaction has no reactants@', lout, char, 1, buff )
                  else if( (nf+nsr+npr) >= 2 ) then
                     call errmes ( 'photo-reaction has two or more reactants@', lout, char, 1, buff )
                  end if
                  if( nf == 1 ) then
                     if( nsp == 0 .and. npp == 0 ) then
                        call errmes( 'fixed species photolysis produces nothing@', lout, char, 1, buff )
                     end if
                     rxno = rxno + 1
!-----------------------------------------------------------------------
!        ... Photolysis of an invariant species; check for products
!-----------------------------------------------------------------------
                     if( nsp /= 0 ) then
!-----------------------------------------------------------------------
!        ... Solution species production from fixed species photolysis
!-----------------------------------------------------------------------
                       prdcnt = prdcnt + 1
                       prdmap(prdcnt,1) = rxno
                       prdmap(prdcnt,2:nsp+1) = rxttab(3,1:nsp)
                     end if
                     if( npp /= 0 ) then
!-----------------------------------------------------------------------
!        ... Pce species production from fixed species photolysis
!-----------------------------------------------------------------------
                       ipl = ipl + 1
                       do k = 1,npp
                         ipcep(1) = ipcep(1) + 1
                         ic = ipcep(1)
                         pcep(ic,1,1) = rxttab(5,k)
                         pcep(ic,2,1) = rxno
                         pcep(ic,3,1) = ipl
                       end do
                     end if
!-----------------------------------------------------------------------
!        ... Set the fixed reactants map
!-----------------------------------------------------------------------
                     fixcnt(1) = fixcnt(1) + 1
                     ic = fixcnt(1)
                     fixmap(ic,1,1) = -rxno
                     fixmap(ic,2,1) = rxttab(1,1)
                     kc = rxttab(1,1)
                     phtsym(rxno) = fixsym(kc)
                  else if( npr == 1 ) then
!-----------------------------------------------------------------------
!        ... Pce species photolysis; put into pce loss map
!-----------------------------------------------------------------------
                     if( npp /= 0 ) then
!-----------------------------------------------------------------------
!        ... Check to see that a pce is not a product species
!-----------------------------------------------------------------------
                         call errmes( 'pce reactants and products in same reaction@', lout, char, 1, buff )
                     end if
                     rxno = rxno + 1
                     ipcel(1) = ipcel(1) + 1
                     ic = ipcel(1)
                     pcel(ic,1,1) = rxttab(4,1)
                     pcel(ic,2,1) = rxno
                     kc = rxttab(4,1)
                     phtsym(rxno) = pcesym(kc)
                     if( nsp /= 0 ) then
                        pcel(ic,3:nsp+2,1) = rxttab(3,1:nsp)
                     end if
!-----------------------------------------------------------------------
!        ... Solution species photolysis
!-----------------------------------------------------------------------
                  else
                     rxno = rxno + 1
                     irc(1) = irc(1) + 1
                     ic = irc(1)
                     rxmap(ic,2,1) = rxttab(2,1)
                     rxmap(ic,1,1) = rxno
                     kc = rxttab(2,1)
                     phtsym(rxno) = solsym(kc)
                     if( nsp /= 0 ) then
!---------------------------------------------------------------------
!        ... Solution species production from solution photolysis
!---------------------------------------------------------------------
                        rxmap(ic,3:2+nsp,1) = rxttab(3,1:nsp)
                     end if
                     if( npp /= 0 ) then
!----------------------------------------------------------------------
!        ... Pce species production from solution photolysis
!----------------------------------------------------------------------
                        ipl = ipl + 1
                        do k = 1,npp
                          ipcep(2) = ipcep(2)+1
                          ic = ipcep(2)
                          pcep(ic,1,2) = rxttab(5,k)
                          pcep(ic,2,2) = rxno
                          pcep(ic,3,2) = rxttab(2,1)
                          pcep(ic,4,2) = ipl
                        end do
                     end if
                  end if
                  if( rxno > rxt_lim ) then
                     call errmes( ' Reaction count exceeds limit@', lout, buff, 1, buff )
                  end if
!-----------------------------------------------------------------------
!        ... Check for non-unity product coefficients
!-----------------------------------------------------------------------
                  if( coeff_flg ) then
                     pcoeff_cnt = pcoeff_cnt + 1
                     pcoeff_ind(rxno) = pcoeff_cnt
                     pcoeff(1:nsp,pcoeff_cnt) = pcoeffs(1:nsp)
                  end if
!-----------------------------------------------------------------------
!        ... Check for photorate alias
!-----------------------------------------------------------------------
                 m = index( loc_rxt_tag, '=' )
                 if( m /= 0 ) then
                    beg_alias = m + 1
                 else
                    m = index( loc_rxt_tag, '->' )
                    if( m > 0 ) then
                       beg_alias = m + 2
                    end if
                 end if
                 if( m > 0 ) then
                    loc_pht_alias(1) = loc_rxt_tag(beg_alias:)
                    loc_rxt_tag(m:)  = ' '
                    m = index( loc_pht_alias(1), ',' )
                    il = 1
                    iu = 2
                    if( m == 0 ) then
                       iu = 1
                    else if( m == 1 ) then
                       il = 2
                       loc_pht_alias(2) = loc_pht_alias(1)(2:)
                    else
                       loc_pht_alias(2) = loc_pht_alias(1)(m+1:)
                       loc_pht_alias(1) = loc_pht_alias(1)(:m-1)
                    end if
                    do ic = il,iu
                       k = index( loc_pht_alias(ic), '*' )
                       if( k == 0 ) then
		          pht_alias(rxno,ic)      = loc_pht_alias(ic)
                       else
		          pht_alias_mult(rxno,ic) = loc_pht_alias(ic)(:k-1)
		          read(pht_alias_mult(rxno,ic),*,iostat=istat) number
                          if( istat /= 0 ) then
                             call errmes ( ' # is not a valid number@', &
                                           lout, &
                                           pht_alias_mult(rxno,ic), &
                                           len_trim(pht_alias_mult(rxno,ic)), &
                                           buff )
                          end if
		          pht_alias(rxno,ic)      = loc_pht_alias(ic)(k+1:)
                       end if
                    end do
                 end if
!-----------------------------------------------------------------------
!        ... Check for duplicate reaction tag
!-----------------------------------------------------------------------
		  do m = 1,rxno-1
		     if( trim( loc_rxt_tag ) /= ' ' ) then
		        if( trim( rxt_tag(m) ) == trim( loc_rxt_tag ) ) then
                           call errmes ( ' # rxtnt alias already in use@', &
                                         lout, &
                                         loc_rxt_tag, &
                                         len_trim(loc_rxt_tag), &
                                         buff )
		        end if
		     end if
		  end do
		  rxt_tag(rxno) = loc_rxt_tag
                  has_cph(rxno)   = cph_flg
!-----------------------------------------------------------------------
!        ... Print the reaction on unit lout
!-----------------------------------------------------------------------
                   call outp( rxparms, nr, np, rxtsym, prdsym, sym_rate, rxno, rate, loc_rxt_tag, lout )
	       end do photolysis_loop
      
	    case( gas_phase )
!=======================================================================
!    	... The chemical reactions
!=======================================================================
               write(lout,*) ' '
               write(lout,260)
gas_phase_rxt_loop : &
	       do
                  call cardin( lin, buff, nchar )
                  buffh = buff
                  call upcase( buffh )
                  if( buffh == 'ENDREACTIONS' ) then
	             cycle keyword_loop
                  end if

                  rxtsym(1:3) = ' '
                  prdsym(1:4) = ' '
		  sym_rate = ' '
                  call rxtprs( nchar, nr, np, rxtsym, prdsym, &
                               rate, pcoeffs, coeff_flg, rxparms, sym_rate, &
                               loc_rxt_tag, cph_flg, .false. )

                  if( nr < 0 ) then
                     call errmes ( 'there are no reactants@', lout, char, 1, buff )
                  end if

                  call mapper( nq,       nfs,      npce,     nr,       np, &
                               rxtsym,   prdsym,   solsym,   fixsym,   pcesym, &
                               nsr,      nsp,      nf,       npr,      npp, &
                               rxttab,   rxtn,     coeff_flg,pcoeffs )

                  if( nsr == 3 ) then
                     call errmes ( ' three solution species reactants@', lout, char, 1, buff )
                  end if
                  if( npr >= 2 ) then
                     call errmes ( ' there are two or more pce reactants@', lout, char, 1, buff )
                  end if
                  if( nsr == 2 .and. npr /= 0 ) then
                        call errmes ( ' there are two solution and one pce reactants@', lout, char, 1, buff )
                  end if
                  if( (nf+nsr+npr) /= nr ) then
                     call errmes ( ' reaction parsing algorithm error@', lout, char, 1, buff )
                  end if
                  rxno = rxno + 1  
                  if( rxno > rxt_lim ) then
                     call errmes( ' Reaction count exceeds limit@', lout, buff, 1, buff )
                  end if
!-----------------------------------------------------------------------
!	... User specified reaction rate ?
!-----------------------------------------------------------------------
                  if( sym_rate(1) /= ' ' ) then
		     if( sym_rate(3) == ' ' ) then
                        rateno = rateno + 1
                        rattab(rateno) = rxno
                        rates(:2,rateno) = rate(:2)
		        sym_rates(:2,rateno) = sym_rate(:2)
		     else if( sym_rate(5) /= ' ' ) then
                        troecnt = troecnt + 1
                        troetab(troecnt) = rxno
                        troe_rxparm(:,troecnt) = rate(:)
		        troe_sym_rates(:,troecnt) = sym_rate(:)
                     end if
                  end if
!-----------------------------------------------------------------------
!        ... Check for duplicate reaction tag
!-----------------------------------------------------------------------
		  do m = 1, rxno-1
		     if( trim( loc_rxt_tag ) /= ' ' ) then
		        if( trim( rxt_tag(m) ) == trim( loc_rxt_tag ) ) then
                           call errmes ( ' # rxtnt alias already in use@', &
                                         lout, &
                                         loc_rxt_tag, &
                                         len_trim(loc_rxt_tag), &
                                         buff )
		        end if
		     end if
		  end do
		  rxt_tag(rxno) = loc_rxt_tag
                  has_cph(rxno)   = cph_flg
                  call outp( rxparms, nr, np, rxtsym, prdsym, sym_rate, rxno-phtcnt, rate, loc_rxt_tag, lout )

                  if( nf /= 0 ) then
                     fixcnt(nf) = fixcnt(nf) + 1  
                     ic = fixcnt(nf)
                     fixmap(ic,1,nf) = rxno
                     fixmap(ic,2:1+nf,nf) = rxttab(1,1:nf)
                     if( nf == nr ) then
!-----------------------------------------------------------------------
!        ... Fixed reactants only
!-----------------------------------------------------------------------
                        prdcnt = prdcnt+1
                        fixmap(ic,1,nf) = -rxno
                        prdmap(prdcnt,1) = rxno
                        if( nsp /= 0 ) then
                           prdmap(prdcnt,2:1+nsp) = rxttab(3,1:nsp)
                        end if
                        if( npp /= 0 ) then
                           ipl = ipl + 1
                           do k = 1,npp
                             ipcep(1) = ipcep(1) + 1
                             ic = ipcep(1)
                             pcep(ic,1,1) = rxttab(5,k)
                             pcep(ic,2,1) = rxno
                             pcep(ic,3,1) = ipl
                           end do
                        end if
                        cycle
                     end if
                  end if

                  if( nsr == nr .or. npr == 0 ) then
!-----------------------------------------------------------------------
!        ... Solution reactants only
!-----------------------------------------------------------------------
                     irc(nsr) = irc(nsr) + 1
                     ic = irc(nsr)
                     rxmap(ic,1,nsr) = rxno
                     rxmap(ic,2:1+nsr,nsr) = rxttab(2,1:nsr)
                     if( nsp /= 0 ) then
                        rxmap(ic,nsr+2:nsr+nsp+1,nsr) = rxttab(3,1:nsp)
                     end if
                     if( npp /= 0 ) then
                        ipl = ipl + 1
                        npl = nsr + 1
                        do k = 1,npp
                          ipcep(npl) = ipcep(npl) + 1
                          ic = ipcep(npl)
                          pcep(ic,1,npl) = rxttab(5,k)
                          pcep(ic,2,npl) = rxno
                          pcep(ic,nsr+3,npl) = ipl
                          pcep(ic,3:2+nsr,npl) = rxttab(2,1:nsr)
                        end do
                     end if
                  else
!-----------------------------------------------------------------------
!        Solution,fixed, and pce reactants possible.
!        There is a pce reactant, and either a solution/fixed
!        reactant or both.  If there is a pce product then
!        terminate with a reaction logic error.
!-----------------------------------------------------------------------
                     if( npp /= 0 ) then
                        call errmes( ' Reaction has both reactant and product pce species@', lout, char, 1, buff )
                     end if
                     npl = nsr + 1
                     ipcel(npl) = ipcel(npl)+1
                     ic = ipcel(npl)
                     pcel(ic,1,npl) = rxttab(4,1)
                     pcel(ic,2,npl) = rxno
                     if( npl == 2 ) then
			pcel(ic,3,2) = rxttab(2,1)
		     end if
                     if( nsp == 0 ) then
			cycle
		     end if
                     pcel(ic,npl+2:npl+nsp+1,npl) = rxttab(3,1:nsp)
                  end if
!-----------------------------------------------------------------------
!        ... Check for non-unity product coefficients
!-----------------------------------------------------------------------
                  if( coeff_flg ) then
                     pcoeff_cnt = pcoeff_cnt + 1
                     pcoeff_ind(rxno) = pcoeff_cnt
                     pcoeff(1:nsp,pcoeff_cnt) = pcoeffs(1:nsp)
                  end if
               end do gas_phase_rxt_loop

	    case( heterogeneous )
!=======================================================================
!        The heterogeneous loss chemistry list
!=======================================================================
               write(lout,*) ' '
               write(lout,7175)
	       allocate( toklen(64) )
	       allocate( tokens(64) )
hetero_loop :  do
                  call cardin( lin, buff, nchar )
                  buffh = buff
                  call upcase( buffh )
                  if( buffh == 'ENDHETEROGENEOUS' ) then
		     deallocate( toklen )
		     deallocate( tokens )
                     cycle keyword_loop
                  end if

		  call gettokens( buff, nchar, ',', 8, tokens, toklen, 64, tokcnt )
		  if( tokcnt <= 0 ) then
		     deallocate( toklen )
		     deallocate( tokens )
                     call errmes ( ' Error in het list@', lout, param, k, buff )
		  end if
               end do hetero_loop

	    case( extraneous )
!=======================================================================
!    	... The extraneous prod/loss chemistry list
!=======================================================================
               write(lout,*) ' '
               write(lout,8175)
               allocate( toklen(64) )
               allocate( ext_tokens(64) )
extfrc_loop :  do
                  call cardin( lin, buff, nchar )
                  buffh = buff
                  call upcase( buffh )
                  if( buffh == 'ENDEXTFORCING' ) then
                     deallocate( toklen )
                     deallocate( ext_tokens )
                     cycle keyword_loop
                  end if

                  call gettokens( buff, nchar, ',', 24, ext_tokens, toklen, 64, tokcnt )
                  if( tokcnt <= 0 ) then
                     deallocate( toklen )
                     deallocate( ext_tokens )
                     call errmes ( ' Error in ext prod list@', lout, param, k, buff )
                  end if
ext_tok_loop :    do j = 1,tokcnt
                     do m = 1,spccnt(1)
                        k = index( ext_tokens(j), '<' ) - 1
                        if( k < 1 ) then
                           found = .false.
                           k = len_trim( ext_tokens(j) )
                        else
                           found = .true.
                        end if
                        if( ext_tokens(j)(:k) == spcsym(m,1) ) then
                           if( usrcnt > 1 ) then
                              if( any(usrmap(:usrcnt) == m ) ) then
                                 call errmes( ' # is already in ext frc list@', lout, ext_tokens(j), toklen(j), buff )
			      end if
			   end if
                           usrcnt = usrcnt + 1
                           if( usrcnt > rxt_lim ) then
                              call errmes( ' Extran reaction count exceeds limit@', lout, buff, 1, buff )
                           end if
                           usrmap(usrcnt) = m
                           !write(lout,'(1x,''('',i2,'')'',3x,a)') usrcnt, trim(spcsym(m,1))

                           if( .not. found ) then
                              write(lout,'(1x,''('',i2,'')'',3x,a)') usrcnt, trim(spcsym(m,1))
                           else
                              frc_from_dataset(usrcnt) = .true.
                              write(lout,'(1x,''('',i2,'')'',3x,a,3x,''(dataset)'')') usrcnt, trim(spcsym(m,1))
                           end if


                           cycle ext_tok_loop
                        end if
                     end do
                     call errmes ( ' # is not in Solution list@', lout, ext_tokens(j), toklen(j), buff )
                  end do ext_tok_loop
               end do extfrc_loop
         end select
      end do keyword_loop
      hetcnt = spccnt(1)
      do j=1,hetcnt
        hetmap(j,1) = j
      enddo
!-----------------------------------------------------------------------
!        ... Check for pce,sol reaction validity
!-----------------------------------------------------------------------
      il = ipcel(2)
      ipp = ipcep(3)
      if( il == 0 .or. ipp == 0 ) then
	 return
      end if
      do i = 1,il
         im = pcel(i,1,2)
         count1(im) = 1
      end do
      do i = 1,ipp
         im = pcep(i,1,3)
         count2(im) = 1
      end do
      do i = 1,npce
         if( count1(i) == 1 .and. count2(i) == 1 ) then
            call errmes ( 'pce species # violates use rules@', lout, pcesym(i), 8, buff )
         end if
         count1(i) = 0
      end do

240   format(5x,'Photolysis')
260   format(5x,'Reactions')
7175  format('Heterogeneous loss species')
7177  format(1x,'(',i2,')',3x,a8)
8175  format('Extraneous prod/loss species')
!-----------------------------------------------------------------------
!        ... End of the chemistry processing code
!-----------------------------------------------------------------------

      end subroutine chem

      subroutine rxtprs( nchar, &
                         rxtcnt, &
                         prdcnt, &
                         rxtsym, &
                         prdsym, &
                         rate, &
                         pcoeffs, &
                         coeff_flg, &
                         prdprms, &
                         sym_rate, &
                         loc_rxt_tag, &
                         cph_flg, &
			 is_photorate )

      use io, only : lin, lout
      use rxt_mod, only : rxtnt_lim, prd_lim

      implicit none

!-----------------------------------------------------------------------
!  	... Rxtprs parses the reaction and places the symbols
!           in the reactant and product character arrays rxtsym &
!           prdsym.  The reactants and products are checked for
!           symbol length violations.  The number of reactants
!           and products are limited to three each.  Reaction rate
!           information is checked for numeric validity.  There
!           must be at least one reactant.  No other error checking
!           is performed in this subroutine.
!-----------------------------------------------------------------------

      integer, parameter ::  symlen = 16

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)  ::     nchar
      integer, intent(out) ::     rxtcnt, prdcnt
      real, intent(out)    ::     rate(*), pcoeffs(prd_lim)
      character(len=*), intent(out)  ::  loc_rxt_tag
      character(len=16), intent(out) ::  rxtsym(rxtnt_lim), prdsym(prd_lim)
      character(len=16), intent(out) :: prdprms(prd_lim)
      character(len=16), intent(out) :: sym_rate(*)
      logical, intent(in)  ::     is_photorate
      logical, intent(out) ::     coeff_flg
      logical, intent(out) ::     cph_flg
      
      
!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::   retcod
      integer  ::   ncharl, tprdcnt
      integer  ::   comma
      integer  ::   k, j, length, ratcnt, start, position
      integer, allocatable  ::   slen(:)
      character(len=320) :: buffl, buffhl
      character(len=16), allocatable  :: rxparms(:)

      rxtcnt = 0
      prdcnt = 0
      coeff_flg = .false.
      cph_flg   = .false.
      
      rate(:5)    = 0.
      pcoeffs(:)  = 1.
      ncharl      = nchar
      allocate( slen(MAX(5,prd_lim)) )
      allocate( rxparms(MAX(5,prd_lim)) )
!-----------------------------------------------------------------------
!	... Check for reaction name alias, cph
!-----------------------------------------------------------------------
      k = index( buff(:ncharl), ']' ) - 1
      if( k > 0 ) then
         j = index( buff(:ncharl), '[' ) + 1
         loc_rxt_tag = buff(j:k)
         comma = index( buff(j:k), ',' )
         if( comma /= 0 ) then
            if( buff(j+comma:k) == 'cph' ) then
               cph_flg = .true.
               loc_rxt_tag = buff(j:j+comma-2)
            end if
         end if
         buff = buff(k+2:ncharl)
	 ncharl = nchar - (k+1)
      else
         loc_rxt_tag = ' '
      end if

      length = index( buff(:ncharl), '=' )
      if( length == 0 ) then
         length = index( buff(:ncharl), '->' )
      end if

      if( length <= 1 ) then
         if( length == 0 ) then
            write(6,102)
         else if( length == 1 ) then
            write(6,100)
         end if
         rxtcnt = -1
         return
      end if
      
!-----------------------------------------------------------------------
!        ... Parse out the reactants
!-----------------------------------------------------------------------
      call gettokens( buff, &
                      length-1, &
                      '+', &
                      symlen, &
                      rxtsym, &
                      slen, &
                      rxtnt_lim, &
                      rxtcnt )
      if( rxtcnt < 0 ) then
         call errmes( 'Too many reactants in reaction@', lout, buff, 1, buff )
      else if ( rxtcnt == 0 ) then
         call errmes( 'Reactant symbol exceeds symbol length@', lout, buff, 1, buff )
      end if

      if( index( buff(:ncharl),'->' ) /= 0 ) then
	 length = length + 1
      end if
      if( length == ncharl ) then       ! reaction has reactants only
         return
      end if
      start = length + 1               ! char after "=" sign
      position = index( buff(start:ncharl), ';' )
      if( position == 0 ) then    ! products, no rates
         length = ncharl - start + 1
      else if( position /= 1 ) then
         length = position - 1
      end if
      
!-----------------------------------------------------------------------
!        ... Parse out the products and multipliers
!-----------------------------------------------------------------------
      if( position /= 1 ) then
         call gettokens( buff(start:), &
                         length, &
                         '+', &
                         16, &
                         prdprms, &
                         slen, &
                         prd_lim, &
                         prdcnt )

         if( prdcnt < 0 ) then
            call errmes( 'Too many products in reaction@', lout, buff, 1, buff )
         else if ( prdcnt == 0 ) then
            call errmes( 'Product symbol exceeds symbol length@', lout, buff, 1, buff )
         end if
!-----------------------------------------------------------------------
!        ... Check each "product" token for an explicit multiplier
!-----------------------------------------------------------------------
         do k = 1,prdcnt
            j = index( prdprms(k)(:slen(k)), '*' )
            if( j == 0 ) then
               prdsym(k) = prdprms(k)(:symlen)
               cycle
            else if( j == 1 .or. j == slen(k) ) then
               call errmes( ' Product & multiplier syntax error@', lout, prdprms(k), 1, prdprms(k) )
            end if
            read(prdprms(k)(:j-1),*,iostat=retcod) pcoeffs(k)
            if( retcod /= 0 ) then
	       call errmes( 'number format error in product multiplier #@', lout, &
                             prdprms(k), slen(k), buff )
	    end if
            prdsym(k) = prdprms(k)(j+1:)
            if( pcoeffs(k) /= 1. ) then
               coeff_flg = .true.
            end if
         end do
      end if

!-----------------------------------------------------------------------
!        ... Set any reaction rate parms
!-----------------------------------------------------------------------
      if( position /= 0 ) then
         start = start + position
         call gettokens( buff(start:), &
                         ncharl-start+1, &
                         ',', &
                         16, &
                         rxparms, &
                         slen, &
                         5, &
                         ratcnt )

         if( ratcnt <= 0 ) then
            call errmes( ' Syntax error in reaction rate parameters@', lout, buff, 1, buff )
         end if
         do k = 1,ratcnt
	    if( rxparms(k) /= ' ' ) then
	       read(rxparms(k)(:slen(k)),*,iostat=retcod) rate(k)
               if( retcod /= 0 ) then
	          call errmes( 'number format error in reaction rate #@', lout, rxparms(k), slen(k), buff )
	       end if
	       sym_rate(k) = rxparms(k)(:slen(k))
	    else
	       rate(k) = 0.
	       sym_rate(k) = ' '
	    end if
         end do
      end if

!-----------------------------------------------------------------------
!        ... New code for extended product lines
!-----------------------------------------------------------------------
      do
         call cardin( lin, buffl, ncharl )
         buffhl = buffl
         call upcase( buffhl )
         if( .not. is_photorate .and. buffhl == 'ENDREACTIONS' ) then
            backspace( lin )
	    exit
         else if( is_photorate .and. buffhl == 'ENDPHOTOLYSIS' ) then
            backspace( lin )
	    exit
         else
            length = index( buffl(:ncharl), '=' )
            if( length == 0 ) then
               length = index( buffl(:ncharl), '->' )
            end if
            if( length /= 0 ) then
               backspace( lin )
	       exit
            end if
         end if
	 if( buffl(1:1) /= '+' ) then
            call errmes ( 'Extended Reactions must start with + delimiter@', lout, buffl, 1, buffl )
	 end if
         if( prdcnt >= prd_lim ) then
            call errmes( 'Too many products in reaction@', lout, buffl, 1, buffl )
	 end if
!-----------------------------------------------------------------------
!        ... Parse out the products and multipliers
!-----------------------------------------------------------------------
         call gettokens( buffl(2:), &
                         ncharl-1, &
                         '+', &
                         16, &
                         prdprms(prdcnt+1), &
                         slen(prdcnt+1), &
                         prd_lim-prdcnt, &
                         tprdcnt )
         if( tprdcnt < 0 ) then
            call errmes( 'Too many products in reaction@', lout, buffl, 1, buffl )
         else if ( tprdcnt == 0 ) then
            call errmes( 'Product symbol exceeds symbol length@', lout, buffl, 1, buffl )
         end if
!-----------------------------------------------------------------------
!        ... Check each "product" token for an explicit multiplier
!-----------------------------------------------------------------------
         do k = prdcnt+1,prdcnt+tprdcnt
            j = index( prdprms(k)(:slen(k)), '*' )
            if( j == 0 ) then
               prdsym(k) = prdprms(k)(:symlen)
               cycle
            else if( j == 1 .or. j == slen(k) ) then
               call errmes( ' Product & multiplier syntax error@', lout, prdprms(k), 1, prdprms(k) )
            end if
            read(prdprms(k)(:j-1),*,iostat=retcod) pcoeffs(k)
            if( retcod /= 0 ) then
	       call errmes( 'number format error in product multiplier #@', lout, &
                             prdprms(k), slen(k), buff )
	    end if
            prdsym(k) = prdprms(k)(j+1:)
            if( pcoeffs(k) /= 1. ) then
               coeff_flg = .true.
            end if
         end do
	 prdcnt = prdcnt + tprdcnt
      end do

      deallocate( slen )
      deallocate( rxparms )

!-----------------------------------------------------------------------
!        ... Formats
!-----------------------------------------------------------------------
100   format('0 **** reaction string has no reactants ****')
102   format('0 **** reaction string has no  =  sign separator ****')

      end subroutine rxtprs

      subroutine mapper( nsol,     nfix,     npce,     rxtcnt, &
                         prdcnt,   rxtsym,   prdsym,   solsym, &
                         fixsym,   pcesym,   nsr,      nsp, &
                         nf,       npr,      npp,      rxttab, &
                         rxtype,   coeff_flg,pcoeffs )

      use io, only : lout
      use rxt_mod, only : rxtnt_lim, prd_lim

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)  :: nsol, nfix, npce, &
                              rxtcnt, prdcnt, rxtype
      integer, intent(out) :: nsr, nsp, nf, npr, npp
      integer, intent(out) :: rxttab(5,*)
      real, intent(inout)  :: pcoeffs(prd_lim)
      logical, intent(inout) :: coeff_flg
      character(len=*), intent(in)  ::  rxtsym(rxtnt_lim), prdsym(prd_lim)
      character(len=*), intent(in)  ::  solsym(:), fixsym(:), pcesym(:)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer  ::  k, l, photo = 0
      logical  ::  local_flag

      nf = 0
      nsp = 0
      nsr = 0
      npr = 0
      npp = 0
      
rxtnt_scan : &
      do k = 1,rxtcnt
         if( rxtype == photo ) then
            if( k == 2 ) then
               if( rxtsym(k) /= 'hv' ) then
                  call errmes( 'Photo reaction has misplaced or missing "hv" operator@', lout, buff, 1, buff )
               end if
               cycle rxtnt_scan
            else if( k > 2 ) then
               call errmes( 'Photo-reaction can have only one reactant@', lout, buff, 1, buff )
            else if( rxtsym(k) == 'hv' ) then
               call errmes( 'Photo-reaction has misplaced "hv" operator@', lout, buff, 1, buff )
            end if
         else if( rxtsym(k) == 'hv' ) then
            call errmes( ' Photolysis operator "hv" is illegal in a non-photolysis reaction@', lout, buff, 1, buff )
         end if
!-----------------------------------------------------------------------
!    	... Parse out fixed reactants
!-----------------------------------------------------------------------
         if( nfix /= 0 ) then
            do l = 1,nfix
               if( rxtsym(k) == fixsym(l) ) then
                  nf = nf + 1
                  rxttab(1,nf) = l
                  cycle rxtnt_scan
               end if
            end do
         end if
!-----------------------------------------------------------------------
!     	... Parse out solution reactants
!-----------------------------------------------------------------------
         do l = 1,nsol
            if( rxtsym(k) == solsym(l) ) then
               nsr = nsr + 1
               rxttab(2,nsr) = l
               cycle rxtnt_scan
            end if
         end do
!-----------------------------------------------------------------------
!     	... Parse out the pce reactants
!-----------------------------------------------------------------------
         if( npce /= 0 .and. (nsr+nf) /= rxtcnt ) then
            do l = 1,npce
               if( rxtsym(k) == pcesym(l) ) then
                  npr = npr + 1
                  rxttab(4,npr) = l
                  cycle rxtnt_scan
               end if
            end do
         end if
         call errmes( ' Reactant "#" is not in sol, pce, or fixed lists@', lout, rxtsym(k), len_trim(rxtsym(k)), buff )
      end do rxtnt_scan

!-----------------------------------------------------------------------
!     	... Parse out solution products
!-----------------------------------------------------------------------
      local_flag = .false.
      do k = 1,prdcnt
         do l = 1,nsol
            if( prdsym(k) == solsym(l) ) then
               nsp = nsp + 1
               rxttab(3,nsp) = l
               if( coeff_flg ) then
                  pcoeffs(nsp) = pcoeffs(k)
                  if( pcoeffs(k) /= 1.e0 ) then
                     local_flag = .true.
                  end if
               end if
               exit
            end if
         end do
      end do 
      coeff_flg = local_flag
!-----------------------------------------------------------------------
!     	... Parse out the pce products
!-----------------------------------------------------------------------
      do k = 1,prdcnt
         do l = 1,npce
            if( prdsym(k) == pcesym(l) ) then
               npp = npp + 1
               rxttab(5,npp) = l
               exit
            end if
         end do
      end do

      end subroutine mapper

      subroutine outp( rxparms, &
                       nr, &
                       np, &
                       rxtsym, &
                       prdsym, &
                       sym_rate, &
                       irxn, &
                       rate, &
                       loc_rxt_tag, &
                       lout )

      use rxt_mod, only : rxtnt_lim, prd_lim
      use var_mod, only : nq, nfs, solsym, fixsym

      implicit none

!-----------------------------------------------------------------------
!        OUTP OUTPuts a single reaction and rate
!
!        Inputs:
!           nr - number of reactants
!           np - number of products
!           rxparms - vector of "full" product terms (including
!                     multipliers)
!           rxtsym - reactant symbol(s)
!           prdsym - product symbol(s)
!           irxn   - reaction number
!           rate   - vector of reaction rate parameters
!           lout   - logical OUTPut unit number
!        Outputs:
!           NONE
!-----------------------------------------------------------------------

      integer, intent(in) ::      nr, np, irxn, lout
      real, intent(in)    ::      rate(:)
      character(len=*), intent(in) :: rxparms(prd_lim)
      character(len=*), intent(in) :: sym_rate(5)
      character(len=*), intent(in)  :: loc_rxt_tag
      character(len=*), intent(in)  :: rxtsym(rxtnt_lim), prdsym(prd_lim)

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::    i, j, k, kl, length, slen, retcod, line_cnt
      integer  ::    buff_pos, arrow_pos
      real     ::    coeff
      character(len=320) :: buff
      character(len=64)  :: rx_piece

!-----------------------------------------------------------------------
!	... function declarations
!-----------------------------------------------------------------------
      integer :: inclist

      buff = ' '
      j = 1

!-----------------------------------------------------------------------
!        ... Form the reactants
!-----------------------------------------------------------------------
      do i = 1,nr
         length = len_trim( rxtsym(i) )
         buff(j:length+j-1) = rxtsym(i)(:length)
         j = length + j + 1
         if( i == nr ) then
            buff(j:) = '->'
            j = j + 3
         else
            buff(j:) = '+'
            j = j + 2
         end if
      end do
      buff_pos = j ; arrow_pos = j

!-----------------------------------------------------------------------
!        ... Form the products
!-----------------------------------------------------------------------
      line_cnt = 1
      if( np /= 0 ) then
         do i = 1,np
	    if( i /= 1 ) then
	       rx_piece = '+'
	    else
	       rx_piece = ' '
	    end if
	    j = 1
            length = index( rxparms(i), '*' )
            if( length /= 0 ) then
               read(rxparms(i)(:length-1),*,iostat=retcod) coeff
	       if( retcod /= 0 ) then
	          call errmes( ' # is not a valid real number@', &
                               lout, &
                               rxparms(i), &
                               length-1, &
                               buff )
	       end if
               if( coeff /= 1. ) then
		  length = length + 1
		  slen = len_trim( rx_piece )
                  rx_piece(slen+2:slen+length+1) = rxparms(i)(:length-1)
                  j = len_trim( rx_piece ) + 1
	       else
                  j = len_trim( rx_piece ) + 2
               end if
	    else
               j = len_trim( rx_piece ) + 2
            end if
	    kl = inclist( trim(prdsym(i)), solsym, nq )
	    if( kl < 1 .and. nfs > 0 ) then
	       kl = inclist( trim(prdsym(i)), fixsym, nfs )
	    end if
            length = len_trim( prdsym(i) )
	    if( kl > 0 ) then
               rx_piece(j:length+j-1) = prdsym(i)(:length)
	    else
               rx_piece(j:length+j+1) = '{' // prdsym(i)(:length) // '}'
	    end if
            length = len_trim( rx_piece )
	    if( (buff_pos + length) <= 69 ) then
	       buff(buff_pos:) = trim( rx_piece )
	       buff_pos = buff_pos + length + 1
	       if( i == np ) then
		  kl = line_cnt
	          do k = kl,max(4,line_cnt)
                     call write_rxt( buff, sym_rate, rate, loc_rxt_tag, irxn, line_cnt )
	             line_cnt = line_cnt + 1
		  end do
	       end if
	    else
               call write_rxt( buff, sym_rate, rate, loc_rxt_tag, irxn, line_cnt )
	       line_cnt = line_cnt + 1
	       buff(arrow_pos+1:) = trim( rx_piece )
	       buff_pos = len_trim( buff ) + 2
	       if( i == np ) then
		  kl = line_cnt
	          do k = kl,max(4,line_cnt)
                     call write_rxt( buff, sym_rate, rate, loc_rxt_tag, irxn, line_cnt )
	             line_cnt = line_cnt + 1
		  end do
	       end if
	    end if
         end do
      else
         buff(j:) = '(No products)'
	 do k = 1,3
            call write_rxt( buff, sym_rate, rate, loc_rxt_tag, irxn, line_cnt )
	    line_cnt = line_cnt + 1
	 end do
      end if

      end subroutine outp

      subroutine write_rxt( buff, sym_rate, rate, loc_rxt_tag, irxn, line_cnt )
!-----------------------------------------------------------------------
!        ... Print the reaction rate
!-----------------------------------------------------------------------

      use io, only : lout
      use rxt_mod, only : phtcnt

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: line_cnt, irxn
      real, intent(in)    :: rate(:)
      character(len=320), intent(inout) :: buff
      character(len=16), intent(in)     :: sym_rate(:)
      character(len=*), intent(in)      :: loc_rxt_tag

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      logical :: troe_rate

      if( line_cnt <= 3 ) then
         if( sym_rate(1) /= ' ' ) then
	    troe_rate = rate(1) /= 0. .and. rate(3) /= 0.
	    if( line_cnt == 1 ) then
	       if( rate(1) == 0. ) then
                  buff(69:) = ' rate = 0.'
                  write(lout,100) loc_rxt_tag, irxn, buff, irxn+phtcnt
	       else if( .not. troe_rate ) then
                  buff(69:) = ' rate = '
                  write(buff(77:),'(1pe8.2)') rate(1)
                  if( rate(2) /= 0. ) then
                     buff(85:) = '*exp('
                     write(buff(90:),'(f8.0)') rate(2)
                     buff(98:) = '/t)'
                  end if
                  write(lout,100) loc_rxt_tag, irxn, buff, irxn+phtcnt
	       else
                  buff(69:) = ' troe : ko='
                  write(buff(80:),'(1pe8.2)') rate(1)
	          if( rate(2) /= 0. ) then
                     buff(88:) = '*(300/t)**'
                     if ( rate(2)>=0. ) then 
                        write(buff(98:),'(f4.2)') rate(2)
                     else
                        write(buff(98:),'(f5.2)') rate(2)
                     endif
                  end if
                  write(lout,110) loc_rxt_tag, irxn, buff, irxn+phtcnt
               end if
	    else if( troe_rate ) then
	       if( line_cnt == 2 ) then
                  buff(69:) = '        ki='
                  write(buff(80:),'(1pe8.2)') rate(3)
	          if( rate(4) /= 0. ) then
	             if( rate(4) /= 1. ) then
                        buff(88:) = '*(300/t)**'
                        if ( rate(4)>=0. ) then 
                           write(buff(98:),'(f4.2)') rate(4)
                        else
                           write(buff(98:),'(f5.2)') rate(4)
                        endif
	             else
                        buff(88:) = '*(300/t)'
	             end if
                  end if
	       else if( line_cnt == 3 ) then
                  buff(69:) = '         f='
                  write(buff(80:),'(f4.2)') rate(5)
	       end if
               write(lout,120) buff
	    else if( buff /= ' ' ) then
               write(lout,120) buff
            end if
         else
	    if( line_cnt == 1 ) then
               buff(69:) = ' rate = ** User defined **'
               write(lout,100) loc_rxt_tag, irxn, buff, irxn+phtcnt
            else if( buff /= ' ' ) then
               write(lout,120) buff
            end if
         end if
      else if( buff /= ' ' ) then
         write(lout,120) buff
      end if
      buff = ' '
      
!-----------------------------------------------------------------------
!        ... Formats
!-----------------------------------------------------------------------
100   format(2x,a16,1x,'(',i3,')',3x,a100,3x,'(',i3,')')
110   format(2x,a16,1x,'(',i3,')',3x,a101,2x,'(',i3,')')
120   format(27x,a103)

      end subroutine write_rxt

      end module mo_chem
