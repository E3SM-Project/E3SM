
subroutine SYMBOL( iout )
  !-----------------------------------------------------------------------
  !	... Parse out the solution, pce, fixed, column, and group variables
  !-----------------------------------------------------------------------

  use IO
  use VAR_MOD, only : indexh2o, spccnt, indexm, grpsym, grpcof, &
       grpmap, grpcnt, colsym, colmap, spcsym, slvdsym, nslvd, &
       colub, relmap, aliases, var_lim

  implicit none

  !-----------------------------------------------------------------------
  !	... Dummy args
  !-----------------------------------------------------------------------
  character(len=80), intent(out) ::  iout(*)

  !-----------------------------------------------------------------------
  !	... Local variables
  !-----------------------------------------------------------------------
  integer, parameter :: symlen = 16

  integer  ::  retcod, parsw(6), nchar
  integer  ::  spclim(5)
  integer  ::  kpar, i
  integer  ::  symbol_len

  integer, allocatable ::      toklen(:)
  integer  ::  j, k, jl, l, ic, m, count,acount
  integer  ::  jeq, no_tokens, relcnt

  character(len=80), allocatable ::  tokens(:)
  character(len=16) ::  param
  character(len=16)  :: spchdr(6) = (/ 'SOLUTION        ', &
                                       'PCE             ', &
                                       'FIXED           ', &
                                       'GROUPS          ', &
                                       'COL-INT         ', &
                                       'NOT-TRANSPORTED '  /)
  character(len=20)  :: spcend(6) = (/ 'ENDSOLUTION         ', &
                                       'ENDPCE              ', &
                                       'ENDFIXED            ', &
                                       'ENDGROUPS           ', &
                                       'ENDCOL-INT          ', &
                                       'ENDNOT-TRANSPORTED  '  /)
  character(len=16)  ::  upname
  character(len=1)  ::  char

  logical  :: found
  spclim(:) = var_lim
  parsw(:) = 0
  buffh = buff
  call UPCASE( buffh )
  if( buffh /= 'SPECIES' ) then
     call ERRMES( '"Species" card missing; run terminated@', lout, char, 1, buff )
  end if
  kpar = 0
  ALLOCATE( tokens(64) )
  ALLOCATE( toklen(64) )
  acount = 0
  do
     call CARDIN( lin, buff, nchar )
     buffh = buff
     call UPCASE( buffh )
     if( buffh == 'ENDSPECIES' ) then
        DEALLOCATE( tokens )
        DEALLOCATE( toklen )
        exit
     end if
     count = 0
     found = .false.
     do kpar = 1,6
        if( buffh == spchdr(kpar) ) then
           found = .true.
           exit
        end if
     end do
     if( .not. found ) then
        call ERRMES( '# is an invalid species header@', lout, buff(:16), 16, buff )
     else if( parsw(kpar) /= 0 ) then
        call ERRMES( '# species class already declared@', lout, spchdr(kpar), LEN_TRIM(spchdr(kpar)), buff )
     else if( kpar >= 4 .and. MAXVAL( spccnt(1:3) ) == 0 ) then
        call ERRMES( 'There are no sol,pce,fixed species; cannot declare # species class@', lout,  &
             spchdr(kpar), LEN_TRIM(spchdr(kpar)), buff )
     else if( kpar == 2 .and. parsw(1) /= 1 ) then
        call ERRMES( ' Solution species must be specified BEFORE pce species@', lout, buff, 1, ' ' )
     end if
     parsw(kpar) = 1

     if( kpar <= 3 ) then
        !-----------------------------------------------------------------------
        !       ... Read the solution, fixed, and/or "relationship" variables
        !-----------------------------------------------------------------------
        do
           call CARDIN( lin, buff, nchar )
           buffh = buff
           call UPCASE( buffh )
           if( buffh == spcend(kpar) ) then
!!$              if( kpar == 1 .and. count == 0 ) then
!!$                 call ERRMES( 'No # variables declared; terminating@', lout, spchdr(1), &
!!$                      LEN_TRIM(spchdr(1)), &
!!$                      buff )
!!$              end if
              spccnt(kpar) = count
              exit
           end if
           if( kpar == 1  .or. kpar == 3 ) then
              symbol_len = 80
           else if( kpar /= 2 ) then
              symbol_len = symlen
           else
              symbol_len = 2*symlen+1
           end if
           call GETTOKENS(  buff,     nchar,    ',',      symbol_len, &
                tokens,   toklen,   64,       no_tokens )
           if( no_tokens <= 0 ) then
              call ERRMES( ' Species input line in error@', lout, buff, 1, ' ' )
           end if
           if( (count + no_tokens) > spclim(kpar) ) then
              call ERRMES( ' Species count exceeds limit@', lout, buff, 1, buff )
           end if
           do j = 1,no_tokens
              count = count + 1
              acount = acount + 1
              if( kpar == 1  .or. kpar == 3 ) then
                 l = INDEX( tokens(j)(:toklen(j)), '->' )
                 if( l == 0 ) then
                    if( toklen(j) <= symlen ) then
                       spcsym(count,kpar) = tokens(j)(:toklen(j))
                    else
                       call ERRMES( ' Species input line in error@', lout, buff, 1, ' ' )
                    end if
                 else
                    l = l - 1
                    if( l <= symlen ) then
                       spcsym(count,kpar) = tokens(j)(:l)
                    else
                       call ERRMES( ' Species input line in error@', lout, buff, 1, ' ' )
                    end if
                    l = l + 3
                    if( toklen(j) - l + 1 <= 64 ) then
                       aliases(acount) = tokens(j)(l:toklen(j))
                    else
                       call ERRMES( ' Species input line in error@', lout, buff, 1, ' ' )
                    end if
                 end if
              else if( kpar == 2 ) then
                 call GETTOKENS(   tokens(j), toklen(j), '/', &
                      symlen, tokens(19), toklen(19), &
                      2, relcnt )
                 if( relcnt /= 2 ) then
                    call ERRMES( ' Relationship syntax error@', lout, buff, 1, ' ' )
                 end if
                 match_loop :         do m = 1,2
                    do l = 1,spccnt(1)
                       if( tokens(18+m) == spcsym(l,1)) then
                          relmap(count,m) = l
                          cycle match_loop
                       end if
                    end do
                    if( relmap(count,m) == 0 ) then
                       call ERRMES( ' Relation member # not in solution list@', lout, &
                            tokens(18+m), toklen(18+m), &
                            buff )
                    end if
                 end do match_loop
              else
                 spcsym(count,kpar) = tokens(j)(:toklen(j))
              end if
           end do
        end do
        !-----------------------------------------------------------------------
        !       ... Check fixed species for atmospheric total density symbol "M"
        !	    and see if water vapor is declard
        !-----------------------------------------------------------------------
        if( kpar == 3 ) then
           found = .false.
           do i = 1,spccnt(3)
              upname = spcsym(i,3)
              call UPCASE( upname )
              if( upname == 'M' ) then
                 indexm = i
                 found = .true.
                 exit
              end if
           end do
           if( .not. found ) then
              call ERRMES( 'There must be a fixed symbol "m"; run terminated@', lout, char, 1, buff )
           end if
           do i = 1,spccnt(3)
              upname = spcsym(i,3)
              call UPCASE( upname )
              if( upname == 'H2O' ) then
                 indexh2o = i
                 exit
              end if
           end do
        end if
     else if( kpar == 4 ) then
        !-----------------------------------------------------------------------
        !     ... The group section
        !-----------------------------------------------------------------------
        do
           call CARDIN( lin, buff, nchar )
           buffh = buff
           call UPCASE( buffh )
           if( buffh == spcend(kpar) ) then
              spccnt(4) = count
              exit
           end if
           count = count + 1
           if( count > spclim(kpar) ) then
              call ERRMES( 'Group count exceeds limit@', lout, buff, 1, buff )
           end if
           jeq = INDEX( buff(:symlen+1),'=' )
           if( jeq == 0 ) then
              call ERRMES( 'Group name exceeds char limit@', lout, buff, 1, buff )
           else if( jeq == 1 ) then
              call ERRMES( 'Group name missing@', lout, buff, 1, buff )
           else if( jeq == nchar ) then
              call ERRMES( 'Group list missing@', lout, buff, 1, buff )
           end if
           grpsym(count) = buff(:jeq-1)
           jl = jeq + 1
           call GETTOKENS(  buff(jl:nchar),     nchar - jl + 1,    '+', 16, &
                tokens,   toklen,   20,       no_tokens )
           if( no_tokens <= 1 ) then
              call ERRMES( 'One or fewer group members@', lout, buff, 1, buff )
           end if
           l = 0
           do j = 1,no_tokens
              k = INDEX( tokens(j)(:toklen(j)),'*' )
              if( k == 1 ) then
                 call ERRMES( 'No group member multiplier specified@', lout, buff, 1, buff )
              else if( k == toklen(j) ) then
                 call ERRMES( 'No group member specified@', lout, buff, 1, buff )
              end if
              param = tokens(j)(k+1:toklen(j))
              found = .false.
              grp_srch_loop :   do ic = 1,2
                 do m = 1,spccnt(ic)
                    if( param == spcsym(m,ic) ) then
                       found = .true.
                       exit grp_srch_loop
                    end if
                 end do
              end do grp_srch_loop
              if( .not. found ) then
                 call ERRMES( 'Group member # not in sol,pce lists@', lout, param, k, buff )
              end if
              l = l + 1
              if( l > 20 ) then
                 call ERRMES( 'Group member count exceeds limit@', lout, buff, 1, buff )
              end if
              if( k /= 0 ) then
                 param = tokens(j)(:k-1)
                 call RELCON( param, k-1, grpcof(l,count), retcod )
                 if( retcod /= 0 ) then
                    call ERRMES( '# is an invalid group member coefficient@', lout, param, k-1, buff )
                 end if
              end if
              grpmap(l,count) = 1000*ic + m
           end do

           if( l == 0 ) then
              call ERRMES( 'Group has no members@', lout, buff, 1, buff )
           end if
           grpcnt(count) = l
           iout(count) = ' '
           iout(count) = buff(:jeq-1) // ' = '
           do j = 1,no_tokens
              jl = LEN_TRIM(iout(count)) + 1
              if( j == 1 ) then
                 iout(count)(jl+1:) = tokens(j)(:toklen(j))
              else
                 iout(count)(jl:) = ' + ' // tokens(j)(:toklen(j))
              end if
           end do
        end do
     else if( kpar == 5 ) then
        !-----------------------------------------------------------------------
        !     ... The column integrals
        !-----------------------------------------------------------------------
        do
           call CARDIN( lin, buff, nchar )
           buffh = buff
           call UPCASE( buffh )
           if( buffh == spcend(kpar) ) then
              spccnt(5) = count
              exit
           end if
           call GETTOKENS(  buff(:nchar),     nchar,    ',', 32, &
                tokens,   toklen,   20,       no_tokens )
           if( no_tokens <= 0 ) then
!              call ERRMES( 'Column integral list improperly specified@', lout, buff, 1, buff )
           end if
           do j = 1,no_tokens
              jeq = INDEX( tokens(j)(:toklen(j)),'=' )
!!$              if( jeq == 0 ) then
!!$                 call ERRMES( 'No assignment operator in column integral list@', lout, tokens(j), toklen(j), buff )
!!$              else if( jeq == 1 ) then
!!$                 call ERRMES( 'No member in column integral list@', lout, tokens(j), toklen(j), buff )
!!$              else if( jeq == toklen(j) ) then
!!$                 call ERRMES( 'No value in column integral list@', lout, tokens(j), toklen(j), buff )
!!$              else if( jeq > (symlen+1) ) then
!!$                 call ERRMES( 'Column name exceeds char limit@', lout, buff, 1, buff )
!!$              end if
              param = tokens(j)(:jeq-1)
              found = .false.
              colub_srch_loop : do ic = 1,3
                 do m = 1,spccnt(ic)
                    if( param == spcsym(m,ic) ) then
                       found = .true.
                       exit colub_srch_loop
                    end if
                 end do
              end do colub_srch_loop
!!$              if( .not. found ) then
!!$                 call ERRMES( 'Column member # not in sol,pce, or fixed lists@', lout, param, k, buff )
!!$              end if
              count = count + 1
              if( count > spclim(kpar) ) then
                 call ERRMES( 'Column count exceeds limit@', lout, buff, 1, buff )
              end if
              param = tokens(j)(jeq+1:toklen(j))
              call RELCON( param, toklen(j)-jeq, colub(count), retcod )
              if( retcod /= 0 ) then
                 call ERRMES( '# is an invalid upper bndy column density@', lout, param, toklen(j)-jeq, buff )
              end if
              colsym(count) = tokens(j)(:jeq-1)
           end do
        end do
     else if( kpar == 6 ) then
        !-----------------------------------------------------------------------
        !     ... The Short Lived Species
        !-----------------------------------------------------------------------
        do

           call CARDIN( lin, buff, nchar )
           buffh = buff
           call UPCASE( buffh )
           if( buffh == spcend(kpar) ) then
              exit
           end if
           call GETTOKENS(  buff(:nchar),     nchar,    ',', 32, &
                tokens,   toklen,   20,       no_tokens )
!!$           if( no_tokens <= 0 ) then
!!$!              call ERRMES( 'Column integral list improperly specified@', lout, buff, 1, buff )
!!$           end if

           do j=1,no_tokens

              found = .false.
              slvd_srch_loop : do ic = 1,2
                 do m = 1,spccnt(ic)
                    if( tokens(j) == spcsym(m,ic) ) then
                       found = .true.
                       exit slvd_srch_loop
                    end if
                 end do
              end do slvd_srch_loop
              if( .not. found ) then
                 call ERRMES( 'Short Lived Species # not in sol,pce lists@', lout, tokens(j), toklen(j), buff )
              end if

              slvdsym(j) = trim(tokens(j))
           enddo

           nslvd = no_tokens


        end do
     
     end if

  end do

end subroutine SYMBOL
