
      subroutine OUTP( rxparms, &
                       nr, &
                       np, &
                       rxtsym, &
                       prdsym, &
                       sym_rate, &
                       irxn, &
                       rate, &
                       loc_rxt_alias, &
                       lout )

      use RXT_MOD, only : rxtnt_lim, prd_lim

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
      character(len=16), intent(in) :: rxparms(prd_lim)
      character(len=16), intent(in) :: sym_rate(5)
      character(len=16), intent(in)  :: loc_rxt_alias
      character(len=16), intent(in)  :: rxtsym(rxtnt_lim), prdsym(prd_lim)

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::    i, j, k, kl, length, retcod, line_cnt
      integer  ::    buff_pos, arrow_pos
      real     ::    coeff
      character(len=320) :: buff
      character(len=64)  :: rx_piece

      buff = ' '
      j = 1

!-----------------------------------------------------------------------
!        ... Form the reactants
!-----------------------------------------------------------------------
      do i = 1,nr
         length = LEN_TRIM( rxtsym(i) )
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
      buff_pos = j ; arrow_pos = j - 1

!-----------------------------------------------------------------------
!        ... Form the products
!-----------------------------------------------------------------------
      line_cnt = 1
      if( np /= 0 ) then
         do i = 1,np
	    rx_piece = ' '
	    j = 1
            length = INDEX( rxparms(i), '*' )
            if( length /= 0 ) then
               read(rxparms(i)(:length-1),*,iostat=retcod) coeff
	       if( retcod /= 0 ) then
	          call ERRMES( ' # is not a valid real number@', &
                               lout, &
                               rxparms(i), &
                               length-1, &
                               buff )
	       end if
               if( coeff /= 1. ) then
		  length = length + 1
                  rx_piece(:length) = rxparms(i)(:length-1) // '*'
                  j = length
               end if
            end if
            length = LEN_TRIM( prdsym(i) )
            rx_piece(j:length+j-1) = prdsym(i)(:length)
            length = LEN_TRIM( rx_piece )
	    if( (buff_pos + length) <= 69 ) then
	       buff(buff_pos:) = TRIM( rx_piece )
	       buff_pos = buff_pos + length + 1
	       if( i /= np ) then
                  buff(buff_pos:buff_pos) = '+'
		  buff_pos = buff_pos + 2
	       else
		  kl = line_cnt
	          do k = kl,3
                     call WRITE_RXT( buff, sym_rate, rate, loc_rxt_alias, irxn, line_cnt )
	             line_cnt = line_cnt + 1
		  end do
	       end if
	    else
               call WRITE_RXT( buff, sym_rate, rate, loc_rxt_alias, irxn, line_cnt )
	       line_cnt = line_cnt + 1
	       if( i /= np ) then
	          buff(arrow_pos:arrow_pos) = '+'
		  buff_pos = arrow_pos + 2
	       else
		  kl = line_cnt
	          do k = kl,3
                     call WRITE_RXT( buff, sym_rate, rate, loc_rxt_alias, irxn, line_cnt )
	             line_cnt = line_cnt + 1
		  end do
	       end if
	    end if
         end do
      else
         buff(j:) = '(No products)'
	 do k = 1,3
            call WRITE_RXT( buff, sym_rate, rate, loc_rxt_alias, irxn, line_cnt )
	    line_cnt = line_cnt + 1
	 end do
      end if

      end subroutine OUTP

      subroutine WRITE_RXT( buff, sym_rate, rate, loc_rxt_alias, irxn, line_cnt )
!-----------------------------------------------------------------------
!        ... Print the reaction rate
!-----------------------------------------------------------------------

      use IO, only : lout
      use RXT_MOD, only : phtcnt

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: line_cnt, irxn
      real, intent(in)    :: rate(:)
      character(len=320), intent(inout) :: buff
      character(len=16), intent(in)     :: sym_rate(:)
      character(len=16), intent(in)      :: loc_rxt_alias

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
                  write(lout,100) loc_rxt_alias, irxn, buff, irxn+phtcnt
	       else if( .not. troe_rate ) then
                  buff(69:) = ' rate = '
                  write(buff(77:),'(1pe8.2)') rate(1)
                  if( rate(2) /= 0. ) then
                     buff(85:) = '*EXP('
                     write(buff(90:),'(f8.0)') rate(2)
                     buff(98:) = '/t)'
                  end if
                  write(lout,100) loc_rxt_alias, irxn, buff, irxn+phtcnt
	       else
                  buff(69:) = ' troe : ko='
                  write(buff(80:),'(1pe8.2)') rate(1)
	          if( rate(2) /= 0. ) then
                     buff(88:) = '*(300/t)**'
                     write(buff(98:),'(f4.2)') rate(2)
                  end if
                  write(lout,110) loc_rxt_alias, irxn, buff, irxn+phtcnt
               end if
	    else if( troe_rate ) then
	       if( line_cnt == 2 ) then
                  buff(69:) = '        ki='
                  write(buff(80:),'(1pe8.2)') rate(3)
	          if( rate(4) /= 0. ) then
	             if( rate(4) /= 1. ) then
                        buff(88:) = '*(300/t)**'
                        write(buff(98:),'(f4.2)') rate(4)
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
               write(lout,100) loc_rxt_alias, irxn, buff, irxn+phtcnt
            end if
         end if
      else if( buff /= ' ' ) then
         write(lout,120) buff
      end if
      buff = ' '
      
!-----------------------------------------------------------------------
!        ... Formats
!-----------------------------------------------------------------------
100   format(2x,a8,1x,'(',i3,')',3x,a100,3x,'(',i3,')')
110   format(2x,a8,1x,'(',i3,')',3x,a101,2x,'(',i3,')')
120   format(19x,a101)

      end subroutine WRITE_RXT
