
      subroutine MAKE_RATE_TAB( rxparm, &
                                rxptab, &
                                rxpcnt )
!-----------------------------------------------------------------------
!        ... Make the code to setup the rate table
!-----------------------------------------------------------------------

      use IO, only : temp_path

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::    rxpcnt
      integer, intent(in) ::    rxptab(*)
      
      real, intent(in)    ::    rxparm(2,*)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer  ::   i, j, k, cnt, pos
      character(len=72) :: line
      logical  ::  lexist
      
      if( rxpcnt == 0 ) then
	 return
      else
!-----------------------------------------------------------------------
!        ... Check for temp dependent rates
!-----------------------------------------------------------------------
         cnt = COUNT( rxparm(2,1:rxpcnt) /= 0. )
         if( cnt == 0 ) then
	    return
         end if
      end if

!-----------------------------------------------------------------------
!        ... First write the table setup routine
!-----------------------------------------------------------------------
      INQUIRE( file = TRIM( temp_path ) // 'rxttab.F', exist = lexist )
      if( lexist ) then
	 call SYSTEM( 'rm ' // TRIM( temp_path ) // 'rxttab.F' )
      end if
      OPEN( unit = 30, file = TRIM( temp_path ) // 'rxttab.F' )

      line = ' '
      write(30,100) line
      line(7:) = 'subroutine RXTTAB( )'
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'implicit none '
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'real      rates'
      write(30,100) line
      line(7:) = 'common / RXTTAB / rates(126,'
      write(line(LEN_TRIM(line)+1:),'(i3,'')'')') cnt
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'real      temp(126)'
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'integer   j'
      write(30,100) line
      line = ' '
      write(30,100) line
      line = '# if defined(EXPHF)'
      write(30,100) line
      line = ' '
      line(7:) = 'real      EXPHF'
      write(30,100) line
      line = 'CDIR$ VFUNCTION EXPHF'
      write(30,100) line
      line = '# endif'
      write(30,100) line
      line = ' '
      write(30,100) line
      
      line(7:) = 'do j = 1,126'
      write(30,100) line
      line = ' '
      line(10:) = 'temp(j) = 1. / (180. + REAL(j-1))'
      write(30,100) line
      line = ' '
      line(7:) = 'end do'
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'do j = 1,126'
      write(30,100) line
      line = ' '
      line(10:) = 'rates(j,   ) ='
      k = 0
      do i = 1,rxpcnt
         if( rxparm(2,i) /= 0.e0 ) then
	    k = k + 1
            write(line(18:20),'(i3)') k
            call R2C( line(25:), rxparm(1,i), 'l' )
	    pos = LEN_TRIM(line)
	    write(30,*) '# if defined(EXPHF)'
            line(pos+1:) = '*EXPHF('
            call R2C( line(LEN_TRIM(line)+1:), rxparm(2,i), 'l' )
            line(LEN_TRIM(line)+1:) = ' * temp(j) )'
            write(30,100) line
	    write(30,*) '# else'
            line(pos+1:) = '*EXP('
            call R2C( line(LEN_TRIM(line)+1:), rxparm(2,i), 'l' )
            line(LEN_TRIM(line)+1:) = ' * temp(j) )'
            write(30,100) line
	    write(30,*) '# endif'
         end if
      end do

      line = ' '
      line(7:) = 'end do'
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'end'
      write(30,100) line
      
      CLOSE(30)
!-----------------------------------------------------------------------
!        ... Finally write the table interpolation routine
!-----------------------------------------------------------------------
      INQUIRE( file = 'setrxt.F', exist = lexist )
      if( lexist ) then
	 call SYSTEM( 'rm setrxt.F' )
      end if
      OPEN( unit = 30, file = 'setrxt.F' )

      line = ' '
      write(30,100) line
      line(7:) = 'subroutine SETRXT( rate,'
      write(30,100) line
      line(6:) = '$                   temp )'
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'implicit none '
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'real      rate(PLNPLV,RXNCNT)'
      write(30,100) line
      line(7:) = 'real      temp(PLNPLV)'
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'real      rates'
      write(30,100) line
      line(7:) = 'common / RXTTAB / rates(126,'
      write(line(LEN_TRIM(line)+1:),'(i3,'')'')') cnt
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'integer   i, ip1, j'
      write(30,100) line
      line(7:) = 'real      del_temp'
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'do j = 1,PLNPLV'
      write(30,100) line
      
      if( cnt /= rxpcnt ) then
         line = ' '
         line(10:) = 'rate(j,   ) ='
         do i = 1,rxpcnt
            if( rxparm(2,i) == 0.e0 ) then
               write(line(17:19),'(i3)') rxptab(i)
               call R2C( line(24:), rxparm(1,i), 'l' )
               write(30,100) line
            end if
         end do
      end if
      line = ' '
      line(7:) = 'end do'
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'do j = 1,PLNPLV'
      write(30,100) line
      line = ' '
      line(10:) = 'i = INT( temp(j) ) - 179'
      write(30,100) line
      line(10:) = 'i = MAX( 1,MIN( 125,i) )'
      write(30,100) line
      line(10:) = 'ip1 = i + 1'
      write(30,100) line
      line(10:) = 'del_temp = temp(j) - AINT(temp(j))'
      write(30,100) line
      line = ' '
      k = 0
      do i = 1,rxpcnt
         if( rxparm(2,i) /= 0.e0 ) then
            line = ' '
            line(10:) = 'rate(j,'
	    k = k + 1
            write(line(LEN_TRIM(line)+1:),'(i3,'') ='')') rxptab(i)
	    line(LEN_TRIM(line)+2:) = 'rates(i,'
            write(line(LEN_TRIM(line)+1:),'(i3,'')'')') k
            write(30,100) line
	    j = INDEX( line,'=' ) + 2
            line(6:) = '$'
            line(j:) = '+ del_temp * (rates(ip1,'
            write(line(LEN_TRIM(line)+1:),'(i3,'')'')') k
	    line(LEN_TRIM(line)+2:) = '- rates(i,'
            write(line(LEN_TRIM(line)+1:),'(i3,''))'')') k
            write(30,100) line
         end if
      end do

      line = ' '
      line(7:) = 'end do'
      write(30,100) line
      line = ' '
      write(30,100) line
      line(7:) = 'end'
      write(30,100) line
      CLOSE(30)
      
100   format(a72)

      end subroutine MAKE_RATE_TAB
