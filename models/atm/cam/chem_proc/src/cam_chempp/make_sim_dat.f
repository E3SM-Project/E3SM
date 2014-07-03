    
      subroutine make_sim_dat( model, sparse )
!-------------------------------------------------------------------
!	... write the simulation data routine; only for CAM
!-------------------------------------------------------------------

      use io,      only : temp_path
      use sp_mods, only : sparsity
      use var_mod, only : clscnt, clsmap, permute, new_nq, new_solsym
      use var_mod, only : nq, newind, mass, c_mass, temp_mass
      use var_mod, only : nfs, fixsym
      use var_mod, only : nslvd, slvdsym
      use rxt_mod, only : cls_rxt_cnt, rxntot
      use rxt_mod, only : rxt_has_tag, rxt_tag
      use rxt_mod, only : phtcnt, pht_alias, pht_alias_mult
      use rxt_mod, only : usrcnt, usrmap, frc_from_dataset

      implicit none

!-------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      character(len=16), intent(in) :: model
      type(sparsity), intent(in)   :: sparse(2)

!-------------------------------------------------------------------
!	... local variables
!-------------------------------------------------------------------
      integer, parameter :: max_len= 132

      integer  ::  i, l, m, m1, n, n1
      integer  ::  lpos
      integer  ::  lstrt
      integer, allocatable   :: ndx(:)
      character(len=max_len) :: line
      character(len=64)      :: frmt
      character(len=24)      :: number
      character(len=12)      :: num12
      character(len=16)      :: rxt_string
      character(len=16)       :: wrk_chr(5)
      logical  ::  flush
      logical  ::  lexist
      integer  ::  numlen

      inquire( file = trim( temp_path ) // 'mo_sim_dat.F', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'mo_sim_dat.F' )
      end if
      open( unit = 30, file = trim( temp_path ) // 'mo_sim_dat.F' )

      line = ' '
      write(30,100) trim(line)
      line(7:) = 'module mo_sim_dat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'private'
      write(30,100) trim(line)
      line(7:) = 'public :: set_sim_dat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'contains'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'subroutine set_sim_dat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:)   = 'use chem_mods,   only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass, crb_mass'
      write(30,100) trim(line)
      if( clscnt(4) > 0 ) then
        line(7:) = 'use chem_mods,   only : diag_map'
        write(30,100) trim(line)
      endif
      line(7:)   = 'use chem_mods,   only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map'
      write(30,100) trim(line)
      line(7:)   = 'use chem_mods,   only : pht_alias_lst, pht_alias_mult'
      write(30,100) trim(line)
      line(7:)   = 'use chem_mods,   only : extfrc_lst, inv_lst, slvd_lst'
      write(30,100) trim(line)
      line(7:)   = 'use abortutils,  only : endrun'
      write(30,100) trim(line)
      line(7:)   = 'use mo_tracname, only : solsym'
      write(30,100) trim(line)
      line(7:)   = 'use chem_mods,   only : frc_from_dataset'
      write(30,100) trim(line)
      line(7:)   = 'use shr_kind_mod,only : r8 => shr_kind_r8'
      write(30,100) trim(line)
      line(7:)   = 'use cam_logfile, only : iulog'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'implicit none '
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------'
      write(30,100) trim(line)
      line = '!      ... local variables'
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------'
      write(30,100) trim(line)
      line = '      integer :: ios'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
!-------------------------------------------------------------------
!	... set the simulation chemical mechanism data
!	    class species count
!-------------------------------------------------------------------
      line = '      clscnt(:) = (/'
      write(line(len_trim(line)+2:),'(5(I6,a))') clscnt(1),',',clscnt(2),',',clscnt(3),',',clscnt(4),',',clscnt(5), ' /)'
      write(30,'(a)') trim(line)
      line = ' '
      write(30,100) trim(line)

!-------------------------------------------------------------------
!	... class reaction count
!-------------------------------------------------------------------
      do i = 1,5
         if( clscnt(i) > 0 ) then
            line = '      cls_rxt_cnt(:,'
            write(line(len_trim(line)+1:),'(i1,") = (/")') i
            m = len_trim(line) + 2
            write(line(m:),'(4(I6,a))') cls_rxt_cnt(1,i),',',cls_rxt_cnt(2,i),',',cls_rxt_cnt(3,i),',',cls_rxt_cnt(4,i),' /)'
            write(30,'(a)') trim(line)
         end if
      end do

!-------------------------------------------------------------------
!	... species symbols
!-------------------------------------------------------------------
      line = ' '
      write(30,100) trim(line)
      write(line,'("      solsym(:",i3,") = (/")') new_nq
      m = len_trim(line) + 2
      do n = 1,new_nq,5
         n1 = min( n+4,new_nq )
         if( n1 /= new_nq ) then
            write(line(m:),'(5("''",a16,"'',")," &")') new_solsym(n:n1)
         else
            if( n1 > n ) then
               write(frmt,'("(",i1)') n1 - n
               frmt(len_trim(frmt)+1:) = '("''",a16,"'',"),"''",a16,"'' /)")'
            else
               frmt = '("''",a16,"'' /)")'
            end if
            write(line(m:),trim(frmt)) new_solsym(n:n1)
         end if
         write(30,'(a)') trim(line)
         line = ' '
      end do

!-------------------------------------------------------------------
!	... species mass
!-------------------------------------------------------------------
      if( nq > 0 ) then
         line = ' '
         write(30,100) trim(line)
         temp_mass(:) = 0.
         do n = 1,nq
            if( newind(n) /= 0 ) then
               temp_mass(newind(n)) = mass(n)
            end if
         end do
         line = '      adv_mass(:'
         write(line(len_trim(line)+1:),'(i3,") = (/")') new_nq
         m     = len_trim(line) + 2
         lstrt = m
         do n = 1,new_nq
            number = ' '
            if ( temp_mass(n) > 1. ) then
               write(num12,'(f12.6)') temp_mass(n)
            else
               write(num12,'(g12.6)') temp_mass(n)
            endif
            numlen = len_trim(num12)
            number(12-numlen+1:12) = num12(1:numlen)
            lpos   = scan( number, '0123456789', back=.true. ) + 1
            if( n < new_nq ) then
               if( mod(n,5) /= 0 ) then
                  number(lpos:) = '_r8,'
                  flush = .false.
               else
                  number(lpos:) = '_r8, &'
                  flush = .true.
               end if
            else
               number(lpos:) = '_r8 /)'
               flush = .true.
            end if
            line(m:) = trim( number )
            if( .not. flush ) then
               m = len_trim(line) + 2
            else
               write(30,'(a)') trim(line)
               line = ' '
               m = lstrt
            end if
         end do
      end if

!-------------------------------------------------------------------
!	... species carbon mass
!-------------------------------------------------------------------
      if( nq > 0 ) then
         line = ' '
         write(30,100) trim(line)
         temp_mass(:) = 0.
         do n = 1,nq
            if( newind(n) /= 0 ) then
               temp_mass(newind(n)) = c_mass(n)
            end if
         end do
         line = '      crb_mass(:'
         write(line(len_trim(line)+1:),'(i3,") = (/")') new_nq
         m     = len_trim(line) + 2
         lstrt = m
         do n = 1,new_nq
            number = ' '
            write(num12,'(f12.6)') temp_mass(n)
            numlen = len_trim(num12)
            number(12-numlen+1:12) = num12(1:numlen)
            lpos   = scan( number, '0123456789', back=.true. ) + 1
            if( n < new_nq ) then
               if( mod(n,5) /= 0 ) then
                  number(lpos:) = '_r8,'
                  flush = .false.
               else
                  number(lpos:) = '_r8, &'
                  flush = .true.
               end if
            else
               number(lpos:) = '_r8 /)'
               flush = .true.
            end if
            line(m:) = trim( number )
            if( .not. flush ) then
               m = len_trim(line) + 2
            else
               write(30,'(a)') trim(line)
               line = ' '
               m = lstrt
            end if
         end do
      end if
!-------------------------------------------------------------------
! 	... fixed species masses
!-------------------------------------------------------------------
      if( nfs > 0 ) then
         line = ' '
         write(30,100) trim(line)
         line = '      fix_mass(:'
         write(line(len_trim(line)+1:),'(i3,") = (/")') nfs
         m     = len_trim(line) + 2
         lstrt = m
         do n = 1,nfs
            number = ' '
            write(number,'(g15.9)') mass(n+new_nq)
            number = adjustl( number )
            lpos   = scan( number, '0123456789', back=.true. ) + 1
            if( n < nfs ) then
               if( mod(n,5) /= 0 ) then
                  number(lpos:) = '_r8,'
                  flush = .false.
               else
                  number(lpos:) = '_r8, &'
                  flush = .true.
               end if
            else
               number(lpos:) = '_r8 /)'
               flush = .true.
            end if
            line(m:) = trim( number )
            if( .not. flush ) then
               m = len_trim(line) + 2
            else
               write(30,'(a)') trim(line)
               line = ' '
               m = lstrt
            end if
         end do
      end if

!-------------------------------------------------------------------
!	... class map
!-------------------------------------------------------------------
      line = ' '
      write(30,100) trim(line)
      do i = 1,5
         if( clscnt(i) > 0 ) then
            write(line,'("      clsmap(:",i3,",",i1,") = (/")') clscnt(i),i
            m = len_trim(line) + 2
            do n = 1,clscnt(i),10
               n1 = min( n+9,clscnt(i) )
               if( n1 /= clscnt(i) ) then
                  write(line(m:),'(10(i4,",")," &")') clsmap(n:n1,i,2)
               else
                  if( n1 > n ) then
                     write(frmt,'("(",i1)') n1 - n
                     frmt(len_trim(frmt)+1:) = '(i4,","),i4," /)")'
                  else
                     frmt = '(i4," /)")'
                  end if
                  write(line(m:),trim(frmt)) clsmap(n:n1,i,2)
               end if
               write(30,'(a)') trim(line)
               line = ' '
            end do
         end if
      end do

!-------------------------------------------------------------------
!	... class permutation map
!-------------------------------------------------------------------
      line = ' '
      write(30,100) trim(line)
      do i = 2,5
         if( clscnt(i) > 0 ) then
            write(line,'("      permute(:",i3,",",i1,") = (/")') clscnt(i),i
            m = len_trim(line) + 2
            do n = 1,clscnt(i),10
               n1 = min( n+9,clscnt(i) )
               if( n1 /= clscnt(i) ) then
                  write(line(m:),'(10(i4,",")," &")') permute(n:n1,i)
               else
                  if( n1 > n ) then
                     write(frmt,'("(",i1)') n1 - n
                     frmt(len_trim(frmt)+1:) = '(i4,","),i4," /)")'
                  else
                     frmt = '(i4," /)")'
                  end if
                  write(line(m:),trim(frmt)) permute(n:n1,i)
               end if
               write(30,'(a)') trim(line)
               line = ' '
            end do
         end if
      end do

!-------------------------------------------------------------------
!	... class diagonal indicies
!-------------------------------------------------------------------
      line = ' '
      write(30,100) trim(line)
      do i = 4,4
         if( clscnt(i) > 0 ) then
            write(line,'("      diag_map(:",i3,") = (/")') clscnt(i)
            m = len_trim(line) + 2
            do n = 1,clscnt(i),10
               n1 = min( n+9,clscnt(i) )
               if( n1 /= clscnt(i) ) then
                  write(line(m:),'(10(i4,",")," &")') sparse(i-3)%diag_map(n:n1)
               else
                  if( n1 > n ) then
                     write(frmt,'("(",i1)') n1 - n
                     frmt(len_trim(frmt)+1:) = '(i4,","),i4," /)")'
                  else
                     frmt = '(i4," /)")'
                  end if
                  write(line(m:),trim(frmt)) sparse(i-3)%diag_map(n:n1)
               end if
               write(30,'(a)') trim(line)
               line = ' '
            end do
         end if
      end do

!-----------------------------------------------------------------------
!        ... Write the ext frcing species
!-----------------------------------------------------------------------
   if( usrcnt > 0 ) then
      line = ' '
      write(30,100) trim(line)
      write(line,'("      extfrc_lst(:",i3,") = (/")') usrcnt
      m = len_trim(line) + 2
      do n = 1,usrcnt,5
	 wrk_chr(:) = ' '
         n1 = min( n+4,usrcnt )
	 do i = 1,n1-n+1 !!n,n1
	    wrk_chr(i) = new_solsym(usrmap(i+n-1))
         end do        
         if( n1 /= usrcnt ) then
            write(line(m:),'(5("''",a16,"'',")," &")') wrk_chr(1:n1-n+1)
         else
            if( n1 > n ) then
               write(frmt,'("(",i1)') n1 - n
               frmt(len_trim(frmt)+1:) = '("''",a16,"'',"),"''",a16,"'' /)")'
            else
               frmt = '("''",a16,"'' /)")'
            end if
            write(line(m:),trim(frmt)) wrk_chr(1:n1-n+1)
         end if
         write(30,'(a)') trim(line)
         line = ' '
      end do
      
      ! frc_from_dataset
      line = ' '
      write(30,100) trim(line)
      write(line,'("      frc_from_dataset(:",i3,") = (/")') usrcnt
      m1 = len_trim(line) + 2
      do n = 1,usrcnt,5
         n1 = min( n+4,usrcnt )
         m = m1
         do l = n,n1
            if( l /= usrcnt ) then
               if( l /= n1 ) then
                  if( frc_from_dataset(l) ) then
                     write(line(m:),'(".true.,")')
                  else
                     write(line(m:),'(".false.,")')
                  end if
               else
                  if( frc_from_dataset(l) ) then
                     write(line(m:),'(".true., &")')
                  else
                     write(line(m:),'(".false., &")')
                  end if
               end if
            else
               if( frc_from_dataset(l) ) then
                  write(line(m:),'(".true. /)")')
               else
                  write(line(m:),'(".false. /)")')
               end if
            end if
            m = len_trim(line) + 2
         end do
         write(30,'(a)') trim(line)
         line = ' '
      end do
   end if

!-------------------------------------------------------------------
!	... fixed species
!-------------------------------------------------------------------
      if( nfs > 0 ) then
         line = ' '
         write(30,100) trim(line)
         write(line,'("      inv_lst(:",i3,") = (/")') nfs
         m1 = len_trim(line) + 2
         do n = 1,nfs,5
            n1 = min( n+4,nfs )
            m = m1
            do l = n,n1
               if( l /= nfs ) then
                  if( l /= n1 ) then
                     write(line(m:),'("''",a16,"'',")') fixsym(l)
                  else
                     write(line(m:),'("''",a16,"'', &")') fixsym(l)
                  end if
               else
                  write(line(m:),'("''",a16,"'' /)")') fixsym(l)
               end if
               m = len_trim(line) + 2
            end do
            write(30,'(a)') trim(line)
            line = ' '
         end do
      end if

!-------------------------------------------------------------------
!	... short lived species
!-------------------------------------------------------------------
      if( nslvd > 0 ) then
         line = ' '
         write(30,100) trim(line)
         write(line,'("      slvd_lst(:",i3,") = (/")') nslvd
         m1 = len_trim(line) + 2
         do n = 1,nslvd,5
            n1 = min( n+4,nslvd )
            m = m1
            do l = n,n1
               if( l /= nslvd ) then
                  if( l /= n1 ) then
                     write(line(m:),'("''",a16,"'',")') slvdsym(l)
                  else
                     write(line(m:),'("''",a16,"'', &")') slvdsym(l)
                  end if
               else
                  write(line(m:),'("''",a16,"'' /)")') slvdsym(l)
               end if
               m = len_trim(line) + 2
            end do
            write(30,'(a)') trim(line)
            line = ' '
         end do
      end if

!-------------------------------------------------------------------
!	... reaction tags
!-------------------------------------------------------------------
      i = count( rxt_has_tag(:rxntot) )
      if( i > 0 ) then
         allocate( ndx(i) )
         l = 0
         do m = 1,rxntot
            if( rxt_has_tag(m) ) then
               l = l + 1
               ndx(l) = m
            end if
         end do 
         line = ' '
         write(30,100) trim(line)
!!$         write(line,'("      rxt_tag_cnt = ",i4)') i
!!$         write(30,100) trim(line)
!!$         line = ' '
         line(7:) = 'if( allocated( rxt_tag_lst ) ) then'
         write(30,100) trim(line)
         line(7:) = '   deallocate( rxt_tag_lst )'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line(7:) = 'allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )'
         write(30,100) trim(line)
         line(7:) = 'if( ios /= 0 ) then'
         write(30,100) trim(line)
         line = ' '
         line(10:) = 'write(iulog,*) ''set_sim_dat: failed to allocate rxt_tag_lst; error = '',ios'
         write(30,100) trim(line)
         line(10:) = 'call endrun'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line(7:) = 'if( allocated( rxt_tag_map ) ) then'
         write(30,100) trim(line)
         line(7:) = '   deallocate( rxt_tag_map )'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line(7:) = 'allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )'
         write(30,100) trim(line)
         line(7:) = 'if( ios /= 0 ) then'
         write(30,100) trim(line)
         line = ' '
         line(10:) = 'write(iulog,*) ''set_sim_dat: failed to allocate rxt_tag_map; error = '',ios'
         write(30,100) trim(line)
         line(10:) = 'call endrun'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line = '      rxt_tag_lst(:rxt_tag_cnt) = (/ '
         m1 = len_trim(line) + 2
         do n = 1,i,4
            n1 = min( n+3,i )
            m = m1
            do l = n,n1
               rxt_string = rxt_tag(ndx(l))
               lpos = index( rxt_string, ',cph' )
               if( lpos > 0 ) then
                  rxt_string = trim( rxt_string(:lpos-1) )
               end if
               if( l /= i ) then
                  if( l /= n1 ) then
                     write(line(m:),'("''",a16,"'',")') rxt_string
                  else
                     write(line(m:),'("''",a16,"'', &")') rxt_string
                  end if
               else
                  write(line(m:),'("''",a16,"'' /)")') rxt_string
               end if
               m = len_trim(line) + 2
            end do
            write(30,'(a)') trim(line)
            line = ' '
         end do

         line = '      rxt_tag_map(:rxt_tag_cnt) = (/'
         m = len_trim(line) + 2
         do n = 1,i,10
            n1 = min( n+9,i )
            if( n1 /= i ) then
               write(line(m:),'(10(i4,",")," &")') ndx(n:n1)
            else
               if( n1 > n ) then
                  write(frmt,'("(",i1)') n1 - n
                  frmt(len_trim(frmt)+1:) = '(i4,","),i4," /)")'
               else
                  frmt = '(i4," /)")'
               end if
               write(line(m:),trim(frmt)) ndx(n:n1)
            end if
            write(30,'(a)') trim(line)
            line = ' '
         end do
         deallocate( ndx )
      end if

!-------------------------------------------------------------------
!	... photoreactions alias
!-------------------------------------------------------------------
      if( phtcnt > 0 ) then
         line = ' '
         line(7:) = 'if( allocated( pht_alias_lst ) ) then'
         write(30,100) trim(line)
         line(7:) = '   deallocate( pht_alias_lst )'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line(7:) = 'allocate( pht_alias_lst(phtcnt,2),stat=ios )'
         write(30,100) trim(line)
         line(7:) = 'if( ios /= 0 ) then'
         write(30,100) trim(line)
         line = ' '
         line(10:) = 'write(iulog,*) ''set_sim_dat: failed to allocate pht_alias_lst; error = '',ios'
         write(30,100) trim(line)
         line(10:) = 'call endrun'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line = ' '
         line(7:) = 'if( allocated( pht_alias_mult ) ) then'
         write(30,100) trim(line)
         line(7:) = '   deallocate( pht_alias_mult )'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line(7:) = 'allocate( pht_alias_mult(phtcnt,2),stat=ios )'
         write(30,100) trim(line)
         line(7:) = 'if( ios /= 0 ) then'
         write(30,100) trim(line)
         line = ' '
         line(10:) = 'write(iulog,*) ''set_sim_dat: failed to allocate pht_alias_mult; error = '',ios'
         write(30,100) trim(line)
         line(10:) = 'call endrun'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         do i = 1,2
            if( i == 1 ) then
               line = '      pht_alias_lst(:,1) = (/ '
            else
               line = '      pht_alias_lst(:,2) = (/ '
            end if
            m1 = len_trim(line) + 2
            do n = 1,phtcnt,4
               n1 = min( n+3,phtcnt )
               m = m1
               do l = n,n1
                  rxt_string = pht_alias(l,i)
                  if( l /= phtcnt ) then
                     if( l /= n1 ) then
                        write(line(m:),'("''",a16,"'',")') rxt_string
                     else
                        write(line(m:),'("''",a16,"'', &")') rxt_string
                     end if
                  else
                     write(line(m:),'("''",a16,"'' /)")') rxt_string
                  end if
                  m = len_trim(line) + 2
               end do
               write(30,'(a)') trim(line)
               line = ' '
            end do
         end do

         do i = 1,2
            if( i == 1 ) then
               line = '      pht_alias_mult(:,1) = (/ '
            else
               line = '      pht_alias_mult(:,2) = (/ '
            end if
            m = len_trim(line) + 2
            do n = 1,phtcnt
               number = ' '
               write(number,'(a)') pht_alias_mult(n,i)
               number = adjustl( number )
               lpos   = scan( number, '0123456789', back=.true. )
               if( lpos == 1 ) then
                  lpos = lpos + 1
                  number(lpos:lpos) = '.'
               end if
               lpos = lpos + 1
               if( n < phtcnt ) then
                  if( mod(n,5) /= 0 ) then
                     number(lpos:) = '_r8,'
                     flush = .false.
                  else
                     number(lpos:) = '_r8, &'
                     flush = .true.
                  end if
               else
                  number(lpos:) = '_r8 /)'
                  flush = .true.
               end if
               line(m:) = trim( number )
               if( .not. flush ) then
                  m = len_trim(line) + 2
               else
                  write(30,'(a)') trim(line)
                  line = ' '
                  m = lstrt
               end if
            end do
         end do
      end if

      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end subroutine set_sim_dat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end module mo_sim_dat'
      write(30,100) trim(line)

100   format(a)

      close( unit = 30 )
      end subroutine make_sim_dat
