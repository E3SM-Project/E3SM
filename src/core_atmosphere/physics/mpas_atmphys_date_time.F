! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_date_time
 use mpas_kind_types

 implicit none
 private
 public:: get_julgmt,             &
          monthly_interp_to_date, &
          monthly_min_max

 character(len=StrKIND),public:: current_date


! MPAS utility module for time management.
! Laura D. Fowler (send comments to laura@ucar.edu).
! 2013-05-01.
!
! subroutines called in mpas_atmphys_date_time:
! ---------------------------------------------
! get_julgmt            : calculates current julian day.
! split_date_char       : used to extract actual date from input string.
! monthly_interp_to_date: interpolates monthly-mean data to current julian day.
! monthly_min_max       : looks for min and max values from monthly-mean data set (greenfrac,...)
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * in subroutine monthly_interp_to_date, change the length of variables day15 and mon from
!   StrKIND to 2 to input correctly the reference date to subroutine get_julgmt_date.
! * in subroutines get_julgmt_date and split_date_char, changed the declaration of date_str
!   from StrKIND to *.
!   Laura D. Fowler (laura@ucar.edu) / 2013-10-18.


 contains


!=================================================================================================================
 subroutine get_julgmt(date_str,julyr,julday,gmt)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: date_str

!output arguments:
 integer,intent(out):: julyr
 integer,intent(out):: julday

 real(kind=RKIND),intent(out) :: gmt

!local variables:
 integer:: ny , nm , nd , nh , ni , ns , nt
 integer:: my1, my2, my3, monss
 integer,dimension(12):: mmd
 data mmd /31,28,31,30,31,30,31,31,30,31,30,31/

!-----------------------------------------------------------------------------------------------------------------

 call split_date_char(date_str,ny,nm,nd,nh,ni,ns,nt)

 gmt = nh + float(ni)/60. + float(ns)/3600.
 my1 = mod(ny,4)
 my2 = mod(ny,100)
 my3 = mod(ny,400)
 if(my1.eq.0.and.my2.ne.0.or.my3.eq.0)mmd(2)=29
 julday=nd
 julyr=ny
 do monss=1,nm-1
    julday=julday+mmd(monss)
 enddo

 end subroutine get_julgmt

!=================================================================================================================
 subroutine split_date_char(date,century_year,month,day,hour,minute,second,ten_thousandth)
!=================================================================================================================
   
!input arguments:
 character(len=*),intent(in):: date

!output arguments:
 integer,intent(out):: century_year,month,day,hour,minute,second,ten_thousandth

!-----------------------------------------------------------------------------------------------------------------

 read(date,fmt='(    I4)') century_year
 read(date,fmt='( 5X,I2)') month
 read(date,fmt='( 8X,I2)') day
 read(date,fmt='(11X,I2)') hour
 read(date,fmt='(14X,I2)') minute
 read(date,fmt='(17X,I2)') second
 read(date,fmt='(20X,I4)') ten_thousandth
   
 end subroutine split_date_char

!=================================================================================================================
 subroutine monthly_interp_to_date(npoints,date_str,field_in,field_out)
!=================================================================================================================

!input arguments:
 character(len=StrKIND),intent(in):: date_str
 integer,intent(in):: npoints
 real(kind=RKIND),intent(in) ,dimension(12,npoints):: field_in

!output arguments:
 real(kind=RKIND),intent(out),dimension(npoints):: field_out

!local variables:
 character(len=2):: day15,mon

 integer:: l,n
 integer:: julyr,julday,int_month,month1,month2
 integer:: target_julyr,target_julday,target_date
 integer,dimension(0:13):: middle

 real(kind=RKIND):: gmt

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine monthly_interp_to_date:'
!write(0,*) '--- current_date  = ',trim(date_str)

 write(day15,fmt='(I2.2)') 15
 do l = 1 , 12
    write(mon,fmt='(I2.2)') l
    call get_julgmt(date_str(1:4)//'-'//mon//'-'//day15//'_'//'00:00:00.0000' , &
                     julyr,julday,gmt)
    middle(l) = julyr*1000 + julday
 enddo

 l = 0
 middle(l) = middle( 1) - 31

 l = 13
 middle(l) = middle(12) + 31

 call get_julgmt(date_str,target_julyr,target_julday,gmt)
 target_date = target_julyr * 1000 + target_julday
!write(0,*) '--- target_julday =',target_julday
!write(0,*) '--- target_date   =',target_date

 find_month : do l = 0 , 12
    if((middle(l) .lt. target_date) .and. (middle(l+1) .ge. target_date)) then
       do n = 1, npoints
          int_month = l
          if((int_month .eq. 0) .or. (int_month .eq. 12)) then
             month1 = 12
             month2 =  1
          else
             month1 = int_month
             month2 = month1 + 1
          endif
          if(n == 1) then
!            write(0,*) '--- month1 =',month1
!            write(0,*) '--- month2 =',month2
          endif
          field_out(n) = ( field_in(month2,n) * (target_date - middle(l))    &
                       +   field_in(month1,n) * (middle(l+1) - target_date)) &
                       / (middle(l+1) - middle(l))
!         if(field_out(n) .ne. 8.) write(0,201) n,field_in(month2,n),field_in(month2,n), &
!                                               field_out(n)
       enddo
       exit find_month
    endif
 enddo find_month
 
!201 format(i6,3(1x,e15.8))

 end subroutine monthly_interp_to_date

!=================================================================================================================
 subroutine monthly_min_max(npoints,field_in,field_min,field_max)
!=================================================================================================================

!input arguments:
 integer,intent(in):: npoints
 real(kind=RKIND),intent(in) ,dimension(12,npoints):: field_in

!output arguments:
 real(kind=RKIND),intent(out),dimension(npoints):: field_min,field_max

!local variables:
 integer:: n,nn
 real(kind=RKIND):: minner,maxxer

!-----------------------------------------------------------------------------------------------------------------
 
 do n = 1, npoints
    minner = field_in(1,n)
    maxxer = field_in(1,n)
    
    do nn = 2, 12
       if(field_in(nn,n) .lt. minner) then
          minner = field_in(nn,n)
       endif
       if(field_in(nn,n) .gt. maxxer) then
          maxxer = field_in(nn,n)
       endif
       field_min(n) = minner
       field_max(n) = maxxer
    enddo

 enddo

 end subroutine monthly_min_max

!=================================================================================================================
 end module mpas_atmphys_date_time
!=================================================================================================================
