CVD$G NOVECTOR
CVD$G NOCONCUR
*deck xerabt
      subroutine xerabt(messg,nmessg)
c***begin prologue  xerabt
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  abort program execution and print error message.
c***description
c
c     abstract
c        ***note*** machine dependent routine
c        xerabt aborts the execution of the program.
c        the error message causing the abort is given in the calling
c        sequence, in case one needs it for printing on a dayfile,
c        for example.
c
c     description of parameters
c        messg and nmessg are as in xerror, except that nmessg may
c        be zero, in which case no message is being supplied.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  1 august 1982
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  (none)
c***end prologue  xerabt
      dimension messg(nmessg)
c***first executable statement  xerabt
      call exit(1)
      end
*deck xerctl
      subroutine xerctl(messg1,nmessg,nerr,level,kontrl)
c***begin prologue  xerctl
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  allow user control over handling of errors.
c***description
c
c     abstract
c        allows user control over handling of individual errors.
c        just after each message is recorded, but before it is
c        processed any further (i.e., before it is printed or
c        a decision to abort is made), a call is made to xerctl.
c        if the user has provided his own version of xerctl, he
c        can then override the value of kontrol used in processing
c        this message by redefining its value.
c        kontrl may be set to any value from -2 to 2.
c        the meanings for kontrl are the same as in xsetf, except
c        that the value of kontrl changes only for this message.
c        if kontrl is set to a value outside the range from -2 to 2,
c        it will be moved back into that range.
c
c     description of parameters
c
c      --input--
c        messg1 - the first word (only) of the error message.
c        nmessg - same as in the call to xerror or xerrwv.
c        nerr   - same as in the call to xerror or xerrwv.
c        level  - same as in the call to xerror or xerrwv.
c        kontrl - the current value of the control flag as set
c                 by a call to xsetf.
c
c      --output--
c        kontrl - the new value of kontrl.  if kontrl is not
c                 defined, it will remain at its original value.
c                 this changed value of control affects only
c                 the current occurrence of the current message.
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  (none)
c***end prologue  xerctl
      character*20 messg1
c***first executable statement  xerctl
      return
      end
*deck xerprt
      subroutine xerprt(messg,nmessg)
c***begin prologue  xerprt
c***date written   790801   (yymmdd)
c***revision date  851213   (yymmdd)
c***category no.  r3
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  print error messages.
c***description
c
c     abstract
c        print the hollerith message in messg, of length nmessg,
c        on each file indicated by xgetua.
c     latest revision ---  1 august 1985
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  i1mach,xgetua
c***end prologue  xerprt
      integer lun(5)
      character*(*) messg
c     obtain unit numbers and write line to each unit
c***first executable statement  xerprt
      call xgetua(lun,nunit)
      lenmes = len(messg)
      do 20 kunit=1,nunit
         iunit = lun(kunit)
         if (iunit.eq.0) iunit = i1mach(4)
         do 10 ichar=1,lenmes,72
            last = min0(ichar+71 , lenmes)
            write (iunit,'(1x,a)') messg(ichar:last)
   10    continue
   20 continue
      return
      end
*deck xerror
      subroutine xerror(messg,nmessg,nerr,level)
c***begin prologue  xerror
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  process an error (diagnostic) message.
c***description
c
c     abstract
c        xerror processes a diagnostic message, in a manner
c        determined by the value of level and the current value
c        of the library error control flag, kontrl.
c        (see subroutine xsetf for details.)
c
c     description of parameters
c      --input--
c        messg - the hollerith message to be processed, containing
c                no more than 72 characters.
c        nmessg- the actual number of characters in messg.
c        nerr  - the error number associated with this message.
c                nerr must not be zero.
c        level - error category.
c                =2 means this is an unconditionally fatal error.
c                =1 means this is a recoverable error.  (i.e., it is
c                   non-fatal if xsetf has been appropriately called.)
c                =0 means this is a warning message only.
c                =-1 means this is a warning message which is to be
c                   printed at most once, regardless of how many
c                   times this call is executed.
c
c     examples
c        call xerror('smooth -- num was zero.',23,1,2)
c        call xerror('integ  -- less than full accuracy achieved.',
c    1                43,2,1)
c        call xerror('rooter -- actual zero of f found before interval f
c    1ully collapsed.',65,3,0)
c        call xerror('exp    -- underflows being set to zero.',39,1,-1)
c
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  xerrwv
c***end prologue  xerror
      character*(*) messg
c***first executable statement  xerror
      call xerrwv(messg,nmessg,nerr,level,0,0,0,0,0.,0.)
      return
      end
*deck xerrwv
      subroutine xerrwv(messg,nmessg,nerr,level,ni,i1,i2,nr,r1,r2)
c***begin prologue  xerrwv
c***date written   800319   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  process an error message allowing 2 integer and 2 real
c            values to be included in the message.
c***description
c
c     abstract
c        xerrwv processes a diagnostic message, in a manner
c        determined by the value of level and the current value
c        of the library error control flag, kontrl.
c        (see subroutine xsetf for details.)
c        in addition, up to two integer values and two real
c        values may be printed along with the message.
c
c     description of parameters
c      --input--
c        messg - the hollerith message to be processed.
c        nmessg- the actual number of characters in messg.
c        nerr  - the error number associated with this message.
c                nerr must not be zero.
c        level - error category.
c                =2 means this is an unconditionally fatal error.
c                =1 means this is a recoverable error.  (i.e., it is
c                   non-fatal if xsetf has been appropriately called.)
c                =0 means this is a warning message only.
c                =-1 means this is a warning message which is to be
c                   printed at most once, regardless of how many
c                   times this call is executed.
c        ni    - number of integer values to be printed. (0 to 2)
c        i1    - first integer value.
c        i2    - second integer value.
c        nr    - number of real values to be printed. (0 to 2)
c        r1    - first real value.
c        r2    - second real value.
c
c     examples
c        call xerrwv('smooth -- num (=i1) was zero.',29,1,2,
c    1   1,num,0,0,0.,0.)
c        call xerrwv('quadxy -- requested error (r1) less than minimum (
c    1r2).,54,77,1,0,0,0,2,errreq,errmin)
c
c     latest revision ---  1 august 1985
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  fdump,i1mach,j4save,xerabt,xerctl,xerprt,xersav,
c                    xgetua
c***end prologue  xerrwv
      character*(*) messg
      character*20 lfirst
      character*37 form
      dimension lun(5)
c     get flags
c***first executable statement  xerrwv
      lkntrl = j4save(2,0,.false.)
      maxmes = j4save(4,0,.false.)
c     check for valid input
      if ((nmessg.gt.0).and.(nerr.ne.0).and.
     1    (level.ge.(-1)).and.(level.le.2)) go to 10
         if (lkntrl.gt.0) call xerprt('fatal error in...',17)
         call xerprt('xerror -- invalid input',23)
c        if (lkntrl.gt.0) call fdump
         if (lkntrl.gt.0) call xerprt('job abort due to fatal error.',
     1  29)
         if (lkntrl.gt.0) call xersav(' ',0,0,0,kdummy)
         call xerabt('xerror -- invalid input',23)
         return
   10 continue
c     record message
      junk = j4save(1,nerr,.true.)
      call xersav(messg,nmessg,nerr,level,kount)
c     let user override
      lfirst = messg
      lmessg = nmessg
      lerr = nerr
      llevel = level
      call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
c     reset to original values
      lmessg = nmessg
      lerr = nerr
      llevel = level
      lkntrl = max0(-2,min0(2,lkntrl))
      mkntrl = iabs(lkntrl)
c     decide whether to print message
      if ((llevel.lt.2).and.(lkntrl.eq.0)) go to 100
      if (((llevel.eq.(-1)).and.(kount.gt.min0(1,maxmes)))
     1.or.((llevel.eq.0)   .and.(kount.gt.maxmes))
     2.or.((llevel.eq.1)   .and.(kount.gt.maxmes).and.(mkntrl.eq.1))
     3.or.((llevel.eq.2)   .and.(kount.gt.max0(1,maxmes)))) go to 100
         if (lkntrl.le.0) go to 20
            call xerprt(' ',1)
c           introduction
            if (llevel.eq.(-1)) call xerprt
     1('warning message...this message will only be printed once.',57)
            if (llevel.eq.0) call xerprt('warning in...',13)
            if (llevel.eq.1) call xerprt
     1      ('recoverable error in...',23)
            if (llevel.eq.2) call xerprt('fatal error in...',17)
   20    continue
c        message
         call xerprt(messg,lmessg)
         call xgetua(lun,nunit)
         isizei = log10(float(i1mach(9))) + 1.0
         isizef = log10(float(i1mach(10))**i1mach(11)) + 1.0
         do 50 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
            do 22 i=1,min(ni,2)
               write (form,21) i,isizei
   21          format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
               if (i.eq.1) write (iunit,form) i1
               if (i.eq.2) write (iunit,form) i2
   22       continue
            do 24 i=1,min(nr,2)
               write (form,23) i,isizef+10,isizef
   23          format ('(11x,21hin above message, r',i1,'=,e',
     1         i2,'.',i2,')')
               if (i.eq.1) write (iunit,form) r1
               if (i.eq.2) write (iunit,form) r2
   24       continue
            if (lkntrl.le.0) go to 40
c              error number
               write (iunit,30) lerr
   30          format (15h error number =,i10)
   40       continue
   50    continue
c        trace-back
c        if (lkntrl.gt.0) call fdump
  100 continue
      ifatal = 0
      if ((llevel.eq.2).or.((llevel.eq.1).and.(mkntrl.eq.2)))
     1ifatal = 1
c     quit here if message is not fatal
      if (ifatal.le.0) return
      if ((lkntrl.le.0).or.(kount.gt.max0(1,maxmes))) go to 120
c        print reason for abort
         if (llevel.eq.1) call xerprt
     1   ('job abort due to unrecovered error.',35)
         if (llevel.eq.2) call xerprt
     1   ('job abort due to fatal error.',29)
c        print error summary
         call xersav(' ',-1,0,0,kdummy)
  120 continue
c     abort
      if ((llevel.eq.2).and.(kount.gt.max0(1,maxmes))) lmessg = 0
      call xerabt(messg,lmessg)
      return
      end
*deck xersav
      subroutine xersav(messg,nmessg,nerr,level,icount)
c***begin prologue  xersav
c***date written   800319   (yymmdd)
c***revision date  851213   (yymmdd)
c***category no.  r3
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  record that an error has occurred.
c***description
c
c     abstract
c        record that this error occurred.
c
c     description of parameters
c     --input--
c       messg, nmessg, nerr, level are as in xerror,
c       except that when nmessg=0 the tables will be
c       dumped and cleared, and when nmessg is less than zero the
c       tables will be dumped and not cleared.
c     --output--
c       icount will be the number of times this message has
c       been seen, or zero if the table has overflowed and
c       does not contain this message specifically.
c       when nmessg=0, icount will not be altered.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  1 august 1985
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  i1mach,xgetua
c***end prologue  xersav
      integer lun(5)
      character*(*) messg
      character*20 mestab(10),mes
      dimension nertab(10),levtab(10),kount(10)
      save mestab,nertab,levtab,kount,kountx
c     next two data statements are necessary to provide a blank
c     error table initially
      data kount(1),kount(2),kount(3),kount(4),kount(5),
     1     kount(6),kount(7),kount(8),kount(9),kount(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      data kountx/0/
c***first executable statement  xersav
      if (nmessg.gt.0) go to 80
c     dump the table
         if (kount(1).eq.0) return
c        print to each unit
         call xgetua(lun,nunit)
         do 60 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
c           print table header
            write (iunit,10)
   10       format (32h0          error message summary/
     1      51h message start             nerr     level     count)
c           print body of table
            do 20 i=1,10
               if (kount(i).eq.0) go to 30
               write (iunit,15) mestab(i),nertab(i),levtab(i),kount(i)
   15          format (1x,a20,3i10)
   20       continue
   30       continue
c           print number of other errors
            if (kountx.ne.0) write (iunit,40) kountx
   40       format (41h0other errors not individually tabulated=,i10)
            write (iunit,50)
   50       format (1x)
   60    continue
         if (nmessg.lt.0) return
c        clear the error tables
         do 70 i=1,10
   70       kount(i) = 0
         kountx = 0
         return
   80 continue
c     process a message...
c     search for this messg, or else an empty slot for this messg,
c     or else determine that the error table is full.
      mes = messg
      do 90 i=1,10
         ii = i
         if (kount(i).eq.0) go to 110
         if (mes.ne.mestab(i)) go to 90
         if (nerr.ne.nertab(i)) go to 90
         if (level.ne.levtab(i)) go to 90
         go to 100
   90 continue
c     three possible cases...
c     table is full
         kountx = kountx+1
         icount = 1
         return
c     message found in table
  100    kount(ii) = kount(ii) + 1
         icount = kount(ii)
         return
c     empty slot found for new message
  110    mestab(ii) = mes
         nertab(ii) = nerr
         levtab(ii) = level
         kount(ii)  = 1
         icount = 1
         return
      end
      subroutine xgetf(kontrl)
c***begin prologue  xgetf
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  return the current value of the error control flag.
c***description
c
c   abstract
c        xgetf returns the current value of the error control flag
c        in kontrl.  see subroutine xsetf for flag value meanings.
c        (kontrl is an output parameter only.)
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  7 june 1978
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xgetf
c***first executable statement  xgetf
      kontrl = j4save(2,0,.false.)
      return
      end
*deck xgetua
      subroutine xgetua(iunita,n)
c***begin prologue  xgetua
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  return unit number(s) to which error messages are being
c            sent.
c***description
c
c     abstract
c        xgetua may be called to determine the unit number or numbers
c        to which error messages are being sent.
c        these unit numbers may have been set by a call to xsetun,
c        or a call to xsetua, or may be a default value.
c
c     description of parameters
c      --output--
c        iunit - an array of one to five unit numbers, depending
c                on the value of n.  a value of zero refers to the
c                default unit, as defined by the i1mach machine
c                constant routine.  only iunit(1),...,iunit(n) are
c                defined by xgetua.  the values of iunit(n+1),...,
c                iunit(5) are not defined (for n .lt. 5) or altered
c                in any way by xgetua.
c        n     - the number of units to which copies of the
c                error messages are being sent.  n will be in the
c                range from 1 to 5.
c
c     latest revision ---  19 mar 1980
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xgetua
      dimension iunita(5)
c***first executable statement  xgetua
      n = j4save(5,0,.false.)
      do 30 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         iunita(i) = j4save(index,0,.false.)
   30 continue
      return
      end
*deck j4save
      function j4save(iwhich,ivalue,iset)
c***begin prologue  j4save
c***refer to  xerror
c***routines called  (none)
c***description
c
c     abstract
c        j4save saves and recalls several global variables needed
c        by the library error handling routines.
c
c     description of parameters
c      --input--
c        iwhich - index of item desired.
c                = 1 refers to current error number.
c                = 2 refers to current error control flag.
c                 = 3 refers to current unit number to which error
c                    messages are to be sent.  (0 means use standard.)
c                 = 4 refers to the maximum number of times any
c                     message is to be printed (as set by xermax).
c                 = 5 refers to the total number of units to which
c                     each error message is to be written.
c                 = 6 refers to the 2nd unit for error messages
c                 = 7 refers to the 3rd unit for error messages
c                 = 8 refers to the 4th unit for error messages
c                 = 9 refers to the 5th unit for error messages
c        ivalue - the value to be set for the iwhich-th parameter,
c                 if iset is .true. .
c        iset   - if iset=.true., the iwhich-th parameter will be
c                 given the value, ivalue.  if iset=.false., the
c                 iwhich-th parameter will be unchanged, and ivalue
c                 is a dummy parameter.
c      --output--
c        the (old) value of the iwhich-th parameter will be returned
c        in the function value, j4save.
c
c     written by ron jones, with slatec common math library subcommittee
c    adapted from bell laboratories port library error handler
c     latest revision ---  1 august 1985
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***end prologue  j4save
      logical iset
      integer iparam(9)
      save iparam
      data iparam(1),iparam(2),iparam(3),iparam(4)/0,2,0,10/
      data iparam(5)/1/
      data iparam(6),iparam(7),iparam(8),iparam(9)/0,0,0,0/
c***first executable statement  j4save
      j4save = iparam(iwhich)
      if (iset) iparam(iwhich) = ivalue
      return
      end
*deck xerclr
      subroutine xerclr
c***begin prologue  xerclr
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  reset current error number to zero.
c***description
c
c     abstract
c        this routine simply resets the current error number to zero.
c        this may be necessary to do in order to determine that
c        a certain error has occurred again since the last time
c        numxer was referenced.
c
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xerclr
c***first executable statement  xerclr
      junk = j4save(1,0,.true.)
      return
      end
      subroutine xerdmp
c***begin prologue  xerdmp
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  print the error tables and then clear them.
c***description
c
c     abstract
c        xerdmp prints the error tables, then clears them.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  7 june 1978
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  xersav
c***end prologue  xerdmp
c***first executable statement  xerdmp
      call xersav(' ',0,0,0,kount)
      return
      end
      subroutine xermax(max)
c***begin prologue  xermax
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  set maximum number of times any error message is to be
c            printed.
c***description
c
c     abstract
c        xermax sets the maximum number of times any message
c        is to be printed.  that is, non-fatal messages are
c        not to be printed after they have occured max times.
c        such non-fatal messages may be printed less than
c        max times even if they occur max times, if error
c        suppression mode (kontrl=0) is ever in effect.
c
c     description of parameter
c      --input--
c        max - the maximum number of times any one message
c              is to be printed.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  7 june 1978
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xermax
c***first executable statement  xermax
      junk = j4save(4,max,.true.)
      return
      end
      subroutine xgetun(iunit)
c***begin prologue  xgetun
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  return the (first) output file to which error messages
c            are being sent.
c***description
c
c     abstract
c        xgetun gets the (first) output file to which error messages
c        are being sent.  to find out if more than one file is being
c        used, one must use the xgetua routine.
c
c     description of parameter
c      --output--
c        iunit - the logical unit number of the  (first) unit to
c                which error messages are being sent.
c                a value of zero means that the default file, as
c                defined by the i1mach routine, is being used.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision --- 23 may 1979
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xgetun
c***first executable statement  xgetun
      iunit = j4save(3,0,.false.)
      return
      end
      subroutine xsetf(kontrl)
c***begin prologue  xsetf
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3a
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  set the error control flag.
c***description
c
c     abstract
c        xsetf sets the error control flag value to kontrl.
c        (kontrl is an input parameter only.)
c        the following table shows how each message is treated,
c        depending on the values of kontrl and level.  (see xerror
c        for description of level.)
c
c        if kontrl is zero or negative, no information other than the
c        message itself (including numeric values, if any) will be
c        printed.  if kontrl is positive, introductory messages,
c        trace-backs, etc., will be printed in addition to the message.
c
c              iabs(kontrl)
c        level        0              1              2
c        value
c          2        fatal          fatal          fatal
c
c          1     not printed      printed         fatal
c
c          0     not printed      printed        printed
c
c         -1     not printed      printed        printed
c                                  only           only
c                                  once           once
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  19 mar 1980
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save,xerrwv
c***end prologue  xsetf
c***first executable statement  xsetf
      if ((kontrl.ge.(-2)).and.(kontrl.le.2)) go to 10
         call xerrwv('xsetf  -- invalid value of kontrl (i1).',33,1,2,
     1  1,kontrl,0,0,0.,0.)
         return
   10 junk = j4save(2,kontrl,.true.)
      return
      end
      subroutine xsetua(iunita,n)
c***begin prologue  xsetua
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3b
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  set logical unit numbers (up to 5) to which error
c            messages are to be sent.
c***description
c
c     abstract
c        xsetua may be called to declare a list of up to five
c        logical units, each of which is to receive a copy of
c        each error message processed by this package.
c        the purpose of xsetua is to allow simultaneous printing
c        of each error message on, say, a main output file,
c        an interactive terminal, and other files such as graphics
c        communication files.
c
c     description of parameters
c      --input--
c        iunit - an array of up to five unit numbers.
c                normally these numbers should all be different
c                (but duplicates are not prohibited.)
c        n     - the number of unit numbers provided in iunit
c                must have 1 .le. n .le. 5.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  19 mar 1980
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save,xerrwv
c***end prologue  xsetua
      dimension iunita(5)
c***first executable statement  xsetua
      if ((n.ge.1).and.(n.le.5)) go to 10
         call xerrwv('xsetua -- invalid value of n (i1).',34,1,2,
     1  1,n,0,0,0.,0.)
         return
   10 continue
      do 20 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         junk = j4save(index,iunita(i),.true.)
   20 continue
      junk = j4save(5,n,.true.)
      return
      end
      subroutine xsetun(iunit)
c***begin prologue  xsetun
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3b
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  set output file to which error messages are to be sent.
c***description
c
c     abstract
c        xsetun sets the output file to which error messages are to
c        be sent.  only one file will be used.  see xsetua for
c        how to declare more than one file.
c
c     description of parameter
c      --input--
c        iunit - an input parameter giving the logical unit number
c                to which error messages are to be sent.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  7 june 1978
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xsetun
c***first executable statement  xsetun
      junk = j4save(3,iunit,.true.)
      junk = j4save(5,1,.true.)
      return
      end
      FUNCTION RAND(R)
C***BEGIN PROLOGUE  RAND
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  L6A21
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),TYPE=SINGLE PRECISION(RAND-S),
C             RANDOM NUMBER,SPECIAL FUNCTIONS,UNIFORM
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Generates a uniformly distributed random number.
C***DESCRIPTION
C
C      This pseudo-random number generator is portable among a wide
C variety of computers.  RAND(R) undoubtedly is not as good as many
C readily available installation dependent versions, and so this
C routine is not recommended for widespread usage.  Its redeeming
C feature is that the exact same random numbers (to within final round-
C off error) can be generated from machine to machine.  Thus, programs
C that make use of random numbers can be easily transported to and
C checked in a new environment.
C      The random numbers are generated by the linear congruential
C method described, e.g., by Knuth in Seminumerical Methods (p.9),
C Addison-Wesley, 1969.  Given the I-th number of a pseudo-random
C sequence, the I+1 -st number is generated from
C             X(I+1) = (A*X(I) + C) MOD M,
C where here M = 2**22 = 4194304, C = 1731 and several suitable values
C of the multiplier A are discussed below.  Both the multiplier A and
C random number X are represented in double precision as two 11-bit
C words.  The constants are chosen so that the period is the maximum
C possible, 4194304.
C      In order that the same numbers be generated from machine to
C machine, it is necessary that 23-bit integers be reducible modulo
C 2**11 exactly, that 23-bit integers be added exactly, and that 11-bit
C integers be multiplied exactly.  Furthermore, if the restart option
C is used (where R is between 0 and 1), then the product R*2**22 =
C R*4194304 must be correct to the nearest integer.
C      The first four random numbers should be .0004127026,
C .6750836372, .1614754200, and .9086198807.  The tenth random number
C is .5527787209, and the hundredth is .3600893021 .  The thousandth
C number should be .2176990509 .
C      In order to generate several effectively independent sequences
C with the same generator, it is necessary to know the random number
C for several widely spaced calls.  The I-th random number times 2**22,
C where I=K*P/8 and P is the period of the sequence (P = 2**22), is
C still of the form L*P/8.  In particular we find the I-th random
C number multiplied by 2**22 is given by
C I   =  0  1*P/8  2*P/8  3*P/8  4*P/8  5*P/8  6*P/8  7*P/8  8*P/8
C RAND=  0  5*P/8  2*P/8  7*P/8  4*P/8  1*P/8  6*P/8  3*P/8  0
C Thus the 4*P/8 = 2097152 random number is 2097152/2**22.
C      Several multipliers have been subjected to the spectral test
C (see Knuth, p. 82).  Four suitable multipliers roughly in order of
C goodness according to the spectral test are
C    3146757 = 1536*2048 + 1029 = 2**21 + 2**20 + 2**10 + 5
C    2098181 = 1024*2048 + 1029 = 2**21 + 2**10 + 5
C    3146245 = 1536*2048 +  517 = 2**21 + 2**20 + 2**9 + 5
C    2776669 = 1355*2048 + 1629 = 5**9 + 7**7 + 1
C
C      In the table below LOG10(NU(I)) gives roughly the number of
C random decimal digits in the random numbers considered I at a time.
C C is the primary measure of goodness.  In both cases bigger is better.
C
C                   LOG10 NU(I)              C(I)
C       A       I=2  I=3  I=4  I=5    I=2  I=3  I=4  I=5
C
C    3146757    3.3  2.0  1.6  1.3    3.1  1.3  4.6  2.6
C    2098181    3.3  2.0  1.6  1.2    3.2  1.3  4.6  1.7
C    3146245    3.3  2.2  1.5  1.1    3.2  4.2  1.1  0.4
C    2776669    3.3  2.1  1.6  1.3    2.5  2.0  1.9  2.6
C   Best
C    Possible   3.3  2.3  1.7  1.4    3.6  5.9  9.7  14.9
C
C             Input Argument --
C R      If R=0., the next random number of the sequence is generated.
C        If R .LT. 0., the last generated number will be returned for
C          possible use in a restart procedure.
C        If R .GT. 0., the sequence of random numbers will start with
C          the seed R mod 1.  This seed is also returned as the value of
C          RAND provided the arithmetic is done exactly.
C
C             Output Value --
C RAND   a pseudo-random number between 0. and 1.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  RAND
      SAVE IA1, IA0, IA1MA0, IC, IX1, IX0
      DATA IA1, IA0, IA1MA0 /1536, 1029, 507/
      DATA IC /1731/
      DATA IX1, IX0 /0, 0/
C***FIRST EXECUTABLE STATEMENT  RAND
      IF (R.LT.0.) GO TO 10
      IF (R.GT.0.) GO TO 20
C
C           A*X = 2**22*IA1*IX1 + 2**11*(IA1*IX1 + (IA1-IA0)*(IX0-IX1)
C                   + IA0*IX0) + IA0*IX0
C
      IY0 = IA0*IX0
      IY1 = IA1*IX1 + IA1MA0*(IX0-IX1) + IY0
      IY0 = IY0 + IC
      IX0 = MOD (IY0, 2048)
      IY1 = IY1 + (IY0-IX0)/2048
      IX1 = MOD (IY1, 2048)
C
 10   RAND = IX1*2048 + IX0
      RAND = RAND / 4194304.
      RETURN
C
 20   IX1 = AMOD(R,1.)*4194304. + 0.5
      IX0 = MOD (IX1, 2048)
      IX1 = (IX1-IX0)/2048
      GO TO 10
C
      END
