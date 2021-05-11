module solver

   use prec

   implicit none

   ! maximum number of months in input dataset
   integer, parameter :: nmax=12*250

CONTAINS
   
   subroutine solvmid (lonindx, latindx, nmon, conv, dt, tmin, tmax, &
                       bbmin, maxiter, a, c, obsmean, ss, icnt, jcnt)

   integer, intent(in) :: lonindx
   integer, intent(in) :: latindx
   integer, intent(in) :: nmon
   integer, intent(in) :: maxiter        ! max number of iterations
   integer, intent(inout) :: icnt, jcnt  ! number of points smoothed
   real(r8), intent(in) :: conv, tmin, tmax, dt, bbmin
   real(r8), intent(inout) :: obsmean(nmon)
   real(r8), intent(in) :: a(nmon), c(nmon)
   real(r8), intent(inout) :: ss(nmon)

   integer i, n, imethod, n1, n2, nn, jj, i1, i2, i3, nnn, jend, &
           kk, j, k, kkk, nm, np
   integer jbeg(nmax)

   real(r8) :: relax, residmax, resid, dxm, dxp, s1, s2, addmax, addmin
   real(r8) :: r(nmax), avg(nmax), aa(nmax), bb(nmax), cc(nmax), add(nmax)

   double precision s(nmax), sum

   imethod = 1
! ???  check following value
   relax = 1.0

   if (nmon > nmax) then
      write(6,*) 'error-- nmax not declared large enough in '
      write(6,*) 'subroutine solvmid'
      stop 999
   end if

! Check for occurance where obs monthly means are consecutively
! at upper and lower limits. If so, smooth data, being careful
! to preserve annual mean.

   do n=1,nmon
      add(n) = 0.0
   end do

   n2 = nmon
   do n=1,nmon
      n1 = n2
      n2 = n
      if ((obsmean(n2)-obsmean(n1)) > (tmax-tmin-2.*dt)) then
          add(n1) = add(n1) + (0.5*(obsmean(n2)-obsmean(n1)-tmax+tmin) + dt)/c(n1)
          add(n2) = add(n2) - (0.5*(obsmean(n2)-obsmean(n1)-tmax+tmin) + dt)/a(n2)
       else if ((obsmean(n1)-obsmean(n2)) > (tmax-tmin-2.*dt)) then
          add(n1) = add(n1) - (0.5*(obsmean(n1)-obsmean(n2)-tmax+tmin) + dt)/c(n1)
          add(n2) = add(n2) + (0.5*(obsmean(n1)-obsmean(n2)-tmax+tmin) + dt)/a(n2)
       end if
    end do
  
    nn = 0
    addmax = 0.0
    addmin = 0.0
    do n=1,nmon
         
       if (add(n) .ne. 0.0) then
!          write(6,*) 'monthly means go from one limit to the other in 1 month'
!          write(6,*) 'alat = ', latindx, ' i = ', lonindx, ' n = ', n, ' add = ', add(n)
          addmax = max(addmax, add(n))
          addmin = min(addmin, add(n))
          obsmean(n) = obsmean(n) + add(n)
          nn = nn + 1
       end if
    end do

    icnt = icnt + nn
    if (nn > 0) then
       jcnt = jcnt + 1
!        if (nmon == 12) then
!          write(6,*) 'Climatology: '
!        endif
       if (jcnt <= 200) then
          write(6,*)  nn, ' monthly values smoothed at lat= ', latindx, ' lon= ', lonindx
          write(6,'(a,1p,e14.7,a,e14.7)') 'max added = ', addmax, '  max subtracted = ', addmin
        end if

        if (jcnt == 200) then
          write(6,*) ' '
          if (nmon == 12) then
            write(6,*) 'No more warnings will be printed concerning smoothing of climatological data'
          else
            write(6,*) 'No more warnings will be printed concerning smoothing of monthly data'
          end if
          write(6,*) ' '
          if (nmon == 12) then
            write(6,*) 'No more warnings will be printed concerning smoothing of climatological data'
          else
            write(6,*) 'No more warnings will be printed concerning smoothing of monthly data'
          endif
          write(6,*) ' '
        endif
      endif

!    check if all are le tmin or all are ge tmax

      if (obsmean(1) <= (tmin+0.01*dt)) then
         do i=2,nmon
            if (obsmean(i) > (tmin+0.01*dt)) go to 99
         end do
         do i=1,nmon
            ss(i) = tmin
         end do
!        if (nmon == 12) write(6,*) 'Climatology: '
!        write(6,*) 'all values were at minimum at this grid cell:'
!        write(6,*) 'latitude = ', latindx, ' longitude = ', lonindx
         return

      else if (obsmean(1) >= (tmax-0.01*dt)) then

         do i=2,nmon
            if (obsmean(i) < (tmax-0.01*dt)) go to 99
         end do
         do i=1,nmon
            ss(i) = tmax
         end do
!        if (nmon == 12) write(6,*) 'Climatology: '
!        write(6,*) 'all values were at maximum at this grid cell:'
!        write(6,*) 'latitude = ', latindx, ' longitude = ', lonindx
        return

      end if

   99 jj = 0
      do i=1,nmon
         i1 = i
         i2 = mod(i,nmon) + 1
         i3 = mod((i+1), nmon) + 1
     
         if ((obsmean(i1) <= tmin+0.01*dt .and. obsmean(i2) <= tmin+0.01*dt .and. &
              obsmean(i3) > tmin+0.01*dt) .or. &
             (obsmean(i1) >= tmax-0.01*dt .and. obsmean(i2) >= tmax-0.01*dt .and. &
              obsmean(i3) < tmax-0.01*dt)) then
            jj = jj + 1
            jbeg(jj) = i2
         end if
      end do

      if (jj == 0) then    ! simple cyclic treatment

! Latest approximation of means (given mid-month values)

         nnn = 0
  105    nnn = nnn + 1
         sum = 0.0
         residmax = 0.0

         do n = 1, nmon
            nm = mod((n+nmon-2), nmon) + 1
            np = mod(n, nmon) + 1
            bb(n) = 0.0
            avg(n) = 0.0
            if (nnn < imethod) then
               call approx (tmin, tmax, a(n), c(n), ss(nm), ss(n), &
                            ss(np), aa(n), bb(n), cc(n), avg(n))

            else
               call numer (conv, tmin, tmax, bbmin, a(n), c(n), ss(nm), &
                           ss(n), ss(np), aa(n), bb(n), cc(n), avg(n))
            end if

            r(n) = obsmean(n) - avg(n) 
            sum = sum + r(n)**2
            residmax = max(residmax, abs(r(n)))
         end do

         resid = dsqrt(sum)/nmon
         if (residmax > conv) then
            if (nnn > maxiter*0.5) then
               write(6,'(a,i8,a,1p,e14.7,a,e14.7)') 'iteration = ', nnn, ' residual = ', &
                        resid, ' maximum residual = ', residmax
            end if

            if (nnn > maxiter*0.9) then
               write(6,*) ' '
               write(6,*) 'latitude = ', latindx, ' longitude = ', lonindx
               do n=1,nmon
                  write(6,'(8(1pe10.2))') obsmean(n), avg(n), r(n), s(n), ss(n), aa(n), bb(n), cc(n)
               end do
            end if

            if (nnn > maxiter) then
               write(6,*) 'latitude = ', latindx, ' longitude = ', lonindx 
               write(6,*) 'does not converge'
!JR uncomment so it barfs if no convergence
               call exit(1)
            end if
              
! Solve for new estimate of mid-month values

            call cyclic (lonindx, latindx, aa, bb, cc, cc(nmon), aa(1), r, s, nmon)
              
            do n=1,nmon
               ss(n) = ss(n) + relax*s(n)
            end do
              
! If ss exceeds tmax or tmin, then it should exceed it no
! more than absolutely necessary:

            do n=1,nmon
               nm = mod((n+nmon-2), nmon) + 1
               np = mod(n, nmon) + 1
               
               if (ss(n) > tmax) then
                  if (ss(nm) <= tmax) then
                     dxm = (ss(n)-tmax)/((ss(n)-ss(nm))*a(n))
                  else 
                     dxm = 0.0
                  end if

                  if (ss(np) <= tmax) then 
                     dxp = (ss(n)-tmax)/((ss(n)-ss(np))*c(n))
                  else
                     dxp = 0.0
                  end if

                  if ((dxm > 0.5) .and. (dxp > 0.5)) then
                     s1 = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                     s2 = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                     ss(n) = min(s1, s2)
                  else if ((dxm == 0.0) .and. (dxp == 0.0)) then
                     ss(n) = tmax
                  else if ((dxp == 0.0) .and. (dxm > 0.5)) then
                     ss(n) = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                  else if ((dxm == 0.0) .and. (dxp > 0.5)) then
                     ss(n) = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                  end if
                   
               else if (ss(n) < tmin) then

                  if (ss(nm) >= tmin) then
                     dxm = (ss(n)-tmin)/((ss(n)-ss(nm))*a(n))
                  else 
                     dxm = 0.0
                  end if
                   
                  if (ss(np) >= tmin) then 
                     dxp = (ss(n)-tmin)/((ss(n)-ss(np))*c(n))
                  else
                     dxp = 0.0
                  endif

                  if ((dxm > 0.5) .and. (dxp > 0.5)) then
                     s1 = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                     s2 = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                     ss(n) = min(s1, s2)
                  else if ((dxm == 0.0) .and. (dxp == 0.0)) then
                     ss(n) = tmin
                  else if ((dxp == 0.0) .and. (dxm > 0.5)) then
                     ss(n) = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                  elseif ((dxm == 0.0) .and. (dxp > 0.5)) then
                     ss(n) = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                  end if
               end if
            end do
            go to 105
         end if
      else
       
! Treat independent segments

         do j=1,jj
            jend = jbeg(j)
150         jend = jend + 1
            i1 = mod((jend-1), nmon) + 1
            i2 = mod(jend, nmon) + 1
            if ((obsmean(i1) <= tmin+0.01*dt .and. obsmean(i2) <= tmin+0.01*dt) .or. &
                (obsmean(i1) >= tmax-0.01*dt .and. obsmean(i2) >= tmax-0.01*dt)) then

! Calculate values for interval jbeg(j) to jend
! latest approximation of means (given mid-month values)

               nnn      = 0
205            nnn      = nnn + 1
               kk       = jend - jbeg(j) + 1
               n        = jbeg(j)
               avg(1)   = obsmean(n)
               r(1)     = 0.0
               n        = mod((jend-1), nmon) + 1
               avg(kk)  = obsmean(n)
               r(kk)    = 0.0
               sum      = 0.0
               residmax = 0.0

               do k = 2, kk-1
                  nm = mod((k+jbeg(j)-3), nmon) + 1
                  n  = mod((k+jbeg(j)-2), nmon) + 1
                  np = mod((k+jbeg(j)-1), nmon) + 1
                  bb(k) = 0.0
                  avg(k) = 0.0

                  if (nnn < imethod) then
                     call approx (tmin, tmax, a(n), c(n), ss(nm), ss(n), &
                                  ss(np), aa(k), bb(k), cc(k), avg(k))
                  else
                     call numer (conv, tmin, tmax, bbmin, a(n), c(n), ss(nm), &
                                 ss(n), ss(np), aa(k), bb(k), cc(k), avg(k))
                  end if

                  r(k) = obsmean(n) - avg(k) 
                  sum = sum + r(k)**2
                  residmax = max(residmax, abs(r(k)))
               end do

               resid = dsqrt(sum)/(kk-2)

               if (residmax > conv) then
                  if (nnn > maxiter*0.5) then
                     write(6,'(a,i8,a,i8,a,1p,e14.7,a,e14.7)') 'iter = ', nnn, ' kk = ', kk, &
                             ' residual = ', resid, ' maximum residual = ', residmax
                  end if

                  if (nnn > maxiter*0.9) then
                     write(6,*) ' '
                     write(6,*) 'latitude = ', latindx, ' longitude = ', lonindx 
                     do k=1,kk
                        n  = mod((k+jbeg(j)-2), nmon) + 1
                        write(6,'(8(1pe10.2))') obsmean(n), avg(k), r(k), s(k), ss(n), &
                                                aa(k), bb(k), cc(k)
                     end do
                  end if

                  if (nnn > maxiter) then
                     write(6,*) 'latitude = ', latindx, ' longitude = ', lonindx 
                     write(6,*) 'does not converge'
!JR uncomment so it barfs if no convergence
                     call exit(1)
                  end if

! Solve for new estimate of mid-month values

                  kkk = kk - 2
                  call tridag(lonindx, latindx, aa(2), bb(2), cc(2), r(2), s(2), kkk)
                  
                  do k=2,kk-1
                     n  = mod((k+jbeg(j)-2), nmon) + 1
                     ss(n) = ss(n) + relax*s(k)
                  end do
                   
! If ss exceeds tmax or tmin, then it should exceed it no
! more than absolutely necessary:

                  n  = mod((jbeg(j)-1), nmon) + 1
                  np  = mod(jbeg(j), nmon) + 1
                   
                  if (obsmean(n) >= (tmax-0.01*dt)) then 
                     ss(n) = max(tmax, (tmax + (tmax-ss(np))*c(n)/(2.-c(n))))
                  else
                     ss(n) = min(tmin, (tmin + (tmin-ss(np))*c(n)/(2.-c(n))))
                  end if

                  nm  = mod((jend+nmon-2), nmon) + 1
                  n  = mod((jend-1), nmon) + 1

                  if (obsmean(n) >= (tmax-0.01*dt)) then 
                     ss(n) = max(tmax, (tmax + (tmax-ss(nm))*a(n)/(2.-a(n))))
                  else
                     ss(n) = min(tmin, (tmin + (tmin-ss(nm))*a(n)/(2.-a(n))))
                  end if

                  do k=2,kk-1
                     nm = mod((k+jbeg(j)+nmon-3), nmon) + 1
                     n  = mod((k+jbeg(j)-2), nmon) + 1
                     np = mod((k+jbeg(j)-1), nmon) + 1
                     
                     if (ss(n) > tmax) then
                        if (ss(nm) <= tmax) then
                           dxm = (ss(n)-tmax)/((ss(n)-ss(nm))*a(n))
                        else 
                           dxm = 0.0
                        end if

                        if (ss(np) <= tmax) then 
                           dxp = (ss(n)-tmax)/((ss(n)-ss(np))*c(n))
                        else
                           dxp = 0.0
                        end if

                        if ((dxm > 0.5) .and. (dxp > 0.5)) then
                           s1 = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                           s2 = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                           ss(n) = min(s1, s2)
                        else if ((dxm == 0.0) .and. (dxp == 0.0)) then
                           ss(n) = tmax
                        else if ((dxp == 0.0) .and. (dxm > 0.5)) then
                           ss(n) = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                        else if ((dxm == 0.0) .and. (dxp > 0.5)) then
                           ss(n) = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                        end if
                        
                     else if (ss(n) < tmin) then

                        if (ss(nm) >= tmin) then
                           dxm = (ss(n)-tmin)/((ss(n)-ss(nm))*a(n))
                        else 
                           dxm = 0.0
                        end if
                        
                        if (ss(np) >= tmin) then 
                           dxp = (ss(n)-tmin)/((ss(n)-ss(np))*c(n))
                        else
                           dxp = 0.0
                        end if
                        
                        if ((dxm > 0.5) .and. (dxp > 0.5)) then
                           s1 = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                           s2 = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                           ss(n) = max(s1, s2)
                        else if ((dxm == 0.0) .and. (dxp == 0.0)) then
                           ss(n) = tmin
                        else if ((dxp == 0.0) .and. (dxm > 0.5)) then
                           ss(n) = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                        else if ((dxm == 0.0) .and. (dxp > 0.5)) then
                           ss(n) = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                        end if
                     end if
                  end do
                  go to 205
               end if
               go to 300
            else
               go to 150
            end if
             
            if ((ss(kk) > 700.0) .or. (ss(kk) < -700.0)) then
               write(6,*) 'exceeds 700 at lat = ', latindx, ' lonindx = ', lonindx
               write(6,*) 'kk = ', kk, ' obs = ', obsmean(kk-1), obsmean(kk), obsmean(kk+1)
               write(6,*) 'kk = ', kk, '  ss = ', ss(kk-1), ss(kk), ss(kk+1)
            end if
300      continue
         end do

! Fill in values where consecutive means are outside limits

         do i=1,nmon
            i1 = mod((i-2+nmon), nmon) + 1
            i2 = mod((i-1), nmon) + 1
            i3 = mod(i, nmon) + 1
            if (obsmean(i1) <= tmin+0.01*dt .and. obsmean(i2) <= tmin+0.01*dt .and. &
                obsmean(i3) <= tmin+0.01*dt) then
               ss(i2) = tmin
            else if (obsmean(i1) >= tmax-0.01*dt .and. obsmean(i2) >= tmax-0.01*dt .and. &
                     obsmean(i3) >= tmax-0.01*dt) then
               ss(i2) = tmax
            end if
         end do
      end if

      return
   end subroutine solvmid

   subroutine numer (conv, tmin, tmax, bbmin, a, &
                     c, ssm, ss, ssp, aa, &
                     bb, cc, avg)
! *********************************************************************

      real(r8), intent(in) :: conv
      real(r8), intent(in) :: tmin
      real(r8), intent(in) :: tmax
      real(r8), intent(in) :: bbmin
      real(r8), intent(in) :: a
      real(r8), intent(in) :: c
      real(r8), intent(in) :: ssm
      real(r8), intent(in) :: ss
      real(r8), intent(in) :: ssp
      real(r8), intent(out) :: aa
      real(r8), intent(out) :: bb
      real(r8), intent(out) :: cc
      real(r8), intent(out) :: avg

      real(r8) :: ssmm
      real(r8) :: ssmp
      real(r8) :: sssm
      real(r8) :: sssp
      real(r8) :: sspm
      real(r8) :: sspp
      real(r8) :: r

      avg = amean (tmin,tmax,a,c,ssm,ss,ssp)

      ssmm = ssm - conv
      ssmp = ssm + conv
      sssm = ss  - conv
      sssp = ss  + conv
      sspm = ssp - conv
      sspp = ssp + conv

      aa = (amean(tmin,tmax,a,c,ssmp,ss,ssp) - &
            amean(tmin,tmax,a,c,ssmm,ss,ssp)) / (2.*conv)

      bb = (amean(tmin,tmax,a,c,ssm,sssp,ssp) - &
            amean(tmin,tmax,a,c,ssm,sssm,ssp)) / (2.*conv)

      cc = (amean(tmin,tmax,a,c,ssm,ss,sspp) - &
            amean(tmin,tmax,a,c,ssm,ss,sspm)) / (2.*conv)

      aa = min(aa, bb)
      cc = min(cc, bb)
      
      if (bb < bbmin) then
         bb = bbmin
         r = 0.2*bbmin
         aa = max(r, aa)
         cc = max(r, cc)
      endif
      
      return
   end subroutine numer

   real(r8) function amean (tmin, tmax, a, c, ssm, ss, ssp)

      real(r8), intent(in) :: tmin, tmax, a, c, ssm, ss, ssp

      real(r8) :: dx, dy, avg

      avg = 0.0

      if (ss <=  tmin) then
         if (ssm <= tmin) then
            avg = avg + tmin*0.5
         else if (ssm >= tmax) then
            dx = (ss-tmin)/((ss-ssm)*a)
            dy = (ss-tmax)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + tmin*0.5
            else if (dy <= 0.5) then
               avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            else
               avg = avg + tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
            end if
         else
            dx = (ss-tmin)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + tmin*0.5
            else
               avg = avg + tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
            end if
         end if

         if (ssp <= tmin) then
            avg = avg + tmin*0.5
         else if (ssp >= tmax) then
            dx = (ss-tmin)/((ss-ssp)*c)
            dy = (ss-tmax)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + tmin*0.5
            else if (dy <= 0.5) then
               avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            else
               avg = avg + tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
            end if
         else
            dx = (ss-tmin)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + tmin*0.5
            else
               avg = avg + tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
            end if
         end if
      else if (ss >= tmax) then
         if (ssm >= tmax) then
            avg = avg + tmax*0.5
         else if (ssm <= tmin) then
            dx = (ss-tmax)/((ss-ssm)*a)
            dy = (ss-tmin)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + tmax*0.5
            else if (dy <= 0.5) then
               avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            else
               avg = avg + tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
            end if
         else
            dx = (ss-tmax)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + tmax*0.5
            else
               avg = avg + tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
            end if
         end if
         
         if (ssp >= tmax) then
            avg = avg + tmax*0.5
         else if (ssp <= tmin) then
            dx = (ss-tmax)/((ss-ssp)*c)
            dy = (ss-tmin)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + tmax*0.5
            else if (dy <= 0.5) then
               avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            else
               avg = avg + tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
            end if
         else
            dx = (ss-tmax)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + tmax*0.5
            else
               avg = avg + tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
            end if
         end if
         
      else

         if (ssm <= tmin) then
            dx = (ss-tmin)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
            else
               avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
            endif
         else if (ssm >= tmax) then
            dx = (ss-tmax)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
            else
               avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
            end if
         else
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
         end if
        
         if (ssp <= tmin) then
            dx = (ss-tmin)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
            else
               avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
            end if
         else if (ssp >= tmax) then
            dx = (ss-tmax)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
            else
               avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
            endif
         else
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
         end if
      end if
      amean = avg

      return
   end function amean

   subroutine approx (tmin, tmax, a, c, ssm, ss, ssp, aa, bb, cc, avg)

      real(r8), intent(in) :: tmin, tmax, a, c, ssm, ss, ssp
      real(r8), intent(inout) :: avg
      real(r8), intent(out) :: aa, bb, cc

      real(r8) :: dx, dy

      if (ss <= tmin) then
         if (ssm <= tmin) then
            avg = avg + tmin*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
         else if (ssm >= tmax) then
            dx = (ss-tmin)/((ss-ssm)*a)
            dy = (ss-tmax)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + tmin*0.5
               aa = a/32.
               bb = bb + 0.125 - a/32.
            else if (dy <= 0.5) then
               avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
               aa = a/16.
               bb = bb + 0.25 - a/16. 
            else
               avg = avg + tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
               aa = a/16.
               bb = bb + 0.25 - a/16.
            end if
         else
            dx = (ss-tmin)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + tmin*0.5
               aa = a/32.
               bb = bb + 0.125 - a/32.
            else
               avg = avg + tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
               aa = a/16.
               bb = bb + 0.25 - a/16.
            end if
         end if

         if (ssp <= tmin) then
            avg = avg + tmin*0.5
            cc = c/32.
            bb = bb + 0.125 - c/32.
         else if (ssp >= tmax) then
            dx = (ss-tmin)/((ss-ssp)*c)
            dy = (ss-tmax)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + tmin*0.5
               cc = c/32.
               bb = bb + 0.125 - c/32.
            elseif (dy <= 0.5) then
               avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
               cc = c/16.
               bb = bb + 0.25 - c/16.
            else
               avg = avg + tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
               cc = c/16.
               bb = bb + 0.25 - c/16.
            end if
         else
            dx = (ss-tmin)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + tmin*0.5
               cc = c/32.
               bb = bb + 0.125 - c/32.
            else
               avg = avg + tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
               cc = c/16.
               bb = bb + 0.25 - c/16.
            end if
         end if
      else if (ss >= tmax) then
         if (ssm >= tmax) then
            avg = avg + tmax*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
         else if (ssm <= tmin) then
            dx = (ss-tmax)/((ss-ssm)*a)
            dy = (ss-tmin)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + tmax*0.5
               aa = a/32.
               bb = bb + 0.125 - a/32.
            else if (dy <= 0.5) then
               avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
               aa = a/16.
               bb = bb + 0.25 - a/16.
            else
               avg = avg + tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
               aa = a/16.
               bb = bb + 0.25 - a/16.
            end if
         else
            dx = (ss-tmax)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + tmax*0.5
               aa = a/32.
               bb = bb + 0.125 - a/32.
            else
               avg = avg + tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
               aa = a/16.
               bb = bb + 0.25 - a/16.
            end if
         end if

         if (ssp >= tmax) then
            avg = avg + tmax*0.5
            cc = c/32.
            bb = bb + 0.125 - c/32.
         else if (ssp <= tmin) then
            dx = (ss-tmax)/((ss-ssp)*c)
            dy = (ss-tmin)/((ss-ssp)*c)
            if (dx >= 0.5) then
               cc = c/32.
               bb = bb + 0.125 - c/32.
               avg = avg + tmax*0.5
            else if (dy <= 0.5) then
               avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
               cc = c/16.
               bb = bb + 0.25 - c/16.
            else
               avg = avg + tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
               cc = c/16.
               bb = bb + 0.25 - c/16.
            end if
         else
            dx = (ss-tmax)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + tmax*0.5
               cc = c/32.
               bb = bb + 0.125 - c/32.
            else
               avg = avg + tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
               cc = c/16.
               bb = bb + 0.25 - c/16.
            end if
         end if
      else
         if (ssm <= tmin) then
            dx = (ss-tmin)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
               aa = a/16.
               bb = bb + 0.25 - a/16.
            else
               avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
               aa = a/16.
               bb = bb + 0.25 - a/16.
            endif
         else if (ssm >= tmax) then
            dx = (ss-tmax)/((ss-ssm)*a)
            if (dx >= 0.5) then
               avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
               aa = a/16.
               bb = bb + 0.25 - a/16.
            else
               avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
               aa = a/16.
               bb = bb + 0.25 - a/16.
            end if
         else
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
            aa = a/8.
            bb = bb + 0.5 - a/8.
         end if

         if (ssp <= tmin) then
            dx = (ss-tmin)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
               cc = c/16.
               bb = bb + 0.25 - c/16.
            else
               avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
               cc = c/16.
               bb = bb + 0.25 - c/16.
            end if
         else if (ssp >= tmax) then
            dx = (ss-tmax)/((ss-ssp)*c)
            if (dx >= 0.5) then
               avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
               cc = c/16.
               bb = bb + 0.25 - c/16.
            else
               avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
               cc = c/16.
               bb = bb + 0.25 - c/16.
            end if
         else
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
            cc = c/8.
            bb = bb + 0.5 - c/8.
         end if
      end if

      return
   end subroutine approx

   subroutine tridag(lonindx,latindx,a,b,c,r,u,n)
      integer, intent(in) :: n
      integer, intent(in) :: lonindx, latindx
      real(r8), intent(in) :: a(n),b(n),c(n),r(n)
      double precision, intent(out) :: u(n)

      integer :: j
      real(r8) :: bet, gam(nmax)

      if (nmax < n) then
         write(6,*) 'Error nmax not declared large enough'
         write(6,*) 'in tridag'
         stop 999
      endif

      if (b(1) == 0.) then
         write(6,*) 'longitude = ', lonindx, '  latitude = ', latindx
!JR make stop if something bad happened
!JR          pause 'tridag: rewrite equations'
         write (6,*) 'tridag: rewrite equations'
         call exit(1)
      end if

      bet = b(1)
      u(1) = r(1)/bet

      if (n > 1) then
         do j=2,n
            gam(j) = c(j-1)/bet
            bet    = b(j)-a(j)*gam(j)
            if (bet == 0.) then
               write(6,*) 'longitude = ', lonindx, '  latitude = ', latindx
!JR make stop if something bad happened
!JR            pause 'tridag failed'
               write (6,*)'tridag failed'
               call exit(1)
            end if
            u(j)=(r(j)-a(j)*u(j-1))/bet
         end do
         do j=n-1,1,-1
            u(j)=u(j)-gam(j+1)*u(j+1)
         end do
      endif
      return
   end subroutine tridag

   subroutine cyclic (lonindx, latindx, a, b, c, alpha, beta, r, x, n)
      integer, intent(in) :: lonindx, latindx
      integer, intent(in) :: n
      real(r8), intent(in) :: alpha,beta,a(n),b(n),c(n),r(n)
      double precision, intent(out) :: x(n)
!u    uses tridag

      integer :: i
      real(r8) :: fact,gamma,bb(nmax),u(nmax)
      double precision z(nmax)

      if (n <= 2) then
         write(6,*) 'n too small in cyclic'
         stop 999
      end if
      if (n > nmax) then
         write(6,*) 'nmax too small in cyclic'
         stop 999
      end if

      gamma = -b(1)
      bb(1) = b(1) - gamma
      bb(n) = b(n) - alpha*beta/gamma
      do i=2,n-1
         bb(i)=b(i)
      end do

      call tridag (lonindx, latindx, a, bb, c, r, x, n)

      u(1) = gamma
      u(n) = alpha
      do i=2,n-1
         u(i) = 0.
      end do

      call tridag (lonindx, latindx, a, bb, c, u, z, n)

      fact = (x(1) + beta*x(n)/gamma)/(1. + z(1) + beta*z(n)/gamma)
      do i=1,n
         x(i)=x(i)-fact*z(i)
      end do

      return
   end subroutine cyclic
end module solver
