!! ******************************************************************
!! The routines listed in this file "adgaquad_mod.F90" are performing
!! Numerical Integrations using some kind of 
!! adaptive Gauss quadrature.
!! They are taken from the Internet (http://www.netlib.org)
!! and parts of different software packages / libraries.
!! ******************************************************************
!! For any restrictions on the use of the routines, please see
!! the original web site. 
!! ******************************************************************
!! Changes: calls to error handler 'xerror()' replaced by
!!          WRITE(7,*) - statements.
!! ******************************************************************
!! list of routines and the libraries they are taken from:
!!  dqag		calling routine, bounded integration interval
!!			QUADPACK; calls: dqage
!!  dqage		the integration routine, bounded interval
!!			QUADPACK; calls: sd1mach,dqk15,dqk21,dqk31,
!!					dqk41,dqk51,dqk61,dqpsrt
!!  dqagi		calling routine, unbounded (semi-infinite or 
!!			infinite) integration interval
!!			QUADPACK; calls: dqagie
!!  dqagie		the integration routine, unbounded interval
!!			QUADPACK; calls: sd1mach,dqelg,dqk15i,dqpsrt
!! ------------------------------------------------------------------
!!  dqk15		QUADPACK; calls: sd1mach
!!  dqk21		QUADPACK; calls: sd1mach
!!  dqk31		QUADPACK; calls: sd1mach
!!  dqk41		QUADPACK; calls: sd1mach
!!  dqk51		QUADPACK; calls: sd1mach
!!  dqk61		QUADPACK; calls: sd1mach
!!  dqpsrt		QUADPACK; calls: none
!!  dqk15i		QUADPACK; calls: sd1mach
!!  dqelg		QUADPACK; calls: sd1mach
!! ------------------------------------------------------------------
!!  xerror		Error handling routine
!!			ALLIANT (/quad); calls: xerrwv
!!  xerrwv		Error handling routine
!!			SODEPACK; calls: none
!!  d1mach		determine machine parameters (accuracies)
!!			BLAS; calls: none
!! ------------------------------------------------------------------

module adgaquad_mod

  use carma_precision_mod
  use adgaquad_types_mod

  implicit none
  
  private
 
  public :: dqag
  public :: dqage
  public :: dqagi
  public :: dqagie
  
  contains 

  !!***begin prologue  dqag
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)  
  !!***category no.  h2a1a1
  !!***keywords  automatic integrator, general-purpose,
  !!             integrand examinator, globally adaptive,
  !!             gauss-kronrod
  !!***author  piessens,robert,appl. math. & progr. div - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !!***purpose  the routine calculates an approximation result to a given
  !!            definite integral i = integral of f over (a,b),
  !!            hopefully satisfying following claim for accuracy
  !!            abs(i-result)le.max(epsabs,epsrel*abs(i)).
  !!***description
  !!
  !!        computation of a definite integral
  !!        standard fortran subroutine
  !!        double precision version
  !!
  !!            fx     - double precision
  !!                     function subprogam defining the integrand
  !!                     function f(x). the actual name for f needs to be
  !!                     declared e x t e r n a l in the driver program.
  !!
  !!            fx_vars- structure containing variables need for integration
  !!                     specific to fractal meanfield scattering code
  !!
  !!            a      - double precision
  !!                     lower limit of integration
  !!
  !!            b      - double precision
  !!                     upper limit of integration
  !!
  !!            epsabs - double precision
  !!                     absolute accoracy requested
  !!            epsrel - double precision
  !!                     relative accuracy requested
  !!                     if  epsabs.le.0
  !!                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
  !!                     the routine will end with ier = 6.
  !!
  !!            key    - integer
  !!                     key for choice of local integration rule
  !!                     a gauss-kronrod pair is used with
  !!                       7 - 15 points if key.lt.2,
  !!                      10 - 21 points if key = 2,
  !!                      15 - 31 points if key = 3,
  !!                      20 - 41 points if key = 4,
  !!                      25 - 51 points if key = 5,
  !!                      30 - 61 points if key.gt.5.
  !!
  !!         on return
  !!            result - double precision
  !!                     approximation to the integral
  !!
  !!            abserr - double precision
  !!                     estimate of the modulus of the absolute error,
  !!                     which should equal or exceed abs(i-result)
  !!
  !!            neval  - integer
  !!                     number of integrand evaluations
  !!
  !!            ier    - integer
  !!                     ier = 0 normal and reliable termination of the
  !!                             routine. it is assumed that the requested
  !!                             accuracy has been achieved.
  !!                     ier.gt.0 abnormal termination of the routine
  !!                             the estimates for result and error are
  !!                             less reliable. it is assumed that the
  !!                             requested accuracy has not been achieved.
  !!                      error messages
  !!                     ier = 1 maximum number of subdivisions allowed
  !!                             has been achieved. one can allow more
  !!                             subdivisions by increasing the value of
  !!                             limit (and taking the according dimension
  !!                             adjustments into account). however, if
  !!                             this yield no improvement it is advised
  !!                             to analyze the integrand in order to
  !!                             determine the integration difficulaties.
  !!                             if the position of a local difficulty can
  !!                             be determined (i.e.singularity,
  !!                             discontinuity within the interval) one
  !!                             will probably gain from splitting up the
  !!                             interval at this point and calling the
  !!                             integrator on the subranges. if possible,
  !!                             an appropriate special-purpose integrator
  !!                             should be used which is designed for
  !!                             handling the type of difficulty involved.
  !!                         = 2 the occurrence of roundoff error is
  !!                             detected, which prevents the requested
  !!                             tolerance from being achieved.
  !!                         = 3 extremely bad integrand behaviour occurs
  !!                             at some points of the integration
  !!                             interval.
  !!                         = 6 the input is invalid, because
  !!                             (epsabs.le.0 and
  !!                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
  !!                             or limit.lt.1 or lenw.lt.limit*4.
  !!                             result, abserr, neval, last are set
  !!                             to zero.
  !!                             except when lenw is invalid, iwork(1),
  !!                             work(limit*2+1) and work(limit*3+1) are
  !!                             set to zero, work(1) is set to a and
  !!                             work(limit+1) to b.
  !!		             = 9 failure in sd1mach determining machine parameters
  !!
  !!         dimensioning parameters
  !!            limit - integer
  !!                    dimensioning parameter for iwork
  !!                    limit determines the maximum number of subintervals
  !!                    in the partition of the given integration interval
  !!                    (a,b), limit.ge.1.
  !!                    if limit.lt.1, the routine will end with ier = 6.
  !!
  !!            lenw  - integer
  !!                    dimensioning parameter for work
  !!                    lenw must be at least limit*4.
  !!                    if lenw.lt.limit*4, the routine will end with
  !!                    ier = 6.
  !!
  !!            last  - integer
  !!                    on return, last equals the number of subintervals
  !!                    produced in the subdiviosion process, which
  !!                    determines the number of significant elements
  !!                    actually in the work arrays.
  !!
  !!         work arrays
  !!            iwork - integer
  !!                    vector of dimension at least limit, the first k
  !!                    elements of which contain pointers to the error
  !!                    estimates over the subintervals, such that
  !!                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
  !!                    form a decreasing sequence with k = last if
  !!                    last.le.(limit/2+2), and k = limit+1-last otherwise
  !!
  !!            work  - double precision
  !!                    vector of dimension at least lenw
  !!                    on return
  !!                    work(1), ..., work(last) contain the left end
  !!                    points of the subintervals in the partition of
  !!                     (a,b),
  !!                    work(limit+1), ..., work(limit+last) contain the
  !!                     right end points,
  !!                    work(limit*2+1), ..., work(limit*2+last) contain
  !!                     the integral approximations over the subintervals,
  !!                    work(limit*3+1), ..., work(limit*3+last) contain
  !!                     the error estimates.
  !!
  !!***references  (none)
  !!***routines called  dqage,xerror
  !!***end prologue  dqag
  subroutine dqag(fx,fx_vars,a,b,epsabs,epsrel,key,result,abserr,neval,ier, &
                  limit,lenw,last,iwork,work)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: a
    real(kind=f) :: b
    real(kind=f) :: epsabs
    real(kind=f) :: epsrel
    integer :: key
    real(kind=f) :: result
    real(kind=f) :: abserr
    integer :: neval
    integer :: ier
    integer :: limit
    integer :: lenw
    integer :: last
    integer :: iwork(limit)
    real(kind=f) :: work(lenw)

    ! Local declarations
    integer :: lvl,l1,l2,l3

    !  check validity of lenw.
    !
    !***first executable statement  dqag
    ier = 6
    neval = 0
    last = 0
    result = 0.0_f
    abserr = 0.0_f
    if(limit.lt.1.or.lenw.lt.limit*4) go to 10

    !  prepare call for dqage.

    l1 = limit+1
    l2 = limit+l1
    l3 = limit+l2

    call dqage(fx,fx_vars,a,b,epsabs,epsrel,key,limit,result,abserr,neval, &
               ier,work(1),work(l1),work(l2),work(l3),iwork,last)
 
    !  call error handler if necessary.
 
    lvl = 0
10  if(ier.eq.6) lvl = 1
    if(ier.ne.0) then
      write(*,*) "ERROR: abnormal return from dqag"
      write(*,*) "       ifail=",ier,"  level=",lvl
    endif
    return
  end subroutine dqag



  !!***begin prologue  dqage
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)
  !!***category no.  h2a1a1
  !!***keywords  automatic integrator, general-purpose,
  !!             integrand examinator, globally adaptive,
  !!             gauss-kronrod
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !!***purpose  the routine calculates an approximation result to a given
  !!            definite integral   i = integral of f over (a,b),
  !!            hopefully satisfying following claim for accuracy
  !!            abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
  !!***description
  !!
  !!        computation of a definite integral
  !!        standard fortran subroutine
  !!        double precision version
  !!
  !!        parameters
  !!         on entry
  !!            fx     - double precision
  !!                     function subprogram defining the integrand
  !!                     function f(x). the actual name for f needs to be
  !!                     declared e x t e r n a l in the driver program.
  !!
  !!            fx_vars- structure containing variables need for integration
  !!                     specific to fractal meanfield scattering code
  !!
  !!            a      - double precision
  !!                     lower limit of integration
  !!
  !!            b      - double precision
  !!                     upper limit of integration
  !!
  !!            epsabs - double precision
  !!                     absolute accuracy requested
  !!            epsrel - double precision
  !!                     relative accuracy requested
  !!                     if  epsabs.le.0
  !!                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
  !!                     the routine will end with ier = 6.
  !!
  !!            key    - integer
  !!                     key for choice of local integration rule
  !!                     a gauss-kronrod pair is used with
  !!                          7 - 15 points if key.lt.2,
  !!                         10 - 21 points if key = 2,
  !!                         15 - 31 points if key = 3,
  !!                         20 - 41 points if key = 4,
  !!                         25 - 51 points if key = 5,
  !!                         30 - 61 points if key.gt.5.
  !!
  !!            limit  - integer
  !!                     gives an upperbound on the number of subintervals
  !!                     in the partition of (a,b), limit.ge.1.
  !!
  !!         on return
  !!            result - double precision
  !!                     approximation to the integral
  !!
  !!            abserr - double precision
  !!                     estimate of the modulus of the absolute error,
  !!                     which should equal or exceed abs(i-result)
  !!
  !!            neval  - integer
  !!                     number of integrand evaluations
  !!
  !!            ier    - integer
  !!                     ier = 0 normal and reliable termination of the
  !!                             routine. it is assumed that the requested
  !!                             accuracy has been achieved.
  !!                     ier.gt.0 abnormal termination of the routine
  !!                             the estimates for result and error are
  !!                             less reliable. it is assumed that the
  !!                             requested accuracy has not been achieved.
  !!            error messages
  !!                     ier = 1 maximum number of subdivisions allowed
  !!                             has been achieved. one can allow more
  !!                             subdivisions by increasing the value
  !!                             of limit.
  !!                             however, if this yields no improvement it
  !!                             is rather advised to analyze the integrand
  !!                             in order to determine the integration
  !!                             difficulties. if the position of a local
  !!                             difficulty can be determined(e.g.
  !!                             singularity, discontinuity within the
  !!                             interval) one will probably gain from
  !!                             splitting up the interval at this point
  !!                             and calling the integrator on the
  !!                             subranges. if possible, an appropriate
  !!                             special-purpose integrator should be used
  !!                             which is designed for handling the type of
  !!                             difficulty involved.
  !!                         = 2 the occurrence of roundoff error is
  !!                             detected, which prevents the requested
  !!                             tolerance from being achieved.
  !!                         = 3 extremely bad integrand behaviour occurs
  !!                             at some points of the integration
  !!                             interval.
  !!                         = 6 the input is invalid, because
  !!                             (epsabs.le.0 and
  !!                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
  !!                             result, abserr, neval, last, rlist(1) ,
  !!                             elist(1) and iord(1) are set to zero.
  !!                             alist(1) and blist(1) are set to a and b
  !!                             respectively.
  !!                         = 9 failure in sd1mach determining machine parameters
  !!
  !!            alist   - double precision
  !!                      vector of dimension at least limit, the first
  !!                       last  elements of which are the left
  !!                      end points of the subintervals in the partition
  !!                      of the given integration range (a,b)
  !!
  !!            blist   - double precision
  !!                      vector of dimension at least limit, the first
  !!                       last  elements of which are the right
  !!                      end points of the subintervals in the partition
  !!                      of the given integration range (a,b)
  !!
  !!            rlist   - double precision
  !!                      vector of dimension at least limit, the first
  !!                       last  elements of which are the
  !!                      integral approximations on the subintervals
  !!
  !!            elist   - double precision
  !!                      vector of dimension at least limit, the first
  !!                       last  elements of which are the moduli of the
  !!                       absolute error estimates on the subintervals
  !!
  !!            iord    - integer
  !!                      vector of dimension at least limit, the first k
  !!                      elements of which are pointers to the
  !!                      error estimates over the subintervals,
  !!                      such that elist(iord(1)), ...,
  !!                      elist(iord(k)) form a decreasing sequence,
  !!                      with k = last if last.le.(limit/2+2), and
  !!                      k = limit+1-last otherwise
  !!
  !!            last    - integer
  !!                      number of subintervals actually produced in the
  !!                      subdivision process
  !!
  !!***references  (none)
  !!***routines called  sd1mach,dqk15,dqk21,dqk31,
  !!                    dqk41,dqk51,dqk61,dqpsrt
  !!***end prologue  dqage
  subroutine dqage(fx,fx_vars,a,b,epsabs,epsrel,key,limit,result,abserr, &
                   neval,ier,alist,blist,rlist,elist,iord,last)
  
    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: a
    real(kind=f) :: b
    real(kind=f) :: epsabs
    real(kind=f) :: epsrel
    integer :: limit
    integer :: key
    real(kind=f) :: result
    real(kind=f) :: abserr
    integer :: neval
    integer :: ier
    real(kind=f) :: alist(limit)
    real(kind=f) :: blist(limit)
    real(kind=f) :: rlist(limit)
    real(kind=f) :: elist(limit)
    integer :: iord(limit)
    integer :: last

    ! Local declarations
    real(kind=f) :: area, area1, area12, area2, a1, a2
    real(kind=f) :: b1, b2, dabs, defabs, defab1, defab2, dmax1, epmach   
    real(kind=f) :: errbnd,errmax,error1,error2,erro12,errsum,resabs,uflow
    integer :: iroff1,iroff2,k,keyf,maxerr,nrmax


    !            list of major variables
    !            -----------------------
    !
    !           alist     - list of left end points of all subintervals
    !                       considered up to now
    !           blist     - list of right end points of all subintervals
    !                       considered up to now
    !           rlist(i)  - approximation to the integral over
    !                      (alist(i),blist(i))
    !           elist(i)  - error estimate applying to rlist(i)
    !           maxerr    - pointer to the interval with largest
    !                       error estimate
    !           errmax    - elist(maxerr)
    !           area      - sum of the integrals over the subintervals
    !           errsum    - sum of the errors over the subintervals
    !           errbnd    - requested accuracy max(epsabs,epsrel*
    !                       abs(result))
    !           *****1    - variable for the left subinterval
    !           *****2    - variable for the right subinterval
    !           last      - index for subdivision
    !
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach  is the largest relative spacing.
    !           uflow  is the smallest positive magnitude.
    !
    !***first executable statement  dqage
    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    call sd1mach(1,uflow,ier)
    if(ier.eq.9) return
 
    !epmach = d1mach(4)
    !uflow = d1mach(1)

    !           test on validity of parameters
    !           ------------------------------
    !
    ier = 0
    neval = 0
    last = 0
    result = 0.0_f
    abserr = 0.0_f
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0_f
    elist(1) = 0.0_f
    iord(1) = 0
    if(epsabs.le.0.0_f.and.epsrel.lt.dmax1(0.5e2_f*epmach,0.5e-28_f)) ier = 6
    if(ier.eq.6) go to 999

    !           first approximation to the integral
    !           -----------------------------------
    !
    keyf = key
    if(key.le.0) keyf = 1
    if(key.ge.7) keyf = 6
    neval = 0
    if(keyf.eq.1) call dqk15(fx,fx_vars,a,b,result,abserr,defabs,resabs,ier)
    if(ier.eq.9) return
    if(keyf.eq.2) call dqk21(fx,fx_vars,a,b,result,abserr,defabs,resabs,ier)
    if(ier.eq.9) return
    if(keyf.eq.3) call dqk31(fx,fx_vars,a,b,result,abserr,defabs,resabs,ier)
    if(ier.eq.9) return
    if(keyf.eq.4) call dqk41(fx,fx_vars,a,b,result,abserr,defabs,resabs,ier)
    if(ier.eq.9) return
    if(keyf.eq.5) call dqk51(fx,fx_vars,a,b,result,abserr,defabs,resabs,ier)
    if(ier.eq.9) return
    if(keyf.eq.6) call dqk61(fx,fx_vars,a,b,result,abserr,defabs,resabs,ier)
    if(ier.eq.9) return
    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1
    !
    !           test on accuracy.
    !
    errbnd = dmax1(epsabs,epsrel*dabs(result))
    if(abserr.le.0.5e2_f*epmach*defabs.and.abserr.gt.errbnd) ier = 2
    if(limit.eq.1) ier = 1
    if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.abserr.eq.0.0d+00) go to 60
    !
    !           initialization
    !           --------------
    !
    errmax = abserr
    maxerr = 1
    area = result
    errsum = abserr
    nrmax = 1
    iroff1 = 0
    iroff2 = 0
    !
    !           main do-loop
    !           ------------
    !
    do last = 2,limit
      !
      !           bisect the subinterval with the largest error estimate.
      !
      a1 = alist(maxerr)
      b1 = 0.5_f*(alist(maxerr)+blist(maxerr))
      a2 = b1
      b2 = blist(maxerr)
      if(keyf.eq.1) call dqk15(fx,fx_vars,a1,b1,area1,error1,resabs,defab1,ier)
      if(ier.eq.9) return
      if(keyf.eq.2) call dqk21(fx,fx_vars,a1,b1,area1,error1,resabs,defab1,ier)
      if(ier.eq.9) return
      if(keyf.eq.3) call dqk31(fx,fx_vars,a1,b1,area1,error1,resabs,defab1,ier)
      if(ier.eq.9) return
      if(keyf.eq.4) call dqk41(fx,fx_vars,a1,b1,area1,error1,resabs,defab1,ier)
      if(ier.eq.9) return
      if(keyf.eq.5) call dqk51(fx,fx_vars,a1,b1,area1,error1,resabs,defab1,ier)
      if(ier.eq.9) return
      if(keyf.eq.6) call dqk61(fx,fx_vars,a1,b1,area1,error1,resabs,defab1,ier)
      if(ier.eq.9) return
      if(keyf.eq.1) call dqk15(fx,fx_vars,a2,b2,area2,error2,resabs,defab2,ier)
      if(ier.eq.9) return
      if(keyf.eq.2) call dqk21(fx,fx_vars,a2,b2,area2,error2,resabs,defab2,ier)
      if(ier.eq.9) return
      if(keyf.eq.3) call dqk31(fx,fx_vars,a2,b2,area2,error2,resabs,defab2,ier)
      if(ier.eq.9) return
      if(keyf.eq.4) call dqk41(fx,fx_vars,a2,b2,area2,error2,resabs,defab2,ier)
      if(ier.eq.9) return
      if(keyf.eq.5) call dqk51(fx,fx_vars,a2,b2,area2,error2,resabs,defab2,ier)
      if(ier.eq.9) return
      if(keyf.eq.6) call dqk61(fx,fx_vars,a2,b2,area2,error2,resabs,defab2,ier)
      if(ier.eq.9) return
      !           improve previous approximations to integral
      !           and error and test for accuracy.
      !
      !       neval = neval+1
      area12 = area1+area2
      erro12 = error1+error2
      errsum = errsum+erro12-errmax
      area = area+area12-rlist(maxerr)
      if(defab1.eq.error1.or.defab2.eq.error2) go to 5
      if(dabs(rlist(maxerr)-area12).le.0.1e-4_f*dabs(area12).and.erro12.ge.0.99_f*errmax) iroff1 = iroff1+1
      if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
    5 rlist(maxerr) = area1
      rlist(last) = area2
      errbnd = dmax1(epsabs,epsrel*dabs(area))
      if(errsum.le.errbnd) go to 8
      !
      !           test for roundoff error and eventually set error flag.
      !
      if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
      !
      !           set error flag in the case that the number of subintervals
      !           equals limit.
      !
      if(last.eq.limit) ier = 1
      !
      !           set error flag in the case of bad integrand behaviour
      !           at a point of the integration range.
      !
      if(dmax1(dabs(a1),dabs(b2)).le.(0.1e1_f+0.1e3_f*epmach)*(dabs(a2)+0.1e4_f*uflow)) ier = 3
      !
      !           append the newly-created intervals to the list.
      !
    8 if(error2.gt.error1) go to 10
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
      go to 20
   10 alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
      !
      !           call subroutine dqpsrt to maintain the descending ordering
      !           in the list of error estimates and select the subinterval
      !           with the largest error estimate (to be bisected next).
      !
   20 call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
      ! ***jump out of do-loop
      if(ier.ne.0.or.errsum.le.errbnd) go to 40
    end do
    !
    !           compute final result.
    !           ---------------------
    !
 40 result = 0.0_f
    do k=1,last
      result = result+rlist(k)
    end do
    abserr = errsum
 60 if(keyf.ne.1) neval = (10*keyf+1)*(2*neval+1)
    if(keyf.eq.1) neval = 30*neval+15
999 return
  end subroutine dqage


  !!***begin prologue  dqagi
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)
  !!***category no.  h2a3a1,h2a4a1
  !!***keywords  automatic integrator, infinite intervals,
  !!             general-purpose, transformation, extrapolation,
  !!             globally adaptive
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div. -k.u.leuven
  !!***purpose  the routine calculates an approximation result to a given
  !!            integral   i = integral of f over (bound,+infinity)
  !!            or i = integral of f over (-infinity,bound)
  !!            or i = integral of f over (-infinity,+infinity)
  !!            hopefully satisfying following claim for accuracy
  !!            abs(i-result).le.max(epsabs,epsrel*abs(i)).
  !!***description
  !!
  !!        integration over infinite intervals
  !!        standard fortran subroutine
  !!
  !!        parameters
  !!         on entry
  !!            fx     - double precision
  !!                     function subprogram defining the integrand
  !!                     function f(x). the actual name for f needs to be
  !!                     declared e x t e r n a l in the driver program.
  !!
  !!            fx_vars- structure containing variables need for integration
  !!                     specific to fractal meanfield scattering code
  !!
  !!            bound  - double precision
  !!                     finite bound of integration range
  !!                     (has no meaning if interval is doubly-infinite)
  !!
  !!            inf    - integer
  !!                     indicating the kind of integration range involved
  !!                     inf = 1 corresponds to  (bound,+infinity),
  !!                     inf = -1            to  (-infinity,bound),
  !!                     inf = 2             to (-infinity,+infinity).
  !!
  !!            epsabs - double precision
  !!                     absolute accuracy requested
  !!            epsrel - double precision
  !!                     relative accuracy requested
  !!                     if  epsabs.le.0
  !!                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
  !!                     the routine will end with ier = 6.
  !!
  !!
  !!         on return
  !!            result - double precision
  !!                     approximation to the integral
  !!
  !!            abserr - double precision
  !!                     estimate of the modulus of the absolute error,
  !!                     which should equal or exceed abs(i-result)
  !!
  !!            neval  - integer
  !!                     number of integrand evaluations
  !!
  !!            ier    - integer
  !!                     ier = 0 normal and reliable termination of the
  !!                             routine. it is assumed that the requested
  !!                             accuracy has been achieved.
  !!                   - ier.gt.0 abnormal termination of the routine. the
  !!                             estimates for result and error are less
  !!                             reliable. it is assumed that the requested
  !!                             accuracy has not been achieved.
  !!            error messages
  !!                     ier = 1 maximum number of subdivisions allowed
  !!                             has been achieved. one can allow more
  !!                             subdivisions by increasing the value of
  !!                             limit (and taking the according dimension
  !!                             adjustments into account). however, if
  !!                             this yields no improvement it is advised
  !!                             to analyze the integrand in order to
  !!                             determine the integration difficulties. if
  !!                             the position of a local difficulty can be
  !!                             determined (e.g. singularity,
  !!                             discontinuity within the interval) one
  !!                             will probably gain from splitting up the
  !!                             interval at this point and calling the
  !!                             integrator on the subranges. if possible,
  !!                             an appropriate special-purpose integrator
  !!                             should be used, which is designed for
  !!                             handling the type of difficulty involved.
  !!                         = 2 the occurrence of roundoff error is
  !!                             detected, which prevents the requested
  !!                             tolerance from being achieved.
  !!                             the error may be under-estimated.
  !!                         = 3 extremely bad integrand behaviour occurs
  !!                             at some points of the integration
  !!                             interval.
  !!                         = 4 the algorithm does not converge.
  !!                             roundoff error is detected in the
  !!                             extrapolation table.
  !!                             it is assumed that the requested tolerance
  !!                             cannot be achieved, and that the returned
  !!                             result is the best which can be obtained.
  !!                         = 5 the integral is probably divergent, or
  !!                             slowly convergent. it must be noted that
  !!                             divergence can occur with any other value
  !!                             of ier.
  !!                         = 6 the input is invalid, because
  !!                             (epsabs.le.0 and
  !!                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
  !!                              or limit.lt.1 or leniw.lt.limit*4.
  !!                             result, abserr, neval, last are set to
  !!                             zero. exept when limit or leniw is
  !!                             invalid, iwork(1), work(limit*2+1) and
  !!                             work(limit*3+1) are set to zero, work(1)
  !!                             is set to a and work(limit+1) to b.
  !!		             = 9 failure in sd1mach determining machine parameters
  !!
  !!          dimensioning parameters
  !!            limit - integer
  !!                    dimensioning parameter for iwork
  !!                    limit determines the maximum number of subintervals
  !!                    in the partition of the given integration interval
  !!                    (a,b), limit.ge.1.
  !!                    if limit.lt.1, the routine will end with ier = 6.
  !!
  !!            lenw  - integer
  !!                    dimensioning parameter for work
  !!                    lenw must be at least limit*4.
  !!                    if lenw.lt.limit*4, the routine will end
  !!                    with ier = 6.
  !!
  !!            last  - integer
  !!                    on return, last equals the number of subintervals
  !!                    produced in the subdivision process, which
  !!                    determines the number of significant elements
  !!                    actually in the work arrays.
  !!
  !!         work arrays
  !!            iwork - integer
  !!                    vector of dimension at least limit, the first
  !!                    k elements of which contain pointers
  !!                    to the error estimates over the subintervals,
  !!                    such that work(limit*3+iwork(1)),... ,
  !!                    work(limit*3+iwork(k)) form a decreasing
  !!                    sequence, with k = last if last.le.(limit/2+2), and
  !!                    k = limit+1-last otherwise
  !!
  !!            work  - double precision
  !!                    vector of dimension at least lenw
  !!                    on return
  !!                    work(1), ..., work(last) contain the left
  !!                     end points of the subintervals in the
  !!                     partition of (a,b),
  !!                    work(limit+1), ..., work(limit+last) contain
  !!                     the right end points,
  !!                    work(limit*2+1), ...,work(limit*2+last) contain the
  !!                     integral approximations over the subintervals,
  !!                    work(limit*3+1), ..., work(limit*3)
  !!                     contain the error estimates.
  !!***references  (none)
  !!***routines called  dqagie,xerror
  !!***end prologue  dqagi
  !!
  subroutine dqagi(fx,fx_vars,bound,inf,epsabs,epsrel,result,abserr,neval, &
                   ier,limit,lenw,last,iwork,work)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: bound
    integer :: inf
    real(kind=f) :: epsabs 
    real(kind=f) :: epsrel
    real(kind=f) :: result
    real(kind=f) :: abserr
    integer :: neval
    integer :: ier
    integer :: limit  
    integer :: lenw  
    integer ::  last
    integer ::  iwork(limit)
    real(kind=f) :: work(lenw)

    ! Local declarations
    integer lvl,l1,l2,l3

    !
    !         check validity of limit and lenw.
    !
    !***first executable statement  dqagi
    ier = 6
    neval = 0
    last = 0
    result = 0.0_f
    abserr = 0.0_f
    if(limit.lt.1.or.lenw.lt.limit*4) go to 10
    !
    !         prepare call for dqagie.
    !
    l1 = limit+1
    l2 = limit+l1
    l3 = limit+l2

    call dqagie(fx,fx_vars,bound,inf,epsabs,epsrel,limit,result,abserr, &
                neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
    !
    !         call error handler if necessary.
    !
    lvl = 0
10  if(ier.eq.6) lvl = 1
    if(ier.ne.0) then
      write(*,*) "ERROR: abnormal return from dqagi"
      write(*,*) "       ifail=",ier,"  level=",lvl
    endif
    return
  end subroutine dqagi


  !!***begin prologue  dqagie
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)
  !!***category no.  h2a3a1,h2a4a1
  !!***keywords  automatic integrator, infinite intervals,
  !!             general-purpose, transformation, extrapolation,
  !!             globally adaptive
  !!***author  piessens,robert,appl. math & progr. div - k.u.leuven
  !!           de doncker,elise,appl. math & progr. div - k.u.leuven
  !!***purpose  the routine calculates an approximation result to a given
  !!            integral   i = integral of f over (bound,+infinity)
  !!            or i = integral of f over (-infinity,bound)
  !!            or i = integral of f over (-infinity,+infinity),
  !!            hopefully satisfying following claim for accuracy
  !!            abs(i-result).le.max(epsabs,epsrel*abs(i))
  !!***description
  !!
  !! integration over infinite intervals
  !! standard fortran subroutine
  !!
  !!            fx     - double precision
  !!                     function subprogram defining the integrand
  !!                     function f(x). the actual name for f needs to be
  !!                     declared e x t e r n a l in the driver program.
  !!
  !!            fx_vars- structure containing variables need for integration
  !!                     specific to fractal meanfield scattering code
  !!
  !!            bound  - double precision
  !!                     finite bound of integration range
  !!                     (has no meaning if interval is doubly-infinite)
  !!
  !!            inf    - double precision
  !!                     indicating the kind of integration range involved
  !!                     inf = 1 corresponds to  (bound,+infinity),
  !!                     inf = -1            to  (-infinity,bound),
  !!                     inf = 2             to (-infinity,+infinity).
  !!
  !!            epsabs - double precision
  !!                     absolute accuracy requested
  !!            epsrel - double precision
  !!                     relative accuracy requested
  !!                     if  epsabs.le.0
  !!                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
  !!                     the routine will end with ier = 6.
  !!
  !!            limit  - integer
  !!                     gives an upper bound on the number of subintervals
  !!                     in the partition of (a,b), limit.ge.1
  !!
  !!         on return
  !!           result - double precision
  !!                     approximation to the integral
  !!
  !!            abserr - double precision
  !!                     estimate of the modulus of the absolute error,
  !!                     which should equal or exceed abs(i-result)
  !!
  !!            neval  - integer
  !!                     number of integrand evaluations
  !!
  !!            ier    - integer
  !!                     ier = 0 normal and reliable termination of the
  !!                             routine. it is assumed that the requested
  !!                             accuracy has been achieved.
  !!                   - ier.gt.0 abnormal termination of the routine. the
  !!                             estimates for result and error are less
  !!                             reliable. it is assumed that the requested
  !!                             accuracy has not been achieved.
  !!            error messages
  !!                     ier = 1 maximum number of subdivisions allowed
  !!                             has been achieved. one can allow more
  !!                             subdivisions by increasing the value of
  !!                             limit (and taking the according dimension
  !!                             adjustments into account). however,if
  !!                             this yields no improvement it is advised
  !!                             to analyze the integrand in order to
  !!                             determine the integration difficulties.
  !!                             if the position of a local difficulty can
  !!                             be determined (e.g. singularity,
  !!                             discontinuity within the interval) one
  !!                             will probably gain from splitting up the
  !!                             interval at this point and calling the
  !!                             integrator on the subranges. if possible,
  !!                             an appropriate special-purpose integrator
  !!                             should be used, which is designed for
  !!                             handling the type of difficulty involved.
  !!                         = 2 the occurrence of roundoff error is
  !!                             detected, which prevents the requested
  !!                             tolerance from being achieved.
  !!                             the error may be under-estimated.
  !!                         = 3 extremely bad integrand behaviour occurs
  !!                             at some points of the integration
  !!                             interval.
  !!                         = 4 the algorithm does not converge.
  !!                             roundoff error is detected in the
  !!                             extrapolation table.
  !!                             it is assumed that the requested tolerance
  !!                             cannot be achieved, and that the returned
  !!                             result is the best which can be obtained.
  !!                         = 5 the integral is probably divergent, or
  !!                             slowly convergent. it must be noted that
  !!                             divergence can occur with any other value
  !!                             of ier.
  !!                         = 6 the input is invalid, because
  !!                             (epsabs.le.0 and
  !!                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
  !!                             result, abserr, neval, last, rlist(1),
  !!                             elist(1) and iord(1) are set to zero.
  !!                             alist(1) and blist(1) are set to 0
  !!                             and 1 respectively.
  !!		             = 9 failure in sd1mach determining machine parameters
  !!
  !!            alist  - double precision
  !!                     vector of dimension at least limit, the first
  !!                      last  elements of which are the left
  !!                     end points of the subintervals in the partition
  !!                     of the transformed integration range (0,1).
  !!
  !!            blist  - double precision
  !!                     vector of dimension at least limit, the first
  !!                      last  elements of which are the right
  !!                     end points of the subintervals in the partition
  !!                     of the transformed integration range (0,1).
  !!
  !!            rlist  - double precision
  !!                     vector of dimension at least limit, the first
  !!                      last  elements of which are the integral
  !!                     approximations on the subintervals
  !!
  !!            elist  - double precision
  !!                     vector of dimension at least limit,  the first
  !!                     last elements of which are the moduli of the
  !!                     absolute error estimates on the subintervals
  !!
  !!            iord   - integer
  !!                     vector of dimension limit, the first k
  !!                     elements of which are pointers to the
  !!                     error estimates over the subintervals,
  !!                     such that elist(iord(1)), ..., elist(iord(k))
  !!                     form a decreasing sequence, with k = last
  !!                     if last.le.(limit/2+2), and k = limit+1-last
  !!                     otherwise
  !!
  !!            last   - integer
  !!                     number of subintervals actually produced
  !!                     in the subdivision process
  !!
  !!***references  (none)
  !!***routines called  sd1mach,dqelg,dqk15i,dqpsrt
  !!***end prologue  dqagie
  subroutine dqagie(fx,fx_vars,bound,inf,epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: bound
    integer :: inf
    real(kind=f) :: epsabs
    real(kind=f) :: epsrel
    integer :: limit
    real(kind=f) :: result
    real(kind=f) :: abserr
    integer :: neval
    integer :: ier
    real(kind=f) :: alist(limit)
    real(kind=f) :: blist(limit)
    real(kind=f) :: rlist(limit)
    real(kind=f) :: elist(limit)
    integer :: iord(limit)
    integer :: last
    
    ! Local declartions
    real(kind=f) :: abseps, area, area1, area12, area2
    real(kind=f) :: a1, a2, b1, b2, correc
    real(kind=f) :: defabs, defab1, defab2
    real(kind=f) :: dmax1, dres, epmach, erlarg, erlast
    real(kind=f) :: errbnd, errmax, error1, error2, erro12, errsum
    real(kind=f) :: ertest, oflow, resabs, reseps, res3la(3), rlist2(52)
    real(kind=f) :: small, uflow, boun
      
    integer :: id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn
    integer :: ktmin, maxerr, nres, nrmax, numrl2
    logical :: extrap, noext


    !
    !            the dimension of rlist2 is determined by the value of
    !            limexp in subroutine dqelg.
    !
    !            list of major variables
    !            -----------------------
    !
    !           alist     - list of left end points of all subintervals
    !                       considered up to now
    !           blist     - list of right end points of all subintervals
    !                       considered up to now
    !           rlist(i)  - approximation to the integral over
    !                       (alist(i),blist(i))
    !           rlist2    - array of dimension at least (limexp+2),
    !                       containing the part of the epsilon table
    !                       wich is still needed for further computations
    !           elist(i)  - error estimate applying to rlist(i)
    !           maxerr    - pointer to the interval with largest error
    !                       estimate
    !           errmax    - elist(maxerr)
    !           erlast    - error on the interval currently subdivided
    !                       (before that subdivision has taken place)
    !           area      - sum of the integrals over the subintervals
    !           errsum    - sum of the errors over the subintervals
    !           errbnd    - requested accuracy max(epsabs,epsrel*
    !                       abs(result))
    !           *****1    - variable for the left subinterval
    !           *****2    - variable for the right subinterval
    !           last      - index for subdivision
    !           nres      - number of calls to the extrapolation routine
    !           numrl2    - number of elements currently in rlist2. if an
    !                       appropriate approximation to the compounded
    !                       integral has been obtained, it is put in
    !                       rlist2(numrl2) after numrl2 has been increased
    !                       by one.
    !           small     - length of the smallest interval considered up
    !                       to now, multiplied by 1.5
    !           erlarg    - sum of the errors over the intervals larger
    !                       than the smallest interval considered up to now
    !           extrap    - logical variable denoting that the routine
    !                       is attempting to perform extrapolation. i.e.
    !                       before subdividing the smallest interval we
    !                       try to decrease the value of erlarg.
    !           noext     - logical variable denoting that extrapolation
    !                       is no longer allowed (true-value)
    !
    !            machine dependent constants
    !            ---------------------------
    !
    !           epmach is the largest relative spacing.
    !           uflow is the smallest positive magnitude.
    !           oflow is the largest positive magnitude.
    !
    !***first executable statement  dqagie

    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    
    !epmach = d1mach(4)
    !
    !           test on validity of parameters
    !           -----------------------------
    !
    ier = 0
    neval = 0
    last = 0
    result = 0.0_f
    abserr = 0.0_f
    alist(1) = 0.0_f
    blist(1) = 0.1e1_f
    rlist(1) = 0.0_f
    elist(1) = 0.0_f
    iord(1) = 0
    if(epsabs.le.0.0_f.and.epsrel.lt.dmax1(0.5e2_f*epmach,0.5e-28_f)) ier = 6
    if(ier.eq.6) go to 999
    !
    !
    !           first approximation to the integral
    !           -----------------------------------
    !
    !           determine the interval to be mapped onto (0,1).
    !           if inf = 2 the integral is computed as i = i1+i2, where
    !           i1 = integral of f over (-infinity,0),
    !           i2 = integral of f over (0,+infinity).
    !
    boun = bound
    if(inf.eq.2) boun = 0.0_f
    call dqk15i(fx,fx_vars,boun,inf,0.0_f,0.1e1_f,result,abserr,defabs,resabs,ier)
    if(ier.eq.9) return
    !
    !           test on accuracy
    !
    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1
    dres = dabs(result)
    errbnd = dmax1(epsabs,epsrel*dres)
    if(abserr.le.1.0e2_f*epmach*defabs.and.abserr.gt.errbnd) ier = 2
    if(limit.eq.1) ier = 1
    if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.abserr.eq.0.0_f) go to 130
    !
    !           initialization
    !           --------------
    !
    call sd1mach(1,uflow,ier)
    if(ier.eq.9) return
    call sd1mach(2,oflow,ier)
    if(ier.eq.9) return

    !uflow = d1mach(1)
    !oflow = d1mach(2)
    rlist2(1) = result
    errmax = abserr
    maxerr = 1
    area = result
    errsum = abserr
    abserr = oflow
    nrmax = 1
    nres = 0
    ktmin = 0
    numrl2 = 2
    extrap = .false.
    noext = .false.
    ierro = 0
    iroff1 = 0
    iroff2 = 0
    iroff3 = 0
    ksgn = -1
    if(dres.ge.(0.1e1_f-0.5e2_f*epmach)*defabs) ksgn = 1
    !
    !           main do-loop
    !           ------------
    !
    do 90 last = 2,limit
      !
      !           bisect the subinterval with nrmax-th largest error estimate.
      !
      a1 = alist(maxerr)
      b1 = 0.5_f*(alist(maxerr)+blist(maxerr))
      a2 = b1
      b2 = blist(maxerr)
      erlast = errmax
      call dqk15i(fx,fx_vars,boun,inf,a1,b1,area1,error1,resabs,defab1,ier)
      if(ier.eq.9) return
      call dqk15i(fx,fx_vars,boun,inf,a2,b2,area2,error2,resabs,defab2,ier)
      if(ier.eq.9) return
      !
      !           improve previous approximations to integral
      !           and error and test for accuracy.
      !
      area12 = area1+area2
      erro12 = error1+error2
      errsum = errsum+erro12-errmax
      area = area+area12-rlist(maxerr)
      if(defab1.eq.error1.or.defab2.eq.error2)go to 15
      if(dabs(rlist(maxerr)-area12).gt.0.1e-4_f*dabs(area12).or.erro12.lt.0.99_f*errmax) go to 10
      if(extrap) iroff2 = iroff2+1
      if(.not.extrap) iroff1 = iroff1+1
 10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
 15   rlist(maxerr) = area1
      rlist(last) = area2
      errbnd = dmax1(epsabs,epsrel*dabs(area))
      !
      !           test for roundoff error and eventually set error flag.
      !
      if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
      if(iroff2.ge.5) ierro = 3
      !
      !           set error flag in the case that the number of
      !           subintervals equals limit.
      !
      if(last.eq.limit) ier = 1
      !
      !           set error flag in the case of bad integrand behaviour
      !           at some points of the integration range.
      !
      if(dmax1(dabs(a1),dabs(b2)).le.(0.1e1_f+0.1e3_f*epmach)*(dabs(a2)+0.1e4_f*uflow)) ier = 4
      !
      !           append the newly-created intervals to the list.
      !
      if(error2.gt.error1) go to 20
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
      go to 30
 20   alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
      !
      !           call subroutine dqpsrt to maintain the descending ordering
      !           in the list of error estimates and select the subinterval
      !           with nrmax-th largest error estimate (to be bisected next).
      !
 30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
      if(errsum.le.errbnd) go to 115
      if(ier.ne.0) go to 100
      if(last.eq.2) go to 80
      if(noext) go to 90
      erlarg = erlarg-erlast
      if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
      if(extrap) go to 40
      !
      !           test whether the interval to be bisected next is the
      !           smallest interval.
      !
      if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
      extrap = .true.
      nrmax = 2
 40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
      !
      !           the smallest interval has the largest error.
      !           before bisecting decrease the sum of the errors over the
      !           larger intervals (erlarg) and perform extrapolation.
      !
      id = nrmax
      jupbnd = last
      if(last.gt.(2+limit/2)) jupbnd = limit+3-last
      do k = id,jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        nrmax = nrmax+1
      end do
      !
      !           perform extrapolation.
      !
 60   numrl2 = numrl2+1
      rlist2(numrl2) = area
      call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres, ier)
      if(ier.eq.9) return
      ktmin = ktmin+1
      if(ktmin.gt.5.and.abserr.lt.0.1e-2_f*errsum) ier = 5
      if(abseps.ge.abserr) go to 70
      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = dmax1(epsabs,epsrel*dabs(reseps))
      if(abserr.le.ertest) go to 100
      !
      !            prepare bisection of the smallest interval.
      !
 70   if(numrl2.eq.1) noext = .true.
      if(ier.eq.5) go to 100
      maxerr = iord(1)
      errmax = elist(maxerr)
      nrmax = 1
      extrap = .false.
      small = small*0.5_f
      erlarg = errsum
      go to 90
 80   small = 0.375_f
      erlarg = errsum
      ertest = errbnd
      rlist2(2) = area
 90 continue
    !
    !           set final result and error estimate.
    !           ------------------------------------
    !
100 if(abserr.eq.oflow) go to 115
    if((ier+ierro).eq.0) go to 110
    if(ierro.eq.3) abserr = abserr+correc
    if(ier.eq.0) ier = 3
    if(result.ne.0.0_f.and.area.ne.0.0_f)go to 105
    if(abserr.gt.errsum)go to 115
    if(area.eq.0.0_f) go to 130
    go to 110
105 if(abserr/dabs(result).gt.errsum/dabs(area))go to 115
    !
    !           test on divergence
    !
110 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.defabs*0.1e-1_f) go to 130
    if(0.1e-1_f.gt.(result/area).or.(result/area).gt.0.1e3_f.or.errsum.gt.dabs(area)) ier = 6
    go to 130
    !
    !           compute global integral sum.
    !
115 result = 0.0_f
    do k = 1,last
      result = result+rlist(k)
    end do
    abserr = errsum
130 neval = 30*last-15
    if(inf.eq.2) neval = 2*neval
    if(ier.gt.2) ier=ier-1
999 return
  end subroutine dqagie


  !!***begin prologue  dqk15
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)
  !!***category no.  h2a1a2
  !!***keywords  15-point gauss-kronrod rules
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div - k.u.leuven
  !!***purpose  to compute i = integral of f over (a,b), with error
  !!                           estimate
  !!                       j = integral of abs(f) over (a,b)
  !!***description
  !!
  !!           integration rules
  !!           standard fortran subroutine
  !!           double precision version
  !!
  !!           parameters
  !!            on entry
  !!              fx     - double precision
  !!                       function subprogram defining the integrand
  !!                       function f(x). the actual name for f needs to be
  !!                       declared e x t e r n a l in the calling program.
  !!
  !!              fx_vars- structure containing variables need for integration
  !!                     specific to fractal meanfield scattering code!
  !!
  !!              a      - double precision
  !!                       lower limit of integration
  !!
  !!              b      - double precision
  !!                       upper limit of integration
  !!
  !!            on return
  !!              result - double precision
  !!                       approximation to the integral i
  !!                       result is computed by applying the 15-point
  !!                       kronrod rule (resk) obtained by optimal addition
  !!                       of abscissae to the7-point gauss rule(resg).
  !!
  !!              abserr - double precision
  !!                       estimate of the modulus of the absolute error,
  !!                       which should not exceed abs(i-result)
  !!
  !!              resabs - double precision
  !!                       approximation to the integral j
  !!
  !!              resasc - double precision
  !!                       approximation to the integral of abs(f-i/(b-a))
  !!                       over (a,b)
  !!
  !!		  ier    - integer
  !!                       ier = 0 normal and reliable termination of the
  !!                               routine. it is assumed that the requested
  !!                               accuracy has been achieved.
  !!                       ier.gt.0 abnormal termination of the routine. the
  !!                                estimates for result and error are less
  !!                                reliable. it is assumed that the requested
  !!                                accuracy has not been achieved.
  !!
  !!***references  (none)
  !!***routines called  sd1mach
  !!***end prologue  dqk15  
  !!
  subroutine dqk15(fx,fx_vars,a,b,result,abserr,resabs,resasc,ier)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: a
    real(kind=f) :: b 
    real(kind=f) :: result
    real(kind=f) :: abserr
    real(kind=f) :: resabs
    real(kind=f) :: resasc
    integer :: ier

    ! Local Declarations
    real(kind=f) :: absc, centr, dabs, dhlgth, dmax2, dmin1
    real(kind=f) :: epmach, fc, fsum, fval1, fval2, fv1(7), fv2(7), hlgth
    real(kind=f) :: resg, resk, reskh, uflow, wg(4), wgk(8), xgk(8)
    integer :: j,jtw,jtwm1

    !
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 15-point kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 7-point
    !                    gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 7-point gauss rule
    !
    !           wgk    - weights of the 15-point kronrod rule
    !
    !           wg     - weights of the 7-point gauss rule
    !
    !
    ! gauss quadrature weights and kronron quadrature abscissae and weights
    ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    ! bell labs, nov. 1981.
    !
    data wg  (  1) / 0.129484966168869693270611432679082_f /
    data wg  (  2) / 0.279705391489276667901467771423780_f /
    data wg  (  3) / 0.381830050505118944950369775488975_f /
    data wg  (  4) / 0.417959183673469387755102040816327_f /

    data xgk (  1) / 0.991455371120812639206854697526329_f /
    data xgk (  2) / 0.949107912342758524526189684047851_f /
    data xgk (  3) / 0.864864423359769072789712788640926_f /
    data xgk (  4) / 0.741531185599394439863864773280788_f /
    data xgk (  5) / 0.586087235467691130294144838258730_f /
    data xgk (  6) / 0.405845151377397166906606412076961_f /
    data xgk (  7) / 0.207784955007898467600689403773245_f /
    data xgk (  8) / 0.000000000000000000000000000000000_f /

    data wgk (  1) / 0.022935322010529224963732008058970_f /
    data wgk (  2) / 0.063092092629978553290700663189204_f /
    data wgk (  3) / 0.104790010322250183839876322541518_f /
    data wgk (  4) / 0.140653259715525918745189590510238_f /
    data wgk (  5) / 0.169004726639267902826583426598550_f /
    data wgk (  6) / 0.190350578064785409913256402421014_f /
    data wgk (  7) / 0.204432940075298892414161999234649_f /
    data wgk (  8) / 0.209482141084727828012999174891714_f /

    !
    !
    !           list of major variables
    !           -----------------------
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 7-point gauss formula
    !           resk   - result of the 15-point kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach is the largest relative spacing.
    !           uflow is the smallest positive magnitude.
    !
    !***first executable statement  dqk15
    !epmach = d1mach(4)
    !uflow = d1mach(1)

    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    call sd1mach(1,uflow,ier)
    if(ier.eq.9) return

    centr = 0.5_f*(a+b)
    hlgth = 0.5_f*(b-a)
    dhlgth = dabs(hlgth)
    !
    !           compute the 15-point kronrod approximation to
    !           the integral, and estimate the absolute error.
    !
    fc = fx(centr,fx_vars)
    resg = fc*wg(4)
    resk = fc*wgk(8)
    resabs = dabs(resk)
    do j=1,3
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = fx(centr-absc,fx_vars)
      fval2 = fx(centr+absc,fx_vars)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1+fval2
      resg = resg+wg(j)*fsum
      resk = resk+wgk(jtw)*fsum
      resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,4
      jtwm1 = j*2-1
      absc = hlgth*xgk(jtwm1)
      fval1 = fx(centr-absc,fx_vars)
      fval2 = fx(centr+absc,fx_vars)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1+fval2
      resk = resk+wgk(jtwm1)*fsum
      resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5_f
    resasc = wgk(8)*dabs(fc-reskh)
    do j=1,7
      resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc.ne.0.0_f.and.abserr.ne.0.0_f) abserr = resasc*dmin1(0.1e1_f,(0.2e3_f*abserr/resasc)**1.5_f)
    if(resabs.gt.uflow/(0.5e2_f*epmach)) abserr = dmax1((epmach*0.5e2_f)*resabs,abserr)
999 return
  end subroutine dqk15

  !!
  !!***begin prologue  dqk21
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)
  !!***category no.  h2a1a2
  !!***keywords  21-point gauss-kronrod rules
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !!***purpose  to compute i = integral of f over (a,b), with error
  !!                           estimate
  !!                       j = integral of abs(f) over (a,b)
  !!***description
  !!
  !!           integration rules
  !!           standard fortran subroutine
  !!           double precision version
  !!
  !!           parameters
  !!            on entry
  !!              fx     - double precision
  !!                       function subprogram defining the integrand
  !!                       function f(x). the actual name for f needs to be
  !!                       declared e x t e r n a l in the driver program.
  !!
  !!            fx_vars- structure containing variables need for integration
  !!                     specific to fractal meanfield scattering code
  !!              a      - double precision
  !!                       lower limit of integration
  !!
  !!              b      - double precision
  !!                       upper limit of integration
  !!
  !!            on return
  !!              result - double precision
  !!                       approximation to the integral i
  !!                       result is computed by applying the 21-point
  !!                       kronrod rule (resk) obtained by optimal addition
  !!                       of abscissae to the 10-point gauss rule (resg).
  !!
  !!              abserr - double precision
  !!                       estimate of the modulus of the absolute error,
  !!                       which should not exceed abs(i-result)
  !!
  !!              resabs - double precision
  !!                       approximation to the integral j
  !!
  !!              resasc - double precision
  !!                       approximation to the integral of abs(f-i/(b-a))
  !!                       over (a,b)
  !!
  !!		  ier    - integer
  !!                       ier = 0 normal and reliable termination of the
  !!                               routine. it is assumed that the requested
  !!                               accuracy has been achieved.
  !!                       ier.gt.0 abnormal termination of the routine. the
  !!                                estimates for result and error are less
  !!                                reliable. it is assumed that the requested
  !!                                accuracy has not been achieved.
  !!
  !!
  !!***references  (none)
  !!***routines called  sd1mach
  !!***end prologue  dqk21
  !!
  subroutine dqk21(fx,fx_vars,a,b,result,abserr,resabs,resasc, ier)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: a
    real(kind=f) :: b
    real(kind=f) :: result
    real(kind=f) :: abserr
    real(kind=f) :: resabs
    real(kind=f) :: resasc
    integer :: ier

    ! Local declarations
    real(kind=f) :: absc, centr, dabs, dhlgth, dmax1, dmin1
    real(kind=f) :: epmach, fc, fsum, fval1, fval2, fv1(10), fv2(10), hlgth
    real(kind=f) :: resg, resk, reskh, uflow, wg(5), wgk(11),xgk(11)
    integer :: j,jtw,jtwm1
 
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 21-point kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 10-point
    !                    gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 10-point gauss rule
    !
    !           wgk    - weights of the 21-point kronrod rule
    !
    !           wg     - weights of the 10-point gauss rule
    !
    !
    ! gauss quadrature weights and kronron quadrature abscissae and weights
    ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    ! bell labs, nov. 1981.
    !
    data wg  (  1) / 0.066671344308688137593568809893332_f /
    data wg  (  2) / 0.149451349150580593145776339657697_f /
    data wg  (  3) / 0.219086362515982043995534934228163_f /
    data wg  (  4) / 0.269266719309996355091226921569469_f /
    data wg  (  5) / 0.295524224714752870173892994651338_f /

    data xgk (  1) / 0.995657163025808080735527280689003_f /
    data xgk (  2) / 0.973906528517171720077964012084452_f /
    data xgk (  3) / 0.930157491355708226001207180059508_f /
    data xgk (  4) / 0.865063366688984510732096688423493_f /
    data xgk (  5) / 0.780817726586416897063717578345042_f /
    data xgk (  6) / 0.679409568299024406234327365114874_f /
    data xgk (  7) / 0.562757134668604683339000099272694_f /
    data xgk (  8) / 0.433395394129247190799265943165784_f /
    data xgk (  9) / 0.294392862701460198131126603103866_f /
    data xgk ( 10) / 0.148874338981631210884826001129720_f /
    data xgk ( 11) / 0.000000000000000000000000000000000_f /

    data wgk (  1) / 0.011694638867371874278064396062192_f /
    data wgk (  2) / 0.032558162307964727478818972459390_f /
    data wgk (  3) / 0.054755896574351996031381300244580_f /
    data wgk (  4) / 0.075039674810919952767043140916190_f /
    data wgk (  5) / 0.093125454583697605535065465083366_f /
    data wgk (  6) / 0.109387158802297641899210590325805_f /
    data wgk (  7) / 0.123491976262065851077958109831074_f /
    data wgk (  8) / 0.134709217311473325928054001771707_f /
    data wgk (  9) / 0.142775938577060080797094273138717_f /
    data wgk ( 10) / 0.147739104901338491374841515972068_f /
    data wgk ( 11) / 0.149445554002916905664936468389821_f /

    !
    !           list of major variables
    !           -----------------------
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 10-point gauss formula
    !           resk   - result of the 21-point kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach is the largest relative spacing.
    !           uflow is the smallest positive magnitude.
    !
    !***first executable statement  dqk21
    !epmach = d1mach(4)
    !uflow = d1mach(1)

    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    call sd1mach(1,uflow,ier)
    if(ier.eq.9) return

    centr = 0.5_f*(a+b)
    hlgth = 0.5_f*(b-a)
    dhlgth = dabs(hlgth)
    !
    !           compute the 21-point kronrod approximation to
    !           the integral, and estimate the absolute error.
    !
    resg = 0.0_f
    fc = fx(centr, fx_vars)
    resk = wgk(11)*fc
    resabs = dabs(resk)
    do j=1,5
      jtw = 2*j
      absc = hlgth*xgk(jtw)
      fval1 = fx(centr-absc, fx_vars)
      fval2 = fx(centr+absc, fx_vars)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1+fval2
      resg = resg+wg(j)*fsum
      resk = resk+wgk(jtw)*fsum
      resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,5
      jtwm1 = 2*j-1
      absc = hlgth*xgk(jtwm1)
      fval1 = fx(centr-absc, fx_vars)
      fval2 = fx(centr+absc, fx_vars)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1+fval2
      resk = resk+wgk(jtwm1)*fsum
      resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5_f
    resasc = wgk(11)*dabs(fc-reskh)
    do j=1,10
      resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc.ne.0.0_f.and.abserr.ne.0.0_f) abserr = resasc*dmin1(0.1e1_f,(0.2e3_f*abserr/resasc)**1.5_f)
    if(resabs.gt.uflow/(0.5e2_f*epmach)) abserr = dmax1((epmach*0.5e2_f)*resabs,abserr)
999 return
  end subroutine dqk21

  !!***begin prologue  dqk31
  !!***date written   800101   (yymmdd)
  !!***revision date  130519   (yymmdd)
  !!***category no.  h2a1a2
  !!***keywords  31-point gauss-kronrod rules
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !!***purpose  to compute i = integral of f over (a,b) with error
  !!                           estimate
  !!                       j = integral of abs(f) over (a,b)
  !!***description
  !!
  !!           integration rules
  !!           standard fortran subroutine
  !!           double precision version
  !!
  !!           parameters
  !!            on entry
  !!              fx     - double precision
  !!                       function subprogram defining the integrand
  !!                       function f(x). the actual name for f needs to be
  !!                       declared e x t e r n a l in the calling program.
  !!
  !!            fx_vars- structure containing variables need for integration
  !!                     specific to fractal meanfield scattering code!
  !!
  !!              a      - double precision
  !!                       lower limit of integration
  !!
  !!              b      - double precision
  !!                       upper limit of integration
  !!
  !!            on return
  !!              result - double precision
  !!                       approximation to the integral i
  !!                       result is computed by applying the 31-point
  !!                       gauss-kronrod rule (resk), obtained by optimal
  !!                       addition of abscissae to the 15-point gauss
  !!                       rule (resg).
  !!
  !!              abserr - double precison
  !!                       estimate of the modulus of the modulus,
  !!                       which should not exceed abs(i-result)
  !!
  !!              resabs - double precision
  !!                       approximation to the integral j
  !!
  !!              resasc - double precision
  !!                       approximation to the integral of abs(f-i/(b-a))
  !!                       over (a,b)
  !!
  !!		  ier    - integer
  !!                       ier = 0 normal and reliable termination of the
  !!                               routine. it is assumed that the requested
  !!                               accuracy has been achieved.
  !!                       ier.gt.0 abnormal termination of the routine. the
  !!                                estimates for result and error are less
  !!                                reliable. it is assumed that the requested
  !!                                accuracy has not been achieved.
  !!
  !!
  !!***references  (none)
  !!***routines called  sd1mach
  !!***end prologue  dqk31
  subroutine dqk31(fx,fx_vars,a,b,result,abserr,resabs,resasc, ier)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: a
    real(kind=f) :: b
    real(kind=f) :: result
    real(kind=f) :: abserr
    real(kind=f) :: resabs
    real(kind=f) :: resasc
    integer :: ier

    ! Local declarations
    real(kind=f) :: absc, centr, dabs, dhlgth, dmax1, dmin1
    real(kind=f) :: epmach, fc, fsum, fval1, fval2, fv1(15), fv2(15), hlgth
    real(kind=f) :: resg, resk, reskh, uflow, wg(8), wgk(16), xgk(16)
    integer :: j,jtw,jtwm1

    !
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 31-point kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 15-point
    !                    gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 15-point gauss rule
    !
    !           wgk    - weights of the 31-point kronrod rule
    !
    !           wg     - weights of the 15-point gauss rule
    !
    !
    ! gauss quadrature weights and kronron quadrature abscissae and weights
    ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    ! bell labs, nov. 1981.
    !
    data wg  (  1) / 0.030753241996117268354628393577204_f /
    data wg  (  2) / 0.070366047488108124709267416450667_f /
    data wg  (  3) / 0.107159220467171935011869546685869_f /
    data wg  (  4) / 0.139570677926154314447804794511028_f /
    data wg  (  5) / 0.166269205816993933553200860481209_f /
    data wg  (  6) / 0.186161000015562211026800561866423_f /
    data wg  (  7) / 0.198431485327111576456118326443839_f /
    data wg  (  8) / 0.202578241925561272880620199967519_f /

    data xgk (  1) / 0.998002298693397060285172840152271_f /
    data xgk (  2) / 0.987992518020485428489565718586613_f /
    data xgk (  3) / 0.967739075679139134257347978784337_f /
    data xgk (  4) / 0.937273392400705904307758947710209_f /
    data xgk (  5) / 0.897264532344081900882509656454496_f /
    data xgk (  6) / 0.848206583410427216200648320774217_f /
    data xgk (  7) / 0.790418501442465932967649294817947_f /
    data xgk (  8) / 0.724417731360170047416186054613938_f /
    data xgk (  9) / 0.650996741297416970533735895313275_f /
    data xgk ( 10) / 0.570972172608538847537226737253911_f /
    data xgk ( 11) / 0.485081863640239680693655740232351_f /
    data xgk ( 12) / 0.394151347077563369897207370981045_f /
    data xgk ( 13) / 0.299180007153168812166780024266389_f /
    data xgk ( 14) / 0.201194093997434522300628303394596_f /
    data xgk ( 15) / 0.101142066918717499027074231447392_f /
    data xgk ( 16) / 0.000000000000000000000000000000000_f /

    data wgk (  1) / 0.005377479872923348987792051430128_f /
    data wgk (  2) / 0.015007947329316122538374763075807_f /
    data wgk (  3) / 0.025460847326715320186874001019653_f /
    data wgk (  4) / 0.035346360791375846222037948478360_f /
    data wgk (  5) / 0.044589751324764876608227299373280_f /
    data wgk (  6) / 0.053481524690928087265343147239430_f /
    data wgk (  7) / 0.062009567800670640285139230960803_f /
    data wgk (  8) / 0.069854121318728258709520077099147_f /
    data wgk (  9) / 0.076849680757720378894432777482659_f /
    data wgk ( 10) / 0.083080502823133021038289247286104_f /
    data wgk ( 11) / 0.088564443056211770647275443693774_f /
    data wgk ( 12) / 0.093126598170825321225486872747346_f /
    data wgk ( 13) / 0.096642726983623678505179907627589_f /
    data wgk ( 14) / 0.099173598721791959332393173484603_f /
    data wgk ( 15) / 0.100769845523875595044946662617570_f /
    data wgk ( 16) / 0.101330007014791549017374792767493_f /
    !
    !
    !           list of major variables
    !           -----------------------
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 15-point gauss formula
    !           resk   - result of the 31-point kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    !           machine dependent constants
    !           ---------------------------
    !           epmach is the largest relative spacing.
    !           uflow is the smallest positive magnitude.
    !***first executable statement  dqk31
    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    call sd1mach(1,uflow,ier)
    if(ier.eq.9) return

    !epmach = d1mach(4)
    !uflow = d1mach(1)

    centr = 0.5_f*(a+b)
    hlgth = 0.5_f*(b-a)
    dhlgth = dabs(hlgth)
    !
    !           compute the 31-point kronrod approximation to
    !           the integral, and estimate the absolute error.
    !
    fc = fx(centr, fx_vars)
    resg = wg(8)*fc
    resk = wgk(16)*fc
    resabs = dabs(resk)
    do j=1,7
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = fx(centr-absc, fx_vars)
      fval2 = fx(centr+absc, fx_vars)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1+fval2
      resg = resg+wg(j)*fsum
      resk = resk+wgk(jtw)*fsum
      resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,8
      jtwm1 = j*2-1
      absc = hlgth*xgk(jtwm1)
      fval1 = fx(centr-absc, fx_vars)
      fval2 = fx(centr+absc, fx_vars)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1+fval2
      resk = resk+wgk(jtwm1)*fsum
      resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5_f
    resasc = wgk(16)*dabs(fc-reskh)
    do j=1,15
      resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc.ne.0.0_f.and.abserr.ne.0.0_f) abserr = resasc*dmin1(0.1e1_f,(0.2e3_f*abserr/resasc)**1.5_f)
    if(resabs.gt.uflow/(0.5e2_f*epmach)) abserr = dmax1((epmach*0.5e2_f)*resabs,abserr)
999 return
  end subroutine dqk31



  !!***begin prologue  dqk41
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)
  !!***category no.  h2a1a2
  !!***keywords  41-point gauss-kronrod rules
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !!***purpose  to compute i = integral of f over (a,b), with error
  !!                           estimate
  !!                       j = integral of abs(f) over (a,b)
  !!***description
  !!
  !!           integration rules
  !!           standard fortran subroutine
  !!           double precision version
  !!
  !!           parameters
  !!            on entry
  !!              fx     - double precision
  !!                       function subprogram defining the integrand
  !!                       function f(x). the actual name for f needs to be
  !!                       declared e x t e r n a l in the calling program.
  !!
  !!            fx_vars- structure containing variables need for integration
  !!                     specific to fractal meanfield scattering code
  !!
  !!              a      - double precision
  !!                       lower limit of integration
  !!
  !!              b      - double precision
  !!                       upper limit of integration
  !!
  !!            on return
  !!              result - double precision
  !!                       approximation to the integral i
  !!                       result is computed by applying the 41-point
  !!                       gauss-kronrod rule (resk) obtained by optimal
  !!                       addition of abscissae to the 20-point gauss
  !!                       rule (resg).
  !!
  !!              abserr - double precision
  !!                       estimate of the modulus of the absolute error,
  !!                       which should not exceed abs(i-result)
  !!
  !!              resabs - double precision
  !!                       approximation to the integral j
  !!
  !!              resasc - double precision
  !!                       approximation to the integal of abs(f-i/(b-a))
  !!                       over (a,b)
  !!
  !!		  ier    - integer
  !!                       ier = 0 normal and reliable termination of the
  !!                               routine. it is assumed that the requested
  !!                               accuracy has been achieved.
  !!                       ier.gt.0 abnormal termination of the routine. the
  !!                                estimates for result and error are less
  !!                                reliable. it is assumed that the requested
  !!                                accuracy has not been achieved.
  !!
  !!
  !!***references  (none)
  !!***routines called  sd1mach
  !!***end prologue  dqk41
  !!
  subroutine dqk41(fx,fx_vars,a,b,result,abserr,resabs,resasc, ier)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: a
    real(kind=f) :: b
    real(kind=f) :: result
    real(kind=f) :: abserr
    real(kind=f) :: resabs
    real(kind=f) :: resasc
    integer :: ier

    ! Local declarations
    real(kind=f) :: absc, centr, dabs, dhlgth, dmax1, dmin1
    real(kind=f) :: epmach, fc, fsum, fval1, fval2, fv1(20), fv2(20), hlgth
    real(kind=f) :: resg, resk, reskh, uflow, wg(10), wgk(21), xgk(21)
    integer :: j, jtw, jtwm1

    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 41-point gauss-kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 20-point
    !                    gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 20-point gauss rule
    !
    !           wgk    - weights of the 41-point gauss-kronrod rule
    !
    !           wg     - weights of the 20-point gauss rule
    !
    !
    ! gauss quadrature weights and kronron quadrature abscissae and weights
    ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    ! bell labs, nov. 1981.
    !
    data wg  (  1) / 0.017614007139152118311861962351853_f /
    data wg  (  2) / 0.040601429800386941331039952274932_f /
    data wg  (  3) / 0.062672048334109063569506535187042_f /
    data wg  (  4) / 0.083276741576704748724758143222046_f /
    data wg  (  5) / 0.101930119817240435036750135480350_f /
    data wg  (  6) / 0.118194531961518417312377377711382_f /
    data wg  (  7) / 0.131688638449176626898494499748163_f /
    data wg  (  8) / 0.142096109318382051329298325067165_f /
    data wg  (  9) / 0.149172986472603746787828737001969_f /
    data wg  ( 10) / 0.152753387130725850698084331955098_f /

    data xgk (  1) / 0.998859031588277663838315576545863_f /
    data xgk (  2) / 0.993128599185094924786122388471320_f /
    data xgk (  3) / 0.981507877450250259193342994720217_f /
    data xgk (  4) / 0.963971927277913791267666131197277_f /
    data xgk (  5) / 0.940822633831754753519982722212443_f /
    data xgk (  6) / 0.912234428251325905867752441203298_f /
    data xgk (  7) / 0.878276811252281976077442995113078_f /
    data xgk (  8) / 0.839116971822218823394529061701521_f /
    data xgk (  9) / 0.795041428837551198350638833272788_f /
    data xgk ( 10) / 0.746331906460150792614305070355642_f /
    data xgk ( 11) / 0.693237656334751384805490711845932_f /
    data xgk ( 12) / 0.636053680726515025452836696226286_f /
    data xgk ( 13) / 0.575140446819710315342946036586425_f /
    data xgk ( 14) / 0.510867001950827098004364050955251_f/
    data xgk ( 15) / 0.443593175238725103199992213492640_f /
    data xgk ( 16) / 0.373706088715419560672548177024927_f /
    data xgk ( 17) / 0.301627868114913004320555356858592_f /
    data xgk ( 18) / 0.227785851141645078080496195368575_f /
    data xgk ( 19) / 0.152605465240922675505220241022678_f /
    data xgk ( 20) / 0.076526521133497333754640409398838_f /
    data xgk ( 21) / 0.000000000000000000000000000000000_f /

    data wgk (  1) / 0.003073583718520531501218293246031_f /
    data wgk (  2) / 0.008600269855642942198661787950102_f /
    data wgk (  3) / 0.014626169256971252983787960308868_f /
    data wgk (  4) / 0.020388373461266523598010231432755_f /
    data wgk (  5) / 0.025882133604951158834505067096153_f /
    data wgk (  6) / 0.031287306777032798958543119323801_f /
    data wgk (  7) / 0.036600169758200798030557240707211_f /
    data wgk (  8) / 0.041668873327973686263788305936895_f /
    data wgk (  9) / 0.046434821867497674720231880926108_f /
    data wgk ( 10) / 0.050944573923728691932707670050345_f /
    data wgk ( 11) / 0.055195105348285994744832372419777_f /
    data wgk ( 12) / 0.059111400880639572374967220648594_f /
    data wgk ( 13) / 0.062653237554781168025870122174255_f /
    data wgk ( 14) / 0.065834597133618422111563556969398_f /
    data wgk ( 15) / 0.068648672928521619345623411885368_f /
    data wgk ( 16) / 0.071054423553444068305790361723210_f /
    data wgk ( 17) / 0.073030690332786667495189417658913_f /
    data wgk ( 18) / 0.074582875400499188986581418362488_f /
    data wgk ( 19) / 0.075704497684556674659542775376617_f /
    data wgk ( 20) / 0.076377867672080736705502835038061_f /
    data wgk ( 21) / 0.076600711917999656445049901530102_f /
   
    !
    !           list of major variables
    !           -----------------------
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 20-point gauss formula
    !           resk   - result of the 41-point kronrod formula
    !           reskh  - approximation to mean value of f over (a,b), i.e.
    !                    to i/(b-a)
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach is the largest relative spacing.
    !           uflow is the smallest positive magnitude.
    !
    !***first executable statement  dqk41
    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    call sd1mach(1,uflow,ier)
    if(ier.eq.9) return

    !epmach = d1mach(4)
    !uflow = d1mach(1)

    centr = 0.5_f*(a+b)
    hlgth = 0.5_f*(b-a)
    dhlgth = dabs(hlgth)
    !
    !           compute the 41-point gauss-kronrod approximation to
    !           the integral, and estimate the absolute error.
    !
    resg = 0.0_f
    fc = fx(centr, fx_vars)
    resk = wgk(21)*fc
    resabs = dabs(resk)
    do j=1,10
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = fx(centr-absc, fx_vars)
      fval2 = fx(centr+absc, fx_vars)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1+fval2
      resg = resg+wg(j)*fsum
      resk = resk+wgk(jtw)*fsum
      resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,10
      jtwm1 = j*2-1
      absc = hlgth*xgk(jtwm1)
      fval1 = fx(centr-absc, fx_vars)
      fval2 = fx(centr+absc, fx_vars)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1+fval2
      resk = resk+wgk(jtwm1)*fsum
      resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5_f
    resasc = wgk(21)*dabs(fc-reskh)
    do j=1,20
      resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc.ne.0.0_f.and.abserr.ne.0._f) abserr = resasc*dmin1(0.1e1_f,(0.2e3_f*abserr/resasc)**1.5_f)
    if(resabs.gt.uflow/(0.5e2_f*epmach)) abserr = dmax1((epmach*0.5e2_f)*resabs,abserr)
999 return
  end subroutine dqk41


  !!***begin prologue  dqk51
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)
  !!***category no.  h2a1a2
  !!***keywords  51-point gauss-kronrod rules
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math & progr. div. - k.u.leuven
  !!***purpose  to compute i = integral of f over (a,b) with error
  !!                           estimate
  !!                       j = integral of abs(f) over (a,b)
  !!***description
  !!
  !!           integration rules
  !!           standard fortran subroutine
  !!           double precision version
  !!
  !!           parameters
  !!            on entry
  !!              fx     - double precision
  !!                       function subroutine defining the integrand
  !!                       function f(x). the actual name for f needs to be
  !!                       declared e x t e r n a l in the calling program.
  !!
  !!              fx_vars- structure containing variables need for integration
  !!                       specific to fractal meanfield scattering code
  !!
  !!              a      - double precision
  !!                       lower limit of integration
  !!
  !!              b      - double precision
  !!                       upper limit of integration
  !!
  !!            on return
  !!              result - double precision
  !!                       approximation to the integral i
  !!                       result is computed by applying the 51-point
  !!                       kronrod rule (resk) obtained by optimal addition
  !!                       of abscissae to the 25-point gauss rule (resg).
  !!
  !!              abserr - double precision
  !!                       estimate of the modulus of the absolute error,
  !!                       which should not exceed abs(i-result)
  !!
  !!              resabs - double precision
  !!                       approximation to the integral j
  !!
  !!              resasc - double precision
  !!                       approximation to the integral of abs(f-i/(b-a))
  !!                       over (a,b)
  !!
  !!		  ier    - integer
  !!                       ier = 0 normal and reliable termination of the
  !!                               routine. it is assumed that the requested
  !!                               accuracy has been achieved.
  !!                       ier.gt.0 abnormal termination of the routine. the
  !!                                estimates for result and error are less
  !!                                reliable. it is assumed that the requested
  !!                                accuracy has not been achieved.
  !!
  !!***references  (none)
  !!***routines called  sd1mach
  !!***end prologue  dqk51
  !!
  subroutine dqk51(fx,fx_vars,a,b,result,abserr,resabs,resasc,ier)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: a
    real(kind=f) :: b
    real(kind=f) :: result
    real(kind=f) :: abserr
    real(kind=f) :: resabs
    real(kind=f) :: resasc
    integer :: ier

    ! Local declarations
    real(kind=f) ::  absc, centr, dabs, dhlgth, dmax1, dmin1
    real(kind=f) ::  epmach, fc, fsum, fval1, fval2, fv1(25), fv2(25), hlgth
    real(kind=f) ::  resg, resk, reskh, uflow, wg(13),wgk(26), xgk(26)
    integer j,jtw,jtwm1

    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 51-point kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 25-point
    !                    gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 25-point gauss rule
    !
    !           wgk    - weights of the 51-point kronrod rule
    !
    !           wg     - weights of the 25-point gauss rule
    !
    !
    ! gauss quadrature weights and kronron quadrature abscissae and weights
    ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    ! bell labs, nov. 1981.
    !
    data wg  (  1) / 0.011393798501026287947902964113235_f /
    data wg  (  2) / 0.026354986615032137261901815295299_f /
    data wg  (  3) / 0.040939156701306312655623487711646_f /
    data wg  (  4) / 0.054904695975835191925936891540473_f /
    data wg  (  5) / 0.068038333812356917207187185656708_f /
    data wg  (  6) / 0.080140700335001018013234959669111_f /
    data wg  (  7) / 0.091028261982963649811497220702892_f /
    data wg  (  8) / 0.100535949067050644202206890392686_f /
    data wg  (  9) / 0.108519624474263653116093957050117_f /
    data wg  ( 10) / 0.114858259145711648339325545869556_f /
    data wg  ( 11) / 0.119455763535784772228178126512901_f /
    data wg  ( 12) / 0.122242442990310041688959518945852_f /
    data wg  ( 13) / 0.123176053726715451203902873079050_f /

    data xgk (  1) / 0.999262104992609834193457486540341_f /
    data xgk (  2) / 0.995556969790498097908784946893902_f /
    data xgk (  3) / 0.988035794534077247637331014577406_f /
    data xgk (  4) / 0.976663921459517511498315386479594_f /
    data xgk (  5) / 0.961614986425842512418130033660167_f /
    data xgk (  6) / 0.942974571228974339414011169658471_f /
    data xgk (  7) / 0.920747115281701561746346084546331_f /
    data xgk (  8) / 0.894991997878275368851042006782805_f /
    data xgk (  9) / 0.865847065293275595448996969588340_f /
    data xgk ( 10) / 0.833442628760834001421021108693570_f /
    data xgk ( 11) / 0.797873797998500059410410904994307_f /
    data xgk ( 12) / 0.759259263037357630577282865204361_f /
    data xgk ( 13) / 0.717766406813084388186654079773298_f /
    data xgk ( 14) / 0.673566368473468364485120633247622_f /
    data xgk ( 15) / 0.626810099010317412788122681624518_f /
    data xgk ( 16) / 0.577662930241222967723689841612654_f /
    data xgk ( 17) / 0.526325284334719182599623778158010_f /
    data xgk ( 18) / 0.473002731445714960522182115009192_f /
    data xgk ( 19) / 0.417885382193037748851814394594572_f /
    data xgk ( 20) / 0.361172305809387837735821730127641_f /
    data xgk ( 21) / 0.303089538931107830167478909980339_f /
    data xgk ( 22) / 0.243866883720988432045190362797452_f /
    data xgk ( 23) / 0.183718939421048892015969888759528_f /
    data xgk ( 24) / 0.122864692610710396387359818808037_f /
    data xgk ( 25) / 0.061544483005685078886546392366797_f /
    data xgk ( 26) / 0.000000000000000000000000000000000_f /

    data wgk (  1) / 0.001987383892330315926507851882843_f /
    data wgk (  2) / 0.005561932135356713758040236901066_f /
    data wgk (  3) / 0.009473973386174151607207710523655_f /
    data wgk (  4) / 0.013236229195571674813656405846976_f /
    data wgk (  5) / 0.016847817709128298231516667536336_f /
    data wgk (  6) / 0.020435371145882835456568292235939_f /
    data wgk (  7) / 0.024009945606953216220092489164881_f /
    data wgk (  8) / 0.027475317587851737802948455517811_f /
    data wgk (  9) / 0.030792300167387488891109020215229_f /
    data wgk ( 10) / 0.034002130274329337836748795229551_f /
    data wgk ( 11) / 0.037116271483415543560330625367620_f /
    data wgk ( 12) / 0.040083825504032382074839284467076_f /
    data wgk ( 13) / 0.042872845020170049476895792439495_f /
    data wgk ( 14) / 0.045502913049921788909870584752660_f /
    data wgk ( 15) / 0.047982537138836713906392255756915_f /
    data wgk ( 16) / 0.050277679080715671963325259433440_f /
    data wgk ( 17) / 0.052362885806407475864366712137873_f /
    data wgk ( 18) / 0.054251129888545490144543370459876_f /
    data wgk ( 19) / 0.055950811220412317308240686382747_f /
    data wgk ( 20) / 0.057437116361567832853582693939506_f /
    data wgk ( 21) / 0.058689680022394207961974175856788_f /
    data wgk ( 22) / 0.059720340324174059979099291932562_f /
    data wgk ( 23) / 0.060539455376045862945360267517565_f /
    data wgk ( 24) / 0.061128509717053048305859030416293_f /
    data wgk ( 25) / 0.061471189871425316661544131965264_f /
    !       note: wgk (26) was calculated from the values of wgk(1..25)
    data wgk ( 26) / 0.061580818067832935078759824240066_f /
   
    !
    !           list of major variables
    !           -----------------------
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 25-point gauss formula
    !           resk   - result of the 51-point kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach is the largest relative spacing.
    !           uflow is the smallest positive magnitude.
    !
    !***first executable statement  dqk51
    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    call sd1mach(1,uflow,ier)
    if(ier.eq.9) return

    !epmach = d1mach(4)
    !uflow = d1mach(1)

    centr = 0.5_f*(a+b)
    hlgth = 0.5_f*(b-a)
    dhlgth = dabs(hlgth)
    !
    !           compute the 51-point kronrod approximation to
    !           the integral, and estimate the absolute error.
    !
    fc = fx(centr, fx_vars)
    resg = wg(13)*fc
    resk = wgk(26)*fc
    resabs = dabs(resk)
    do j=1,12
      jtw = j*2
      absc = hlgth*xgk(jtw)
      fval1 = fx(centr-absc, fx_vars)
      fval2 = fx(centr+absc, fx_vars)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1+fval2
      resg = resg+wg(j)*fsum
      resk = resk+wgk(jtw)*fsum
      resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,13
      jtwm1 = j*2-1
      absc = hlgth*xgk(jtwm1)
      fval1 = fx(centr-absc, fx_vars)
      fval2 = fx(centr+absc, fx_vars)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1+fval2
      resk = resk+wgk(jtwm1)*fsum
      resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5_f
    resasc = wgk(26)*dabs(fc-reskh)
    do j=1,25
      resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc.ne.0.0_f.and.abserr.ne.0.0_f) abserr = resasc*dmin1(0.1e1_f,(0.2e3_f*abserr/resasc)**1.5_f)
    if(resabs.gt.uflow/(0.5e2_f*epmach)) abserr = dmax1((epmach*0.5e2_f)*resabs,abserr)
999 return
  end subroutine dqk51

  !!***begin prologue  dqk61
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)
  !!***category no.  h2a1a2
  !!***keywords  61-point gauss-kronrod rules
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !!***purpose  to compute i = integral of f over (a,b) with error
  !!                           estimate
  !!                       j = integral of dabs(f) over (a,b)
  !!***description
  !!
  !!        integration rule
  !!        standard fortran subroutine
  !!        double precision version
  !!
  !!
  !!        parameters
  !!         on entry
  !!           fx     - double precision
  !!                    function subprogram defining the integrand
  !!                    function f(x). the actual name for f needs to be
  !!                    declared e x t e r n a l in the calling program.
  !!
  !!            fx_vars- structure containing variables need for integration
  !!                     specific to fractal meanfield scattering code
  !!
  !!           a      - double precision
  !!                    lower limit of integration
  !!
  !!           b      - double precision
  !!                    upper limit of integration
  !!
  !!         on return
  !!           result - double precision
  !!                    approximation to the integral i
  !!                    result is computed by applying the 61-point
  !!                    kronrod rule (resk) obtained by optimal addition of
  !!                    abscissae to the 30-point gauss rule (resg).
  !!
  !!           abserr - double precision
  !!                    estimate of the modulus of the absolute error,
  !!                    which should equal or exceed dabs(i-result)
  !!
  !!           resabs - double precision
  !!                    approximation to the integral j
  !!
  !!           resasc - double precision
  !!                    approximation to the integral of dabs(f-i/(b-a))
  !!
  !!		  ier    - integer
  !!                       ier = 0 normal and reliable termination of the
  !!                               routine. it is assumed that the requested
  !!                               accuracy has been achieved.
  !!                       ier.gt.0 abnormal termination of the routine. the
  !!                                estimates for result and error are less
  !!                                reliable. it is assumed that the requested
  !!                                accuracy has not been achieved.
  !!
  !!
  !!***references  (none)
  !!***routines called  sd1mach
  !!***end prologue  dqk61
  !!
  subroutine dqk61(fx,fx_vars,a,b,result,abserr,resabs,resasc, ier)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: a
    real(kind=f) :: b
    real(kind=f) :: result
    real(kind=f) :: abserr
    real(kind=f) :: resabs
    real(kind=f) :: resasc
    integer :: ier

    ! Local declartions
    real(kind=f) :: dabsc, centr, dabs, dhlgth, dmax1, dmin1
    real(kind=f) :: epmach, fc, fsum, fval1, fval2, fv1(30), fv2(30), hlgth
    real(kind=f) :: resg, resk, reskh, uflow, wg(15) ,wgk(31), xgk(31)
    integer :: j, jtw, jtwm1

    !
    !           the abscissae and weights are given for the
    !           interval (-1,1). because of symmetry only the positive
    !           abscissae and their corresponding weights are given.
    !
    !           xgk   - abscissae of the 61-point kronrod rule
    !                   xgk(2), xgk(4)  ... abscissae of the 30-point
    !                   gauss rule
    !                   xgk(1), xgk(3)  ... optimally added abscissae
    !                   to the 30-point gauss rule
    !
    !           wgk   - weights of the 61-point kronrod rule
    !
    !           wg    - weigths of the 30-point gauss rule
    !
    !
    ! gauss quadrature weights and kronron quadrature abscissae and weights
    ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    ! bell labs, nov. 1981.
    !
    data wg  (  1) / 0.007968192496166605615465883474674_f /
    data wg  (  2) / 0.018466468311090959142302131912047_f /
    data wg  (  3) / 0.028784707883323369349719179611292_f /
    data wg  (  4) / 0.038799192569627049596801936446348_f /
    data wg  (  5) / 0.048402672830594052902938140422808_f /
    data wg  (  6) / 0.057493156217619066481721689402056_f /
    data wg  (  7) / 0.065974229882180495128128515115962_f /
    data wg  (  8) / 0.073755974737705206268243850022191_f /
    data wg  (  9) / 0.080755895229420215354694938460530_f /
    data wg  ( 10) / 0.086899787201082979802387530715126_f /
    data wg  ( 11) / 0.092122522237786128717632707087619_f /
    data wg  ( 12) / 0.096368737174644259639468626351810_f /
    data wg  ( 13) / 0.099593420586795267062780282103569_f /
    data wg  ( 14) / 0.101762389748405504596428952168554_f /
    data wg  ( 15) / 0.102852652893558840341285636705415_f /

    data xgk (  1) / 0.999484410050490637571325895705811_f /
    data xgk (  2) / 0.996893484074649540271630050918695_f /
    data xgk (  3) / 0.991630996870404594858628366109486_f /
    data xgk (  4) / 0.983668123279747209970032581605663_f /
    data xgk (  5) / 0.973116322501126268374693868423707_f /
    data xgk (  6) / 0.960021864968307512216871025581798_f /
    data xgk (  7) / 0.944374444748559979415831324037439_f /
    data xgk (  8) / 0.926200047429274325879324277080474_f /
    data xgk (  9) / 0.905573307699907798546522558925958_f /
    data xgk ( 10) / 0.882560535792052681543116462530226_f /
    data xgk ( 11) / 0.857205233546061098958658510658944_f /
    data xgk ( 12) / 0.829565762382768397442898119732502_f /
    data xgk ( 13) / 0.799727835821839083013668942322683_f /
    data xgk ( 14) / 0.767777432104826194917977340974503_f /
    data xgk ( 15) / 0.733790062453226804726171131369528_f /
    data xgk ( 16) / 0.697850494793315796932292388026640_f /
    data xgk ( 17) / 0.660061064126626961370053668149271_f /
    data xgk ( 18) / 0.620526182989242861140477556431189_f /
    data xgk ( 19) / 0.579345235826361691756024932172540_f /
    data xgk ( 20) / 0.536624148142019899264169793311073_f /
    data xgk ( 21) / 0.492480467861778574993693061207709_f /
    data xgk ( 22) / 0.447033769538089176780609900322854_f /
    data xgk ( 23) / 0.400401254830394392535476211542661_f /
    data xgk ( 24) / 0.352704725530878113471037207089374_f /
    data xgk ( 25) / 0.304073202273625077372677107199257_f /
    data xgk ( 26) / 0.254636926167889846439805129817805_f /
    data xgk ( 27) / 0.204525116682309891438957671002025_f /
    data xgk ( 28) / 0.153869913608583546963794672743256_f /
    data xgk ( 29) / 0.102806937966737030147096751318001_f /
    data xgk ( 30) / 0.051471842555317695833025213166723_f /
    data xgk ( 31) / 0.000000000000000000000000000000000_f /

    data wgk (  1) / 0.001389013698677007624551591226760_f /
    data wgk (  2) / 0.003890461127099884051267201844516_f /
    data wgk (  3) / 0.006630703915931292173319826369750_f /
    data wgk (  4) / 0.009273279659517763428441146892024_f /
    data wgk (  5) / 0.011823015253496341742232898853251_f /
    data wgk (  6) / 0.014369729507045804812451432443580_f /
    data wgk (  7) / 0.016920889189053272627572289420322_f /
    data wgk (  8) / 0.019414141193942381173408951050128_f /
    data wgk (  9) / 0.021828035821609192297167485738339_f /
    data wgk ( 10) / 0.024191162078080601365686370725232_f /
    data wgk ( 11) / 0.026509954882333101610601709335075_f /
    data wgk ( 12) / 0.028754048765041292843978785354334_f /
    data wgk ( 13) / 0.030907257562387762472884252943092_f /
    data wgk ( 14) / 0.032981447057483726031814191016854_f /
    data wgk ( 15) / 0.034979338028060024137499670731468_f /
    data wgk ( 16) / 0.036882364651821229223911065617136_f /
    data wgk ( 17) / 0.038678945624727592950348651532281_f /
    data wgk ( 18) / 0.040374538951535959111995279752468_f /
    data wgk ( 19) / 0.041969810215164246147147541285970_f /
    data wgk ( 20) / 0.043452539701356069316831728117073_f /
    data wgk ( 21) / 0.044814800133162663192355551616723_f /
    data wgk ( 22) / 0.046059238271006988116271735559374_f /
    data wgk ( 23) / 0.047185546569299153945261478181099_f /
    data wgk ( 24) / 0.048185861757087129140779492298305_f /
    data wgk ( 25) / 0.049055434555029778887528165367238_f /
    data wgk ( 26) / 0.049795683427074206357811569379942_f /
    data wgk ( 27) / 0.050405921402782346840893085653585_f /
    data wgk ( 28) / 0.050881795898749606492297473049805_f /
    data wgk ( 29) / 0.051221547849258772170656282604944_f /
    data wgk ( 30) / 0.051426128537459025933862879215781_f /
    data wgk ( 31) / 0.051494729429451567558340433647099_f /
    
    !           list of major variables
    !           -----------------------
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           dabsc  - abscissa
    !           fval*  - function value
    !           resg   - result of the 30-point gauss rule
    !           resk   - result of the 61-point kronrod rule
    !           reskh  - approximation to the mean value of f
    !                    over (a,b), i.e. to i/(b-a)
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach is the largest relative spacing.
    ! a        uflow is the smallest positive magnitude.
    !
    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    call sd1mach(1,uflow,ier)
    if(ier.eq.9) return

    !epmach = d1mach(4)
    !uflow = d1mach(1)

    centr = 0.5_f*(b+a)
    hlgth = 0.5_f*(b-a)
    dhlgth = dabs(hlgth)
    !
    !           compute the 61-point kronrod approximation to the
    !           integral, and estimate the absolute error.
    !
    !***first executable statement  dqk61
    resg = 0.0_f
    fc = fx(centr, fx_vars)
    resk = wgk(31)*fc
    resabs = dabs(resk)
    do j=1,15
      jtw = j*2
      dabsc = hlgth*xgk(jtw)
      fval1 = fx(centr-dabsc, fx_vars)
      fval2 = fx(centr+dabsc, fx_vars)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1+fval2
      resg = resg+wg(j)*fsum
      resk = resk+wgk(jtw)*fsum
      resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j=1,15
      jtwm1 = j*2-1
      dabsc = hlgth*xgk(jtwm1)
      fval1 = fx(centr-dabsc, fx_vars)
      fval2 = fx(centr+dabsc, fx_vars)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1+fval2
      resk = resk+wgk(jtwm1)*fsum
      resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5_f
    resasc = wgk(31)*dabs(fc-reskh)
    do j=1,30
      resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc.ne.0.0_f.and.abserr.ne.0.0_f) abserr = resasc*dmin1(0.1e1_f,(0.2e3_f*abserr/resasc)**1.5_f)
    if(resabs.gt.uflow/(0.5e2_f*epmach)) abserr = dmax1((epmach*0.5e2_f)*resabs,abserr)
999 return
  end subroutine dqk61

  !!
  !!***begin prologue  dqpsrt
  !!***refer to  dqage,dqagie,dqagpe,dqawse
  !!***routines called  (none)
  !!***revision date  130319   (yymmdd)
  !!***keywords  sequential sorting
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !!***purpose  this routine maintains the descending ordering in the
  !!            list of the local error estimated resulting from the
  !!            interval subdivision process. at each call two error
  !!            estimates are inserted using the sequential search
  !!            method, top-down for the largest error estimate and
  !!            bottom-up for the smallest error estimate.
  !!***description
  !!
  !!           ordering routine
  !!           standard fortran subroutine
  !!           double precision version
  !!
  !!           parameters (meaning at output)
  !!              limit  - integer
  !!                       maximum number of error estimates the list
  !!                       can contain
  !!
  !!              last   - integer
  !!                       number of error estimates currently in the list
  !!
  !!              maxerr - integer
  !!                       maxerr points to the nrmax-th largest error
  !!                       estimate currently in the list
  !!
  !!              ermax  - double precision
  !!                       nrmax-th largest error estimate
  !!                       ermax = elist(maxerr)
  !!
  !!              elist  - double precision
  !!                       vector of dimension last containing
  !!                       the error estimates
  !!
  !!              iord   - integer
  !!                       vector of dimension last, the first k elements
  !!                       of which contain pointers to the error
  !!                       estimates, such that
  !!                       elist(iord(1)),...,  elist(iord(k))
  !!                       form a decreasing sequence, with
  !!                       k = last if last.le.(limit/2+2), and
  !!                       k = limit+1-last otherwise
  !!
  !!              nrmax  - integer
  !!                       maxerr = iord(nrmax)
  !!
  !!***end prologue  dqpsrt
  !!
  subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)

    ! Arguments
    integer :: limit
    integer :: last
    integer :: maxerr
    real(kind=f) :: ermax
    real(kind=f) :: elist(last)
    integer :: iord(last)
    integer :: nrmax

    ! Local declarations
    real(kind=f) :: errmax, errmin
    integer :: i, ibeg, ido, isucc, j, jbnd, jupbn, k

    !
    !           check whether the list contains more than
    !           two error estimates.
    !
    !***first executable statement  dqpsrt
    if(last.gt.2) go to 10
    iord(1) = 1
    iord(2) = 2
    go to 90
    !
    !           this part of the routine is only executed if, due to a
    !           difficult integrand, subdivision increased the error
    !           estimate. in the normal case the insert procedure should
    !           start after the nrmax-th largest error estimate.
    !
 10 errmax = elist(maxerr)
    if(nrmax.eq.1) go to 30
    ido = nrmax-1
    do i = 1,ido
      isucc = iord(nrmax-1)
      ! ***jump out of do-loop
      if(errmax.le.elist(isucc)) go to 30
      iord(nrmax) = isucc
      nrmax = nrmax-1
    end do
    !
    !           compute the number of elements in the list to be maintained
    !           in descending order. this number depends on the number of
    !           subdivisions still allowed.
    !
 30 jupbn = last
    if(last.gt.(limit/2+2)) jupbn = limit+3-last
    errmin = elist(last)
    !
    !           insert errmax by traversing the list top-down,
    !           starting comparison from the element elist(iord(nrmax+1)).
    !
    jbnd = jupbn-1
    ibeg = nrmax+1
    if(ibeg.gt.jbnd) go to 50
    do i=ibeg,jbnd
      isucc = iord(i)
      ! ***jump out of do-loop
      if(errmax.ge.elist(isucc)) go to 60
      iord(i-1) = isucc
    end do
 50 iord(jbnd) = maxerr
    iord(jupbn) = last
    go to 90
    !
    !           insert errmin by traversing the list bottom-up.
    !
 60 iord(i-1) = maxerr
    k = jbnd
    do j=i,jbnd
      isucc = iord(k)
      ! ***jump out of do-loop
      if(errmin.lt.elist(isucc)) go to 80
      iord(k+1) = isucc
      k = k-1
    end do
    iord(i) = last
    go to 90
 80 iord(k+1) = last
    !
    !           set maxerr and ermax.
    !
 90 maxerr = iord(nrmax)
    ermax = elist(maxerr)
    return
  end subroutine dqpsrt

  !!
  !!***begin prologue  dqk15i
  !!***date written   800101   (yymmdd)
  !!***revision date  130319   (yymmdd)
  !!***category no.  h2a3a2,h2a4a2
  !!***keywords  15-point transformed gauss-kronrod rules
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !!***purpose  the original (infinite integration range is mapped
  !!            onto the interval (0,1) and (a,b) is a part of (0,1).
  !!            it is the purpose to compute
  !!            i = integral of transformed integrand over (a,b),
  !!            j = integral of abs(transformed integrand) over (a,b).
  !!***description
  !!
  !!           integration rule
  !!           standard fortran subroutine
  !!           double precision version
  !!
  !!           parameters
  !!            on entry
  !!              fx     - double precision
  !!                       fuction subprogram defining the integrand
  !!                       function f(x). the actual name for f needs to be
  !!                       declared e x t e r n a l in the calling program.
  !!
  !!              fx_vars- structure containing variables need for integration
  !!                       specific to fractal meanfield scattering code!
  !!
  !!              boun   - double precision
  !!                       finite bound of original integration
  !!                       range (set to zero if inf = +2)
  !!
  !!              inf    - integer
  !!                       if inf = -1, the original interval is
  !!                                   (-infinity,bound),
  !!                       if inf = +1, the original interval is
  !!                                   (bound,+infinity),
  !!                       if inf = +2, the original interval is
  !!                                   (-infinity,+infinity) and
  !!                       the integral is computed as the sum of two
  !!                       integrals, one over (-infinity,0) and one over
  !!                       (0,+infinity).
  !!
  !!              a      - double precision
  !!                       lower limit for integration over subrange
  !!                       of (0,1)
  !!
  !!              b      - double precision
  !!                       upper limit for integration over subrange
  !!                       of (0,1)
  !!
  !!            on return
  !!              result - double precision
  !!                       approximation to the integral i
  !!                       result is computed by applying the 15-point
  !!                       kronrod rule(resk) obtained by optimal addition
  !!                       of abscissae to the 7-point gauss rule(resg).
  !!
  !!              abserr - double precision
  !!                       estimate of the modulus of the absolute error,
  !!                       which should equal or exceed abs(i-result)
  !!
  !!              resabs - double precision
  !!                       approximation to the integral j
  !!
  !!              resasc - double precision
  !!                       approximation to the integral of
  !!                       abs((transformed integrand)-i/(b-a)) over (a,b)
  !!
  !!		  ier    - integer
  !!                       ier = 0 normal and reliable termination of the
  !!                               routine. it is assumed that the requested
  !!                               accuracy has been achieved.
  !!                       ier.gt.0 abnormal termination of the routine. the
  !!                                estimates for result and error are less
  !!                                reliable. it is assumed that the requested
  !!                                accuracy has not been achieved.
  !!
  !!***references  (none)
  !!***routines called  sd1mach
  !!***end prologue  dqk15i
  subroutine dqk15i(fx,fx_vars,boun,inf,a,b,result,abserr,resabs,resasc, ier)

    ! Arguments
    interface
      function fx(centr, vars)
        use carma_precision_mod, only : f
        use adgaquad_types_mod
        real(kind=f), intent(in) :: centr
        type(adgaquad_vars_type), intent(inout) :: vars
        real(kind=f)             :: fx
      end function fx
    end interface
    type(adgaquad_vars_type) :: fx_vars
    real(kind=f) :: boun
    integer :: inf
    real(kind=f) :: a 
    real(kind=f) :: b
    real(kind=f) :: result
    real(kind=f) :: abserr
    real(kind=f) :: resabs
    real(kind=f) :: resasc
    integer :: ier

    ! Local declarations
    real(kind=f) :: absc, absc1, absc2, centr, dabs, dinf
    real(kind=f) :: dmax1, dmin1, epmach, fc, fsum, fval1, fval2, fv1(7) ,fv2(7), hlgth
    real(kind=f) :: resg, resk, reskh, tabsc1, tabsc2, uflow, wg(8), wgk(8), xgk(8)
    integer :: j

    !
    !           the abscissae and weights are supplied for the interval
    !           (-1,1).  because of symmetry only the positive abscissae and
    !           their corresponding weights are given.
    !
    !           xgk    - abscissae of the 15-point kronrod rule
    !                    xgk(2), xgk(4), ... abscissae of the 7-point
    !                    gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 7-point gauss rule
    !
    !           wgk    - weights of the 15-point kronrod rule
    !
    !           wg     - weights of the 7-point gauss rule, corresponding
    !                    to the abscissae xgk(2), xgk(4), ...
    !                    wg(1), wg(3), ... are set to zero.
    !
    data wg(1) / 0.0_f /
    data wg(2) / 0.129484966168869693270611432679082_f /
    data wg(3) / 0.0_f /
    data wg(4) / 0.279705391489276667901467771423780_f /
    data wg(5) / 0.0_f /
    data wg(6) / 0.381830050505118944950369775488975_f /
    data wg(7) / 0.0_f /
    data wg(8) / 0.417959183673469387755102040816327_f /

    data xgk(1) / 0.991455371120812639206854697526329_f /
    data xgk(2) / 0.949107912342758524526189684047851_f /
    data xgk(3) / 0.864864423359769072789712788640926_f /
    data xgk(4) / 0.741531185599394439863864773280788_f /
    data xgk(5) / 0.586087235467691130294144838258730_f /
    data xgk(6) / 0.405845151377397166906606412076961_f /
    data xgk(7) / 0.207784955007898467600689403773245_f /
    data xgk(8) / 0.000000000000000000000000000000000_f /

    data wgk(1) / 0.022935322010529224963732008058970_f /
    data wgk(2) / 0.063092092629978553290700663189204_f /
    data wgk(3) / 0.104790010322250183839876322541518_f /
    data wgk(4) / 0.140653259715525918745189590510238_f /
    data wgk(5) / 0.169004726639267902826583426598550_f /
    data wgk(6) / 0.190350578064785409913256402421014_f /
    data wgk(7) / 0.204432940075298892414161999234649_f /
    data wgk(8) / 0.209482141084727828012999174891714_f /
    !
    !
    !           list of major variables
    !           -----------------------
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc*  - abscissa
    !           tabsc* - transformed abscissa
    !           fval*  - function value
    !           resg   - result of the 7-point gauss formula
    !           resk   - result of the 15-point kronrod formula
    !           reskh  - approximation to the mean value of the transformed
    !                    integrand over (a,b), i.e. to i/(b-a)
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach is the largest relative spacing.
    !           uflow is the smallest positive magnitude.
    !
    !*** first executable statement  dqk15i
    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    call sd1mach(1,uflow,ier)
    if(ier.eq.9) return

    !epmach = d1mach(4)
    !uflow = d1mach(1)
    dinf = min0(1,inf)

    centr = 0.5_f*(a+b)
    hlgth = 0.5_f*(b-a)
    tabsc1 = boun+dinf*(0.1e1_f-centr)/centr
    fval1 = fx(tabsc1, fx_vars)
    if(inf.eq.2) fval1 = fval1+fx(-tabsc1, fx_vars)
    fc = (fval1/centr)/centr
    !
    !           compute the 15-point kronrod approximation to
    !           the integral, and estimate the error.
    !
    resg = wg(8)*fc
    resk = wgk(8)*fc
    resabs = dabs(resk)
    do j=1,7
      absc = hlgth*xgk(j)
      absc1 = centr-absc
      absc2 = centr+absc
      tabsc1 = boun+dinf*(0.1d+01-absc1)/absc1
      tabsc2 = boun+dinf*(0.1d+01-absc2)/absc2
      fval1 = fx(tabsc1, fx_vars)
      fval2 = fx(tabsc2, fx_vars)
      if(inf.eq.2) fval1 = fval1+fx(-tabsc1, fx_vars)
      if(inf.eq.2) fval2 = fval2+fx(-tabsc2, fx_vars)
      fval1 = (fval1/absc1)/absc1
      fval2 = (fval2/absc2)/absc2
      fv1(j) = fval1
      fv2(j) = fval2
      fsum = fval1+fval2
      resg = resg+wg(j)*fsum
      resk = resk+wgk(j)*fsum
      resabs = resabs+wgk(j)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5_f
    resasc = wgk(8)*dabs(fc-reskh)
    do j=1,7
      resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resasc = resasc*hlgth
    resabs = resabs*hlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc.ne.0.0_f.and.abserr.ne.0._f) abserr = resasc*dmin1(0.1e1_f,(0.2e3_f*abserr/resasc)**1.5_f)
    if(resabs.gt.uflow/(0.5e2_f*epmach)) abserr = dmax1((epmach*0.5e2_f)*resabs,abserr)
999 return
  end subroutine dqk15i

  !!
  !!***begin prologue  dqelg
  !!***refer to  dqagie,dqagoe,dqagpe,dqagse
  !!***routines called  sd1mach
  !!***revision date  130319   (yymmdd)
  !!***keywords  epsilon algorithm, convergence acceleration,
  !!             extrapolation
  !!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !!           de doncker,elise,appl. math & progr. div. - k.u.leuven
  !!***purpose  the routine determines the limit of a given sequence of
  !!            approximations, by means of the epsilon algorithm of
  !!            p.wynn. an estimate of the absolute error is also given.
  !!            the condensed epsilon table is computed. only those
  !!            elements needed for the computation of the next diagonal
  !!            are preserved.
  !!***description
  !!
  !!           epsilon algorithm
  !!           standard fortran subroutine
  !!           double precision version
  !!
  !!           parameters
  !!              n      - integer
  !!                       epstab(n) contains the new element in the
  !!                       first column of the epsilon table.
  !!
  !!              epstab - double precision
  !!                       vector of dimension 52 containing the elements
  !!                       of the two lower diagonals of the triangular
  !!                       epsilon table. the elements are numbered
  !!                       starting at the right-hand corner of the
  !!                       triangle.
  !!
  !!              result - double precision
  !!                       resulting approximation to the integral
  !!
  !!              abserr - double precision
  !!                       estimate of the absolute error computed from
  !!                       result and the 3 previous results
  !!
  !!              res3la - double precision
  !!                       vector of dimension 3 containing the last 3
  !!                       results
  !!
  !!              nres   - integer
  !!                       number of calls to the routine
  !!                       (should be zero at first call)
  !!
  !!		  ier    - integer
  !!                       ier = 0 normal and reliable termination of the
  !!                               routine. it is assumed that the requested
  !!                               accuracy has been achieved.
  !!                       ier.gt.0 abnormal termination of the routine. the
  !!                                estimates for result and error are less
  !!                                reliable. it is assumed that the requested
  !!                                accuracy has not been achieved.
  !!
  !!***end prologue  dqelg
  !!
  subroutine dqelg(n,epstab,result,abserr,res3la,nres,ier)

    ! Arguments
    integer :: n
    real(kind=f) :: epstab(52)
    real(kind=f) :: result
    real(kind=f) :: abserr
    real(kind=f) :: res3la(3)
    integer :: nres  
    integer :: ier

    ! Local declarations
    real(kind=f) :: dabs, delta1, delta2, delta3, dmax1
    real(kind=f) :: epmach, epsinf, error, err1, err2, err3, e0, e1, e1abs, e2, e3
    real(kind=f) :: oflow, res, ss, tol1, tol2, tol3
    integer :: i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num
    !
    !           list of major variables
    !           -----------------------
    !
    !           e0     - the 4 elements on which the computation of a new
    !           e1       element in the epsilon table is based
    !           e2
    !           e3                 e0
    !                        e3    e1    new
    !                              e2
    !           newelm - number of elements to be computed in the new
    !                    diagonal
    !           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
    !           result - the element in the new diagonal with least value
    !                    of error
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach is the largest relative spacing.
    !           oflow is the largest positive magnitude.
    !           limexp is the maximum number of elements the epsilon
    !           table can contain. if this number is reached, the upper
    !           diagonal of the epsilon table is deleted.
    !
    !***first executable statement  dqelg
    call sd1mach(4,epmach,ier)
    if(ier.eq.9) return
    call sd1mach(2,oflow,ier)
    if(ier.eq.9) return

    !epmach = d1mach(4)
    !oflow = d1mach(2)
    nres = nres+1
    abserr = oflow
    result = epstab(n)
    if(n.lt.3) go to 100
    limexp = 50
    epstab(n+2) = epstab(n)
    newelm = (n-1)/2
    epstab(n) = oflow
    num = n
    k1 = n
    do 40 i = 1,newelm
      k2 = k1-1
      k3 = k1-2
      res = epstab(k1+2)
      e0 = epstab(k3)
      e1 = epstab(k2)
      e2 = res
      e1abs = dabs(e1)
      delta2 = e2-e1
      err2 = dabs(delta2)
      tol2 = dmax1(dabs(e2),e1abs)*epmach
      delta3 = e1-e0
      err3 = dabs(delta3)
      tol3 = dmax1(e1abs,dabs(e0))*epmach
      if(err2.gt.tol2.or.err3.gt.tol3) go to 10
      !
      !           if e0, e1 and e2 are equal to within machine
      !           accuracy, convergence is assumed.
      !           result = e2
      !           abserr = abs(e1-e0)+abs(e2-e1)
      !
      result = res
      abserr = err2+err3
      ! ***jump out of do-loop
      go to 100
 10   e3 = epstab(k1)
      epstab(k1) = e1
      delta1 = e1-e3
      err1 = dabs(delta1)
      tol1 = dmax1(e1abs,dabs(e3))*epmach
      !
      !           if two elements are very close to each other, omit
      !           a part of the table by adjusting the value of n
      !
      if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
      ss = 0.1e1_f/delta1+0.1e1_f/delta2-0.1e1_f/delta3
      epsinf = dabs(ss*e1)
      !
      !           test to detect irregular behaviour in the table, and
      !           eventually omit a part of the table adjusting the value
      !           of n.
      !
      if(epsinf.gt.0.1e-3_f) go to 30
 20   n = i+i-1
      ! ***jump out of do-loop
      go to 50
      !
      !           compute a new element and eventually adjust
      !           the value of result.
      !
 30   res = e1+0.1e1_f/ss
      epstab(k1) = res
      k1 = k1-2
      error = err2+dabs(res-e2)+err3
      if(error.gt.abserr) go to 40
      abserr = error
      result = res
 40 continue
    !
    !           shift the table.
    !
 50 if(n.eq.limexp) n = 2*(limexp/2)-1
    ib = 1
    if((num/2)*2.eq.num) ib = 2
    ie = newelm+1
    do 60 i=1,ie
      ib2 = ib+2
      epstab(ib) = epstab(ib2)
      ib = ib2
 60 continue
    if(num.eq.n) go to 80
    indx = num-n+1
    do 70 i = 1,n
      epstab(i)= epstab(indx)
      indx = indx+1
 70 continue
 80 if(nres.ge.4) go to 90
    res3la(nres) = result
    abserr = oflow
    go to 100
    !
    !           compute error estimate
    !
 90 abserr = dabs(result-res3la(3))+dabs(result-res3la(2))+dabs(result-res3la(1))
    res3la(1) = res3la(2)
    res3la(2) = res3la(3)
    res3la(3) = result
100 abserr = dmax1(abserr,0.5e1_f*epmach*dabs(result))
999 return
  end subroutine dqelg
  
  !!
  !! *********************************************************
  !! taken from BLAS library
  !!  (http://netlib.bell-labs.com/netlib/blas)
  !! *********************************************************
  SUBROUTINE SD1MACH(I,D1MACH_OUT,IER)
    INTEGER, INTENT(in) :: I
    REAL(kind=f), INTENT(out) :: D1MACH_OUT
    INTEGER, INTENT(out) :: IER
  !
  !  DOUBLE-PRECISION MACHINE CONSTANTS
  !  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
  !  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
  !  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
  !  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
  !  D1MACH( 5) = LOG10(B)
  !
    INTEGER :: SMALL(2)
    INTEGER :: LARGE(2)
    INTEGER :: RIGHT(2)
    INTEGER :: DIVER(2)
    INTEGER :: LOG10(2)
    INTEGER :: SC, CRAY1(38), J
    SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
    REAL(kind=f) :: DMACH(5)
    EQUIVALENCE (DMACH(1),SMALL(1))
    EQUIVALENCE (DMACH(2),LARGE(1))
    EQUIVALENCE (DMACH(3),RIGHT(1))
    EQUIVALENCE (DMACH(4),DIVER(1))
    EQUIVALENCE (DMACH(5),LOG10(1))
    !  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
    !  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
    !  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
    !  MANY MACHINES YET.
    !  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
    !  ON THE NEXT LINE
    DATA SC/0/
    !  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
    !  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
    !          mail netlib@research.bell-labs.com
    !          send old1mach from blas
    !  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
    !
    !     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
    !      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
    !      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
    !      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
    !      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
    !      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
    !
    !     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
    !     32-BIT INTEGERS.
    !      DATA SMALL(1),SMALL(2) /    8388608,           0 /
    !      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
    !      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
    !      DATA DIVER(1),DIVER(2) /  620756992,           0 /
    !      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
    !
    !  MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
    !      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
    !      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
    !      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
    !      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
    !      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
    !
    !     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
    IER = 0
    IF (SC .NE. 987) THEN
      DMACH(1) = 1.e13_f
      IF (      SMALL(1) .EQ. 1117925532 .AND. SMALL(2) .EQ. -448790528) THEN
      !           *** IEEE BIG ENDIAN ***
        SMALL(1) = 1048576
        SMALL(2) = 0
        LARGE(1) = 2146435071
        LARGE(2) = -1
        RIGHT(1) = 1017118720
        RIGHT(2) = 0
        DIVER(1) = 1018167296
        DIVER(2) = 0
        LOG10(1) = 1070810131
        LOG10(2) = 1352628735
      ELSE IF ( SMALL(2) .EQ. 1117925532 .AND. SMALL(1) .EQ. -448790528) THEN
      !           *** IEEE LITTLE ENDIAN ***
        SMALL(2) = 1048576
        SMALL(1) = 0
        LARGE(2) = 2146435071
        LARGE(1) = -1
        RIGHT(2) = 1017118720
        RIGHT(1) = 0
        DIVER(2) = 1018167296
        DIVER(1) = 0
        LOG10(2) = 1070810131
        LOG10(1) = 1352628735
      ELSE IF ( SMALL(1) .EQ. -2065213935 .AND. SMALL(2) .EQ. 10752) THEN
      !               *** VAX WITH D_FLOATING ***
        SMALL(1) = 128
        SMALL(2) = 0
        LARGE(1) = -32769
        LARGE(2) = -1
        RIGHT(1) = 9344
        RIGHT(2) = 0
        DIVER(1) = 9472
        DIVER(2) = 0
        LOG10(1) = 546979738
        LOG10(2) = -805796613
      ELSE IF ( SMALL(1) .EQ. 1267827943 .AND. SMALL(2) .EQ. 704643072) THEN
      !               *** IBM MAINFRAME ***
        SMALL(1) = 1048576
        SMALL(2) = 0
        LARGE(1) = 2147483647
        LARGE(2) = -1
        RIGHT(1) = 856686592
        RIGHT(2) = 0
        DIVER(1) = 873463808
        DIVER(2) = 0
        LOG10(1) = 1091781651
        LOG10(2) = 1352628735
      ELSE IF ( SMALL(1) .EQ. 1120022684 .AND. SMALL(2) .EQ. -448790528) THEN
      !           *** CONVEX C-1 ***
        SMALL(1) = 1048576
        SMALL(2) = 0
        LARGE(1) = 2147483647
        LARGE(2) = -1
        RIGHT(1) = 1019215872
        RIGHT(2) = 0
        DIVER(1) = 1020264448
        DIVER(2) = 0
        LOG10(1) = 1072907283
        LOG10(2) = 1352628735
      ELSE IF ( SMALL(1) .EQ. 815547074 .AND. SMALL(2) .EQ. 58688) THEN
      !           *** VAX G-FLOATING ***
        SMALL(1) = 16
        SMALL(2) = 0
        LARGE(1) = -32769
        LARGE(2) = -1
        RIGHT(1) = 15552
        RIGHT(2) = 0
        DIVER(1) = 15568
        DIVER(2) = 0
        LOG10(1) = 1142112243
        LOG10(2) = 2046775455
      ELSE
        DMACH(2) = 1.e27_f + 1
        DMACH(3) = 1.e27_f
        LARGE(2) = LARGE(2) - RIGHT(2)
        IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
          CRAY1(1) = 67291416
          DO J = 1, 20
            CRAY1(J+1) = CRAY1(J) + CRAY1(J)
          END DO
          CRAY1(22) = CRAY1(21) + 321322
          DO J = 22, 37
            CRAY1(J+1) = CRAY1(J) + CRAY1(J)
          END DO
          IF (CRAY1(38) .EQ. SMALL(1)) THEN
          !                  *** CRAY ***
            CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
            SMALL(2) = 0
            CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
            CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
            CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
            RIGHT(2) = 0
            CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
            DIVER(2) = 0
            CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
            CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
          ELSE
            IER=9
          END IF
          ELSE
            IER=9
          END IF
        END IF
        SC = 987
      END IF
      !    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) IER=9
      IF (I .LT. 1 .OR. I .GT. 5) THEN
        IER=9
      END IF
      D1MACH_OUT = DMACH(I)
      RETURN
  END SUBROUTINE SD1MACH

  SUBROUTINE I1MCRY(A, A1, B, C, D)
    !*** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
    INTEGER A, A1, B, C, D
    A1 = 16777216*B + C
    A = 16777216*A1 + D
  END SUBROUTINE I1MCRY



end module adgaquad_mod
