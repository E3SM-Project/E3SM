!! This module (fractal_meanfield_mod.F90) contains the main routines 
!! necessary to calculate the solution of the mean field approximation 
!! for a dry fractal particle composed of identical spherical monomers.
!! This is used to generate optical properties for these paticles in CARMA.
!!
!! See Botet et al. 1997 "Mean-field approximation of Mie
!! scattering by fractal aggregates of identical spheres."
!! Applied Optics 36(33) 8791-8797
!!
!! Original code from P. Rannou and R. Botet.
!! Translated to F90 and ported into CARMA by E. Wolf
!!
!! master: fractal_meanfield calling: cmie,ludcmpc,lubksbc,dqagi
!!
!! calculating the monomer Mie scattering
!!   - SUBROUTINE cmie()      calling: intmie()
!!   - SUBROUTINE intmie()    calling: intmie()
!!
!! calculating the matrix elements
!!   - FUNCTION funa()        calling: dqag,fpl
!!   - FUNCTION fpl()         calling: plgndr
!!   - FUNCTION plgndr()
!!   - FUNCTION funb_n()
!!   - FUNCTION funs_n()      calling: dq2agi,xfreal_n,xfimag_n
!!   - FUNCTION xfreal_n()    calling: besseljy,phi
!!   - FUNCTION xfimag_n()    calling: besseljy,phi
!!   - FUNCTION BESSELJY()
!! 
!! Routines to calculate the scattered wave
!! of monomer:
!!   - FUNCTION fpi()         calling: plgndr()
!!   - FUNCTION ftau()        calling: plgndr()
!! of agglomerate/cluster:
!!   - FUNCTION fp1()
!!
!! Routines related to the probability distribution:
!!   - FUNCTION anorm()       calling: dqdagi,fdval
!!   - FUNCTION fdval()
!!   - FUNCTION phi()         calling: fdval
!!   - FUNCTION fco()         calling: fdval
!!
!! @author P. Rannou, R. Botet, Eric Wolf
!! version March 2013
module fractal_meanfield_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carma_mod
  
  use adgaquad_types_mod
  use adgaquad_mod
  use lusolvec_mod

  implicit none
  
  private

  public :: fractal_meanfield

  ! Private module varibles: Moved from COMMON blocks
  integer, parameter :: nmi=40
  integer, parameter :: n2m = 2*nmi
  
  contains

  !!
  !! Generate optical properties for CARMA fractal particles.
  !!
  !! See Botet et al. 1997 "Mean-field approximation of Mie
  !! scattering by fractal aggregates of identical spheres."
  !! Applied Optics 36(33) 8791-8797
  !!
  !! @author  P.Rannou, R.Botet, Eric Wolf
  !! @version March 2013
  subroutine fractal_meanfield(carma, xl_in, xk_in, xn_in, nb_in, alpha_in, &
    df_in, rmon,xv, ang, Qext, Qsca, gfac, rc)

    ! some of these may be included in carma object
    type(carma_type), intent(in)       :: carma         !! the carma object
    real(kind=f),intent(in)            :: xl_in         !! Wavelength [microns]
    real(kind=f),intent(in)            :: xk_in         !! imaginary index of refraction
    real(kind=f),intent(in)            :: xn_in         !! real index of refraction
    real(kind=f),intent(in)            :: nb_in         !! number of monomers
    real(kind=f),intent(in)            :: alpha_in      !! Packing coefficient
    real(kind=f),intent(in)            :: df_in         !! Fractal dimension
    real(kind=f),intent(in)            :: rmon          !! monomer size [microns]
    real(kind=f),intent(in)            :: xv            !! set to 1
    real(kind=f),intent(in)            :: ang           !! angle set to zero
    real(kind=f),intent(out)           :: Qext          !! EFFICIENCY FACTOR FOR EXTINCTION
    real(kind=f),intent(out)           :: Qsca          !! EFFICIENCY FACTOR FOR SCATTERING
    real(kind=f),intent(out)           :: gfac          !! asymmetry factor
    integer,intent(inout)              :: rc            !! return code, negative indicates failure

    ! Local declarations
    integer, parameter                 :: nth = 10001
    integer, parameter                 :: maxsub = 50000
    integer, parameter                 :: lenw = 200000
    integer                            :: nstop         ! index with the last mie-coefficient
    integer                            :: n1stop
    integer                            :: pp,tt,mm
    real(kind=f)                       :: krg,rg        ! Particle structure
    real(kind=f)                       :: sigmas,sigmae,nc1
    real(kind=f)                       :: sigext        ! extinction cross section
    real(kind=f)                       :: sigsca        ! scattering cross section
    real(kind=f)                       :: sigabs        ! absorption cross section
    real(kind=f)                       :: totg          ! asymmetry parameter 
    real(kind=f)                       :: sigext2,sigext3
    real(kind=f)                       :: rems          ! radius of equivalent mass sphere
    real(kind=f)                       :: gems          ! geometric cross-section of equivalent mass sphere
    real(kind=f)                       :: dthetar,angler,weight
    real(kind=f)                       :: sumsca
    real(kind=f)                       :: lbd,beta      ! optical characteristics
    real(kind=f)                       :: sstest(0:nf)
    real(kind=f)                       :: setest(0:nf)
    real(kind=f)                       :: xl(39)        ! place holder for  wavelength
    real(kind=f)                       :: xn(39)        ! place holder for real index of refraction
    real(kind=f)                       :: xk(39)        ! place holder for imaginary index of refraction
    real(kind=f)                       :: val           
    real(kind=f)                       :: funca(nmi,nmi,0:n2m)       ! for storage of funa(nu,n;p)  
    complex(kind=f)                    :: res           ! for storage of funs_n
    complex(kind=f)                    :: funcs(0:n2m)
    real(kind=f)                       :: s11(0:nth-1)
    real(kind=f)                       :: s11_n(0:nth-1)
    real(kind=f)                       :: xint(0:nth-1)
    real(kind=f)                       :: wom
    real(kind=f)                       :: pol(0:nth-1)
    complex(kind=f)                    :: s01,s02
    complex(kind=f)                    :: s1(0:nth-1)
    complex(kind=f)                    :: s2(0:nth-1)
    complex(kind=f)                    :: ajt
    complex(kind=f)                    :: an(nf)   
    complex(kind=f)                    :: bn(nf)
    complex(kind=f)                    :: ni,i,id,onec,zeroc
    complex(kind=f)                    :: d1(nmi)
    complex(kind=f)                    :: d2(nmi)
    complex(kind=f)                    :: Ap1(nmi,nmi)
    complex(kind=f)                    :: Bp1(nmi,nmi)
    complex(kind=f)                    :: dvec(n2m)  ! For matrix eqn of order 2N
    complex(kind=f)                    :: cvec(n2m)
    complex(kind=f)                    :: EpABC(n2m,n2m)
    integer                            :: luindx(n2m)  ! For LU decomposition
    real(kind=f)                       :: dlu
    integer                            :: ifail 
    integer                            :: iwork(maxsub)
    integer                            :: neval,nsubin
    real(kind=f)                       :: work(lenw)

    ! Previously these were implicitly defined
    real(kind=f)                       :: angle, pi, rn, ri
    real(kind=f)                       :: deltas, deltae, xfact
    real(kind=f)                       :: bound, errrel, p1, dp1
    real(kind=f)                       :: errabs, total
    integer                            :: n2stop, n3stop, ntheta, ii, kk, nn, jj, iy, ir, q, interv
    real(kind=f)                       :: a0, c0, a1, c1, a2, c2
    real(kind=f)                       :: qabs
    
    ! Previously these were globals, which wouldn't be thread safe.
    type(adgaquad_vars_type)           :: fx_vars

    ! Set the return code to default to okay.
    rc = RC_OK

    ! *** Set from input arguments
    fx_vars%nb = nb_in
    fx_vars%df = df_in
    fx_vars%alpha = alpha_in
    xl(1) = xl_in
    xk(1) = xk_in
    xn(1) = xn_in    

    ! *** Complex constants 1, 1, identity(1,1), zero(0,0) :
    i     = cmplx(0._f,1._f,kind=f)
    onec  = cmplx(1._f,0._f,kind=f)
    id    = cmplx(1._f,1._f,kind=f)
    zeroc = cmplx(0._f,0._f,kind=f)

    ! Other initializations
    funca(:,:,:) = 0.0_f
    fx_vars%a = rmon *1.e-2_f           ! a = r_monomer in m
    beta=ang*(3.1415926_f / 180._f)     ! =0 when ang=0  
    Ap1(:,:) = zeroc
    Bp1(:,:) = zeroc
    sstest(:) = 0.0_f
    setest(:) = 0.0_f

    ! **************************************************************** 
    ! *** Definition and calculation  of factorials 0 - nf   
    ! *** (nf set in adgaqaud_types_mod.F90)
    ! *** and storage   [ real*8   fact()   (double prec.) ] 
    ! ****************************************************************
  
    fx_vars%fact(0)=1._f      ! factorials fact(n)=n!
    do ii=1,nf
      fx_vars%fact(ii) = fx_vars%fact(ii-1)*ii*1._f
    end do

    pi=4._f*atan(1._f)       ! 3.1415926535
    fx_vars%coeff=anorm(carma,fx_vars,rc)
    if (rc < 0) return

    ! ****************************************************************   
    ! anorm() integrated INT_0^inf[ x**(df-1.)*exp(-x**df/2._f)  dx ]
    ! and occupied 
    !    anorm := 4 pi * INT_0^inf[ x**(df-1.)*exp(-x**df/2._f)  dx ]
    !          == geometric scalingfactor Eq.(10) in [Botet et al, 1995]
    !    c     := 0.5
    ! ****************************************************************  
 
    kk=1 
    ni=xn(kk)*1._f+i*xk(kk)*xv*1._f      ! ni := complex index of refraction of monomer
                                         ! (xv := 1 ; input parameter in file "calpha")
    lbd=xl(kk)*1.e-6_f                   ! lbd := wavelength in m
                                         !        (in matrix medium / material !)
    fx_vars%k=2._f*pi/lbd                ! k   := abs.val. of wavevector in m^-1 
                                         !        (in matrix medium / material !)
  
    ! *** ******************************************************************
    ! *** Calculation of Mie coefficients for monomer scattering             
    ! *** up to a maximum order of nf=50                                        
    ! *** ******************************************************************

    do ii=1,nf
      an(ii) = zeroc
      bn(ii) = zeroc
    end do

    rn=xn(kk)        ! Re(relative_n_complex,monomer)
    ri=xk(kk)*xv     ! Im(relative_n_complex,monomer)
                      ! xv should be set to 1 (sse above)              
                      ! a = monomer sphere radius                      
                      ! lbd = wavelength in matrix medium              

    ! Call Mie routine
    call cmie(lbd,rn,ri,fx_vars%a,an,bn,nstop)
    
    do ii=1,nf
      if (an(ii).ne.0._f) nstop=ii
    end do

    ! nstop is now the index with the last mie-coefficient
    ! (highest index i) an(i) not equal zero.
    ! since all the an were set to zero before calling
    ! cmie(), nstop is the termination index used in cmie()
    ! or in intmie(). Usually, a termination index
    ! nstop = INT( 2 + x + 4 x^(1/3) ) is used; in intmie(),
    ! however, a value of
    ! nstop := MAX( INT(...), |m*x| )+15 is used !???

    sigmas=0._f
    sigmae=0._f

    do nn=1,nstop
      nc1=abs(an(nn))**2._f+abs(bn(nn))**2._f
      nc1=nc1*(2._f*nn+1._f)
      sigmas=sigmas+nc1*(2._f*3.14159265_f)/(fx_vars%k**2._f)
      nc1=real(an(nn)+bn(nn))
      nc1=nc1*(2._f*nn+1._f)
      sigmae=sigmae+nc1*(2._f*3.14159265_f)/(fx_vars%k**2._f)
      sstest(nn)=sigmas
      setest(nn)=sigmae
      deltas=abs(sstest(nn-1)-sstest(nn))/sstest(nn)
      deltae=abs(setest(nn-1)-setest(nn))/setest(nn)
      if(deltas.gt.1.e-6_f) n2stop=nn
      if(deltae.gt.1.e-6_f) n3stop=nn
    end do

    n1stop=n2stop
    if (n3stop.gt.n2stop) n1stop=n3stop
    ! The order of the set of linear equations is chosen
    ! as the number of mie coefficients where the sum yielding
    ! the monomer ext./scatt. cross sections do not change more
    ! than 1.D-3 compared to the values with one summand less.

    rg=fx_vars%alpha*fx_vars%nb**(1._f/fx_vars%df)*fx_vars%a    ! rg := radius of gyration
    krg=fx_vars%k*rg
    ntheta=7800                 !180-int(krg**.5*28*log10(dxk(kk)))               
    if (ntheta*0.5_f .eq. (ntheta/2)*1._f) ntheta=ntheta+1

    ! *** ******************************************************************
    ! *** FIRST PART:  Solution of self consistent mean field equation,
    ! ***              i.e. the set of linear equations (SLE) defining the
    ! ***              mean field coefficients d^1_1,n and d^2_1,n
    ! ***              according to Eq.(12) of (Botet 1997)
    ! ***              To do so,
    ! ***               - matrix elements A^1,nu_1,n and B^1,nu_1,n
    ! ***                 are calculated with Eq.(13), using Eqns.(14)-(16)
    ! ***               - the set of lin. Eqns. is solved yielding the d's
    ! *** ******************************************************************
    ! *** Eq.(12) of Botet et al 1997 defines a matrix eqn. of order 2N :
    ! *** (since N=n1stop, 2N = 2 * n1stop = order of SLE)
    ! ***                                                                       
    ! ***      EpABC * dvec = cvec     with                                     
    ! ***                                                                       
    ! *** dvec and cvec  being the 2N-vectors                                   
    ! ***                                                                       
    ! ***              ( d^(1)_1,1 )                 ( a_1 )                    
    ! ***              ( d^(1)_1,2 )                 ( a_2 )                    
    ! ***              (   ...     )                 ( ... )                    
    ! ***              ( d^(1)_1,n )                 ( a_n )                    
    ! ***      dvec := ( d^(2)_1,1 )   and   cvec := ( b_1 )    and further     
    ! ***              ( d^(2)_1,2 )                 ( b_2 )                    
    ! ***              (   ...     )                 ( ... )                    
    ! ***              ( d^(2)_1,n )                 ( b_n )                    
    ! ***                                                        
    ! ***                                                                       
    ! ***    EpABC :=  1 + AB * C    where AB, 1, C are the 2N*2N - matrices    
    ! ***                                                                       
    ! ***          ( a_1 0  0 ...          ... 0 )                              
    ! ***          (  0 a_2 0 ...          ... 0 )                              
    ! ***    AB := (  0  0  ...  0  0  0   ... 0 )   1 := 2N*2N unity matrix    
    ! ***          (  0  0  ... a_n 0  0   ... 0 )                              
    ! ***          (  0  0  ...  0 b_1 0   ... 0 )                              
    ! ***          (  ...             ...  ... 0 )                              
    ! ***          (  ...          ... 0 b_n-1 0 )                              
    ! ***          (  0  0  ...    ... 0   0  b_n)                              
    ! ***                                                                       
    ! *** and      (  A   B  )                                                  
    ! ***    C  := (         )  where A and B are the two N*N matrices          
    ! ***          (  B   A  )  given by Eq.(13),                               
    ! ***                       including the factor (N_monomers - 1):          
    ! ***                                                                       
    ! ***    A_n,nu :=  (N_m-1) * A_(1,n)^(1,nu)  and                           
    ! ***    B_n,nu :=  (N_m-1) * B_(1,n)^(1,nu)                                
    ! ***                                                                       
    ! ***    (A_(1,n... and B_(1,n... according to Eq.(13) of Botet 1997)       
    ! *** ******************************************************************  

    n2stop = 2 * n1stop       ! n2stop = order of SLE

    ! *** ******************************************************************  
    ! *** Error handling moved from xfreal_n, xfimag_n.  Calculations fail 
    ! *** in integration package when n2stop > 48. n2stop is related to the 
    ! *** number of complex mie scattering coefficients used in teh calculation
    ! *** which is in turn related to the size parameter of monomers.
    ! *** If nstop>48 end calculation here instead of continuing.

    if (n2stop.gt.48) then
      if (carma%f_do_print) then
         write(carma%f_LUNOPRT, *) "fractal_meanfield_mod::n2stop greater &
              &than 48. Size parameter (2*pi*rmon/lambda): ", &
              2._f*3.14159265_f*fx_vars%a/lbd, "Monomer Size parameter &
              &must be less than ~17."
      end if
      rc = RC_ERROR
      return
    endif


    ! *** ******************************************************************
    do ii=1,n1stop
      cvec(ii)        = an(ii)       ! right hand side vector
      cvec(n1stop+ii) = bn(ii)
      dvec(ii)        = zeroc        ! solution vector d
      dvec(n1stop+ii) = zeroc
    end do
  
    do pp=0,n2stop       !variable p
      res = funs_n(carma,fx_vars,pp,rc)   ! Eq.(16)    S_p(k R_g)
      if (rc < 0) return
      funcs(pp) = res
    end do
    
    ! Calculate terms A and B

    ! *** loops over indices nu,n,p :                                              
    do ii=1,n1stop       !variable n
      do jj=1,n1stop     !variable nu
        mm = IABS(ii-jj)
        tt = ii+jj

        ! *** ******************************************************************       
        !        calculation of A_(1,n)^(1,nu) according to Eq.(13)                    
        ! *** ******************************************************************       
        do pp=mm,tt

          funca(jj,ii,pp) = funa(carma,fx_vars,jj,ii,pp,rc)       ! Eq.(14)    a(nu,n;p)a
          if (rc < 0) return

          Ap1(ii,jj) = Ap1(ii,jj) + &
                       ( onec * (ii*(ii+1)+jj*(jj+1)-pp*(pp+1)) ) &
                                 * funca(jj,ii,pp) * funcs(pp)
        end do  ! loop over pp (variable p)                            

        ! scaling factors of eq.(13), factor (N_mon-1) from eq.(12)
        Ap1(ii,jj) = Ap1(ii,jj) * (2._f*jj+1._f)/(jj*(jj*1._f+1._f))
        Ap1(ii,jj) = Ap1(ii,jj) * (fx_vars%nb-1._f) / (ii*(ii*1._f+1._f))

        ! *** ******************************************************************
        !        calculation of B_(1,n)^(1,nu) according to Eq.(13)
        ! *** ******************************************************************
        do pp=mm,tt
          Bp1(ii,jj) = Bp1(ii,jj) + funb_n(jj,ii,pp,funca) * funcs(pp)
        end do      ! loop over pp (variable p)

        ! scaling factors of eq.(13), factor (N_mon-1) from eq.(12)
        Bp1(ii,jj) = Bp1(ii,jj) * (2._f*jj+1._f)/(jj*(jj*1._f+1._f))
        Bp1(ii,jj) = Bp1(ii,jj) * (fx_vars%nb-1._f) * 2._f/(ii*(ii*1._f+1._f))
      end do  ! loop over jj=1,n1stop (variable nu)                         
    end do ! loop over ii=1,n1stop (variable n) 

    ! *** ******************************************************************      
    ! End of Calculation of terms A and B

    ! *** ******************************************************************       
    ! *** Setup and solution of matrix equation of order 2N ( = n2stop )           
    ! *** constituted by eq.(12)                                                   
    ! *** ******************************************************************       
    ! *** matrix product (AB * C)  (definitions see above)                         
    do ii=1,n1stop
      do jj=1,n1stop
        EpABC(ii,jj)               = an(ii) * Ap1(ii,jj)
        EpABC(ii,jj+n1stop)        = an(ii) * Bp1(ii,jj)
        EpABC(ii+n1stop,jj)        = bn(ii) * Bp1(ii,jj)
        EpABC(ii+n1stop,jj+n1stop) = bn(ii) * Ap1(ii,jj)
      end do
    end do
    
    ! *** ******************************************************************       
    ! *** add 2N*2N unity matrix                                                   
    do ii=1,n1stop
      EpABC(ii,ii)               = EpABC(ii,ii) + onec
      EpABC(ii+n1stop,ii+n1stop) = EpABC(ii+n1stop,ii+n1stop) + onec
    end do
    
    ! ======================================================================       
    ! *** solve matrix equation using external routines (LU decomposition)         
    CALL LUDCMPC(EpABC,n2stop,n2m,luindx,dlu)
    CALL LUBKSBC(EpABC,n2stop,n2m,luindx,cvec)
    do ii=1,n1stop
      d1(ii) = cvec(ii)
      d2(ii) = cvec(ii+n1stop)
    end do

    ! *** ******************************************************************       
    ! *** SECOND PART:  Recomposition of the total wave scattered by               
    ! ***               the entire agglomerate/cluster by adding the               
    ! ***               waves scattered by each monomer taking into                
    ! ***               account the respective phase of the single waves.          
    ! *** ****************************************************************** 

    ! *** ******************************************************************       
    ! ----------------------------------------------------------------------       
    !  1) Calculate the amplitude functions |S1^j(th)|  et |S2^j(th)|              
    !     of one monomer of the agglomerate/cluster:                               
    !     ( see e.g. Bohren, Huffman (1983) p.112, Eq.(4.74) with the              
    !       substitutions a_n -> d^1_1,n  and b_n -> d^2_1,n                       
    !       or   Rannou (1999) Eq.(1)-(6)                        )                 
    ! ----------------------------------------------------------------------       
    ! *** ******************************************************************  

    do iy=0,ntheta-1,1          ! loop over angles
      angle=iy*180._f/(ntheta-1)
      s1(iy)=0._f
      s2(iy)=0._f
      wom=cos(angle*3.1415926353_f/180._f)

      do ir=1,n1stop             ! loop over Mie - indices
        xfact=2._f*(2._f*ir+1._f)/(ir*1._f*(ir*1._f+1._f))
        ajt=d1(ir)*fpi(ir,wom,fx_vars)+d2(ir)*tau(ir,wom,fx_vars)
        s1(iy)=s1(iy)+xfact*ajt
        ajt=d1(ir)*tau(ir,wom,fx_vars)+d2(ir)*fpi(ir,wom,fx_vars)
        s2(iy)=s2(iy)+xfact*ajt
      end do

      s11(iy)=abs(s1(iy))**2._f+abs(s2(iy))**2._f
      pol(iy)=abs(s1(iy))**2._f-abs(s2(iy))**2._f
      pol(iy)=pol(iy)/(abs(s1(iy))**2_f+abs(s2(iy))**2._f)
      ! ***    S_11(theta) = 1/2 * ( |S_1|^2 + |S_2|^2 )                             
      ! ***    above, s1(theta) = 2 * S_1(theta)                                     
      ! ***  =>S_11(theta) = 1/2 * ( |1/2*s1|^2 + |1/2*s2|^2 )                       
      !                     = 1/8 * ( |s1|^2 + |s2|^2 )                              
      s11_n(iy)=.125_f*(abs(s1(iy))**2._f+abs(s2(iy))**2._f)
    end do

    s01=s1(0)
    s02=s2(0)

    ! *** Extinction cross section sigext( d^1_1,n , d^2_1,n ) ***                 
    sigext=0._f
    do ir=1,n1stop                ! loop (sum) over Mie-indices
      sigext=sigext+(2._f*ir+1._f)*REAL(d1(ir)+d2(ir))
    end do
    sigext  = fx_vars%nb * 2._f*pi/fx_vars%k**2._f * sigext    ! Eq.(27)

    ! *** Alternatively (in a test, all values agreed with rel.acc. 1e-6),         
    ! *** Extinction cross section sigext( S(0 deg) ) (optical theorem) ***        
    ! *** (see e.g. Bohren, Huffman (1983), Eq. (4.76))                            
    ! *** S(0)=S_1(0)=S_2(0);  sigma_ext = 4 pi / k^2 * Re(S(0))                   
    ! *** above, s1(theta) = 2 * S_1(theta) (factor 2 in 'xfact')                  
    !     sigext2 = nb * 4._f*pi/k**2._f * 0.5_f*REAL(s01)                         
    !     sigext3 = nb * 4._f*pi/k**2._f * 0.5_f*REAL(s02)                         
    ! *** ******************************************************************    

    ! *** ******************************************************************       
    ! ----------------------------------------------------------------------       
    ! 2) Calculate the phase integral in Eq.(26) with P(r) already                 
    !    substituted ( compare Eq.(10) and (Botet 1995) ) :                        
    !         INT(0;infinity)[ sin(2XuZ) u^(d-2) f_co(u) du ]                      
    !    taking into account the different phases of the single                    
    !    scattered waves.                                                          
    ! ----------------------------------------------------------------------       
    ! *** ******************************************************************   
    do q=0,ntheta-1,1
      angle=q*180._f/(ntheta-1)
      if (angle .eq.   0._f) angle=0.001_f
      if (angle .eq. 180._f) angle=179.999_f
      fx_vars%zed=sin(angle*3.1415928353_f/180._f/2._f)

      bound=0._f
      interv=1
      errrel=1e-5_f
      p1=0._f
      dp1=0._f

      !======================================================================       
      !---    Version using the QUADPACK - routine :                                
      !----------------------------------------------------------------------       
      ifail = 0
      CALL dqagi(fp1,fx_vars,bound,interv,errabs,errrel,p1,dp1,neval,ifail,maxsub,lenw,nsubin,iwork,work)
      if(ifail.ne.0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "fractal_meanfield_mod::ifail=",ifail," returned by dqagi()!"
        rc = RC_ERROR
        return
      endif
      !======================================================================       

      p1=2._f*pi * (fx_vars%nb-1._f) / (fx_vars%coeff*fx_vars%zed*krg)*p1 + 1._f
      xint(q)=p1
    end do

    ! *** now, xint(theta) contains the square bracket terms in                    
    ! *** Botet (1997) Eq.(26)  or Rannou (1999) Eq.(1)            

    ! *** ******************************************************************       
    ! ----------------------------------------------------------------------       
    ! 3) Calculation of the phase function, calculation of the optical             
    !    properties (asymmetrie factor g, scatt. cross section sigma_s)            
    !    by angular integration:   INT_0^180[ ... d_theta ]                        
    ! ----------------------------------------------------------------------       
    ! *** ******************************************************************  

    total=0._f
    totg=0._f

    do q=1,ntheta-2,2
      angle=(q-1)*180._f/(ntheta-1)          ! angle in deg 
      a0=fx_vars%nb*xint(q-1)*s11(q-1)*sin(angle*3.1415926353_f/180._f)
      c0=cos(angle*3.1415926353_f/180._f)

      angle=q*180._f/(ntheta-1)
      a1=fx_vars%nb*xint(q)*s11(q)*sin(angle*3.1415926353_f/180._f)
      c1=cos(angle*3.1415926353_f/180._f)

      angle=(q+1)*180._f/(ntheta-1)
      a2=fx_vars%nb*xint(q+1)*s11(q+1)*sin(angle*3.1415926353_f/180._f)
      c2=cos(angle*3.1415926353_f/180._f)

      total=total+2._f/6._f*3.1415926353_f/(ntheta-1)*(a0+4._f*a1+a2)
      totg=totg+2._f/6._f*3.1415926353_f/(ntheta-1)*(a0*c0+4._f*a1*c1+a2*c2)
    end do
    totg=totg/total

    ! *** ******************************************************************       
    ! *** angular integration of I(theta) according to                             
    ! *** Botet (1997) Eq.(26)  or  Rannou (1999) Eq.(1)                           
    ! ***      I(theta)  =  N 2pi/k^2 * S(theta) * [ phase integral ]              
    ! *** with                                                                     
    ! ***      S(theta)      =  s11_n(i)                                           
    ! ***      [ phase i. ]  =  xint(i)                                            
    ! *** Perfom integration using the following rule:                             
    ! ***  Integral_0^pi[ I(theta) sin(theta) d_theta ]                            
    ! ***                                                                          
    ! ***   = Sum_q=1^ntheta-1{ Integral_th_(i-1)^th_i [                           
    ! ***                                                                          
    ! ***       1/2(I(th_(i-1))+I(th_i)) * sin(th) d_th ] }                        
    ! ***                                                                          
    ! ***   = sin(delta_theta/2) * Sum_q=1^ntheta-1{                               
    ! ***                                                                          
    ! ***       ( I(th_(i-1)) + I(th_i) ) * sin(th_middle)  }                      
    ! ***                                                                          
    ! *** ****************************************************************** 

    !dthetad = 180._f / (ntheta-1)     ! angular interval in deg              
    dthetar =   pi   / (ntheta-1)     ! angular interval in rad              
    sumsca = 0._f
    do q=1,ntheta-1,1
      angler = (DBLE(q)-.5_f)*dthetar ! middle of interval in rad           
      weight = SIN(angler)           ! integration weight                  
      val    = s11_n(q-1)*xint(q-1) + s11_n(q)*xint(q)
      sumsca = sumsca + val*weight
    end do
    
    sumsca = sin(.5_f*dthetar) * sumsca  ! interval width factor            
    ! *** Scattering cross section                                                 
    sigsca = 2._f * pi / fx_vars%k**2._f * DBLE(fx_vars%nb) * sumsca
    ! Warning! sigabs is well computed using this approximation                    
    sigabs=fx_vars%nb*(sigmae-sigmas)
    ! sigext=sigabs+sigsca is better than the mean-field value
    ! previously defined. This is used hereafter. (P.Rannou)

    ! *** Radius of equivalent mass sphere
    rems = fx_vars%a * fx_vars%nb**(1._f/3._f)

    ! *** reference area in definition of efficiencies is the geometrical
    ! *** cross section of equivalent mass sphere
    gems = pi * rems**2._f

    ! *** Extinction and scattering efficiencies:                                  
    qsca = sigsca / gems
    qabs = sigabs / gems
    qext = qabs + qsca

    gfac = totg
    
  end subroutine fractal_meanfield

  !!
  !! Mie-scattering routine calling interface
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  subroutine cmie(lambda,xn,xk,rad,an,bn,nstop)
    
    ! Arguments
    real(kind=f), intent(in) :: lambda      !! wavelength (microns)
    real(kind=f), intent(in) :: xn          !! real index of refraction
    real(kind=f), intent(in) :: xk          !! imaginary index of refraction
    real(kind=f), intent(in) :: rad         !! monomer radius (meters)
    complex(kind=f), intent(out) :: an(50)  !! Mie wave coefficient an
    complex(kind=f), intent(out) :: bn(50)  !! Mie wave coefficient bn
    integer, intent(out) :: nstop           !! index of last mie-coefficent

    ! Local declarations
    integer, parameter :: nang = 451    ! number of angles
    complex(kind=f) :: refrel           ! complex index of refraction
    real(kind=f) :: theta(10000)        
    real(kind=f) :: x,dang

    refrel=cmplx(xn,xk,kind=f)
    x=2._f*3.14159265_f*rad/lambda      ! size parameter of monomer
    dang=1.570796327_f/real(nang-1,kind=f)

    call intmie(x,refrel,nang,an,bn,nstop)

    return
  end subroutine cmie

  !!
  !! Mie scattering calculations
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  SUBROUTINE intmie(x,refrel,nang,an,bn,nstop)

    ! Arguments
    real(kind=f), intent(in) :: x               !! size parameter of monomer 
    complex(kind=f), intent(in) :: refrel       !! complex index of refraction
    integer, intent(in) :: nang                 !! number of angles
    complex(kind=f), intent(out) :: an(nf)      !! Mie wave coefficient an
    complex(kind=f), intent(out) :: bn(nf)      !! Mie wave coefficient an
    integer, intent(out) :: nstop               !! index of last mie-coefficent

    ! Local declarations
    real(kind=f) :: amu(10000),pi(10000)
    real(kind=f) :: pi0(10000),pi1(10000)
    complex(kind=f) :: d(3000),y,xi,xi0,xi1
    complex(kind=f) ::  s1(2000),s2(2000)
    real(kind=f) psi0,psi1,psi,dn,dx
    integer :: nmx,nn,n,j
    real(kind=f) :: rn, xstop, dang, ymod, chi0, chi1, apsi0, apsi1, fn, chi, apsi

    dx=x            
    y=x*refrel

    xstop=x+4._f*x**.3333_f+2._f
    nstop=xstop
    ymod=abs(y)
    nmx=dmax1(xstop,ymod)+15
    dang=1.570796327_f/real(nang-1,kind=f)

    ! Initializations
    pi0(:) = 0._f
    pi1(:) = 0._f
    s1(:) = cmplx(0._f,0._f,kind=f)
    s2(:) = cmplx(0._f,0._f,kind=f)
    amu(:) = 0.0_f
    pi(:) = 0.0_f

    d(:) = cmplx(0._f,0._f,kind=f)
    nn=nmx-1

    do n=1,nn
      rn=nmx-n+1
      d(nmx-n)=(rn/y)-(1._f/(d(nmx-n+1)+rn/y))
    end do

    do j=1,nang
      pi0(j)=0._f       ! Legendre functions
      pi1(j)=1._f
    end do

    nn=2*nang-1

    do j=1,nn
      s1(j)=cmplx(0._f,0._f,kind=f)
      s2(j)=cmplx(0._f,0._f,kind=f)
    end do

    psi0=cos(dx)      ! Initialize Bessel functions
    psi1=sin(dx)
    chi0=-sin(x)
    chi1=cos(x)

    apsi0=psi0        
    apsi1=psi1

    xi0=cmplx(apsi0,-chi0,kind=f)
    xi1=cmplx(apsi1,-chi1,kind=f)

    n=1

    !   ************* iterate over index n ************* 
200 dn=n
    rn=n
    fn=(2._f*rn+1._f)/(rn*(rn+1._f))
 
    psi=(2._f*dn-1._f)*psi1/dx-psi0     ! calculate Bessel functions
    chi=(2._f*rn-1._f)*chi1/x-chi0      
    apsi=psi
    xi=cmplx(apsi,-chi,kind=f)
 
    an(n)=(d(n)/refrel+rn/x)*apsi-apsi1
    an(n)=an(n)/((d(n)/refrel+rn/x)*xi-xi1)
    bn(n)=(refrel*d(n)+rn/x)*apsi-apsi1
    bn(n)=bn(n)/((refrel*d(n)+rn/x)*xi-xi1)
 
    psi0=psi1
    psi1=psi
    apsi1=psi1           
 
    chi0=chi1
    chi1=chi
    xi1=cmplx(apsi1,-chi1,kind=f)

    n=n+1
    rn=n
 
    do 999 j=1,nang
      pi1(j)=((2._f*rn-1._f)/(rn-1._f))*amu(j)*pi(j)
      pi1(j)=pi1(j)-rn*pi0(j)/(rn-1._f)
999   pi0(j)=pi(j)

    if (n-1-nstop) 200,300,300
300 continue

    return
  END SUBROUTINE intmie

  !!
  !!
  !! CALLS:  FUNCTION dqag/dqdag/DADAPT_() Integration
  !!      FUNCTION fpl()           Integrand
  !!
  !! Integral in eq. 14, Botet et al. 1997
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  function funa(carma,fx_vars,nu,n,p,rc)
    type(carma_type), intent(in)    :: carma            !! the carma object
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated
    integer, intent(in)             :: n                !! indices
    integer, intent(in)             :: nu               !! indices
    integer, intent(in)             :: p                !! indices 
    integer, intent(inout)          :: rc               !! return code
    real(kind=f)                    :: funa             !! 

    ! Local declarations
    integer, parameter              :: maxsub=1000   
    real(kind=f)                    :: r,xa,xb,era,erl
    integer                         :: interv      
    integer                         :: ifail
    integer, parameter              :: lenw=4000                         ! .ge. 4*maxsub
    integer                         :: iwork(maxsub),neval,nsubin        ! nsubin=last
    real(kind=f)                    :: work(lenw)
    real(kind=f)                    :: bound, rres, rerr 

    ! Set return code assuming success.
    rc = RC_OK
    
    ! Initializations
    funa=0._f
    fx_vars%u1=n
    fx_vars%u2=1
    fx_vars%u3=nu
    fx_vars%u4=1
    fx_vars%u5=p
    fx_vars%u6=0
    xa=-1._f
    xb=1._f
    bound=0._f
    interv=1
    era=0._f
    erl=1.e-4_f
    rres=0._f
    rerr=0._f

    !======================================================================
    !--- Version using the QUADPACK - routine :
    !----------------------------------------------------------------------
    ifail = 0
    
    call dqag(fpl,fx_vars,xa,xb,era,erl,3,rres,rerr,neval,ifail,maxsub,lenw,nsubin,iwork,work)

    if (ifail.ne.0) then
      if (carma%f_do_print) then
         write(carma%f_LUNOPRT, *) "funa::ifail=",ifail, &
              " returned by dqag() during call of funa(",nu,",",n,",",p,")"
      end if
      rc = RC_ERROR
      return
    endif

    rres=rres-2._f          ! ceci est un artifice pour eviter que
                            ! la routine se plante quand la fonction
                            ! est paire (res=0.;err=1.d-3 impossible
                            ! a atteindre!! j'ai fpl'=fpl+1....d'ou
                            ! int(fpl)=int(fpl')-2. integr de -1 a 1!

    r = (2._f*p+1._f)/2._f
    funa = r * rres

    return
  END FUNCTION funa

  !!
  !! CALLS:  FUNCTION plgndr() Legendre-Functions
  !!
  !! Used in funa.  Integrand of eq. 14, Botet et al. 1997
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  FUNCTION fpl(x, fx_vars)

    ! Arguments    
    real(kind=f),intent(in) :: x                        !!
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated

    ! Local declarations
    real(kind=f) :: fpl
    integer :: m,n,mu,nu,p,pmu
    real(kind=f) :: c1,c2,c3

    c1=plgndr(fx_vars%u1,fx_vars%u2,x,fx_vars)
    c2=plgndr(fx_vars%u3,fx_vars%u4,x,fx_vars)
    c3=plgndr(fx_vars%u5,fx_vars%u6,x,fx_vars)

    fpl=c1*c2*c3+1._f        !this is a trick!

    return
  END FUNCTION fpl

  !!
  !! Adapted from FUNCTION plgndr() in: Press, Teukolsky, Vetterling, Flannery
  !! "Numerical Recipes in ???" (e.g. Num.Rec.in C, 2nd Ed., Cambridge Univ.Press, 1992, page 254)
  !!
  !! Calculate Legendre Polynomials, used in eq. 14 Botet et al. 1997
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  FUNCTION plgndr(l,m,x,fx_vars)

    ! Arguments
    integer, intent(in) :: l                         !! indices                          
    integer, intent(in) :: m                         !! indices
    real(kind=f), intent(in) :: x                    !! return result
    type(adgaquad_vars_type), intent(in) :: fx_vars  !! variables for functions being integrated

    ! Local declarations
    real(kind=f) :: plgndr
    integer ::lbl  
    real(kind=f) :: pll, pmm, somx2, pmmp1
    integer :: i, ll
    real(kind=f) :: fact1
    integer :: mstar

    mstar=m

    lbl=0
    plgndr=0._f

    if (mstar.lt.0)then
      mstar=-m
      lbl=1
    endif

    if (mstar.gt.l) then
      pll=0._f
      plgndr=0._f
      return          ! si m>l, Pl,m=0 !
    endif

    pmm=1._f

    if(mstar.gt.0) then
      somx2=sqrt((1._f-x)*(1._f+x))
      fact1=1._f
      do i=1,mstar
        pmm=+pmm*fact1*somx2    !cghmt - en + !!
        fact1=fact1+2._f
      end do
    endif

    if(l.eq.mstar) then
      plgndr=pmm
    else
      pmmp1=x*(2*mstar+1)*pmm

      if(l.eq.mstar+1) then
        plgndr=pmmp1
      else
        do ll=mstar+2,l
          pll=(x*(2*ll-1)*pmmp1-(ll+mstar-1)*pmm)/(ll-mstar)
          pmm=pmmp1
          pmmp1=pll
        end do
        plgndr=pll
      endif
    endif

    if (lbl.eq.1) then
      plgndr=(-1)**mstar*(fx_vars%fact(l-mstar)/fx_vars%fact(l+mstar))*plgndr
      mstar=-m         !restitution du parametre m!!!!!
    endif

    return
  END FUNCTION plgndr

  !! 
  !! replaces funb(nu,n,p) in original code, 
  !! saving n*n re-calculations of funa(nu,n,p).
  !!
  !! Calculates eq. 15, Botet et al. 1997
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  FUNCTION funb_n(nu,n,p,funca)
    
    ! Arguments
    integer, intent(in) :: nu                            !! indices
    integer, intent(in) :: n                             !! indices
    integer, intent(in) :: p                             !! indices  
    real(kind=f), intent(in) :: funca(nmi,nmi,0:n2m)     !! return result

    ! Local Declarations
    real(kind=f)  :: funb_n
    integer :: i, l, j
    real(kind=f) :: var
  
    funb_n = 0._f
    i = int((p*1._f-1._f-abs(n*1._f-nu*1._f))*1._f/2._f)
    !print*,nu,n,p,i
 
    do l=0,i
      j = p-2*l-1 

      ! omit j = -1 (when nu=n and p=l=i=0)
      IF (j .GE. 0) THEN

        var = funca(nu,n,j)    ! in main, a(nu,n,p) was stored in
                               ! funca(nu,n;p) 
        funb_n = funb_n + var
      ENDIF

    end do
    funb_n = (2._f*p+1._f) * funb_n
    return
  END FUNCTION funb_n

  !! Replaces funs(pp,k) in original code
  !!
  !! CALLS:
  !!    FUNCTION dqagi/dq2agi/DADAPT_() Integration
  !!    FUNCTION xfreal_n()      Integrand
  !!    FUNCTION xfimag_n()      Integrand
  !!
  !! Calculates eq. 16 , Botet et al. 1997
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version Mar 2013
  function funs_n(carma,fx_vars,pp,rc)

    ! Arguments
    type(carma_type), intent(in)    :: carma            !! the carma object
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated
    integer, intent(in)             :: pp               !! indices
    integer, intent(inout)          :: rc               !! return code

    ! Local Declarations
    integer, parameter :: maxsub=50000
    complex(kind=f) :: rcomplex,funs_n
    real(kind=f) :: rres,ires,rerr,ierr,afun
    real(kind=f) :: xa,xb
    integer :: ifail
    integer, parameter :: lenw=200000            ! .ge. 4*maxsub
    integer :: iwork(maxsub),neval,nsubin         ! nsubin=last
    real(kind=f) ::  work(lenw)
    real(kind=f) :: rg, bound, errabs, errrel
    integer :: interv
    
    rc = RC_OK

    rg=fx_vars%alpha*fx_vars%nb**(1._f/fx_vars%df)*fx_vars%a

    afun=(2._f*3.1415926_f)/(fx_vars%k**3._f)


    fx_vars%pbes=pp
    fx_vars%kbes=fx_vars%k

    bound=0._f
    interv=1

    errabs=0._f
    errrel=1.e-3_f

    rres=0._f
    !trres=0._f
    rerr=0._f
    !trerr=0._f
    xa=0._f
    xb=5._f*fx_vars%k*rg

    !======================================================================
    !--- Version using the QUADPACK - routine :
    !----------------------------------------------------------------------
    ifail = 0
    CALL dqagi(xfreal_n,fx_vars,bound,interv,errabs,errrel,rres,rerr,neval,ifail,maxsub,lenw,nsubin,iwork,work)
    if (ifail.ne.0) then
      if (carma%f_do_print) then
         write(carma%f_LUNOPRT, *) "funs_n::ifail=",ifail, &
              " returned by dqag() during call of funs(",pp, &
              ") while integrating xfreal_n()"
      end if
      rc = RC_ERROR
      return
    endif

    bound=0._f
    interv=1

    ires=0._f
    ierr=0._f
    xa=0._f
    xb=5._f*fx_vars%k*rg

    !======================================================================
    !--- Version using the QUADPACK - routine :
    !----------------------------------------------------------------------
    ifail = 0
    CALL dqagi(xfimag_n,fx_vars,bound,interv,errabs,errrel,ires,ierr,neval,ifail,maxsub,lenw,nsubin,iwork,work)
    if(ifail.ne.0) then
      if (carma%f_do_print) then
         write(carma%f_LUNOPRT, *) "funs_n::ifail=",ifail, &
              " returned by dqagi() during call of funs(",pp, &
              ") while integrating xfimag_n()"
      end if
      rc = RC_ERROR
      return
    endif

    rcomplex = cmplx(1._f,0._f,kind=f)*rres + cmplx(0._f,1._f,kind=f)*ires

    funs_n = afun * rcomplex

    continue
    return
  END FUNCTION funs_n

  !!
  !! replaces xfreel(xx) in original code
  !! CALLS:  FUNCTION BESSELJY() Spherical Bessel functions
  !!         FUNCTION phi()             Probability distrib.
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version Mar 2013
  FUNCTION xfreal_n(xx, fx_vars)

    ! Arguments
    real(kind=f), intent(in) :: xx
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated

    ! Local Declarations
    complex(kind=f) :: z,xj(0:nf),xjp(0:nf),xy(0:nf),xyp(0:nf)
    complex(kind=f)  :: jsol,ysol,hsol,hpsol
    real(kind=f) :: x,r,xfreal_n
    integer :: ifail,p,pc
    real(kind=f) :: rg

    ifail = 0
    rg = fx_vars%alpha*fx_vars%nb**(1._f/fx_vars%df)*fx_vars%a
    x = xx
    if (x.GT.3000._f)  x=3000._f
    z = x*cmplx(1._f,0._f,kind=f)

    pc = fx_vars%pbes
    if( fx_vars%pbes .eq. 0 ) pc = fx_vars%pbes + 1

    CALL BESSELJY(z,pc,xj,xjp,xy,xyp,ifail)

    r=x/fx_vars%kbes

    xfreal_n = real( z*z*xj(fx_vars%pbes)*xj(fx_vars%pbes) * phi(r,fx_vars) )

    return
  END FUNCTION xfreal_n

  !!
  !! replaces xfima(xx) in original code
  !! CALLS:  FUNCTION BESSELJY() Spherical Bessel functions
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version Mar 2013
  FUNCTION xfimag_n(xx, fx_vars)
    
    ! Arguments
    real(kind=f), intent(in) :: xx
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! variables for functions being integrated

    ! Local Declarations
    complex(kind=f) :: z,xj(0:nf),xjp(0:nf),xy(0:nf),xyp(0:nf)
    real(kind=f) :: x,r,xfimag_n
    integer :: ifail,p,pc
    real(kind=f) :: rg

    rg = fx_vars%alpha*fx_vars%nb**(1._f/fx_vars%df)*fx_vars%a
    x = xx
    if (x.gt.3000._f)  x=3000._f
    ifail = 0
    z = x*cmplx(1._f,0._f,kind=f)

    pc = fx_vars%pbes
    if( fx_vars%pbes .eq. 0 ) pc = fx_vars%pbes + 1

    CALL BESSELJY(z,pc,xj,xjp,xy,xyp,ifail)

    r=x/fx_vars%kbes

    xfimag_n = real( z*z*xj(fx_vars%pbes)*xy(fx_vars%pbes) * phi(r,fx_vars) )

    return
  END FUNCTION xfimag_n


  !! Spherical Bessel functions j_n(z) and y_n(z) of complex
  !! argument to desired accuracy,
  !! and their derivatives, up to a maximal order n=LMAX. 
  !! j_n(z) = SQRT(pi/2 / z) * J_(n + 1/2)(z) 
  !! y_n(z) = SQRT(pi/2 / z) * Y_(n + 1/2)(z) 
  !! Adapted from:
  !!         I.J.Thompson, A.R.Barnett
  !!        "Modified Bessel Funkctions I_v(z) and K_v(z)
  !!         of Real Order and Complex Argument, to Selected
  !!         Accuracy"
  !!         COMP.PHYS.COMMUN. 47 (1987) 245-57
  !!         (Source code printed on page 249)
  !! ******************************************************************
  !! INPUTS:
  !!   X  argument z, dble cmplx
  !!            z in the upper half plane, Im(z) > -3
  !!   LMAX       largest desired order of Bessel functions, int
  !!            j_n,y_n,j_n',y_n' are calculated for n=0 to n=LMAX
  !!            Dimension of arrays xj,xjp,xy,xyp at least (0:LMAX)
  !!   XJ(M)    Spher. Bessel function j_m(z), dble cmplx
  !!   XJP(M)   Derivative of Spher. Bessel function d/dz [ j_m(z) ],
  !!            dble cmplx
  !!   XY(M)    Spher. Bessel function y_m(z), dble cmplx
  !!   XYP(M)   Derivative of Spher. Bessel function d/dz [ y_m(z) ],
  !!            dble cmplx
  !!   IFAIL    error flag, int
  !!             =   0  if all results are satisfactory
  !!             =  -1  for arguments out of range
  !!             = > 0  for results ok up to and including the
  !!                    function of order LMAX-IFAIL
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version Mar 2013
  SUBROUTINE BESSELJY (X, LMAX, XJ, XJP, XY, XYP, IFAIL)
    
    ! Arguments
    complex(kind=f), intent(in) :: X
    integer, intent(in) :: LMAX
    complex(kind=f), intent(out) :: XJ(0:LMAX)
    complex(kind=f), intent(out) :: XJP(0:LMAX)
    complex(kind=f), intent(out) :: XY(0:LMAX)
    complex(kind=f), intent(out) :: XYP(0:LMAX)
    integer, intent(out) :: IFAIL   

    ! Local Declarations
    INTEGER, PARAMETER :: LIMIT = 20000
    REAL(kind=f),parameter :: ZERO  = 0._f
    REAL(kind=f),parameter :: ONE   = 1._f
    REAL(kind=f),parameter :: ACCUR = 1e-12_f
    REAL(kind=f),parameter :: TM30  = 1e-30_f
    COMPLEX(kind=f), parameter :: CI = (0._f, 1._f)
    complex(kind=f) :: XI, W, PL, B, D, FF, DEL, C, XJ0, XH1, XH1P, XTEMP
    integer :: L

    IF (ABS(X).LT.ACCUR .OR. AIMAG(X) .LT. -3.d0) THEN
      IFAIL=-1
      GOTO 5
    END IF

    ! *** Lentz - Algorithmus (?) :
    XI = ONE/X
    W = XI + XI
    PL = LMAX * XI
    FF = PL + XI
    B = FF + FF + XI
    D = ZERO
    C = FF
    DO 1 L=1,LIMIT
      D = B - D
      C = B - ONE/C
      IF(ABS(D).LT. TM30) D = TM30
      IF(ABS(C).LT. TM30) C = TM30
      D = ONE / D
      DEL = D * C
      FF = FF * DEL
      B = B + W
1     IF(ABS(DEL-ONE).LT.ACCUR) GOTO 2
    IFAIL = -2
    GOTO 5

2   XJ(LMAX) = TM30
    XJP(LMAX) = FF * XJ(LMAX)

    ! *** Abwaertsrekursion
    DO 3 L = LMAX-1,0,-1
      XJ(L) = PL * XJ(L+1) + XJP(L+1)
      XJP(L) = PL * XJ(L) - XJ(L+1)
3     PL = PL - XI
   
    ! *** Calculate the l=0 Besselfunktionen
    XJ0 = XI * SIN(X)
    XY(0) = - XI * COS(X)
    XH1 = EXP(CI * X) * XI * (-CI)
    XH1P = XH1 * (CI - XI)
    B = XH1P

    ! *** Rescale XJ, XJP, converting to spherical Bessels
    ! *** Recur   XH1,XH1P   as sperical Bessels
    W = ONE / XJ(0)
    PL = XI
    DO 4 L = 0,LMAX
      XJ(L) = XJ0 * (W*XJ(L))
      XJP(L) = XJ0 * (W*XJP(L)) - XI * XJ(L)
      IF (L.EQ.0) GOTO 4
      XTEMP = XH1
      XH1 = (PL-XI) * XTEMP - XH1P
      PL = PL + XI
      XH1P = - PL * XH1 + XTEMP
      XY(L) = CI * (XJ(L) - XH1)    ! y_n = i * ( j_n - h^1_n )
      XYP(L) = CI * (XJP(L) - XH1P) ! und dito fuer Ableitungen
4   CONTINUE
    XYP(0) = CI * (XJP(0) - B)
    RETURN

5   WRITE(*,10) IFAIL
10  FORMAT( 'ERROR in SUBR BESSELJY() : IFAIL = ', I4)
    RETURN
  END SUBROUTINE BESSELJY

  !! 
  !! Angular function pi_l( x=cos(theta) )
  !! e.g. Bohren,Huffman (1983) 
  !!      pp.94 ff Eq.(4.46)-(4.49)
  !!      p.112
  !! CALLS:  FUNCTION plgndr() Legendre-Functions
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version Mar 2013
  FUNCTION fpi(l,x,fx_vars)
  
    ! Arguments
    integer, intent(in) :: l
    real(kind=f), intent(in) :: x
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated
    
    ! Local declarations
    real(kind=f) :: fpi
    real(kind=f) :: y 
    real(kind=f) :: flag

    y=x
    if (x.eq.1._f) y=1._f-1.e-6_f
          ! alternatively, one could use Bohren,Huffman
          ! p.112:   pi_n(1)=tau_n(1)= 1/2 * n * (n+1) !!!
    flag=plgndr(l,1,y,fx_vars)
    fpi=(1._f-y**2._f)**(-0.5_f)*flag
    return
  END FUNCTION fpi

  !!
  !! Angular function tau_l( x=cos(theta) )
  !! e.g. Bohren,Huffman (1983) 
  !!      pp.94 ff Eq.(4.46)-(4.49)
  !!      p.112
  !! CALLS:  FUNCTION plgndr() Legendre-Functions
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  FUNCTION tau(l,x,fx_vars)

    ! Arguments
    integer, intent(in) :: l
    real(kind=f), intent(in) :: x
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated

    ! Local Declarations
    real(kind=f) :: fp
    real(kind=f) :: tau
    real(kind=f) :: flag
    real(kind=f) :: y

    y=x
    if (x.eq.1._f) y=1._f-1.e-6_f
        ! alternatively, one could use Bohren,Huffman
  ! p.112:   pi_n(1)=tau_n(1)= 1/2 * n * (n+1) !!!
    flag=plgndr(l,0,y,fx_vars)
    fp=fpi(l,y,fx_vars)
    tau=-y*fp+l*(l*1._f+1._f)*flag
    return
  END FUNCTION tau

  !!
  !!
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  function fp1(u, fx_vars)

    ! Arguments
    real(kind=f), intent(in)                  :: u            !!
    type(adgaquad_vars_type), intent(inout)   :: fx_vars      !! varaibles for functions being integrated
    real(kind=f)                              :: fp1          !! returns

    ! Local Declarations
    real(kind=f) :: krg,s1,s2,s3,rg
   
    rg=fx_vars%alpha*fx_vars%a*fx_vars%nb**(1._f/fx_vars%df)
    krg=fx_vars%k*rg
    s1=sin(2._f*krg*fx_vars%zed*u)
    s2=u**(fx_vars%df-2._f)
    s3=fco(u, fx_vars)
    fp1=s1*s2*s3
    
    return
  END FUNCTION fp1

  !!
  !! CALLS:  FUNCTION dqagi/dqdagi/DADAPT_() Integration
  !!         FUNCTION fdval()     Integrand
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  FUNCTION anorm(carma, fx_vars, rc)

    ! arguments
    type(carma_type), intent(in)    :: carma         !! the carma object
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated
    integer, intent(inout)          :: rc            !! return code

    ! Local Declarations
    real(kind=f) :: anorm
    integer :: interv
    integer, parameter :: maxsub=50000
    integer :: ifail
    integer, parameter :: lenw=200000                   ! .ge. 4*maxsub
    integer :: iwork(maxsub),neval,nsubin        ! nsubin=last
    real(kind=f) :: work(lenw)
    real(kind=f) :: bound,errrel,errabs,b,db,c
    
    rc = RC_OK

    bound=0._f
    interv=1
    errrel=1.e-3_f
    errabs=0._f
    b=0._f
    db=0._f

    !======================================================================
    !--- Version using the QUADPACK - routine :
    !----------------------------------------------------------------------
    ifail = 0
    CALL dqagi(fdval,fx_vars,bound,interv,errabs,errrel,b,db,neval,ifail,maxsub,lenw,nsubin,iwork,work)
    if(ifail.ne.0) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "anorm::ifail=",ifail," returned by dqagi() during call of anorm"
      rc = RC_ERROR
      return
    endif
    
    c=0.5_f
    anorm=b*4._f*3.1415926_f
    return
  END FUNCTION anorm

  !!
  !! Probability distribution of monomer location within cluster
  !! CALLS:  FUNCTION fdval()
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  FUNCTION phi(x,fx_vars)   

    ! Arguments
    real(kind=f), intent(in) :: x
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated

    ! Local Declarations
    real(kind=f) :: fval,pref,phi, rg, z

    rg=fx_vars%alpha*fx_vars%nb**(1._f/fx_vars%df)*fx_vars%a
    z=x/rg
    pref=(x/rg)**(fx_vars%df-3._f)/(fx_vars%coeff*rg**3._f)
    fval=z**(1._f-fx_vars%df)*fdval(z, fx_vars)
    phi=pref*fval
    continue
    return
  END FUNCTION phi

  !!
  !! Probability distribution of monomer location within cluster
  !! CALLS:  FUNCTION fdval()
  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  FUNCTION fco(z, fx_vars)
    
    ! Arguments
    real(kind=f), intent(in) :: z
    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated

    ! Local Declarations
    real(kind=f) :: fco
    real(kind=f) :: fval

    fval=z**(1._f-fx_vars%df)*fdval(z, fx_vars)
    fco=fval
    continue
    return
  END FUNCTION fco

  !!
  !! @author P. Rannou, R. Botet, Eric Wolf
  !! @version March 2013
  FUNCTION fdval(x, fx_vars)

    type(adgaquad_vars_type), intent(inout) :: fx_vars  !! varaibles for functions being integrated
   
    ! Arguments
    real(kind=f), intent(in) :: x
   
    ! Local Declarations
    real(kind=f) :: fdval

    fdval=x**(fx_vars%df-1._f)*exp(-x**fx_vars%df/2._f)
    return
  END FUNCTION fdval

end module


