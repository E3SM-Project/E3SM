! $Revision: 88 $, $Date: 2013-11-13 07:08:38 -0700 (Wed, 13 Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/predict_mom07.F90 $
        subroutine predict_mom07(m2,tc,n,m)
        
        ! subroutine to predict nth moment m given m2, tc, n
        ! if m2=-9999 then the routine will predict m2 given m,n,tc
            implicit none
        
            real*8 :: a1,a2,a3,b1,b2,b3,c1,c2,c3
            real*8 :: m2,tc,n,m,a_,b_,c_,A,B,C,n2
        
            a1=      13.6078
            a2=     -7.76062
            a3=     0.478694
            b1=   -0.0360722
            b2=    0.0150830
            b3=   0.00149453
            c1=     0.806856
            c2=   0.00581022
            c3=    0.0456723
        
            n2=n*n
            a_=a1+a2*n+a3*n2
            b_=b1+b2*n+b3*n2
            c_=c1+c2*n+c3*n2
        
            A=exp(a_)
            B=b_
            C=c_
            
        ! predict m from m2 and tc
                if(m2.ne.-9999) then 
                m=A*exp(B*tc)*m2**C
                endif
        ! get m2 if mass-dimension relationship not proportional to D**2
                if(m2.eq.-9999) then 
                m2=(m/(A*exp(B*tc)))**(1.0/C)    
                endif
        
                return
        end subroutine predict_mom07
