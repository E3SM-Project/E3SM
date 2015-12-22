c $Revision: 23 $, $Date: 2011-03-31 07:41:37 -0600 (Thu, 31 Mar 2011) $
c $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/quickbeam/predict_psd07.f $
      program predict_psd07

c 14 Feb 2007

c This procedure uses the second moment and in-cloud temperature
c to predict the ice particle size distribution. 
c The moment relations use the 2nd moment as a reference moment.


c For more details please read 'Snow size distribution parameterization for midlatitude and tropical ice clouds' by Field et al JAS2007

c !!!!!!!!!!!variables!!!!!!!!!!!!!!!!!!!!
c Input
c tc in C
c iwc in g m^-3
c x,y prefactor and exponent in mass-size relation m=xD^y
c for the UM    x=0.069, y=2.0    (SI units)
c regime: 'T' for tropical, 'M' for midlatitude

c Output
c D particle size in microns
c dN_dD particle size distribution: dN/dD m^-4
c L23 ratio of M3/M2 - mean mass weighted size if y=2


       real D(500),dN_dD(500),phi(500)
       real x,y,tc,iwc,iwcsi,My,m2,n,m3,xx
       character regime


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c    read in input file
       open(unit=2,file="input.txt")
       read(2,*) iwc,tc,x,y,regime
       close(2)


c incloud temperature check

      if (tc.gt.0) then
         print*,'in-cloud temperature too warm (>0C)'
         goto 1000
      endif  
      if (tc.lt.-70) then 
         print*,'Warning: in-cloud temp colder than in original data'
      endif  

c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



c Calculate 2nd moment
      iwcsi=iwc/1e3  !kg m^-3

      My=iwcsi/x  !yth moment

c use moment relations to get 2nd moment required for predicting psd

      m2=My 
c default for UM y=2 

c if y ne 2 then need to find second moment before proceeding
      if(y.ne.2.0) then 
         m2=-9999
         call predict_mom07(m2,tc,y,My)
      endif
c !!M2 range check (1e-5 3e-2 m^-1)!!!!!!!!!!!!!!!!!!!!!

      if (m2.lt.0.0) then
         print*,'Negative M2'
         goto 1000
      endif
      if (m2.lt.1e-5.or.m2.gt.0.1) then
        print*,'Warning: M2 outside of range in original data'
      endif 
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c Use 2nd and 3rd moment to predict psd;;;;;;;;;;;;;;;;;;;;


       n=3
       call predict_mom07(m2,tc,n,m3)
        

c !!!define universal psd!!!!
      do 10 i=1,500
        xx=i/40. !dimensionless size
        if (regime.eq.'T') then 
               phi(i)=152.0*exp(-12.4*xx)+3.28*xx**(-0.78)*exp(-1.94*xx)
        endif

        if (regime.eq.'M') then 
               phi(i)=141.0*exp(-16.8*xx)+102.0*xx**(2.07)*exp(-4.82*xx)
        endif
 	
        D(i)=1e6*xx*(m3/m2)
        dN_dD(i)=phi(i)*m2**4/m3**3
 10   continue




c !write output file

      open(unit=2,file="output.txt")
      write(2,*) 'Input values'
      write(2,*) 'IWC [g/m3]=',iwc
      write(2,*) 'T [C] =',tc
      write(2,*) 'x , y [SI]=', x, y
      write(2,*),'Regime =',regime
      write(2,*) ' '
      write(2,*) 'Output values'
      write(2,*) 'L32 [microns]=',1e6*m3/m2
      write(2,*) 'D (microns) dN/dD (m^-4)'
      do 20 i=1,500
        write(2,*),D(i),dN_dD(i)
 20   continue
      close(2)	       
 

 1000   continue

        stop 
        end


        subroutine predict_mom07(m2,tc,n,m)

c subroutine to predict nth moment m given m2, tc, n
c if m2=-9999 then the routine will predict m2 given m,n,tc
 
        real a1,a2,a3,b1,b2,b3,c1,c2,c3
        real m2,tc,n,m,a_,b_,c_,A,B,C

       a1=      13.6078
       a2=     -7.76062
       a3=     0.478694
       b1=   -0.0360722
       b2=    0.0150830
       b3=   0.00149453
       c1=     0.806856
       c2=   0.00581022
       c3=    0.0456723


       a_=a1+a2*n+a3*n**2.0
       b_=b1+b2*n+b3*n**2.0
       c_=c1+c2*n+c3*n**2.0

       A=exp(a_)
       B=b_
       C=c_
       
c predict m from m2 and tc
        if(m2.ne.-9999) then 
           m=A*exp(B*tc)*m2**C
        endif
c get m2 if mass-dimension relationship not proportional to D**2
        if(m2.eq.-9999) then 
           m2=(m/(A*exp(B*tc)))**(1.0/C)    
        endif
        

       return
       end






