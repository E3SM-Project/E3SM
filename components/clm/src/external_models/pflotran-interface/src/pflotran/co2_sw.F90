module co2_sw_module

  ! module contains only 1 interface with other part, read as:

  ! co2_sw_prop(p,t,rho,dddt,dddp,fg,dfgdp,dfgdt,
  !                 eng,ent,dhdt,dhdp,visc,dvdt,dvdp)

  use PFLOTRAN_Constants_module

      implicit none

#include "petsc/finclude/petscsys.h"

      private
      save
      PetscReal,  public :: co2_sw_c_t0_tab = -5.d0+273.15D0, co2_sw_c_p0_tab = 0.5d0
      PetscReal,  public :: co2_sw_c_t1_tab = 5.d2+273.15D0, co2_sw_c_p1_tab = 250.d0
      PetscReal,  public :: co2_sw_f_t0_tab = -5.d0+273.15D0, co2_sw_f_p0_tab = 0.5d0
  !   PetscReal,  public :: co2_sw_f_t1_tab = 1.d2+273.15D0, co2_sw_f_p1_tab = 15.d0
      PetscReal,  public :: co2_sw_f_t1_tab = 1.d2+273.15D0, co2_sw_f_p1_tab = 30.d0

      PetscInt :: ntab_t_c = 50, ntab_p_c = 100
     ! PetscInt :: ntab_t_f = 105, ntab_p_f = 145
      PetscInt :: ntab_t_f = 105, ntab_p_f = 295
       
      PetscReal :: dt_tab_c = 10.d0, dp_tab_c = 2.5d0
      PetscReal :: dt_tab_f = 1.d0, dp_tab_f = 0.1d0
  
      PetscReal ,private, allocatable :: co2_prop_sw_c(:,:,:)
      PetscReal ,private, allocatable :: co2_prop_sw_f(:,:,:)
      PetscInt var_index(4)
      logical ifinetable

      PetscReal, parameter, private ::  denc = 467.6d0, tc = 304.1282d0,&
                                     rg = 0.1889241d0, pc = 7.3773d0

      public initialize_sw_interp, co2_sw_interp, init_span_wagner  
 contains

     
! ************************************************************************** !

subroutine init_span_wagner(option)
  ! 
  ! init_span_wagner
  ! 
  ! Author: Chuan Lu
  ! Date: 5/13/08
  ! 
  use Option_module
  use co2_span_wagner_module
  use co2_span_wagner_spline_module

  implicit none
  type(option_type) :: option
  PetscMPIInt :: myrank

  if (option%co2eos == EOS_SPAN_WAGNER) then
    select case(option%itable)
       case(0,1,2)
         call initialize_span_wagner(option%itable, &
                                     option%myrank, &
                                     option)
       case(4,5)
         myrank = option%myrank
         call initialize_span_wagner(ZERO_INTEGER,myrank, &
                                     option)
         call initialize_sw_interp(option%itable,myrank)
       case(3)
         call sw_spline_read
       case default
         print *, 'Wrong table option : STOP'
         stop
    end select
  endif
end subroutine init_span_wagner

! ************************************************************************** !

subroutine initialize_sw_interp(itable,myrank)
  ! 
  ! prepare data table for interpolation
  ! 
      use co2_span_wagner_module, only: co2_span_wagner, vappr
        
        implicit none
        PetscInt itable
        PetscMPIInt myrank
        
        PetscInt i, j
        PetscInt :: iflag = 1
        
        PetscReal p,t, ts,tmp2, tmp, dum1, dum2
!       PetscReal dtmp,dddt,dddp
        PetscReal dtemp,dpres
        PetscReal rhodp,rhodt,fgdp,fgdt,engdp,engdt,entdp,entdt,vdp,vdt

        character*3 :: q
        character*1 :: tab
      
      
      
      

        tab = char(9)
         q = '","'
        
        allocate( co2_prop_sw_c(0:ntab_p_c,0:ntab_t_c,1:15)) 
        allocate( co2_prop_sw_f(0:ntab_p_f,0:ntab_t_f,1:15))
      
        var_index(1) = 3; var_index(2) = 6; var_index(3) = 10; var_index(4) = 13;  
        
        
        select case(itable)
        case(5)
! prepare the coarse table  
         tmp2=0D0
         do j = 0, ntab_t_c
            tmp=tmp2
            t = co2_sw_c_t0_tab + dt_tab_c * real(j)
           do i = 0, ntab_p_c
             p = co2_sw_c_p0_tab + dp_tab_c * real(i)
             co2_prop_sw_c(i, j, 1) = p
             co2_prop_sw_c(i, j, 2) = t
             co2_prop_sw_c(i,j,3) = tmp
             call co2_span_wagner(p,t,co2_prop_sw_c(i,j,3),co2_prop_sw_c(i,j,4),&
                co2_prop_sw_c(i,j,5),co2_prop_sw_c(i,j,6),co2_prop_sw_c(i,j,7),&
                co2_prop_sw_c(i,j,8),co2_prop_sw_c(i,j,9),co2_prop_sw_c(i,j,10),&
                co2_prop_sw_c(i,j,11),co2_prop_sw_c(i,j,12),co2_prop_sw_c(i,j,13),&
                co2_prop_sw_c(i,j,14),co2_prop_sw_c(i,j,15),iflag)
             
               ts = 1000.D0
               if (p .le.pc .and. t .le.tc) then
                 call vappr (ts,p,dum1,dum2,TWELVE_INTEGER)
                endif

             
             dpres = 1.e-3
             if (t >ts) dpres = - dpres
              call co2_span_wagner(p+dpres,t,rhodp,co2_prop_sw_c(i,j,4),&
                      co2_prop_sw_c(i,j,5),fgdp,co2_prop_sw_c(i,j,7),&
                      co2_prop_sw_c(i,j,8),engdp,entdp,&
                      co2_prop_sw_c(i,j,11),co2_prop_sw_c(i,j,12),vdp,&
                      co2_prop_sw_c(i,j,14),co2_prop_sw_c(i,j,15),iflag)

        
               dtemp = 1.e-6
               if (t  < ts) dtemp = -dtemp
               call co2_span_wagner(p,t+dtemp,rhodt,co2_prop_sw_c(i,j,4),&
                 co2_prop_sw_c(i,j,5),fgdt,co2_prop_sw_c(i,j,7),&
                 co2_prop_sw_c(i,j,8),engdt,entdt,&
                 co2_prop_sw_c(i,j,11),co2_prop_sw_c(i,j,12),vdt,&
                 co2_prop_sw_c(i,j,14),co2_prop_sw_c(i,j,15),iflag)

               tmp = co2_prop_sw_c(i,j,3)
               if (j==0) tmp2 = co2_prop_sw_c(i,j,3)
      
                
                   !dddt
          co2_prop_sw_c(i,j,4) = (rhodt-co2_prop_sw_c(i,j,3))/dtemp
          
          !dfgdt
          co2_prop_sw_c(i,j,8) = (fgdt-co2_prop_sw_c(i,j,6))/dtemp
          
          !dhdt
          co2_prop_sw_c(i,j,11) = (entdt-co2_prop_sw_c(i,j,10))/dtemp
          
          !dvdt
          co2_prop_sw_c(i,j,14) = (vdt-co2_prop_sw_c(i,j,13))/dtemp
!       endif
!       if (i>0) then
          !dddp
          co2_prop_sw_c(i,j,5) = (rhodp-co2_prop_sw_c(i,j,3))/dpres
          
          !dfgdp
          co2_prop_sw_c(i,j,7) = (fgdp-co2_prop_sw_c(i,j,6))/dpres
          
          !dhdp
          co2_prop_sw_c(i,j,12) = (entdp-co2_prop_sw_c(i,j,10))/dpres
          
          !dvdp
          co2_prop_sw_c(i,j,15) = (vdp-co2_prop_sw_c(i,j,13))/dpres
          
            enddo
!            print *,'co2_sw: ', i,j, p,t, co2_prop_sw_c(i,j,3)
         enddo    


 ! Prepare fine table
         tmp2=0D0
         do j = 0, ntab_t_f
            tmp=tmp2
            t = co2_sw_f_t0_tab + dt_tab_f * real(j)
           do i = 0, ntab_p_f
             p = co2_sw_f_p0_tab + dp_tab_f * real(i)
             co2_prop_sw_f(i, j, 1) = p
             co2_prop_sw_f(i, j, 2) = t
             co2_prop_sw_f(i,j,3) = tmp
             call co2_span_wagner(p,t,co2_prop_sw_f(i,j,3),co2_prop_sw_f(i,j,4),&
                co2_prop_sw_f(i,j,5),co2_prop_sw_f(i,j,6),co2_prop_sw_f(i,j,7),&
                co2_prop_sw_f(i,j,8),co2_prop_sw_f(i,j,9),co2_prop_sw_f(i,j,10),&
                co2_prop_sw_f(i,j,11),co2_prop_sw_f(i,j,12),co2_prop_sw_f(i,j,13),&
                co2_prop_sw_f(i,j,14),co2_prop_sw_f(i,j,15),iflag)
             
               ts = 1000.D0
               if (p .le.pc .and. t .le.tc) then
                 call vappr (ts,p,dum1,dum2,TWELVE_INTEGER)
                endif

             
             dpres = 1.e-3
             if (t >ts) dpres = - dpres
              call co2_span_wagner(p+dpres,t,rhodp,co2_prop_sw_f(i,j,4),&
                      co2_prop_sw_f(i,j,5),fgdp,co2_prop_sw_f(i,j,7),&
                      co2_prop_sw_f(i,j,8),engdp,entdp,&
                      co2_prop_sw_f(i,j,11),co2_prop_sw_f(i,j,12),vdp,&
                      co2_prop_sw_f(i,j,14),co2_prop_sw_f(i,j,15),iflag)

        
               dtemp = 1.e-6
               if (t < ts) dtemp = -dtemp
               call co2_span_wagner(p,t+dtemp,rhodt,co2_prop_sw_f(i,j,4),&
                 co2_prop_sw_f(i,j,5),fgdt,co2_prop_sw_f(i,j,7),&
                 co2_prop_sw_f(i,j,8),engdt,entdt,&
                 co2_prop_sw_f(i,j,11),co2_prop_sw_f(i,j,12),vdt,&
                 co2_prop_sw_f(i,j,14),co2_prop_sw_f(i,j,15),iflag)

               tmp = co2_prop_sw_f(i,j,3)
               if (j==0) tmp2 = co2_prop_sw_f(i,j,3)
      
                
                   !dddt
          co2_prop_sw_f(i,j,4) = (rhodt-co2_prop_sw_f(i,j,3))/dtemp
          
          !dfgdt
          co2_prop_sw_f(i,j,8) = (fgdt-co2_prop_sw_f(i,j,6))/dtemp
          
          !dhdt
          co2_prop_sw_f(i,j,11) = (entdt-co2_prop_sw_f(i,j,10))/dtemp
          
          !dvdt
          co2_prop_sw_f(i,j,14) = (vdt-co2_prop_sw_f(i,j,13))/dtemp
!       endif
!       if (i>0) then
          !dddp
          co2_prop_sw_f(i,j,5) = (rhodp-co2_prop_sw_f(i,j,3))/dpres
          
          !dfgdp
          co2_prop_sw_f(i,j,7) = (fgdp-co2_prop_sw_f(i,j,6))/dpres
          
          !dhdp
          co2_prop_sw_f(i,j,12) = (entdp-co2_prop_sw_f(i,j,10))/dpres
          
          !dvdp
          co2_prop_sw_f(i,j,15) = (vdp-co2_prop_sw_f(i,j,13))/dpres
           
           enddo
!          print *,' co2_sw: ', i,j, p,t,co2_prop_sw_f(i,j,3)
         enddo    


      
      if (myrank == 0) then
        print *,'Writing Table lookup file: Coarse...'
            if (myrank==0) print *,'--> open co2data_c.dat'
           open(unit=122,file='co2data_c.dat',status='unknown')
             write(122,'(''TITLE= "'',''co2data_c.dat'',''"'')')
             write(122,'(''VARIABLES= "'',a6,100(a3,a6))') &
          'p',q,'T',q,'d',q,'dddT',q,'dddp',q,'fg',q,'dfgdp',q,'dfgdT',q, &
          'u',q,'h',q,'dhdT',q,'dhdp',q,'vis',q,'dvdT',q,'dvdp','"'
             write(122,'(''ZONE T= "'',''",'','' I='',i4,'' , J='',i4)') ntab_t_c+1,ntab_p_c+1
            do i = 0, ntab_p_c
             
             do j = 0, ntab_t_c
                 write(122,'(1p15e14.6)') co2_prop_sw_c(i,j,1:15)
             enddo
          enddo
        close (122)
      
         print *,'Writing Table lookup file: fine...'
            if (myrank==0) print *,'--> open co2data_f.dat'
           open(unit=122,file='co2data_f.dat',status='unknown')
             write(122,'(''TITLE= "'',''co2data_c.dat'',''"'')')
             write(122,'(''VARIABLES= "'',a6,100(a3,a6))') &
          'p',q,'T',q,'d',q,'dddT',q,'dddp',q,'fg',q,'dfgdp',q,'dfgdT',q, &
          'u',q,'h',q,'dhdT',q,'dhdp',q,'vis',q,'dvdT',q,'dvdp','"'
             write(122,'(''ZONE T= "'',''",'','' I='',i4,'' , J='',i4)') ntab_t_f+1,ntab_p_f+1
            do i = 0, ntab_p_f
                do j = 0, ntab_t_f
               write(122,'(1p15e14.6)') co2_prop_sw_f(i,j,1:15)
             enddo
          enddo
        close (122)
     
      endif
 
       
      case(4)
        if (myrank==0) print *,'Reading Table ...'
        if (myrank==0) print *,'--> open co2data0.dat'
        open(unit=122,file='co2data_c0.dat',status='unknown')
        read(122,*)
        read(122,*)
        read(122,*)
        do i = 0, ntab_p_c
          do j = 0, ntab_t_c
            read(122,'(1p15e14.6)') co2_prop_sw_c(i,j,1:15)
          enddo
        enddo
        close (122)

        open(unit=122,file='co2data_f0.dat',status='unknown')
        read(122,*)
        read(122,*)
        read(122,*)
        do i = 0, ntab_p_f
          do j = 0, ntab_t_f
            read(122,'(1p15e14.6)') co2_prop_sw_f(i,j,1:15)
          enddo
        enddo
        close (122)
     end select
      
end subroutine initialize_sw_interp 

! ************************************************************************** !

PetscReal function co2_prop_spwag(ip,it,iv)
     implicit none 
 !    PetscReal co2_prop_spwag
     PetscInt ip,it,iv
      
      
      if (ifinetable)then
         co2_prop_spwag = co2_prop_sw_f(ip,it,iv) 
       else
          co2_prop_spwag = co2_prop_sw_c(ip,it,iv) 
      endif 
     
  end function co2_prop_spwag

! ************************************************************************** !

subroutine interp(x1,x2,y)
  ! 
  ! 2-d function interpolation
  ! 

      use co2_span_wagner_module, only: vappr, co2_span_wagner
      implicit none 
  
#include "petsc/finclude/petscsys.h"
  
      PetscReal :: x1,x2,y(15)   
      
      PetscReal factor(1:4), fac(2,2) , iindex, jindex
      PetscReal funv(1:2,1:2,0:2), funi(1:2,1:2,0:2)
             
      PetscReal ps, tmp, tmp2, ntab_t, ntab_p, dt_tab, dp_tab, p0_tab, t0_tab 
      PetscInt isucc, i1,i2,j1,j2, icross, i,j
      PetscInt :: iflag = 1
      
      ifinetable = PETSC_FALSE
      if (x2 <= co2_sw_f_t1_tab .and. x1<= co2_sw_f_p1_tab) then
       ! within the fine grid table 
        ifinetable = PETSC_TRUE
      endif   


    isucc= 0 
    if (ifinetable)then
      ntab_t = ntab_t_f
      ntab_p = ntab_p_f
      dt_tab = dt_tab_f 
      dp_tab = dp_tab_f
      p0_tab =  co2_sw_f_p0_tab
      t0_tab =  co2_sw_f_t0_tab
!     print *, 'using fine table', x1,x2
    else
      ntab_t = ntab_t_c
      ntab_p = ntab_p_c
      dt_tab = dt_tab_c 
      dp_tab = dp_tab_c
      p0_tab =  co2_sw_c_p0_tab
      t0_tab =  co2_sw_c_t0_tab
    endif


    tmp = (x1 - p0_tab) / dp_tab; i1 = floor(tmp); i2 = i1+1; iindex=tmp 
    tmp = (x2 - t0_tab) / dt_tab; j1 = floor(tmp); j2 = j1+1; jindex=tmp 

    isucc=1

 ! Check whether the table block covers the saturation line, missed special case.
    icross =0
    if (ifinetable) then
      call vappr(co2_prop_spwag(i1,j1,TWO_INTEGER),ps,tmp,tmp2,ELEVEN_INTEGER)
      if ((ps - co2_prop_spwag(i1,j1,ONE_INTEGER)) * (ps - co2_prop_spwag(i1,j2,ONE_INTEGER)) <0.D0)then
        icross = 1; isucc=0
      else
        call vappr(co2_prop_spwag(i2,j1,TWO_INTEGER),ps,tmp,tmp,ELEVEN_INTEGER)
        if ((ps - co2_prop_spwag(i2,j1,ONE_INTEGER)) * (ps - co2_prop_spwag(i2,j2,ONE_INTEGER)) <0.D0)then
          icross = 1; isucc=0
    !     else
    !     call vappr(0.5D0 * (co2_prop_spwag(i1,j1,TWO_INTEGER)+co2_prop_spwag(i2,j1,TWO_INTEGER)),ps, tmp,tmp2,TWELVE_INTEGER)
    !     if ((ts - co2_prop_spwag(i2,j1,TWO_INTEGER)) * (ts - co2_prop_spwag(i2,j2,TWO_INTEGER)) <0.D0)then
    !       icross = 1; isucc=0
    !     endif
        endif
      endif   
    endif

    if (icross == 1) print *,'co2_sw: cross sat line'

    if (iindex > ntab_p .or. iindex < 0.d0 .or. jindex < 0.d0 .or. jindex > ntab_t) then
      print  *,' Out of Table Bounds (interp): ', 'p=',x1,' t=',x2,' i=',iindex,' j=',jindex
      isucc=0
    endif


    if (isucc>0 .and. icross ==0 )then

      factor(1)= (iindex-i2) * (jindex-j2)
      factor(2)= -(iindex-i1) * (jindex-j2)
      factor(3)= -(iindex-i2) * (jindex-j1)
      factor(4)= (iindex-i1) * (jindex-j1)

      fac(1,1) = (iindex-i1 ); fac(2,1) = 1.D0 -  fac(1,1) 
      fac(1,2) = (jindex-j1 ); fac(2,2) = 1.D0 -  fac(1,2)
 
       !  print *,  icross,isucc, factor, fac 

      i=1
      y(1) = factor(1)*co2_prop_spwag(i1,j1,i) + factor(2)*co2_prop_spwag(i2,j1,i) &
           + factor(3)*co2_prop_spwag(i1,j2,i) + factor(4)*co2_prop_spwag(i2,j2,i)
      if (dabs(y(1)-x1)>1D-10 ) then
        print *,' Error in intropolate::P',x1,iindex,factor;isucc=0
      endif
     !print *, 'Table: P ',iindex,jindex, factor,i  
      i=i+1
      y(2) = factor(1)*co2_prop_spwag(i1,j1,i) + factor(2)*co2_prop_spwag(i2,j1,i) &
           + factor(3)*co2_prop_spwag(i1,j2,i) + factor(4)*co2_prop_spwag(i2,j2,i)
      if (dabs(y(2)-x2)>1D-10 ) then
        print *,' Error in intropolate:;T', x2,jindex,factor; isucc=0
      endif
    endif

  ! print *, 'Table: T',iindex,jindex,factor,i,isucc,itable


    if (isucc==1)then
       
      do i =3,15
      ! if (i==var_index(1) .or. i==var_index(2) .or.i==var_index(3) .or.i==var_index(4)) cycle
        y(i) = factor(1)*co2_prop_spwag(i1,j1,i) + factor(2)*co2_prop_spwag(i2,j1,i) &
           + factor(3)*co2_prop_spwag(i1,j2,i) + factor(4)*co2_prop_spwag(i2,j2,i)
      enddo 
       
!#if 0      
      do j=  1, 4
       i= var_index(j)
         funv(1,1,0) = co2_prop_spwag(i1,j1,i) 
         funv(1,2,0) = co2_prop_spwag(i1,j2,i)
         funv(2,1,0) = co2_prop_spwag(i2,j1,i)
         funv(2,2,0) = co2_prop_spwag(i2,j2,i)
          
       select case(i)
       case(3,10,13)     
         funv(1,1,1) = co2_prop_spwag(i1,j1,i+2)   
         funv(1,2,1) = co2_prop_spwag(i1,j2,i+2)
         funv(2,1,1) = co2_prop_spwag(i2,j1,i+2)
         funv(2,2,1) = co2_prop_spwag(i2,j2,i+2)
         funv(1,1,2) = co2_prop_spwag(i1,j1,i+1)   
         funv(1,2,2) = co2_prop_spwag(i1,j2,i+1)
         funv(2,1,2) = co2_prop_spwag(i2,j1,i+1)
         funv(2,2,2) = co2_prop_spwag(i2,j2,i+1)
        case(6)
         funv(1,1,1) = co2_prop_spwag(i1,j1,i+1)   
         funv(1,2,1) = co2_prop_spwag(i1,j2,i+1)
         funv(2,1,1) = co2_prop_spwag(i2,j1,i+1)
         funv(2,2,1) = co2_prop_spwag(i2,j2,i+1)
         funv(1,1,2) = co2_prop_spwag(i1,j1,i+2)   
         funv(1,2,2) = co2_prop_spwag(i1,j2,i+2)
         funv(2,1,2) = co2_prop_spwag(i2,j1,i+2)
         funv(2,2,2) = co2_prop_spwag(i2,j2,i+2)
       end select          

    
        funi(:,:,0)= funv(:,:,0)   
     
        funi(1,1,1) = funv(1,1,1) - (funv(2,1,0) -funv(1,1,0)) /dp_tab
        funi(1,1,2) = funv(1,1,2) - (funv(1,2,0) -funv(1,1,0)) /dt_tab
        funi(2,1,1) = funv(2,1,1) - (funv(2,1,0) -funv(1,1,0)) /dp_tab
        funi(2,1,2) = funv(2,1,2) - (funv(2,2,0) -funv(2,1,0)) /dt_tab
        funi(1,2,1) = funv(1,2,1) - (funv(2,2,0) -funv(1,2,0)) /dp_tab
        funi(1,2,2) = funv(1,2,2) - (funv(1,2,0) -funv(1,1,0)) /dt_tab
        funi(2,2,1) = funv(2,2,1) - (funv(2,2,0) -funv(1,2,0)) /dp_tab
        funi(2,2,2) = funv(2,2,2) - (funv(2,2,0) -funv(2,1,0)) /dt_tab
  
         y(i) = fac(2,2) *( fac(2,1)*(funi(1,1,0) + fac(1,1)* fac(2,1)* funi(1,1,1)+ fac(1,2)* fac(2,2)* funi(1,1,2))&    
                           +fac(1,1)*(funi(2,1,0) - fac(1,1)* fac(2,1)* funi(2,1,1)- fac(1,2)* fac(2,2)* funi(2,1,2)))&  
              + fac(1,2) *( fac(2,1)*(funi(1,2,0) + fac(1,1)* fac(2,1)* funi(1,2,1)+ fac(1,2)* fac(2,2)* funi(1,2,2))&    
                           +fac(1,1)*(funi(2,2,0) - fac(1,1)* fac(2,1)* funi(2,2,1)- fac(1,2)* fac(2,2)* funi(2,2,2)))  
        enddo
!#endif
       endif
       
      if ((icross ==1 .and. ifinetable).or. isucc < 1)then
      ! print *, ' Exit table looking', icross, isucc, ifinetable
       call co2_span_wagner(x1,x2,y(3),y(4),y(5),y(6),y(7),y(8), &
        y(9),y(10),y(11),y(12),y(13),y(14),y(15),ZERO_INTEGER,iflag)
       
      endif  
end subroutine 

! ************************************************************************** !

subroutine co2_sw_interp(p,tc,rho,dddt,dddp,fg,dfgdp,dfgdt, &
        eng,ent,dhdt,dhdp,visc,dvdt,dvdp,itable)
         
       implicit none
       
      PetscReal :: p,tc,rho,eng,ent,dhdt,dhdp,dddt,dddp,visc,dvdt,dvdp
      PetscReal :: fg,dfgdp,dfgdt
      PetscInt itable
       
      PetscReal :: t, prop(15)
     
  
      t=tc+273.15D0
      prop(1)=p
      prop(2)=t
       
      call  interp(p,t,prop)       
      
      rho = prop(3)
      dddt = prop(4)
      dddp = prop(5)
      fg = prop(6)
      dfgdp = prop(7)
      dfgdt = prop(8)
      eng = prop(9)
      ent = prop(10)
      dhdt = prop(11)
      dhdp = prop(12)
      visc = prop(13)
      dvdt = prop(14)
      dvdp = prop(15)
  
!  contains
end subroutine co2_sw_interp 

      
end module co2_sw_module     
