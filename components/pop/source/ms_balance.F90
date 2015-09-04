 module ms_balance
 
!BOP
! !MODULE: ms_balance

! !DESCRIPTION:
!  This module contains routines necessary for the balancing of evaporation,
!  precipitation, melt, runoff, and salt in marginal-seas regions of the
!  ccsm coupled ocean model
!
! !REVISION HISTORY:
!  SVN:$Id: ms_balance.F90 38648 2012-07-12 19:51:27Z mlevy@ucar.edu $
!
! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod

   use kinds_mod
   use domain_size
   use domain
   use global_reductions
   use gather_scatter
   use grid   
   use communicate
   use io_tools
   use ice
   use constants
   use time_management
 
   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_ms_balance, &
             ms_balancing
 
! !PUBLIC DATA MEMBERS:

   real (r8), dimension(:,:,:,:), allocatable :: MASK_FRAC

   real (r8), dimension(max_regions) :: &
     annual_depth                      ,&
     monthly_depth

!EOP
!BOC
!EOC
!***********************************************************************

 contains
 
 
!***********************************************************************
!BOP
! !IROUTINE:  init_ms_balance
! !INTERFACE:

   subroutine init_ms_balance
 
! !DESCRIPTION:
! Provides initialization for marginal-sea balancing
!
!  For a given marginal sea, determine the corresponding active-ocean 
!  distribution region for balancing evaporation, precipitation, melt, and runoff

! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: max_points = 5000, &
                                    max_iter   =  16  
 
   logical (log_kind)           :: & ! Is this point a:
     duplicate                    ,& !   duplicate point?
     marginal_sea                 ,& !   marginal sea point?
     land                         ,& !   land point?   
     debug1 = .false.             ,& ! Print level-1 debug statements 
     debug2 = .false.                ! Print level-2 debug statements 
 
   integer (int_kind)           :: &
     n                            ,& ! region number
     i,j                          ,& ! global do-loop indices
     npts                         ,& ! number of redistribution points
     iter, nn, pt                 ,& ! do-loop indices
     region_number                ,&  
     mask_index = 0                 
 
   integer (int_kind), dimension (max_points)   :: &
     ipts                         ,& ! storage for i&j points in 
     jpts                            !  distribution region
 
 
   integer (int_kind), allocatable, dimension(:,:) :: &
     REGION_MASK_G                   ! global region mask
 
   real (r8)                    :: &
     mblat, mblon                 ,& ! center of lat & lon search
     lat_reach, lon_reach         ,& ! search reach in lat & lon
     area                         ,& ! search area
     sum_frac                        ! sum of distribution fractions
 
   real (r8),dimension (max_points) :: &
     areas                         ,& ! distribution areas     
     fracs                            ! distribution fractions 
 
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
     WORK
 
   real (r8), allocatable, dimension(:,:) :: &
     TLAT_G                        ,& ! global latitude  of cell center
     TLON_G                        ,& ! global longitude of cell center
     AREAT_G                       ,& ! global areas centered at T points
     MASK_G                           ! global redistribution mask
 
 
   if (my_task == master_task) then
     write(stdout,delim_fmt)
     write(stdout,blank_fmt)
     write(stdout,'(a)') ' Marginal-sea balancing information'
     write(stdout,blank_fmt)
     write(stdout,delim_fmt)
     call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
   endif

   allocate (MASK_FRAC (nx_block,ny_block,max_blocks_clinic,max_ms) )

   annual_depth  = c0
   monthly_depth = c0

 
!-----------------------------------------------------------------------
!     create global arrays of TLON, TLAT, AREAT_G, and REGION_MASK
!-----------------------------------------------------------------------
 
   allocate (REGION_MASK_G(nx_global,ny_global))
   allocate (TLON_G  (nx_global,ny_global), &
             TLAT_G  (nx_global,ny_global), &
             AREAT_G (nx_global,ny_global), &
             MASK_G  (nx_global,ny_global)  )


   call gather_global(TLON_G, TLOND, master_task,distrb_clinic) !global TLON in degrees
   

   call gather_global(TLAT_G, TLATD, master_task,distrb_clinic) !global TLAT in degrees
   
   WORK = DXT*DYT
   call gather_global(AREAT_G , WORK,  master_task,distrb_clinic)
 
   call gather_global(REGION_MASK_G,REGION_MASK,master_task,distrb_clinic)

   if (debug1) call print_regions (REGION_MASK_G)
 
 
!-----------------------------------------------------------------------
!     for each marginal sea, determine associated distribution region:
!-----------------------------------------------------------------------

   marginal_seas: do n=1,num_regions
 
   if (region_info(n)%marginal_sea) then

   mask_index = mask_index + 1
 
   MASK_G = c0
!-----------------------------------------------------------------------
!     search globally for associated regions 
!-----------------------------------------------------------------------
   if (my_task == master_task) then
     npts  =  0
     ipts  =  0
     jpts  =  0
     areas = c0
     fracs = c0
     area  = c0
 
     mblat =  region_info(n)%ms_bal%lat
     mblon =  region_info(n)%ms_bal%lon
 
 
!-----------------------------------------------------------------------
!      set initial extent of search in latitude and longitude
!-----------------------------------------------------------------------
     lat_reach = c1
     lon_reach = c1
 
 
     iter_loop: do iter = 1, max_iter
          
!-----------------------------------------------------------------------
!      select active-ocean points which lie within search area
!-----------------------------------------------------------------------
       j_loop: do j=1,ny_global
       i_loop: do i=1,nx_global

         if (mblat-lat_reach <= TLAT_G(i,j)    .and. &
             mblat+lat_reach >= TLAT_G(i,j)    .and. &
             mblon-lon_reach <= TLON_G(i,j)    .and. &
             mblon+lon_reach >= TLON_G(i,j)  ) then
               
             !-----------------------------------------------------
             ! is this point identical to a previously selected one?
             !-----------------------------------------------------
             duplicate = .false.
 
             dup_loop: do nn=1,npts
                if (ipts(nn) == i  .and. jpts(nn) == j) then
                    duplicate = .true.
                    exit dup_loop
                endif
             end do dup_loop
               
             !-----------------------------------------------------
             ! is this point a marginal-sea point?
             !-----------------------------------------------------
             if (REGION_MASK_G(i,j) < 0) then
                 marginal_sea = .true.
             else
                 marginal_sea = .false.
             endif
               
             !-----------------------------------------------------
             ! is this point a land point?
             !-----------------------------------------------------
             if (REGION_MASK_G(i,j) == 0) then
                 land = .true.
             else
                 land = .false.
             endif

             !--------------------------------------------------
             ! reject duplicate points
             !--------------------------------------------------
             if (duplicate) then
                if(debug2) call pt_print('reject duplicate point',i,j) 
 
             !--------------------------------------------------
             ! reject marginal-sea points
             !--------------------------------------------------
             else if (marginal_sea) then
                if(debug2) &
                   call pt_print('reject marginal sea point',i,j) 
    
             !--------------------------------------------------
             ! reject land points
             !--------------------------------------------------
             else if (land) then             
                if(debug2) call pt_print('reject land point',i,j) 
 
             !--------------------------------------------------
             ! select unique active-ocean points
             !--------------------------------------------------
             else   
             !--------------------------------------------------
             !   has maximum number of distribution points been
             !   selected?
             !--------------------------------------------------
                if (npts ==  max_points) then
                    exit iter_loop     
                endif
 
                npts        = npts + 1
                ipts (npts) = i
                jpts (npts) = j
                areas(npts) = AREAT_G(i,j)
                area        = area + areas(npts)
 
                !---------------------------------------------------
                !  is search area acceptable size?
                !---------------------------------------------------
                if( area >= region_info(n)%ms_bal%area)then
                    if (debug1) &
                    write(stdout,1002)'(init_ms_balance) ', &
                      'search area = ', region_info(n)%ms_bal%area
                    exit iter_loop   
                endif
             endif
         endif 
       end do i_loop
       end do j_loop
 
         
       !------------------------------------------------------------
       !   increase search area 
       !------------------------------------------------------------
       lat_reach = lat_reach + c1 
       lon_reach = lon_reach + c1
 
     end do iter_loop
           
       
!-----------------------------------------------------------------------
!     end of master_task region
!-----------------------------------------------------------------------
    endif
 
!-----------------------------------------------------------------------
!     spread the news
!-----------------------------------------------------------------------
   call broadcast_scalar(area , master_task)
   call broadcast_scalar(npts , master_task)
   call broadcast_array (ipts , master_task)
   call broadcast_array (jpts , master_task)
   call broadcast_array (areas, master_task)
 
 
!-----------------------------------------------------------------------
!     distribution points are selected; now test for reasonableness
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!     has a non-zero distribution area been selected?
!-----------------------------------------------------------------------
     if (npts <= 0) then
         call document ('init_ms_balance', 'marginal sea ' /&
                     &/  trim(region_info(n)%name) )
         call document ('init_ms_balance', &
                        'no points selected for distribution area')
         call exit_POP (sigAbort,  &
                  'must select at least one set of active points')
     endif
 
 
!-----------------------------------------------------------------------
!     might there be other points within the distribution area
!     that were not selected because the search loop exited
!     upon reaching the maximum number of points limit?
!-----------------------------------------------------------------------
     if (npts ==  max_points) then
         call document ('init_ms_balance', 'marginal sea ' /&
                    &/   trim(region_info(n)%name) )
         call document ('init_ms_balance', &
                        'warning: an increase of max_points may be necessary')
     endif
         
 
!-----------------------------------------------------------------------
!     are all of the distribution points within one region?
!-----------------------------------------------------------------------
    if (my_task == master_task) then
     region_number = REGION_MASK_G(ipts(1),jpts(1))
     do pt=1,npts
       if(region_number /= REGION_MASK_G(ipts(pt),jpts(pt)))then
          write(stdout,1000) '(init_ms_balance) ','marginal sea '/&
                          &/ trim(region_info(n)%name)
        write(stdout,1000) '(init_ms_balance)', &
                           ' WARNING: distribution points span two regions' 
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
          if (pt > 1) then
            write(stdout,*)'(init_ms_balance)', ' region1 : ' ,&
                             REGION_MASK_G(ipts(pt-1),jpts(pt-1)) 
            write(stdout,*) '(init_ms_balance)', ' region2 : ', &
                             REGION_MASK_G(ipts(pt),jpts(pt)) 
            call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
          endif
        endif
     enddo
    endif
 
!-----------------------------------------------------------------------
!     define distribution fractions; ensure sum = c1
!-----------------------------------------------------------------------
    sum_frac = c0
    do pt=1,npts
       fracs(pt) = areas(pt)/area
       sum_frac   = sum_frac + fracs(pt)
    enddo
 
    fracs(npts) = fracs(npts) + (c1 - sum_frac)
    
 
!-----------------------------------------------------------------------
!     define the distribution mask (globally, then locally)
!-----------------------------------------------------------------------
   if (my_task == master_task) then
    do pt=1,npts
      j_loop2: do j=1,ny_global
        do i=1,nx_global
         if (i == ipts(pt)  .and. j == jpts(pt) ) then
             MASK_G(i,j) = fracs(pt)
             exit j_loop2
         endif
        enddo
      end do j_loop2
    enddo
   endif
       
   call scatter_global(MASK_FRAC(:,:,:,mask_index),MASK_G,master_task, &
                       distrb_clinic, &
                       field_loc_center, field_type_scalar)

 
   region_info(n)%ms_bal%area       = area
   region_info(n)%ms_bal%mask_index = mask_index
 

!-----------------------------------------------------------------------
!     document the distribution points
!-----------------------------------------------------------------------
   if (debug1) then
     if (my_task == master_task) then
      write(stdout,1002)'(init_ms_balance) ','actual area = ',area 
      write(stdout,*) trim(region_info(n)%name)
      write(stdout,*)'pt   ipts   jpts   TLAT  TLON  REGION_MASK frac'
      do pt=1,npts
         write(stdout,*) pt, ipts(pt), jpts(pt),      & 
                    TLAT_G(ipts(pt),jpts(pt)),        &
                    TLON_G(ipts(pt),jpts(pt)),        &
                    REGION_MASK_G(ipts(pt),jpts(pt)), &
                    fracs(pt)
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      enddo
     endif

     call document('init_ms_balance', 'area       = ', area       )
     call document('init_ms_balance', 'npts       = ', npts       )
     call document('init_ms_balance', 'mask_index = ', mask_index )
 
    endif
 
!-----------------------------------------------------------------------
!   end  region_info(n)%marginal_sea if-block
!-----------------------------------------------------------------------
    end if 
 
!-----------------------------------------------------------------------
!   end marginal-seas do-loop
!-----------------------------------------------------------------------
    end do marginal_seas
 
   deallocate (REGION_MASK_G)
   deallocate (TLON_G,  &
               TLAT_G,  &
               AREAT_G, &
               MASK_G   )

1000  format(5x,'(',a,')', a )
1002  format(5x,'(',a,')  ', a ,1x, 1pe15.5)
 
!-----------------------------------------------------------------------
!EOC
   end subroutine init_ms_balance

 
!***********************************************************************
!BOP
! !IROUTINE:  ms_balancing
! !INTERFACE:

   subroutine ms_balancing (STF2, EVAP_F, PREC_F, MELT_F,ROFF_F, IOFF_F, &
                            SALT_F, QFLUX, flux_type, ICEOCN_F)
 
! !DESCRIPTION:
!
!    The total excess or deficit of freshwater (kg/s) for each marginal sea 
!      is transported to or from its associated active-ocean region 
!
!    The transport term, T, is computed in kg/s of freshwater:
!
!EOP
!        T = [ max(0,QFLUX)*c_q +(E+P+M+R+I)*c_f +S*c_s ]*DXT*DYT*c_a
!        T == 0  over marginal seas
!
!        c_f converts E+P+M+R+I kg/m^2/s freshwater to kg/m^2/s freshwater
!        c_f = c1
!
!        c_s converts S       kg/m^2/s salt       to kg/m^2/s freshwater
!        c_s = -(1.0e3_r8/ocn_ref_salinity)*rho_fw/rho_sw
!                       salinity = salt(g)/saltwater(kg)
!
!        c_q converts QFLUX   W/m^2    heat       to kg/m^2/s freshwater
!        c_q = -(1.0e4_r8/latent_heat_fusion)*
!               (ocn_ref_salinity - sea_ice_salinity)/ocn_ref_salinity
!
!        c_a converts T       kg/m^2/s * cm^2     to kg/s     freshwater
!        c_a = 1.0e-4_r8
!
!    The total transport is the sum of T kg/s of freshwater. The transport
!      is distributed over the designated active-ocean regions
!      in freshwater kg/s/cm^2, and must also be converted to the same
!      units as STF(:,:,2)
!
!        c_t converts flux kg/s/cm^2 to the same units as STF(:,:,2)
!           if S(:,:,2) is a salt flux in msu*cm/s
!               c_t = -ocn_ref_salinity/rho_fw 
!           if S(:,:,2) is a freshwater flux in g/cm^2/s
!               c_t = 1.0e3_r8 
!-----------------------------------------------------------------------
!BOP
 
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
 
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic),intent(in)  :: &
     EVAP_F                     ,& ! evaporation   flux kg/m^2/s  fw
     PREC_F                     ,& ! precipitation flux kg/m^2/s  fw
     MELT_F                     ,& ! snow&ice melt flux kg/m^2/s  fw
     ROFF_F                     ,& ! river runoff  flux kg/m^2/s  fw
     IOFF_F                     ,& ! ice   runoff  flux kg/m^2/s  fw
     SALT_F                     ,& ! salt          flux kg/m^2/s  salt
     QFLUX

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), optional, &
     intent(in) :: ICEOCN_F

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic),intent(inout):: &
     STF2                          ! contains STF(:,:,2) from forcing_coupled

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real (r8)  :: &
     c_q,        &
     c_f,        &
     c_s,        &
     c_a,        &
     c_t

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic)  :: &
      WORK,                        &! temporary work array
      TRANSPORT                     ! pointwise transport term

   character (*)  ::  &
      flux_type        ! flux type of STF(:,:,2) (freshwater or salt )

   integer (int_kind) ::  &
      n,                  &
      iblock

   logical (log_kind)          :: &
      print_transport_daily    = .false. ,&
      print_transport_monthly  = .true.  ,&
      print_transport_annually = .false. ,&
      iceocn_f_present

 
   iceocn_f_present = present(ICEOCN_F)

   c_q = -(1.0e4_r8/latent_heat_fusion)* &
          (ocn_ref_salinity - sea_ice_salinity)/ocn_ref_salinity
   c_f = c1
   c_s = -(1.0e3_r8/ocn_ref_salinity)*(rho_fw/rho_sw)
   c_a = 1.0e-4_r8
 
   select case (flux_type) 
     case ('salt') 
       c_t = -ocn_ref_salinity
     case ('freshwater') 
       c_t = 1.0e3_r8
     case default
          call exit_POP (sigAbort, &
              'ms_balancing -- must set flux_type to either salt or freshwater')
    end select
             
    if (eoy) annual_depth  = c0
    if (eom) monthly_depth = c0
 
    WORK = EVAP_F + PREC_F + MELT_F + ROFF_F + IOFF_F
    if ( iceocn_f_present )  WORK = WORK + ICEOCN_F
 
    do n=1,num_regions
     
     TRANSPORT = c0
 
!-----------------------------------------------------------------------
!     in each marginal sea, determine transport term and distribute to
!        or from the associated active-ocean region
!-----------------------------------------------------------------------
     if (region_info(n)%marginal_sea) then
 
        !-----------------------------------------------------------
        ! accumulate total transport, in kg/s,  for each marginal sea 
        !
        !   T = [( max(0,QFLUX)*c_q + (E+P+M+R)*c_f +S*c_s )]*area*c_a 
        !
        ! and then set STF2 to zero or the opposite of ice formation
        ! flux there
        !-----------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock = 1,nblocks_clinic

         where (REGION_MASK(:,:,iblock) == region_info(n)%number)
           TRANSPORT(:,:,iblock) =                                    &
                                 ( max(c0,QFLUX(:,:,iblock))*c_q      &
                                 + WORK(:,:,iblock)*c_f               &
                                 + SALT_F(:,:,iblock)*c_s )*          &
                                   DXT(:,:,iblock)*DYT(:,:,iblock)*c_a  
           STF2(:,:,iblock) = -max(c0,QFLUX(:,:,iblock))*c_q*c_t*1.0e-4_r8
         end where

         enddo ! iblock
         !$OMP END PARALLEL DO


         region_info(n)%ms_bal%transport = &
            global_sum(TRANSPORT,distrb_clinic,field_loc_center)
 
 
        !-----------------------------------------------------------
        ! accumulate the average depth of freshwater transported
        !-----------------------------------------------------------
         annual_depth(n) = annual_depth(n) + &
          (region_info(n)%ms_bal%transport/region_info(n)%ms_bal%area) &
          *1.0e3_r8*seconds_in_day
 
         monthly_depth(n) = monthly_depth(n) + &
          (region_info(n)%ms_bal%transport/region_info(n)%ms_bal%area) &
          *1.0e3_r8*seconds_in_day
 
           
        !----------------------------------------------------------------
        ! transport excess/deficit to/from associated active-ocean region
        !----------------------------------------------------------------
         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock = 1,nblocks_clinic

         STF2(:,:,iblock) = STF2(:,:,iblock)  &
                     + MASK_FRAC(:,:,iblock,region_info(n)%ms_bal%mask_index) &
                     * region_info(n)%ms_bal%transport*c_t/  &
                      (DXT(:,:,iblock)*DYT(:,:,iblock))
         enddo ! iblock
         !$OMP END PARALLEL DO

     endif
 
    enddo ! num_regions

 
 
    !-----------------------------------------------------------
    ! record the average daily depth of freshwater transported
    !   (total transport/total distribution area)
    !-----------------------------------------------------------
 
    if (print_transport_daily .and. my_task == master_task)  then
      write(stdout,*) ' '   
      call int_to_char (4,iyear ,cyear )
      call int_to_char (2,imonth,cmonth)
      call int_to_char (2,iday  ,cday  )
      write(stdout,*) cyear/&
                            &/  '/'  /&
                                      &/cmonth/&
                                               &/  '/' /&
                                                        &/cday
      do n=1,num_regions
       if (region_info(n)%marginal_sea) then
           write(stdout,1100) &
                 'average depth of freshwater transported to ', &
                 trim(region_info(n)%name) /&
                                            &/ ': ', &
                (region_info(n)%ms_bal%transport/region_info(n)%ms_bal%area)* &
                 1.0e3_r8*seconds_in_day, ' (cm/day)'   
       endif
      enddo
      write(stdout,*) ' '
    endif
 
 
 
    !-----------------------------------------------------------
    ! record the total monthly depth of freshwater transported
    !-----------------------------------------------------------
 
    if (print_transport_monthly  .and. eom_next &
                                 .and. my_task == master_task) then
      write(stdout,*) ' '
      do n=1,num_regions
       if (region_info(n)%marginal_sea) then
           write(stdout,1102) month3_all(imonth) /&
                                                  &/ ' ' /&
                                                          &/ &
          'total monthly depth of freshwater transported to ' , &
           trim(region_info(n)%name) /&
                                      &/ ': ',monthly_depth(n),' (cm)'
       endif
      enddo
      write(stdout,*) ' '
    endif
 
    !-----------------------------------------------------------
    ! record the total annual depth of freshwater transported
    !-----------------------------------------------------------
 
    if (print_transport_annually .and. eom_next .and. imonth_next == 1 &
       .and. my_task == master_task) then
      write(stdout,*) ' '
      do n=1,num_regions
       if (region_info(n)%marginal_sea) then
           write(stdout,1101) &
           iyear,  ' total annual depth of freshwater transported to ' , &
           trim(region_info(n)%name) /&
                                      &/ ': ',annual_depth(n),' (cm)'
       endif
      enddo
      write(stdout,*) ' '
    endif
 
 
1100  format (1x, a45, a20, 1pe25.15, a)
1101  format (1x, i5, a45, a20, 1pe25.15, a)
1102  format (1x, a55, a20, 1pe25.15, a)
 
!-----------------------------------------------------------------------
!EOC
   end subroutine ms_balancing


!***********************************************************************
!BOP
! !IROUTINE:  pt_print
! !INTERFACE:

   subroutine pt_print (string, i, j)


   character (*)              :: string
   integer   (int_kind)       :: i,j

   write(stdout,1000) trim(string), i,j
1000  format (1x,a, 1x, 2i4)

!-----------------------------------------------------------------------
!EOC
   end subroutine pt_print

 
!***********************************************************************
!BOP
! !IROUTINE:  print_regions
! !INTERFACE:

   subroutine print_regions (REGION_MASK_G)
 
! !DESCRIPTION:
!  Print global REGION_MASK
 
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), dimension(nx_global,ny_global) :: REGION_MASK_G 

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (4), dimension(25) :: line
   integer (int_kind)           :: i1, i2, i,j, nn
 
 
   i1 = 1   
   i2 = 25
 
   if (my_task == master_task) then
 
   do nn=1,4
     write(stdout,2000) 'j ', (i, i=i1,i2)
     do j=ny_global,1,-1
        do i=i1,i2
         if (REGION_MASK_G(i,j) == 0) then
           write(line(i-i1+1),'(a4)') '-'
         else
           write(line(i-i1+1),'(i4)') REGION_MASK_G(i,j)
         endif
        enddo ! i
      write (stdout,2001) j, (line(i-i1+1), i=i1,i2)
     enddo ! j
 
     i1 = i2+1;  i2 = i1+24
     write(stdout,*) ' '
     write(stdout,*) ' '
   enddo
 
   endif
 
2001  format(1x, i4,2x, 30a4) 
2000  format(2x, a3,2x, 30i4)
 
!-----------------------------------------------------------------------
!EOC

   end subroutine print_regions
 
 end module ms_balance
