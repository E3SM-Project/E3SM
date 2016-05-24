
      module mo_setcld

      use shr_kind_mod, only : r8 => shr_kind_r8
      use cam_logfile,  only : iulog

      save

      integer,  parameter :: nblk   = 10, &                ! maximum dimension of cloud blks (sing+mu
                             mreg   = 16, &                ! maximum dimension of regions 
                             nmax   = 3                    ! Maximum number of the blocks
      real(r8), parameter :: wden   = 1.0_r8   *  1.e6_r8  ! g/m3  (1 m3 water = 1e6 g water)
      real(r8), parameter :: re     = 10.0_r8  *  1.e-6_r8 ! assuming cloud drop radius = 10 um to M
      real(r8), parameter :: cldmin = 0.01_r8              ! minimum cld cover

      integer :: r1(1,2) = reshape( (/ 0, 1 /),(/1,2/) ) ! 1=blk, 2=binary prob (regions)
      integer :: r2(2,4) = reshape( (/ 0,0, 0,1, 1,1, 1,0/), (/2,4/) )
      integer :: r3(3,8) = reshape( (/ 0,0,0, 0,0,1, 0,1,0, 0,1,1,  &
                                       1,0,0, 1,0,1, 1,1,0, 1,1,1/), (/3,8/) )
      integer :: r4(4,16) = reshape( (/ 0,0,0,0, 0,0,1,0, 0,0,1,1, 0,0,0,1, &
                                        0,1,0,0, 0,1,1,0, 0,1,1,1, 0,1,0,1, &
                                        1,0,0,0, 1,0,1,0, 1,0,1,1, 1,0,0,1, &
                                        1,1,0,0, 1,1,1,0, 1,1,1,1, 1,1,0,1/),(/4,16/) )
      real(r8) :: a1(1,2), b1(1,2)
      real(r8) :: a2(2,4), b2(2,4)
      real(r8) :: a3(3,8), b3(3,8)
      real(r8) :: a4(4,16), b4(4,16)

      private
      public :: setcld, setcld_inti

      contains

      subroutine setcld_inti
!-----------------------------------------------------------------------------
!	... Initialize setcld module
!-----------------------------------------------------------------------------
      
      implicit none

      a1(:,:) = r1(:,:)
      a2(:,:) = r2(:,:)
      a3(:,:) = r3(:,:)
      a4(:,:) = r4(:,:)

      b1(:,:) = (-1)**r1(:,:)
      b2(:,:) = (-1)**r2(:,:)
      b3(:,:) = (-1)**r3(:,:)
      b4(:,:) = (-1)**r4(:,:)


      end subroutine setcld_inti

      subroutine setcld( z, xlwc, cldfrc, nreg, fp, optr)
!-----------------------------------------------------------------------------
!   PURPOSE:                                                                 
!   Set up an altitude profile of ozone, and corresponding absorption        
!   optical depths.  Subroutine includes a shape-conserving scaling method   
!   that allows scaling of the entire profile to a given overhead ozone      
!   column amount.                                                           
!-----------------------------------------------------------------------------
!   PARAMETERS:                                                              
!   PARAMETERS:                                                              
!   NZ      - INTEGER, number of specified altitude levels in the working (I)
!             grid                                                           
!   Z       - REAL, specified altitude working grid (km)                  (I)
!   XLWC     Cloud water content g/M3                                     (I)
!   CLDFRC   Cloud fraction
! 
!   NREG     INTEGER regions #                                            (O)
!                                                                            
!   FP       REAL, probability at each region                             (O)
!                                                                            
!   optr     REAL, optical depth at each region                           (O)
!                                                                            
!-----------------------------------------------------------------------------
!     VERTICAL DOMAIN is from bottom(1)  to TOP (TOP=plevp)
!        CCM from top(1) to bottom(plevp)
!-----------------------------------------------------------------------------
!     Calculate the UV flux at short wave-length 
!       with cloud overlap conditions
!       References - Stubenrauch, J. Climate, 273-1997
!
!               regions
!       n=1    n=2  n=3  n=4(nreg)   
!      ------------------------
!      !    !      !   !      !
!      !    !==========!      !    cld(3),opt(3)  nblk(2)  nki=1  sing-layer (maximum cover=cld(3)
!      !    !      !   !      ! 
!      !    !      ========   !    cld(2),opt(2)
!      !    !      ===========!    cld(1),opt(1)  nblk(1), nki=2  mult-layer (maximum cover=cld(1)
!      !    !      !   !      !
!------------------------------
!
!     PROCESSURES:
!
!     (1) To find cloud layers
!
!     (2) To find single layer cloud         =======
!
!     (3) To find mult adjcent cloud layers  ========
!                                            ===========
!
!     (4) To deteming # (NBLK) of cloud layers (Single + mult)
!
!     (5) To find maximum cloud cover at each cloud block
!
!     (6) To assume in mult-layer, all layer cld cover=maximum cover
!            and cloud opt will be adjusted to be maximum cover
!               opt = opt * cld/cld-max
!     
!     (7) To determine # of regions (NREG)
!                nreg = 2**nblk
!
!     (8) To determine the binary probality of NREG
!           NBLK=1   r1   (nreg=2)
!           NBLK=2   r2   (nreg=4)
!           NBLK=3   r3   (nreg=8)
!           NBLK=4   r4   (nreg=16)
!
!         The maximum nblk in each horizontal grid = 4, nblk<=4
!         to elimanite the nblk>4, find smaller opt, and let opt(small=0)
!           
!     (9) To calculate the probality at each region
!             
!               clear sky = 1-cld
!            so cld-cov = 1-{PRODUCT[1 to nblk] (1-cld)}
!
!     (10) To get opt (wei opt) at each vertical level at each hor grid
!
!     (11) input opt in rtlink for flux at each region
!             need a do-loop for call rtlink
!             to get flxreg(nreg=1,2,...)
!      
!     (12) Adding each flux at each region weighted by
!             probability of each region
!    
!             FLUX = flxreg(1)*prob(1) + flxreg(2)*prob(2) + ....
!----------------------------------------------------------------------------

      use mo_params
      use mo_waveo3
      use ppgrid, only: pver, pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      real(r8), intent(in)  :: z(pverp)
      real(r8), intent(in)  :: xlwc(pverp)
      real(r8), intent(in)  :: cldfrc(pverp)
      real(r8), intent(out) :: optr(pver,mreg)             ! cld opt (z dependent) at each region
      real(r8), intent(out) :: fp(mreg)                    ! probality at each region

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: ilat, nn, m, mm
      integer  :: k, l, n                               ! vertical index for each blk
      integer  :: nb, &                                 ! # blk index found
                  nreg                                  ! # regions   found
      integer  :: indx(1)                               ! wrk array
      integer  :: ki(nblk), ki2(nblk)                   ! level count in each cloud block
      integer  :: kl(nblk), ku(nblk)                    ! top, bot level index in each cloud block
      integer  :: kl2(nblk), ku2(nblk)                  ! top, bot level index in each cloud block

      real(r8) :: small
      real(r8) :: wrk
      real(r8) :: cld(pverp)
      real(r8) :: optin(pver)
      real(r8) :: opt(pver)                             ! cld opt   where cloud is found

      real(r8) :: cldadj(pver,nblk), &                  ! cld cover  at each block nblk=1,2,..     
                  zkadj(pver,nblk),  &                  ! blk height at each block nblk=1,2,.. 
                  optadj(pver,nblk), &                  ! cld opt    at each block nblk=1,2,.. 
                  optwe(pver,nblk),  &                  ! cld opt(weighted) at each block nblk=1,2,.. 
                  opta(nblk),        &                  ! total cld opt(weighted) at each block nblk=1,2,.. 
                  zbk(nblk),  &                         ! wrk array
                  maxcld(nblk)                          ! maximum cloud cover at each blk

      real(r8) :: cldadj2(pver,nblk), &                 ! cld cover  at each block nblk=1,2,..     
                  zkadj2(pver,nblk),  &                 ! blk height at each block nblk=1,2,.. 
                  optadj2(pver,nblk)                    ! cld opt    at each block nblk=1,2,.. 

      logical  :: mask(nblk)                            ! wrk array

!-----------------------------------------------------------------------------    
!     	... adjust cloud fraction for cloud liquid water content
!-----------------------------------------------------------------------------    
      do k = 1,pverp
         if( xlwc(k) <= .01_r8  .and. cldfrc(k) /=  0._r8 ) then
            cld(k) = 0._r8
	 else
            cld(k) = cldfrc(k)
         end if
      end do

Have_clouds : &
      if( count( cld(:pver) > cldmin ) > 0 ) then
!-----------------------------------------------------------------------------    
!     	... Find cloud layers
!-----------------------------------------------------------------------------    
!-----------------------------------------------------------------------------    
!  	... calculate  cloud optical depth T
!           following Liao et al. JGR, 104, 23697, 1999
!-----------------------------------------------------------------------------    
         do k = 1,pver
            optin(k)  = max( 1.5_r8 * xlwc(k)*(z(k+1) - z(k))*1.e3_r8/ (wden * re),0._r8 )
         end do
#ifdef DEBUG
	 write(iulog,*) ' '
	 write(iulog,*) '--------------------------------------------------------'
	 write(iulog,*) 'k, z(k), cld(k),xlwc(k),optin(k)'
         do k = 1,pver
            write(iulog,981) k,z(k),cld(k),xlwc(k),optin(k)
         end do
	 write(iulog,*) '--------------------------------------------------------'
#endif
!--------------------------------------------------------
!     	... find cloud layer and make adjcent layer block
!--------------------------------------------------------
         nb = 1
         mm = 0
         do k = 1,pver
            if( cld(k) > cldmin ) then
               mm = mm + 1
               cldadj(mm,nb) = cld(k) 
               zkadj(mm,nb)  = z(k)
               optadj(mm,nb) = optin(k)
               ki(nb) = mm
               ku(nb) = k
	       if( mm == 1 ) then
	          kl(nb) = k
	       end if
               if( cld(k+1) <= cldmin) then
                  mm = 0
		  nb = nb + 1
               end if
            end if
         end do

	 if( mm == 0 ) then
	    nb = nb - 1
	 end if

!--------------------------------------------------------
!     	...Maximum overlap for adjacent cloud layers
!          (1) to find maximum cld cover
!--------------------------------------------------------
         do l = 1,nb
            maxcld(l) = maxval( cldadj(1:ki(l),l) )  !maxcld= maxmum cloud in one block
         end do

#ifdef DEBUG
	 write(iulog,*) ' '
	 write(iulog,*) 'setcld: has ',nb,' cloud layers'
	 write(iulog,*) '        ki,kl,ku'
	 write(iulog,*) ki(:nb)
	 write(iulog,*) kl(:nb)
	 write(iulog,*) ku(:nb)
	 write(iulog,'(1p,10g12.5)') maxcld(:nb)
	 write(iulog,*) ' '
#endif

!--------------------------------------------------------
!     	... limit the total block <= nmax
!--------------------------------------------------------
         if( nb > nmax ) then
#ifdef DEBUG
	    write(iulog,*) ' '
	    write(iulog,*) 'setcld: has ',nb,' cloud layers'
	    write(iulog,*) '        ki,kl,ku'
	    write(iulog,*) ki(:nb)
	    write(iulog,*) kl(:nb)
	    write(iulog,*) ku(:nb)
	    write(iulog,'(1p,10g12.5)') maxcld(:nb)
#endif
            do l = 1,nb
               opta(l) = sum( optadj(1:ki(l),l) )/maxcld(l) 
               zbk(l) = maxcld(l) 
               ki2(l) = ki(l)
               kl2(l) = kl(l)
               ku2(l) = ku(l)
               do k = 1,ki(l)
                  cldadj2(k,l) = cldadj(k,l)
                  zkadj2(k,l)  = zkadj(k,l)
                  optadj2(k,l) = optadj(k,l)
               end do
            end do     
#ifdef DEBUG
	    write(iulog,'(1p,10g12.5)') opta(:nb)
	    write(iulog,*) ' '
#endif
	    mask(:nb) = .true.
	    do l = 1,nmax
	       indx(:)       = maxloc( opta(:nb),mask=mask(:nb) )
	       mask(indx(1)) = .false.
	    end do
            nn = 0
            do l = 1,nb
               if( .not. mask(l) ) then
                  nn = nn + 1
                  ki(nn)     = ki2(l)
                  kl(nn)     = kl2(l)
                  ku(nn)     = ku2(l)
                  maxcld(nn) = zbk(l) 
                  do k = 1,ki(nn)
                     cldadj(k,nn) = cldadj2(k,l)
                     zkadj (k,nn) = zkadj2(k,l)
                     optadj(k,nn) = optadj2(k,l)
                  end do
               end if          
            end do   
            nb = nmax
#ifdef DEBUG
	    write(iulog,*) ' '
	    write(iulog,*) 'setcld: has ',nb,' cloud layers'
	    write(iulog,*) '        ki,kl,ku'
	    write(iulog,*) ki(:nb)
	    write(iulog,*) kl(:nb)
	    write(iulog,*) ku(:nb)
	    write(iulog,'(1p,10g12.5)') maxcld(:nb)
	    write(iulog,*) ' '
#endif
         end if
#ifdef DEBUG
	 write(iulog,*) ' '
	 write(iulog,*) '--------------------------------------------------------'
	 write(iulog,*) 'nb, l, k, zkadj, cldadj, optadj'
         do l = 1,nb
            do k = 1,ki(l)
               write(iulog,992) nb,l,k,zkadj(k,l),cldadj(k,l),optadj(k,l)            
            end do
         end do
	 write(iulog,*) '--------------------------------------------------------'
#endif

!--------------------------------------------------------
!     (2) calculate the weighted cloud optical depth
!--------------------------------------------------------
         do l = 1,nb
            if( maxcld(l) > 0._r8 ) then
	       wrk = 1._r8/maxcld(l)
               do k = 1,ki(l)
                  optwe(k,l) = cldadj(k,l)*optadj(k,l)*wrk
               end do
            else 
               do k = 1,ki(l)
                  optwe(k,l) = 0._r8
               end do
            end if
         end do
!--------------------------------------------------------
!     (3) calculate the probability
!--------------------------------------------------------
         nreg = 2**nb    ! each column of grid has nregions; nreg probobilities
         zbk(:nb) = 1._r8 - maxcld(:nb) 
         if( nb == 1 ) then
            do n = 1,nreg
               fp(n) = product( a1(:nb,n) + b1(:nb,n)*zbk(:nb) )
            end do
         else if( nb == 2 ) then
            do n = 1,nreg
               fp(n) = product( a2(:nb,n) + b2(:nb,n)*zbk(:nb) )
            end do
         else if( nb == 3 ) then
            do n = 1,nreg
               fp(n) = product( a3(:nb,n) + b3(:nb,n)*zbk(:nb) )
            end do
         else if( nb == 4 ) then
            do n = 1,nreg
               fp(n) = product( a4(:nb,n) + b4(:nb,n)*zbk(:nb) )
            end do
         else if( nb == 0 ) then
            do n = 1,nreg
               fp(n) = 1._r8
            end do
         end if
#ifdef DEBUG
	 write(iulog,*) ' '
	 write(iulog,*) '--------------------------------------------------------'
	 write(iulog,*) 'nreg,fp'
         write(iulog,993) nreg,(fp(n),n=1,nreg)
	 write(iulog,*) 'nreg,k,zkadj,optwe'
         do l = 1,nb
            do k = 1,ki(l)
               write(iulog,994) nreg,k,zkadj(k,l),optwe(k,l)
            end do
         end do
	 write(iulog,*) '--------------------------------------------------------'
#endif

!-------------------------------------------------
!     (4) convert to weight optical depth at each region
!-------------------------------------------------
         do n = 1,nreg      
            optr(:,n) = 0._r8
         end do
         if( nb == 1) then
            do n = 1,nreg
               do l = 1,nb
		  if( r1(l,n) /= 0 ) then
                     do k = kl(l),ku(l)
                        optr(k,n) = optwe(k-kl(l)+1,l)
                     end do
		  end if
               end do
            end do
         else if( nb == 2 ) then
            do n = 1,nreg
               do l = 1,nb
                  if( r2(l,n) /= 0 ) then
                     do k = kl(l),ku(l)
                        optr(k,n) = optwe(k-kl(l)+1,l)
                     end do
                  end if
               end do
            end do
         else if( nb == 3 ) then
            do n = 1,nreg
               do l = 1,nb
                  if( r3(l,n) /= 0 ) then
                     do k = kl(l),ku(l)
                        optr(k,n) = optwe(k-kl(l)+1,l)
                     end do
                  end if
               end do
            end do
         else if( nb == 4 ) then
            do n = 1,nreg
               do l = 1,nb
                  if( r4(l,n) /= 0 ) then
                     do k = kl(l),ku(l)
                        optr(k,n) = optwe(k-kl(l)+1,l)
                     end do
                  end if
               end do
            end do
         else if( nb == 0 ) then
            do k = 1,pver
               optr(k,1) = 0._r8
            end do
         end if
      else
!-----------------------------------------------------------------------------    
!     	... no clouds are found
!-----------------------------------------------------------------------------    
         nreg = 1
         do k = 1,pver
           optr(k,1) = 0._r8     ! clear sky
         end do
         fp(1) = 1._r8
      end if Have_clouds

#ifdef DEBUG
      if( nreg > 1 ) then
         write(iulog,*) ' '
         write(iulog,*) '--------------------------------------------------------'
         write(iulog,*) 'fp'
         write(iulog,'(1p,10g12.5)') fp(1:nreg)
         do n = 1,nreg
            write(iulog,*) 'z,optr for nreg = ',n
	    do k = 1,pver
               write(iulog,'(1p,2g15.5)') z(k),optr(k,n)
	    end do
            write(iulog,*) ' '
         end do
         write(iulog,*) '--------------------------------------------------------'
      end if
#endif

 981  format(i10,4g10.3)
 992  format(3i10,3g13.3)
 993  format(i10,10g13.3)
 994  format(2i10,10g13.3)

      end subroutine setcld

      end module mo_setcld
