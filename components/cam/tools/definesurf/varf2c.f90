subroutine varf2c(flon   ,flat  ,nflon   ,nflat   ,fine   , &
                  clon   ,clat  ,nclon   ,nclat   ,cmean  , &
                  cvar  )

  use shr_kind_mod, only: r8 => shr_kind_r8

!-----------------------------------------------------------------------------
!       Bin going from a fine grid to a coarse grid.
!       A schematic for the coarse and fine grid systems is shown in
!       Figure 1.  This code assumes that each data point is represent
!       it's surrounding area, called a cell.  The first grid data point
!       for both grids is assumed to be located at 0E (GM).  This
!       implies that the 1st cell for both the fine and the coarse grids
!       strattles the Greenwich Meridian (GM).  This code also assumes
!       that there is no data wraparound (last data value is located at
!       360-dx).
!
!       FIGURE 1:  Overview of the coarse (X) and fine (@) grids
!                  longitudinal structure where:
!                  X = location of each coarse grid data point
!                  @ = location of each fine   grid data point
!
!           Greenwich                                     Greenwich
!              0          Coarse cells                       360
!              :          v                                   :
!       clon(1):  clon(2) v clon(3)                clon(nclon):
!         v    :    v     v   v                        v      :
!         xxxxxxxxxxxxxxxxxxxxxxxxxxxx..xxxxxxxxxxxxxxxx      :
!         x         x         x              x         x      :
!         x         x         x              x         x      :
!         x  c(1)   x  c(2)   x              x c(nclon)x      :
!         x    X    x    X    x              x    X    x      :
!         x   ___ ___ ___ ___ ___ ___    ___ ___ ___ ___ ___  :
!         x  |   |   |   |   |   |   |  |   |   |   |   |   | :
!         x  | @ | @ | @ | @ | @ | @ |..| @ | @ | @ | @ | @ | :
!         xxx|___|___|___|___|___|___|  |___|___|___|___|___| :
!              v       v       v                      v   v   :
!          flon(1)   flon(3)   v            flon(nflon-1) flon(nflon)
!              :               v                              :
!              :               Fine cells                     :
!              0                                             360
!
!    The Longitude/Latitude search:
!    ------------------------------
!
!      Given a coarse grid cell with west and east boundaries of cWest
!      and cEast and south and north boundaries of cSouth and cNorth
!      (outlined by "x" in figure 2), find the indices of the fine grid
!      points which are contained within the coarse grid cell. imin and 
!      imax are the indices fine grid points which overlap the western
!      and eastern boundary of the coarse cell.  jmin and jmax are the
!      corresponding indices in the S-N direction.  Bin these overlapping
!      values to generate coarse(n), the coarse grid data values.
!
!        FIGURE 2: Detail of Coarse and Fine cell overlap.
!                   @ = fine   grid data point
!                   X = coarse grid data point
!
!                             cWest             cEast
!                   |       |   x   |       |     x |
!                  -@-------@---x---@-------@-----x-@-
!                   |       |  x*xxxxxxxxxxxxxxxxx*x|xx cNorth
!                   |       |   x   |       |     x |
!                   |       |   x   |       |     x |
!                   @-------@---x---@-------@-----x-@- jmax
!                   |       |   x   |  c(n) |     x |
!                   |   @   |   x   |       |     x |
!                   |       |   x   |       |     x |
!                   @-------@---x---@-------@-----x-@- jmin
!                   |       |   x   |       |     x |
!                   |   @   |  x*xxxxxxx@xxxxxxxxx*x|xx cSouth
!                   |       |   x   |       |     x |
!                  -@-------@---x---@-------@-----x-@-  
!                   |             imin    imax      |
!
!
!      When a cell coarse cell strattles the Greenwich Meridian
!      ---------------------------------------------------------
!
!      The first coarse grid cell strattles the GM, so when the western
!      boundary of the coarse cell is < 0, an additional search is carried out. 
!      It ASSUMES that the easternmost fine grid point overlaps and searches
!      westward from nflon, looking for a grid point west of clon(1)
!      This generates a second set of longitudinal indices, imin1 and imax1.   
!      See Figure 3.
!
!       Figure 3:  Detail of Coarse cell strattling GM:
!       -----------------------------------------------
!
!              Greenwich                                   Greenwich
!                 0                                           360
!          cWest  :  cEast                            cWest    :
!          clon(1):  clon(2)               clon(nclon+1)=clon(1)
!            v    :    v                                v      :
!            xxxxxxxxxxxxxxxxxxxxxxx ... xxxxxxxxxxxxxxxx      :
!            x         x        x            x          x      :
!            x         x        x            x          x      :
!            x  c(1)   x        x            x  c(nclon)x      :
!            x    X    x        x            x     X    x      :
!            x   ___ ___ ___ _                    ___ ___ ___  :
!            x  |   |   |   |                        |   |   | :
!            x  | @ | @ | @ |                      @ | @ | @ | :
!            xxx|___|___|___|_                    ___|___|___| :
!                 ^ : ^   ^                            ^   ^   :
!            flon(1): ^ flon(3)              flon(nflon-1) ^   :
!                 ^ : ^                                ^   ^   :
!                 ^ :flon(2)                           ^ flon(nflon)
!                 ^ : ^                                ^   ^   :
!              imin : imax                         imin1 imax1 :
!                   :                                          :
!
!
!         In this case, imin=1, imax=2, imin1=nflon-1 and imax1=nflon.
!         because the last two cells of the fine grid will have some
!         contribution the the 1st cell of the coarse grid.
!
!-----------------------------------------------------------------------
  implicit none
!-----------------------------Arguments---------------------------------

  integer nflon           ! Input: number of fine longitude points
  integer nflat           ! Input: number of fine latitude points
  integer nclon           ! Input: number of coarse longitude points
  integer nclat           ! Input: number of coarse latitude points

  real(r8) flon(nflon)        ! Input: fine grid lons, centers (deg)
  real(r8) flat(nflat)        ! Input: fine grid lats, centers (deg)
  real(r8) fine(nflon,nflat)  ! Input: Fine grid data array
  real(r8) clon(nclon+1,nclat) ! Input: coarse grid cell lons, west  edge (deg)
  real(r8) clat(nclat+1)      ! Input: coarse grid cell lat,  south edge (deg)
  real(r8) cmean(nclon,nclat) ! Input: mean of fine points over coarse grid cell
  real(r8) cvar (nclon,nclat) ! Output:variance of fine points over coarse cell

!--------------------------Local variables------------------------------

  real(r8) cWest            ! Coarse cell longitude, west edge (deg)
  real(r8) cEast            ! Coarse cell longitude, east edge (deg)
  real(r8) cSouth           ! Coarse cell latitude, south edge (deg)
  real(r8) cNorth           ! Coarse cell latitude, notrh edge (deg)
  real(r8) sum              ! coarse tmp value
      
  integer i,j             ! Indices
  integer imin ,imax      ! Max/Min E-W indices of intersecting fine cell.
  integer imin1,imax1     ! fine E-W indices when coarse cell strattles GM
  integer jmin ,jmax      ! Max/Min N-S indices of intersecting fine cell.
  integer iclon,jclat     ! coarse grid indices
  integer num             ! increment 

!-----------------------------------------------------------------------------

  do jclat= 1,nclat         ! loop over coarse latitudes
    cSouth = clat(jclat)
    cNorth = clat(jclat+1)

    do iclon=1,nclon       ! loop over coarse longitudes
      cWest  = clon(iclon,jclat)
      cEAST  = clon(iclon+1,jclat)

!  1. Normal longitude search:  Find imin and imax

      imin = 0
      imax = 0
      do i=1,nflon-1                     ! loop over fine lons, W -> E
        if (flon(i) .gt. cEast) goto 10 ! fine grid point is E of coarse box
        if (flon(i) .ge. cWest .and. imin.eq.0) imin=i
        imax=i
      enddo

!  2. If cWest < 0, then coarse cell strattles GM.  Hunt westward
!     from the end to find indices of any overlapping fine grid cells:
!     imin1 and imax1.

10    imin1 = 0           ! borders for cWest, cEast
      imax1 = -1          ! borders for cWest, cEast
      if (cWest .lt. 0) then
        cWest = cWest + 360.
        imax1 = nflon
        do i=nflon,1,-1                     ! loop over fine lons, E -> W
          imin1=i
          if (flon(i) .le. cWest) goto 20 ! fine grid point is W of coarse box 
        enddo
      endif

!  3. Do the latitude search S -> N for jmin and jmax
      
20    jmin = 0
      jmax = 0
      do  j=1,nflat                       ! loop over fine lats, S -> N
        if (flat(j) .gt. cNorth) goto 30 ! fine grid point is N of coarse box
        if (flat(j) .ge. cSouth .and. jmin.eq.0) jmin=j
        jmax=j
      enddo
30    continue

!  4. Sdv

      sum = 0.                        ! Initialize coarse data value
      num  = 0

      do j=jmin,jmax                   ! loop over fine lats, S -> N
        do i=imin,imax                ! loop over fine lons, W -> E
          sum = sum + (fine(i,j) - cmean(iclon,jclat))**2
          num  = num + 1
        enddo
        do i=imin1,imax1              ! If coarse cell strattles GM
          sum = sum + (fine(i,j) - cmean(iclon,jclat))**2
          num  = num + 1 
        enddo
      enddo

      if (num .gt. 0) then
        cvar(iclon,jclat) = sum/num
      else
        cvar(iclon,jclat) = 1.e30
      endif
    end do
  end do
  return
end subroutine varf2c
