module fvm_overlap_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8

  real (kind=r8),parameter, private  :: bignum = 1.0e20_r8
  real (kind=r8),parameter, private  :: tiny   = 1.0e-12_r8
  real (kind=r8),parameter, private  :: fuzzy_width = 10.0_r8*tiny

  public:: compute_weights_cell

  private
  integer, parameter :: max_cross = 10
contains
  subroutine compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
       xcell_in,ycell_in,jx,jy,nreconstruction,xgno,ygno,igno_min,igno_max,&
       jx_min, jx_max, jy_min, jy_max,&
       ngauss,gauss_weights,abscissae,weights,weights_eul_index,jcollect,jmax_segments)

    implicit none
    integer , intent(in) :: nvertex
    logical, intent(in) :: lexact_horizontal_line_integrals
    integer , intent(in):: nreconstruction, jx,jy,ngauss,jmax_segments
    !
    ! dimension(nvertex)
    !
    real (kind=r8)   ,  dimension(4), intent(in):: xcell_in,ycell_in
    !
    integer , intent(in)               :: jx_min, jy_min, jx_max, jy_max,igno_min,igno_max
    !
    ! dimension(-ihalo:nc+2+ihalo)
    !
    real (kind=r8), dimension(igno_min:igno_max), intent(in) :: xgno, ygno
    !
    ! for Gaussian quadrature
    !
    real (kind=r8), dimension(:), intent(in) :: gauss_weights, abscissae !dimension(ngauss)
    !
    ! Number of Eulerian sub-cell integrals for the cell in question
    !
    integer , intent(out)       :: jcollect
    !
    ! local workspace
    !
    !
    ! max number of line segments is:
    !
    ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
    !
    real (kind=r8)   ,  &
         dimension(jmax_segments,nreconstruction), intent(out) :: weights
    integer ,  &
         dimension(jmax_segments,2), intent(out)      :: weights_eul_index

    integer :: jsegment
    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer :: jcross_lat
    !
    ! max. crossings per side is 2*ihalo
    !
    real (kind=r8), &
         dimension(max_cross,2) :: r_cross_lat
    integer , &
         dimension(max_cross,2) :: cross_lat_eul_index
    real (kind=r8)   ,  dimension(nvertex) :: xcell,ycell

    xcell = xcell_in(1:nvertex)
    ycell = ycell_in(1:nvertex)

    jsegment   = 0
    weights    = 0.0_r8
    jcross_lat = 0

    call side_integral(lexact_horizontal_line_integrals,xcell,ycell,nvertex,jsegment,jmax_segments,&
         weights,weights_eul_index,nreconstruction,jx,jy,xgno,ygno,igno_min,igno_max,jx_min, jx_max, jy_min, jy_max,&
         ngauss,gauss_weights,abscissae,&
         jcross_lat,r_cross_lat,cross_lat_eul_index)
    !
    !**********************
    !
    ! Do inner integrals
    !
    !**********************
    !
    call compute_inner_line_integrals_lat(lexact_horizontal_line_integrals,&
         r_cross_lat,cross_lat_eul_index,&
         jcross_lat,jsegment,xgno,igno_min,igno_max,jx_min, jx_max, jy_min, jy_max,&
         weights,weights_eul_index,&
         nreconstruction,ngauss,gauss_weights,abscissae)

    IF (ABS((jcross_lat/2)-DBLE(jcross_lat)/2.0_r8)>tiny) then
      WRITE(*,*) "number of latitude crossings are not even: ABORT",jcross_lat,jx,jy
      STOP
    END IF

    !
    ! collect line-segment that reside in the same Eulerian cell
    !
    if (jsegment>0) then
      call collect(weights,weights_eul_index,nreconstruction,jcollect,jsegment,jmax_segments)
    else
      jcollect = 0
    end if
  end subroutine compute_weights_cell
  !
  !****************************************************************************
  !
  ! organize data and store it
  !
  !****************************************************************************
  !
  subroutine collect(weights,weights_eul_index,nreconstruction,jcollect,jsegment,jmax_segments)
    implicit none
    integer ,                                  INTENT(IN   ) :: jsegment,jmax_segments
    integer , intent(in)    :: nreconstruction
    !
    real (kind=r8)  , dimension(:,:), intent(inout) :: weights !dimension(jmax_segments,nreconstruction)
    integer , dimension(:,:), intent(inout) :: weights_eul_index !dimension(jmax_segments,2)
    integer ,                                  INTENT(OUT  ) :: jcollect
    !
    ! local workspace
    !
    integer :: imin, imax, jmin, jmax, i,j,k,h
    logical                 :: ltmp

    real (kind=r8)  , dimension(jmax_segments,nreconstruction) :: weights_out
    integer , dimension(jmax_segments,2     ) :: weights_eul_index_out

    weights_out           = 0.0_r8
    weights_eul_index_out = -100

    imin = MINVAL(weights_eul_index(1:jsegment,1))
    imax = MAXVAL(weights_eul_index(1:jsegment,1))
    jmin = MINVAL(weights_eul_index(1:jsegment,2))
    jmax = MAXVAL(weights_eul_index(1:jsegment,2))

    ltmp = .FALSE.

    jcollect = 1

    do j=jmin,jmax
       do i=imin,imax
          do k=1,jsegment
             if (weights_eul_index(k,1)==i.AND.weights_eul_index(k,2)==j) then
                weights_out(jcollect,1:nreconstruction) = &
                weights_out(jcollect,1:nreconstruction) + weights(k,1:nreconstruction)
                ltmp = .TRUE.
                h = k
             endif
          enddo
          if (ltmp) then
             weights_eul_index_out(jcollect,:) = weights_eul_index(h,:)
             jcollect = jcollect+1
          endif
          ltmp = .FALSE.
       enddo
    enddo
    jcollect = jcollect-1
    weights           = weights_out
    weights_eul_index = weights_eul_index_out
  end subroutine collect
  !
  !*****************************************************************************************
  !
  ! compute crossings with Eulerian latitudes and longitudes
  !
  !*****************************************************************************************
  !
  subroutine compute_inner_line_integrals_lat(lexact_horizontal_line_integrals,r_cross_lat,&
       cross_lat_eul_index,&
       jcross_lat,jsegment,xgno,igno_min,igno_max,jx_min,jx_max,jy_min, jy_max,weights,weights_eul_index,&
       nreconstruction,ngauss,gauss_weights,abscissae)
    implicit none
    logical, intent(in) :: lexact_horizontal_line_integrals
    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer ,         intent(in):: jcross_lat, nreconstruction,ngauss,igno_min,igno_max
    integer ,         intent(inout):: jsegment
    !
    ! for Gaussian quadrature
    !
    real (kind=r8), dimension(ngauss), intent(in) :: gauss_weights, abscissae
    !
    ! max. crossings per side is 2*ihalo
    !

    real (kind=r8)  , dimension(:,:), intent(in):: r_cross_lat ! dimension(8*ihalo,2)
    integer , dimension(:,:), intent(in):: cross_lat_eul_index ! ! dimension(8*ihalo,2)
    integer , intent(in):: jx_min, jx_max, jy_min, jy_max

    real (kind=r8), dimension(igno_min:igno_max), intent(in)  :: xgno !dimension(-ihalo:nc+2+ihalo)
    !
    ! dimension(jmax_segments,nreconstruction)
    !
    real (kind=r8), dimension(:,:), intent(inout) :: weights
    !
    ! dimension(jmax_segments,2)
    !
    integer , dimension(:,:), intent(inout)               :: weights_eul_index

    real (kind=r8)   , dimension(nreconstruction)         :: weights_tmp
    integer :: imin,imax,i,j,k,h
    real (kind=r8), dimension(2)  :: rstart,rend,rend_tmp
    real (kind=r8), dimension(2)  :: xseg, yseg
    5   FORMAT(10e14.6)
    if (jcross_lat>0) then
       do i=MINVAL(cross_lat_eul_index(1:jcross_lat,2)),MAXVAL(cross_lat_eul_index(1:jcross_lat,2))
          !
          ! find "first" crossing with Eulerian cell i
          !
          do k=1,jcross_lat
             if (cross_lat_eul_index(k,2)==i) exit
          enddo
          do j=k+1,jcross_lat
             !
             ! find "second" crossing with Eulerian cell i
             !
             if (cross_lat_eul_index(j,2)==i) then
                if (r_cross_lat(k,1)<r_cross_lat(j,1)) then
                   rstart = r_cross_lat(k,1:2)
                   rend   = r_cross_lat(j,1:2)
                   imin   = cross_lat_eul_index(k,1)
                   imax   = cross_lat_eul_index(j,1)
                else
                   rstart = r_cross_lat(j,1:2)
                   rend   = r_cross_lat(k,1:2)
                   imin   = cross_lat_eul_index(j,1)
                   imax   = cross_lat_eul_index(k,1)
                endif
                do h=imin,imax
                   if (h==imax) then
                      rend_tmp = rend
                   else
                      rend_tmp(1) = xgno(h+1)
                      rend_tmp(2) = r_cross_lat(k,2)
                   endif
                   xseg(1) = rstart(1)
                   xseg(2) = rend_tmp(1)
                   yseg(1) = rstart(2)
                   yseg(2) = rend_tmp(2)
                   call get_weights_exact(lexact_horizontal_line_integrals, weights_tmp,xseg,yseg,&
                        nreconstruction,ngauss,gauss_weights,abscissae)

                   if (i.LE.jy_max-1.AND.i.GE.jy_min.AND.h.LE.jx_max-1.AND.h.GE.jx_min) then
                      jsegment=jsegment+1
                      weights_eul_index(jsegment,1) = h
                      weights_eul_index(jsegment,2) = i
                      weights(jsegment,1:nreconstruction) = -weights_tmp
                   endif
                   !
                   ! subtract the same weights on the west side of the line
                   !
                   if (i.LE.jy_max.AND.i.GE.jy_min+1.AND.h.LE.jx_max-1.AND.h.GE.jx_min) then
                      jsegment = jsegment+1
                      weights_eul_index(jsegment,1) = h
                      weights_eul_index(jsegment,2) = i-1
                      weights(jsegment,1:nreconstruction) = weights_tmp
                   endif
                   !
                   ! prepare for next iteration
                   !
                   rstart = rend_tmp
                enddo
             endif
          enddo
       enddo
    endif
  end subroutine compute_inner_line_integrals_lat

  !
  ! line integral from (a1_in,a2_in) to (b1_in,b2_in)
  ! If line is coniciding with an Eulerian longitude or latitude the routine
  ! needs to know where an adjacent side is located to determine which
  ! reconstruction must be used. therefore (c1,c2) is passed to the routine
  !
  !

  subroutine side_integral(lexact_horizontal_line_integrals,&
       x_in,y_in,nvertex,jsegment,jmax_segments,&
       weights,weights_eul_index,nreconstruction,jx,jy,xgno,ygno,igno_min,igno_max,&
       jx_min,jx_max,jy_min,jy_max,&
       ngauss,gauss_weights,abscissae,&!)!phl add jx_min etc.
       jcross_lat,r_cross_lat,cross_lat_eul_index)
    implicit none


    logical, intent(in) :: lexact_horizontal_line_integrals
    integer ,            intent(in)    :: nreconstruction,jx,jy,jmax_segments,ngauss
    integer , intent(in)               :: nvertex,igno_min,igno_max
    !
    ! for Gaussian quadrature
    !
    real (kind=r8), dimension(:), intent(in) :: gauss_weights, abscissae !dimension(ngauss)
    real (kind=r8), dimension(:)        , intent(in)    :: x_in,y_in !dimension(1:nvertex)

    integer , intent(in)               :: jx_min, jy_min, jx_max, jy_max
    real (kind=r8), dimension(igno_min:igno_max), intent(in) :: xgno, ygno !dimension(-ihalo:nc+2+ihalo)
    integer ,            intent(inout) :: jsegment
!    integer ,dimension(0:2),intent(in)    :: jx_eul_in, jy_eul_in
    real (kind=r8)   , dimension(:,:), intent(out) :: weights !dimension(jmax_segments,nreconstruction)
    integer ,  &
         dimension(jmax_segments,2), intent(out) :: weights_eul_index

    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer ,         intent(inout):: jcross_lat
    !
    ! max. crossings per side is 2*ihalo
    !
    real (kind=r8), &
         dimension(max_cross,2), intent(inout):: r_cross_lat
    integer , &
         dimension(max_cross,2), intent(inout):: cross_lat_eul_index
    !
    ! local variables
    !
    real (kind=r8), dimension(2) :: xseg,yseg
    real (kind=r8), dimension(0:3) :: x,y
    real (kind=r8)               :: xeul,yeul,xcross,ycross,slope
    integer ::    jx_eul_tmp,jy_eul_tmp
    integer :: xsgn1,ysgn1,xsgn2,ysgn2
    integer :: iter
    logical :: lcontinue, lsame_cell_x, lsame_cell_y

    integer :: jx_eul, jy_eul, side_count
    real (kind=r8), dimension(0:nvertex+2)  :: xcell,ycell


5   FORMAT(10e14.6)
    !
    !***********************************************
    !
    ! find jx_eul and jy_eul for (x(1),y(1))
    !
    !***********************************************
    !
    jx_eul = jx; jy_eul = jy
    xcell(1:nvertex)=x_in; ycell(1:nvertex)=y_in
    DO iter=1,nvertex
      CALL truncate_vertex(xcell(iter),jx_eul,xgno,igno_min,igno_max)
      CALL truncate_vertex(ycell(iter),jy_eul,ygno,igno_min,igno_max)
    END DO
    xcell(0) = xcell(nvertex); xcell(nvertex+1)=xcell(1); xcell(nvertex+2)=xcell(2);
    ycell(0) = ycell(nvertex); ycell(nvertex+1)=ycell(1); ycell(nvertex+2)=ycell(2);


    IF ((&
         MAXVAL(xcell).LE.xgno(jx_min).OR.MINVAL(xcell).GE.xgno(jx_max).OR.&
         MAXVAL(ycell).LE.ygno(jy_min).OR.MINVAL(ycell).GE.ygno(jy_max))) THEN
      !
      ! entire cell off panel
      !
    ELSE
      jx_eul = MIN(MAX(jx,jx_min),jx_max-1)
      jy_eul = MIN(MAX(jy,jy_min),jy_max-1)
      CALL which_eul_cell(xcell(1:3),jx_eul,xgno,igno_min,igno_max)
      CALL which_eul_cell(ycell(1:3),jy_eul,ygno,igno_min,igno_max)

      side_count = 1
      DO WHILE (side_count<nvertex+1)
        iter = 0
        lcontinue = .TRUE.
        x(0:3) = xcell(side_count-1:side_count+2); y(0:3) = ycell(side_count-1:side_count+2);
        DO while (lcontinue)
          iter = iter+1
          IF (iter>10) THEN
            WRITE(*,*) "search not converging",iter
            STOP
          END IF
          lsame_cell_x = (x(2).GE.xgno(jx_eul).AND.x(2).LE.xgno(jx_eul+1))
          lsame_cell_y = (y(2).GE.ygno(jy_eul).AND.y(2).LE.ygno(jy_eul+1))
          IF (lsame_cell_x.AND.lsame_cell_y) THEN
            !
            !****************************
            !
            ! same cell integral
            !
            !****************************
            !
            xseg(1) = x(1); yseg(1) = y(1); xseg(2) = x(2); yseg(2) = y(2)
            jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul;
            lcontinue = .FALSE.
            !
            ! prepare for next side if (x(2),y(2)) is on a grid line
            !
            IF (x(2).EQ.xgno(jx_eul+1).AND.x(3)>xgno(jx_eul+1)) THEN
              !
              ! cross longitude jx_eul+1
              !
              jx_eul=jx_eul+1
            ELSE IF (x(2).EQ.xgno(jx_eul  ).AND.x(3)<xgno(jx_eul)) THEN
              !
              ! cross longitude jx_eul
              !
              jx_eul=jx_eul-1
            END IF
            IF (y(2).EQ.ygno(jy_eul+1).AND.y(3)>ygno(jy_eul+1)) THEN
              !
              ! register crossing with latitude: line-segments point Northward
              !
              jcross_lat = jcross_lat + 1
              jy_eul     = jy_eul     + 1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              cross_lat_eul_index(jcross_lat,2) = jy_eul
              r_cross_lat(jcross_lat,1) = x(2)
              r_cross_lat(jcross_lat,2) = y(2)
!              write(*,*) "A register crossing with latitude",x(2),y(2),jx_eul,jy_eul
            ELSE IF (y(2).EQ.ygno(jy_eul  ).AND.y(3)<ygno(jy_eul)) THEN
              !
              ! register crossing with latitude: line-segments point Southward
              !
              jcross_lat = jcross_lat+1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              cross_lat_eul_index(jcross_lat,2) = jy_eul
              r_cross_lat(jcross_lat,1) = x(2)
              r_cross_lat(jcross_lat,2) = y(2)
!              write(*,*) "B register crossing with latitude",x(2),y(2),jx_eul,jy_eul
              !
              jy_eul=jy_eul-1
            END IF
            lcontinue=.FALSE.
          ELSE
            !
            !****************************
            !
            ! not same cell integral
            !
            !****************************
            !
            IF (lsame_cell_x) THEN
              ysgn1 = (1+INT(SIGN(1.0_r8,y(2)-y(1))))/2 !"1" if y(2)>y(1) else "0"
              ysgn2 = INT(SIGN(1.0_r8,y(2)-y(1)))       !"1" if y(2)>y(1) else "-1"
              !
              !*******************************************************************************
              !
              ! there is at least one crossing with latitudes but no crossing with longitudes
              !
              !*******************************************************************************
              !
              yeul   = ygno(jy_eul+ysgn1)
              IF (x(1).EQ.x(2)) THEN
                !
                ! line segment is parallel to longitude (infinite slope)
                !
                xcross = x(1)
              ELSE
                slope  = (y(2)-y(1))/(x(2)-x(1))
                xcross = x_cross_eul_lat(x(1),y(1),yeul,slope)
                !
                ! constrain crossing to be "physically" possible
                !
                xcross = MIN(MAX(xcross,xgno(jx_eul)),xgno(jx_eul+1))
                !
                ! debugging
                !
                IF (xcross.GT.xgno(jx_eul+1).OR.xcross.LT.xgno(jx_eul)) THEN
                  WRITE(*,*) "xcross is out of range",jx,jy
                  WRITE(*,*) "xcross-xgno(jx_eul+1), xcross-xgno(jx_eul))",&
                       xcross-xgno(jx_eul+1), xcross-ygno(jx_eul)
                  STOP
                END IF
              END IF
              xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xcross; yseg(2) = yeul
              jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul;
              !
              ! prepare for next iteration
              !
              x(0) = x(1); y(0) = y(1); x(1) = xcross; y(1) = yeul; jy_eul = jy_eul+ysgn2
              !
              ! register crossing with latitude
              !
              jcross_lat = jcross_lat+1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              if (ysgn2>0) then
                cross_lat_eul_index(jcross_lat,2) = jy_eul
              else
                cross_lat_eul_index(jcross_lat,2) = jy_eul+1
              end if
              r_cross_lat(jcross_lat,1) = xcross
              r_cross_lat(jcross_lat,2) = yeul
            ELSE IF (lsame_cell_y) THEN
              !
              !*******************************************************************************
              !
              ! there is at least one crossing with longitudes but no crossing with latitudes
              !
              !*******************************************************************************
              !
              xsgn1 = (1+INT(SIGN(1.0_r8,x(2)-x(1))))/2 !"1" if x(2)>x(1) else "0"
              xsgn2 = INT(SIGN(1.0_r8,x(2)-x(1))) !"1" if x(2)>x(1) else "-1"
              xeul   = xgno(jx_eul+xsgn1)
              IF (ABS(x(2)-x(1))<fuzzy_width) THEN
                ! fuzzy crossing
                ycross = 0.5_r8*(y(2)-y(1))
              ELSE
                slope  = (y(2)-y(1))/(x(2)-x(1))
                ycross = y_cross_eul_lon(x(1),y(1),xeul,slope)
              END IF
              !
              ! constrain crossing to be "physically" possible
              !
              ycross = MIN(MAX(ycross,ygno(jy_eul)),ygno(jy_eul+1))
              !
              ! debugging
              !
              IF (ycross.GT.ygno(jy_eul+1).OR.ycross.LT.ygno(jy_eul)) THEN
                WRITE(*,*) "ycross is out of range"
                WRITE(*,*) "jx,jy,jx_eul,jy_eul",jx,jy,jx_eul,jy_eul
                WRITE(*,*) "ycross-ygno(jy_eul+1), ycross-ygno(jy_eul))",&
                     ycross-ygno(jy_eul+1), ycross-ygno(jy_eul)
                STOP
              END IF
              xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xeul; yseg(2) = ycross
              jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul;
              !
              ! prepare for next iteration
              !
              x(0) = x(1); y(0) = y(1); x(1) = xeul; y(1) = ycross; jx_eul = jx_eul+xsgn2
            ELSE
              !
              !*******************************************************************************
              !
              ! there are crossings with longitude(s) and latitude(s)
              !
              !*******************************************************************************
              !
              xsgn1 = (1+INT(SIGN(1.0_r8,x(2)-x(1))))/2 !"1" if x(2)>x(1) else "0"
              xsgn2 = (INT(SIGN(1.0_r8,x(2)-x(1)))) !"1" if x(2)>x(1) else "0"
              xeul   = xgno(jx_eul+xsgn1)
              ysgn1 = (1+INT(SIGN(1.0_r8,y(2)-y(1))))/2 !"1" if y(2)>y(1) else "0"
              ysgn2 = INT(SIGN(1.0_r8,y(2)-y(1)))       !"1" if y(2)>y(1) else "-1"
              yeul   = ygno(jy_eul+ysgn1)

              slope  = (y(2)-y(1))/(x(2)-x(1))
              IF (ABS(x(2)-x(1))<fuzzy_width) THEN
                ycross = 0.5_r8*(y(2)-y(1))
              ELSE
                ycross = y_cross_eul_lon(x(1),y(1),xeul,slope)
              END IF
              xcross = x_cross_eul_lat(x(1),y(1),yeul,slope)


              IF ((xsgn2>0.AND.xcross.LE.xeul).OR.(xsgn2<0.AND.xcross.GE.xeul)) THEN
                !
                ! cross latitude
                !
                xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xcross; yseg(2) = yeul
                jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul;
                !
                ! prepare for next iteration
                !
                x(0) = x(1); y(0) = y(1); x(1) = xcross; y(1) = yeul; jy_eul = jy_eul+ysgn2
                !
                ! register crossing with latitude
                !
                jcross_lat = jcross_lat+1
                cross_lat_eul_index(jcross_lat,1) = jx_eul
                if (ysgn2>0) then
                  cross_lat_eul_index(jcross_lat,2) = jy_eul
                else
                  cross_lat_eul_index(jcross_lat,2) = jy_eul+1
                end if
                r_cross_lat(jcross_lat,1) = xcross
                r_cross_lat(jcross_lat,2) = yeul
!              write(*,*) "D register crossing with latitude",xcross,yeul,jx_eul,cross_lat_eul_index(jcross_lat,2)
              ELSE
                !
                ! cross longitude
                !
                xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xeul; yseg(2) = ycross
                jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul;
                !
                ! prepare for next iteration
                !
                x(0) = x(1); y(0) = y(1); x(1) = xeul; y(1) = ycross; jx_eul = jx_eul+xsgn2
              END IF

            END IF
          END IF
          !
          ! register line-segment (don't register line-segment if outside of panel)
          !
          if (jx_eul_tmp>=jx_min.AND.jy_eul_tmp>=jy_min.AND.&
               jx_eul_tmp<=jx_max-1.AND.jy_eul_tmp<=jy_max-1) then
            jsegment=jsegment+1
            weights_eul_index(jsegment,1) = jx_eul_tmp
            weights_eul_index(jsegment,2) = jy_eul_tmp

            call get_weights_exact(lexact_horizontal_line_integrals.AND.ABS(yseg(2)-yseg(1))<tiny,&
                 weights(jsegment,:),&
                 xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
!old            call get_weights_gauss(weights(jsegment,1:nreconstruction),&
!old                 xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
          ELSE
            !
            ! segment outside of panel
            !
          END IF

        END DO
        side_count = side_count+1
      END DO
    END IF
  end subroutine side_integral


  real (kind=r8) function y_cross_eul_lon(x,y,xeul,slope)
    implicit none
    real (kind=r8), intent(in) :: x,y
    real (kind=r8)              , intent(in) :: xeul,slope
    !
    ! line: y=a*x+b
    !
    real (kind=r8) :: b

    b = y-slope*x
    y_cross_eul_lon = slope*xeul+b
  end function y_cross_eul_lon

  real (kind=r8) function x_cross_eul_lat(x,y,yeul,slope)
    implicit none
    real (kind=r8), intent(in) :: x,y
    real (kind=r8)              , intent(in) :: yeul,slope

    if (fuzzy(ABS(slope),fuzzy_width)>0) THEN
        x_cross_eul_lat = x+(yeul-y)/slope
    ELSE
      x_cross_eul_lat = bignum
    END IF
  end function x_cross_eul_lat

  subroutine get_weights_exact(lexact_horizontal_line_integrals,weights,xseg,yseg,nreconstruction,&
       ngauss,gauss_weights,abscissae)
    use fvm_analytic_mod, only: I_00, I_10, I_01, I_20, I_02, I_11
    implicit none
    logical, intent(in) :: lexact_horizontal_line_integrals
    integer , intent(in) :: nreconstruction, ngauss
    real (kind=r8), intent(out) :: weights(:)
    real (kind=r8), dimension(:), intent(in) :: gauss_weights, abscissae !dimension(ngauss)


    real (kind=r8), dimension(:), intent(in) :: xseg,yseg !dimension(2)
    !
    ! compute weights
    !
    if(lexact_horizontal_line_integrals) then
      weights(1) = ((I_00(xseg(2),yseg(2))-I_00(xseg(1),yseg(1))))
      if (ABS(weights(1))>1.0_r8) THEN
        WRITE(*,*) "1 exact weights(jsegment)",weights(1),xseg,yseg
        stop
      end if
      if (nreconstruction>1) then
         weights(2) = ((I_10(xseg(2),yseg(2))-I_10(xseg(1),yseg(1))))
         weights(3) = ((I_01(xseg(2),yseg(2))-I_01(xseg(1),yseg(1))))
      endif
      if (nreconstruction>3) then
         weights(4) = ((I_20(xseg(2),yseg(2))-I_20(xseg(1),yseg(1))))
         weights(5) = ((I_02(xseg(2),yseg(2))-I_02(xseg(1),yseg(1))))
         weights(6) = ((I_11(xseg(2),yseg(2))-I_11(xseg(1),yseg(1))))
      endif
    else
      call get_weights_gauss(weights,xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
    endif
  end subroutine get_weights_exact



  subroutine get_weights_gauss(weights,xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
    use fvm_analytic_mod, only: F_00, F_10, F_01, F_20, F_02, F_11
    implicit none
    integer , intent(in) :: nreconstruction,ngauss
    real (kind=r8), intent(out) :: weights(:)
    real (kind=r8), dimension(2     ), intent(in) :: xseg,yseg
    real (kind=r8) :: slope
    !
    ! compute weights
    !
    !
    ! for Gaussian quadrature
    !
    real (kind=r8), dimension(ngauss), intent(in) :: gauss_weights, abscissae

    ! if line-segment parallel to x or y use exact formulaes else use qudrature
    !
    real (kind=r8) :: b,integral,dx2,xc,x,y
    integer :: i

!    if (fuzzy(abs(xseg(1) -xseg(2)),fuzzy_width)==0)then
    if (xseg(1).EQ.xseg(2))then
      weights = 0.0_r8
    else
      slope    = (yseg(2)-yseg(1))/(xseg(2)-xseg(1))
      b        = yseg(1)-slope*xseg(1)
      dx2      = 0.5_r8*(xseg(2)-xseg(1))
      xc       = 0.5_r8*(xseg(1)+xseg(2))
      integral = 0.0_r8
      do i=1,ngauss
        x        = xc+abscissae(i)*dx2
        y        = slope*x+b
        integral = integral+gauss_weights(i)*F_00(x,y)
      enddo
      weights(1) = integral*dx2
      if (nreconstruction>1) then
        integral = 0.0_r8
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_10(x,y)
        enddo
        weights(2) = integral*dx2
        integral = 0.0_r8
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_01(x,y)
        enddo
        weights(3) = integral*dx2
      endif
      if (nreconstruction>3) then
        integral = 0.0_r8
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_20(x,y)
        enddo
        weights(4) = integral*dx2
        integral = 0.0_r8
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_02(x,y)
        enddo
        weights(5) = integral*dx2
        integral = 0.0_r8
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_11(x,y)
        enddo
        weights(6) = integral*dx2
      endif
    end if
  end subroutine get_weights_gauss

  subroutine truncate_vertex(x,j_eul,gno,igno_min,igno_max)
    implicit none
    integer , intent(inout) :: j_eul
    integer , intent(in)    :: igno_min,igno_max

    real (kind=r8)                    , intent(inout)    :: x
    real (kind=r8), dimension(igno_min:igno_max), intent(in)    :: gno !dimension(-ihalo:nc+2+ihalo)
!    real (kind=r8), intent(in)    :: eps

    logical                 :: lcontinue
    integer :: iter, xsgn
    real (kind=r8) :: dist,dist_new,tmp

    lcontinue = .TRUE.
    iter = 0
    dist = bignum

    xsgn     = INT(SIGN(1.0_r8,x-gno(j_eul)))

    DO WHILE (lcontinue)
      if ((j_eul<igno_min) .or. (j_eul>igno_max)) then
        write(*,*) 'something is wrong', j_eul,igno_min,igno_max, iter
        stop
      endif
      iter     = iter+1
      tmp      = x-gno(j_eul)
      dist_new = ABS(tmp)
      IF (dist_new>dist) THEN
        lcontinue = .FALSE.
      ELSE IF (ABS(tmp)<1.0E-9_r8) THEN
        x = gno(j_eul)
        lcontinue = .FALSE.
      ELSE
        j_eul = j_eul+xsgn
        dist = dist_new
      END IF
      IF (iter>100) THEN
        WRITE(*,*) "truncate vertex not converging"
        STOP
      END IF
    END DO
  END subroutine truncate_vertex

  subroutine which_eul_cell(x,j_eul,gno,igno_min,igno_max)
    implicit none
    integer , intent(inout) :: j_eul
    integer , intent(in)    :: igno_min,igno_max
    real (kind=r8), dimension(:)                    , intent(in):: x !dimension(3)
    real (kind=r8), dimension(igno_min:igno_max), intent(in)    :: gno ! dimension(-ihalo:nc+2+ihalo)

    logical :: lcontinue
    integer :: iter

    lcontinue = .TRUE.
    iter = 0

    DO WHILE (lcontinue)
      iter = iter+1
      IF (x(1).GE.gno(j_eul).AND.x(1).LT.gno(j_eul+1)) THEN
        lcontinue = .FALSE.
        !
        ! special case when x(1) is on top of grid line
        !
        IF (x(1).EQ.gno(j_eul)) THEN
          !
          ! x(1) is on top of gno(J_eul)
          !
          IF (x(2).GT.gno(j_eul)) THEN
            j_eul = j_eul
          ELSE IF (x(2).LT.gno(j_eul)) THEN
            j_eul = j_eul-1
          ELSE
            !
            ! x(2) is on gno(j_eul) grid line; need x(3) to determine Eulerian cell
            !
            IF (x(3).GT.gno(j_eul)) THEN
              !
              ! x(3) to the right
              !
              j_eul = j_eul
            ELSE IF (x(3).LT.gno(j_eul)) THEN
              !
              ! x(3) to the left
              !
              j_eul = j_eul-1
            ELSE
              WRITE(*,*) "inconsistent cell: x(1)=x(2)=x(3)",x(1),x(2),x(3)
              STOP
            END IF
          END IF
        END IF
      ELSE
        !
        ! searching - prepare for next iteration
        !
        IF (x(1).GE.gno(j_eul+1)) THEN
          j_eul = j_eul + 1
        ELSE
          !
          ! x(1).LT.gno(j_eul)
          !
          j_eul = j_eul - 1
        END IF
      END IF
      IF (iter>1000.OR.j_eul<igno_min.OR.j_eul>igno_max) THEN
        WRITE(*,*) "search is which_eul_cell not converging!", iter, j_eul,igno_min,igno_max
        WRITE(*,*) "gno", gno(igno_min), gno(igno_max)
        write(*,*) gno
        STOP
      END IF
    END DO
  END subroutine which_eul_cell


  function fuzzy(x,epsilon)
    implicit none

    integer :: fuzzy
    real (kind=r8), intent(in) :: epsilon
    real (kind=r8) :: x

    IF (ABS(x)<epsilon) THEN
      fuzzy = 0
    ELSE IF (x >epsilon) THEN
      fuzzy = 1
    ELSE !IF (x < fuzzy_width) THEN
      fuzzy = -1
    ENDIF
  end function

end module fvm_overlap_mod
