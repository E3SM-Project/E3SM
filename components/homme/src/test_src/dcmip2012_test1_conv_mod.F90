module dcmip2012_test1_conv_mod

  ! Based on DCMIP 2012 tests 1-1,2,3.

  use parallel_mod,       only: abortmp
  ! Use physical constants consistent with HOMME
  use physical_constants, only: a => rearth0, Rd => Rgas, g, cp, pi => dd_pi, p0

  implicit none
  private

  integer, parameter :: rt = 8

  real(rt), parameter :: &
       tau     = 12.d0 * 86400.d0, & ! period of motion 12 days
       T0      = 300.d0,           & ! temperature (K)
       ztop    = 12000.d0,         & ! model top (m)
       H       = Rd * T0 / g         ! scale height

  ! For tracers.
  real(rt), parameter :: qlon1 = 5.d0*(pi/6.d0), qlat1 = 0, &
       &                 qlon2 = -qlon1, qlat2 = 0, &
       &                 qsc_min = 0.1d0

  real(rt), parameter :: zero = 0.d0, one = 1.d0, two = 2.d0

  public :: test1_conv_advection, test1_conv_print_results

contains

  subroutine get_nondiv2d_uv(time, lon, lat, u, v)
    ! Classic 2D nondivergent flow field.

    real(rt), intent(in ) :: time, lon, lat
    real(rt), intent(out) :: u, v

    real(rt), parameter :: &
         u0      = (2.d0*pi*a)/tau,    &  ! 2 pi a / 12 days
         k0      = (10.d0*a)/tau          ! velocity magnitude

    real(rt) :: lonp

    ! translational longitude
    lonp = lon - 2.d0*pi*time/tau
    ! zonal velocity
    u = k0*sin(lonp)*sin(lonp)*sin(2.d0*lat)*cos(pi*time/tau) + u0*cos(lat)
    ! meridional velocity
    v = k0*sin(2.d0*lonp)*cos(lat)*cos(pi*time/tau)
  end subroutine get_nondiv2d_uv

  subroutine get_nondiv3d_uv(bs, pbot, ptop, zbot, ztop, ztaper, time, lon, lat, p, z, u, v, w)
    real(rt), intent(in ) :: bs, pbot, ptop, zbot, ztop, ztaper, time, lon, lat, p, z
    real(rt), intent(out) :: u, v, w

    real(rt), parameter :: omega0  = (2*23000.d0*pi)/tau

    real(rt) :: s, s_p, lonp, ud, c, arg
    
    ! This is essentially the test 1-1 flow. The key difference in this flow
    ! removes the factor of 2 in ud and w cos-time factors. The 2 in the
    ! original code makes trajectories not return to their initial points.

    ! Shape function in p.
    if (p >= pbot .or. p <= ptop) then
       s = 0
       s_p = 0
    else
       c = 0.3d0
       arg = pi*(p - ptop)/(pbot - ptop)
       s = c*sin(arg)**3
       s_p = (3*c*pi/(pbot - ptop))*sin(arg)**2*cos(arg)
    end if
    ! Translational longitude.
    lonp = lon - 2.d0*pi*time/tau
    ! Nondivergent 2D flow.
    call get_nondiv2d_uv(time, lon, lat, u, v)
    ! Taper the 2D nondiv (u,v) flow in the z direction. This does not induce
    ! any w, and the 2D field remains nondivergent at each z.
    u = u*ztaper
    v = v*ztaper
    ! Divergent flow.
    ud = (omega0*a)*cos(lonp)*(cos(lat)**2.0)*cos(pi*time/tau)*s_p
    u = u + ud
    w = -((Rd*T0)/(g*p))*omega0*sin(lonp)*cos(lat)*cos(pi*time/tau)*s
  end subroutine get_nondiv3d_uv

  function get_2d_cinf_tracer(lon, lat) result(q)
    real(rt), intent(in) :: lon, lat

    real(rt) :: q

    real(rt) :: x, y, zeta

    x = cos(lat)*cos(lon)
    y = cos(lat)*sin(lon)
    zeta = sin(lat)
    q = 1.5d0*(1 + sin(pi*x)*sin(pi*y)*sin(pi*zeta))
  end function get_2d_cinf_tracer

  subroutine ll2xyz(lon, lat, x, y, z)
    ! Unit sphere.

    real(rt), intent(in) :: lon, lat
    real(rt), intent(out) :: x, y, z

    real(rt) :: sinl, cosl

    sinl = sin(lat)
    cosl = cos(lat)
    x = cos(lon)*cosl
    y = sin(lon)*cosl
    z = sinl
  end subroutine ll2xyz

  function great_circle_dist(lon1, lat1, lon2, lat2) result(d)
    ! Unit sphere.
    
    real(rt), intent(in) :: lon1, lat1, lon2, lat2
    real(rt) :: d
    
    real(rt) xA, yA, zA, xB, yB, zB, cp1, cp2, cp3, cpnorm, dotprod
    
    call ll2xyz(lon1, lat1, xA, yA, zA)
    call ll2xyz(lon2, lat2, xB, yB, zB)
    cp1 = yA*zB - yB*zA
    cp2 = xB*zA - xA*zB
    cp3 = xA*yB - xB*yA
    cpnorm = sqrt(cp1*cp1 + cp2*cp2 + cp3*cp3)
    dotprod = xA*xB + yA*yB + zA*zB
    d = atan2(cpnorm, dotprod)
  end function great_circle_dist

  function q_gh(x, y, z, xi, yi, zi) result(q)
    real(rt), intent(in) :: x, y, z, xi, yi, zi
    real(rt) :: q

    real(rt), parameter :: h_max = 0.95d0, b = 5.d0
    real(rt) :: r2
    
    r2 = (x - xi)**2 + (y - yi)**2 + (z - zi)**2
    q = h_max*exp(-b*r2)
  end function q_gh

  function q_cb(r, ri) result(q)
    real(rt), intent(in) :: r, ri
    real(rt) :: q

    real(rt), parameter :: h_max = one
    
    q = 0.5d0*h_max*(1 + cos(pi*ri/r))
  end function q_cb

  function q_sc(clon_in, clat, lon, lat, up_slot) result(q)
    real(rt), intent(in) :: clon_in, clat, lon, lat
    logical, intent(in) :: up_slot
    real(rt) :: q

    real(rt), parameter :: b = qsc_min, c = one, r = 0.5d0, &
         &                 lon_thr = r/6.d0, lat_thr = 5*(r/12.d0)

    real(rt) :: clon, ri

    clon = clon_in
    if (clon < zero) clon = clon + two*pi

    ri = great_circle_dist(lon, lat, clon, clat)
    q = b
    if (ri <= r) then
       if (abs(lon - clon) >= lon_thr) then
          q = c
          return
       else
          if (up_slot) then
             if (lat - clat < -lat_thr) then
                q = c
                return
             end if
          else
             if (lat - clat >  lat_thr) then
                q = c
                return
             end if
          end if
       end if
    end if
  end function q_sc

  function get_2d_gaussian_hills(lon, lat) result(q)
    real(rt), intent(in) :: lon, lat
    real(rt) :: q

    real(rt) :: x1, y1, z1, x2, y2, z2, x, y, z

    call ll2xyz(qlon1, qlat1, x1, y1, z1)
    call ll2xyz(qlon2, qlat2, x2, y2, z2)
    call ll2xyz(lon, lat, x, y, z)
    q = q_gh(x, y, z, x1, y1, z1) + q_gh(x, y, z, x2, y2, z2)
  end function get_2d_gaussian_hills

  function get_2d_cosine_bells(lon, lat) result(q)
    real(rt), intent(in) :: lon, lat
    real(rt) :: q

    real(rt), parameter :: r = 0.5d0, b = 0.1d0, c = 0.9d0
    real(rt) :: h, ri

    h = 0
    ri = great_circle_dist(lon, lat, qlon1, qlat1)
    if (ri < r) then
       h = q_cb(r, ri)
    else
       ri = great_circle_dist(lon, lat, qlon2, qlat2)
       if (ri < r) h = q_cb(r, ri)
    end if
    q = b + c*h
  end function get_2d_cosine_bells

  function get_2d_correlated_cosine_bells(lon, lat) result(q)
    real(rt), intent(in) :: lon, lat
    real(rt) :: q

    real(rt), parameter :: a = -0.8d0, b = 0.9d0

    q = get_2d_cosine_bells(lon, lat)
    q = a*q + b
  end function get_2d_correlated_cosine_bells

  function get_2d_slotted_cylinders(lon, lat) result(q)
    real(rt), intent(in) :: lon, lat
    real(rt) :: q

    q = q_sc(qlon1, qlat1, lon, lat, .true.)
    if (q < 0.5d0) q = q_sc(qlon2, qlat2, lon, lat, .false.)
  end function get_2d_slotted_cylinders

  subroutine test1_conv_advection_orography( &
       test_minor,time,lon,lat,p,z,zcoords,cfv,hybrid_eta,hya,hyb,u,v,w,t,phis,ps,rho,q1,q2,q3,q4)

    character(len=1), intent(in) :: test_minor ! a, b, c, d, or e
    real(rt), intent(in)  :: time         ! simulation time (s)
    real(rt), intent(in)  :: lon          ! Longitude (radians)
    real(rt), intent(in)  :: lat          ! Latitude (radians)
    real(rt), intent(in)  :: hya          ! A coefficient for hybrid-eta coordinate
    real(rt), intent(in)  :: hyb          ! B coefficient for hybrid-eta coordinate

    logical, intent(in)  :: hybrid_eta    ! flag to indicate whether the hybrid sigma-p (eta) coordinate is used

    real(rt), intent(out)  :: p           ! Pressure  (Pa)
    real(rt), intent(out)  :: z           ! Height (m)

    integer , intent(in)     :: zcoords   ! 0 or 1 see below
    integer , intent(in)     :: cfv       ! 0, 1 or 2 see below
    real(rt), intent(out)    :: u         ! Zonal wind (m s^-1)
    real(rt), intent(out)    :: v         ! Meridional wind (m s^-1)
    real(rt), intent(out)    :: w         ! Vertical Velocity (m s^-1)
    real(rt), intent(out)    :: t         ! Temperature (K)
    real(rt), intent(out)    :: phis      ! Surface Geopotential (m^2 s^-2)
    real(rt), intent(out)    :: ps        ! Surface Pressure (Pa)
    real(rt), intent(out)    :: rho       ! density (kg m^-3)
    real(rt), intent(out)    :: q1        ! Tracer q1 (kg/kg)
    real(rt), intent(out)    :: q2        ! Tracer q2 (kg/kg)
    real(rt), intent(out)    :: q3        ! Tracer q3 (kg/kg)
    real(rt), intent(out)    :: q4        ! Tracer q4 (kg/kg)

    real(rt), parameter :: &
         u0      = 2.d0*pi*a/tau,       &  ! Velocity Magnitude (m/s)
         alpha   = pi/6.d0,             &  ! rotation angle (radians), 30 degrees
         lambdam = 3.d0*pi/2.d0,        &  ! mountain longitude center point (radians)
         phim    = zero,                &  ! mountain latitude center point (radians)
         h0      = 2000.d0,             &  ! peak height of the mountain range (m)
         Rm      = 3.d0*pi/4.d0,        &  ! mountain radius (radians)
         ztop_t  = 2000.d0,             &  ! transition layer
         zbot_q  = ztop_t + 500.d0,     &  ! bottom of tracers; below, all q = 0
         lon_offset = 0.5d0*pi,         &  ! longitudinal translation of std 2d test flow and qs
         ! For Hadley-like flow. Multiply w and tracer vertical extent by (ztop
         ! - ztop_t)/ztop to compensate for smaller domain.
         tau_h   = 86400.d0,            &  ! period of motion 1 day (in s)
         z1_h    = ztop_t + 1000.d0,    &  ! position of lower tracer bound (m)
         z2_h    = z1_h + 6000.d0,      &  ! position of upper tracer bound (m)
         z0_h    = 0.5d0*(z1_h+z2_h),   &  ! midpoint (m)
         u0_h    = 250.d0,              &  ! Zonal velocity magnitude (m/s)
         ! w0_h is the main parameter to modify to make the test easier (smaller
         ! w0_h) or harder (larger).
         w0_h    = 0.05d0,              &  ! Vertical velocity magnitude (m/s)
         ! For 3D deformational flow.
         bs_a    = 1.0d0                   ! shape function smoothness

    real(rt) :: r, height, zs, zetam, ztaper, rho0, z_q_shape, ptop, ptop_t, &
         &      c0, fl, fl_lat, gz, gz_z, fz, fz_z, delta, lambdam_t, u_topo_fac, &
         &      u0_topo, tau_topo
    logical :: ps_timedep

    if (cfv /= 0)         call abortmp('test1_conv_advection_orography does not support cfv != 0')
    if (.not. hybrid_eta) call abortmp('test1_conv_advection_orography does not support !hybrid_eta')
    if (zcoords /= 0)     call abortmp('test1_conv_advection_orography does not support zcoords != 0')

    ! Mountain oscillation half-width (radians).
    zetam = pi/14.d0
    ! Smooth mountains for very less resource-intensive convergence testing.
    if (test_minor == 'c') zetam = pi/2.d0
    ! Smoother than default but still fairly rough.
    if (test_minor == 'd' .or. test_minor == 'f') zetam = pi/6.d0

    ps_timedep = test_minor == 'e' .or. test_minor == 'f'
    lambdam_t = lambdam
    if (ps_timedep) then
       ! Move the topography to make ps depend on time.
       u0_topo = u0
       tau_topo = tau
       if (test_minor == 'e') then
          u0_topo = u0_h
          tau_topo = tau_h
       end if
       u_topo_fac = -u0_topo/two
       lambdam_t = lambdam_t + &
            &      sin(pi*time/tau_topo)*(tau_topo/pi)*u_topo_fac & ! integral of u at lat = 0
            &      /a ! to radians
    end if
    r = great_circle_dist(lambdam_t, phim, lon, lat)
    if (r .lt. Rm) then
       zs = (h0/2.d0)*(one+cos(pi*r/Rm))*cos(pi*r/zetam)**2.d0
    else
       zs = zero
    endif
    if (test_minor == 'a') zs = zero
    zs = -zs ! holes instead of mountains
    phis = g*zs
    ps = p0 * exp(-zs/H)

    p = hya*p0 + hyb*ps
    height = H * log(p0/p)
    z = height

    T = T0

    rho = p/(Rd*T)
    rho0 = p0/(Rd*T)

    if (z <= 0) then
       ztaper = 0
    elseif (z >= ztop_t) then
       ztaper = 1
    else
       ztaper = (1 + cos(pi*(1 + z/ztop_t)))/2
    end if

    w = zero

    select case(test_minor)
    case('z') ! currently unused
       ! Solid body rotation
       ! Zonal Velocity
       u = u0*(cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))
       ! Meridional Velocity
       v = -u0*(sin(lon)*sin(alpha))
       u = u*ztaper
       v = v*ztaper
    case('b')
       ! 2D nondiv flow in each layer.
       call get_nondiv2d_uv(time, lon + lon_offset, lat, u, v)
       u = u*ztaper
       v = v*ztaper
    case('a', 'c', 'd', 'f')
       ! 3D nondiv flow.
       ptop_t = p0*exp(-ztop_t/H)
       ptop = p0*exp(-ztop/H)
       call get_nondiv3d_uv(bs_a, ptop_t, ptop, ztop_t, ztop, ztaper, &
            &               time, lon + lon_offset, lat, p, z, u, v, w)
    case('e')
       ! Similar to Hadley-like flow but with more smoothness in derivatives.
       u = u0_h*cos(lat)*cos(pi*time/tau_h)*ztaper
       fl = cos(lat)**2
       fl_lat = -2*cos(lat)*sin(lat)
       if (z <= 0) then
          fz = 0
          fz_z = 0
       else
          gz = pi*z/ztop
          gz_z = pi/ztop
          fz = -sin(gz)**3
          fz_z = -3*sin(gz)**2*cos(gz)*gz_z
       end if
       c0 = w0_h*(rho0/rho)*cos(pi*time/tau_h)
       w =    c0*(cos(lat)*fl_lat - 2*sin(lat)*fl)*fz
       v = -a*c0*(cos(lat)*fl                    )*fz_z
    case default
       call abortmp('test1_conv_advection_orography: invalid case')
    end select

    if (ps_timedep) then
       ! Low-level solid-body rotational wind for consistency with the moving ps
       ! field.
       u = u + cos(pi*time/tau_topo)*u_topo_fac*(1 - ztaper)*cos(lat)
    end if

    if (time > 0) then
       q1 = 0; q2 = 0; q3 = 0; q4 = 0
       return
    end if

    z_q_shape = 0.5d0*(1 - cos(2*pi*(z - zbot_q)/(ztop - zbot_q)))
    if (z < zbot_q .or. z > ztop) z_q_shape = zero

    select case(test_minor)
    case('e')
       if (height < z2_h .and. height > z1_h) then
          q1 = 0.5d0 * (one + cos(2.d0*pi*(z-z0_h)/(z2_h-z1_h)))
       else
          q1 = zero
       end if
       q2 = q1 * get_2d_cinf_tracer(lon, lat)
       q3 = q1 * get_2d_gaussian_hills(lon - lon_offset, lat)
       q4 = q1 * get_2d_cosine_bells(lon - lon_offset, lat)

    case default
       q1 = z_q_shape * get_2d_gaussian_hills(lon - lon_offset, lat)
       q2 = z_q_shape * get_2d_cosine_bells(lon - lon_offset, lat)
       q4 = z_q_shape * get_2d_correlated_cosine_bells(lon - lon_offset, lat)
       ! Tracer discontinuous in 3D.
       q3 = qsc_min
       delta = z2_h - z1_h
       if ( (z >= z1_h                .and. z <= z1_h + 0.25d0*delta) .or. &
            (z >= z1_h + 0.4d0 *delta .and. z <= z2_h - 0.4d0 *delta) .or. &
            (z <= z2_h                .and. z >= z2_h - 0.25d0*delta)) then
          q3 = get_2d_slotted_cylinders(lon - lon_offset, lat)
       end if
    end select
  end subroutine test1_conv_advection_orography

  subroutine test1_conv_advection(test_case,time,lon,lat,hya,hyb,p,z,u,v,w,use_w,t,phis,ps,rho,q)
    character(len=*), intent(in) :: test_case  ! dcmip2012_test1_{3a-f}_conv
    real(rt), intent(in)     :: time       ! simulation time (s)
    real(rt), intent(in)     :: lon, lat   ! Longitude, latitude (radians)
    real(rt), intent(in)     :: hya, hyb   ! Hybrid a, b coefficients
    real(rt), intent(inout)  :: z          ! Height (m)
    real(rt), intent(inout)  :: p          ! Pressure  (Pa)
    real(rt), intent(out)    :: u          ! Zonal wind (m s^-1)
    real(rt), intent(out)    :: v          ! Meridional wind (m s^-1)
    real(rt), intent(out)    :: w          ! Vertical Velocity (m s^-1)
    logical , intent(out)    :: use_w      ! Should caller use w or instead div(u,v)?
    real(rt), intent(out)    :: T          ! Temperature (K)
    real(rt), intent(out)    :: phis       ! Surface Geopotential (m^2 s^-2)
    real(rt), intent(out)    :: ps         ! Surface Pressure (Pa)
    real(rt), intent(out)    :: rho        ! density (kg m^-3)
    real(rt), intent(out)    :: q(5)       ! Tracer q1 (kg/kg)

    integer, parameter :: cfv = 0, zcoords = 0
    logical, parameter :: use_eta = .true.

    character(len=1) :: test_major, test_minor

    test_major = test_case(17:17)
    if (test_major == '3') test_minor = test_case(18:18)

    use_w = .false.
    select case(test_major)
    case('3')
       call test1_conv_advection_orography( &
            test_minor,time,lon,lat,p,z,zcoords,cfv,use_eta,hya,hyb,u,v,w,t,phis,ps,rho, &
            q(1),q(2),q(3),q(4))
    end select
  end subroutine test1_conv_advection

  subroutine test1_conv_print_results(test_case, elem, tl, hvcoord, par, subnum)
    use element_mod, only: element_t
    use time_mod, only: timelevel_t
    use hybvcoord_mod, only: hvcoord_t
    use parallel_mod, only: parallel_t, pmax_1d
    use dimensions_mod, only: nelemd, nlev, qsize, np
    use parallel_mod, only: global_shared_buf, global_shared_sum
    use global_norms_mod, only: wrap_repro_sum
    use physical_constants, only: Rd => Rgas, p0

    character(len=*), intent(in) :: test_case
    type(element_t), intent(in) :: elem(:)
    type(timelevel_t), intent(in) :: tl
    type(hvcoord_t), intent(in) :: hvcoord
    type(parallel_t), intent(in) :: par
    integer, intent(in) :: subnum

    real(rt) :: q(np,np,5), lon, lat, z, p, phis, u, v, w, T, phis_ps, ps, rho, time, &
         hya, hyb, a, b, reldif, linf_num(qsize), linf_den(qsize)
    integer :: ie, k, iq, i, j
    logical :: use_w

    ! Set time to 0 to get the initial conditions.
    time = 0._rt

    linf_num = 0
    linf_den = 0
    do ie = 1,nelemd
       global_shared_buf(ie,:2*qsize) = 0._rt
       do k = 1,nlev
          ! test1_conv_advection_orography uses these:
          hya = hvcoord%hyam(k)
          hyb = hvcoord%hybm(k)
          ! test1_advection_deformation uses these, in which ps = p0:
          p = p0 * hvcoord%etam(k)
          z = H * log(1.0d0/hvcoord%etam(k))

          ! Normwise relative errors. We weight the horizontal direction by
          ! sphereme but do not weight the vertical direction; each vertical
          ! level in a column has equal weight.
          
          do j = 1,np
             do i = 1,np
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                select case(subnum)
                case (1)
                   call test1_conv_advection( &
                        test_case,time,lon,lat,hya,hyb,p,z,u,v,w,use_w,T,phis,ps,rho,q(i,j,:))
                end select
             end do
          end do
          
          do iq = 1,qsize
             global_shared_buf(ie,2*iq-1) = global_shared_buf(ie,2*iq-1) + &
                  sum(elem(ie)%spheremp*(elem(ie)%state%Q(:,:,k,iq) - q(:,:,iq))**2)
             global_shared_buf(ie,2*iq) = global_shared_buf(ie,2*iq) + &
                  sum(elem(ie)%spheremp*q(:,:,iq)**2)
             linf_num(iq) = max(linf_num(iq), &
                  maxval(abs(elem(ie)%state%Q(:,:,k,iq) - q(:,:,iq))))
             linf_den(iq) = max(linf_den(iq), &
                  maxval(abs(q(:,:,iq))))
          end do
       end do
    end do

    call wrap_repro_sum(nvars=2*qsize, comm=par%comm)
    do iq = 1, qsize
       linf_num(iq) = pmax_1d(linf_num(iq:iq), par)
       linf_den(iq) = pmax_1d(linf_den(iq:iq), par)
    end do
    
    if (par%masterproc) then
       print '(a)', 'test1_conv>                          l2                    linf'
       do iq = 1,qsize
          a = global_shared_sum(2*iq-1)
          b = global_shared_sum(2*iq)
          reldif = sqrt(a/b)
          print '(a,i2,es24.16,es24.16)', 'test1_conv> Q', &
               iq, reldif, linf_num(iq)/linf_den(iq)
       end do
    end if
  end subroutine test1_conv_print_results

end module dcmip2012_test1_conv_mod
