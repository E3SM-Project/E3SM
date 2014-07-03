
module sulchem

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Contains the sulfur cycle chemistry codes.
   ! 
   ! Author: 
   ! Chemistry code is from Mary Barth.
   ! Module coded by B. Eaton.
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols, pver
   use abortutils,   only: endrun
   use physconst,    only: tmelt
   use perf_mod
   use cam_logfile,  only: iulog

   implicit none
   save
   private
   public :: &
      inisulchem,   &! Initialize sulfur cycle chemistry.
      sulf_chemdr,  &
      cldychmini, dicor

   ! private module data

   integer, parameter ::&
      nz=18,              &
      nza=8,              &
      nalb=2,             &
      nz1=nz-1,           &
      nza1=nza-1,         &
      nz2=51,             &
      nsiz =nz *nza*nalb, &
      nsiz2=nz2*nza*nalb

   real(r8) ::&
      jper2(nsiz2),          &! j values interpolated to new ht grid
      delz(nz1),             &! Difference between tabl. z's
      delang(nza1),          &! Difference between tabl. sec(zen.ang)
      delalb,                &! Difference in tabl. albedo = .5-.2
      pie,                   &! 3.14159...
      dayspy                  ! Number of days per 1 year

   real(r8), parameter :: mwa2so2 = 28.966_r8/64._r8 ! ratio molecular wt of air to so2
   real(r8), parameter :: mwa2h2o2 = 28.966_r8/34._r8 ! ratio molecular wt of air to h2o2


   real(r8), parameter, dimension(nsiz) :: jper =    &! Tabulated j values
(/-11.63010025_r8,-11.53059959_r8,-11.47570038_r8,-11.44130039_r8,-11.40859985_r8,-11.40690041_r8,-11.42399979_r8,-11.46020031_r8,&
  -11.48509979_r8,-11.49100018_r8,-11.46759987_r8,-11.41129971_r8,-11.32050037_r8,-11.20450020_r8,-10.96660042_r8,-10.65079975_r8,&
  -10.31389999_r8,-10.15680027_r8,-11.93789959_r8,-11.81239986_r8,-11.73900032_r8,-11.68900013_r8,-11.62989998_r8,-11.60620022_r8,&
  -11.60490036_r8,-11.61820030_r8,-11.62650013_r8,-11.62040043_r8,-11.58790016_r8,-11.52509975_r8,-11.43029976_r8,-11.31110001_r8,&
  -11.06599998_r8,-10.75179958_r8,-10.39190006_r8,-10.20240021_r8,-12.20409966_r8,-12.06019974_r8,-11.97229958_r8,-11.90939999_r8,&
  -11.82750034_r8,-11.78409958_r8,-11.76560020_r8,-11.75730038_r8,-11.74950027_r8,-11.73139954_r8,-11.68959999_r8,-11.62010002_r8,&
  -11.52070045_r8,-11.39859962_r8,-11.14780045_r8,-10.83170033_r8,-10.46100044_r8,-10.24470043_r8,-12.51459980_r8,-12.35309982_r8,&
  -12.24969959_r8,-12.17240047_r8,-12.06379986_r8,-11.99619961_r8,-11.95629978_r8,-11.92070007_r8,-11.89220047_r8,-11.85890007_r8,&
  -11.80539989_r8,-11.72710037_r8,-11.62139988_r8,-11.49520016_r8,-11.23849964_r8,-10.91720009_r8,-10.54150009_r8,-10.29699993_r8,&
  -13.14570045_r8,-12.95919991_r8,-12.83010006_r8,-12.72649956_r8,-12.56420040_r8,-12.44439983_r8,-12.35589981_r8,-12.25739956_r8,&
  -12.18089962_r8,-12.11240005_r8,-12.03289986_r8,-11.93480015_r8,-11.81389999_r8,-11.67720032_r8,-11.40919971_r8,-11.07330036_r8,&
  -10.69830036_r8,-10.41100025_r8,-14.34010029_r8,-14.15550041_r8,-14.01509953_r8,-13.88949966_r8,-13.65480042_r8,-13.43579960_r8,&
  -13.23670006_r8,-12.97910023_r8,-12.77519989_r8,-12.61489964_r8,-12.47119999_r8,-12.32730007_r8,-12.17119980_r8,-12.00699997_r8,&
  -11.70969963_r8,-11.35060024_r8,-10.96660042_r8,-10.65600014_r8,-15.41930008_r8,-15.26990032_r8,-15.16590023_r8,-15.07689953_r8,&
  -14.89920044_r8,-14.68239975_r8,-14.41380024_r8,-13.96300030_r8,-13.55609989_r8,-13.23820019_r8,-12.98649979_r8,-12.77019978_r8,&
  -12.56400013_r8,-12.36279964_r8,-12.01889992_r8,-11.63490009_r8,-11.22630024_r8,-10.92059994_r8,-17.23080063_r8,-17.09090042_r8,&
  -16.99939919_r8,-16.92959976_r8,-16.82430077_r8,-16.74169922_r8,-16.66230011_r8,-16.49410057_r8,-16.09469986_r8,-15.37609959_r8,&
  -14.65410042_r8,-14.07880020_r8,-13.62609959_r8,-13.26119995_r8,-12.75020027_r8,-12.25749969_r8,-11.79109955_r8,-11.45209980_r8,&
  -11.17109966_r8,-11.14850044_r8,-11.14080048_r8,-11.14089966_r8,-11.15730000_r8,-11.18900013_r8,-11.23019981_r8,-11.29240036_r8,&
  -11.33609962_r8,-11.35649967_r8,-11.34700012_r8,-11.30399990_r8,-11.22620010_r8,-11.12189960_r8,-10.90159988_r8,-10.60239983_r8,&
  -10.27840042_r8,-10.12619972_r8,-11.52420044_r8,-11.47599983_r8,-11.44919968_r8,-11.43290043_r8,-11.42099953_r8,-11.42879963_r8,&
  -11.44970036_r8,-11.48649979_r8,-11.51080036_r8,-11.51659966_r8,-11.49520016_r8,-11.44289970_r8,-11.35809994_r8,-11.24800014_r8,&
  -11.01669979_r8,-10.71510029_r8,-10.36550045_r8,-10.18050003_r8,-11.81890011_r8,-11.75230026_r8,-11.71059990_r8,-11.68099976_r8,&
  -11.64519978_r8,-11.63220024_r8,-11.63479996_r8,-11.64830017_r8,-11.65499973_r8,-11.64729977_r8,-11.61489964_r8,-11.55410004_r8,&
  -11.46290016_r8,-11.34809971_r8,-11.10849953_r8,-10.80249977_r8,-10.44029999_r8,-10.22789955_r8,-12.15250015_r8,-12.06820011_r8,&
  -12.01109982_r8,-11.96700001_r8,-11.90400028_r8,-11.86610031_r8,-11.84650040_r8,-11.83150005_r8,-11.81620026_r8,-11.79199982_r8,&
  -11.74639988_r8,-11.67520046_r8,-11.57610035_r8,-11.45580006_r8,-11.20800018_r8,-10.89459991_r8,-10.52560043_r8,-10.28439999_r8,&
  -12.80440044_r8,-12.69659996_r8,-12.61540031_r8,-12.54609966_r8,-12.43089962_r8,-12.34109974_r8,-12.27270031_r8,-12.19359970_r8,&
  -12.12870026_r8,-12.06789970_r8,-11.99429989_r8,-11.90139961_r8,-11.78509998_r8,-11.65240002_r8,-11.39009953_r8,-11.05939960_r8,&
  -10.68850040_r8,-10.40359974_r8,-13.98359966_r8,-13.87979984_r8,-13.79139996_r8,-13.70489979_r8,-13.52700043_r8,-13.34549999_r8,&
  -13.17129993_r8,-12.93599987_r8,-12.74400043_r8,-12.59039974_r8,-12.45110035_r8,-12.31060028_r8,-12.15719986_r8,-11.99520016_r8,&
  -11.70090008_r8,-11.34430027_r8,-10.96230030_r8,-10.65270042_r8,-15.04539967_r8,-14.97089958_r8,-14.91469955_r8,-14.86209965_r8,&
  -14.74139977_r8,-14.57019997_r8,-14.33740044_r8,-13.92109966_r8,-13.53120041_r8,-13.22140026_r8,-12.97410011_r8,-12.76060009_r8,&
  -12.55630016_r8,-12.35649967_r8,-12.01439953_r8,-11.63179970_r8,-11.22420025_r8,-10.91909981_r8,-16.85479927_r8,-16.78750038_r8,&
  -16.74150085_r8,-16.70509911_r8,-16.64769936_r8,-16.59880066_r8,-16.54520035_r8,-16.40889931_r8,-16.04310036_r8,-15.35249996_r8,&
  -14.64319992_r8,-14.07279968_r8,-13.62240028_r8,-13.25860023_r8,-12.74860001_r8,-12.25650024_r8,-11.79049969_r8,-11.45160007_r8/)

   real(r8), parameter, dimension(nz) :: zz =     &! Tabulated heights
      (/ 0._r8, 1._r8, 2._r8, 3._r8, 5._r8, 7._r8, 9._r8, 12._r8, 15._r8, 18._r8, 21._r8, 24._r8, &
         27._r8, 30._r8, 35._r8, 40._r8, 45._r8, 50._r8 /)

   real(r8), parameter, dimension(nza) :: zasec = &! Tabulated sec(zen.ang)
      (/ 1._r8, 1.3_r8, 1.6_r8, 2._r8, 3._r8, 6._r8, 12._r8, 50._r8 /)

   real(r8), parameter, dimension(nalb) :: alb =  &! Tabulated albedo = 0.2 and 0.5
      (/ 0.2_r8, 0.5_r8 /)

   real(r8), parameter, dimension(121) :: yield = &! yield of DMS+OH rxn going to SO2 
      (/ 0.76905_r8,  0.76953_r8,  0.77002_r8,  0.77050_r8,  0.77098_r8,  0.77146_r8, &
         0.77195_r8,  0.77243_r8,  0.77291_r8,  0.77340_r8,  0.77388_r8,  0.77436_r8, &
         0.77485_r8,  0.77533_r8,  0.77581_r8,  0.77629_r8,  0.77678_r8,  0.77726_r8, &
         0.77774_r8,  0.77822_r8,  0.77870_r8,  0.77918_r8,  0.77967_r8,  0.78015_r8, &
         0.78063_r8,  0.78111_r8,  0.78159_r8,  0.78207_r8,  0.78255_r8,  0.78303_r8, &
         0.78351_r8,  0.78399_r8,  0.78447_r8,  0.78495_r8,  0.78543_r8,  0.78591_r8, &
         0.78640_r8,  0.78688_r8,  0.78737_r8,  0.78786_r8,  0.78835_r8,  0.78884_r8, &
         0.78934_r8,  0.78984_r8,  0.79035_r8,  0.79086_r8,  0.79137_r8,  0.79190_r8, &
         0.79243_r8,  0.79297_r8,  0.79353_r8,  0.79409_r8,  0.79467_r8,  0.79526_r8, &
         0.79587_r8,  0.79650_r8,  0.79715_r8,  0.79783_r8,  0.79853_r8,  0.79927_r8, &
         0.80004_r8,  0.80084_r8,  0.80169_r8,  0.80258_r8,  0.80352_r8,  0.80452_r8, &
         0.80558_r8,  0.80671_r8,  0.80790_r8,  0.80918_r8,  0.81055_r8,  0.81200_r8, &
         0.81356_r8,  0.81522_r8,  0.81700_r8,  0.81890_r8,  0.82093_r8,  0.82309_r8, &
         0.82540_r8,  0.82786_r8,  0.83047_r8,  0.83325_r8,  0.83618_r8,  0.83928_r8, &
         0.84255_r8,  0.84598_r8,  0.84957_r8,  0.85332_r8,  0.85722_r8,  0.86126_r8, &
         0.86543_r8,  0.86973_r8,  0.87412_r8,  0.87861_r8,  0.88316_r8,  0.88777_r8, &
         0.89240_r8,  0.89705_r8,  0.90169_r8,  0.90631_r8,  0.91087_r8,  0.91537_r8, &
         0.91978_r8,  0.92409_r8,  0.92829_r8,  0.93236_r8,  0.93630_r8,  0.94009_r8, &
         0.94373_r8,  0.94721_r8,  0.95054_r8,  0.95370_r8,  0.95670_r8,  0.95954_r8, &
         0.96223_r8,  0.96476_r8,  0.96714_r8,  0.96938_r8,  0.97148_r8,  0.97344_r8, &
         0.97528_r8 /)

   real(r8), parameter, dimension(121) :: dmsrate = &! rate coeff of DMS+OH 
      (/ 0.21261E-10_r8,0.21112E-10_r8,0.20966E-10_r8,0.20823E-10_r8,0.20683E-10_r8, &
         0.20546E-10_r8,0.20411E-10_r8,0.20280E-10_r8,0.20151E-10_r8,0.20024E-10_r8, &
         0.19900E-10_r8,0.19779E-10_r8,0.19660E-10_r8,0.19543E-10_r8,0.19428E-10_r8, &
         0.19316E-10_r8,0.19206E-10_r8,0.19098E-10_r8,0.18991E-10_r8,0.18887E-10_r8, &
         0.18785E-10_r8,0.18684E-10_r8,0.18585E-10_r8,0.18488E-10_r8,0.18393E-10_r8, &
         0.18299E-10_r8,0.18207E-10_r8,0.18116E-10_r8,0.18027E-10_r8,0.17939E-10_r8, &
         0.17852E-10_r8,0.17766E-10_r8,0.17682E-10_r8,0.17598E-10_r8,0.17516E-10_r8, &
         0.17434E-10_r8,0.17354E-10_r8,0.17274E-10_r8,0.17194E-10_r8,0.17115E-10_r8, &
         0.17036E-10_r8,0.16958E-10_r8,0.16880E-10_r8,0.16801E-10_r8,0.16722E-10_r8, &
         0.16643E-10_r8,0.16563E-10_r8,0.16482E-10_r8,0.16401E-10_r8,0.16317E-10_r8, &
         0.16233E-10_r8,0.16146E-10_r8,0.16057E-10_r8,0.15966E-10_r8,0.15871E-10_r8, &
         0.15774E-10_r8,0.15673E-10_r8,0.15567E-10_r8,0.15458E-10_r8,0.15343E-10_r8, &
         0.15223E-10_r8,0.15097E-10_r8,0.14965E-10_r8,0.14826E-10_r8,0.14680E-10_r8, &
         0.14526E-10_r8,0.14365E-10_r8,0.14195E-10_r8,0.14016E-10_r8,0.13828E-10_r8, &
         0.13632E-10_r8,0.13426E-10_r8,0.13211E-10_r8,0.12987E-10_r8,0.12754E-10_r8, &
         0.12514E-10_r8,0.12265E-10_r8,0.12009E-10_r8,0.11747E-10_r8,0.11479E-10_r8, &
         0.11207E-10_r8,0.10932E-10_r8,0.10655E-10_r8,0.10376E-10_r8,0.10098E-10_r8, &
         0.98217E-11_r8,0.95482E-11_r8,0.92787E-11_r8,0.90145E-11_r8,0.87565E-11_r8, &
         0.85057E-11_r8,0.82630E-11_r8,0.80289E-11_r8,0.78042E-11_r8,0.75893E-11_r8, &
         0.73845E-11_r8,0.71899E-11_r8,0.70058E-11_r8,0.68321E-11_r8,0.66687E-11_r8, &
         0.65155E-11_r8,0.63722E-11_r8,0.62384E-11_r8,0.61140E-11_r8,0.59985E-11_r8, &
         0.58915E-11_r8,0.57926E-11_r8,0.57013E-11_r8,0.56173E-11_r8,0.55402E-11_r8, &
         0.54694E-11_r8,0.54047E-11_r8,0.53456E-11_r8,0.52917E-11_r8,0.52426E-11_r8, &
         0.51981E-11_r8,0.51578E-11_r8,0.51214E-11_r8,0.50886E-11_r8,0.50591E-11_r8, &
         0.50327E-11_r8 /)


!##############################################################################
contains
!##############################################################################


  subroutine inisulchem()

    use cam_history,    only: addfld, add_default, phys_decomp

    !----------------------------------------------------------------------- 
    ! Purpose: 
    ! Initialize sulfur cycle chemistry module.
    ! 
    ! Author: B. Eaton
    !----------------------------------------------------------------------- 

    implicit none

    ! Local variables
    integer ::&
         i, l, ik, k

    !-----------------------------------------------------------------------
    !     initialize variables for jh2o2 interpolation
    !     interpolate jval data onto 1 km height grid
    do l = 1, nza*nalb
       do ik = 1, nz2-1
          k = 2
11        continue
          if( zz(k) .gt. ik-1 ) then
             jper2((l-1)*nz2+ik) = jper((l-1)*nz+k-1) &
                  + ((ik-1)-zz(k-1))/(zz(k)-zz(k-1)) &
                  * (jper((l-1)*nz+k)-jper((l-1)*nz+k-1))
          else
             k = k + 1
             go to 11
          endif
       end do
       jper2(l*nz2) = jper(l*nz) 
    end do

    do i = 1,nza1
       delang(i) = 1._r8 / (zasec(i+1) - zasec(i))
    end do
    delalb  = 1._r8 / (alb(2) - alb(1))

    !     Constants from radiation routines needed for dicor()
    dayspy  =  365._r8
    pie     =  4._r8*atan(1._r8)

    call addfld('DMSSNK  ','kg/kg/s' ,pver, 'A', 'DMSsink ',phys_decomp)
    call add_default ('DMSSNK', 1, ' ')
    call addfld('SO2SRCG ','kg/kg/s' ,pver, 'A', 'SO2 Src G   ',phys_decomp)
    call add_default ('SO2SRCG', 1, ' ')
    call addfld('SO2SRCG2','kg/kg/s' ,pver, 'A', 'SO2 Src G2  ',phys_decomp)
    call add_default ('SO2SRCG2', 1, ' ')
    call addfld('SO4SRC  ','kg/kg/s' ,pver, 'A', 'SO4 Src Tot ',phys_decomp)
    call add_default ('SO4SRC', 1, ' ')
    call addfld('SO4SRCG ','kg/kg/s' ,pver, 'A', 'SO4 Src G   ',phys_decomp)
    call add_default ('SO4SRCG', 1, ' ')
    call addfld('H2O2SRC ','kg/kg/s' ,pver, 'A', 'H2O2Src ',phys_decomp)
    call add_default ('H2O2SRC', 1, ' ')
    call addfld('H2O2SNKA','kg/kg/s' ,pver, 'A', 'H2O2SnkA',phys_decomp)
    call add_default ('H2O2SNKA', 1, ' ')
    call addfld('H2O2SNKG','kg/kg/s' ,pver, 'A', 'H2O2SnkG',phys_decomp)
    call add_default ('H2O2SNKG', 1, ' ')

    call addfld ('PH','pH',pver, 'A','Cloud Water pH',phys_decomp)
    call add_default ('PH', 1, ' ')

  end subroutine inisulchem

!##############################################################################
   subroutine sulf_chemdr( t, pmid, cldwat, rain, cldf, cldv, &
                           icwmr1, icwmr2, indcp, ncldypts, hion, &
                           so2, so4, dms,  h2o2, ekso2, ekh2o2, &
                           o3, h2o23d, oh, no3, ho2, q,  &
                           zm, dtime, calday, clat, clon, ncol, lchnk  )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Sulfur chemistry driver separated from chemwdepdr.  The prognostic species 
      ! are DMS, SO2, SO4, and H2O2.  Needed a chem driver without the wet 
      ! deposition for MOZART chemistry
      ! 
      ! Author: F Vitt
      !----------------------------------------------------------------------- 

      use cam_history, only: outfld

      real(r8), intent(in) ::&
         clat(pcols),              &
         clon(pcols),              &
         calday,                   &
         dtime,                    &! time step
         pmid(pcols,pver),         &! pressure at layer midpoints
         zm(pcols,pver),           &! height of layer midpoints
         t(pcols,pver),            &! air temperature (K)
         q(pcols,pver),            &! specific humidity (kg/kg)
         cldwat(pcols,pver),       &! cloud water mixing ratio
         cldf(pcols,pver),         &! total cloud fraction
         cldv(pcols,pver),         &! cloudy volume which is occupied by rain or cloud water 
         rain(pcols,pver),         &! total precip mixing ratio
         icwmr1(pcols,pver),       &! in cloud water mixing ratio for zhang+hack scheme
         icwmr2(pcols,pver)         ! in cloud water mixing ratio for zhang+hack scheme

      real(r8), dimension(pcols,pver), intent(inout) ::&
         so2,         &! so2 array
         so4,         &! so4 array
         dms,         &! dms
         h2o2          ! h2o2 

      real(r8), dimension(pcols,pver), intent(in) ::&
         o3,          &
         no3,         &
         ho2,         &
         oh,          &
         h2o23d

      integer, intent(in) ::    ncol
      integer, intent(in) ::    lchnk

      real(r8), intent(in) ::&
         ekh2o2(pcols,pver)    ! H2O2  Henry's Law coefficient
      real(r8), intent(inout) ::&
         ekso2(pcols,pver)     ! SO2  effective Henry's Law coefficient

      integer, intent(in) ::&
         indcp(pcols,pver),    &! Indices of cloudy grid points
         ncldypts(pver)         ! Number of cloudy grid points

      real(r8), dimension(pcols,pver), intent(inout) ::&
         hion        ! Hydrogen ion concentration in drops

      ! Local variables:

      real(r8), dimension(pcols,pver) ::&
         h2o2new,     &! h2o2
         h2o2tmp,     &! temporary h2o2
         so2new,      &! so2 array
         so2tmp,      &! temporary so2 array
         so4new,      &! so4 array
         so4tmp,      &! temporary so4 array
         dmsnew,      &! dms
         dmssnk,      &! total sink of dms
         so2srcg,     &! total source of so2  
         so2srcg2,    &! dms + oh source of so2  
         so4src,      &! total source of so4  
         so4srcg,     &! gas source of so4  
         so4srca,     &! aq. source of so4  
         h2o2src,     &! total source of h2o2  
         h2o2snkg,    &! gas sink of h2o2  
         h2o2snka      ! aqueous sink of h2o2  

      real(r8), dimension(pcols,pver) ::&
         jh2o2         ! photolysis rate of H2O2

      integer :: i, k

!----------------------------------------------------------------------

     call getj( clat, clon, calday, zm, ncol, jh2o2 )

     call aqchem( so2, so4, h2o2, t, pmid, &
          cldwat, rain, so2new, so4new, h2o2new, &
          o3, cldf, so4srca, hion, h2o2snka, &
          ekso2, ekh2o2, dtime, indcp, ncldypts, &
          cldv, icwmr1, icwmr2, ncol )

     so2tmp(:ncol,:) = so2new(:ncol,:)
     so4tmp(:ncol,:) = so4new(:ncol,:)
     h2o2tmp(:ncol,:) = h2o2new(:ncol,:)

     call gaschem( so2tmp, so4tmp, dms, h2o2tmp, t, &
          pmid, so2new, so4new, dmsnew, h2o2new, &
          h2o23d, oh, no3, ho2, q, &
          jh2o2, so4srcg, dmssnk, so2srcg, so2srcg2, &
          h2o2src, h2o2snkg, dtime, ncol )

#ifdef SCYC_MASSBGT
         call endrun ('sulf_chemdr: gamwet expects scalar for actual argument lat')
         call gamwet( 'dmssnk', lat, pdel, dmssnk )
         call gamwet( 'so2srcg', lat, pdel, so2srcg )
         call gamwet( 'so4srca', lat, pdel, so4srca )
         call gamwet( 'so4srcg', lat, pdel, so4srcg )
#endif

      so2( :ncol,:) = so2new( :ncol,:)
      so4( :ncol,:) = so4new( :ncol,:)
      dms( :ncol,:) = dmsnew( :ncol,:)
      h2o2(:ncol,:) = h2o2new(:ncol,:)

      do k = 1, pver
         do i = 1, ncol
            so4src(i,k) = so4srca(i,k) + so4srcg(i,k)
         end do
      end do

      call outfld( 'DMSSNK  ', dmssnk, pcols, lchnk )
      call outfld( 'SO2SRCG ', so2srcg, pcols, lchnk )
      call outfld( 'SO2SRCG2', so2srcg2, pcols, lchnk )
      call outfld( 'SO4SRC  ', so4src, pcols, lchnk )
      call outfld( 'SO4SRCG ', so4srcg, pcols, lchnk )
      call outfld( 'H2O2SRC ', h2o2src, pcols, lchnk )
      call outfld( 'H2O2SNKA', h2o2snka, pcols, lchnk )
      call outfld( 'H2O2SNKG', h2o2snkg, pcols, lchnk )

   endsubroutine sulf_chemdr


   subroutine aqchem( so2, so4, h2o2, tair, pres, &
                      qc, qr, so2new, so4new, h2o2new, &
                      o3, cldf, so4srca, hion, h2o2snka, &
                      ekso2, ekh2o2, ztodt, ind, ncpts, &
                      cldv, icwmr1, icwmr2, ncol )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate the sources and sinks from the aqueous chemistry:
      !          S(IV) + H2O2 --> SO4 
      !          S(IV) + O3   --> SO4 
      !
      !  where concentrations are concentrations in aqueous phase.
      !
      ! Author: M. Barth
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::&
!         lat(pcols),          &! latitude indices
!         lon(pcols),          &! longititude indices
         ind(pcols,pver),     &! indices for cloudy grid points
         ncpts(pver),         &! number of cloudy grid points
         ncol                  ! number of atmospheric columns used in chunk
      real(r8), intent(in) ::&
         so2(pcols,pver),     &! so2 array
         so4(pcols,pver),     &! so4 from aqueous chemistry 
         h2o2(pcols,pver),    &! h2o2 in kg/kg
         o3(pcols,pver),      &! o3
         tair(pcols,pver),    &! air temperature (K)
         pres(pcols,pver),    &! midpoint pressures
         qc(pcols,pver),      &! cloud water mixing ratio
         qr(pcols,pver),      &! rain mixing ratio
         ztodt,               &! time step
         cldf(pcols,pver)     ! fraction cloud cover

      real(r8), intent(in) ::&
         cldv(pcols,pver),    &! estimate of local volume occupied by clouds
         icwmr1(pcols,pver),  &! in cloud water mixing ration for zhang scheme
         icwmr2(pcols,pver)    ! in cloud water mixing ration for hack  scheme

      real(r8), intent(inout) ::&
         hion(pcols,pver)      ! hydrogen ion concentration (mol/L)

      real(r8), intent(in) ::&
         ekh2o2(pcols,pver)    ! H2O2  Henry's Law coefficient
      real(r8), intent(inout) ::&
         ekso2(pcols,pver)     ! SO2  effective Henry's Law coefficient

      real(r8), intent(out) ::&
         so2new(pcols,pver),  &! so2 array
         so4new(pcols,pver),  &! so4 from aqueous chemistry 
         h2o2new(pcols,pver), &! h2o2 from aqueous chemistry 
         so4srca(pcols,pver), &! so4 production rate from achem
         h2o2snka(pcols,pver)  ! aqueous sink of H2O2

      ! Local variables:

       real(r8)         &
         theff(pcols),  & ! ekso2 at cloudy grid points
         sheff            ! temp array for theff

      integer :: i, k, n
      integer :: i2         ! indice for cloudy grid points
      integer :: ilw        ! number of cloudy grid points
      integer, parameter :: ncyc = 20 ! number of times to loop through aqchem
                                      ! perform aqchem on smaller dt

      real(r8), dimension(pcols,pver) ::&
         cwat,       &! cloud water
         rwat,       &! rain
         twat         ! total liquid water = cwat + rwat
      real(r8) ::&
         h2o2rate,   &! rate of S(IV) + H2O2 reaction 
         o3rate,     &! first order rate coefficient for S(IV) + O3
         molarso4,   &! total so4 in units mol SO4/L H2O
         ahso3        ! (molar)**2 concentration of HSO3-

      real(r8), dimension(pcols) ::&
         tso2,       &! so2 at cloudy grid points
         tso4,       &! so4 at cloudy grid points
         th2o2,      &! h2o2 at cloudy grid points
         to3,        &! o3 at cloudy grid points
         ttwat,      &! twat at cloudy grid points
         temp,       &! tair at cloudy grid points
         tprs,       &! pres at cloudy grid points
         tso2n,      &! so2new at cloudy grid points
         tso4n,      &! so4new at cloudy grid points
         th2o2n,     &! h2o2new at cloudy grid points
         tso4src,    &! so4srca at cloudy grid points
         th2o2snk,   &! h2o2snk at cloudy grid points
         thion,      &! hion at cloudy grid points
         tkh2o2       ! ekh2o2 at cloudy grid points

      real(r8) ::&
         molso2,     &! total so2 in units mol SO2/mol air
         molh2o2,    &! total h2o2 in units mol H2O2/mol air
         molarh2o2,  &! total h2o2 in units mol H2O2/L H2O
         molo3,      &! total o3 in units mol O3/mol air
         eko3,       &! O3   Henry's Law coefficient
         so4src,     &! source of SO4 from SO2 + H2O2
         so2snk,     &! sink of SO2 from SO2 + H2O2
         h2o2snk,    &! sink of H2O2 from SO2 + H2O2
         siv,        &! concentration of S(IV) (mol/L)
         a1, a2,     &! fraction of S(IV) that is HSO3, SO3
         totsb,      &! sum of sulfur at beg of aqchem for one point
         totse,      &! sum of sulfur at end of aqchem for one point
         weight,     &
         tc           ! temperature in deg C

      real(r8) ::&
         dtchem,     &! chemistry timestep
         patm,       &! pressure (atm)
         fa1,        &!  = Khs*K1s/[H+]
         fa2,        &!  = Khs*K1s*K2s/[H+]**2
         khs,        &! Henry's Law constant for SO2
         khsk1s,     &! Khs * K1s
         kh12s,      &! Khs * K1s * K2s
         kh2o2,      &! rate constant for H2O2 + HSO3-
         ko3hs,      &! rate constant for O3 + HSO3-
         ko3so3,     &! rate constant for O3 + SO3=
         qso2,       &! temp array for so2
         qso4,       &! temp array for so4
         qh2o2,      &! temp array for h2o2
         ahion,      &! temp array for ahion
         molrate,    &! rate of rxn in mol/mol-s
         tmprate,    &! temp rate of rxn
         sumrate,    &! h2o2rate + o3rate
         omsm

      ! Function:
!-drb       real(r8) e298, dhr, tt, ek
      real(r8) e298, dhr, tt
      real (r8) ek
      ek(e298, dhr, tt) = e298*exp(dhr*(1._r8/tt - 1._r8/298._r8))
      !----------------------------------------------------------------------

!      call t_startf ('aqchem')

      dtchem = ztodt/ncyc
      omsm = 1._r8 - 1.e-10_r8

      do k=1,pver
         do i=1,ncol

            ! Determine cloud water and rain water mixing ratios
            tc     = tair(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice

            cwat(i,k) = (qc(i,k)/max(cldf(i,k), 0.01_r8) + &
            ! add suspended water from convection and do aqchem in only the liquid phase
                        icwmr1(i,k) + icwmr2(i,k)) * (1._r8-weight)

            if(tair(i,k) .gt. tmelt) then
               weight = 0._r8                    ! fraction of condensate that is ice
            else
               weight = 1._r8
            endif

            rwat(i,k) = qr(i,k)/max(cldv(i,k), 0.01_r8) &
                        *(1._r8-weight) ! aqchem in only the liquid phase

            ! Sum cwat and rwat 
            twat(i,k) = cwat(i,k) + rwat(i,k)

            !        if(twat(i,k) .gt. 0.5e-3) then
            !         tc     = tair(i,k) - tmelt
            !         weight = max(0._r8,min(-tc*0.05,1.0_r8)) 
            !         write(iulog,'(a,3i4,5e12.5)') 'TWAT',i,lat,k, twat(i,k), cwat(i,k),
            !     _    rwat(i,k), icwmr1(i,k)*(1.-weight), icwmr2(i,k)*(1.-weight)
            !         write(iulog,*) ' rwat stuff, qr, cldv, cldf ', qr(i,k), cldv(i,k),
            !     $        cldf(i,k)
            !        endif

         end do
      end do

      ! Initialize arrays
      do k=1,pver
         do i=1,ncol
            so2new(i,k) = so2(i,k)       !keep concentration
            so4new(i,k) = so4(i,k)       
            h2o2new(i,k) = h2o2(i,k)       
            so4srca(i,k) = 0._r8
            h2o2snka(i,k) = 0._r8
         end do
      end do

      ! Initialize cloudy grid point arrays
      !      write(iulog,*) ' Beginning of aqchem'
      do k=1,pver
         ilw = ncpts(k)

         do i2=1,ilw
            tso2(i2)   = so2(ind(i2,k),k)
            tso4(i2)   = so4(ind(i2,k),k)
            th2o2(i2)  = h2o2(ind(i2,k),k)
            to3(i2)    = o3(ind(i2,k),k)
            ttwat(i2)  = twat(ind(i2,k),k)
            temp(i2)   = tair(ind(i2,k),k)
            tprs(i2)   = pres(ind(i2,k),k)
            thion(i2)  = hion(ind(i2,k),k)
            theff(i2)  = ekso2(ind(i2,k),k)
            tkh2o2(i2) = ekh2o2(ind(i2,k),k)
            tso4src(i2) = 0._r8
            th2o2snk(i2) = 0._r8
            !        write(iulog,'(3i5,4e12.5,f8.3,a)') ind(i2,k), lat, k, ttwat(i2), &
            !        thion(i2), tso4(i2), tso2(i2), -log10(thion(i2)), '  begin'
         end do
                    
         ! Aqueous chemistry calculations begin here
         do i2=1,ilw
            patm = tprs(i2)/101325._r8
            khs    = ek(1.23_r8,3120._r8,temp(i2))
            khsk1s = khs * ek(1.3e-2_r8,2015._r8,temp(i2))
            kh12s = khsk1s * ek(6.31e-8_r8,1505._r8,temp(i2)) 
            kh2o2  = 6.2534e14_r8*exp(-4751._r8/temp(i2)) 
            molo3 = to3(i2) 
            eko3  = ek(1.15e-2_r8,2560._r8,temp(i2))
            ko3hs = 4.28e13_r8*exp(-5533._r8/temp(i2)) 
            ko3so3= 7.436e16_r8*exp(-5280._r8/temp(i2)) 

            !        if(lat .eq. 8 .and. k .eq. 14 .and. i2 .gt. 8) then
            !         write(iulog,'(i3,7e11.4)') i2, tso2(i2), th2o2(i2), thion(i2), &
            !          theff(i2), temp(i2), patm, molo3
            !         write(iulog,'(3x,e11.4,f8.1)') ttwat(i2), dtchem
            !        endif
!cdir expand=ncyc
            do n=1,ncyc                     !Do aqchem on a smaller dt   5/5/97

               if(n .eq. 1) then
                  qso2 = tso2(i2)
                  qso4 = tso4(i2)
                  qh2o2 = th2o2(i2)
                  ahion = thion(i2)
                  sheff = theff(i2)
               endif

               fa1 = khsk1s/max(ahion,1.e-30_r8)
               fa2 = kh12s/max((ahion*ahion),1.e-30_r8)

               ! Calculate rate of HSO3 + H2O2 reaction
               molso2  = qso2  * mwa2so2
               molh2o2 = qh2o2 * mwa2h2o2
               ahso3 = fa1 * patm*molso2
               molarh2o2= tkh2o2(i2)*patm*molh2o2

               h2o2rate = kh2o2 * ahion * ahso3 * molarh2o2/(1._r8 + 13._r8*ahion)  

               ! Convert to mol/(mol-s) 
               molrate = h2o2rate*ttwat(i2)*28.97e-3_r8
               ! check for conservation in kg/(kg-s)
               tmprate = mwa2h2o2*min(molrate/mwa2h2o2, qh2o2/dtchem)
               tmprate = mwa2so2*min(tmprate/mwa2so2, qso2/dtchem)

               h2o2rate = tmprate/(max(ttwat(i2), 1.e-30_r8)*28.97e-3_r8)

               ! Calculate rate of o3 oxidation
               siv = sheff * molso2 * patm         ![S(IV)]
               a1 = fa1 / sheff
               a2 = fa2 / sheff

               o3rate  = (ko3hs*a1 + ko3so3*a2) * eko3*molo3 * patm*siv

               ! Convert to mol/(mol-s) 
               molrate = o3rate*ttwat(i2)*28.97e-3_r8
               ! check for conservation in kg/(kg-s)
               tmprate = mwa2so2*min(molrate/mwa2so2, qso2/dtchem)

               o3rate = tmprate/(max(ttwat(i2), 1.e-30_r8)*28.97e-3_r8)

               ! check total rate for conservation
               sumrate = h2o2rate + o3rate
               !         if(lat .eq. 8 .and. k .eq. 14 .and. i2 .gt. 8) then
               !          write(iulog,'(i3,5e11.4)') n, sumrate, h2o2rate, o3rate, &
               !            qso2/dtchem, qh2o2/dtchem
               !         endif
               molrate = sumrate*ttwat(i2)*28.97e-3_r8
               molrate = mwa2so2*min(molrate/mwa2so2, qso2/dtchem)
               tmprate = molrate/(max(ttwat(i2), 1.e-30_r8)*28.97e-3_r8)
               h2o2rate = h2o2rate/max(sumrate, 1.e-30_r8) * tmprate
               o3rate = o3rate/max(sumrate, 1.e-30_r8) * tmprate

               ! Update so2, so4, h2o2
               so2snk = (h2o2rate + o3rate)*ttwat(i2)*1.e-3_r8 * 64._r8
               h2o2snk = h2o2rate*ttwat(i2)*1.e-3_r8 * 34._r8
               !         if(qh2o2 .lt. omsm*h2o2snk*dtchem) then
               !          write(iulog,'(a,2e12.5,3i4)') 'AQCHEM qh2o2 < h2o2snk*dt ', qh2o2, &
               !           h2o2snk*dtchem, lat, ind(i2,k), k
               !         endif
               !         if(qso2 .lt. omsm*so2snk*dtchem) then
               !          write(iulog,'(a,2e12.5,3i4)') 'AQCHEM qso2 < so2snk*dt ', qso2, &
               !           so2snk*dtchem, lat, ind(i2,k), k
               !         endif
               so4src = 1.5_r8*so2snk

               qso2 = max(qso2 - so2snk*dtchem, 0._r8)
               qso4 = qso4 + so4src*dtchem
               qh2o2 = max(qh2o2 - h2o2snk*dtchem, 0._r8)

               tso4src(i2) = tso4src(i2) + so4src
               th2o2snk(i2) = th2o2snk(i2) + h2o2snk

               ! calculate new ahion, sheff
               molso2 = qso2 * mwa2so2
               molarso4 = qso4/(max(ttwat(i2), 1.e-30_r8)*96.e-3_r8) 
               ahso3 = khsk1s * molso2 * patm
!++pjr
! rewrote this expression so we dont take the sqrt(max(a,1.e-60))
!                               but rather max(sqrt(a),1.e-30)
!               ahion = (molarso4 + sqrt(max(molarso4*molarso4 + 4.*ahso3, &
!                       1.e-60_r8)))*0.5
!--pjr
               ahion = (molarso4 + sqrt(max(molarso4*molarso4 + 4._r8*ahso3, &
                       0._r8)))*0.5_r8
               ! Limit to keep pH > 0.0 (ion mixing ratio < 1)
               ahion = min(max(ahion,1.e-30_r8),0.999_r8)

               sheff = khs + khsk1s/ahion + kh12s/max((ahion*ahion),1.e-30_r8)
            end do    !ncyc

            tso2n(i2) = qso2
            tso4n(i2) = qso4
            th2o2n(i2)= qh2o2
            thion(i2) = ahion
            theff(i2) = sheff
            tso4src(i2) = tso4src(i2)/ncyc
            th2o2snk(i2) = th2o2snk(i2)/ncyc

         end do


         do i2=1,ilw
            ! so2new is modified for just the cloudy portion of the grid cell.  Now
            !  adjust it for the entire grid cell.
            so2new(ind(i2,k),k) = cldv(ind(i2,k),k)*tso2n(i2) + &
                                  (1._r8-cldv(ind(i2,k),k))*tso2(i2)
            so4new(ind(i2,k),k) = cldv(ind(i2,k),k)*tso4n(i2) + &
                                  (1._r8-cldv(ind(i2,k),k))*tso4(i2)
            so4srca(ind(i2,k),k) = cldv(ind(i2,k),k)*tso4src(i2) 
            h2o2snka(ind(i2,k),k)= cldv(ind(i2,k),k)*th2o2snk(i2) 
            h2o2new(ind(i2,k),k) = cldv(ind(i2,k),k)*th2o2n(i2) + &
                                   (1._r8-cldv(ind(i2,k),k))*th2o2(i2)
            hion(ind(i2,k),k) = thion(i2)  
            ekso2(ind(i2,k),k) = theff(i2)  
            !         write(iulog,'(3i5,4e12.5,f8.3,a)') ind(i2,k), lat, k, ttwat(i2), &
            !         thion(i2), tso4(i2), tso2(i2), -log10(thion(i2)), '  end'
         end do


         do i=1,ncol
            totsb = 0.5_r8*so2(i,k) + 1._r8/3._r8 * so4(i,k)
            totse = 0.5_r8*so2new(i,k) + 1._r8/3._r8 * so4new(i,k)
            !        if(abs(totsb-totse) .gt. .001*totsb) then
            !         write(iulog,*) 'Total sulfur not conserved in aqueous chemistry'
            !         write(iulog,*) i, k, totsb, totse
            !         stop
            !        endif

         end do

      end do

!      call t_stopf ('aqchem')
      return
   end subroutine aqchem

!##############################################################################

   subroutine cldychmini( tair, pres, qc, qr, cldf, &
                          cldv, icwmr1, icwmr2, so2, so4, &
                          ph, hion, ekso2, ekh2o2, ind, &
                          ncpts, ncol ) 

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate the pH, solubility constants of SO2 and H2O2, and the
      !  location of the cloudy grid points
      ! 
      ! Author: M. Barth
      !----------------------------------------------------------------------- 

      implicit none

      real(r8), intent(in) ::&
         tair(pcols,pver),    &! air temperature (K)
         pres(pcols,pver),    &! midpoint pressures
         qc(pcols,pver),      &! cloud water mixing ratio
         qr(pcols,pver),      &! rain mixing ratio
         cldf(pcols,pver),   &! fraction cloud cover
         cldv(pcols,pver),    &! cloudy volume which is occupied by rain or cloud water
         so2(pcols,pver),     &! so2 array
         so4(pcols,pver),     &! so4 array
         icwmr1(pcols,pver),  &!  in cloud water mixing ration for zhang scheme
         icwmr2(pcols,pver)    !  in cloud water mixing ration for hack  scheme

      real(r8), intent(out) ::&
         ph(pcols,pver),      &! pH of drops
         hion(pcols,pver)      ! hydrogen ion concentration (mol/L)

      real(r8), intent(out) ::         &
         ekh2o2(pcols,pver),  &! H2O2 Henry's Law coefficient
         ekso2(pcols,pver)     ! SO2  effective Henry's Law coefficient

      integer, intent(in) ::  &
         ncol                  ! number of atmospheric columns used in chunk
      integer, intent(out) ::&
         ind(pcols,pver),     &! indices for cloudy grid points
         ncpts(pver)           ! number of cloudy grid points

      ! Local variables:

      real(r8) :: &
         theff(pcols)          ! ekso2 at cloudy grid points

      integer :: &
         i, k,                &
         i2,                  &! index for cloudy grid points
         ilw                   ! number of cloudy grid points

      real(r8), dimension(pcols,pver) ::&
         cwat,       &! cloud water
         rwat,       &! rain
         twat         ! total liquid water = cwat + rwat

      real(r8), dimension(pcols) ::&
         tso2,       &! so2 at cloudy grid points
         tso4,       &! so4 at cloudy grid points
         ttwat,      &! twat at cloudy grid points
         temp,       &! tair at cloudy grid points
         tprs,       &! pres at cloudy grid points
         tph,        &! ph at cloudy grid points
         thion        ! hion at cloudy grid points
      real(r8)  ::   &
         molarso4,   &! total so4 in units mol SO4/L H2O
         ahso3,      &! (molar)**2 concentration of HSO3-
         molso2,     &! so2 in mol/mol units
         weight,     &
         tc           ! temperature in deg C

      ! Function:
!-drb      real(r8) e298, dhr, tt, ek
      real(r8) e298, dhr, tt
      real(r8) ek
      ek(e298, dhr, tt) = e298*exp(dhr*(1._r8/tt - 1._r8/298._r8))
      !----------------------------------------------------------------------

      do k=1,pver
         do i=1,ncol

            ekh2o2(i,k) = ek(7.4e4_r8, 6621._r8, tair(i,k))
            ekso2(i,k)  = ek(1.23_r8,  3120._r8, tair(i,k)) 
            ! Determine cloud water and rain water mixing ratios
            tc     = tair(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice

            cwat(i,k) = (qc(i,k)/max(cldf(i,k), 0.00001_r8) + &
               ! add suspended water from convection and do aqchem in only the liquid phase
                        icwmr1(i,k) + icwmr2(i,k)) * (1._r8-weight)

            if(tair(i,k) .gt. tmelt) then
               weight = 0._r8                    ! fraction of condensate that is ice
            else
               weight = 1._r8
            endif

            rwat(i,k) = qr(i,k)/max(cldv(i,k), 0.00001_r8) &
                        *(1._r8-weight) ! aqchem in only the liquid phase

            ! Sum cwat and rwat 
            twat(i,k) = cwat(i,k) + rwat(i,k)

         end do
      end do

      ! Initialize arrays
      do k=1,pver
         do i=1,ncol
            ph(i,k) = 99._r8
            hion(i,k) = 0._r8
            ind(i,k) = 0
         end do
      end do

      ! Find which grid points have liquid water
      do k=1,pver
         ncpts(k) = 0
         ilw = 0
         do i=1,ncol
            if(cldv(i,k) .ge. 0.02_r8 .and. twat(i,k) .ge. 1.e-12_r8) then
               ilw = ilw + 1
               ind(ilw,k) = i
            endif
         end do

         ncpts(k) = ilw       !assign number of cloudy grid points to array

         do i2=1,ilw
            tso2(i2) = so2(ind(i2,k),k)
            tso4(i2) = so4(ind(i2,k),k)
            ttwat(i2) = twat(ind(i2,k),k)
            temp(i2) = tair(ind(i2,k),k)
            tprs(i2) = pres(ind(i2,k),k)

            ! Set pH as a function of HSO3- and SO4=
            !  Say NH4+ = 1.5 SO4= and NO3- = 0.5 SO4=
            !  Thus, H+ = HSO3- + SO4=  as done in ECHAM
            molso2 = tso2(i2)  * mwa2so2

            molarso4 = tso4(i2)/(max(ttwat(i2), 1.e-30_r8)*96.e-3_r8) 
            ahso3 = ek(1.23_r8,3120._r8,temp(i2)) * &
                    ek(1.3e-2_r8,2015._r8,temp(i2)) * &
                    molso2 * tprs(i2)/101325._r8
!++pjr
! rewrote this expression so we dont take the sqrt(max(a,1.e-50))
!                               but rather max(sqrt(a),1.e-25)
!            thion(i2) = (molarso4 + sqrt(max(molarso4*molarso4 + 4.*ahso3, &
!                         1.e-50_r8)))*0.5
!--pjr
            thion(i2) = (molarso4 + &
                         sqrt(max(molarso4*molarso4 + 4._r8*ahso3,0._r8)))*0.5_r8


!++ag Limit to keep pH > 0.0 (ion mixing ratio < 1)
!           thion(i2) = max(thion(i2),1.e-25_r8)
            thion(i2) = min(max(thion(i2),1.e-25_r8),0.999_r8)
!--ag
            tph(i2) = -log10(thion(i2))

            theff(i2) = ek(1.23_r8, 3120._r8, temp(i2)) * (1._r8 + &
                        ek(1.3e-2_r8, 2015._r8, temp(i2))/thion(i2) * (1._r8 + &
                        ek(6.31e-8_r8, 1505._r8, temp(i2))/thion(i2)))

         end do


         do i2=1,ilw
            hion(ind(i2,k),k) = thion(i2)  
            ph(ind(i2,k),k) = tph(i2)  
            ekso2(ind(i2,k),k) = theff(i2)

         end do
      end do

   end subroutine cldychmini


!##############################################################################

   subroutine coszen( calday, dodiavg, clat, clon, ncol, coszrs )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Compute solar distance factor and cosine solar zenith angle usi
      ! day value where a round day (such as 213.0) refers to 0z at
      ! Greenwich longitude.
      !
      ! Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
      ! Processes in Meterology and Climatology, Elsevier Scientific
      ! Publishing Company, New York  p. 57, p. 62,63.
      ! 
      ! Author: CCM2
      !----------------------------------------------------------------------- 

      implicit none

      real(r8), intent(in) ::&
         calday,              &! Calendar day, including fraction
         clat(pcols),         &! latitudes (radians)
         clon(pcols)

      integer, intent(in) ::&
!         lat(pcols),          &! latitude indices
!         lon(pcols),          &! longitude indices
         ncol                  ! number of atmospheric columns used in chunk

      logical, intent(in) ::&
         dodiavg               ! true => do diurnal averaging

      real(r8), intent(out) ::&
         coszrs(pcols)        ! Cosine solar zenith angle

      ! Local variables
      integer :: i     ! Longitude loop index
      real(r8) ::&
         phi,     &! Greenwich calendar day + local time + long offset
         theta,   &! Earth orbit seasonal angle in radians
         delta,   &! Solar declination angle  in radians
         sinc,    &! Sine   of latitude
         cosc,    &! Cosine of latitude
         sind,    &! Sine   of declination
         cosd      ! Cosine of declination
      real(r8) ::&
         frac,    &! Daylight fraction
         arg,     &!
         tsun,    &! temporary term in diurnal averaging
         coszrsu   ! uniform cosine zenith solar angle 
      !-----------------------------------------------------------------------

      theta = 2._r8*pie*calday/dayspy

      ! Solar declination in radians:
      delta = .006918_r8 - .399912_r8*cos(theta) + .070257_r8*sin(theta) - &
              .006758_r8*cos(2._r8*theta) + .000907_r8*sin(2._r8*theta) - &
              .002697_r8*cos(3._r8*theta) + .001480_r8*sin(3._r8*theta)

      do i=1,ncol
         ! Compute local cosine solar zenith angle,
         sinc = sin(clat(i))
         sind = sin(delta)
         cosc = cos(clat(i))
         cosd = cos(delta)
   
         ! If using diurnal averaging, then compute the average local cosine solar 
         ! zenith angle using formulas from paltridge and platt 1976  p. 57, p. 62,63.
         if (dodiavg) then
            arg = -(sinc/cosc)*(sind/cosd)
            if (arg .lt. -1._r8) then
               frac = 1.0_r8
            else if (arg .gt. 1._r8) then
               frac = 0.0_r8
            else
               frac = (1._r8/pie)*acos(arg)
            endif
            tsun = pie*frac
            if (tsun .gt. 0._r8) then
               coszrsu =  sinc*sind + (cosc*cosd*sin(tsun))/tsun
            else
               coszrsu = 0.0_r8
            endif
            coszrs(i) = coszrsu
         else                       ! No diurnal averaging
   
            ! Calday is the calender day for Greenwich, including fraction
            ! of day; the fraction of the day represents a local time at
            ! Greenwich; to adjust this to produce a true instantaneous time
            ! For other longitudes, we must correct for the local time change:
            ! local time based on the longitude and day of year
            ! then compute the local cosine solar zenith angle
	     coszrs(i) = sinc*sind - cosc*cosd*cos(2._r8*pie*calday + clon(i))
         end if
      end do

   end subroutine coszen

!##############################################################################

   subroutine dicor( calday, clat, ncol, corr )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate a correction factor for the ho2 reaction rate to account for
      ! the diurnal variation of ho2
      !
      ! Author: Mary Barth
      !----------------------------------------------------------------------- 

      implicit none

      real(r8), intent(in) :: calday                ! Calendar day, including fraction
      real(r8), intent(in) :: clat(pcols)           ! latitudes (radians)

      integer, intent(in) ::  ncol                  ! number of atmospheric columns used in chunk

      real(r8), intent(out) :: corr(pcols)          ! correction factor

      !---------------------------Local variables-----------------------------

      integer ic        ! index in chunk
      integer i         ! Longitude loop index
      real(r8) phi      ! Greenwich calendar day + local time + long offset
      real(r8) theta    ! Earth orbit seasonal angle in radians
      real(r8) delta    ! Solar declination angle  in radians
      real(r8) sinc     ! Sine   of latitude
      real(r8) cosc     ! Cosine of latitude
      real(r8) sind     ! Sine   of declination
      real(r8) cosd     ! Cosine of declination
      real(r8) coszrs       ! Cosine solar zenith angle
      real(r8) xxx      ! Work space
      real(r8) yyy      ! Work space
      integer nc        ! counter

      real(r8) ::     dicor_corr          ! h2o reaction rate correction
      integer ntimes    ! number of points to evaluate the diurnal average

      !-----------------------------------------------------------------------
      !

      do ic=1,ncol

         ! Compute solar distance factor and cosine solar zenith angle using
         ! day value where a round day (such as 213.0) refers to 0z at
         ! Greenwich longitude.
         !
         ! Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
         ! Processes in Meterology and Climatology, Elsevier Scientific
         ! Publishing Company, New York  p. 57, p. 62,63.

         theta = 2._r8*pie*calday/dayspy
         !
         ! Solar declination in radians:
         !
         delta = .006918_r8 - .399912_r8*cos(theta) + .070257_r8*sin(theta) - &
              .006758_r8*cos(2._r8*theta) + .000907_r8*sin(2._r8*theta) - &
              .002697_r8*cos(3._r8*theta) + .001480_r8*sin(3._r8*theta)

         ! now calculate the solar zenith angle 
         ! and some useful quantities for the whole day
         !
         ! Compute local cosine solar zenith angle,
         !
         sinc = sin(clat(ic))
         sind = sin(delta)
         cosc = cos(clat(ic))
         cosd = cos(delta)
         !
         ! Calday is the calender day for Greenwich, including fraction
         ! of day; 
         ! since all we care about in this calculation is the diurnal variation
         ! of some quantities we just increment calday over one day.
         !
         xxx = 0._r8
         yyy = 0._r8
         nc  = 0
         ntimes = 48
         do i=1,ntimes
            phi       = calday + (real(i,r8)-1)/real(ntimes,r8)
            coszrs = sinc*sind - cosc*cosd*cos(2._r8*pie*phi)
            if (coszrs.gt.0) then
               nc  = nc + 1
               xxx = xxx + coszrs**2
               yyy = yyy + coszrs
            endif
         end do

         if (yyy.gt.0._r8) then
            dicor_corr = ntimes*xxx/(yyy**2)
         else
            dicor_corr = 1._r8
         endif


         !
         ! Return correction
         !
         corr(ic) = dicor_corr
      end do

    end subroutine dicor

!##############################################################################

   subroutine gaschem( so2, so4, dms, h2o2, tair, &
                       pres, so2new, so4new, dmsnew, h2o2new, &
                       h2o23d, oh, no3, ho2, h2o, &
                       jh2o2, so4srcg, dmssnk, so2srcg, so2srcg2, &
                       h2o2src, h2o2snkg, ztodt, ncol )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate the sources and sinks from the gas chemistry:
      !          DMS + OH     --> a*SO2 + (1-a)*MSA
      !          DMS + NO3    --> SO2 
      !          SO2 + OH + M --> SO4 + M
      !          HO2 + HO2    --> H2O2    
      !          H2O2+ OH     --> HO2 + H2O   
      !          H2O2+ hv     --> 2OH
      !
      !  where a is the yield of SO2 from DMS oxidation, MSA is methane sulfonic
      !   acid which is not carried in this model, and M is an air molecule and
      !   is accounted for in the rate coefficient expression.
      !
      ! Author: M. Barth
      !----------------------------------------------------------------------- 

      implicit none

      real(r8), intent(in) :: &
         so2(pcols,pver),     &! so2 array
         so4(pcols,pver),     &! so4 from gas chemistry 
         dms(pcols,pver),     &! dms
         h2o2(pcols,pver),    &! h2o2
         oh(pcols,pver),      &! oh
         no3(pcols,pver),     &! no3
         ho2(pcols,pver),     &! ho2
         h2o(pcols,pver),     &! h2o
         jh2o2(pcols,pver),   &! photolysis rate of H2O2
         tair(pcols,pver),    &! air temperature (K)
         pres(pcols,pver),    &! midpoint pressures
         ztodt,               &! time step
         h2o23d(pcols,pver)    ! h2o2 -- 3d field from images

      integer, intent(in) ::  &
!         lat(pcols),          &! latitude indices
         ncol                  ! number of atmospheric columns used in chunk

      real(r8), intent(out) ::  &
         so4new(pcols,pver),    &! so4 from gas chemistry 
         dmsnew(pcols,pver),    &! dms
         h2o2new(pcols,pver),   &! h2o2
         so4srcg(pcols,pver),   &
         so2srcg2(pcols,pver)

      real(r8), intent(inout) ::  &
         so2new(pcols,pver),      &! so2 array
         dmssnk(pcols,pver),      &
         so2srcg(pcols,pver),     &
         h2o2src(pcols,pver),     &
         h2o2snkg(pcols,pver)

      ! Local variables:

      integer ::&
         i, k, ik
      logical :: found
      real(r8), parameter ::&         !  from NASA(1992)
         rl300=3.e-31_r8,    &! rate constant at low pressure
         an=3.3_r8,          &! constant:  (300/Tair)**an
         rh300=1.5e-12_r8,   &! rate constant at high pressure
         am=0._r8             ! constant:  (300/Tair)**am

      real(r8) ::&
         air,             &! concentration of air (molec/cm3)
         rlt,             &! used in calculation of SO2 rate coefficient
         rht,             &! used in calculation of SO2 rate coefficient
         ax, bx,          &! used in calculation of SO2 rate coefficient
         psi,             &! used in calculation of SO2 rate coefficient
         so2rate,         &! first order rate coefficient for SO2 + OH
         so2snk,          &
         dmsrate1,        &! first order rate coefficient for DMS + OH
         dmsno3,          &! first order rate coefficient for DMS + NO3
         dmstot            ! first order rate coefficient for both rxns
      real(r8) rk2         ! rate coeff. in HO2 + HO2 rxn
      real(r8) rk3         ! rate coeff. in HO2 + HO2 rxn
      real(r8) fh2o        ! rate coeff. in HO2 + HO2 rxn
      real(r8) ho2rate     ! total rate coeff for HO2 + HO2 rxn
      real(r8) h2o2des     ! rate for H2O2*OH
      real(r8) h2o2phot    ! rate jh2o2*H2O2
      !      real(r8) tauh2o2         !relaxation time for H2O2 generation
      !      real(r8) xh2o2           !zonally avgd H2O2 in kg/kg
      !----------------------------------------------------------------------c
      
      !  Assume a 1.5 day relaxation time for H2O2 generation
      !      tauh2o2 = 129600.    !seconds

!      call t_startf ('gaschem')

      found = .false.
      do k=1,pver
         do i=1,ncol
            ! SO2 + OH rate coefficient:
            if(tair(i,k) .le. 0._r8) then
               found = .true.
               exit
            endif
         end do
      end do
      if (found) then
        do k=1,pver
           do i=1,ncol
              ! SO2 + OH rate coefficient:
              if(tair(i,k) .le. 0._r8) then
                 write(iulog,*) 'TAIR bad ', i, k, tair(i,k)
                 call endrun ('GASCHEM')
              endif
           end do
        end do
      endif
      
      do k=1,pver
         do i=1,ncol

            so2rate = 0._r8

            ! determine concentration of air (molec/cm3) = rho(kg/m3)*Na/MWair
            air = pres(i,k)/(287._r8*tair(i,k)) * 1.e-3_r8/28.966_r8 * 6.022e23_r8
            !         air = 0.8*air                   !N2 concentration
 
            ! SO2 + OH rate coefficient:
            rlt = rl300 * (300._r8/tair(i,k))**an
            rht = rh300 * (300._r8/tair(i,k))**am
            ax = rlt*air/rht 
            psi = 1._r8/(1._r8 + (log10(ax))**2)
            bx = rlt*air/(1._r8 + ax)
            so2rate = bx * 0.6_r8**psi             !rate coef. in cm3/molec-s
            ! rewrite so2rate as first order rate coefficient
            so2rate = so2rate * oh(i,k)

            so2snk = so2(i,k)*(1._r8-exp(-so2rate*ztodt))
            ! Update sulfur species
            so2new(i,k) = so2(i,k)  - so2snk

            so4new(i,k) = so4(i,k) + 1.5_r8*so2snk
            so4srcg(i,k) = 1.5_r8*so2snk/ztodt

            !   Note: 1.5 factor accounts for different molecular weights of so4 and so2
            !    mw(so4)=96;  mw(so2)=64;  mw(so4)/mw(so2)=3/2=1.5
         end do
      end do

      do k=1,pver
         do i=1,ncol

            ! DMS oxidation
            !c DMS + OH rate coefficient and SO2 yield:
            dmstot = 0._r8
            dmsrate1 = 0._r8
            dmsno3 = 0._r8
            ! Get dmsrate from tabulated data
            ik = int(tair(i,k)-195._r8)
            ik = min(max(ik,1), 121)

            if(dms(i,k) .ge. 1.e-30_r8) then
               dmsrate1 = dmsrate(ik)*oh(i,k)

               ! DMS + NO3
               dmsno3 = 1.e-12_r8*exp(-500._r8*(0.0033557_r8 - 1._r8/tair(i,k))) * &
                        no3(i,k)

               ! Total first order rate coeff for DMS
               dmstot = dmsrate1 + dmsno3
            endif

            ! Update sulfur species
            dmssnk(i,k) = dms(i,k)*(1._r8-exp(-dmstot*ztodt)) / ztodt
            dmsnew(i,k) = dms(i,k) - dmssnk(i,k)*ztodt
            so2srcg(i,k) = 1.032258_r8 * (yield(ik) * dms(i,k)*(1._r8-exp(-dmsrate1*ztodt)) + &
                           dms(i,k)*(1._r8-exp(-dmsno3*ztodt)) ) / ztodt
            so2new(i,k) = so2new(i,k) + so2srcg(i,k)*ztodt
            so2srcg2(i,k) = 1.032258_r8 / ztodt * &
                            yield(ik) * dms(i,k)*(1._r8-exp(-dmsrate1*ztodt)) 
            !   Note: 1.032258 = 64/62 accounts for different MW of so2 and dms
         end do
      end do

      do k=1,pver
         do i=1,ncol
            !-----------------------------------------------------------------------
            !c Generate H2O2 via gas phase processes -- parameterizing chemistry
            !        xh2o2 = h2o23d(i,k) * (287.*tair(i,k)/pres(i,k))*34./6.022e20
            !        h2o2new(i,k) = h2o2(i,k) + ztodt/tauh2o2 * (xh2o2 - h2o2(i,k))
            !
            ! Generate H2O2 via the HO2 + HO2 reaction
            !  HO2 concentrations are prescribed
            !  water vapor concentrations taken from q(i,k,1)
            !
            ! determine concentration of air (molec/cm3) = rho(kg/m3)*Na/MWair
            air = pres(i,k)/(287._r8*tair(i,k)) * 1.e-3_r8/28.966_r8 * 6.022e23_r8
            rk2 = 1.7e-12_r8*exp( 600._r8*(1._r8/tair(i,k) - 0.0033557_r8))
            rk3 = 4.9e-32_r8*exp(1000._r8*(1._r8/tair(i,k) - 0.0033557_r8))
            fh2o = 1.4e-21_r8*exp(2200._r8/tair(i,k))
            ho2rate = (rk2 + rk3*air)*(1._r8 + fh2o*h2o(i,k)) * &
                      ho2(i,k)*ho2(i,k)
            !        xh2o2 = h2o23d(i,k) * (287.*tair(i,k)/pres(i,k))*34./6.022e20

            ! Gas-phase destruction of H2O2 occurs via:
            !      H2O2 + hv --> 2OH
            !      H2O2 + OH --> HO2 + H2O
            !  These reactions seem to have the same reaction rate.  Therefore
            !   parameterize this chemistry as doubling the rate of reaction of
            !   H2O2 + OH.
      
            ! h2o2*air*28.966/34. = h2o2 in molec/cm3

            h2o2des = 1.7e-12_r8*exp(-200._r8*(1._r8/tair(i,k) - 0.0033557_r8)) * &
                      oh(i,k) * h2o2(i,k)*air*mwa2h2o2

            h2o2phot = jh2o2(i,k) * h2o2(i,k)*air*mwa2h2o2

            h2o2src(i,k) = ho2rate/(air*mwa2h2o2)
            h2o2snkg(i,k) = (h2o2des + h2o2phot)/(air*mwa2h2o2)
!           h2o2new(i,k) = h2o2(i,k) + (ho2rate - h2o2des - h2o2phot) * &
!                           ztodt/(air*mwa2h2o2)
            h2o2new(i,k) = h2o2(i,k) + (h2o2src(i,k) - h2o2snkg(i,k)) * ztodt

            !c        h2o2new(i,k) = min(h2o2(i,k) + 
            !c     _      ho2rate * ztodt*(287.*tair(i,k)/pres(i,k))*34./6.022e20,
            !c     _      xh2o2)

         end do
      end do

!      call t_stopf ('gaschem')
      
      return
   end subroutine gaschem

!##############################################################################

   subroutine getj( clat, clon, calday, zm, ncol, jout )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Set h2o2 photolysis rates.
      ! 
      ! Author: M. Barth
      !----------------------------------------------------------------------- 

      implicit none

      real(r8), intent(in) ::&
         clat(pcols),           &! latitudes (radians)
         clon(pcols),           &! longitude (radians)
         calday,                &! Current calendar day = julian day + fraction
         zm(pcols,pver)          ! height of midpoint of layer

      integer, intent(in) ::&
!         lat(pcols),            &! latitude indices
!         lon(pcols),            &! longitude indices
         ncol                    ! number of atmospheric columns used in chunk

      real(r8), intent(out) ::&
         jout(pcols,pver)        ! photolysis rate for location

      ! Local variables
      logical dodiavg                 ! when true, diurnally avg zenith angle
      real(r8) cosza(pcols)               ! cosine of zenith angle
      real(r8) dicosza                    ! diurnally avgd cosine of zenith angle
      real(r8) disecza                    ! diurnally avgd cosine of zenith angle

      integer i, is, k
      integer ind(3)
      integer mz
      real(r8) psum
      real(r8) ratza                      !ratio of zenith angles
      real(r8) ratalb                     !ratio of albedos = (0.3-0.2)/(0.5-0.2)
      real(r8) sumrat                     !sum of the ratios
      real(r8) w1                         !1-sumrat
      !-----------------------------------------------------------------------c
      dodiavg = .true.
      call coszen( calday, dodiavg, clat, clon, ncol, cosza )
! TBH:  Loops are not in optimal order for memory access...  
      do i=1,ncol
         dicosza = cosza(i)
   
         if( dicosza .lt. 0.02_r8 ) then
            do k=1,pver
               jout(i,k) = 0._r8
            end do
         else
            disecza = min( 1._r8/max( dicosza, 1.e-7_r8 ), 50._r8 )
      
            ! Interpolate data onto current level and zenith angle for alb=0.3
      
            ! zenith angle:
            do is = 1, 8                            !zasec = za in table
               if( zasec(is) .gt. disecza ) go to 5 !secza = za at gridpt
            end do
            is = is - 1
5           is = max( is-1, 1 )
            ratza = (disecza - zasec(is)) * delang(is)
      
            !     albedo -- set to 0.3 (data is for 0.2 and 0.5)
            ratalb = (0.3_r8 - alb(1)) * delalb
      
            sumrat = ratza + ratalb
            w1 = 1._r8 - sumrat
            
            ! Set indices
            ! note albedo index is either 1 or 2; therefore only one jh2o2 will
            ! have ialb=2
            ind(1) = nz2*(is-1) 
            ind(2) = nz2*is 
            ind(3) = nz2*nza + nz2*(is-1) 
      
            do k = 1, pver
               mz = nint( min( 50._r8, zm(i,k)*0.001_r8 + 1._r8 ) )
               psum = w1*jper2(ind(1)+mz) + ratza*jper2(ind(2)+mz) + &
                      ratalb*jper2(ind(3)+mz)
               jout(i,k) = exp( psum )
            end do
         endif
      end do

   end subroutine getj

!##############################################################################

end module sulchem

