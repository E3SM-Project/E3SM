module CLMVICMapMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs  the mapping from CLM layers to VIC layers
  ! Specifically, 10 (or 23 when more_vertlayers == .true.) 
  ! CLM hydrologically active soil layers are mapped to three VIC layers
  ! by assigning the first nlvic(1) layers to VIC layer 1
  !              the next nlvic(2) layers  to VIC alyer 2
  !              and the remaining to VIC layer 3
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: initCLMVICMap   ! map layer/node fractions
  public  :: CLMVICMap          ! map from VIC to CLM layers
  private :: linear_interp      ! function for linear interperation 
  !
  ! !REVISION HISTORY:
  ! Created by Aihui Wang, 2008
  ! Revised by Maoyi Huang, 02/12/2010
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initCLMVICMap(c)

    ! !DESCRIPTION:
    ! This subroutine calculates mapping between CLM and VIC layers
    ! added by AWang
    ! modified by M.Huang for CLM4 
    !
    ! !USES:
    use clmtype
    use clm_varcon  , only : denh2o, denice, pondmx
    use clm_varpar  , only : nlevsoi, nlayer, nlayert, nlevgrnd 
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: c
    !
    ! !REVISION HISTORY:
    ! Created by Maoyi Huang
    ! 11/13/2012, Maoyi Huang: rewrite the mapping modules in CLM4VIC 
    !
    real(r8) :: sum_frac(1:nlayer)                  ! sum of fraction for each layer
    real(r8) :: deltal(1:nlayer+1)                  ! temporary
    real(r8) :: zsum                                ! temporary
    real(r8) :: lsum                                ! temporary
    real(r8) :: temp                                ! temporary

    ! other local variables
   
    integer :: i, j, fc
    ! note: in CLM h2osoil_liq unit is kg/m2, in VIC moist is mm
    ! h2osoi_ice is actually water equavlent ice content.
    !-----------------------------------------------------------------------
    
   associate(& 
   dz            =>    cps%dz             , & ! Input:  [real(r8) (:,:)]  layer depth (m)                       
   zi            =>    cps%zi             , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m) 
   z             =>    cps%z              , & ! Input:  [real(r8) (:,:)]  layer thickness (m)                   
   depth         =>    cps%depth          , & ! Input:  [real(r8) (:,:)]  layer depth of VIC (m)                
   vic_clm_fract =>    cps%vic_clm_fract    & ! Input:  [real(r8) (:,:,:)]  fraction of VIC layers in clm layers
   )
   !************************************************************************  

   !  set fraction of VIC layer in each CLM layer
 
   lsum = 0._r8
   do i = 1, nlayer
      deltal(i) = depth(c,i)
   end do
   do i = 1, nlayer
      zsum = 0._r8
      sum_frac(i) = 0._r8
      do j = 1, nlevsoi
         if( (zsum < lsum) .and. (zsum + dz(c,j) >= lsum ))  then
             call linear_interp(lsum, temp, zsum, zsum + dz(c,j), 0._r8, 1._r8)
             vic_clm_fract(c,i,j) = 1._r8 - temp
             if(lsum + deltal(i) < zsum + dz(c,j)) then
                   call linear_interp(lsum + deltal(i), temp, zsum, zsum + dz(c,j), 1._r8, 0._r8)
                   vic_clm_fract(c,i,j) = vic_clm_fract(c,i,j) - temp
             end if
         else if( (zsum < lsum + deltal(i)) .and. (zsum + dz(c,j) >= lsum + deltal(i)) ) then
             call linear_interp(lsum + deltal(i), temp, zsum, zsum + dz(c,j), 0._r8, 1._r8)
             vic_clm_fract(c,i,j) = temp
             if(zsum<=lsum) then
                 call linear_interp(lsum, temp, zsum, zsum + dz(c,j), 0._r8, 1._r8)
                 vic_clm_fract(c,i,j) = vic_clm_fract(c,i,j) - temp
             end if
         else if( (zsum >= lsum .and. zsum + dz(c,j) <= lsum + deltal(i)) )  then
              vic_clm_fract(c,i,j) = 1._r8
         else
              vic_clm_fract(c,i,j) = 0._r8
         end if
         zsum = zsum + dz(c,j)
         sum_frac(i) = sum_frac(i) + vic_clm_fract(c,i,j)
     end do                           ! end CLM layer calculation
     lsum = lsum + deltal(i)
   end do                             ! end VIC layer calcultion 
   
    end associate 
 end subroutine initCLMVICMap

 !-----------------------------------------------------------------------
 subroutine CLMVICMap(bounds, numf, filter)
   !
   ! !DESCRIPTION:
   ! mapping from VIC to CLM layers, M.Huang
   !
   ! !USES:
   use clmtype
   use clm_varcon  , only : denh2o, denice, pondmx, watmin
   use clm_varpar  , only : nlevsoi, nlayer, nlayert, nlevgrnd 
   use decompMod   , only : bounds_type
   !
   ! !REVISION HISTORY:
   ! Created by Maoyi Huang
   ! 11/13/2012, Maoyi Huang: rewrite the mapping modules in CLM4VIC 
   !
   ! !ARGUMENTS:
   implicit none
   type(bounds_type), intent(in) :: bounds  ! bounds
   integer , intent(in)  :: numf                   ! number of column soil points in column filter
   integer , intent(in)  :: filter(:)      ! column filter for soil points
   !
   ! !LOCAL VARIABLES
   real(r8) :: ice0(1:nlayer)            ! last step ice lens (mm)  (new)
   real(r8) :: moist0(1:nlayer)          ! last step soil water (mm)  (new)
   integer  :: i, j, c, fc
   ! note: in CLM3 h2osoil_liq unit is kg/m2, in VIC moist is mm
   ! h2osoi_ice is actually water equavlent ice content.
   !-----------------------------------------------------------------------

   associate(& 
   dz             => cps%dz           , & ! Input:  [real(r8) (:,:)] layer depth (m)                        
   zi             => cps%zi           , & ! Input:  [real(r8) (:,:)] interface level below a "z" level (m)  
   z              => cps%z            , & ! Input:  [real(r8) (:,:)] layer thickness (m)                    
   h2osoi_liq     => cws%h2osoi_liq   , & ! Input:  [real(r8) (:,:)] liquid water (kg/m2)                   
   h2osoi_ice     => cws%h2osoi_ice   , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2)                       
   moist          => cws%moist        , & ! Input:  [real(r8) (:,:)] liquid water (mm)                      
   ice            => cws%ice          , & ! Input:  [real(r8) (:,:)] ice lens (mm)                          
   h2osoi_vol     => cws%h2osoi_vol   , & ! Input:  [real(r8) (:,:)] volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
   moist_vol      => cws%moist_vol    , & ! Input:  [real(r8) (:,:)] volumetric soil moisture for VIC soil layers
   porosity       => cps%porosity     , & ! Input:  [real(r8) (:,:)] soil porisity (1-bulk_density/soil_density)
   depth          => cps%depth        , & ! Input:  [real(r8) (:,:)] layer depth of upper layer (m)         
   max_moist      => cps%max_moist    , & ! Input:  [real(r8) (:,:)] max layer moist + ice (mm)             
   vic_clm_fract  => cps%vic_clm_fract  & ! Input:  [real(r8) (:,:,:)] fraction of VIC layers in each CLM layer
   )

   ! map CLM to VIC
   do fc = 1, numf
      c = filter(fc)
      do i = 1, nlayer
         ice0(i) = ice(c,i)
         moist0(i) = moist(c,i)
         ice(c,i) = 0._r8
         moist(c,i) = 0._r8
         do j = 1, nlevsoi
            ice(c,i) = ice(c,i) + h2osoi_ice(c,j) * vic_clm_fract(c,i,j)
            moist(c,i) = moist(c,i) + h2osoi_liq(c,j) * vic_clm_fract(c,i,j)
         end do
         ice(c,i) = min((moist0(i) + ice0(i)), ice(c,i))
         ice(c,i) = max(0._r8, ice(c,i))
         moist(c,i) =max(watmin, moist(c,i))
         moist(c,i) =min(max_moist(c,i)-ice(c,i), moist(c,i))
         moist_vol(c,i) = moist(c,i)/(depth(c,i)*denice) &
                      + ice(c,i)/(depth(c,i)*denh2o)
         moist_vol(c,i) = min(porosity(c,i), moist_vol(c,i))
         moist_vol(c,i) = max(0.01_r8, moist_vol(c,i))
      end do

     ! hydrologic inactive layers
     ice(c, nlayer+1:nlayert) = h2osoi_ice(c, nlevsoi+1:nlevgrnd)
     moist(c, nlayer+1:nlayert) = h2osoi_liq(c, nlevsoi+1:nlevgrnd)
     moist_vol(c, nlayer+1:nlayert) = h2osoi_vol(c, nlevsoi+1:nlevgrnd)
  end do

    end associate 
 end subroutine CLMVICMap

 !-------------------------------------------------------------------
 subroutine linear_interp(x,y, x0, x1, y0, y1)
   !
   ! !DESCRIPTION:
   ! This subroutine provides linear interpolation
   !
   ! !USES:  
   !
   ! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: x, x0, y0, x1, y1
   real(r8), intent(out) :: y
   !-------------------------------------------------------------------


   y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)

 end subroutine linear_interp

end module CLMVICMapMod
