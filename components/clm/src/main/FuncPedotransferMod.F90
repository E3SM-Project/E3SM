module FuncPedotransferMod
!
!DESCRIPTIONS:
!module contains different pedotransfer functions to
!compute the mineral soil hydraulic properties.
!currenty, only the Clapp-Hornberg formulation is used.
!HISTORY:
!created by Jinyun Tang, Mar.1st, 2014
implicit none
  private
  public :: pedotransf
  public :: get_ipedof
  public :: init_pedof

  integer, parameter :: cosby_1984_table5 = 0    !by default uses this form
  integer, parameter :: cosby_1984_table4 = 1
  integer, parameter :: noilhan_lacarrere_1995 = 2
  integer :: ipedof0
contains

   subroutine init_pedof()
   !
   !DESCRIPTIONS
   !initialize the default pedotransfer function
   implicit none
   
   
   ipedof0 = cosby_1984_table5         !the default pedotransfer function
   end subroutine init_pedof
   
   subroutine pedotransf(ipedof, sand, clay, watsat, bsw, sucsat, xksat)
   !pedotransfer function to compute hydraulic properties of mineral soil
   !based on input soil texture
   
   use shr_kind_mod         , only : r8 => shr_kind_r8
   use abortutils    , only : endrun         
   implicit none
   integer,  intent(in) :: ipedof !type of pedotransfer function, use the default pedotransfer function  
   real(r8), intent(in) :: sand   !% sand
   real(r8), intent(in) :: clay   !% clay
   real(r8), intent(out):: watsat !v/v saturate moisture
   real(r8), intent(out):: bsw    !b shape parameter
   real(r8), intent(out):: sucsat !mm, soil matric potential
   real(r8), intent(out):: xksat  !mm/s, saturated hydraulic conductivity
   
   character(len=32) :: subname = 'pedotransf'  ! subroutine name 
   select case (ipedof)
   case (cosby_1984_table4)
      call pedotransf_cosby1984_table4(sand, clay, watsat, bsw, sucsat, xksat)
   case (noilhan_lacarrere_1995)
      call pedotransf_noilhan_lacarrere1995(sand, clay, watsat, bsw, sucsat, xksat) 
   case (cosby_1984_table5)  
      call pedotransf_cosby1984_table5(sand, clay, watsat, bsw, sucsat, xksat)
   case default
      call endrun(subname // ':: a pedotransfer function must be specified!')  
   end select

   end subroutine pedotransf
   
!------------------------------------------------------------------------------------------
   subroutine pedotransf_cosby1984_table4(sand, clay, watsat, bsw, sucsat, xksat)
   !
   !DESCRIPTIONS
   !compute hydraulic properties based on functions derived from Table 4 in cosby et al, 1984
   use shr_kind_mod         , only : r8 => shr_kind_r8   
   implicit none
   real(r8), intent(in) :: sand   !% sand
   real(r8), intent(in) :: clay   !% clay
   real(r8), intent(out):: watsat !v/v saturate moisture
   real(r8), intent(out):: bsw    !b shape parameter
   real(r8), intent(out):: sucsat !mm, soil matric potential
   real(r8), intent(out):: xksat  !mm/s, saturated hydraulic conductivity
   
   !Cosby et al. Table 4
   watsat = 0.505_r8-0.00142_r8*sand-0.00037*clay
   bsw = 3.10+0.157*clay-0.003*sand            
   sucsat  = 10._r8 * ( 10._r8**(1.54_r8-0.0095_r8*sand+0.0063*(100._r8-sand-clay)))            
   xksat         = 0.0070556 *(10.**(-0.60+0.0126*sand-0.0064*clay) )     !mm/s now use table 4.
      
   end subroutine pedotransf_cosby1984_table4
   
!------------------------------------------------------------------------------------------
   subroutine pedotransf_cosby1984_table5(sand, clay, watsat, bsw, sucsat, xksat)
   !
   !DESCRIPTIONS
   !compute hydraulic properties based on functions derived from Table 5 in cosby et al, 1984
   
   use shr_kind_mod         , only : r8 => shr_kind_r8   
   implicit none
   real(r8), intent(in) :: sand   !% sand
   real(r8), intent(in) :: clay   !% clay
   real(r8), intent(out):: watsat !v/v saturate moisture
   real(r8), intent(out):: bsw    !b shape parameter
   real(r8), intent(out):: sucsat !mm, soil matric potential
   real(r8), intent(out):: xksat  !mm/s, saturated hydraulic conductivity
   
   !Cosby et al. Table 5     
   watsat = 0.489_r8 - 0.00126_r8*sand
   bsw    = 2.91 + 0.159*clay
   sucsat = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )            
   xksat         = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s, from table 5 
      
   end subroutine pedotransf_cosby1984_table5
   
!------------------------------------------------------------------------------------------
   subroutine pedotransf_noilhan_lacarrere1995(sand, clay, watsat, bsw, sucsat, xksat)
   !
   !DESCRIPTIONS
   !compute hydraulic properties based on functions derived from Noilhan and Lacarrere, 1995
   
   use shr_kind_mod         , only : r8 => shr_kind_r8   
   implicit none
   real(r8), intent(in) :: sand   !% sand
   real(r8), intent(in) :: clay   !% clay
   real(r8), intent(out):: watsat !v/v saturate moisture
   real(r8), intent(out):: bsw    !b shape parameter
   real(r8), intent(out):: sucsat !mm, soil matric potential
   real(r8), intent(out):: xksat  !mm/s, saturated hydraulic conductivity
   
   !Noilhan and Lacarrere, 1995
   watsat = -0.00108*sand+0.494305
   bsw = 0.137*clay + 3.501
   sucsat = 10._r8**(-0.0088*sand+2.85)
   xksat = 10._r8**(-0.0582*clay-0.00091*sand+0.000529*clay**2._r8-0.0001203*sand**2._r8-1.38)
   end subroutine pedotransf_noilhan_lacarrere1995
!------------------------------------------------------------------------------------------
   function get_ipedof(soil_order)result(ipedof)
   !
   ! DESCRIPTION
   ! select the pedotransfer function to be used
   implicit none
   integer, intent(in) :: soil_order
   
   integer :: ipedof
   
   if(soil_order==0)then
      ipedof=ipedof0
   endif
   
   end function get_ipedof   
end module FuncpedotransferMod
