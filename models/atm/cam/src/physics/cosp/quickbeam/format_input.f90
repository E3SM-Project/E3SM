! FORMAT_INPUT: Procedures to prepare data for input to the simulator
! Compiled/Modified:
!   08/28/2006  John Haynes (haynes@atmos.colostate.edu)
!
! irreg_to_grid (subroutine)
! order_data (subroutine)

  module format_input

  contains

! ----------------------------------------------------------------------------
! SUBROUTINE IRREG_TO_GRID
! ----------------------------------------------------------------------------
  subroutine irreg_to_grid(hgt_matrix,t_matrix,p_matrix,rh_matrix, &
    env_hgt_matrix,env_t_matrix,env_p_matrix,env_rh_matrix)
  use array_lib
  implicit none

! Purpose:
!   Linearly interpolate sounding-level data to the hydrometeor-level
!   resolution
!
! Inputs:
!   [hgt_matrix]       hydrometeor-level heights
!   [env_hgt_matrix]   sounding-level heights
!   [env_t_matrix]     sounding-level temperatures
!   [env_p_matrix]     sounding-level pressures
!   [env_rh_matrix]    sounding-level relative humidities
!
! Outputs:
!   [t_matrix]         hydrometeor-level temperatures
!   [p_matrix]         hydrometeor-level pressures
!   [rh_matrix]        hydrometeor-level relative humidities
!
! Created:
!   08/28/2006  John Haynes (haynes@atmos.colostate.edu)

! ----- INPUTS -----
  real*8, dimension(:,:), intent(in) :: &
    hgt_matrix,env_hgt_matrix,env_t_matrix,env_p_matrix,env_rh_matrix

! ----- OUTPUTS -----
  real*8, dimension(:,:), intent(out) :: &
    t_matrix,p_matrix,rh_matrix

! ----- INTERNAL -----
  integer :: nprof, i
  integer,parameter :: KR8 = selected_real_kind(15,300)

  nprof = size(hgt_matrix,1)
  do i=1,nprof
    call lin_interpolate(env_t_matrix(i,:),env_hgt_matrix(i,:), &
      t_matrix(i,:),hgt_matrix(i,:),1000._KR8)
    call lin_interpolate(env_p_matrix(i,:),env_hgt_matrix(i,:), &
      p_matrix(i,:),hgt_matrix(i,:),1000._KR8)
    call lin_interpolate(env_rh_matrix(i,:),env_hgt_matrix(i,:), &
      rh_matrix(i,:),hgt_matrix(i,:),1000._KR8)
  enddo

  end subroutine irreg_to_grid

! ----------------------------------------------------------------------------
! SUBROUTINE ORDER_DATA
! ----------------------------------------------------------------------------
  subroutine order_data(hgt_matrix,hm_matrix,p_matrix,t_matrix, &
    rh_matrix,sfc_radar,hgt_reversed)
  implicit none

! Purpose:
!   Ensure that input data is in top-down order/bottom-up order,
!   for space-based/surface based radars, respectively
!
! Inputs:
!   [hgt_matrix]   heights
!   [hm_matrix]    mixing ratios
!   [t_matrix]     temperatures
!   [p_matrix]     pressures
!   [rh_matrix]    relative humidities
!   [sfc_radar]    1=surface radar, 0=spaceborne
!
! Outputs:
!   [hgt_matrix],[hm_matrix],[p_matrix,[t_matrix],[rh_matrix] in proper
!   order for input to the radar simulator routine
!   [hgt_reversed]   T=heights were reordered,F=heights were not reordered
!
! Note:
!   The order for all profiles is assumed to the same as the first profile.
!
! Created:
!   08/28/2006  John Haynes (haynes@atmos.colostate.edu)

! ----- INPUTS -----
  integer, intent(in) :: sfc_radar

! ----- OUTPUTS -----
  real*8, dimension(:,:), intent(inout) :: &
    hgt_matrix,p_matrix,t_matrix,rh_matrix
  real*8, dimension(:,:,:), intent(inout) :: &
    hm_matrix
  logical, intent(out) :: hgt_reversed

! ----- INTERNAL -----
  integer :: ngate
  logical :: hgt_descending
  

  ngate = size(hgt_matrix,2)
  hgt_descending = hgt_matrix(1,1) > hgt_matrix(1,ngate)
      
! :: surface: heights must be ascending
! :: space-based: heights must be descending
  if ( &
     (sfc_radar == 1 .and. hgt_descending) .or.  &
     (sfc_radar == 0 .and. (.not. hgt_descending)) &
     ) &
  then

    hgt_matrix(:,:) = hgt_matrix(:,ngate:1:-1)
    hm_matrix(:,:,:) = hm_matrix(:,:,ngate:1:-1)
    p_matrix(:,:) = p_matrix(:,ngate:1:-1)
    t_matrix(:,:) = t_matrix(:,ngate:1:-1)
    rh_matrix(:,:) = rh_matrix(:,ngate:1:-1) 

    hgt_reversed = .true.
  else
    hgt_reversed = .false.
  endif

  end subroutine order_data

  end module format_input
