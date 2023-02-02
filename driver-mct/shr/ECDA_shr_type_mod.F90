!BOP
!   !MODULE: ecda_shr_type_mod
!   !INTERFACE:
module ECDA_shr_type_mod
!   !DESCRIPTION:
!   ! Define the state/obs type
!   !
!   !REVISION HISTORY:
!    Jul  2019 Y. Liu <liu6@tamu.edu> initial release
!   !USES:

  use shr_kind_mod, only : r8=> shr_kind_r8

  integer, public, parameter :: filename_len      = 256
  integer, public, parameter :: type_name_len     = 72
  integer, public, parameter :: var_name_len      = 72

  type, public ::  ECDA_state_type
       character(len=type_name_len)                ::   name              ! state variable  name,ex : ocn_temp,ocn_salt...
       integer                                     ::   state_type = 1    ! state type for mpas grid: 1 cells (default); 2,edgers; 3 Vertex 
       integer                                     ::   state_dim         ! dimension of state variables (1, 2,3 for now) 
       logical                                     ::   required          ! if required for DA by obsope or update 
       logical                                     ::   update            ! if  final analysis results to updated in model, IAU 
       logical                                     ::   averaged          ! if the observations averaged
       real(r8)                                    ::   r2p_fac           ! the relaxiation to prior factor, current use Zhang F. scheme
       real(r8), allocatable,dimension(:)          ::   KKlev             ! the vertical coordinator for model state
       real(r8)                                    ::   mx_value          ! Max values of the state 
       real(r8)                                    ::   mn_value          ! Min values of the state 
       
       != model original grid
       real(r8), allocatable,dimension(:,:)        ::   h_da              ! the water depth for the grid
       real(r8), allocatable,dimension(:,:)        ::   lon_da            ! longitude 
       real(r8), allocatable,dimension(:,:)        ::   lat_da            ! latitude  

       real(r8), allocatable,dimension(:)        ::   h_da_1d              ! the water depth for the grid
       real(r8), allocatable,dimension(:)        ::   lon_da_1d            ! longitude 
       real(r8), allocatable,dimension(:)        ::   lat_da_1d            ! latitude  

       real(r8), allocatable,dimension(:,:,:,:)    ::   prior             ! the original prior value, ensemble
       real(r8), allocatable,dimension(:,:,:)      ::   prior2d           ! the original prior value for 2D, ensemble
       real(r8), allocatable,dimension(:,:)        ::   prior1d           ! the original prior value for 2D, ensemble
       real(r8), allocatable,dimension(:,:,:,:)    ::   posterior         ! the original posterior value before final inflation,ensemble
       real(r8), allocatable,dimension(:,:,:)      ::   posterior2d       ! the original posterior value for 2D, ensemble
       real(r8), allocatable,dimension(:,:)        ::   posterior1d       ! the original posterior value for 1D, ensemble

       real(r8), allocatable,dimension(:,:,:,:)    ::   prior_ave         ! the original prior ave value ,ensemble
       real(r8), allocatable,dimension(:,:,:)      ::   prior2d_ave       ! the original prior ave value for 2D, ensemble
       real(r8), allocatable,dimension(:,:)        ::   prior1d_ave       ! the original prior ave value for 1D, ensemble

       != model original grid
       integer                                     ::   ave_num           ! the total average number
       real(r8), allocatable,dimension(:,:,:)      ::   ave_value         !  average values for average observations
       real(r8), allocatable,dimension(:,:)        ::   ave_value2d       !  average values for average observations (2D) 
       real(r8), allocatable,dimension(:)          ::   ave_value1d       !  average values for average observations (21D) 
       real(r8), allocatable,dimension(:,:,:)      ::   IAU_increm        ! IAU increment at current PE
       real(r8), allocatable,dimension(:,:)        ::   IAU_increm2d      ! IAU increment at current PE
       real(r8), allocatable,dimension(:)          ::   IAU_increm1d      ! IAU increment at current PE
  end type ECDA_state_type


  type, public ::  ECDA_obs_type
       character(len=type_name_len)                ::   obs_name          ! observation name, use for identify observations, it is the observation namelist name
       integer                                     ::   obs_type          ! observation type :  ex:TEMP_ID,SALT_ID, use for obs_ope
       character(len=filename_len)                 ::   filename          ! the obs file name or part
       character(len=var_name_len)                 ::   varname           ! the obs variable name in ncfile 
       character(len=var_name_len)                 ::   uncert_name           ! the obs variable name in ncfile 
       character(len=var_name_len)                 ::   lonname           ! the lon coordinator name in ncfile 
       character(len=var_name_len)                 ::   latname           ! the lat coordinator name in ncfile 
       character(len=var_name_len)                 ::   vertical_name     ! the vertical coordinator name in ncfile 
       character(len=var_name_len)                 ::   dm_name           ! use for profile observation 
       character(len=10)                           ::   init_obs_YMDH     ! the inital time to observations, use to calculate record
       integer                                     ::   uncert_type       ! variance 0, std 1
       integer                                     ::   init_obs_YR       ! ex:2016
       integer                                     ::   init_obs_MO       ! ex: 1~12
       integer                                     ::   init_obs_DY       ! ex:1~31
       integer                                     ::   init_obs_HR       ! ex:0~23
       integer                                     ::   grid_type         ! netcdf file, 1,2 are grid file, 1: model grid, 2 other resolution, need read coordinator information, 3: station data, (lon,lat,lev, value,error) 
       integer                                     ::   obs_dm            ! data dimension in file, only used for grid file
       integer                                     ::   start_layer       ! data start from (vertical layer) for 3D,only used for grid file
       integer                                     ::   nlayer            ! total vertical layer for 3D,only used for grid file
       integer                                     ::   ilayer            ! pick vertical layer ,only used for grid file
       logical                                     ::   assimilate        ! if this observation avaliable 
       integer                                     ::   inst_type         ! instrument types are defined by platform class (e.g. MOORING, DROP) and instrument type (XBT, CTD, ...)
       integer                                     ::   nobs              ! the number of observation record at current PE
       logical, allocatable,dimension(:)           ::   required_states   ! the required state variables for obsope 
       logical, allocatable,dimension(:)           ::   update_states     ! the updated  state variables for current variable 
       real(r8),allocatable,dimension(:)           ::   skip_sfc          ! skip surface for ssh lik observation 
       real(r8), allocatable,dimension(:)          ::   rx, ry,rz         ! obs location in grids in real 
       real(r8), allocatable,dimension(:)          ::   lon,lat,depth     ! obs location in longitude and latitude and depth (ocn:negtive, atm:positive)
       real(r8), allocatable,dimension(:)          ::   values            ! the observation values
       real(r8), allocatable,dimension(:)          ::   err_var           ! the observation error variance. 
       real(r8)                                    ::   err_var0          ! specify constant observation error variance (>0). When <=0, input from observation file.
       real(r8)                                    ::   adj_scale = 1.0   ! some variable need adjust scale
       logical,  allocatable,dimension(:)          ::   valid             ! partial DA switch, to be developed
       integer, allocatable,dimension(:)           ::   qc_flag           ! observation quality flag
       real(r8)                                    ::   obs_range_min     ! observation need larger than the obs_range_min 
       real(r8)                                    ::   obs_range_max     ! observation need smaller than obs_range_max
       real(r8)                                    ::   RR_xy,RR_z        ! half of localization cutoff for spatial and vertical
       real(r8)                                    ::   Mx_PO_dist        ! maxium distance allowed between observation and forecast
       real(r8)                                    ::   Mx_misfit         ! maxium rate allowed for inflate  observation uncertainty

       integer                                     ::   gl_nobs                  ! global observation record for letkf
       real(r8), allocatable,dimension(:)          ::   gl_lon,gl_lat,gl_depth   ! global obs location in longitude and latitude and depth (ocn:-, atm:+)
       real(r8), allocatable,dimension(:,:)        ::   fcst_obs                 ! local forecast observations


  end type ECDA_obs_type

end module ECDA_shr_type_mod

