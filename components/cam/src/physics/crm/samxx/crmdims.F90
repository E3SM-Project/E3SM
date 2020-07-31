module crmdims
    use params, only: crm_rknd, crm_iknd
    implicit none

    integer(crm_iknd), parameter :: crm_nx=CRM_NX
    integer(crm_iknd), parameter :: crm_ny=CRM_NY
    integer(crm_iknd), parameter :: crm_nz=CRM_NZ

    integer(crm_iknd), parameter :: crm_nx_rad=CRM_NX_RAD
    integer(crm_iknd), parameter :: crm_ny_rad=CRM_NY_RAD

    real(crm_rknd), parameter :: crm_dx=CRM_DX
    real(crm_rknd), parameter :: crm_dy=crm_dx
    real(crm_rknd), parameter :: crm_dt=CRM_DT

end module crmdims
