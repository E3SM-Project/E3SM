module test_mkpftMod
! Module for testing mkpftMod

  use mkpftMod
  use test_mod
  use shr_kind_mod, only : r8 => shr_kind_r8
  
  implicit none
  private
  
  public :: test_mkpft_normalize
  
  character(len=*), parameter :: modname = 'test_mkpftMod'

contains
   
!------------------------------------------------------------------------------
  subroutine test_mkpft_normalize

    use mkvarctl, only : numpft

    implicit none
    
    character(len=128) :: testname

    real(r8), allocatable :: pctpft_full(:)
    real(r8), allocatable :: pct_pft(:)
    real(r8), allocatable :: pct_cft(:)
    real(r8), allocatable :: pct_pft_t(:)
    real(r8), allocatable :: pct_cft_t(:)
    real(r8) :: pct_special, pct_natveg, pct_crop, pct_natveg_t, pct_crop_t
    
    real(r8), parameter :: eps = 1.e-12_r8

    character(len=*), parameter :: subname = 'test_mkpft_normalize'


    ! TESTS WITH NO SPECIFIC CROPS

    ! Set some module-level variables in mkpftMod
    numpft = 16
    num_natpft = numpft
    num_cft = 0
    natpft_lb = 0
    natpft_ub = num_natpft
    cft_lb = num_natpft+1
    cft_ub = cft_lb + num_cft - 1
    allocate(pctpft_full(0:numpft), pct_pft(natpft_lb:natpft_ub), pct_cft(cft_lb:cft_ub))
    allocate(pct_pft_t(natpft_lb:natpft_ub), pct_cft_t(cft_lb:cft_ub))

    testname='no specific crop, sum(pctpft_full) = 0, pct_special = 100'
    pctpft_full = (/0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,&
                    0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8/)
    pct_special = 0._r8
    pct_natveg_t = 0._r8
    pct_crop_t = 0._r8
    pct_pft_t  = (/100._r8, 0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,&
                   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8/)
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.false.)
    
    ! same, but pct_special = 25
    testname='no specific crop, sum(pctpft_full) = 0, pct_special = 25'
    pct_special = 25._r8
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.false.)
    
    ! same, but pct_special = 0
    testname='no specific crop, sum(pctpft_full) = 0, pct_special = 0'
    pct_special = 0._r8
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.false.)

    testname='no specific crop, sum(pctpft_full) = 100, pct_special = 0'
    pctpft_full = (/10._r8,  0._r8,   30._r8,  0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,&
                    0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   60._r8/)
    pct_special = 0._r8
    pct_natveg_t = 100._r8
    pct_crop_t = 0._r8
    pct_pft_t  = pctpft_full(:)
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.false.)

    testname='no specific crop, sum(pctpft_full) = 80, pct_special = 0'
    pctpft_full = (/10._r8,  0._r8,   10._r8,  0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,&
                    0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   60._r8/)
    pct_special = 0._r8
    pct_natveg_t = 80._r8
    pct_crop_t = 0._r8
    pct_pft_t  = pctpft_full(:) * 100._r8 / 80._r8
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.false.)

    testname='no specific crop, sum(pctpft_full) = 80, pct_special = 25'
    pctpft_full = (/0._r8,   10._r8,  0._r8,   30._r8,  0._r8,   0._r8,   0._r8,   0._r8,   0._r8,&
                    0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   40._r8,  0._r8/)
    pct_special = 25._r8
    pct_natveg_t = 80._r8 * 0.75_r8
    pct_crop_t = 0._r8
    pct_pft_t  = pctpft_full(:) * 100._r8 / 80._r8
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.false.)
    
    testname='no specific crop, sum(pctpft_full) = 80, pct_special = 100'
    pctpft_full = (/0._r8,   80._r8, 0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,&
                    0._r8,   0._r8,  0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8/)
    pct_special = 100._r8
    pct_natveg_t = 0._r8
    pct_crop_t = 0._r8
    pct_pft_t  = pctpft_full(:) * 100._r8 / 80._r8
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.false.)

    deallocate(pctpft_full, pct_pft, pct_cft, pct_pft_t, pct_cft_t)

    ! TESTS WITH SPECIFIC CROPS
    ! Set some module-level variables in mkpftMod
    numpft = 24
    num_natpft = 14
    num_cft = 10
    natpft_lb = 0
    natpft_ub = num_natpft
    cft_lb = num_natpft+1
    cft_ub = cft_lb + num_cft - 1
    allocate(pctpft_full(0:numpft), pct_pft(natpft_lb:natpft_ub), pct_cft(cft_lb:cft_ub))
    allocate(pct_pft_t(natpft_lb:natpft_ub), pct_cft_t(cft_lb:cft_ub))


    testname='specific crops, sum(pctpft_full) = 0, pct_special = 25'
    pctpft_full = (/0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8, &  ! natural
                    0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   &         ! more natural
                    0._r8,   0._r8, &  ! generic crops
                    0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8/)    ! specific crops
    pct_special = 0._r8
    pct_natveg_t = 0._r8
    pct_crop_t = 0._r8
    pct_pft_t  = (/100._r8, 0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8, &
                   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8 /)
    pct_cft_t  = (/100._r8, 0._r8,   &
                   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8/)
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.true.) 

    testname='specific crops, sum(pctpft_full) = 80 (all natveg), pct_special = 25'
    pctpft_full = (/0._r8,   40._r8,  0._r8,   0._r8,   30._r8,  0._r8,   0._r8,   0._r8, &  ! natural
                    0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   10._r8,  &         ! more natural
                    0._r8,   0._r8, &  ! generic crops
                    0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8/)    ! specific crops
    pct_special = 25._r8
    pct_natveg_t = 80._r8 * 0.75_r8
    pct_crop_t = 0._r8
    pct_pft_t  = (/0._r8,   50._r8,  0._r8,   0._r8,   37.5_r8, 0._r8,   0._r8,   0._r8, &
                   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   12.5_r8 /)
    pct_cft_t  = (/100._r8, 0._r8,   &
                   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8/)
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.true.)

    testname='specific crops, sum(pctpft_full) = 80 (all crop), pct_special = 25'
    pctpft_full = (/0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8, &  ! natural
                    0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   &         ! more natural
                    0._r8,   0._r8, &  ! generic crops
                    0._r8,   40._r8,  0._r8,   30._r8,  0._r8,   0._r8,   0._r8,   10._r8/)    ! specific crops
    pct_special = 25._r8
    pct_natveg_t = 0._r8
    pct_crop_t = 80._r8 * 0.75_r8
    pct_pft_t  = (/100._r8, 0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8, &
                   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8 /)
    pct_cft_t  = (/0._r8,   0._r8,   &
                   0._r8,   50._r8,  0._r8,   37.5_r8, 0._r8,   0._r8,   0._r8,   12.5_r8/)
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.true.)

    ! sum(pctpft_full) = 80 (60 natveg, 20 crop), pct_special = 25
    testname='specific crops, sum(pctpft_full) = 80 (60 natveg, 20 crop), pct_special = 25'
    pctpft_full = (/30._r8,  15._r8,  0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8, &  ! natural
                    0._r8,   15._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   &         ! more natural
                    5._r8,   10._r8, &  ! generic crops
                    0._r8,   0._r8,   5._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8/)    ! specific crops
    pct_special = 25._r8
    pct_natveg_t = 60._r8 * 0.75_r8
    pct_crop_t = 20._r8 * 0.75_r8
    pct_pft_t  = (/50._r8,  25._r8,  0._r8,   0._r8,   0._r8,   0._r8,   0._r8,   0._r8, &
                   0._r8,   25._r8,  0._r8,   0._r8,   0._r8,   0._r8,   0._r8 /)
    pct_cft_t  = (/25._r8,  50._r8,   &
                   0._r8,   0._r8,   25._r8,  0._r8,   0._r8,   0._r8,   0._r8,   0._r8/)
    call mkpft_normalize(pctpft_full, pct_special, pct_natveg, pct_crop, pct_pft, pct_cft)
    call check_results(.true.)

    deallocate(pctpft_full, pct_pft, pct_cft, pct_pft_t, pct_cft_t)

  contains
    subroutine check_results(include_crop)
      logical, intent(in) :: include_crop

      call test_close(pct_natveg, pct_natveg_t, eps, modname//' -- '//subname//' -- '//&
           trim(testname)//' -- pct_natveg')
      call test_close(pct_crop, pct_crop_t, eps, modname//' -- '//subname//' -- '//&
           trim(testname)//' -- pct_crop')
      call test_close(pct_pft, pct_pft_t, eps, modname//' -- '//subname//' -- '//&
           trim(testname)//' -- pct_pft')

      if (include_crop) then
         call test_close(pct_cft, pct_cft_t, eps, modname//' -- '//subname//' -- '//&
              trim(testname)//' -- pct_cft', rel_diff = .true.)
      end if
    end subroutine check_results

  end subroutine test_mkpft_normalize
!------------------------------------------------------------------------------

end module test_mkpftMod
