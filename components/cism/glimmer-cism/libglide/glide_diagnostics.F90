!TODO - May want to eliminate calculations of iarea, iareaf, and areag in calc_iareaf_iareag() and glide_set_mask().  
!       Instead just use the calculations made here.  These should be saved to the model derived type 
!        (model%geometry%iarea, etc.) for output.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_diagnostics.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glide_diagnostics

  ! subroutines for computing various useful diagnostics
  ! Author: William Lipscomb, LANL 
 
  use glimmer_global, only: dp
  use glimmer_log
  use glide_types

  implicit none

contains

  subroutine glide_write_diagnostics (model,  time,    &
                                      tstep_count,     &
                                      minthick_in)

    ! Short driver subroutine to decide whether it's time to write diagnostics.
    ! If so, it calls glide_write_diag.   

    ! input/output arguments

    type(glide_global_type), intent(in) :: model    ! model instance
    real(dp), intent(in)   :: time                  ! current time in years

    integer, intent(in), optional :: tstep_count    ! current timestep

    real(dp), intent(in), optional :: &
       minthick_in       ! ice thickness threshold (m) for including in diagnostics

    ! local arguments

    real(dp) :: minthick ! ice thickness threshold (m) for including in diagnostics
                         ! defaults to eps (a small number) if not passed in

    real(dp), parameter ::   &
       eps = 1.0d-11

    real(dp) ::   &
       quotient, nint_quotient

    if (present(minthick_in)) then
       minthick = minthick_in
    else
       minthick = eps  
    endif
 
! debug
!      print*, '	'
!      print*, 'In glide_write_diagnostics'
!      print*, 'time =', time
!      print*, 'dt_diag =', model%numerics%dt_diag
!      print*, 'ndiag =', model%numerics%ndiag
!      print*, 'tstep_count =', tstep_count

    !TODO - Make this method more robust (i.e., less prone to accumulated roundoff errors).
    !       Might want to derive ndiag from dt_diag at initialization.
    !       Then we would work with integers (tstep_count and ndiag) and avoid roundoff errors.

    if (model%numerics%dt_diag > 0.d0) then                            ! usual case

!!       if (mod(time,model%numerics%dt_diag)) < eps) then  ! not robust because of roundoff error

       quotient = time/model%numerics%dt_diag
       nint_quotient = nint(quotient)
       if (abs(quotient - real(nint_quotient,dp)) < eps) then  ! time to write

          call glide_write_diag(model,                 &
                                time,                  &
                                minthick,              &
                                model%numerics%idiag,  &
                                model%numerics%jdiag)
       endif

    elseif (present(tstep_count) .and. model%numerics%ndiag > 0) then  ! decide based on ndiag

       if (mod(tstep_count, model%numerics%ndiag) == 0)  then          ! time to write
          call glide_write_diag(model,                 &
                                time,                  &
                                minthick,              &
                                model%numerics%idiag,  &
                                model%numerics%jdiag)
       endif

    endif    ! dt_diag > 0

  end subroutine glide_write_diagnostics
 
!--------------------------------------------------------------------------

  subroutine glide_write_diag (model,       time,         &
                               minthick,                  &
                               idiag,       jdiag)

    ! Write global diagnostics
    ! Optionally, write local diagnostics for a selected grid cell
 
    use parallel

    !TODO - Remove scaling
    use glimmer_paramets, only: thk0, len0, vel0, tim0
    use glimmer_physcon, only: scyr, rhoi, shci
 
    implicit none
 
    ! input/output arguments

    type(glide_global_type), intent(in) :: model    ! model instance
    real(dp),  intent(in) :: time                   ! current time in years
    real(dp), intent(in)  :: &
         minthick          ! ice thickness threshold (m) for including in diagnostics

    integer, intent(in), optional :: &
         idiag, jdiag         ! i,j indices for diagnostics (on full grid)
                              ! indices will generally be different on local processor
 
    ! local arguments

    real(dp) ::                         &
         tot_area,                      &    ! total ice area (km^2)
         tot_volume,                    &    ! total ice volume (km^3)
         tot_energy,                    &    ! total ice energy (J)
         mean_thck,                     &    ! mean ice thickness (m)
         mean_temp,                     &    ! mean ice temperature (deg C)
         mean_acab,                     &    ! mean surface accumulation/ablation rate (m/yr)
         mean_bmlt,                     &    ! mean basal melt (m/yr)
         max_thck, max_thck_global,     &    ! max ice thickness (m)
         max_temp, max_temp_global,     &    ! max ice temperature (deg C)
         min_temp, min_temp_global,     &    ! min ice temperature (deg C)
         max_spd_sfc, max_spd_sfc_global,   &    ! max surface ice speed (m/yr)
         max_spd_bas, max_spd_bas_global,   &    ! max basal ice speed (m/yr)
         thck,                          &    ! thickness
         spd,                           &    ! speed
         thck_diag, usrf_diag,          &    ! local column diagnostics
         topg_diag, relx_diag,          &    
         artm_diag, acab_diag,          &
         bmlt_diag, bwat_diag,          &
         bheatflx_diag, level

    real(dp), dimension(model%general%upn) ::  &
         temp_diag,                     &    ! Note: sfc temp not included if temps are staggered
                                             !       (use artm instead)
         spd_diag

    real(dp), dimension(model%lithot%nlayer) ::  &
         lithtemp_diag                       ! lithosphere column diagnostics

    integer :: i, j, k, ktop, kbed,               &
               imax, imin,                        &
               jmax, jmin,                        &
               kmax, kmin,                        &
               imax_global, imin_global,          &
               jmax_global, jmin_global,          &
               kmax_global, kmin_global,          &
               ewn, nsn, upn,                     &    ! model%numerics%ewn, etc.
               nlith,                             &    ! model%lithot%nlayer
               velo_ew_ubound, velo_ns_ubound          ! upper bounds for velocity variables
 
    character(len=100) :: message
 
    real(dp), parameter ::   &
       eps = 1.0d-11,           &! small number
       unphys_val = -999999.d0   ! unphysical negative number
 
    integer ::   &
         global_row, global_col,    &! global row and column indices for diagnostic point
         idiag_local, jdiag_local,  &! local indices of diagnostic point 
         rdiag_local                 ! this_rank for diagnostic point

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn

    nlith = model%lithot%nlayer

    if (uhalo > 0) then
       velo_ns_ubound = nsn-uhalo
       velo_ew_ubound = ewn-uhalo
    else
       ! for uhalo==0 (as is the case for the glide dycore), the velocity grid has one less
       ! point than the main grid, so we need to subtract one to avoid out-of-bounds problems
       velo_ns_ubound = nsn-uhalo-1
       velo_ew_ubound = ewn-uhalo-1
    end if

!NOTE: Some of the global reductions below may seem unnecessary.
!      But at present, subroutine write_log permits writes only from main_task,
!       and the broadcast subroutines allow broadcasts only from main_task,
!       not from other processors.
!      So the way we get info to main_task is by parallel reductions.
!TODO: Support broadcasting from tasks other than main.

    !-----------------------------------------------------------------
    ! Determine whether global diagnostic point is on this processor.
    ! If so, communicate this information to the main processor.
    !-----------------------------------------------------------------

    if (present(idiag) .and. present(jdiag)) then

       rdiag_local = -999
       idiag_local = -999
       jdiag_local = -999

       if (idiag >= 1 .and. idiag <= global_ewn  &
                      .and.                             &
           jdiag >= 1 .and. jdiag <= global_nsn) then

          ! loop over gridcells owned by this processor
          do j = lhalo+1, nsn-uhalo
          do i = lhalo+1, ewn-uhalo
             global_row = (j - lhalo) + global_row_offset
             global_col = (i - lhalo) + global_col_offset
             if (global_col == idiag .and.   &
                 global_row == jdiag) then   ! diag point lives on this processor
                rdiag_local = this_rank
                idiag_local = i 
                jdiag_local = j 
             endif
          enddo   ! i
          enddo   ! j

       else
          call write_log('Error, global diagnostic point (idiag, jdiag) is out of bounds', GM_FATAL)
       endif      ! diagnostic point in bounds

       rdiag_local = parallel_reduce_max(rdiag_local)
       idiag_local = parallel_reduce_max(idiag_local)
       jdiag_local = parallel_reduce_max(jdiag_local)

    endif         ! present(idiag, jdiag)

    !-----------------------------------------------------------------
    ! Compute and write global diagnostics
    !-----------------------------------------------------------------
 
    call write_log('----------------------------------------------------------')
    call write_log(' ')
    write(message,'(a25,f24.16)') 'Diagnostic output, time =', time
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    call write_log(' ')

    ! total ice area (m^2)
 
    tot_area = 0.d0
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) * thk0 > minthick) then
             tot_area = tot_area + model%numerics%dew * model%numerics%dns
          endif
       enddo
    enddo
    tot_area = tot_area * len0**2
    tot_area = parallel_reduce_sum(tot_area)

    ! total ice volume (m^3)
 
    tot_volume = 0.d0
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) * thk0 > minthick) then
             tot_volume = tot_volume + model%geometry%thck(i,j)  &
                                     * model%numerics%dew * model%numerics%dns
          endif
       enddo
    enddo
    tot_volume = tot_volume * thk0 * len0**2
    tot_volume = parallel_reduce_sum(tot_volume)

    ! total ice energy relative to T = 0 deg C (J)
 
    tot_energy = 0.d0
    if (size(model%temper%temp,1) == upn+1) then  ! temps are staggered in vertical
       do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) * thk0 > minthick) then
             do k = 1, upn-1
                tot_energy = tot_energy +   &
                             model%geometry%thck(i,j) * model%temper%temp(k,i,j)    &
                            * model%numerics%dew * model%numerics%dns               &
                            *(model%numerics%sigma(k+1) - model%numerics%sigma(k))
             enddo
          endif
       enddo
       enddo
    
    else   ! temps are unstaggered in vertical
       do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) * thk0 > minthick) then
             ! upper half-layer, T = upper sfc temp
             tot_energy = tot_energy +   &
                          model%geometry%thck(i,j) * model%temper%temp(1,i,j)  &
                         * model%numerics%dew * model%numerics%dns             &
                         * 0.5d0 * model%numerics%sigma(2)
             do k = 2, upn-1
                tot_energy = tot_energy +   &
                             model%geometry%thck(i,j) * model%temper%temp(k,i,j) &
                           * model%numerics%dew * model%numerics%dns             &
                           * 0.5d0*(model%numerics%sigma(k+1) - model%numerics%sigma(k-1))
             enddo
             ! lower half-layer, T = lower sfc temp
             tot_energy = tot_energy +   &
                          model%geometry%thck(i,j) * model%temper%temp(upn,i,j)  &
                         * model%numerics%dew * model%numerics%dns               &
                         * 0.5d0 * (1.0d0 - model%numerics%sigma(upn-1))
          endif
       enddo
       enddo
    endif

    tot_energy = tot_energy * thk0 * len0**2 * rhoi * shci
    tot_energy = parallel_reduce_sum(tot_energy)

    ! mean thickness

    if (tot_area > eps) then
       mean_thck = tot_volume/tot_area
    else
       mean_thck = 0.d0
    endif

    ! mean temperature
 
    if (tot_volume > eps) then
       mean_temp = tot_energy/ (rhoi*shci*tot_volume)
    else
       mean_temp = 0.d0
    endif
 
    ! mean surface accumulation/ablation rate (m/yr)
 
    mean_acab = 0.d0
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) * thk0 > minthick) then
             mean_acab = mean_acab + model%climate%acab(i,j)  &
                                   * model%numerics%dew * model%numerics%dns
          endif
       enddo
    enddo
    mean_acab = mean_acab * scyr * thk0 / tim0 * len0**2   ! convert to m^3/yr
    mean_acab = parallel_reduce_sum(mean_acab)

    if (tot_area > eps) then
       mean_acab = mean_acab/tot_area    ! divide by total area to get m/yr
    else
       mean_acab = 0.d0
    endif

    ! mean basal melting rate (positive for ice loss)
 
    mean_bmlt = 0.d0
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) * thk0 > minthick) then
             mean_bmlt = mean_bmlt + model%temper%bmlt(i,j)  &
                                   * model%numerics%dew * model%numerics%dns
          endif
       enddo
    enddo

    mean_bmlt = mean_bmlt * scyr * thk0 / tim0 * len0**2   ! convert to m^3/yr
    mean_bmlt = parallel_reduce_sum(mean_bmlt)

    if (tot_area > eps) then
       mean_bmlt = mean_bmlt/tot_area    ! divide by total area to get m/yr
    else
       mean_bmlt = 0.d0
    endif

    ! write global sums and means

    write(message,'(a25,e24.16)') 'Total ice area (km^2)    ',   &
                                   tot_area*1.0d-6    ! convert to km^2
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Total ice volume (km^3)  ',   &
                                   tot_volume*1.0d-9  ! convert to km^3
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Total ice energy (J)     ', tot_energy
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,f24.16)') 'Mean thickness (m)       ', mean_thck
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,f24.16)') 'Mean temperature (C)     ', mean_temp
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Mean accum/ablat (m/yr)  ', mean_acab     
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Mean basal melt (m/yr)   ', mean_bmlt
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    
    ! Find various global maxes and mins

    ! max thickness

    imax = 0
    jmax = 0
    max_thck = -999.d0
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) > max_thck) then
             max_thck = model%geometry%thck(i,j)
             imax = i
             jmax = j
          endif
       enddo
    enddo

    imax_global = 0
    jmax_global = 0
    max_thck_global = parallel_reduce_max(max_thck)
    if (max_thck == max_thck_global) then  ! max_thck lives on this processor
       imax_global = (imax - lhalo) + global_col_offset
       jmax_global = (jmax - lhalo) + global_row_offset
    endif
    imax_global = parallel_reduce_max(imax_global)
    jmax_global = parallel_reduce_max(jmax_global)

    write(message,'(a25,f24.16,2i4)') 'Max thickness (m), i, j  ',   &
                                       max_thck_global*thk0, imax_global, jmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max temperature

    ktop = lbound(model%temper%temp,1)
    kbed = ubound(model%temper%temp,1)

    imax = 0
    jmax = 0
    kmax = 0
    max_temp =  -999.d0
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) * thk0 > minthick) then
             do k = ktop, kbed
                if (model%temper%temp(k,i,j) > max_temp) then
                   max_temp = model%temper%temp(k,i,j)
                   imax = i
                   jmax = j
                   kmax = k
                endif
             enddo
          endif
       enddo
    enddo

    imax_global = 0
    jmax_global = 0
    kmax_global = 0
    max_temp_global = parallel_reduce_max(max_temp)
    if (max_temp == max_temp_global) then  ! max_temp lives on this processor
       imax_global = (imax - lhalo) + global_col_offset
       jmax_global = (jmax - lhalo) + global_row_offset
       kmax_global = kmax
    endif
    imax_global = parallel_reduce_max(imax_global)
    jmax_global = parallel_reduce_max(jmax_global)
    kmax_global = parallel_reduce_max(kmax_global)

    ! min temperature

    imin = 0
    jmin = 0
    kmin = 0
    min_temp =  999.d0
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) * thk0 > minthick) then
             do k = ktop, kbed
                if (model%temper%temp(k,i,j) < min_temp) then
                   min_temp = model%temper%temp(k,i,j)
                   imin = i
                   jmin = j
                   kmin = k
                endif
             enddo
          endif
       enddo
    enddo

    imin_global = 0
    jmin_global = 0
    kmin_global = 0
    min_temp_global = -parallel_reduce_max(-min_temp)
    if (min_temp == min_temp_global) then  ! min_temp lives on this processor
       imin_global = (imin - lhalo) + global_col_offset
       jmin_global = (jmin - lhalo) + global_row_offset
       kmin_global = kmin
    endif
    imin_global = parallel_reduce_max(imin_global)
    jmin_global = parallel_reduce_max(jmin_global)
    kmin_global = parallel_reduce_max(kmin_global)
 
    write(message,'(a25,f24.16,3i4)') 'Max temperature, i, j, k ',   &
                    max_temp_global, imax_global, jmax_global, kmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    write(message,'(a25,f24.16,3i4)') 'Min temperature, i, j, k ',   &
                    min_temp_global, imin_global, jmin_global, kmin_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max surface speed

    imax = 0
    jmax = 0
    max_spd_sfc = -999.d0

    do j = lhalo+1, velo_ns_ubound
       do i = lhalo+1, velo_ew_ubound
          spd = sqrt(model%velocity%uvel(1,i,j)**2   &
                   + model%velocity%vvel(1,i,j)**2)
          if (model%geometry%thck(i,j) * thk0 > minthick .and. spd > max_spd_sfc) then
             max_spd_sfc = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo

    imax_global = 0
    jmax_global = 0
    max_spd_sfc_global = parallel_reduce_max(max_spd_sfc)
    if (max_spd_sfc == max_spd_sfc_global) then  ! max_spd_sfc lives on this processor
       imax_global = (imax - lhalo) + global_col_offset
       jmax_global = (jmax - lhalo) + global_row_offset
    endif
    imax_global = parallel_reduce_max(imax_global)
    jmax_global = parallel_reduce_max(jmax_global)

    write(message,'(a25,f24.16,2i4)') 'Max sfc spd (m/yr), i, j ',   &
                    max_spd_sfc_global*vel0*scyr, imax_global, jmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max basal speed

    imax = 0
    jmax = 0
    max_spd_bas = -999.d0
    do j = lhalo+1, velo_ns_ubound
       do i = lhalo+1, velo_ew_ubound
          spd = sqrt(model%velocity%uvel(upn,i,j)**2   &
                   + model%velocity%vvel(upn,i,j)**2)
          if (model%geometry%thck(i,j) * thk0 > minthick  .and. spd > max_spd_bas) then
             max_spd_bas = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo

    imax_global = 0
    jmax_global = 0
    max_spd_bas_global = parallel_reduce_max(max_spd_bas)
    if (max_spd_bas == max_spd_bas_global) then  ! max_spd_bas lives on this processor
       imax_global = (imax - lhalo) + global_col_offset
       jmax_global = (jmax - lhalo) + global_row_offset
    endif
    imax_global = parallel_reduce_max(imax_global)
    jmax_global = parallel_reduce_max(jmax_global)
 
    write(message,'(a25,f24.16,2i4)') 'Max base spd (m/yr), i, j',   &
                    max_spd_bas_global*vel0*scyr, imax_global, jmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! local diagnostics

    ! initialize to unphysical negative values
    usrf_diag     = unphys_val
    thck_diag     = unphys_val
    topg_diag     = unphys_val
    relx_diag     = unphys_val
    artm_diag     = unphys_val
    acab_diag     = unphys_val
    bmlt_diag     = unphys_val
    bwat_diag     = unphys_val
    bheatflx_diag = unphys_val
    temp_diag(:)  = unphys_val
    spd_diag (:)  = unphys_val
    lithtemp_diag(:) = unphys_val    

    if (present(idiag) .and. present(jdiag)) then

       ! Set local diagnostic values, and communicate them to main_task
       
       if (this_rank == rdiag_local) then

          i = idiag_local
          j = jdiag_local
          usrf_diag = model%geometry%usrf(i,j)*thk0
          thck_diag = model%geometry%thck(i,j)*thk0
          topg_diag = model%geometry%topg(i,j)*thk0
          relx_diag = model%isostasy%relx(i,j)*thk0
          artm_diag = model%climate%artm(i,j)
          acab_diag = model%climate%acab(i,j) * thk0*scyr/tim0
          bmlt_diag = model%temper%bmlt(i,j) * thk0*scyr/tim0
          bwat_diag = model%temper%bwat(i,j) * thk0
          bheatflx_diag = model%temper%bheatflx(i,j)
  
          temp_diag(:) = model%temper%temp(1:upn,i,j)          
          spd_diag(:) = sqrt(model%velocity%uvel(1:upn,i,j)**2   &
                           + model%velocity%vvel(1:upn,i,j)**2) * vel0*scyr
          if (model%options%gthf == GTHF_COMPUTE) &
             lithtemp_diag(:) = model%lithot%temp(i,j,:)
       endif

       usrf_diag = parallel_reduce_max(usrf_diag)
       thck_diag = parallel_reduce_max(thck_diag)
       topg_diag = parallel_reduce_max(topg_diag)
       relx_diag = parallel_reduce_max(relx_diag)
       artm_diag = parallel_reduce_max(artm_diag)
       acab_diag = parallel_reduce_max(acab_diag)
       bmlt_diag = parallel_reduce_max(bmlt_diag)
       bwat_diag = parallel_reduce_max(bwat_diag)
       bheatflx_diag = parallel_reduce_max(bheatflx_diag)

       do k = 1, upn
          temp_diag(k) = parallel_reduce_max(temp_diag(k))
          spd_diag(k)  = parallel_reduce_max(spd_diag(k))
       enddo

       do k = 1, nlith
          lithtemp_diag(k) = parallel_reduce_max(lithtemp_diag(k))
       enddo

       call write_log(' ')
       write(message,'(a39,2i4)')  &
            'Grid point diagnostics: (i,j) =', idiag, jdiag
       call write_log(trim(message), type = GM_DIAGNOSTIC)
       write(message,'(a39,3i4)')  &
            '                  Local (i,j,rank) =', idiag_local, jdiag_local, rdiag_local
       call write_log(trim(message), type = GM_DIAGNOSTIC)
       call write_log(' ')
 
       write(message,'(a25,f24.16)') 'Upper surface (m)        ', usrf_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Thickness (m)            ', thck_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Bedrock topo (m)         ', topg_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       if (model%options%isostasy == ISOSTASY_COMPUTE) then
          write(message,'(a25,f24.16)') 'Relaxed bedrock (m)   ', relx_diag
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       write(message,'(a25,f24.16)') 'Sfc mass balance (m/yr)  ', acab_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Basal melt rate (m/yr)   ', bmlt_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Basal water depth (m)    ', bwat_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Basal heat flux (W/m^2)  ', bheatflx_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)
 
       ! Vertical profile of ice speed and temperature

       call write_log(' ')
       write(message,'(a55)') ' Sigma       Ice speed (m/yr)       Ice temperature (C)'
       call write_log(trim(message), type = GM_DIAGNOSTIC)
 
       if (size(model%temper%temp,1) == upn+1) then   ! temperatures staggered in vertical
                                                      ! (at layer midpoints)

           ! upper surface 
           write (message,'(f6.3,2f24.16)') model%numerics%sigma(1), spd_diag(1), artm_diag
           call write_log(trim(message), type = GM_DIAGNOSTIC)

           ! internal
           do k = 1, upn-1

              ! speed at top of layer
              if (k > 1) then
                 write (message,'(f6.3,f24.16)') model%numerics%sigma(k), spd_diag(k)
                 call write_log(trim(message), type = GM_DIAGNOSTIC)
              endif

              ! temp at layer midpoint
              write (message,'(f6.3,24x,f24.16)') model%numerics%stagsigma(k), temp_diag(k)
              call write_log(trim(message), type = GM_DIAGNOSTIC)

           enddo

           ! lower surface
           write (message,'(f6.3,2f24.16)') model%numerics%sigma(upn), spd_diag(upn), temp_diag(upn)
           call write_log(trim(message), type = GM_DIAGNOSTIC)
           
       else    ! temperatures unstaggered in vertical (at layer interfaces)
 
           do k = 1, upn
             write (message,'(f6.3,2f24.16)') model%numerics%sigma(k), spd_diag(k), temp_diag(k)
             call write_log(trim(message), type = GM_DIAGNOSTIC)
          enddo

       endif  ! temps staggered

       ! Vertical profile of upper lithosphere temperature

       if (model%options%gthf == GTHF_COMPUTE) then

          call write_log(' ')
          write(message,'(a41)') '  Level (m)          Lithosphere temp (C)'
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          level = 0.d0
          do k = 1, nlith
             level = level + model%lithot%deltaz(nlith)
             write (message,'(f10.0,6x,f24.16)') level, lithtemp_diag(k)
             call write_log(trim(message), type = GM_DIAGNOSTIC)
          enddo

       endif  ! gthf_compute

    endif     ! idiag and jdiag present

    call write_log(' ')

  end subroutine glide_write_diag
     
!==============================================================

end module glide_diagnostics
