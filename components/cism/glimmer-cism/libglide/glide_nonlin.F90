!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_nonlin.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!Contains helper functions for nonlinear iteration, both to embed in the
!iteration loop and to serialize the data into the vector format that these
!functions require.
!Currently only unstable manifold correction is implemented.

module glide_nonlin

    use glimmer_global, only: dp
    use glimmer_physcon, only: pi
    implicit none

contains

    subroutine check_vector_size(start, veclen, ni, nj, nk)
       use glimmer_log
       integer :: start, veclen, ni, nj, nk
       character(256) :: message
       if (ni*nj*nk > veclen - start + 1) then
           write(message, *) "Need ",ni*nj*nk," elements in vector, starting from element ",start," only have ",veclen - start+1
           call write_log(message, GM_FATAL)
       end if
    end subroutine


    subroutine linearize_3d(vector, start, field)
        use glimmer_paramets, only: GLC_DEBUG
        real(dp), dimension(:) :: vector
        integer :: start
        real(dp), dimension(:,:,:) :: field
        integer :: ni, nj, nk
        integer :: i,j,k
        
        ni = size(field, 1)
        nj = size(field, 2)
        nk = size(field, 3)
        if (GLC_DEBUG) then
           call check_vector_size(start, size(vector), ni, nj, nk)
        end if
        do i=1,ni
            do j=1,nj
                do k=1,nk
                    vector(start) = field(i,j,k)
                    start = start + 1
                end do
            end do
        end do
    end subroutine

    subroutine linearize_2d(vector, start, field)
        use glimmer_paramets, only: GLC_DEBUG
        real(dp), dimension(:) :: vector
        integer :: start
        real(dp), dimension(:,:) :: field
        integer :: ni, nj
        integer :: i,j

        ni = size(field, 1)
        nj = size(field, 2)
        if (GLC_DEBUG) then
           call check_vector_size(start, size(vector), ni, nj, 1)
        end if
        do i=1,ni
            do j=1,nj
                vector(start) = field(i,j)
                start = start + 1
            end do
        end do
    end subroutine
    
    subroutine delinearize_3d(vector, start, field)
        real(dp), dimension(:) :: vector
        integer :: start
        real(dp), dimension(:,:,:) :: field
        integer :: ni, nj, nk
        integer :: i,j,k

        ni = size(field, 1)
        nj = size(field, 2)
        nk = size(field, 3)
        
        do i=1,ni
            do j=1,nj
                do k=1,nk
                    field(i,j,k) = vector(start)
                    start = start + 1
                end do
            end do
        end do
    end subroutine

    subroutine delinearize_2d(vector, start, field)
        real(dp), dimension(:) :: vector
        integer :: start
        real(dp), dimension(:,:) :: field
        integer :: ni, nj
        integer :: i,j

        ni = size(field, 1)
        nj = size(field, 2)
        
        do i=1,ni
            do j=1,nj
                field(i,j) = vector(start)
                start = start + 1
            end do
        end do
    end subroutine

    function picard_iterate(vec_new, vec_old, vec_size, toler, tot_out)
        logical :: picard_iterate

        real(dp), dimension(:), intent(in) :: vec_new
        real(dp), dimension(:), intent(inout) :: vec_old
        integer :: vec_size
        real(dp) :: toler
        real(dp), optional, intent(out) :: tot_out

        real(dp) :: err, norm1, norm2

        norm1 = sqrt(sum(vec_new**2))
        norm2 = sqrt(sum((vec_new-vec_old)**2))

        err = norm2/(norm1 + 1d-10)
        picard_iterate = err >= toler

        vec_old = vec_new

        if (present(tot_out)) then
            tot_out = err
        end if
    end function picard_iterate

    function unstable_manifold_correction(vec_new, vec_old, vec_correction, &
                                          vec_size, toler, tot_out, theta_out)
        logical :: unstable_manifold_correction
        
        real(dp), dimension(:), intent(in) :: vec_new
        real(dp), dimension(:), intent(inout) :: vec_old
        real(dp), dimension(:), intent(inout) :: vec_correction
        integer :: vec_size
        real(dp) :: toler
        real(dp), optional, intent(out) :: tot_out
        real(dp), optional, intent(out) :: theta_out

        real(dp) :: norm1, norm2, norm3, norm4, norm5
        real(dp) :: tot
        real(dp) :: theta
        real(dp) :: alpha
        integer :: i
        real(dp) :: vmean
        real(dp) :: vstd

        real(dp), dimension(vec_size) :: vec_correction_new

        !Assume we need to iterate again until proven otherwise
        unstable_manifold_correction = .true.

        norm1 = 0.d0
        norm2 = 0.d0
        norm3 = 0.d0
        norm4 = 0.d0
        norm5 = 0.d0

        vec_correction_new = vec_new(1:vec_size) - vec_old(1:vec_size)

        do i = 1, vec_size
            vmean = vmean + abs(vec_correction_new(i))
        end do
        vmean = vmean / vec_size

        do i = 1, vec_size
            vstd = vstd + (vec_correction_new(i) - vmean)**2
        end do
        vstd = sqrt(vstd/vec_size)

        do i = 1,vec_size
            norm1 = norm1 + (vec_correction_new(i) - vec_correction(i)) ** 2
            norm2 = norm2 + vec_correction(i) ** 2
            !if (abs(vec_correction_new(i)) > vmean * 4. * vstd) then
            !else
            norm3 = norm3 + vec_correction_new(i) ** 2
            !endif
            norm4 = norm4 + vec_correction(i) * vec_correction_new(i)
            norm5 = norm5 + vec_new(i) ** 2
        end do

        !TODO - Change PI to pi.
        !Compute the angle between successive correction vectors
        if ((abs(norm2) < 1d-10) .or. (abs(norm3) < 1d-10)) then
            theta=PI/2.
        else
            theta=acos(norm4/sqrt(norm2*norm3))
        endif

        if ( (theta <= (5.*PI/6.) ) ) then
            !We've requested unstable manifold correction, and the angle is
            !small (less than 5pi/6, a value identified by Hindmarsh and Payne
            !to work well).   If this is the case, we compute and apply
            !a correction vector.
            
            !Compute the error between the last two *correction vectors* (not
            !the last two iteration values!)  (See (51) in Pattyn's paper)
            if (abs(norm2) > 0.) then !We're just avoiding a divide by 0 here
                alpha=sqrt(norm1/norm2)
            else
                alpha=1.
            endif
            
            if (alpha < 1.e-6) then
                !If the correction vector didn't change much, we're done
                unstable_manifold_correction = .false.
            else
                !Update the previous guess of the velocity with the correction
                !vector.  This throws out the current iteration's computed
                !velocity, and instead uses the computed correction vector.
                vec_old = vec_old + vec_correction_new / alpha
                vec_correction = vec_correction_new

            endif
        else
            !Copy this iteration's new values to the old values
            !for the next iteration - because the angle between correction
            !vectors is large we do not want to apply a correction, so
            !we just go Picard instead
            vec_old = vec_new
            vec_correction = vec_correction_new
        endif
        
        tot=sqrt(norm3/(norm5+1d-10)) !Regularize the denominator so we don't get NAN with simple geometries

        if (present(tot_out)) then
            tot_out = tot
        end if
        
        if (present(theta_out)) then
            theta_out = theta * 180 / pi
        end if

        if (tot < toler) unstable_manifold_correction = .false. 
    end function unstable_manifold_correction
end module glide_nonlin
