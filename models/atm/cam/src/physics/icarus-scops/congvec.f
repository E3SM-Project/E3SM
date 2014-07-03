
! *****************************COPYRIGHT****************************
! (c) British Crown Copyright 2009, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the
! following conditions are met:
! 
!     * Redistributions of source code must retain the above 
!       copyright  notice, this list of conditions and the following 
!       disclaimer.
!     * Redistributions in binary form must reproduce the above 
!       copyright notice, this list of conditions and the following 
!       disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its 
!       contributors may be used to endorse or promote products
!       derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
! 
! *****************************COPYRIGHT*******************************
! *****************************COPYRIGHT*******************************

      do irand = 1, npoints
          ! Marsaglia CONG algorithm
          seed(irand)=transfer(69069_i8*seed(irand)+1234567,1)
          ! mod 32 bit overflow
          seed(irand)=mod(seed(irand),2**30)   
          ran(irand)=seed(irand)*0.931322574615479E-09
      enddo

      ! convert to range 0-1 (32 bit only)
      overflow_32=transfer(i2_16*int(i2_16,i8),overflow_32)
      if ( overflow_32 .le. huge32 ) then
          do irand = 1, npoints
              ran(irand)=ran(irand)+1
              ran(irand)=(ran(irand))-int(ran(irand))
          enddo
      endif


