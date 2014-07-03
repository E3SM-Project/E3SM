! (c) 2008, Lawrence Livermore National Security Limited Liability Corporation.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list 
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Lawrence Livermore National Security Limited Liability Corporation 
!       nor the names of its contributors may be used to endorse or promote products derived from 
!       this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
      
      subroutine pf_to_mr(npoints,nlev,ncol,rain_ls,snow_ls,grpl_ls,
     &                    rain_cv,snow_cv,prec_frac,
     &                    p,t,mx_rain_ls,mx_snow_ls,mx_grpl_ls,
     &                    mx_rain_cv,mx_snow_cv)


      implicit none

      INTEGER npoints       !  number of model points in the horizontal
      INTEGER nlev          !  number of model levels in column
      INTEGER ncol          !  number of subcolumns

      INTEGER j,ilev,ibox
      
      REAL rain_ls(npoints,nlev),snow_ls(npoints,nlev) ! large-scale precip. flux
      REAL grpl_ls(npoints,nlev)
      REAL rain_cv(npoints,nlev),snow_cv(npoints,nlev) ! convective precip. flux

      REAL prec_frac(npoints,ncol,nlev) ! 0 -> clear sky
                                        ! 1 -> LS precipitation
                                        ! 2 -> CONV precipitation
                                        ! 3 -> both
      REAL mx_rain_ls(npoints,ncol,nlev),mx_snow_ls(npoints,ncol,nlev)
      REAL mx_grpl_ls(npoints,ncol,nlev)
      REAL mx_rain_cv(npoints,ncol,nlev),mx_snow_cv(npoints,ncol,nlev)
      REAL p(npoints,nlev),t(npoints,nlev)
      REAL ar,as,ag,br,bs,bg,nr,ns,ng,rho0,rhor,rhos,rhog,rho
      REAL term1r,term1s,term1g,term2r,term2s,term2g,term3
      REAL term4r_ls,term4s_ls,term4g_ls,term4r_cv,term4s_cv
      REAL term1x2r,term1x2s,term1x2g,t123r,t123s,t123g
      
      ! method from Khairoutdinov and Randall (2003 JAS)

      ! --- List of constants from Appendix B
      ! Constant in fall speed formula
      ar=842.
      as=4.84
      ag=94.5
      ! Exponent in fall speed formula
      br=0.8
      bs=0.25
      bg=0.5
      ! Intercept parameter
      nr=8.*1000.*1000.
      ns=3.*1000.*1000.
      ng=4.*1000.*1000.
      ! Densities for air and hydrometeors
      rho0=1.29
      rhor=1000.
      rhos=100.
      rhog=400.
      ! Term 1 of Eq. (A19).
      term1r=ar*17.8379/6.
      term1s=as*8.28508/6.
      term1g=ag*11.6317/6.
      ! Term 2 of Eq. (A19).
      term2r=(3.14159265*rhor*nr)**(-br/4.)
      term2s=(3.14159265*rhos*ns)**(-bs/4.)
      term2g=(3.14159265*rhog*ng)**(-bg/4.)
      
      term1x2r=term1r*term2r
      term1x2s=term1s*term2s
      term1x2g=term1g*term2g
      do ilev=1,nlev
        do j=1,npoints
            rho=p(j,ilev)/(287.05*t(j,ilev))
            term3=(rho0/rho)**0.5
            ! Term 4 of Eq. (A19).
            t123r=term1x2r*term3
            t123s=term1x2s*term3
            t123g=term1x2g*term3
            term4r_ls=rain_ls(j,ilev)/(t123r)
            term4s_ls=snow_ls(j,ilev)/(t123s)
            term4g_ls=grpl_ls(j,ilev)/(t123g)
            term4r_cv=rain_cv(j,ilev)/(t123r)
            term4s_cv=snow_cv(j,ilev)/(t123s)
            do ibox=1,ncol
                mx_rain_ls(j,ibox,ilev)=0.
                mx_snow_ls(j,ibox,ilev)=0.
                mx_grpl_ls(j,ibox,ilev)=0.
                mx_rain_cv(j,ibox,ilev)=0.
                mx_snow_cv(j,ibox,ilev)=0.
                if ((prec_frac(j,ibox,ilev) .eq. 1.) .or.
     &              (prec_frac(j,ibox,ilev) .eq. 3.)) then 
                    mx_rain_ls(j,ibox,ilev)=
     &                     (term4r_ls**(1./(1.+br/4.)))/rho
                    mx_snow_ls(j,ibox,ilev)=
     &                     (term4s_ls**(1./(1.+bs/4.)))/rho
                    mx_grpl_ls(j,ibox,ilev)=
     &                     (term4g_ls**(1./(1.+bg/4.)))/rho
                endif
                if ((prec_frac(j,ibox,ilev) .eq. 2.) .or.
     &              (prec_frac(j,ibox,ilev) .eq. 3.)) then 
                    mx_rain_cv(j,ibox,ilev)=
     &                     (term4r_cv**(1./(1.+br/4.)))/rho
                    mx_snow_cv(j,ibox,ilev)=
     &                     (term4s_cv**(1./(1.+bs/4.)))/rho
                endif
            enddo ! loop over ncol
        enddo ! loop over npoints
      enddo ! loop over nlev
  
      end

