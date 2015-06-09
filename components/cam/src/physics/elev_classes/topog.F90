MODULE topog

  use shr_kind_mod,   only: r8 => shr_kind_r8, SHR_KIND_CL
  use physics_types,  only: physics_state, physics_ptend

  implicit none
  private
  save

  public :: topog_tend

CONTAINS

  subroutine topog_tend(state, ptend, dt)

!---- determine the height rise of an air parcel in transition from grid
!---- cell mean elevation to elevation of different classes. interpolate 
!---- thermodynamic variables to sigma coordinates at different classes.
!---- calculate orographic tendency from profiles.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM3.6.7
! Author:            Steve Ghan, April 2000
! Applied to ACME    Steve Ghan, June 2015
!
!-----------------------------------------------------------------------
    use constituents,   only: pcnst, orotendname
    use physconst,      only: rair, gravit, cpair, latvap, cappa, zvir
    use geopotential,   only: geopotential_dse
    use dimensions_mod, only: nelemd
    use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
    use physics_types,  only: physics_ptend_init, physics_update
    use elev_classes,   only: max_elevation_classes

!
!     i/o fields
!

!     subgrid fields

    real(r8),            intent(in)    :: dt
    type(physics_state), intent(inout) :: state
    type(physics_ptend), intent(out)   :: ptend      

!     local fields
!     grid cell mean fields

    type(physics_state), pointer :: stateg

!     grid cell mean

    real(r8) :: gthx(pverp) ! potential temperature
    real(r8) :: gbvf(pver) ! brunt-vaisala frequency
    real(r8) :: gpmid(pver) ! pressure at model levels
    real(r8) :: gpint(pverp)   ! Interface pressures
    real(r8) :: glnpint(pverp) ! ln(pressure) at model interfaces
    real(r8) :: gpdel(pver) ! pdel(k)   = pint  (k+1)-pint  (k)
    real(r8) :: gzs ! surface elevation
    real(r8) :: gzint(pverp) ! height at model interfaces
    real(r8) :: gzmid(pver)  ! height at model levels
    real(r8) :: ghmax ! maximum height rise (m)
    real(r8) :: gbvfmin ! minimum brunt-vaisala frequency
    real(r8) :: gwspeed ! wind speed
    real(r8) :: gsvht
    real(r8) :: gcmass ! column mass

!     subgrid

    real(r8), allocatable :: s(:,:) ! dry static energy at model levels
    real(r8), allocatable :: pmid(:,:) ! pressure at model levels
    real(r8), allocatable :: lnpmid(:,:) ! ln(pressure) at model levels
    real(r8), allocatable :: pint(:,:)   ! Interface pressures
    real(r8), allocatable :: lnpint(:,:) ! ln(pressure) at model interfaces
    real(r8), allocatable :: pdel(:,:) ! pdel(k)   = pint  (k+1)-pint  (k)
    real(r8), allocatable :: rpdel(:,:) ! rpdel(k)   = 1/pdel(k)
    real(r8), allocatable :: zi(:,:),zm(:,:)
    real(r8) :: zs ! surface elevation
    real(r8) :: htdiff
    real(r8) :: zblk ! depth of blocking layer
    real(r8) :: hrise ! height rise of air parcels
    real(r8) :: zgrid(pver)
    real(r8) :: zint(pverp) ! height at model interfaces
    real(r8) :: exner(pver) ! conversion from potential temp to temp
    real(r8) :: zeclass(max_elevation_classes,pver)  ! height at model levels
    integer  :: kgh(pver) ! interpolation index
    integer  :: kg ! kgh
    real(r8) :: fact(pver) ! interpolation factor
    real(r8) :: uf(max_elevation_classes,pver) ! orographic profile of u wind
    real(r8) :: vf(max_elevation_classes,pver) ! orographic profile of v wind
    real(r8) :: tf(max_elevation_classes,pver) ! orographic profile of temperature
    real(r8) :: qf(max_elevation_classes,pver,pcnst) ! orographic profile of constituents
    real(r8) :: normt(pver) ! normalization factor
    real(r8) :: normq(pver,pcnst) ! normalization factor

    integer  :: k, kk  ! layer index
    integer  :: ig ! grid cell index
    integer  :: ih, index ! column index
    integer  :: n  ! constituent index
    integer  :: istat
    real(r8) :: rovg  ! rair/gravit
    real(r8) :: work(pcols,pver)

    real(r8), parameter :: frc = 1.0_r8 ! critical froude number
    real(r8), parameter :: timescale = 10._r8 ! orographic adjustment timescale (hours)
    real(r8), parameter :: rate = 1._r8 / (timescale*3600._r8) ! orographic adjustment rate
    real(r8) :: qold
    logical  :: lq(pcnst)

    nullify(stateg)

    allocate (s(pcols,pver), stat=istat)
    allocate (pmid(pcols,pver), stat=istat)
    allocate (lnpmid(pcols,pver), stat=istat)
    allocate (pint(pcols,pverp), stat=istat)
    allocate (lnpint(pcols,pverp), stat=istat)
    allocate (pdel(pcols,pver), stat=istat)
    allocate (rpdel(pcols,pver), stat=istat)
    allocate (zi(pcols,pverp), stat=istat)
    allocate (zm(pcols,pver), stat=istat)

    lq(:) = .TRUE.
    lq(1) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'topog', ls=.true., lq=lq)! initialize local ptend type

!!XXgoldyXX!! v Fix this!!!
#if 0
    do ie = 1, nelemd  ! Loop over elements 

      index = 0
      do ig = 1, ncol ! Loop over physics columns in an element

        stateg => state%grid_phys_state(ig)
        if (associated(stateg)) then
!        if(phys_columns(i, ie)%max_elevation_classes.eq.1)then

!        no subgrid classes. set subgrid values to grid mean values

          index = index + 1

          state%t(index,:,ie) = stateg%t(ig,:,ie)
          state%u(index,:v) = stateg%u(ig,:,ie)
          state%v(index,:,ie) = stateg%v(ig,:,ie)
          state%q(index,:,:,ie) = stateg%q(ig,:,:,ie)

        else

          do k = 1, pverp
            gpint(k) = hyai(k)*ps0 + hybi(k)*stateg%ps(ig,ie)
            glnpint(k) = log(gpint(k))
          end do
          do k = 1, pver
            gpdel(k) = gpint(k+1) - gpint(k)
            gpmid(k) = 0.5_r8*(gpint(k) + gpint(k+1))
          end do
!
!---- determine the brunt-vaisala frequency
!
          do k = 1, pver
            gthx(k)=stateg%t(ig,k,ie)*(stateg%ps(ig,ie)/gpmid(k))**cappa
          end do
!
          gthx(pver+1) = stateg%t(ig,pver,ie)
!
          gzs = stateg%phis(ig,ie) / gravit
          gzint(pver+1) = gzs
!
          rovg=rair/gravit
          do k = pver, 1, -1
            gzint(k) = gzint(k+1)+                             &
                 rovg*stateg%t(ig,k,ie)*(1. + zvir * stateg%q(ig,k,1,ie))*(glnpint(k+1)-glnpint(k))
          end do
!
          do k = 1, pver
            gzmid(k) = 0.5_r8*(gzint(k)+gzint(k+1))
          end do
          gzmid(pver+1) = gzint(pver+1)
!
!---- maximum height rise of an air parcel
!
          gbvfmin = 1.e5_r8
          do k = pver-5, pver-1
            gbvf(k) = gravit*(gthx(k)-gthx(k+1))/         &
                 (0.5_r8*(gthx(k)+gthx(k+1))*(gzmid(k)-gzmid(k+1)))
            gbvfmin=min(gbvf(k),gbvfmin)
          end do
          if(gbvfmin <= 0.0_r8) then
            ghmax = 1.e8_r8
          else
            gsvht = 0.0_r8
            gcmass = 0.0_r8
            do k = pver-5, pver-1
              gwspeed=stateg%u(ig,k)*stateg%u(ig,k,ie)+stateg%v(ig,k,ie)*stateg%v(ig,k,ie)
              gsvht=gsvht+gpdel(ig,k,ie)*sqrt(gwspeed/gbvf(k))
              gcmass=gcmass+gpdel(k)
            end do
            gsvht=gsvht*frc/gcmass
            ghmax=gsvht
          end if

          normt(:) = 0.0_r8
          normq(:,:) = 0.0_r8
          indexsave = index

          do ic = 1, phys_columns(ig, ie)%max_elevation_classes

            index = index + 1
            ih=index
!
!         height parcel rises to, for each elevation class
!
            zs = state%phis(ih)/gravit
            htdiff = zs-gzs
            zblk = htdiff-ghmax+gzs
            if(zblk <= gzs)then
              do k = 1, pver
                zgrid(k) = gzmid(k)+htdiff
              end do
            else
              do k  = 1, pver
                if(gzmid(k) < zblk) then
                  hrise = ghmax(ig)*(gzmid(k)-gzs)/(zblk-gzs))
                else
                  hrise = ghmax
                end if
                zgrid(k)=gzmid(k)+min(hrise,htdiff)
              end do
            end if

!        now calculate z for all sigma levels in each class
!        using grid-cell mean as first guess for the potential temperature profile

            zint(pver+1) = zs(ih)
!
            do k = 1, pverp
              pint(ih,k) = hyai(k)*ps0 + hybi(k)*state%ps(index,ie)
              lnpint(ih,k) = log(pint(ih,k))
            end do
            do k = 1, pver
              pdel(ih,k) = pint(ih,k+1) - pint(ih,k)
              rpdel(ih,k) = 1._r8 / pdel(ih,k)
              pmid(ih,k) = 0.5_r8*(pint(ih,k) + pint(ih,k+1))
              lnpmid(ih,k) = log(pmid(ih,k))
            end do

            do k = pver, 1, -1
              exner(k)=(pmid(ih,k)/stateg%ps(ig,ie))**cappa
              zint(k)=zint(k+1)+rovg*gthx(k)*exner(k)*(1. + zvir * stateg%q(ig,k,1,ie))    &
                   *(lnpint(ih,k+1)-lnpint(ih,k))
            end do
!
            do k = 1, pver
              zeclass(ic,k)=0.5_r8 * (zint(k)+zint(k+1))
            end do
!
!---- determine the k values for interpolating from grid to class
!
            kgh = -1
            do k = 1, pver
              do kk = 2, pver
                if(zgrid(kk) < zeclass(ic,k)) then
                  kgh(k) = kk
                  exit
                end if
              end do
              if (kgh(k) < 1) then
                kgh(k) = pver
              end if
            end do
!
!---- interpolate variables to subgrid levels
!
            do k = 1, pver
              kg = kgh(k)
              fact(k) = (zeclass(k)-zgrid(kg))/(zgrid(kg-1)-zgrid(kg))
              fact(ih,k) = max(fact(k),0._r8)
              fact(ih,k) = min(fact(k),1._r8)
              uf(ic,k) = stateg%u(ig,k,ie) ! ensures energy conservation
              vf(ic,k) = stateg%v(ig,k,ie) ! ensures energy conservation
              tf(ic,k) = (fact(k)*(gthx(kg-1) - gthx(kg)) + gthx(kg)) *exner(k)
            end do
!
!
!       transformation conserves grid cell means
!

            normt(:) = normt(:) + phys_columns(ig, ie)%area(ih) * pdel(ih,:) * tf(ic,:)
            do n = 1, pcnst
              do k = 1, pver
                qf(ic,k,n) = fact(k)*(stateg%q(ig,kg-1,n)-stateg%q(ig,kg,n)) + stateg%q(ig,kg,n)
                normq(k,n) = normq(k,n)+ phys_columns(ig, ie)%area(ih) * pdel(ih,k) * qf(ic,k,n)
              end do
            end do

          end do  ! ic

!     Add orographic change
          index = indexsave
          do ic = 1, phys_columns(ig, ie)%max_elevation_classes
            index = index + 1
            ih = index

            tf(ic,:) = tf(ic,:) * gpdel(:) * stateg%t(ig,:,ie) / normt(:)

            do n = 1, pcnst
              do k = 1, pver
                if(normq(k,n) > 1.e-20_r8)then
                  qf(ic,k,n) = qf(ic,k,n) * gpdel(k) * stateg%q(ig,k,n,ie) / normq(k,n)
                end if
                ptend%q(ih,k,n) = (qf(ic,k,n)-state%q(ih,k,n,ie)) * rate
                state%q(ih,k,n) = state%q(ih,k,n,ie) + ptend%q(ih,k,n) * dt
              end do
            end do
            state%u(ih,:) = uf(ih,:)
            state%v(ih,:) = vf(ih,:)
            ptend%s(ih,k) = gravit*zeclass(ic,k) + cpair*(tf(ic,:)-state%t(ih,:,ie)) * rate
            s(ih,:) = cpair*state%t(ih,:,ie) + ptend%s(ih,k) * dt
      
          end do

        end if

      end do   ! ig

    end do ! ie

    work(:ncol,:) = ptend%s(:ncol,:) / cpair
    call outfld('OROTTEND', work, pcols, lchnk   )
    do m = 1, pcnst
      call outfld(orotendname(m), ptend%q(1,1,m), pcols, lchnk   )
    end do
#endif
!!XXgoldyXX!! ^ Fix this!!!

    deallocate (s)
    deallocate (pmid)
    deallocate (lnpmid)
    deallocate (pint)
    deallocate (lnpint)
    deallocate (pdel)
    deallocate (rpdel)
    deallocate (zi)
    deallocate (zm)

    return
  end subroutine topog_tend

END MODULE topog

