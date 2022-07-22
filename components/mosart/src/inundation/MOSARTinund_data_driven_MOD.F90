MODULE MOSARTinund_data_driven_MOD
!------------------------------------------------------------------------------!
! DESCRIPTION: MOSART inundation with linear log relationship.
!              FF = a * log(Vtotal) + b
! 
!------------------------------------------------------------------------------!

    use shr_kind_mod, only: r8 => shr_kind_r8
    use shr_sys_mod, only: shr_sys_abort
    use RunoffMod, only: rtmCTL, Tctl, TUnit, TRunoff

    implicit none
    private

    public inundation_run

    contains

    subroutine inundation_run(iunit)

        implicit none
        integer, intent(in)  :: iunit
        real(r8)             :: Vtotal

        if (TUnit%linear_a(iunit) > -999._r8) then
            Vtotal = TRunoff%wr(iunit,1) + TRunoff%wf_ini(iunit)
            TRunoff%ff_unit(iunit) = calculate_flooded_fraction(Vtotal,TUnit%linear_a(iunit),TUnit%linear_b(iunit),TUnit%linear_vcri(iunit))
            if (TRunoff%ff_unit(iunit) > TUnit%a_chnl(iunit) ) then
                TRunoff%ff_fp(iunit) = TRunoff%ff_unit(iunit) - TUnit%a_chnl(iunit) 
            else
                TRunoff%ff_fp(iunit) = 0._r8
            endif
        else
            TRunoff%ff_fp(iunit) = 0._r8
            TRunoff%ff_fp(iunit) = 0._r8
        endif

        if (TRunoff%ff_fp(iunit) > 0._r8) then
            call calculate_volume_exchange(TRunoff%wr(iunit,1),TRunoff%wf_ini(iunit), &
                                           TUnit%wr_bf(iunit),TRunoff%ff_unit(iunit), &
                                           TRunoff%ff_fp(iunit))
            TRunoff%yr(iunit,1) = TRunoff%wr(iunit,1) / TUnit%rwidth(iunit) / TUnit%rlen(iunit)
        endif

        TRunoff%ff_ini(iunit)     = TRunoff%ff_fp(iunit)
        TRunoff%ffunit_ini(iunit) = TRunoff%ff_unit(iunit)

    end subroutine inundation_run

    real (r8) function calculate_flooded_fraction(Vtotal,a,b,Vcri)

        implicit none
        real(r8), intent(in) :: Vtotal ! total volume [m^{3}], Vchannel + Vfloodplain
        real(r8), intent(in) :: a, b   ! FF = a * log(Vtotal) + b
        real(r8), intent(in) :: Vcri   ! critical volume, FF = 0 when Vtotal <= Vcri
        character( len = * ), parameter :: subname = '(calculate_flooded_fraction)'

        if (Vtotal < 1._r8 .or. Vtotal <= Vcri) then
            calculate_flooded_fraction = 0._r8
        else
            calculate_flooded_fraction = a*log(Vtotal) + b
            if (calculate_flooded_fraction > 1._r8) then
                calculate_flooded_fraction= 1._r8
            elseif (calculate_flooded_fraction < 0._r8) then
                calculate_flooded_fraction = 0._r8
            endif
        endif

        return

    end function calculate_flooded_fraction

    subroutine calculate_volume_exchange(wr,wf,wr_bf,ff_unit,ff_fp)

        implicit none
        real(r8), intent(inout) :: wr, wf
        real(r8), intent(in)    :: wr_bf, ff_unit, ff_fp
        real(r8)                :: w_over     ! channel storage + floodplain storage - channel storage capacity (m^3).
        real(r8)                :: hbf_excess ! bankfull excess water depth
        character( len = * ), parameter :: subname = '(calculate_volume_exchange)'

        w_over  = wr + wf - wr_bf
        if (w_over > 0) then
            hbf_excess = w_over / ff_unit
            wr         = wr_bf + hbf_excess * (ff_unit - ff_fp)
            wf         = hbf_excess * ff_fp
        else
            wr = wr + wf
            wf = 0._r8
        endif

    end subroutine

end MODULE MOSARTinund_data_driven_MOD