!---------------------------------------------------------------
!	... tgcm upper bndy values
!---------------------------------------------------------------

      module mo_tgcm_ubc

        use ppgrid,       only : pver
        use shr_kind_mod, only : r8 => shr_kind_r8
        use constituents, only : pcnst, cnst_fixed_ubc

        use abortutils,   only: endrun
        use cam_logfile,  only: iulog

        use tracer_data,  only : trfld,trfile,MAXTRCRS
        use cam_history,  only : addfld, phys_decomp

        implicit none

        private
        public  :: tgcm_ubc_inti, set_tgcm_ubc, tgcm_timestep_init

        save

        type(trfld), pointer :: fields(:)
        type(trfile)         :: file

        integer :: ub_nspecies
        character(len=16) :: ubc_name(MAXTRCRS)
        integer :: map(MAXTRCRS)

        logical :: ubc_from_tgcm(pcnst)  = .false.

      contains

      subroutine tgcm_ubc_inti( tgcm_ubc_file, tgcm_ubc_data_type, tgcm_ubc_cycle_yr, tgcm_ubc_fixed_ymd, tgcm_ubc_fixed_tod)
        !------------------------------------------------------------------
        !	... initialize upper boundary values
        !------------------------------------------------------------------
        use tracer_data, only : trcdata_init

        use ppgrid,          only : pcols, begchunk, endchunk
        use constituents,    only : cnst_get_ind, cnst_name
        use physics_buffer, only : physics_buffer_desc

        !------------------------------------------------------------------
        !	... dummy args
        !------------------------------------------------------------------
        character(len=*),   intent(in)     :: tgcm_ubc_file
        integer,            intent(in)     :: tgcm_ubc_cycle_yr
        integer,            intent(in)     :: tgcm_ubc_fixed_ymd
        integer,            intent(in)     :: tgcm_ubc_fixed_tod
        character(len=32),  intent(in)     :: tgcm_ubc_data_type


        ! local vars
        integer :: vid, i,ii 

        character(len=256), parameter :: filelist = ' '
        character(len=256), parameter :: datapath = ' '
        logical,            parameter :: rmv_file = .false.
        integer,            parameter :: nubc = 3
        character(len=4),   parameter :: specifier(nubc) = (/'CO  ','CO2 ','H2  '/)

        ii = 0

        do i = 1,nubc

           call cnst_get_ind( specifier(i), vid, abort=.false. )
           if( vid > 0 ) then
              if( cnst_fixed_ubc(vid) ) then
                 ii = ii+1
                 ubc_from_tgcm(vid) = .true.
                 map(ii) = vid
                 ubc_name(ii) = trim(specifier(i))//'_tgcm'
                 call addfld( ubc_name(ii), 'kg/kg', 1, 'I', 'upper boundary mmr', phys_decomp )
              end if
           end if
        enddo

        ub_nspecies = count( ubc_from_tgcm )

        if (ub_nspecies > 0) then
           file%top_bndry = .true.
           allocate(file%in_pbuf(size(specifier)))
           file%in_pbuf(:) = .false.
           call trcdata_init( specifier, tgcm_ubc_file, filelist, datapath, fields, file, &
                              rmv_file, tgcm_ubc_cycle_yr, tgcm_ubc_fixed_ymd, tgcm_ubc_fixed_tod, tgcm_ubc_data_type)
        endif

      end subroutine tgcm_ubc_inti

      subroutine tgcm_timestep_init(pbuf2d, state )

        use tracer_data,  only : advance_trcdata
        use physics_types,only : physics_state
        use ppgrid,       only : begchunk, endchunk
        use physics_buffer, only : physics_buffer_desc

        !--------------------------------------------------------------------
        !	... Advance ub values
        !--------------------------------------------------------------------
        implicit none

        ! args
        type(physics_state), intent(in):: state(begchunk:endchunk)
        type(physics_buffer_desc), pointer :: pbuf2d(:,:)

        if (ub_nspecies > 0) then
           call advance_trcdata( fields, file, state, pbuf2d )
        endif

      end subroutine tgcm_timestep_init

      subroutine set_tgcm_ubc( lchunk, ncol, mmr, mw_dry )
        !--------------------------------------------------------------------
        !	... Set the upper boundary values h2o, h2, and h
        !--------------------------------------------------------------------
        
        use ppgrid,       only : pcols
        use constituents, only : cnst_get_ind, cnst_mw

        use cam_history,  only : outfld

        implicit none

        !--------------------------------------------------------------------
        !	... dummy args
        !--------------------------------------------------------------------
        integer,  intent(in)    :: lchunk            ! chunk id
        integer,  intent(in)    :: ncol              ! columns in chunk
        real(r8), intent(in)    :: mw_dry(pcols)     ! mean mass at top model level
        real(r8), intent(inout) :: mmr(pcols,pcnst)

        !--------------------------------------------------------------------
        !	... local variables
        !--------------------------------------------------------------------
        real(r8), parameter ::  h2o_ubc_vmr = 2.e-8_r8            ! fixed ub h2o concentration (kg/kg)
        real(r8), parameter ::  ch4_ubc_vmr = 2.e-10_r8           ! fixed ub ch4 concentration (kg/kg)

        integer  :: m,n,i

        if (ub_nspecies > 0) then
           do m = 1,ub_nspecies
!---------------------------------------------------------------
!	... tgcm upper bndy values
!---------------------------------------------------------------

              n = map(m)
              mmr(:ncol,n) = fields(m)%data(:ncol,1,lchunk)
              call outfld( ubc_name(m), mmr(:ncol,n), ncol, lchunk )
           enddo
        endif

        !--------------------------------------------------------
        !	... special section to set h2o and ch4 ub concentrations
        !--------------------------------------------------------
        mmr(:ncol,1) = cnst_mw(1)*h2o_ubc_vmr/mw_dry(:ncol)
        call cnst_get_ind( 'CH4', m, abort=.false. )
        if( m > 0 ) then
           mmr(:ncol,m) = cnst_mw(m)*ch4_ubc_vmr/mw_dry(:ncol)
        end if

#ifdef TGCM_DIAGS
        call cnst_get_ind( 'H2', m, abort=.false. )
        if( m > 0 ) then
           write(iulog,*) 'set_ub_vals: diagnostics for chunk = ',lchunk
           write(iulog,*) 'last,next,dels = ',last,next,dels
           write(iulog,*) 'h2 mmr at level ',k
           write(iulog,'(1x,1p,10g12.5)') mmr(:ncol,m))
        end if
#endif

      end subroutine set_tgcm_ubc

      end module mo_tgcm_ubc
