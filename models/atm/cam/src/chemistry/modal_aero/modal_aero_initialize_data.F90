module modal_aero_initialize_data
  use cam_logfile,           only : iulog
  use abortutils,            only: endrun
  use spmd_utils,            only: masterproc, iam
  use ppgrid,                only: pcols, pver, begchunk, endchunk
  use modal_aero_data
  use time_manager,          only: is_first_step

  implicit none
  private

  public :: modal_aero_register
  public :: modal_aero_initialize
  public :: modal_aero_initialize_q

contains

  subroutine modal_aero_register
    use constituents,only: pcnst, cnst_name
    use physics_buffer, only : pbuf_add_field, dtype_r8

    character(len=8)  :: &
         xname_massptr(maxd_aspectype,ntot_amode), &
         xname_massptrcw(maxd_aspectype,ntot_amode)
    character(len=10) :: xname_spectype(maxd_aspectype,ntot_amode)


    !   input species to hold interstitial & activated number
#if ( defined MODAL_AERO_7MODE )
    character(len=*), parameter :: xname_numptr(ntot_amode)   = (/ 'num_a1  ', 'num_a2  ', 'num_a3  ', &
         'num_a4  ', 'num_a5  ', 'num_a6  ', 'num_a7  ' /)
    character(len=*), parameter ::     xname_numptrcw(ntot_amode) = (/ 'num_c1  ', 'num_c2  ', 'num_c3  ', &
         'num_c4  ', 'num_c5  ', 'num_c6  ', 'num_c7  ' /)
#elif ( defined MODAL_AERO_3MODE )
    character(len=*), parameter ::     xname_numptr(ntot_amode)   = (/ 'num_a1  ', 'num_a2  ', &
         'num_a3  ' /)
    character(len=*), parameter ::     xname_numptrcw(ntot_amode) = (/ 'num_c1  ', 'num_c2  ', &
         'num_c3  ' /)
#endif



    integer :: m, l, iptr
    character(len=3) :: trnum       ! used to hold mode number (as characters)

       !   input species to hold aerosol water and "kohler-c"
       !     xname_waterptr(:ntot_amode)   = (/ 'wat_a1  ', 'wat_a2  ', 'wat_a3  ', &
       !                                        'wat_a4  ', 'wat_a5  ', 'wat_a6  ', 'wat_a7  ' /)
       !   input chemical species for the mode
       ! mode 1 (accumulation) species
       ! JPE 02022011: These could also be parameters but a bug in the pathscale compiler prevents
       !               parameter initialization of 2D variables
#if ( defined MODAL_AERO_7MODE )
       xname_massptr(:nspec_amode(1),1)   = (/ 'so4_a1  ', 'nh4_a1  ', &
            'pom_a1  ', 'soa_a1  ', 'bc_a1   ', 'ncl_a1  ' /)
       xname_massptrcw(:nspec_amode(1),1) = (/ 'so4_c1  ', 'nh4_c1  ', &
            'pom_c1  ', 'soa_c1  ', 'bc_c1   ', 'ncl_c1  ' /)
       xname_spectype(:nspec_amode(1),1)  = (/ 'sulfate   ', 'ammonium  ', &
            'p-organic ', 's-organic ', 'black-c   ', 'seasalt   ' /)
#elif ( defined MODAL_AERO_3MODE )
       xname_massptr(:nspec_amode(1),1)   = (/ 'so4_a1  ', &
            'pom_a1  ', 'soa_a1  ', 'bc_a1   ', &
            'dst_a1  ', 'ncl_a1  ' /)
       xname_massptrcw(:nspec_amode(1),1) = (/ 'so4_c1  ', &
            'pom_c1  ', 'soa_c1  ', 'bc_c1   ', &
            'dst_c1  ', 'ncl_c1  ' /)
       xname_spectype(:nspec_amode(1),1)  = (/ 'sulfate   ', &
            'p-organic ', 's-organic ', 'black-c   ', &
            'dust      ', 'seasalt   ' /)
#endif

       ! mode 2 (aitken) species
#if ( defined MODAL_AERO_7MODE )
       xname_massptr(:nspec_amode(2),2)   = (/ 'so4_a2  ', 'nh4_a2  ', &
            'soa_a2  ', 'ncl_a2  ' /)
       xname_massptrcw(:nspec_amode(2),2) = (/ 'so4_c2  ', 'nh4_c2  ', &
            'soa_c2  ', 'ncl_c2  ' /)
       xname_spectype(:nspec_amode(2),2)  = (/ 'sulfate   ', 'ammonium  ', &
            's-organic ', 'seasalt   ' /)
#elif ( defined MODAL_AERO_3MODE )
       xname_massptr(:nspec_amode(2),2)   = (/ 'so4_a2  ', &
            'soa_a2  ', 'ncl_a2  ' /)
       xname_massptrcw(:nspec_amode(2),2) = (/ 'so4_c2  ', &
            'soa_c2  ', 'ncl_c2  ' /)
       xname_spectype(:nspec_amode(2),2)  = (/ 'sulfate   ', &
            's-organic ', 'seasalt   ' /)
#endif

#if ( defined MODAL_AERO_7MODE )
       ! mode 3 (primary carbon) species
       xname_massptr(:nspec_amode(3),3)   = (/ 'pom_a3  ', 'bc_a3   ' /)
       xname_massptrcw(:nspec_amode(3),3) = (/ 'pom_c3  ', 'bc_c3   ' /)
       xname_spectype(:nspec_amode(3),3)  = (/ 'p-organic ', 'black-c   ' /)
#elif ( defined MODAL_AERO_3MODE )
       ! mode 3 (coarse dust & seasalt) species
       xname_massptr(:nspec_amode(3),3)   = (/ 'dst_a3  ', 'ncl_a3  ', 'so4_a3  ' /)
       xname_massptrcw(:nspec_amode(3),3) = (/ 'dst_c3  ', 'ncl_c3  ', 'so4_c3  ' /)
       xname_spectype(:nspec_amode(3),3)  = (/ 'dust      ', 'seasalt   ', 'sulfate   ' /)
#endif


#if ( defined MODAL_AERO_7MODE )
       ! mode 4 (fine seasalt) species
       xname_massptr(:nspec_amode(4),4)   = (/ 'ncl_a4  ', 'so4_a4  ', 'nh4_a4  ' /)
       xname_massptrcw(:nspec_amode(4),4) = (/ 'ncl_c4  ', 'so4_c4  ', 'nh4_c4  ' /)
       xname_spectype(:nspec_amode(4),4)  = (/ 'seasalt   ', 'sulfate   ', 'ammonium  ' /)

       ! mode 5 (fine dust) species
       xname_massptr(:nspec_amode(5),5)   = (/ 'dst_a5  ', 'so4_a5  ', 'nh4_a5  ' /)
       xname_massptrcw(:nspec_amode(5),5) = (/ 'dst_c5  ', 'so4_c5  ', 'nh4_c5  ' /)
       xname_spectype(:nspec_amode(5),5)  = (/ 'dust      ', 'sulfate   ', 'ammonium  ' /)

       ! mode 6 (coarse seasalt) species
       xname_massptr(:nspec_amode(6),6)   = (/ 'ncl_a6  ', 'so4_a6  ', 'nh4_a6  ' /)
       xname_massptrcw(:nspec_amode(6),6) = (/ 'ncl_c6  ', 'so4_c6  ', 'nh4_c6  ' /)
       xname_spectype(:nspec_amode(6),6)  = (/ 'seasalt   ', 'sulfate   ', 'ammonium  ' /)

       ! mode 7 (coarse dust) species
       xname_massptr(:nspec_amode(7),7)   = (/ 'dst_a7  ', 'so4_a7  ', 'nh4_a7  ' /)
       xname_massptrcw(:nspec_amode(7),7) = (/ 'dst_c7  ', 'so4_c7  ', 'nh4_c7  ' /)
       xname_spectype(:nspec_amode(7),7)  = (/ 'dust      ', 'sulfate   ', 'ammonium  ' /)
#endif

    do m = 1, ntot_amode

       if (masterproc) then
          write(iulog,9231) m, modename_amode(m)
          write(iulog,9232)                                          &
               'nspec                       ',                         &
               nspec_amode(m)
          write(iulog,9232)                                          &
               'mprognum, mdiagnum, mprogsfc',                         &
               mprognum_amode(m), mdiagnum_amode(m), mprogsfc_amode(m)
          write(iulog,9232)                                          &
               'mcalcwater                  ',                         &
               mcalcwater_amode(m)
       endif

       !    define species to hold interstitial & activated number
       call search_list_of_names(                                      &
            xname_numptr(m), numptr_amode(m), cnst_name, pcnst )
       if (numptr_amode(m) .le. 0) then
          write(iulog,9061) 'xname_numptr', xname_numptr(m), m
          call endrun()
       end if
       if (numptr_amode(m) .gt. pcnst) then
          write(iulog,9061) 'numptr_amode', numptr_amode(m), m
          write(iulog,9061) 'xname_numptr', xname_numptr(m), m
          call endrun()
       end if

       species_class(numptr_amode(m)) = spec_class_aerosol


       numptrcw_amode(m) = numptr_amode(m)  !use the same index for Q and QQCW arrays
       if (numptrcw_amode(m) .le. 0) then
          write(iulog,9061) 'xname_numptrcw', xname_numptrcw(m), m
          call endrun()
       end if
       if (numptrcw_amode(m) .gt. pcnst) then
          write(iulog,9061) 'numptrcw_amode', numptrcw_amode(m), m
          write(iulog,9061) 'xname_numptrcw', xname_numptrcw(m), m
          call endrun()
       end if
       species_class(numptrcw_amode(m)) = spec_class_aerosol

       call pbuf_add_field(xname_numptrcw(m),'global',dtype_r8,(/pcols,pver/),iptr)
       call qqcw_set_ptr(numptrcw_amode(m),iptr)

       !   output mode information
       if ( masterproc ) then
          write(iulog,9233) 'numptr         ',                           &
               numptr_amode(m), xname_numptr(m)
          write(iulog,9233) 'numptrcw       ',                           &
               numptrcw_amode(m), xname_numptrcw(m)
       end if


       !   define the chemical species for the mode
       do l = 1, nspec_amode(m)

          call search_list_of_names(                                  &
               xname_spectype(l,m), lspectype_amode(l,m),              &
               specname_amode, ntot_aspectype )
          if (lspectype_amode(l,m) .le. 0) then
             write(iulog,9062) 'xname_spectype', xname_spectype(l,m), l, m
             call endrun()
          end if

          call search_list_of_names(                                  &
               xname_massptr(l,m), lmassptr_amode(l,m), cnst_name, pcnst )
          if (lmassptr_amode(l,m) .le. 0) then
             write(iulog,9062) 'xname_massptr', xname_massptr(l,m), l, m
             call endrun()
          end if
          species_class(lmassptr_amode(l,m)) = spec_class_aerosol

          lmassptrcw_amode(l,m) = lmassptr_amode(l,m)  !use the same index for Q and QQCW arrays
          if (lmassptrcw_amode(l,m) .le. 0) then
             write(iulog,9062) 'xname_massptrcw', xname_massptrcw(l,m), l, m
             call endrun()
          end if
          call pbuf_add_field(xname_massptrcw(l,m),'global',dtype_r8,(/pcols,pver/),iptr)
          call qqcw_set_ptr(lmassptrcw_amode(l,m), iptr)
          species_class(lmassptrcw_amode(l,m)) = spec_class_aerosol

          if ( masterproc ) then
             write(iulog,9236) 'spec, spectype ', l,                    &
                  lspectype_amode(l,m), xname_spectype(l,m)
             write(iulog,9236) 'spec, massptr  ', l,                    &
                  lmassptr_amode(l,m), xname_massptr(l,m)
             write(iulog,9236) 'spec, massptrcw', l,                    &
                  lmassptrcw_amode(l,m), xname_massptrcw(l,m)
          end if

       enddo

       if ( masterproc ) write(iulog,*)


       !   set names for aodvis and ssavis
       write(unit=trnum,fmt='(i3)') m+100
       aodvisname(m) = 'AODVIS'//trnum(2:3)
       aodvislongname(m) = 'Aerosol optical depth for mode '//trnum(2:3)
       ssavisname(m) = 'SSAVIS'//trnum(2:3)
       ssavislongname(m) = 'Single-scatter albedo for mode '//trnum(2:3)
       fnactname(m) = 'FNACT'//trnum(2:3)
       fnactlongname(m) = 'Number faction activated for mode '//trnum(2:3)
       fmactname(m) = 'FMACT'//trnum(2:3)
       fmactlongname(m) = 'Fraction mass activated for mode'//trnum(2:3)
    end do

       if (masterproc) write(iulog,9230)
9230   format( // '*** init_aer_modes mode definitions' )
9231   format( 'mode = ', i4, ' = "', a, '"' )
9232   format( 4x, a, 4(1x, i5 ) )
9233   format( 4x, a15, 4x, i7, '="', a, '"' )
9236   format( 4x, a15, i4, i7, '="', a, '"' )
9061   format( '*** subr init_aer_modes - bad ', a /                   &
            5x, 'name, m =  ', a, 5x, i5 )
9062   format( '*** subr init_aer_modesaeromodeinit - bad ', a /                       &
            5x, 'name, l, m =  ', a, 5x, 2i5 )

  end subroutine modal_aero_register


  !==============================================================
  subroutine modal_aero_initialize(pbuf2d)

       use constituents,          only: pcnst
       use physconst,             only: rhoh2o, mwh2o
       use modal_aero_calcsize,   only: modal_aero_calcsize_init
       use modal_aero_coag,       only: modal_aero_coag_init
       use modal_aero_deposition, only: modal_aero_deposition_init
       use modal_aero_gasaerexch, only: modal_aero_gasaerexch_init
       use modal_aero_newnuc,     only: modal_aero_newnuc_init
       use modal_aero_rename,     only: modal_aero_rename_init
       use rad_constituents,      only: rad_cnst_get_info, rad_cnst_get_aer_props, &
                                        rad_cnst_get_mode_props
       use aerodep_flx,           only: aerodep_flx_prescribed
       use physics_buffer,        only: physics_buffer_desc, pbuf_get_chunk

       type(physics_buffer_desc), pointer :: pbuf2d(:,:)

       !--------------------------------------------------------------
       ! ... local variables
       !--------------------------------------------------------------
       integer :: l, m, i, lchnk
       integer :: m_idx, s_idx

       character(len=3) :: trnum       ! used to hold mode number (as characters)
       integer :: iaerosol, ibulk
       integer  :: numaerosols     ! number of bulk aerosols in climate list
       character(len=20) :: bulkname
       real(r8) :: pi
       complex(r8), pointer  :: refindex_aer_sw(:), &
            refindex_aer_lw(:)
       real(r8), pointer :: qqcw(:,:)
       real(r8), parameter :: huge_r8 = huge(1._r8)
       character(len=*), parameter :: routine='modal_aero_initialize'
       !-----------------------------------------------------------------------

       pi = 4._r8*atan(1._r8)    

       ! safety check on modal_aero, and modal_aero_3mode, modal_aero_7mode
#if ( defined MODAL_AERO_3MODE ) && ( defined MODAL_AERO_7MODE )
       call endrun( 'Error - when modal_aero defined, just 1 of modal_aero_3/7mode must be defined'
#elif ( ! ( defined MODAL_AERO_3MODE ) ) && ( ! ( defined MODAL_AERO_7MODE ) )
       call endrun( 'Error - when modal_aero defined, at least 1 of modal_aero_3/7mode must be defined'
#endif

       ! Mode specific properties.
       do m = 1, ntot_amode
          call rad_cnst_get_mode_props(0, m, &
             sigmag=sigmag_amode(m), dgnum=dgnum_amode(m), dgnumlo=dgnumlo_amode(m), &
             dgnumhi=dgnumhi_amode(m), rhcrystal=rhcrystal_amode(m), rhdeliques=rhdeliques_amode(m))

          !   compute frequently used parameters: ln(sigmag),
          !   volume-to-number and volume-to-surface conversions, ...
          alnsg_amode(m) = log( sigmag_amode(m) )

          voltonumb_amode(m) = 1._r8 / ( (pi/6._r8)*                            &
             (dgnum_amode(m)**3._r8)*exp(4.5_r8*alnsg_amode(m)**2._r8) )
          voltonumblo_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
             (dgnumlo_amode(m)**3._r8)*exp(4.5_r8*alnsg_amode(m)**2._r8) )
          voltonumbhi_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
             (dgnumhi_amode(m)**3._r8)*exp(4.5_r8*alnsg_amode(m)**2._r8) )

          alnv2n_amode(m)   = log( voltonumb_amode(m) )
          alnv2nlo_amode(m) = log( voltonumblo_amode(m) )
          alnv2nhi_amode(m) = log( voltonumbhi_amode(m) )
       end do


       ! Properties of mode specie types.

       !     values from Koepke, Hess, Schult and Shettle, Global Aerosol Data Set 
       !     Report #243, Max-Planck Institute for Meteorology, 1997a
       !     See also Hess, Koepke and Schult, Optical Properties of Aerosols and Clouds (OPAC)
       !     BAMS, 1998.

       !      specrefndxsw(:ntot_aspectype)     = (/ (1.53,  0.01),   (1.53,  0.01),  (1.53,  0.01), &
       !                                           (1.55,  0.01),   (1.55,  0.01),  (1.90, 0.60), &
       !                                           (1.50, 1.0e-8), (1.50, 0.005) /)
       !      specrefndxlw(:ntot_aspectype)   = (/ (2.0, 0.5),   (2.0, 0.5), (2.0, 0.5), &
       !                                           (1.7, 0.5),   (1.7, 0.5), (2.22, 0.73), &
       !                                           (1.50, 0.02), (2.6, 0.6) /)
       !     get refractive indices from phys_prop files

       ! The following use of the rad_constituent interfaces makes the assumption that the
       ! prognostic modes are used in the mode climate (index 0) list.
       do l = 1, ntot_aspectype

          ! specname_amode is the species type.  This info call will return the mode and species
          ! indices of the first occurance of the species type.
          call rad_cnst_get_info(0, specname_amode(l), mode_idx=m_idx, spec_idx=s_idx)

          if (m_idx > 0 .and. s_idx > 0) then

             call rad_cnst_get_aer_props(0, m_idx, s_idx, &
                refindex_aer_sw=refindex_aer_sw, &
                refindex_aer_lw=refindex_aer_lw, &
                density_aer=specdens_amode(l), &
                hygro_aer=spechygro(l))

             specrefndxsw(:nswbands,l) = refindex_aer_sw(:nswbands)
             specrefndxlw(:nlwbands,l) = refindex_aer_lw(:nlwbands)

          else
             if (masterproc) then
                write(iulog,*) routine//': INFO: props not found for species type: ',trim(specname_amode(l))
             end if
             specdens_amode(l)         = huge_r8
             spechygro(l)              = huge_r8
             specrefndxsw(:nswbands,l) = (huge_r8, huge_r8)
             specrefndxlw(:nlwbands,l) = (huge_r8, huge_r8)
          endif

       end do


       if (masterproc) write(iulog,9210)
       do l = 1, ntot_aspectype
          !            spechygro(l) = specnu(l)*specphi(l)*specsolfrac(l)*mwh2o*specdens_amode(l) / &
          !	               (rhoh2o*specmw_amode(l))
          if (masterproc) then
             write(iulog,9211) l
             write(iulog,9212) 'name            ', specname_amode(l)
             write(iulog,9213) 'density, MW     ',                  &
                  specdens_amode(l), specmw_amode(l)
             write(iulog,9213) 'hygro', spechygro(l)
             do i=1,nswbands
                write(iulog,9213) 'ref index sw    ', (specrefndxsw(i,l))
             end do
             do i=1,nlwbands
                write(iulog,9213) 'ref index ir    ', (specrefndxlw(i,l))
             end do
          end if
       end do

9210   format( // '*** init_aer_modes aerosol species-types' )
9211   format( 'spectype =', i4)
9212   format( 4x, a, 3x, '"', a, '"' )
9213   format( 4x, a, 5(1pe14.5) )





       do i = 1, pcnst
          species_class(i) = spec_class_undefined
       end do



       !   set cnst_name_cw
       call initaermodes_set_cnstnamecw()


       !
       !   set the lptr_so4_a_amode(m), lptr_so4_cw_amode(m), ...
       !
       call initaermodes_setspecptrs

       !
       !   set threshold for reporting negatives from subr qneg3
       !   for aerosol number species set this to
       !      1e3 #/kg ~= 1e-3 #/cm3 for accum, aitken, pcarbon, ufine modes
       !      3e1 #/kg ~= 3e-5 #/cm3 for fineseas and finedust modes 
       !      1e0 #/kg ~= 1e-6 #/cm3 for other modes which are coarse
       !   for other species, set this to zero so that it will be ignored
       !      by qneg3
       !
       if ( masterproc ) write(iulog,'(/a)') &
            'mode, modename_amode, qneg3_worst_thresh_amode'
       qneg3_worst_thresh_amode(:) = 0.0_r8
       do m = 1, ntot_amode
          l = numptr_amode(m)
          if ((l <= 0) .or. (l > pcnst)) cycle

          if      (m == modeptr_accum) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8
          else if (m == modeptr_aitken) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8
          else if (m == modeptr_pcarbon) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8
          else if (m == modeptr_ufine) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8

          else if (m == modeptr_fineseas) then
             qneg3_worst_thresh_amode(l) = 3.0e1_r8
          else if (m == modeptr_finedust) then
             qneg3_worst_thresh_amode(l) = 3.0e1_r8

          else
             qneg3_worst_thresh_amode(l) = 1.0e0_r8
          end if

          if ( masterproc ) write(iulog,'(i3,2x,a,1p,e12.3)') &
               m, modename_amode(m), qneg3_worst_thresh_amode(l)
       end do


       !
       !   call other initialization routines
       !
       call modal_aero_rename_init
       !   calcsize call must follow rename call
       call modal_aero_calcsize_init( pbuf2d )
       call modal_aero_gasaerexch_init
       !   coag call must follow gasaerexch call
       call modal_aero_coag_init
       call modal_aero_newnuc_init

       ! call modal_aero_deposition_init only if the user has not specified 
       ! prescribed aerosol deposition fluxes
       if (.not.aerodep_flx_prescribed()) then
          call modal_aero_deposition_init
       endif

       if (is_first_step()) then
          ! initialize cloud bourne constituents in physics buffer

          do i = 1, pcnst
             do lchnk = begchunk, endchunk
                qqcw => qqcw_get_field(pbuf_get_chunk(pbuf2d,lchnk), i, lchnk, .true.)
                if (associated(qqcw)) then
                   qqcw = 1.e-38_r8
                end if
             end do
          end do
       end if

       return
     end subroutine modal_aero_initialize


     !==============================================================
     subroutine search_list_of_names(                                &
          name_to_find, name_id, list_of_names, list_length )
       !
       !   searches for a name in a list of names
       !
       !   name_to_find - the name to be found in the list  [input]
       !   name_id - the position of "name_to_find" in the "list_of_names".
       !       If the name is not found in the list, then name_id=0.  [output]
       !   list_of_names - the list of names to be searched  [input]
       !   list_length - the number of names in the list  [input]
       !
       character(len=*), intent(in):: name_to_find, list_of_names(:)
       integer, intent(in) :: list_length
       integer, intent(out) :: name_id
       
       integer :: i
       name_id = -999888777
       if (name_to_find .ne. ' ') then
          do i = 1, list_length
             if (name_to_find .eq. list_of_names(i)) then
                name_id = i
                exit
             end if
          end do
       end if
     end subroutine search_list_of_names


     !==============================================================
     subroutine initaermodes_setspecptrs
       !
       !   sets the lptr_so4_a_amode(m), lptr_so4_cw_amode(m), ...
       !       and writes them to iulog
       !   ALSO sets the mode-pointers:  modeptr_accum, modeptr_aitken, ...
       !       and writes them to iulog
       !   ALSO sets values of specdens_XX_amode and specmw_XX_amode
       !       (XX = so4, om, bc, dust, seasalt)
       !
       implicit none

       !   local variables
       integer l, l2, m
       character*8 dumname
       integer, parameter :: init_val=-999888777

       !   all processes set the pointers

       modeptr_accum = init_val
       modeptr_aitken = init_val
       modeptr_ufine = init_val
       modeptr_coarse = init_val
       modeptr_pcarbon = init_val
       modeptr_fineseas = init_val
       modeptr_finedust = init_val
       modeptr_coarseas = init_val
       modeptr_coardust = init_val
       do m = 1, ntot_amode
          if (modename_amode(m) .eq. 'accum') then
             modeptr_accum = m
          else if (modename_amode(m) .eq. 'aitken') then
             modeptr_aitken = m
          else if (modename_amode(m) .eq. 'ufine') then
             modeptr_ufine = m
          else if (modename_amode(m) .eq. 'coarse') then
             modeptr_coarse = m
          else if (modename_amode(m) .eq. 'primary_carbon') then
             modeptr_pcarbon = m
          else if (modename_amode(m) .eq. 'fine_seasalt') then
             modeptr_fineseas = m
          else if (modename_amode(m) .eq. 'fine_dust') then
             modeptr_finedust = m
          else if (modename_amode(m) .eq. 'coarse_seasalt') then
             modeptr_coarseas = m
          else if (modename_amode(m) .eq. 'coarse_dust') then
             modeptr_coardust = m
          end if
       end do

       do m = 1, ntot_amode
          lptr_so4_a_amode(m)   = init_val
          lptr_so4_cw_amode(m)  = init_val
          lptr_msa_a_amode(m)   = init_val
          lptr_msa_cw_amode(m)  = init_val
          lptr_nh4_a_amode(m)   = init_val
          lptr_nh4_cw_amode(m)  = init_val
          lptr_no3_a_amode(m)   = init_val
          lptr_no3_cw_amode(m)  = init_val
          lptr_pom_a_amode(m)   = init_val
          lptr_pom_cw_amode(m)  = init_val
          lptr_soa_a_amode(m)   = init_val
          lptr_soa_cw_amode(m)  = init_val
          lptr_bc_a_amode(m)    = init_val
          lptr_bc_cw_amode(m)   = init_val
          lptr_nacl_a_amode(m)  = init_val
          lptr_nacl_cw_amode(m) = init_val
          lptr_dust_a_amode(m)  = init_val
          lptr_dust_cw_amode(m) = init_val
          do l = 1, nspec_amode(m)
             l2 = lspectype_amode(l,m)
             if ( (specname_amode(l2) .eq. 'sulfate') .and.  &
                  (lptr_so4_a_amode(m) .le. 0) ) then
                lptr_so4_a_amode(m)  = lmassptr_amode(l,m)
                lptr_so4_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'msa') .and.      &
                  (lptr_msa_a_amode(m) .le. 0) ) then
                lptr_msa_a_amode(m)  = lmassptr_amode(l,m)
                lptr_msa_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'ammonium') .and.  &
                  (lptr_nh4_a_amode(m) .le. 0) ) then
                lptr_nh4_a_amode(m)  = lmassptr_amode(l,m)
                lptr_nh4_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'nitrate') .and.  &
                  (lptr_no3_a_amode(m) .le. 0) ) then
                lptr_no3_a_amode(m)  = lmassptr_amode(l,m)
                lptr_no3_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'p-organic') .and.   &
                  (lptr_pom_a_amode(m) .le. 0) ) then
                lptr_pom_a_amode(m)  = lmassptr_amode(l,m)
                lptr_pom_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 's-organic') .and.   &
                  (lptr_soa_a_amode(m) .le. 0) ) then
                lptr_soa_a_amode(m)  = lmassptr_amode(l,m)
                lptr_soa_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'black-c') .and.  &
                  (lptr_bc_a_amode(m) .le. 0) ) then
                lptr_bc_a_amode(m)  = lmassptr_amode(l,m)
                lptr_bc_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'seasalt') .and.  &
                  (lptr_nacl_a_amode(m) .le. 0) ) then
                lptr_nacl_a_amode(m)  = lmassptr_amode(l,m)
                lptr_nacl_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
             if ( (specname_amode(l2) .eq. 'dust') .and.     &
                  (lptr_dust_a_amode(m) .le. 0) ) then
                lptr_dust_a_amode(m)  = lmassptr_amode(l,m)
                lptr_dust_cw_amode(m) = lmassptrcw_amode(l,m)
             end if
          end do
       end do

       !   all processes set values of specdens_XX_amode and specmw_XX_amode
       specdens_so4_amode = 2.0_r8
       specdens_nh4_amode = 2.0_r8
       specdens_no3_amode = 2.0_r8
       specdens_pom_amode = 2.0_r8
       specdens_soa_amode = 2.0_r8
       specdens_bc_amode = 2.0_r8
       specdens_dust_amode = 2.0_r8
       specdens_seasalt_amode = 2.0_r8
       specmw_so4_amode = 1.0_r8
       specmw_nh4_amode = 1.0_r8
       specmw_no3_amode = 1.0_r8
       specmw_pom_amode = 1.0_r8
       specmw_soa_amode = 1.0_r8
       specmw_bc_amode = 1.0_r8
       specmw_dust_amode = 1.0_r8
       specmw_seasalt_amode = 1.0_r8
       do m = 1, ntot_aspectype
          if      (specname_amode(m).eq.'sulfate   ') then
             specdens_so4_amode = specdens_amode(m)
             specmw_so4_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'ammonium  ') then
             specdens_nh4_amode = specdens_amode(m)
             specmw_nh4_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'nitrate   ') then
             specdens_no3_amode = specdens_amode(m)
             specmw_no3_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'p-organic ') then
             specdens_pom_amode = specdens_amode(m)
             specmw_pom_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'s-organic ') then
             specdens_soa_amode = specdens_amode(m)
             specmw_soa_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'black-c   ') then
             specdens_bc_amode = specdens_amode(m)
             specmw_bc_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'dust      ') then
             specdens_dust_amode = specdens_amode(m)
             specmw_dust_amode = specmw_amode(m)
          else if (specname_amode(m).eq.'seasalt   ') then
             specdens_seasalt_amode = specdens_amode(m)
             specmw_seasalt_amode = specmw_amode(m)
          end if
       enddo

       !   masterproc writes out the pointers
       if ( .not. ( masterproc ) ) return

       write(iulog,9230)
       write(iulog,*) 'modeptr_accum    =', modeptr_accum
       write(iulog,*) 'modeptr_aitken   =', modeptr_aitken
       write(iulog,*) 'modeptr_ufine    =', modeptr_ufine
       write(iulog,*) 'modeptr_coarse   =', modeptr_coarse
       write(iulog,*) 'modeptr_pcarbon  =', modeptr_pcarbon
       write(iulog,*) 'modeptr_fineseas =', modeptr_fineseas
       write(iulog,*) 'modeptr_finedust =', modeptr_finedust
       write(iulog,*) 'modeptr_coarseas =', modeptr_coarseas
       write(iulog,*) 'modeptr_coardust =', modeptr_coardust

       dumname = 'none'
       write(iulog,9240)
       write(iulog,9000) 'sulfate    '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_so4_a_amode(m), lptr_so4_cw_amode(m),  'so4' )
       end do

       write(iulog,9000) 'msa        '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_msa_a_amode(m), lptr_msa_cw_amode(m),  'msa' )
       end do

       write(iulog,9000) 'ammonium   '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_nh4_a_amode(m), lptr_nh4_cw_amode(m),  'nh4' )
       end do

       write(iulog,9000) 'nitrate    '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_no3_a_amode(m), lptr_no3_cw_amode(m),  'no3' )
       end do

       write(iulog,9000) 'p-organic  '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_pom_a_amode(m), lptr_pom_cw_amode(m),  'pom' )
       end do

       write(iulog,9000) 's-organic  '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_soa_a_amode(m), lptr_soa_cw_amode(m),  'soa' )
       end do

       write(iulog,9000) 'black-c    '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_bc_a_amode(m), lptr_bc_cw_amode(m),  'bc' )
       end do

       write(iulog,9000) 'seasalt   '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_nacl_a_amode(m), lptr_nacl_cw_amode(m),  'nacl' )
       end do

       write(iulog,9000) 'dust       '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_dust_a_amode(m), lptr_dust_cw_amode(m),  'dust' )
       end do

9000   format( a )
9230   format(                                                         &
            / 'mode-pointer output from subr initaermodes_setspecptrs' )
9240   format(                                                         &
            / 'species-pointer output from subr initaermodes_setspecptrs' / &
            'mode', 12x, 'id  name_a  ', 12x, 'id  name_cw' )

       return
     end subroutine initaermodes_setspecptrs


     !==============================================================
     subroutine initaermodes_setspecptrs_write2(                     &
          m, laptr, lcptr, txtdum )
       !
       !   does some output for initaermodes_setspecptrs

       use constituents, only: pcnst, cnst_name

       implicit none

       !   subr arguments
       integer m, laptr, lcptr
       character*(*) txtdum

       !   local variables
       character*8 dumnamea, dumnamec

       dumnamea = 'none'
       dumnamec = 'none'
       if (laptr .gt. 0) dumnamea = cnst_name(laptr)
       if (lcptr .gt. 0) dumnamec = cnst_name(lcptr)
       write(iulog,9241) m, laptr, dumnamea, lcptr, dumnamec, txtdum

9241   format( i4, 2( 2x, i12, 2x, a ),                                &
            4x, 'lptr_', a, '_a/cw_amode' )

       return
     end subroutine initaermodes_setspecptrs_write2


     !==============================================================
     subroutine initaermodes_set_cnstnamecw
       !
       !   sets the cnst_name_cw
       !
       use constituents, only: pcnst, cnst_name
       implicit none

       !   subr arguments (none)

       !   local variables
       integer j, l, la, lc, ll, m

       !   set cnst_name_cw
       cnst_name_cw = ' '
       do m = 1, ntot_amode
          do ll = 0, nspec_amode(m)
             if (ll == 0) then
                la = numptr_amode(m)
                lc = numptrcw_amode(m)
             else
                la = lmassptr_amode(ll,m)
                lc = lmassptrcw_amode(ll,m)
             end if
             if ((la < 1) .or. (la > pcnst) .or.   &
                  (lc < 1) .or. (lc > pcnst)) then
                write(*,'(/2a/a,5(1x,i10))')   &
                     '*** initaermodes_set_cnstnamecw error',   &
                     ' -- bad la or lc',   &
                     '    m, ll, la, lc, pcnst =', m, ll, la, lc, pcnst
                call endrun( '*** initaermodes_set_cnstnamecw error' )
             end if
             do j = 2, len( cnst_name(la) ) - 1
                if (cnst_name(la)(j:j+1) == '_a') then
                   cnst_name_cw(lc) = cnst_name(la)
                   cnst_name_cw(lc)(j:j+1) = '_c'
                   exit
                else if (cnst_name(la)(j:j+1) == '_A') then
                   cnst_name_cw(lc) = cnst_name(la)
                   cnst_name_cw(lc)(j:j+1) = '_C'
                   exit
                end if
             end do
             if (cnst_name_cw(lc) == ' ') then
                write(*,'(/2a/a,3(1x,i10),2x,a)')   &
                     '*** initaermodes_set_cnstnamecw error',   &
                     ' -- bad cnst_name(la)',   &
                     '    m, ll, la, cnst_name(la) =',   &
                     m, ll, la, cnst_name(la)
                call endrun( '*** initaermodes_set_cnstnamecw error' )
             end if
          end do   ! ll = 0, nspec_amode(m)
       end do   ! m = 1, ntot_amode

       if ( masterproc ) then
          write(*,'(/a)') 'l, cnst_name(l), cnst_name_cw(l)'
          do l = 1, pcnst
             write(*,'(i4,2(2x,a))') l, cnst_name(l), cnst_name_cw(l)
          end do
       end if

       return
     end subroutine initaermodes_set_cnstnamecw


     !==============================================================
     subroutine modal_aero_initialize_q( name, q )
       !
       ! this routine is for initial testing of the modal aerosol cam3
       !
       ! it initializes several gas and aerosol species to 
       !    "low background" values, so that very short (e.g., 1 day)
       !    test runs are working with non-zero values
       !
       use constituents, only: pcnst, cnst_name
       use pmgrid,      only: plat, plon, plev

       implicit none

       !--------------------------------------------------------------
       ! ... arguments
       !--------------------------------------------------------------
       character(len=*), intent(in) :: name                   !  constituent name
       real(r8), intent(inout) :: q(plon,plev,plat)           !  mass mixing ratio

       !--------------------------------------------------------------
       ! ... local variables
       !--------------------------------------------------------------
       integer k, l
       real(r8) duma, dumb, dumz


       !
       ! to deactivate this routine, just return here
       !
       !     return


       if ( masterproc ) then
          write( *, '(2a)' )   &
               '*** modal_aero_initialize_q - name = ', name
          if (name == 'H2O2'   ) write( *, '(2a)' ) '    doing ', name
          if (name == 'SO2'    ) write( *, '(2a)' ) '    doing ', name
          if (name == 'H2SO4'  ) write( *, '(2a)' ) '    doing ', name
          if (name == 'DMS'    ) write( *, '(2a)' ) '    doing ', name
          if (name == 'NH3'    ) write( *, '(2a)' ) '    doing ', name
          if (name == 'so4_a1' ) write( *, '(2a)' ) '    doing ', name
          if (name == 'so4_a2' ) write( *, '(2a)' ) '    doing ', name
          if (name == 'pom_a3' ) write( *, '(2a)' ) '    doing ', name
          if (name == 'ncl_a4' ) write( *, '(2a)' ) '    doing ', name
          if (name == 'dst_a5' ) write( *, '(2a)' ) '    doing ', name
          if (name == 'ncl_a6' ) write( *, '(2a)' ) '    doing ', name
          if (name == 'dst_a7' ) write( *, '(2a)' ) '    doing ', name
       end if

       do k = 1, plev

          ! init gases
          dumz = (k+1.0e-5_r8)/(plev+1.0e-5_r8)
          dumb = dumz*1.0e-9_r8/28.966_r8
          if (name == 'H2O2'   ) q(:,k,:) = dumb*34.0_r8*1.0_r8
          if (name == 'SO2'    ) q(:,k,:) = dumb*64.0_r8*0.1_r8
          if (name == 'H2SO4'  ) q(:,k,:) = dumb*98.0_r8*0.001_r8
          if (name == 'DMS'    ) q(:,k,:) = dumb*62.0_r8*0.01_r8
          if (name == 'NH3'    ) q(:,k,:) = dumb*17.0_r8*0.1_r8

          ! init first mass species of each aerosol mode
          duma = dumz*1.0e-10_r8
          if (name == 'so4_a1' ) q(:,k,:) = duma*1.0_r8
          if (name == 'so4_a2' ) q(:,k,:) = duma*0.002_r8
          if (name == 'pom_a3' ) q(:,k,:) = duma*0.3_r8
          if (name == 'ncl_a4' ) q(:,k,:) = duma*0.4_r8
          if (name == 'dst_a5' ) q(:,k,:) = duma*0.5_r8
          if (name == 'ncl_a6' ) q(:,k,:) = duma*0.6_r8
          if (name == 'dst_a7' ) q(:,k,:) = duma*0.7_r8

          ! init aerosol number
          !
          ! at k=plev, duma = 1e-10 kgaero/kgair = 0.1 ugaero/kgair
          !            dumb = duma/(2000 kgaero/m3aero)
          duma = dumz*1.0e-10_r8
          dumb = duma/2.0e3_r8
          ! following produces number 1000X too small, and Dp 10X too big
          !        dumb = dumb*1.0e-3
          ! following produces number 1000X too big, and Dp 10X too small
          !        dumb = dumb*1.0e3
          if (name == 'num_a1' ) q(:,k,:) = dumb*1.0_r8  *3.0e20_r8
          if (name == 'num_a2' ) q(:,k,:) = dumb*0.002_r8*4.0e22_r8
          if (name == 'num_a3' ) q(:,k,:) = dumb*0.3_r8  *5.7e21_r8
          if (name == 'num_a4' ) q(:,k,:) = dumb*0.4_r8  *2.7e19_r8
          if (name == 'num_a5' ) q(:,k,:) = dumb*0.5_r8  *4.0e20_r8
          if (name == 'num_a6' ) q(:,k,:) = dumb*0.6_r8  *2.7e16_r8
          if (name == 'num_a7' ) q(:,k,:) = dumb*0.7_r8  *4.0e17_r8

          !*** modal_aero_calcsize_sub - ntot_amode    7
          !mode, dgn, dp*, v2n, v2nhi, v2nlo    1  1.100E-07  1.847E-07  3.031E+20  4.736E+18  2.635E+21
          !mode, dgn, dp*, v2n, v2nhi, v2nlo    2  2.600E-08  3.621E-08  4.021E+22  5.027E+21  1.073E+24
          !mode, dgn, dp*, v2n, v2nhi, v2nlo    3  5.000E-08  6.964E-08  5.654E+21  7.068E+20  7.068E+23
          !mode, dgn, dp*, v2n, v2nhi, v2nlo    4  2.000E-07  4.112E-07  2.748E+19  2.198E+17  1.758E+21
          !mode, dgn, dp*, v2n, v2nhi, v2nlo    5  1.000E-07  1.679E-07  4.035E+20  3.228E+18  3.228E+21
          !mode, dgn, dp*, v2n, v2nhi, v2nlo    6  2.000E-06  4.112E-06  2.748E+16  3.434E+15  2.198E+17
          !mode, dgn, dp*, v2n, v2nhi, v2nlo    7  1.000E-06  1.679E-06  4.035E+17  5.043E+16  3.228E+18

       end do   ! k

       if ( masterproc ) then
          write( *, '(7x,a,1p,10e10.2)' )   &
               name, (q(1,k,1), k=plev,1,-5) 
       end if

       if (plev > 0) return


       if ( masterproc ) then
          write( *, '(/a,i5)' )   &
               '*** modal_aero_initialize_q - ntot_amode', ntot_amode
          do k = 1, ntot_amode
             write( *, '(/a)' ) 'mode, dgn, v2n',   &
                  k, dgnum_amode(k), voltonumb_amode(k)
          end do
       end if

       return
     end subroutine modal_aero_initialize_q



     !==============================================================
   end module modal_aero_initialize_data

