module q3d_restart
   use shr_kind_mod, only: r8 => SHR_KIND_R8, CS => SHR_KIND_CS
   use ppgrid,       only: pcols, pver
   use pio,          only: var_desc_t

   implicit none
   private
   save

   ! Public interfaces
   public :: init_q3d_restart     ! Write the q3d restart attributes out
   public :: write_q3d_restart    ! Write the q3d restart info out
   public :: read_q3d_restart     ! Read the q3d restart info in

   ! Private variables
   type(var_desc_t) :: TH3D_desc
   type(var_desc_t) :: QV3D_desc
   type(var_desc_t) :: QC3D_desc
   type(var_desc_t) :: QI3D_desc
   type(var_desc_t) :: QR3D_desc
   type(var_desc_t) :: QS3D_desc
   type(var_desc_t) :: QG3D_desc
   type(var_desc_t) :: QT3D_desc
   type(var_desc_t) :: Z3DX_desc
   type(var_desc_t) :: Z3DY_desc
   type(var_desc_t) :: Z3DZ_desc
   type(var_desc_t) :: U3DX_desc
   type(var_desc_t) :: U3DY_desc
   type(var_desc_t) :: U3DX_CO_desc
   type(var_desc_t) :: U3DY_CO_desc
   type(var_desc_t) :: W3D_desc
   type(var_desc_t) :: W3DNM1_desc
   type(var_desc_t) :: PSI_desc
   type(var_desc_t) :: PSINM1_desc
   type(var_desc_t) :: CHI_desc
   type(var_desc_t) :: CHINM1_desc
   type(var_desc_t) :: TG_desc
   type(var_desc_t) :: ZROUGH_desc
   type(var_desc_t) :: GWET_desc
   type(var_desc_t) :: UW_desc
   type(var_desc_t) :: WV_desc
   type(var_desc_t) :: UW_CON_desc
   type(var_desc_t) :: WV_CON_desc
   type(var_desc_t) :: WTH_desc
   type(var_desc_t) :: WQV_desc
   type(var_desc_t) :: SPREC_desc

#ifdef JUNG_RESTART  
!  Since mapping is recalculated, these variables are also recalculated (No reading) 
   type(var_desc_t) :: COEFX_desc
   type(var_desc_t) :: COEFY_desc
   type(var_desc_t) :: COEFA_desc
   type(var_desc_t) :: COEFX_K_desc
   type(var_desc_t) :: COEFY_K_desc
   type(var_desc_t) :: COEFA_K_desc
   type(var_desc_t) :: COEFX_E_desc
   type(var_desc_t) :: COEFY_E_desc
   type(var_desc_t) :: COEFA_E_desc
   type(var_desc_t) :: COEFX_Z_desc
   type(var_desc_t) :: COEFY_Z_desc
   type(var_desc_t) :: COEFA_Z_desc
   type(var_desc_t) :: AGAU_CO_desc
   type(var_desc_t) :: BGAU_CO_desc
   type(var_desc_t) :: CGAU_CO_desc
   type(var_desc_t) :: RX2D_C1_desc
   type(var_desc_t) :: RX2D_C2_desc
#endif

   ! Variables to determine iodesc needs
   integer          :: ldof_ndims       =  0 ! Number of buffer dims in iodesc
   integer          :: ldof_lbound(4,4) =  0 ! Array lower bounds in iodesc
   integer          :: ldof_ubound(4,4) = -1 ! Array upper bounds in iodesc
   logical          :: need_decomp = .true.
   logical          :: have_decomp = .false.

   ! Private interfaces
   private q3d_chan_bounds
   private q3d_crm_decomp
   private q3d_free_decomp
   private crm_to_buffer
   private alloc_crm
   private buffer_to_crm
   private def_var
   private write_var
   private read_var

   interface q3d_chan_bounds
      module procedure q3d_chan_bounds_2d
      module procedure q3d_chan_bounds_3d
      module procedure q3d_chan_bounds_4d
   end interface q3d_chan_bounds

   interface crm_to_buffer
      module procedure crm_to_buffer_2d
      module procedure crm_to_buffer_3d
      module procedure crm_to_buffer_4d
   end interface crm_to_buffer

   interface buffer_to_crm
      module procedure buffer_to_crm_2d
      module procedure buffer_to_crm_3d
      module procedure buffer_to_crm_4d
   end interface buffer_to_crm

   interface def_var
      module procedure def_var_2d
      module procedure def_var_3d
      module procedure def_var_4d
   end interface def_var

   interface write_var
      module procedure write_var_2d
      module procedure write_var_3d
      module procedure write_var_4d
   end interface write_var

   interface read_var
      module procedure read_var_2d
      module procedure read_var_3d
      module procedure read_var_4d
   end interface read_var

!=======================================================================
CONTAINS
!=======================================================================

   subroutine init_q3d_restart(File)
      use cam_pio_utils,       only: cam_pio_def_dim, cam_pio_def_var
      use q3d_runtime,         only: q3d_total_channels
      use q3d_runtime,         only: q3d_begchan, q3d_endchan
      use vvm_data_types,      only: channel_t, channel_data
      use parmsld,             only: ntracer
      use pio,                 only: file_desc_t, PIO_DOUBLE
      use pio,                 only: pio_seterrorhandling, PIO_BCAST_ERROR

      !
      ! Dummy argument
      !
      type(file_desc_t), intent(inout) :: file
      !
      ! Local variables
      !
      integer                          :: numchanid, chanlenid
      integer                          :: hdimcnt, ierr, i
      integer                          :: dimids(4)
      integer, allocatable             :: hdimids(:)
      integer                          :: err_handling
      character(len=4)                 :: num
      real(r8),        pointer         :: buff1_2d(:,:)
      real(r8),        pointer         :: buff2_2d(:,:)
      real(r8),        pointer         :: buff3_2d(:,:)
      real(r8),        pointer         :: buff4_2d(:,:)
      real(r8),        pointer         :: buff1_3d(:,:,:)
      real(r8),        pointer         :: buff2_3d(:,:,:)
      real(r8),        pointer         :: buff3_3d(:,:,:)
      real(r8),        pointer         :: buff4_3d(:,:,:)
      real(r8),        pointer         :: buff1_4d(:,:,:,:)
      real(r8),        pointer         :: buff2_4d(:,:,:,:)
      real(r8),        pointer         :: buff3_4d(:,:,:,:)
      real(r8),        pointer         :: buff4_4d(:,:,:,:)
      type(channel_t), pointer         :: channel

      call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
      ! Do we need to save any variables from CONSTLD.F90 (should be no)?

      ! We need a dimension for the number of channels
      call cam_pio_def_dim(File, 'q3d_num_channels', q3d_total_channels,      &
           numchanid, existOK=.false.)

      nullify(buff1_2d)
      nullify(buff2_2d)
      nullify(buff3_2d)
      nullify(buff4_2d)
      nullify(buff1_3d)
      nullify(buff2_3d)
      nullify(buff3_3d)
      nullify(buff4_3d)
      nullify(buff1_4d)
      nullify(buff2_4d)
      nullify(buff3_4d)
      nullify(buff4_4d)
      ! Write attributes for CRM variables
      if (q3d_begchan <= q3d_endchan) then
         channel => channel_data(q3d_begchan)
      end if
      ! 2-D fields, isize, jsize
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WTH
         buff2_2d => channel%seg(2)%WTH
         buff3_2d => channel%seg(3)%WTH
         buff4_2d => channel%seg(4)%WTH
      end if
      call def_var(File, 'WTH', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, WTH_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WQV
         buff2_2d => channel%seg(2)%WQV
         buff3_2d => channel%seg(3)%WQV
         buff4_2d => channel%seg(4)%WQV
      end if
      call def_var(File, 'WQV', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, WQV_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%SPREC
         buff2_2d => channel%seg(2)%SPREC
         buff3_2d => channel%seg(3)%SPREC
         buff4_2d => channel%seg(4)%SPREC
      end if
      call def_var(File, 'SPREC', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, SPREC_desc)

#ifdef JUNG_RESTART       
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX
         buff2_2d => channel%seg(2)%COEFX
         buff3_2d => channel%seg(3)%COEFX
         buff4_2d => channel%seg(4)%COEFX
      end if
      call def_var(File, 'COEFX', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFX_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY
         buff2_2d => channel%seg(2)%COEFY
         buff3_2d => channel%seg(3)%COEFY
         buff4_2d => channel%seg(4)%COEFY
      end if
      call def_var(File, 'COEFY', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFY_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA
         buff2_2d => channel%seg(2)%COEFA
         buff3_2d => channel%seg(3)%COEFA
         buff4_2d => channel%seg(4)%COEFA
      end if
      call def_var(File, 'COEFA', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFA_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX_K
         buff2_2d => channel%seg(2)%COEFX_K
         buff3_2d => channel%seg(3)%COEFX_K
         buff4_2d => channel%seg(4)%COEFX_K
      end if
      call def_var(File, 'COEFX_K', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFX_K_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY_K
         buff2_2d => channel%seg(2)%COEFY_K
         buff3_2d => channel%seg(3)%COEFY_K
         buff4_2d => channel%seg(4)%COEFY_K
      end if
      call def_var(File, 'COEFY_K', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFY_K_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA_K
         buff2_2d => channel%seg(2)%COEFA_K
         buff3_2d => channel%seg(3)%COEFA_K
         buff4_2d => channel%seg(4)%COEFA_K
      end if
      call def_var(File, 'COEFA_K', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFA_K_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX_E
         buff2_2d => channel%seg(2)%COEFX_E
         buff3_2d => channel%seg(3)%COEFX_E
         buff4_2d => channel%seg(4)%COEFX_E
      end if
      call def_var(File, 'COEFX_E', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFX_E_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY_E
         buff2_2d => channel%seg(2)%COEFY_E
         buff3_2d => channel%seg(3)%COEFY_E
         buff4_2d => channel%seg(4)%COEFY_E
      end if
      call def_var(File, 'COEFY_E', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFY_E_desc)
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA_E
         buff2_2d => channel%seg(2)%COEFA_E
         buff3_2d => channel%seg(3)%COEFA_E
         buff4_2d => channel%seg(4)%COEFA_E
      end if
      call def_var(File, 'COEFA_E', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFA_E_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX_Z
         buff2_2d => channel%seg(2)%COEFX_Z
         buff3_2d => channel%seg(3)%COEFX_Z
         buff4_2d => channel%seg(4)%COEFX_Z
      end if
      call def_var(File, 'COEFX_Z', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFX_Z_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY_Z
         buff2_2d => channel%seg(2)%COEFY_Z
         buff3_2d => channel%seg(3)%COEFY_Z
         buff4_2d => channel%seg(4)%COEFY_Z
      end if
      call def_var(File, 'COEFY_Z', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFY_Z_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA_Z
         buff2_2d => channel%seg(2)%COEFA_Z
         buff3_2d => channel%seg(3)%COEFA_Z
         buff4_2d => channel%seg(4)%COEFA_Z
      end if
      call def_var(File, 'COEFA_Z', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, COEFA_Z_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%RX2D_C1
         buff2_2d => channel%seg(2)%RX2D_C1
         buff3_2d => channel%seg(3)%RX2D_C1
         buff4_2d => channel%seg(4)%RX2D_C1
      end if
      call def_var(File, 'RX2D_C1', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, RX2D_C1_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%RX2D_C2
         buff2_2d => channel%seg(2)%RX2D_C2
         buff3_2d => channel%seg(3)%RX2D_C2
         buff4_2d => channel%seg(4)%RX2D_C2
      end if
      call def_var(File, 'RX2D_C2', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, RX2D_C2_desc)
#endif
      
      ! 2-D fields, isize+1, jsize+1
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%UW
         buff2_2d => channel%seg(2)%UW
         buff3_2d => channel%seg(3)%UW
         buff4_2d => channel%seg(4)%UW
      end if
      call def_var(File, 'UW', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, UW_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WV
         buff2_2d => channel%seg(2)%WV
         buff3_2d => channel%seg(3)%WV
         buff4_2d => channel%seg(4)%WV
      end if
      call def_var(File, 'WV', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, WV_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%UW_CON
         buff2_2d => channel%seg(2)%UW_CON
         buff3_2d => channel%seg(3)%UW_CON
         buff4_2d => channel%seg(4)%UW_CON
      end if
      call def_var(File, 'UW_CON', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, UW_CON_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WV_CON
         buff2_2d => channel%seg(2)%WV_CON
         buff3_2d => channel%seg(3)%WV_CON
         buff4_2d => channel%seg(4)%WV_CON
      end if
      call def_var(File, 'WV_CON', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, WV_CON_desc)
      
      ! 2-D fields, full channel length
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%PSI
         buff2_2d => channel%seg(2)%PSI
         buff3_2d => channel%seg(3)%PSI
         buff4_2d => channel%seg(4)%PSI
      end if
      call def_var(File, 'PSI', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, PSI_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%PSINM1
         buff2_2d => channel%seg(2)%PSINM1
         buff3_2d => channel%seg(3)%PSINM1
         buff4_2d => channel%seg(4)%PSINM1
      end if
      call def_var(File, 'PSINM1', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, PSINM1_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%CHI
         buff2_2d => channel%seg(2)%CHI
         buff3_2d => channel%seg(3)%CHI
         buff4_2d => channel%seg(4)%CHI
      end if
      call def_var(File, 'CHI', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, CHI_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%CHINM1
         buff2_2d => channel%seg(2)%CHINM1
         buff3_2d => channel%seg(3)%CHINM1
         buff4_2d => channel%seg(4)%CHINM1
      end if
      call def_var(File, 'CHINM1', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, CHINM1_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%TG
         buff2_2d => channel%seg(2)%TG
         buff3_2d => channel%seg(3)%TG
         buff4_2d => channel%seg(4)%TG
      end if
      call def_var(File, 'TG', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, TG_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%ZROUGH
         buff2_2d => channel%seg(2)%ZROUGH
         buff3_2d => channel%seg(3)%ZROUGH
         buff4_2d => channel%seg(4)%ZROUGH
      end if
      call def_var(File, 'ZROUGH', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, ZROUGH_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%GWET
         buff2_2d => channel%seg(2)%GWET
         buff3_2d => channel%seg(3)%GWET
         buff4_2d => channel%seg(4)%GWET
      end if
      call def_var(File, 'GWET', buff1_2d, buff2_2d, buff3_2d, buff4_2d, numchanid, GWET_desc)
      
#ifdef JUNG_RESTART    
      ! 3-D, isize, jsize, nk1-1   
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%AGAU_CO
         buff2_3d => channel%seg(2)%AGAU_CO
         buff3_3d => channel%seg(3)%AGAU_CO
         buff4_3d => channel%seg(4)%AGAU_CO
      end if
      call def_var(File, 'AGAU_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, AGAU_CO_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%BGAU_CO
         buff2_3d => channel%seg(2)%BGAU_CO
         buff3_3d => channel%seg(3)%BGAU_CO
         buff4_3d => channel%seg(4)%BGAU_CO
      end if
      call def_var(File, 'BGAU_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, BGAU_CO_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%CGAU_CO
         buff2_3d => channel%seg(2)%CGAU_CO
         buff3_3d => channel%seg(3)%CGAU_CO
         buff4_3d => channel%seg(4)%CGAU_CO
      end if
      call def_var(File, 'CGAU_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, CGAU_CO_desc)
#endif      
      
      ! 3-D, nk2
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%Z3DX
         buff2_3d => channel%seg(2)%Z3DX
         buff3_3d => channel%seg(3)%Z3DX
         buff4_3d => channel%seg(4)%Z3DX
      end if
      call def_var(File, 'Z3DX', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, Z3DX_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%Z3DY
         buff2_3d => channel%seg(2)%Z3DY
         buff3_3d => channel%seg(3)%Z3DY
         buff4_3d => channel%seg(4)%Z3DY
      end if
      call def_var(File, 'Z3DY', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, Z3DY_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%W3D
         buff2_3d => channel%seg(2)%W3D
         buff3_3d => channel%seg(3)%W3D
         buff4_3d => channel%seg(4)%W3D
      end if
      call def_var(File, 'W3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, W3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%W3DNM1
         buff2_3d => channel%seg(2)%W3DNM1
         buff3_3d => channel%seg(3)%W3DNM1
         buff4_3d => channel%seg(4)%W3DNM1
      end if
      call def_var(File, 'W3DNM1', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, W3DNM1_desc)
      
      ! 3-D, nk3
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%TH3d
         buff2_3d => channel%seg(2)%TH3d
         buff3_3d => channel%seg(3)%TH3d
         buff4_3d => channel%seg(4)%TH3d
      end if
      call def_var(File, 'TH3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, TH3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QV3D
         buff2_3d => channel%seg(2)%QV3D
         buff3_3d => channel%seg(3)%QV3D
         buff4_3d => channel%seg(4)%QV3D
      end if
      call def_var(File, 'QV3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, QV3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QC3D
         buff2_3d => channel%seg(2)%QC3D
         buff3_3d => channel%seg(3)%QC3D
         buff4_3d => channel%seg(4)%QC3D
      end if
      call def_var(File, 'QC3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, QC3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QI3D
         buff2_3d => channel%seg(2)%QI3D
         buff3_3d => channel%seg(3)%QI3D
         buff4_3d => channel%seg(4)%QI3D
      end if
      call def_var(File, 'QI3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, QI3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QR3D
         buff2_3d => channel%seg(2)%QR3D
         buff3_3d => channel%seg(3)%QR3D
         buff4_3d => channel%seg(4)%QR3D
      end if
      call def_var(File, 'QR3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, QR3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QS3D
         buff2_3d => channel%seg(2)%QS3D
         buff3_3d => channel%seg(3)%QS3D
         buff4_3d => channel%seg(4)%QS3D
      end if
      call def_var(File, 'QS3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, QS3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QG3D
         buff2_3d => channel%seg(2)%QG3D
         buff3_3d => channel%seg(3)%QG3D
         buff4_3d => channel%seg(4)%QG3D
      end if
      call def_var(File, 'QG3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, QG3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%Z3DZ
         buff2_3d => channel%seg(2)%Z3DZ
         buff3_3d => channel%seg(3)%Z3DZ
         buff4_3d => channel%seg(4)%Z3DZ
      end if
      call def_var(File, 'Z3DZ', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, Z3DZ_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DX
         buff2_3d => channel%seg(2)%U3DX
         buff3_3d => channel%seg(3)%U3DX
         buff4_3d => channel%seg(4)%U3DX
      end if
      call def_var(File, 'U3DX', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, U3DX_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DY
         buff2_3d => channel%seg(2)%U3DY
         buff3_3d => channel%seg(3)%U3DY
         buff4_3d => channel%seg(4)%U3DY
      end if
      call def_var(File, 'U3DY', buff1_3d, buff2_3d, buff3_3d, buff4_3d,numchanid, U3DY_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DX_CO
         buff2_3d => channel%seg(2)%U3DX_CO
         buff3_3d => channel%seg(3)%U3DX_CO
         buff4_3d => channel%seg(4)%U3DX_CO
      end if
      call def_var(File, 'U3DX_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, U3DX_CO_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DY_CO
         buff2_3d => channel%seg(2)%U3DY_CO
         buff3_3d => channel%seg(3)%U3DY_CO
         buff4_3d => channel%seg(4)%U3DY_CO
      end if
      call def_var(File, 'U3DY_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, numchanid, U3DY_CO_desc)
      
      ! 4-D Field
      if (q3d_begchan <= q3d_endchan) then
         buff1_4d => channel%seg(1)%QT3D
         buff2_4d => channel%seg(2)%QT3D
         buff3_4d => channel%seg(3)%QT3D
         buff4_4d => channel%seg(4)%QT3D
      end if
      if (ntracer > 0) then
         call def_var(File, 'QT3D', buff1_4d, buff2_4d, buff3_4d, buff4_4d, numchanid, QT3D_desc)
      end if

      call pio_seterrorhandling(File, err_handling)

   end subroutine init_q3d_restart

   subroutine write_q3d_restart(File)
      !-----------------------------------------------------------------------
      use cam_history_support, only: fillvalue
      use spmd_utils,          only: iam
      use pio,                 only: file_desc_t, io_desc_t
      use vvm_data_types,      only: channel_t, channel_data
      use q3d_runtime,         only: q3d_begchan, q3d_endchan
      use parmsld,             only: ntracer
      !
      ! Dummy argument
      !
      type(file_desc_t), intent(inout) :: File
      !
      ! Local variables
      !
      type(io_desc_t)                  :: iodesc
      type(channel_t), pointer         :: channel
      real(r8),        pointer         :: buff1_2d(:,:)
      real(r8),        pointer         :: buff2_2d(:,:)
      real(r8),        pointer         :: buff3_2d(:,:)
      real(r8),        pointer         :: buff4_2d(:,:)
      real(r8),        pointer         :: buff1_3d(:,:,:)
      real(r8),        pointer         :: buff2_3d(:,:,:)
      real(r8),        pointer         :: buff3_3d(:,:,:)
      real(r8),        pointer         :: buff4_3d(:,:,:)
      real(r8),        pointer         :: buff1_4d(:,:,:,:)
      real(r8),        pointer         :: buff2_4d(:,:,:,:)
      real(r8),        pointer         :: buff3_4d(:,:,:,:)
      real(r8),        pointer         :: buff4_4d(:,:,:,:)

      nullify(buff1_2d)
      nullify(buff2_2d)
      nullify(buff3_2d)
      nullify(buff4_2d)
      nullify(buff1_3d)
      nullify(buff2_3d)
      nullify(buff3_3d)
      nullify(buff4_3d)
      nullify(buff1_4d)
      nullify(buff2_4d)
      nullify(buff3_4d)
      nullify(buff4_4d)
      if (q3d_begchan <= q3d_endchan) then
         channel => channel_data(q3d_begchan)
      end if
      ! Reset any decomps
      ldof_lbound(:,:) =  0
      ldof_ubound(:,:) = -1
      need_decomp = .true.
      ! Write CRM values to restart file, sort by levels and length for efficiency
      
      ! 2-D fields, isize, jsize
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WTH
         buff2_2d => channel%seg(2)%WTH
         buff3_2d => channel%seg(3)%WTH
         buff4_2d => channel%seg(4)%WTH
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, WTH_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WQV
         buff2_2d => channel%seg(2)%WQV
         buff3_2d => channel%seg(3)%WQV
         buff4_2d => channel%seg(4)%WQV
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, WQV_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%SPREC
         buff2_2d => channel%seg(2)%SPREC
         buff3_2d => channel%seg(3)%SPREC
         buff4_2d => channel%seg(4)%SPREC
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, SPREC_desc)
      
#ifdef JUNG_RESTART      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX
         buff2_2d => channel%seg(2)%COEFX
         buff3_2d => channel%seg(3)%COEFX
         buff4_2d => channel%seg(4)%COEFX
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFX_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY
         buff2_2d => channel%seg(2)%COEFY
         buff3_2d => channel%seg(3)%COEFY
         buff4_2d => channel%seg(4)%COEFY
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFY_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA
         buff2_2d => channel%seg(2)%COEFA
         buff3_2d => channel%seg(3)%COEFA
         buff4_2d => channel%seg(4)%COEFA
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFA_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX_K
         buff2_2d => channel%seg(2)%COEFX_K
         buff3_2d => channel%seg(3)%COEFX_K
         buff4_2d => channel%seg(4)%COEFX_K
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFX_K_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY_K
         buff2_2d => channel%seg(2)%COEFY_K
         buff3_2d => channel%seg(3)%COEFY_K
         buff4_2d => channel%seg(4)%COEFY_K
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFY_K_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA_K
         buff2_2d => channel%seg(2)%COEFA_K
         buff3_2d => channel%seg(3)%COEFA_K
         buff4_2d => channel%seg(4)%COEFA_K
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFA_K_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX_E
         buff2_2d => channel%seg(2)%COEFX_E
         buff3_2d => channel%seg(3)%COEFX_E
         buff4_2d => channel%seg(4)%COEFX_E
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFX_E_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY_E
         buff2_2d => channel%seg(2)%COEFY_E
         buff3_2d => channel%seg(3)%COEFY_E
         buff4_2d => channel%seg(4)%COEFY_E
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFY_E_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA_E
         buff2_2d => channel%seg(2)%COEFA_E
         buff3_2d => channel%seg(3)%COEFA_E
         buff4_2d => channel%seg(4)%COEFA_E
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFA_E_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX_Z
         buff2_2d => channel%seg(2)%COEFX_Z
         buff3_2d => channel%seg(3)%COEFX_Z
         buff4_2d => channel%seg(4)%COEFX_Z
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFX_Z_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY_Z
         buff2_2d => channel%seg(2)%COEFY_Z
         buff3_2d => channel%seg(3)%COEFY_Z
         buff4_2d => channel%seg(4)%COEFY_Z
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFY_Z_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA_Z
         buff2_2d => channel%seg(2)%COEFA_Z
         buff3_2d => channel%seg(3)%COEFA_Z
         buff4_2d => channel%seg(4)%COEFA_Z
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, COEFA_Z_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%RX2D_C1
         buff2_2d => channel%seg(2)%RX2D_C1
         buff3_2d => channel%seg(3)%RX2D_C1
         buff4_2d => channel%seg(4)%RX2D_C1
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, RX2D_C1_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%RX2D_C2
         buff2_2d => channel%seg(2)%RX2D_C2
         buff3_2d => channel%seg(3)%RX2D_C2
         buff4_2d => channel%seg(4)%RX2D_C2
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, RX2D_C2_desc)  
#endif
      
      call q3d_free_decomp(File, iodesc)
      
      ! 2-D fields, isize+1, jsize+1
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%UW
         buff2_2d => channel%seg(2)%UW
         buff3_2d => channel%seg(3)%UW
         buff4_2d => channel%seg(4)%UW
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, UW_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WV
         buff2_2d => channel%seg(2)%WV
         buff3_2d => channel%seg(3)%WV
         buff4_2d => channel%seg(4)%WV
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, WV_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%UW_CON
         buff2_2d => channel%seg(2)%UW_CON
         buff3_2d => channel%seg(3)%UW_CON
         buff4_2d => channel%seg(4)%UW_CON
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, UW_CON_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WV_CON
         buff2_2d => channel%seg(2)%WV_CON
         buff3_2d => channel%seg(3)%WV_CON
         buff4_2d => channel%seg(4)%WV_CON
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, WV_CON_desc)
      
      call q3d_free_decomp(File, iodesc)
      
      ! 2-D fields, full channel length
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%PSI
         buff2_2d => channel%seg(2)%PSI
         buff3_2d => channel%seg(3)%PSI
         buff4_2d => channel%seg(4)%PSI
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, PSI_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%PSINM1
         buff2_2d => channel%seg(2)%PSINM1
         buff3_2d => channel%seg(3)%PSINM1
         buff4_2d => channel%seg(4)%PSINM1
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, PSINM1_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%CHI
         buff2_2d => channel%seg(2)%CHI
         buff3_2d => channel%seg(3)%CHI
         buff4_2d => channel%seg(4)%CHI
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, CHI_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%CHINM1
         buff2_2d => channel%seg(2)%CHINM1
         buff3_2d => channel%seg(3)%CHINM1
         buff4_2d => channel%seg(4)%CHINM1
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, CHINM1_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%TG
         buff2_2d => channel%seg(2)%TG
         buff3_2d => channel%seg(3)%TG
         buff4_2d => channel%seg(4)%TG
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, TG_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%ZROUGH
         buff2_2d => channel%seg(2)%ZROUGH
         buff3_2d => channel%seg(3)%ZROUGH
         buff4_2d => channel%seg(4)%ZROUGH
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, ZROUGH_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%GWET
         buff2_2d => channel%seg(2)%GWET
         buff3_2d => channel%seg(3)%GWET
         buff4_2d => channel%seg(4)%GWET
      end if
      call write_var(File, buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc, GWET_desc)
      
      call q3d_free_decomp(File, iodesc)
      
#ifdef JUNG_RESTART   
      ! 3-D, isize, jsize, nk1-1    
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%AGAU_CO
         buff2_3d => channel%seg(2)%AGAU_CO
         buff3_3d => channel%seg(3)%AGAU_CO
         buff4_3d => channel%seg(4)%AGAU_CO
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, AGAU_CO_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%BGAU_CO
         buff2_3d => channel%seg(2)%BGAU_CO
         buff3_3d => channel%seg(3)%BGAU_CO
         buff4_3d => channel%seg(4)%BGAU_CO
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, BGAU_CO_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%CGAU_CO
         buff2_3d => channel%seg(2)%CGAU_CO
         buff3_3d => channel%seg(3)%CGAU_CO
         buff4_3d => channel%seg(4)%CGAU_CO
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, CGAU_CO_desc)
      
      call q3d_free_decomp(File, iodesc)
#endif      
      
      ! 3-D, nk2
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%Z3DX
         buff2_3d => channel%seg(2)%Z3DX
         buff3_3d => channel%seg(3)%Z3DX
         buff4_3d => channel%seg(4)%Z3DX
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, Z3DX_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%Z3DY
         buff2_3d => channel%seg(2)%Z3DY
         buff3_3d => channel%seg(3)%Z3DY
         buff4_3d => channel%seg(4)%Z3DY
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, Z3DY_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%W3D
         buff2_3d => channel%seg(2)%W3D
         buff3_3d => channel%seg(3)%W3D
         buff4_3d => channel%seg(4)%W3D
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, W3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%W3DNM1
         buff2_3d => channel%seg(2)%W3DNM1
         buff3_3d => channel%seg(3)%W3DNM1
         buff4_3d => channel%seg(4)%W3DNM1
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, W3DNM1_desc)
      
      call q3d_free_decomp(File, iodesc)
      
      ! 3-D, nk3
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%TH3d
         buff2_3d => channel%seg(2)%TH3d
         buff3_3d => channel%seg(3)%TH3d
         buff4_3d => channel%seg(4)%TH3d
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, TH3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QV3D
         buff2_3d => channel%seg(2)%QV3D
         buff3_3d => channel%seg(3)%QV3D
         buff4_3d => channel%seg(4)%QV3D
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, QV3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QC3D
         buff2_3d => channel%seg(2)%QC3D
         buff3_3d => channel%seg(3)%QC3D
         buff4_3d => channel%seg(4)%QC3D
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, QC3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QI3D
         buff2_3d => channel%seg(2)%QI3D
         buff3_3d => channel%seg(3)%QI3D
         buff4_3d => channel%seg(4)%QI3D
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, QI3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QR3D
         buff2_3d => channel%seg(2)%QR3D
         buff3_3d => channel%seg(3)%QR3D
         buff4_3d => channel%seg(4)%QR3D
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, QR3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QS3D
         buff2_3d => channel%seg(2)%QS3D
         buff3_3d => channel%seg(3)%QS3D
         buff4_3d => channel%seg(4)%QS3D
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, QS3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QG3D
         buff2_3d => channel%seg(2)%QG3D
         buff3_3d => channel%seg(3)%QG3D
         buff4_3d => channel%seg(4)%QG3D
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, QG3D_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%Z3DZ
         buff2_3d => channel%seg(2)%Z3DZ
         buff3_3d => channel%seg(3)%Z3DZ
         buff4_3d => channel%seg(4)%Z3DZ
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, Z3DZ_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DX
         buff2_3d => channel%seg(2)%U3DX
         buff3_3d => channel%seg(3)%U3DX
         buff4_3d => channel%seg(4)%U3DX
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, U3DX_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DY
         buff2_3d => channel%seg(2)%U3DY
         buff3_3d => channel%seg(3)%U3DY
         buff4_3d => channel%seg(4)%U3DY
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, U3DY_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DX_CO
         buff2_3d => channel%seg(2)%U3DX_CO
         buff3_3d => channel%seg(3)%U3DX_CO
         buff4_3d => channel%seg(4)%U3DX_CO
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, U3DX_CO_desc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DY_CO
         buff2_3d => channel%seg(2)%U3DY_CO
         buff3_3d => channel%seg(3)%U3DY_CO
         buff4_3d => channel%seg(4)%U3DY_CO
      end if
      call write_var(File, buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc, U3DY_CO_desc)
      
      call q3d_free_decomp(File, iodesc)
      
      ! 4-D Field
      if (ntracer > 0) then
         if (q3d_begchan <= q3d_endchan) then
            buff1_4d => channel%seg(1)%QT3D
            buff2_4d => channel%seg(2)%QT3D
            buff3_4d => channel%seg(3)%QT3D
            buff4_4d => channel%seg(4)%QT3D
         end if
         call write_var(File, buff1_4d, buff2_4d, buff3_4d, buff4_4d, iodesc, QT3D_desc)
      end if

      ! Cleanup
      call q3d_free_decomp(File, iodesc)

   end subroutine write_q3d_restart

   subroutine read_q3d_restart(File)
      !-----------------------------------------------------------------------
      use vvm_data_types,      only: channel_t, channel_data
      use pio,                 only: file_desc_t, io_desc_t
      use q3d_runtime,         only: q3d_begchan, q3d_endchan
      use parmsld,             only: ntracer
      !
      ! Dummy argument
      !
      type(file_desc_t),   intent(inout) :: File
      !
      ! Local variables
      !
      type(io_desc_t)                  :: iodesc
      type(channel_t), pointer         :: channel
      real(r8),        pointer         :: buff1_2d(:,:)
      real(r8),        pointer         :: buff2_2d(:,:)
      real(r8),        pointer         :: buff3_2d(:,:)
      real(r8),        pointer         :: buff4_2d(:,:)
      real(r8),        pointer         :: buff1_3d(:,:,:)
      real(r8),        pointer         :: buff2_3d(:,:,:)
      real(r8),        pointer         :: buff3_3d(:,:,:)
      real(r8),        pointer         :: buff4_3d(:,:,:)
      real(r8),        pointer         :: buff1_4d(:,:,:,:)
      real(r8),        pointer         :: buff2_4d(:,:,:,:)
      real(r8),        pointer         :: buff3_4d(:,:,:,:)
      real(r8),        pointer         :: buff4_4d(:,:,:,:)

      nullify(buff1_2d)
      nullify(buff2_2d)
      nullify(buff3_2d)
      nullify(buff4_2d)
      nullify(buff1_3d)
      nullify(buff2_3d)
      nullify(buff3_3d)
      nullify(buff4_3d)
      nullify(buff1_4d)
      nullify(buff2_4d)
      nullify(buff3_4d)
      nullify(buff4_4d)
      
      if (q3d_begchan <= q3d_endchan) then
         channel => channel_data(q3d_begchan)
      end if
      
      ! Reset any decomps
      ldof_lbound(:,:) =  0
      ldof_ubound(:,:) = -1
      
      ! Write CRM values to restart file, sort by levels and length for efficiency
      
      ! 2-D fields, isize, jsize
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WTH
         buff2_2d => channel%seg(2)%WTH
         buff3_2d => channel%seg(3)%WTH
         buff4_2d => channel%seg(4)%WTH
      end if
      call read_var(File, 'WTH', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WQV
         buff2_2d => channel%seg(2)%WQV
         buff3_2d => channel%seg(3)%WQV
         buff4_2d => channel%seg(4)%WQV
      end if
      call read_var(File, 'WQV', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%SPREC
         buff2_2d => channel%seg(2)%SPREC
         buff3_2d => channel%seg(3)%SPREC
         buff4_2d => channel%seg(4)%SPREC
      end if
      call read_var(File, 'SPREC', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
#ifdef JUNG_RESTART        
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX
         buff2_2d => channel%seg(2)%COEFX
         buff3_2d => channel%seg(3)%COEFX
         buff4_2d => channel%seg(4)%COEFX
      end if
      call read_var(File, 'COEFX', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY
         buff2_2d => channel%seg(2)%COEFY
         buff3_2d => channel%seg(3)%COEFY
         buff4_2d => channel%seg(4)%COEFY
      end if
      call read_var(File, 'COEFY', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA
         buff2_2d => channel%seg(2)%COEFA
         buff3_2d => channel%seg(3)%COEFA
         buff4_2d => channel%seg(4)%COEFA
      end if
      call read_var(File, 'COEFA', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX_K
         buff2_2d => channel%seg(2)%COEFX_K
         buff3_2d => channel%seg(3)%COEFX_K
         buff4_2d => channel%seg(4)%COEFX_K
      end if
      call read_var(File, 'COEFX_K', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY_K
         buff2_2d => channel%seg(2)%COEFY_K
         buff3_2d => channel%seg(3)%COEFY_K
         buff4_2d => channel%seg(4)%COEFY_K
      end if
      call read_var(File, 'COEFY_K', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA_K
         buff2_2d => channel%seg(2)%COEFA_K
         buff3_2d => channel%seg(3)%COEFA_K
         buff4_2d => channel%seg(4)%COEFA_K
      end if
      call read_var(File, 'COEFA_K', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX_E
         buff2_2d => channel%seg(2)%COEFX_E
         buff3_2d => channel%seg(3)%COEFX_E
         buff4_2d => channel%seg(4)%COEFX_E
      end if
      call read_var(File, 'COEFX_E', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY_E
         buff2_2d => channel%seg(2)%COEFY_E
         buff3_2d => channel%seg(3)%COEFY_E
         buff4_2d => channel%seg(4)%COEFY_E
      end if
      call read_var(File, 'COEFY_E', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA_E
         buff2_2d => channel%seg(2)%COEFA_E
         buff3_2d => channel%seg(3)%COEFA_E
         buff4_2d => channel%seg(4)%COEFA_E
      end if
      call read_var(File, 'COEFA_E', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFX_Z
         buff2_2d => channel%seg(2)%COEFX_Z
         buff3_2d => channel%seg(3)%COEFX_Z
         buff4_2d => channel%seg(4)%COEFX_Z
      end if
      call read_var(File, 'COEFX_Z', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFY_Z
         buff2_2d => channel%seg(2)%COEFY_Z
         buff3_2d => channel%seg(3)%COEFY_Z
         buff4_2d => channel%seg(4)%COEFY_Z
      end if
      call read_var(File, 'COEFY_Z', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%COEFA_Z
         buff2_2d => channel%seg(2)%COEFA_Z
         buff3_2d => channel%seg(3)%COEFA_Z
         buff4_2d => channel%seg(4)%COEFA_Z
      end if
      call read_var(File, 'COEFA_Z', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%RX2D_C1
         buff2_2d => channel%seg(2)%RX2D_C1
         buff3_2d => channel%seg(3)%RX2D_C1
         buff4_2d => channel%seg(4)%RX2D_C1
      end if
      call read_var(File, 'RX2D_C1', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%RX2D_C2
         buff2_2d => channel%seg(2)%RX2D_C2
         buff3_2d => channel%seg(3)%RX2D_C2
         buff4_2d => channel%seg(4)%RX2D_C2
      end if
      call read_var(File, 'RX2D_C2', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
#endif      
      
      call q3d_free_decomp(File, iodesc)
      
      ! 2-D fields, isize+1, jsize+1
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%UW
         buff2_2d => channel%seg(2)%UW
         buff3_2d => channel%seg(3)%UW
         buff4_2d => channel%seg(4)%UW
      end if
      call read_var(File, 'UW', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WV
         buff2_2d => channel%seg(2)%WV
         buff3_2d => channel%seg(3)%WV
         buff4_2d => channel%seg(4)%WV
      end if
      call read_var(File, 'WV', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%UW_CON
         buff2_2d => channel%seg(2)%UW_CON
         buff3_2d => channel%seg(3)%UW_CON
         buff4_2d => channel%seg(4)%UW_CON
      end if
      call read_var(File, 'UW_CON', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%WV_CON
         buff2_2d => channel%seg(2)%WV_CON
         buff3_2d => channel%seg(3)%WV_CON
         buff4_2d => channel%seg(4)%WV_CON
      end if
      call read_var(File, 'WV_CON', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      call q3d_free_decomp(File, iodesc)
      
      ! 2-D fields, full channel length
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%PSI
         buff2_2d => channel%seg(2)%PSI
         buff3_2d => channel%seg(3)%PSI
         buff4_2d => channel%seg(4)%PSI
      end if
      call read_var(File, 'PSI', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%PSINM1
         buff2_2d => channel%seg(2)%PSINM1
         buff3_2d => channel%seg(3)%PSINM1
         buff4_2d => channel%seg(4)%PSINM1
      end if
      call read_var(File, 'PSINM1', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%CHI
         buff2_2d => channel%seg(2)%CHI
         buff3_2d => channel%seg(3)%CHI
         buff4_2d => channel%seg(4)%CHI
      end if
      call read_var(File, 'CHI', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%CHINM1
         buff2_2d => channel%seg(2)%CHINM1
         buff3_2d => channel%seg(3)%CHINM1
         buff4_2d => channel%seg(4)%CHINM1
      end if
      call read_var(File, 'CHINM1', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%TG
         buff2_2d => channel%seg(2)%TG
         buff3_2d => channel%seg(3)%TG
         buff4_2d => channel%seg(4)%TG
      end if
      call read_var(File, 'TG', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%ZROUGH
         buff2_2d => channel%seg(2)%ZROUGH
         buff3_2d => channel%seg(3)%ZROUGH
         buff4_2d => channel%seg(4)%ZROUGH
      end if
      call read_var(File, 'ZROUGH', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_2d => channel%seg(1)%GWET
         buff2_2d => channel%seg(2)%GWET
         buff3_2d => channel%seg(3)%GWET
         buff4_2d => channel%seg(4)%GWET
      end if
      call read_var(File, 'GWET', buff1_2d, buff2_2d, buff3_2d, buff4_2d, iodesc)
      
      call q3d_free_decomp(File, iodesc)

#ifdef JUNG_RESTART        
      ! 3-D, isize, jsize, nk1-1
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%AGAU_CO
         buff2_3d => channel%seg(2)%AGAU_CO
         buff3_3d => channel%seg(3)%AGAU_CO
         buff4_3d => channel%seg(4)%AGAU_CO
      end if
      call read_var(File, 'AGAU_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%BGAU_CO
         buff2_3d => channel%seg(2)%BGAU_CO
         buff3_3d => channel%seg(3)%BGAU_CO
         buff4_3d => channel%seg(4)%BGAU_CO
      end if
      call read_var(File, 'BGAU_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%CGAU_CO
         buff2_3d => channel%seg(2)%CGAU_CO
         buff3_3d => channel%seg(3)%CGAU_CO
         buff4_3d => channel%seg(4)%CGAU_CO
      end if
      call read_var(File, 'CGAU_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      call q3d_free_decomp(File, iodesc)
#endif
      
      ! 3-D, nk2
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%Z3DX
         buff2_3d => channel%seg(2)%Z3DX
         buff3_3d => channel%seg(3)%Z3DX
         buff4_3d => channel%seg(4)%Z3DX
      end if
      call read_var(File, 'Z3DX', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%Z3DY
         buff2_3d => channel%seg(2)%Z3DY
         buff3_3d => channel%seg(3)%Z3DY
         buff4_3d => channel%seg(4)%Z3DY
      end if
      call read_var(File, 'Z3DY', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%W3D
         buff2_3d => channel%seg(2)%W3D
         buff3_3d => channel%seg(3)%W3D
         buff4_3d => channel%seg(4)%W3D
      end if
      call read_var(File, 'W3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%W3DNM1
         buff2_3d => channel%seg(2)%W3DNM1
         buff3_3d => channel%seg(3)%W3DNM1
         buff4_3d => channel%seg(4)%W3DNM1
      end if
      call read_var(File, 'W3DNM1', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      call q3d_free_decomp(File, iodesc)
      ! 3-D, nk3
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%TH3d
         buff2_3d => channel%seg(2)%TH3d
         buff3_3d => channel%seg(3)%TH3d
         buff4_3d => channel%seg(4)%TH3d
      end if
      call read_var(File, 'TH3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QV3D
         buff2_3d => channel%seg(2)%QV3D
         buff3_3d => channel%seg(3)%QV3D
         buff4_3d => channel%seg(4)%QV3D
      end if
      call read_var(File, 'QV3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QC3D
         buff2_3d => channel%seg(2)%QC3D
         buff3_3d => channel%seg(3)%QC3D
         buff4_3d => channel%seg(4)%QC3D
      end if
      call read_var(File, 'QC3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QI3D
         buff2_3d => channel%seg(2)%QI3D
         buff3_3d => channel%seg(3)%QI3D
         buff4_3d => channel%seg(4)%QI3D
      end if
      call read_var(File, 'QI3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QR3D
         buff2_3d => channel%seg(2)%QR3D
         buff3_3d => channel%seg(3)%QR3D
         buff4_3d => channel%seg(4)%QR3D
      end if
      call read_var(File, 'QR3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QS3D
         buff2_3d => channel%seg(2)%QS3D
         buff3_3d => channel%seg(3)%QS3D
         buff4_3d => channel%seg(4)%QS3D
      end if
      call read_var(File, 'QS3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%QG3D
         buff2_3d => channel%seg(2)%QG3D
         buff3_3d => channel%seg(3)%QG3D
         buff4_3d => channel%seg(4)%QG3D
      end if
      call read_var(File, 'QG3D', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%Z3DZ
         buff2_3d => channel%seg(2)%Z3DZ
         buff3_3d => channel%seg(3)%Z3DZ
         buff4_3d => channel%seg(4)%Z3DZ
      end if
      call read_var(File, 'Z3DZ', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DX
         buff2_3d => channel%seg(2)%U3DX
         buff3_3d => channel%seg(3)%U3DX
         buff4_3d => channel%seg(4)%U3DX
      end if
      call read_var(File, 'U3DX', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DY
         buff2_3d => channel%seg(2)%U3DY
         buff3_3d => channel%seg(3)%U3DY
         buff4_3d => channel%seg(4)%U3DY
      end if
      call read_var(File, 'U3DY', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DX_CO
         buff2_3d => channel%seg(2)%U3DX_CO
         buff3_3d => channel%seg(3)%U3DX_CO
         buff4_3d => channel%seg(4)%U3DX_CO
      end if
      call read_var(File, 'U3DX_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      if (q3d_begchan <= q3d_endchan) then
         buff1_3d => channel%seg(1)%U3DY_CO
         buff2_3d => channel%seg(2)%U3DY_CO
         buff3_3d => channel%seg(3)%U3DY_CO
         buff4_3d => channel%seg(4)%U3DY_CO
      end if
      call read_var(File, 'U3DY_CO', buff1_3d, buff2_3d, buff3_3d, buff4_3d, iodesc)
      
      call q3d_free_decomp(File, iodesc)
      
      ! 4-D Field
      if (ntracer > 0) then
         if (q3d_begchan <= q3d_endchan) then
            buff1_4d => channel%seg(1)%QT3D
            buff2_4d => channel%seg(2)%QT3D
            buff3_4d => channel%seg(3)%QT3D
            buff4_4d => channel%seg(4)%QT3D
         end if
         call read_var(File, 'QT3D', buff1_4d, buff2_4d, buff3_4d, buff4_4d, iodesc)
      end if

      ! Cleanup
      call q3d_free_decomp(File, iodesc)

   end subroutine read_q3d_restart

   subroutine def_var_2d(File, name, buff1, buff2, buff3, buff4, numchanid, vardesc)
      use cam_pio_utils,       only: cam_pio_def_dim, cam_pio_def_var
      use pio,                 only: file_desc_t, PIO_DOUBLE

      type(file_desc_t), intent(inout) :: file
      character(len=*),  intent(in)    :: name
      real(r8), pointer                :: buff1(:,:)
      real(r8), pointer                :: buff2(:,:)
      real(r8), pointer                :: buff3(:,:)
      real(r8), pointer                :: buff4(:,:)
      integer                          :: numchanid
      type(var_desc_t),  intent(out)   :: vardesc

      integer                          :: lbounds(4,4)
      integer                          :: ubounds(4,4)
      integer                          :: chan_len ! Channel buffer size
      integer                          :: clenid   ! Channel size dimid
      character(len=16)                :: chan_dimname

      ! Make sure we have a dimension for the channel length
      call q3d_chan_bounds(buff1, buff2, buff3, buff4,                        &
           lbounds, ubounds, chan_size=chan_len)
      write(chan_dimname, '("q3d_chan_",i0)') chan_len
      call cam_pio_def_dim(File, trim(chan_dimname), chan_len,                &
           clenid, existOK=.true.)

      call cam_pio_def_var(File, name, PIO_DOUBLE, (/ clenid, numchanid /),   &
           vardesc)

   end subroutine def_var_2d

   subroutine def_var_3d(File, name, buff1, buff2, buff3, buff4, numchanid, vardesc)
      use cam_pio_utils,       only: cam_pio_def_dim, cam_pio_def_var
      use pio,                 only: file_desc_t, PIO_DOUBLE

      type(file_desc_t), intent(inout) :: file
      character(len=*),  intent(in)    :: name
      real(r8), pointer                :: buff1(:,:,:)
      real(r8), pointer                :: buff2(:,:,:)
      real(r8), pointer                :: buff3(:,:,:)
      real(r8), pointer                :: buff4(:,:,:)
      integer                          :: numchanid
      type(var_desc_t),  intent(out)   :: vardesc

      integer                          :: lbounds(4,4)
      integer                          :: ubounds(4,4)
      integer                          :: chan_len ! Channel buffer size
      integer                          :: clenid   ! Channel size dimid
      character(len=16)                :: chan_dimname

      ! Make sure we have a dimension for the channel length
      call q3d_chan_bounds(buff1, buff2, buff3, buff4,                        &
           lbounds, ubounds, chan_size=chan_len)
      write(chan_dimname, '("q3d_chan_",i0)') chan_len
      call cam_pio_def_dim(File, trim(chan_dimname), chan_len,                &
           clenid, existOK=.true.)

      call cam_pio_def_var(File, name, PIO_DOUBLE, (/ clenid, numchanid /),   &
           vardesc)

   end subroutine def_var_3d

   subroutine def_var_4d(File, name, buff1, buff2, buff3, buff4, numchanid, vardesc)
      use cam_pio_utils,       only: cam_pio_def_dim, cam_pio_def_var
      use pio,                 only: file_desc_t, PIO_DOUBLE

      type(file_desc_t), intent(inout) :: file
      character(len=*),  intent(in)    :: name
      real(r8), pointer                :: buff1(:,:,:,:)
      real(r8), pointer                :: buff2(:,:,:,:)
      real(r8), pointer                :: buff3(:,:,:,:)
      real(r8), pointer                :: buff4(:,:,:,:)
      integer                          :: numchanid
      type(var_desc_t),  intent(out)   :: vardesc

      integer                          :: lbounds(4,4)
      integer                          :: ubounds(4,4)
      integer                          :: chan_len ! Channel buffer size
      integer                          :: clenid   ! Channel size dimid
      character(len=16)                :: chan_dimname

      ! Make sure we have a dimension for the channel length
      call q3d_chan_bounds(buff1, buff2, buff3, buff4,                        &
           lbounds, ubounds, chan_size=chan_len)
      write(chan_dimname, '("q3d_chan_",i0)') chan_len
      call cam_pio_def_dim(File, trim(chan_dimname), chan_len,                &
           clenid, existOK=.true.)

      call cam_pio_def_var(File, name, PIO_DOUBLE, (/ clenid, numchanid /),   &
           vardesc)

   end subroutine def_var_4d

   subroutine write_var_2d(File, buff1, buff2, buff3, buff4, iodesc, vardesc)
      use pio,                 only: file_desc_t, var_desc_t, io_desc_t
      use pio,                 only: pio_setframe, pio_write_darray
      use pio,                 only: pio_offset_kind

      type(file_desc_t), intent(inout)      :: file
      real(r8),                 pointer     :: buff1(:,:)
      real(r8),                 pointer     :: buff2(:,:)
      real(r8),                 pointer     :: buff3(:,:)
      real(r8),                 pointer     :: buff4(:,:)
      type(io_desc_t),   intent(inout)      :: iodesc
      type(var_desc_t),  intent(inout)      :: vardesc

      ! local variables
      integer                               :: ierr
      integer                               :: lbounds(4,4)
      integer                               :: ubounds(4,4)
      integer                               :: chan_size ! Channel buffer size
      real(r8),                 allocatable :: tmpfield(:,:)
      character(len=16)                     :: chan_dimname
      integer(pio_offset_kind), parameter   :: t_idx = 1

      ! Make sure we have a dimension for the channel length
      call q3d_chan_bounds(buff1, buff2, buff3, buff4,                        &
           lbounds, ubounds, chan_size=chan_size)
      ! Get a (possibly updated) decomp
      call q3d_crm_decomp(File, lbounds, ubounds, chan_size, iodesc)
      ! Collect the data from the buffers
      call crm_to_buffer(buff1, buff2, buff3, buff4, lbounds, ubounds, &
           chan_size, tmpfield)
      ! Write the data
      call PIO_Setframe(file, vardesc, t_idx)
      call PIO_Write_Darray(file, vardesc, iodesc, tmpfield, ierr)

   end subroutine write_var_2d

   subroutine write_var_3d(File, buff1, buff2, buff3, buff4, iodesc, vardesc)
      use pio,                 only: file_desc_t, var_desc_t, io_desc_t
      use pio,                 only: pio_setframe, pio_write_darray
      use pio,                 only: pio_offset_kind

      type(file_desc_t), intent(inout)      :: file
      real(r8),                 pointer     :: buff1(:,:,:)
      real(r8),                 pointer     :: buff2(:,:,:)
      real(r8),                 pointer     :: buff3(:,:,:)
      real(r8),                 pointer     :: buff4(:,:,:)
      type(io_desc_t),   intent(inout)      :: iodesc
      type(var_desc_t),  intent(inout)      :: vardesc

      ! local variables
      integer                               :: ierr
      integer                               :: lbounds(4,4)
      integer                               :: ubounds(4,4)
      integer                               :: chan_size ! Channel buffer size
      real(r8),                 allocatable :: tmpfield(:,:)
      character(len=16)                     :: chan_dimname
      integer(pio_offset_kind), parameter   :: t_idx = 1

      ! Make sure we have a dimension for the channel length
      call q3d_chan_bounds(buff1, buff2, buff3, buff4,                        &
           lbounds, ubounds, chan_size=chan_size)
      ! Get a (possibly updated) decomp
      call q3d_crm_decomp(File, lbounds, ubounds, chan_size, iodesc)
      ! Collect the data from the buffers
      call crm_to_buffer(buff1, buff2, buff3, buff4, lbounds, ubounds, &
           chan_size, tmpfield)
      ! Write the data
      call PIO_Setframe(file, vardesc, t_idx)
      call PIO_Write_Darray(file, vardesc, iodesc, tmpfield, ierr)

   end subroutine write_var_3d

   subroutine write_var_4d(File, buff1, buff2, buff3, buff4, iodesc, vardesc)
      use pio,                 only: file_desc_t, var_desc_t, io_desc_t
      use pio,                 only: pio_setframe, pio_write_darray
      use pio,                 only: pio_offset_kind

      type(file_desc_t), intent(inout)      :: file
      real(r8),                 pointer     :: buff1(:,:,:,:)
      real(r8),                 pointer     :: buff2(:,:,:,:)
      real(r8),                 pointer     :: buff3(:,:,:,:)
      real(r8),                 pointer     :: buff4(:,:,:,:)
      type(io_desc_t),   intent(inout)      :: iodesc
      type(var_desc_t),  intent(inout)      :: vardesc

      ! local variables
      integer                               :: ierr
      integer                               :: lbounds(4,4)
      integer                               :: ubounds(4,4)
      integer                               :: chan_size ! Channel buffer size
      real(r8),                 allocatable :: tmpfield(:,:)
      character(len=16)                     :: chan_dimname
      integer(pio_offset_kind), parameter   :: t_idx = 1

      ! Make sure we have a dimension for the channel length
      call q3d_chan_bounds(buff1, buff2, buff3, buff4,                        &
           lbounds, ubounds, chan_size=chan_size)
      ! Get a (possibly updated) decomp
      call q3d_crm_decomp(File, lbounds, ubounds, chan_size, iodesc)
      ! Collect the data from the buffers
      call crm_to_buffer(buff1, buff2, buff3, buff4, lbounds, ubounds, &
           chan_size, tmpfield)
      ! Write the data
      call pio_setframe(file, vardesc, t_idx)
      call pio_write_darray(file, vardesc, iodesc, tmpfield, ierr)

   end subroutine write_var_4d

   subroutine read_var_2d(File, name, buff1, buff2, buff3, buff4, iodesc)
      use cam_pio_utils,       only: cam_pio_handle_error
      use pio,                 only: file_desc_t, var_desc_t, io_desc_t
      use pio,                 only: pio_setframe, pio_read_darray
      use pio,                 only: pio_offset_kind, pio_inq_varid

      type(file_desc_t), intent(inout)      :: file
      character(len=*),  intent(in)         :: name
      real(r8),                 pointer     :: buff1(:,:)
      real(r8),                 pointer     :: buff2(:,:)
      real(r8),                 pointer     :: buff3(:,:)
      real(r8),                 pointer     :: buff4(:,:)
      type(io_desc_t),   intent(inout)      :: iodesc

      ! local variables
      integer                               :: ierr
      integer                               :: lbounds(4,4)
      integer                               :: ubounds(4,4)
      integer                               :: chan_size ! Channel buffer size
      real(r8),                 allocatable :: tmpfield(:,:)
      character(len=16)                     :: chan_dimname
      type(var_desc_t)                      :: vardesc
      integer(pio_offset_kind), parameter   :: t_idx = 1
      character(len=*),         parameter   :: sub = 'read_var_2d'

      ! Make sure we have a dimension for the channel length
      call q3d_chan_bounds(buff1, buff2, buff3, buff4,                        &
           lbounds, ubounds, chan_size=chan_size)
      ! Get a (possibly updated) decomp
      call q3d_crm_decomp(File, lbounds, ubounds, chan_size, iodesc)
      ! Find the variable in the file
      ierr = pio_inq_varid(File, trim(name), vardesc)
      call cam_pio_handle_error(ierr, sub//': cannot find '//trim(name))
      ! Read the data
      call pio_setframe(file, vardesc, t_idx)
      call pio_read_darray(file, vardesc, iodesc, tmpfield, ierr)
      call cam_pio_handle_error(ierr, sub//': reading '//trim(name))
      ! Disburse the data to the CRM segment buffers
      call buffer_to_crm(tmpfield, lbounds, ubounds, buff1, buff2, buff3, buff4)

   end subroutine read_var_2d

   subroutine read_var_3d(File, name, buff1, buff2, buff3, buff4, iodesc)
      use cam_pio_utils,       only: cam_pio_handle_error
      use pio,                 only: file_desc_t, var_desc_t, io_desc_t
      use pio,                 only: pio_setframe, pio_read_darray
      use pio,                 only: pio_offset_kind, pio_inq_varid

      type(file_desc_t), intent(inout)      :: file
      character(len=*),  intent(in)         :: name
      real(r8),                 pointer     :: buff1(:,:,:)
      real(r8),                 pointer     :: buff2(:,:,:)
      real(r8),                 pointer     :: buff3(:,:,:)
      real(r8),                 pointer     :: buff4(:,:,:)
      type(io_desc_t),   intent(inout)      :: iodesc

      ! local variables
      integer                               :: ierr
      integer                               :: lbounds(4,4)
      integer                               :: ubounds(4,4)
      integer                               :: chan_size ! Channel buffer size
      real(r8),                 allocatable :: tmpfield(:,:)
      character(len=16)                     :: chan_dimname
      type(var_desc_t)                      :: vardesc
      integer(pio_offset_kind), parameter   :: t_idx = 1
      character(len=*),         parameter   :: sub = 'read_var_3d'

      ! Make sure we have a dimension for the channel length
      call q3d_chan_bounds(buff1, buff2, buff3, buff4,                        &
           lbounds, ubounds, chan_size=chan_size)
      ! Get a (possibly updated) decomp
      call q3d_crm_decomp(File, lbounds, ubounds, chan_size, iodesc)
      ! Find the variable in the file
      ierr = pio_inq_varid(File, trim(name), vardesc)
      call cam_pio_handle_error(ierr, sub//': cannot find '//trim(name))
      ! Read the data
      call pio_setframe(file, vardesc, t_idx)
      call pio_read_darray(file, vardesc, iodesc, tmpfield, ierr)
      call cam_pio_handle_error(ierr, sub//': reading '//trim(name))
      ! Disburse the data to the CRM segment buffers
      call buffer_to_crm(tmpfield, lbounds, ubounds, buff1, buff2, buff3, buff4)

   end subroutine read_var_3d

   subroutine read_var_4d(File, name, buff1, buff2, buff3, buff4, iodesc)
      use cam_pio_utils,       only: cam_pio_handle_error
      use pio,                 only: file_desc_t, var_desc_t, io_desc_t
      use pio,                 only: pio_setframe, pio_read_darray
      use pio,                 only: pio_offset_kind, pio_inq_varid

      type(file_desc_t), intent(inout)      :: file
      character(len=*),  intent(in)         :: name
      real(r8),                 pointer     :: buff1(:,:,:,:)
      real(r8),                 pointer     :: buff2(:,:,:,:)
      real(r8),                 pointer     :: buff3(:,:,:,:)
      real(r8),                 pointer     :: buff4(:,:,:,:)
      type(io_desc_t),   intent(inout)      :: iodesc

      ! local variables
      integer                               :: ierr
      integer                               :: lbounds(4,4)
      integer                               :: ubounds(4,4)
      integer                               :: chan_size ! Channel buffer size
      real(r8),                 allocatable :: tmpfield(:,:)
      character(len=16)                     :: chan_dimname
      type(var_desc_t)                      :: vardesc
      integer(pio_offset_kind), parameter   :: t_idx = 1
      character(len=*),         parameter   :: sub = 'read_var_4d'

      ! Make sure we have a dimension for the channel length
      call q3d_chan_bounds(buff1, buff2, buff3, buff4,                        &
           lbounds, ubounds, chan_size=chan_size)
      ! Get a (possibly updated) decomp
      call q3d_crm_decomp(File, lbounds, ubounds, chan_size, iodesc)
      ! Find the variable in the file
      ierr = pio_inq_varid(File, trim(name), vardesc)
      call cam_pio_handle_error(ierr, sub//': cannot find '//trim(name))
      ! Read the data
      call pio_setframe(file, vardesc, t_idx)
      call pio_read_darray(file, vardesc, iodesc, tmpfield, ierr)
      call cam_pio_handle_error(ierr, sub//': reading '//trim(name))
      ! Disburse the data to the CRM segment buffers
      call buffer_to_crm(tmpfield, lbounds, ubounds, buff1, buff2, buff3, buff4)

   end subroutine read_var_4d

   subroutine q3d_chan_bounds_2d(buffer1, buffer2, buffer3, buffer4,          &
        lbounds, ubounds, chan_size, buffer_rank)
      use spmd_utils, only: iam, mpicom, MPI_INTEGER

      real(r8), pointer               :: buffer1(:,:)
      real(r8), pointer               :: buffer2(:,:)
      real(r8), pointer               :: buffer3(:,:)
      real(r8), pointer               :: buffer4(:,:)
      integer,            intent(out) :: lbounds(4,4)
      integer,            intent(out) :: ubounds(4,4)
      integer,  optional, intent(out) :: chan_size
      integer,  optional, intent(out) :: buffer_rank

      integer                         :: seg, dim
      integer                         :: ierr
      integer                         :: bsize

      ! First, get buffer bounds to find out if we need to compute a new iodesc
      lbounds = 0
      ubounds = -1
      if (associated(buffer1)) then
         lbounds(1:2,1) = LBOUND(buffer1)
         ubounds(1:2,1) = UBOUND(buffer1)
      end if
      if (associated(buffer2)) then
         lbounds(1:2,2) = LBOUND(buffer2)
         ubounds(1:2,2) = UBOUND(buffer2)
      end if
      if (associated(buffer3)) then
         lbounds(1:2,3) = LBOUND(buffer3)
         ubounds(1:2,3) = UBOUND(buffer3)
      end if
      if (associated(buffer4)) then
         lbounds(1:2,4) = LBOUND(buffer4)
         ubounds(1:2,4) = UBOUND(buffer4)
      end if

      if (present(chan_size)) then
         ! We need a consistent size so do on root and broadcast
         if (iam == 0) then
            chan_size = 0
            do seg = 1, 4
               bsize = 1
               do dim = 1, 2
                  bsize = bsize * (ubounds(dim, seg) - lbounds(dim, seg) + 1)
               end do
               chan_size = chan_size + bsize
            end do
         end if
      end if
      call MPI_bcast(chan_size, 1, MPI_INTEGER, 0, mpicom, ierr)
      if (present(buffer_rank)) then
         buffer_rank = 2
      end if

   end subroutine q3d_chan_bounds_2d

   subroutine q3d_chan_bounds_3d(buffer1, buffer2, buffer3, buffer4,          &
        lbounds, ubounds, chan_size, buffer_rank)
      use spmd_utils, only: iam, mpicom, MPI_INTEGER

      real(r8), pointer               :: buffer1(:,:,:)
      real(r8), pointer               :: buffer2(:,:,:)
      real(r8), pointer               :: buffer3(:,:,:)
      real(r8), pointer               :: buffer4(:,:,:)
      integer,            intent(out) :: lbounds(4,4)
      integer,            intent(out) :: ubounds(4,4)
      integer,  optional, intent(out) :: chan_size
      integer,  optional, intent(out) :: buffer_rank

      integer                         :: seg, dim
      integer                         :: ierr
      integer                         :: bsize

      ! First, get buffer bounds to find out if we need to compute a new iodesc
      lbounds = 0
      ubounds = -1
      if (associated(buffer1)) then
         lbounds(1:3,1) = LBOUND(buffer1)
         ubounds(1:3,1) = UBOUND(buffer1)
      end if
      if (associated(buffer2)) then
         lbounds(1:3,2) = LBOUND(buffer2)
         ubounds(1:3,2) = UBOUND(buffer2)
      end if
      if (associated(buffer3)) then
         lbounds(1:3,3) = LBOUND(buffer3)
         ubounds(1:3,3) = UBOUND(buffer3)
      end if
      if (associated(buffer4)) then
         lbounds(1:3,4) = LBOUND(buffer4)
         ubounds(1:3,4) = UBOUND(buffer4)
      end if

      if (present(chan_size)) then
         ! We need a consistent size so do on root and broadcast
         if (iam == 0) then
            chan_size = 0
            do seg = 1, 4
               bsize = 1
               do dim = 1, 3
                  bsize = bsize * (ubounds(dim, seg) - lbounds(dim, seg) + 1)
               end do
               chan_size = chan_size + bsize
            end do
         end if
      end if
      call MPI_bcast(chan_size, 1, MPI_INTEGER, 0, mpicom, ierr)
      if (present(buffer_rank)) then
         buffer_rank = 3
      end if

   end subroutine q3d_chan_bounds_3d

   subroutine q3d_chan_bounds_4d(buffer1, buffer2, buffer3, buffer4,          &
        lbounds, ubounds, chan_size, buffer_rank)
      use spmd_utils, only: iam, mpicom, MPI_INTEGER

      real(r8), pointer               :: buffer1(:,:,:,:)
      real(r8), pointer               :: buffer2(:,:,:,:)
      real(r8), pointer               :: buffer3(:,:,:,:)
      real(r8), pointer               :: buffer4(:,:,:,:)
      integer,            intent(out) :: lbounds(4,4)
      integer,            intent(out) :: ubounds(4,4)
      integer,  optional, intent(out) :: chan_size
      integer,  optional, intent(out) :: buffer_rank

      integer                         :: seg, dim
      integer                         :: ierr
      integer                         :: bsize

      ! First, get buffer bounds to find out if we need to compute a new iodesc
      lbounds = 0
      ubounds = -1
      if (associated(buffer1)) then
         lbounds(1:4,1) = LBOUND(buffer1)
         ubounds(1:4,1) = UBOUND(buffer1)
      end if
      if (associated(buffer2)) then
         lbounds(1:4,2) = LBOUND(buffer2)
         ubounds(1:4,2) = UBOUND(buffer2)
      end if
      if (associated(buffer3)) then
         lbounds(1:4,3) = LBOUND(buffer3)
         ubounds(1:4,3) = UBOUND(buffer3)
      end if
      if (associated(buffer4)) then
         lbounds(1:4,4) = LBOUND(buffer4)
         ubounds(1:4,4) = UBOUND(buffer4)
      end if

      if (present(chan_size)) then
         ! We need a consistent size so do on root and broadcast
         if (iam == 0) then
            chan_size = 0
            do seg = 1, 4
               bsize = 1
               do dim = 1, 4
                  bsize = bsize * (ubounds(dim, seg) - lbounds(dim, seg) + 1)
               end do
               chan_size = chan_size + bsize
            end do
         end if
      end if
      call MPI_bcast(chan_size, 1, MPI_INTEGER, 0, mpicom, ierr)
      if (present(buffer_rank)) then
         buffer_rank = 4
      end if

   end subroutine q3d_chan_bounds_4d

   subroutine q3d_crm_decomp(File, lbounds, ubounds, chan_size, iodesc)
      use cam_pio_utils,  only: pio_subsystem
      use pio,            only: io_desc_t, file_desc_t, PIO_DOUBLE
      use pio,            only: pio_initdecomp
      use cam_abortutils, only: endrun
      use cam_map_utils,  only: iMap
      use q3d_runtime,    only: q3d_begchan, q3d_endchan, q3d_total_channels

      type(file_desc_t), intent(inout) :: File
      integer,           intent(in)    :: lbounds(4,4)
      integer,           intent(in)    :: ubounds(4,4)
      integer,           intent(in)    :: chan_size
      type(io_desc_t),   intent(inout) :: iodesc

      integer                          :: i, j, k, m, seg, lind, chan
      integer                          :: num_chans
      integer(iMap)                    :: lpos
      integer(iMap),   allocatable     :: ldof(:)

      if (need_decomp) then
         num_chans = q3d_endchan - q3d_begchan + 1
         allocate(ldof(chan_size * num_chans))
         lind = 0 ! The index for ldof
         lpos = chan_size * (q3d_begchan - 1) - 1
         do chan = q3d_begchan, q3d_endchan
            do seg = 1, 4
               do m = lbounds(4, seg), max(ubounds(4, seg),lbounds(4, seg))
                  do k = lbounds(3, seg), max(ubounds(3, seg),lbounds(3, seg))
                     do j = lbounds(2, seg), ubounds(2, seg)
                        do i = lbounds(1, seg), ubounds(1, seg)
                           lind = lind + 1
                           lpos = lpos + 1
                           ldof(lind) = lpos
                        end do
                     end do
                  end do
               end do
            end do
         end do
         call pio_initdecomp(pio_subsystem, pio_double,                      &
              (/ chan_size, q3d_total_channels /), ldof, iodesc)
         deallocate(ldof)
         have_decomp = .true.
         ldof_lbound = lbounds
         ldof_ubound = ubounds
         need_decomp = .false.
      else if (ANY(lbounds /= ldof_lbound) .or. ANY(ubounds /= ldof_ubound)) then
         call endrun("q3d_crm_decomp: bad bounds")
      end if

   end subroutine q3d_crm_decomp

   subroutine q3d_free_decomp(File, iodesc)
      use pio,            only: file_desc_t, io_desc_t
      use pio,            only: pio_freedecomp

      type(file_desc_t), intent(inout) :: File
      type(io_desc_t),   intent(inout) :: iodesc

      if (have_decomp) then
         call pio_freedecomp(File, iodesc)
         have_decomp = .false.
      end if
      need_decomp = .true.

   end subroutine q3d_free_decomp

   subroutine alloc_crm(chan_size, num_chans, buffer)
      integer,               intent(in)    :: chan_size
      integer,               intent(in)    :: num_chans
      real(r8), allocatable, intent(inout) :: buffer(:,:)

      if (allocated(buffer)) then
         if ((size(buffer,1) /= chan_size) .or. (size(buffer,2) /= num_chans)) then
            deallocate(buffer)
         end if
      end if
      if (.not. allocated(buffer)) then
         allocate(buffer(chan_size, num_chans))
      end if
   end subroutine alloc_crm

   subroutine crm_to_buffer_2d(crm_f1, crm_f2, crm_f3, crm_f4, lbounds, ubounds, chan_size, buffer)
      use cam_abortutils, only: endrun
      use q3d_runtime,    only: q3d_begchan, q3d_endchan

      real(r8), pointer                  :: crm_f1(:,:)
      real(r8), pointer                  :: crm_f2(:,:)
      real(r8), pointer                  :: crm_f3(:,:)
      real(r8), pointer                  :: crm_f4(:,:)
      integer,               intent(in)  :: lbounds(4,4)
      integer,               intent(in)  :: ubounds(4,4)
      integer,               intent(in)  :: chan_size
      real(r8), allocatable, intent(out) :: buffer(:,:)

      integer                                    :: i, j, chan, seg, bind
      character(len=*), parameter                :: subname = 'crm_to_buffer_2d'

      if (ANY(ubounds(1:2, :) < lbounds(1:2,:))) then
         call alloc_crm(0, q3d_endchan - q3d_begchan + 1, buffer)
      else
         call alloc_crm(chan_size, q3d_endchan - q3d_begchan + 1, buffer)
      end if
      do chan = q3d_begchan, q3d_endchan
         bind = 0
         do seg = 1, 4
            do j = lbounds(2, seg), ubounds(2, seg)
               do i = lbounds(1, seg), ubounds(1, seg)
                  bind = bind + 1
                  select case(seg)
                  case(1)
                     buffer(bind, chan - q3d_begchan + 1) = crm_f1(i, j)
                  case(2)
                     buffer(bind, chan - q3d_begchan + 1) = crm_f2(i, j)
                  case(3)
                     buffer(bind, chan - q3d_begchan + 1) = crm_f3(i, j)
                  case(4)
                     buffer(bind, chan - q3d_begchan + 1) = crm_f4(i, j)
                  end select
               end do
            end do
         end do
      end do
   end subroutine crm_to_buffer_2d

   subroutine crm_to_buffer_3d(crm_f1, crm_f2, crm_f3, crm_f4, lbounds, ubounds, chan_size, buffer)
      use cam_abortutils, only: endrun
      use q3d_runtime,    only: q3d_begchan, q3d_endchan

      real(r8), pointer                  :: crm_f1(:,:,:)
      real(r8), pointer                  :: crm_f2(:,:,:)
      real(r8), pointer                  :: crm_f3(:,:,:)
      real(r8), pointer                  :: crm_f4(:,:,:)
      integer,               intent(in)  :: lbounds(4,4)
      integer,               intent(in)  :: ubounds(4,4)
      integer,               intent(in)  :: chan_size
      real(r8), allocatable, intent(out) :: buffer(:,:)

      integer                                    :: i, j, k, chan, seg, bind
      character(len=*), parameter                :: subname = 'crm_to_buffer_3d'

      if (ANY(ubounds(1:3, :) < lbounds(1:3,:))) then
         call alloc_crm(0, q3d_endchan - q3d_begchan + 1, buffer)
      else
         call alloc_crm(chan_size, q3d_endchan - q3d_begchan + 1, buffer)
      end if
      do chan = q3d_begchan, q3d_endchan
         bind = 0
         do seg = 1, 4
            do k = lbounds(3, seg), ubounds(3, seg)
               do j = lbounds(2, seg), ubounds(2, seg)
                  do i = lbounds(1, seg), ubounds(1, seg)
                     bind = bind + 1
                     select case(seg)
                     case(1)
                        buffer(bind, chan - q3d_begchan + 1) = crm_f1(i, j, k)
                     case(2)
                        buffer(bind, chan - q3d_begchan + 1) = crm_f2(i, j, k)
                     case(3)
                        buffer(bind, chan - q3d_begchan + 1) = crm_f3(i, j, k)
                     case(4)
                        buffer(bind, chan - q3d_begchan + 1) = crm_f4(i, j, k)
                     end select
                  end do
               end do
            end do
         end do
      end do
   end subroutine crm_to_buffer_3d

   subroutine crm_to_buffer_4d(crm_f1, crm_f2, crm_f3, crm_f4, lbounds, ubounds, chan_size, buffer)
      use cam_abortutils, only: endrun
      use q3d_runtime,    only: q3d_begchan, q3d_endchan

      real(r8), pointer                  :: crm_f1(:,:,:,:)
      real(r8), pointer                  :: crm_f2(:,:,:,:)
      real(r8), pointer                  :: crm_f3(:,:,:,:)
      real(r8), pointer                  :: crm_f4(:,:,:,:)
      integer,               intent(in)  :: lbounds(4,4)
      integer,               intent(in)  :: ubounds(4,4)
      integer,               intent(in)  :: chan_size
      real(r8), allocatable, intent(out) :: buffer(:,:)

      integer                                    :: i, j, k, m, chan, seg, bind
      character(len=*), parameter                :: subname = 'crm_to_buffer_4d'

      ! Sanity check on crm fields
      if (ANY(ubounds(1:4, :) < lbounds(1:4,:))) then
         call alloc_crm(0, q3d_endchan - q3d_begchan + 1, buffer)
      else
         call alloc_crm(chan_size, q3d_endchan - q3d_begchan + 1, buffer)
      end if
      do chan = q3d_begchan, q3d_endchan
         bind = 0
         do seg = 1, 4
            do m = lbounds(4, seg), ubounds(4, seg)
               do k = lbounds(3, seg), ubounds(3, seg)
                  do j = lbounds(2, seg), ubounds(2, seg)
                     do i = lbounds(1, seg), ubounds(1, seg)
                        bind = bind + 1
                        select case(seg)
                        case(1)
                           buffer(bind, chan - q3d_begchan + 1) = crm_f1(i, j, k, m)
                        case(2)
                           buffer(bind, chan - q3d_begchan + 1) = crm_f2(i, j, k, m)
                        case(3)
                           buffer(bind, chan - q3d_begchan + 1) = crm_f3(i, j, k, m)
                        case(4)
                           buffer(bind, chan - q3d_begchan + 1) = crm_f4(i, j, k, m)
                        end select
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine crm_to_buffer_4d

   subroutine buffer_to_crm_2d(buffer, lbounds, ubounds, crm_f1, crm_f2, crm_f3, crm_f4)
      use cam_abortutils, only: endrun
      use q3d_runtime,    only: q3d_begchan, q3d_endchan

      real(r8), intent(in)        :: buffer(:,:)
      integer,  intent(in)        :: lbounds(4,4)
      integer,  intent(in)        :: ubounds(4,4)
      real(r8),         pointer   :: crm_f1(:,:)
      real(r8),         pointer   :: crm_f2(:,:)
      real(r8),         pointer   :: crm_f3(:,:)
      real(r8),         pointer   :: crm_f4(:,:)

      integer                     :: i, j, chan, seg, bind
      character(len=*), parameter :: subname = 'buffer_to_crm_2d'

      ! Sanity check on crm fields
      if (ANY(ubounds(1:2, :) < lbounds(1:2,:))) then
         call endrun(subname//': crm bounds must be rank two')
      end if
      if (ANY(ubounds(3:4, :) >= lbounds(3:4,:))) then
         call endrun(subname//': crm bounds must be rank two')
      end if
      do chan = q3d_begchan, q3d_endchan
         bind = 0
         do seg = 1, 4
            do j = lbounds(2, seg), ubounds(2, seg)
               do i = lbounds(1, seg), ubounds(1, seg)
                  bind = bind + 1
                  select case(seg)
                  case(1)
                     crm_f1(i, j) = buffer(bind, chan - q3d_begchan + 1)
                  case(2)
                     crm_f2(i, j) = buffer(bind, chan - q3d_begchan + 1)
                  case(3)
                     crm_f3(i, j) = buffer(bind, chan - q3d_begchan + 1)
                  case(4)
                     crm_f4(i, j) = buffer(bind, chan - q3d_begchan + 1)
                  end select
               end do
            end do
         end do
      end do
   end subroutine buffer_to_crm_2d

   subroutine buffer_to_crm_3d(buffer, lbounds, ubounds, crm_f1, crm_f2, crm_f3, crm_f4)
      use cam_abortutils, only: endrun
      use q3d_runtime,    only: q3d_begchan, q3d_endchan

      real(r8), intent(in)        :: buffer(:,:)
      integer,  intent(in)        :: lbounds(4,4)
      integer,  intent(in)        :: ubounds(4,4)
      real(r8),         pointer   :: crm_f1(:,:,:)
      real(r8),         pointer   :: crm_f2(:,:,:)
      real(r8),         pointer   :: crm_f3(:,:,:)
      real(r8),         pointer   :: crm_f4(:,:,:)

      integer                     :: i, j, k, chan, seg, bind
      character(len=*), parameter :: subname = 'buffer_to_crm_3d'

      ! Sanity check on crm fields
      if (ANY(ubounds(1:3, :) < lbounds(1:3,:))) then
         call endrun(subname//': crm bounds must be rank three')
      end if
      if (ANY(ubounds(4, :) >= lbounds(4,:))) then
         call endrun(subname//': crm bounds must be rank three')
      end if
      do chan = q3d_begchan, q3d_endchan
         bind = 0
         do seg = 1, 4
            do k = lbounds(3, seg), ubounds(3, seg)
               do j = lbounds(2, seg), ubounds(2, seg)
                  do i = lbounds(1, seg), ubounds(1, seg)
                     bind = bind + 1
                     select case(seg)
                     case(1)
                        crm_f1(i, j, k) = buffer(bind, chan - q3d_begchan + 1)
                     case(2)
                        crm_f2(i, j, k) = buffer(bind, chan - q3d_begchan + 1)
                     case(3)
                        crm_f3(i, j, k) = buffer(bind, chan - q3d_begchan + 1)
                     case(4)
                        crm_f4(i, j, k) = buffer(bind, chan - q3d_begchan + 1)
                     end select
                  end do
               end do
            end do
         end do
      end do
   end subroutine buffer_to_crm_3d

   subroutine buffer_to_crm_4d(buffer, lbounds, ubounds, crm_f1, crm_f2, crm_f3, crm_f4)
      use cam_abortutils, only: endrun
      use q3d_runtime,    only: q3d_begchan, q3d_endchan

      real(r8), intent(in)        :: buffer(:,:)
      integer,  intent(in)        :: lbounds(4,4)
      integer,  intent(in)        :: ubounds(4,4)
      real(r8),         pointer   :: crm_f1(:,:,:,:)
      real(r8),         pointer   :: crm_f2(:,:,:,:)
      real(r8),         pointer   :: crm_f3(:,:,:,:)
      real(r8),         pointer   :: crm_f4(:,:,:,:)

      integer                     :: i, j, k, m, chan, seg, bind
      character(len=*), parameter :: subname = 'buffer_to_crm_4d'

      ! Sanity check on crm fields
      if (ANY(ubounds(1:4, :) < lbounds(1:4,:))) then
         call endrun(subname//': crm bounds must be rank four')
      end if
      do chan = q3d_begchan, q3d_endchan
         bind = 0
         do seg = 1, 4
            do m = lbounds(4, seg), ubounds(4, seg)
               do k = lbounds(3, seg), ubounds(3, seg)
                  do j = lbounds(2, seg), ubounds(2, seg)
                     do i = lbounds(1, seg), ubounds(1, seg)
                        bind = bind + 1
                        select case(seg)
                        case(1)
                           crm_f1(i, j, k, m) = buffer(bind, chan - q3d_begchan + 1)
                        case(2)
                           crm_f2(i, j, k, m) = buffer(bind, chan - q3d_begchan + 1)
                        case(3)
                           crm_f3(i, j, k, m) = buffer(bind, chan - q3d_begchan + 1)
                        case(4)
                           crm_f4(i, j, k, m) = buffer(bind, chan - q3d_begchan + 1)
                        end select
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine buffer_to_crm_4d

end module q3d_restart
