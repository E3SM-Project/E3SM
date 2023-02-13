module conditional_diag_output_utils
!----------------------------------------------------------------------
! Utility subroutines used for handling model output 
! for the conditional sampling and budget analysis tool. 
! Currently included are a few small subroutines that construct output
! variable names, and a subroutine that makes addfld and add_default 
! calls to register output variables during model initialization.
! The outfld calls are in module conditional_diag_main.
!
! History:
!  First version by Hui Wan, PNNL, March-May 2021
!----------------------------------------------------------------------
  use cam_abortutils, only: endrun
  use conditional_diag, only: cnd_diag_info_t

  implicit none
  public

contains

subroutine cnd_diag_output_init(pver, cnd_diag_info)
!----------------------------------------------------------------------------------- 
! 
! Purpose: Register variables related to conditional sampling and budget analysis
!          for history output
!
! Method: (1) Add variables to the master field list by doing addfld calls
!         (2) Add the whole suite of supported output variables to 
!             user-specified history tapes by doing add_default calls.
!         These two things are done for each sampling condition. 
!         Registered output variables include, for each condition,
!          - the metric used for evaluating the user-specified sampling condition,
!          - the flag field resulting from evaluating the sampling condition,
!          - the various QoIs (and their increments if requested) to which 
!            conditional sampling is applied. Per user's choice, 
!            these QoIs and increments might be monitored at various checkpoints
!            in each time step and might be multiplied with the pressure layer
!            thickness.
!-----------------------------------------------------------------------------------
  use cam_history_support, only: max_fieldname_len, horiz_only
  use cam_history,         only: addfld, add_default
  use conditional_diag,    only: FILLVALUE

  integer,intent(in) :: pver
  type(cnd_diag_info_t), intent(in) :: cnd_diag_info

  integer          :: ntape,itape
  integer          :: icnd, iqoi, ichkpt, ii
  character(len=4) :: val_inc_suff(2), suff
  logical          :: l_output(2)

  character(len=max_fieldname_len) :: output_fld_name
  character(len=max_fieldname_len) :: output_fld_name2

  character(len=256) :: fld_long_name

  character(len=*),parameter :: subname = 'cnd_diag_output_init'

  if (cnd_diag_info%ncnd==0) return

  ntape = cnd_diag_info%ntape

  ! Loop through all sampling conditions. Each of them will have
  ! its own set of output variables distinguished by the prefix cndxx_
  ! where xx is a two-digit integer.

  do icnd = 1,cnd_diag_info%ncnd

     !------------------------------------------------------------------------------
     ! Register the metric of the sampling condition and the associated flag field
     !------------------------------------------------------------------------------
     call get_metric_and_flag_names_for_output( icnd, cnd_diag_info, &!in
                                                output_fld_name,     &!out
                                                output_fld_name2    ) !out

     ! Add the 2 variables to the master list of possible output variables.
     ! The arguments of an addfld call are: (1) variable name, 
     ! (2) vertical dimension name, (3) avg. flag, (4) unit, (5) long_name.
     ! Units are set to blank right now; we could add a namelist variable
     ! to let the user provide the info. 

     if (cnd_diag_info%metric_nver(icnd)==1) then

       call addfld(trim(output_fld_name ), horiz_only, 'A',' ','Metric used in conditional sampling')
       call addfld(trim(output_fld_name2), horiz_only, 'A',' ','Flags  used in conditional sampling')

     elseif(cnd_diag_info%metric_nver(icnd)==pver) then

       call addfld(trim(output_fld_name ), (/'lev'/),  'A',' ','Metric used in conditional sampling')
       call addfld(trim(output_fld_name2), (/'lev'/),  'A',' ','Flags  used in conditional sampling')

     elseif(cnd_diag_info%metric_nver(icnd)==pver+1) then

       call addfld(trim(output_fld_name ), (/'ilev'/), 'A',' ','Metric used in conditional sampling')
       call addfld(trim(output_fld_name2), (/'ilev'/), 'A',' ','Flags  used in conditional sampling')

     else 
       call endrun(subname//': invalid number of vertical levels for metric '//trim(cnd_diag_info%metric_name(icnd)))
     end if 

     ! Add the variable to user-specified history tapes.
     ! The 3 arguments of an add_default call are (1) variable name, (2) hist. tape index,
     ! and (3) hist. averaging flag

     do itape = 1,ntape
        call add_default(trim(output_fld_name ), cnd_diag_info%hist_tape_with_all_output(itape), ' ')  
        call add_default(trim(output_fld_name2), cnd_diag_info%hist_tape_with_all_output(itape), ' ')
     end do

     !-------------------------------------------------------------------
     ! Register the diagnostics fields and their increments 
     !-------------------------------------------------------------------
     ! Put the on/off switches for value and increment output into one
     ! array so that we can deal with both using the ii-loop below

     l_output = (/cnd_diag_info%l_output_state, cnd_diag_info%l_output_incrm/)

     ! In terms of variable names in the output files, the values and 
     ! increments of a field are distinguished by different suffixes

     val_inc_suff = (/"","_inc"/)

     do ii = 1,2 ! field value (state) or increment

        if (.not.l_output(ii)) cycle
        suff = val_inc_suff(ii)

        do iqoi  = 1,cnd_diag_info%nqoi
        do ichkpt = 1,cnd_diag_info%nchkpt

           call get_fld_name_for_output( suff, cnd_diag_info, &! in
                                         icnd, iqoi, ichkpt,  &! in
                                         output_fld_name      )! out

           call get_fld_longname_for_output( suff, cnd_diag_info, &! in
                                             icnd, iqoi, ichkpt,  &! in
                                             fld_long_name        )! out

           ! Add the variable to the master list of possible output variables.
           ! The arguments of an addfld call are: (1) variable name, 
           ! (2) vertical dimension name, (3) avg. flag, (4) unit, (5) long_name.
           ! Units are set to blank right now; we could add a namelist variable
           ! to let the user provide the info. 

           if (cnd_diag_info%qoi_nver_save(iqoi)==1) then
              call addfld(trim(output_fld_name), horiz_only, 'A',' ',trim(fld_long_name)) 

           elseif (cnd_diag_info%qoi_nver_save(iqoi)==pver) then
              call addfld(trim(output_fld_name), (/'lev'/),  'A',' ',trim(fld_long_name)) 

           elseif (cnd_diag_info%qoi_nver_save(iqoi)==pver+1) then
              call addfld(trim(output_fld_name), (/'ilev'/), 'A',' ',trim(fld_long_name)) 
           else
              call endrun(subname//': invalid number of vertical levels for '//cnd_diag_info%qoi_name(iqoi))
           end if

           ! Add the variable to user-specified history tapes.
           ! The 3 arguments of an add_default call are (1) variable name, (2) hist. tape index,
           ! and (3) hist. averaging flag

           do itape = 1,ntape
              call add_default(trim(output_fld_name ), cnd_diag_info%hist_tape_with_all_output(itape), ' ')  
           end do

        end do ! ichkpt
        end do ! iqoi

     end do ! ii = 1,2, field value (state) or tendency

  end do ! icnd = 1,ncnd

end subroutine cnd_diag_output_init

!======================================================
subroutine get_metric_and_flag_names_for_output( icnd, cnd_diag_info,   &! in
                                               metric_name_in_output, &! out
                                               flag_name_in_output    )! out

   use cam_history_support, only: max_fieldname_len

   integer,               intent(in)  :: icnd
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info

   character(len=max_fieldname_len) :: metric_name_in_output
   character(len=max_fieldname_len) :: flag_name_in_output

   character(len=2) :: icnd_str ! condition index as a string

   write(icnd_str,'(i2.2)') icnd
   metric_name_in_output = 'cnd'//icnd_str//'_'//trim(cnd_diag_info%metric_name(icnd))
     flag_name_in_output = 'cnd'//icnd_str//'_'//trim(cnd_diag_info%metric_name(icnd))//'_flag'

end subroutine get_metric_and_flag_names_for_output

!======================================================
subroutine get_fld_name_for_output( suff, cnd_diag_info,    &!in
                                    icnd, iqoi, ichkpt,     &!in
                                    fld_name_in_output      )!out

   use cam_history_support, only: max_fieldname_len
   use conditional_diag,    only: NODP

   character(len=*),      intent(in)  :: suff
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info
   integer,               intent(in)  :: icnd, iqoi, ichkpt

   character(len=max_fieldname_len),intent(out) :: fld_name_in_output 

   character(len=2) :: icnd_str   ! condition index as a string
   character(len=2) :: vint       ! suffix

   ! If the field will be vertically integrated, append "_v" to the QoI name

   if (cnd_diag_info%x_dp(iqoi,ichkpt)/=NODP) then
      vint = '_v'
   else
      vint = ''
   end if

   ! Now construct the full name of the variable in output file

   write(icnd_str,'(i2.2)') icnd

   fld_name_in_output = 'cnd'//icnd_str//'_'// &
                        trim(cnd_diag_info%qoi_name(iqoi))//trim(vint)//'_'// &
                        trim(cnd_diag_info%qoi_chkpt(ichkpt))//suff

end subroutine get_fld_name_for_output 

!======================================================
subroutine get_fld_longname_for_output( suff, cnd_diag_info,    &!in
                                        icnd, iqoi, ichkpt,      &!in
                                        fld_long_name_in_output )!out

   use cam_history_support, only: max_fieldname_len
   use conditional_diag,    only: NODP

   character(len=*),      intent(in)  :: suff
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info
   integer,               intent(in)  :: icnd, iqoi, ichkpt

   character(len=256),intent(out)     :: fld_long_name_in_output 

   character(len=2) :: icnd_str ! condition index as a string
   character(len=2) :: vint       ! suffix

   ! If the field will be vertically integrated, append "_v" to the QoI name

   if (cnd_diag_info%x_dp(iqoi,ichkpt)/=NODP) then
      vint = '_v'
   else
      vint = ''
   end if

   write(icnd_str,'(i2.2)') icnd

   fld_long_name_in_output = trim(cnd_diag_info%qoi_name(iqoi))//trim(vint)//trim(suff)// &
                             ' at '//trim(cnd_diag_info%qoi_chkpt(ichkpt))// &
                             ' sampled under condition '//icnd_str// &
                             ' ('//trim(cnd_diag_info%metric_name(icnd))//')' 

end subroutine get_fld_longname_for_output 

end module conditional_diag_output_utils
