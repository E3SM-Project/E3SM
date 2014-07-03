module esmfshr_attribute_mod

   use seq_flds_mod
   use seq_comm_mct
   use ESMF

   implicit none

   private

   public :: esmfshr_attribute_appl_init, esmfshr_attribute_fields_init

!----------------------------------------------------------------------------
   contains
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
! Application attribute initialization
!----------------------------------------------------------------------------

   subroutine esmfshr_attribute_appl_init(comp, rc)

      implicit none

      type(ESMF_CplComp), intent(inout)  :: comp
      integer, intent(out)               :: rc

      integer :: localrc
      character(ESMF_MAXSTR) :: convCIM, purpComp

      convCIM  = 'CIM'
      purpComp = 'Model Component Simulation Description'

      call ESMF_AttributeAdd(comp,  &
                             convention=convCIM, purpose=purpComp, rc=localrc)

      call ESMF_AttributeSet(comp, &
                             'ShortName', 'CESM',  &
                             convention=convCIM, purpose=purpComp, rc=localrc)
      call ESMF_AttributeSet(comp, &
                             'LongName', &
                             'Community Earth System Model', &
                             convention=convCIM, purpose=purpComp, rc=localrc)
      call ESMF_AttributeSet(comp, &
                             'Description',  &
                             "The Community Earth System Model (CESM) is " // &
                             "a coupled climate model for simulating the " // &
                             "earth's climate system. Composed of four " // &
                             "separate models simultaneously simulating " // &
                             "the earth's atmosphere, ocean, land surface " // &
                             "and sea-ice, and one central coupler " // &
                             "component, the CESM allows researchers to " // &
                             "conduct fundamental research into the " // &
                             "earth's past, present and future climate " // &
                             "states.",  &
                             convention=convCIM, purpose=purpComp, rc=localrc)
      call ESMF_AttributeSet(comp, &
                             'ReleaseDate', '2010',  &
                             convention=convCIM, purpose=purpComp, rc=localrc)
      call ESMF_AttributeSet(comp, &
                             'ModelType', 'Earth System',  &
                             convention=convCIM, purpose=purpComp, rc=localrc)
      call ESMF_AttributeSet(comp, &
                             'URL', 'www.cesm.ucar.edu',  &
                             convention=convCIM, purpose=purpComp, rc=localrc)

      rc = localrc

   end subroutine


!----------------------------------------------------------------------------
! Field attribute initialization
!----------------------------------------------------------------------------

   subroutine esmfshr_attribute_fields_init(cesmState, rc)

      implicit none

      type(ESMF_State), intent(inout)  :: cesmState
      integer, intent(out)             :: rc

      ! Local variables
      integer :: localrc

!      call ESMF_LogWrite('a2x fields = >'//trim(seq_flds_a2x_fields)//'<', &
!                         ESMF_LOG_INFO, rc=localrc)

      call esmfshr_attribute_field_init( &
              seq_flds_a2x_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_x2a_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_i2x_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_x2i_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_l2x_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_x2l_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_o2x_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_x2o_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_g2x_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_x2g_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_s2x_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_x2s_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_xao_fields, cesmState, rc=localrc)
      call esmfshr_attribute_field_init( &
              seq_flds_r2x_fields, cesmState, rc=localrc)

      rc = localrc

   end subroutine


   subroutine esmfshr_attribute_field_init(fieldList, state, rc)

      implicit none

      character(len=*), intent(in)     :: fieldList
      type(ESMF_State), intent(inout)  :: state
      integer, intent(out)             :: rc

      ! Local variables
      character(ESMF_MAXSTR)              :: convCIM, purpField
      character(ESMF_MAXSTR),dimension(4) :: attrList

      integer :: i
      integer :: localrc
      integer :: nFields
      integer :: curIndex, preIndex

      character(len=80) :: shortname, longname, stdname, units
      character(len=80) :: str

      type(ESMF_Field), allocatable, dimension(:) :: fields
      type(ESMF_FieldBundle) :: fbundle
      !type(ESMF_LOG) :: alog

      ! Set up the attribute descriptions
      convCIM   = 'CIM'
      purpField = 'Inputs Description'

      !--- make sure field list is not empty ---
      if (trim(fieldList) /= ':') then

         !--- get number of fields in the list ---
         nFields = 0
         do i = 1, len_trim(fieldList)
            if (fieldList(i:i) == ':') nFields = nFields+1
         end do
         nFields = nFields +1

         !--- write log ---
!         call ESMF_LogWrite('Field list = >'//trim(fieldList)//'<', &
!                            ESMF_LOG_INFO, rc=localrc)
!         write(str, fmt='("Field count = ", I5)') nFields
!         call ESMF_LogWrite(str, ESMF_LOG_INFO, rc=localrc)

         !--- allocate ESMF_Field array ---
         if (allocated(fields)) deallocate(fields)
         allocate(fields(nFields))

         !--- set initial values of indices ---
         curIndex = index(fieldList, ':')
         preIndex = 1

         !--- set attribute list ---
         attrList(1) = 'ShortName'
         attrList(2) = 'StandardName'
         attrList(3) = 'LongName'
         attrList(4) = 'Units'
!
! KDS: Do I need to add these attributes for each of the fields, or for just
!      specific fields???
!
         !attrList(5) = 'InputType'
         !attrList(6) = 'InputTargetComponent'
         !attrList(7) = 'InputSpatialRegriddingMethod'
         !attrList(8) = 'InputSpatialRegriddingType'
         !attrList(9) = 'InputFrequency'
         !attrList(10) = 'InputTimeTransformationType'

         ! Go through the list of fields in fieldList
         i = 0
         do while (curIndex > preIndex)

            i = i + 1
            shortname = fieldList(preIndex:curIndex - 1)
!            call ESMF_LogWrite('KDS: shortname => ' //trim(shortname)// ' <', &
!                               ESMF_LOG_INFO, rc=localrc)

            !--- create empty ESMF_Field ---
            fields(i) = ESMF_FieldEmptyCreate(name=trim(shortname), rc=localrc)

            !--- get field metadata ---
            call seq_flds_esmf_metadata_get(shortname, longname, stdname, units)

!            call ESMF_LogWrite(trim(shortname)//' '// &
!                               trim(longname)//' '// &
!                               trim(stdname)//' '// &
!                               trim(units), ESMF_LOG_INFO, rc=localrc)

            !--- add attributes ---
            call ESMF_AttributeAdd(fields(i), &
                                   convention=convCIM, purpose=purpField, &
                                   rc=localrc)

            call ESMF_AttributeSet(fields(i), &
                                   name=attrList(1), value=trim(shortname), &
                                   convention=convCIM, purpose=purpField, &
                                   rc=localrc)
            call ESMF_AttributeSet(fields(i), &
                                   name=attrList(2), value=trim(stdname), &
                                   convention=convCIM, purpose=purpField, &
                                   rc=localrc)
            call ESMF_AttributeSet(fields(i), &
                                   name=attrList(3), value=trim(longname), &
                                   convention=convCIM, purpose=purpField, &
                                   rc=localrc)
            call ESMF_AttributeSet(fields(i), &
                                   name=attrList(4), value=trim(units), &
                                   convention=convCIM, purpose=purpField, &
                                   rc=localrc)

            !--- set current and previous index ---
            preIndex = curIndex + 1
            curIndex = curIndex + index(fieldList(preIndex:), ':')

         end do

!
! KDS: There is a smarter way to do this than to just repeat this stuff...
!      maybe use a for loop instead of a do loop... is there a for loop in
!      fortran???
!
         i = i + 1
         shortname = fieldList(preIndex:curIndex - 1)

         !--- create empty ESMF_Field ---
         fields(i) = ESMF_FieldEmptyCreate(name=trim(shortname), rc=localrc)

         !--- get field metadata ---
         call seq_flds_esmf_metadata_get(shortname, longname, stdname, units)

         !--- add attributes ---
         call ESMF_AttributeAdd(fields(i), &
                                convention=convCIM, purpose=purpField, &
                                rc=localrc)

         call ESMF_AttributeSet(fields(i), &
                                name=attrList(1), value=trim(shortname), &
                                convention=convCIM, purpose=purpField, &
                                rc=localrc)
         call ESMF_AttributeSet(fields(i), &
                                name=attrList(2), value=trim(stdname), &
                                convention=convCIM, purpose=purpField, &
                                rc=localrc)
         call ESMF_AttributeSet(fields(i), &
                                name=attrList(3), value=trim(longname), &
                                convention=convCIM, purpose=purpField, &
                                rc=localrc)
         call ESMF_AttributeSet(fields(i), &
                                name=attrList(4), value=trim(units), &
                                convention=convCIM, purpose=purpField, &
                                rc=localrc)

         !--- create ESMF FieldBundle ---
         fbundle = ESMF_FieldBundleCreate(name="fbundle", rc=localrc)

         !--- add the Fields to the FieldBundle ---
         nFields = i
         do i = 1, nFields
            call ESMF_FieldBundleAdd(fbundle, (/Fields(i)/), rc=localrc)
         end do

         !--- create ESMF State using ESMF FieldBundle ---
         call ESMF_StateAdd(state, (/fbundle/), rc=localrc)

      end if

      rc = localrc

   end subroutine

end module esmfshr_attribute_mod
