module rrlw_ncpar
	use parkind ,only : im => kind_im, rb => kind_rb

	implicit none
        save
	
        real(kind=rb), parameter :: cpdair = 1003.5  ! Specific heat capacity of dry air
                                        		 ! at constant pressure at 273 K
                                        		 ! (J kg-1 K-1)

	
	integer(kind=im), parameter :: maxAbsorberNameLength = 5, &
                          Absorber              = 12
    character(len = maxAbsorberNameLength), dimension(Absorber), parameter :: &
    AbsorberNames = (/        &
     				'N2   ',  &
     				'CCL4 ',  &
     				'CFC11',  &
     				'CFC12',  &
     				'CFC22',  &
     				'H2O  ',  &
     				'CO2  ',  &
     				'O3   ',  &
     				'N2O  ',  & 
     				'CO   ',  &
     				'CH4  ',  &
     				'O2   '  /)
	
	integer(kind=im), dimension(40) :: status
	integer(kind=im) :: i
	integer(kind=im), parameter :: keylower  = 9,   &
						  keyupper  = 5,   &
						  Tdiff     = 5,   &
						  ps        = 59,  &
						  plower    = 13,  &
						  pupper    = 47,  &
						  Tself     = 10,  &
						  Tforeign  = 4,   &
						  pforeign  = 4,   &
						  T         = 19,  &
						  Tplanck   = 181, &
						  band      = 16,  &
						  GPoint    = 16,  &
						  GPointSet = 2
						  
	contains 
	
	subroutine getAbsorberIndex(AbsorberName,AbsorberIndex)
		character(len = *), intent(in) :: AbsorberName
		integer(kind=im), intent(out)           :: AbsorberIndex
		
		integer(kind=im) :: m
	
		AbsorberIndex = -1
		do m = 1, Absorber
			if (trim(AbsorberNames(m)) == trim(AbsorberName)) then
				AbsorberIndex = m
			end if
		end do
		
		if (AbsorberIndex == -1) then
			print*, "Absorber name index lookup failed."
		end if
	end subroutine getAbsorberIndex

end module rrlw_ncpar
