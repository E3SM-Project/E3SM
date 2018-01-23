module rrsw_ncpar
	use parkind ,only : im => kind_im, rb => kind_rb

	implicit none
    save
	
    real(kind=rb), parameter :: cpdair = 1003.5  ! Specific heat capacity of dry air
                                        		 ! at constant pressure at 273 K
                                        		 ! (J kg-1 K-1)

	integer(kind=im), dimension(50) :: status
	integer(kind=im) :: i
	integer(kind=im), parameter :: keylower      = 9,  &
						  keyupper      = 5,  &
						  Tdiff         = 5,  &
						  ps            = 59, &
						  plower        = 13, &
						  pupper        = 47, &
						  Tself         = 10, &
						  Tforeignlower = 3,  &
						  Tforeignupper = 2,  &
						  pforeign      = 4,  &
						  T             = 19, &
						  band          = 14, &
						  GPoint        = 16, &
						  GPointSet     = 2
	
	integer(kind=im), parameter :: maxAbsorberNameLength   = 5,  &
                          Absorber                = 12, &
                          maxKeySpeciesNameLength = 3,  &
                          maxKeySpeciesNames      = 2
                          
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
	
	character(len = maxKeySpeciesNameLength), dimension(band,maxKeySpeciesNames), parameter :: &
    KeySpeciesNamesLower = RESHAPE( SOURCE = (/ 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', &
							  					'H2O', 'H2O', 'H2O', '   ', 'O3 ', 'O3 ', 'H2O', &
							  					'CH4', 'CO2', 'CH4', 'CO2', '   ', 'CO2', 'O2 ', &
							  					'   ', 'O2 ', '   ', '   ', '   ', 'O2 ', '   '  /), &
							  		SHAPE = (/ band, maxKeySpeciesNames /) )
							  
	character(len = maxKeySpeciesNameLength), dimension(band,maxKeySpeciesNames), parameter :: &
    KeySpeciesNamesUpper = RESHAPE( SOURCE = (/ 'CH4', 'H2O', 'CH4', 'CO2', 'H2O', 'H2O', 'O2 ', &
							 					'   ', 'O2 ', '   ', '   ', 'O3 ', 'O3 ', 'CO2', &
							  					'   ', 'CO2', '   ', '   ', '   ', 'CO2', '   ', &
							  					'   ', '   ', '   ', '   ', '   ', 'O2 ', '   '  /), &
							  		SHAPE = (/ band, maxKeySpeciesNames /) )
							
	integer(kind=im), dimension(band)     :: BandNums = (/ 16, 17, 18, 19, 20, 21, 22, &
										  		        23, 24, 25, 26, 27, 28, 29 /)
										      
	real(kind=rb), dimension(keylower) :: KeySpeciesLower = (/ 1.0, 0.125, 0.25, 0.375, &
											      			   0.50, 0.625, 0.75, 0.875, 1.0 /)
											          
	real(kind=rb), dimension(keyupper) :: KeySpeciesUpper = (/ 0.0, 0.25, 0.50, 0.75, 1.0 /)
		
	real(kind=rb), dimension(Tdiff)    :: TempDiffs = (/ -30, -15, 0, 15, 30 /)
										      
	real(kind=rb), dimension(Tself)    :: TempSelf = (/ 245.6,252.8,260.0,267.2,274.4, &
														281.6,288.8,296.0,303.2,310.4 /)		
	
	real(kind=rb), dimension(Tforeignlower) :: TempForeignlower = (/ 296, 260, 224 /)
	
	real(kind=rb), dimension(Tforeignupper) :: TempForeignupper = (/ 224, 260 /)
	
	real(kind=rb), dimension(pforeign) :: PressForeign = (/ 970, 475, 219, 3 /)
			
	real(kind=rb), dimension(T)        :: Temp = (/188.0, 195.2, 202.4, 209.6, 216.8, 224.0, &
								   				   231.2, 238.4, 245.6, 252.8, 260.0, 267.2, &
								   				   274.4, 281.6, 288.8, 296.0, 303.2, 310.4, 317.6 /)
	
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

end module rrsw_ncpar
