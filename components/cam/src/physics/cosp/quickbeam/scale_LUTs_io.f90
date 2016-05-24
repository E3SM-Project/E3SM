  ! scale_LUT_io:  Contains subroutines to load and save scaling Look Up Tables (LUTs) to a file
  ! 
  ! June 2010   Written by Roj Marchand
  
  module scale_LUTs_io
  implicit none

  contains

  subroutine load_scale_LUTs(hp)
  
    use radar_simulator_types

    type(class_param), intent(inout) :: hp

    logical :: LUT_file_exists
    integer :: i,j,k,ind
    
    !
    ! load scale LUT from file 
    !
    inquire(file=trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat', &
        exist=LUT_file_exists)

    if(.not.LUT_file_exists) then
    
        write(*,*) '*************************************************'
        write(*,*) 'Warning: Could NOT FIND radar LUT file: ', &
        trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'        
        write(*,*) 'Will calculated LUT values as needed'
        write(*,*) '*************************************************'
        
        return
    else

        OPEN(unit=12,file=trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat',&
        form='unformatted', &
        err= 89, &
            access='DIRECT',&
            recl=28)
         
            write(*,*) 'Loading radar LUT file: ', &
        trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    
            do i=1,maxhclass
            do j=1,mt_ntt
            do k=1,nRe_types
        
            ind = i+(j-1)*maxhclass+(k-1)*(nRe_types*mt_ntt)
            
            read(12,rec=ind) hp%Z_scale_flag(i,j,k), &
                    hp%Ze_scaled(i,j,k), &
                    hp%Zr_scaled(i,j,k), &
                    hp%kr_scaled(i,j,k)
                                                    
                ! if(ind==1482329) then
                !   write (*,*) ind, hp%Z_scale_flag(i,j,k), &
                !   hp%Ze_scaled(i,j,k), &
                !   hp%Zr_scaled(i,j,k), &
                !   hp%kr_scaled(i,j,k)
                !endif
     
            enddo
            enddo
            enddo
        
            close(unit=12)
        return 
    endif
    
  89    write(*,*) 'Error: Found but could NOT READ radar LUT file: ', &
        trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    stop
    
  end subroutine load_scale_LUTs
  
  subroutine save_scale_LUTs(hp)

    use radar_simulator_types
  
    type(class_param), intent(inout) :: hp
    
    logical :: LUT_file_exists
    integer :: i,j,k,ind
    
    inquire(file=trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat', &
        exist=LUT_file_exists)

    OPEN(unit=12,file=trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat',&
        form='unformatted', &
        err= 99, &
            access='DIRECT',&
            recl=28)
         
        write(*,*) 'Creating or Updating radar LUT file: ', &
        trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    
        do i=1,maxhclass
        do j=1,mt_ntt
        do k=1,nRe_types
        
            ind = i+(j-1)*maxhclass+(k-1)*(nRe_types*mt_ntt)
            
            if(.not.LUT_file_exists .or. hp%Z_scale_added_flag(i,j,k)) then
            
                hp%Z_scale_added_flag(i,j,k)=.false.
            
                write(12,rec=ind) hp%Z_scale_flag(i,j,k), &
                    hp%Ze_scaled(i,j,k), &
                    hp%Zr_scaled(i,j,k), &
                    hp%kr_scaled(i,j,k)
                     
                !  1482329 T  0.170626345026495        0.00000000000000       1.827402935860823E-003
            
                !if(ind==1482329) then
                !   write (*,*) ind, hp%Z_scale_flag(i,j,k), &
                !   hp%Ze_scaled(i,j,k), &
                !   hp%Zr_scaled(i,j,k), &
                !   hp%kr_scaled(i,j,k)
                !endif
            endif
        enddo
        enddo
        enddo
        
        close(unit=12)
    return 
    
  99    write(*,*) 'Error: Unable to create/update radar LUT file: ', &
        trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    return  
    
  end subroutine save_scale_LUTs


  end module scale_LUTs_io
  
