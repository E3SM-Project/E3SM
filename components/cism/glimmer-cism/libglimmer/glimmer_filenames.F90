!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_filenames.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Module to handle setting a working directory for glimmer
!!
!! \author Ian Rutt
!! \date May 2007
module glimmer_filenames

  use glimmer_global,only: dirsep,fname_length

  implicit none

  character(fname_length) :: workingdir = '' !< Working directory for all file operations. Absolute paths are unaffected
  character(fname_length) :: configdir = ''  !< the directory where the config file lives and possibly other input files

contains

  !> initialise the config directory
  !!
  !! \author Magnus Hagdorn
  !! \date September 2009
  subroutine filenames_init(cname)
    implicit none
    character(len=*), intent(in) :: cname !< the configuration file name include path

    ! local variables
    integer pos
    
    ! find the last directory separator, the remaining bit is the filename
    pos = scan(trim(cname),dirsep,back=.true.)
    if (pos > 0) then
       configdir = cname(:pos)
    end if

  end subroutine filenames_init

  !> prepend path to filename
  !!
  !! \author Magnus Hagdorn
  !! \date September 2009
  !!
  !! first check if name starts with a dir sparator if so don't change name, 
  !! then check if file exists in present working directory if so do not modify file. if it doesn't exist
  !! prepend config dir
  !! \return modified file name
  function filenames_inputname(infile)
    implicit none
    character(len=*), intent(in) :: infile
    character(len=fname_length) :: filenames_inputname

    logical :: fexist

    filenames_inputname = trim(infile)

    ! check if configdir exists
    if (len(trim(configdir))  ==  0) then
       return
    end if
    ! check if path is absolute
    !! \todo figure out absolute paths for windows    
    if (infile(1:1) == dirsep) then
       return
    else
       inquire(file=infile,exist=fexist)
       ! check if the file exists in the local directory
       if (fexist) then
          return
       else
          filenames_inputname = trim(configdir)//trim(infile)
       end if
    end if
  end function filenames_inputname


  !> set the working directory
  subroutine glimmer_set_path(path)

    use glimmer_log

    character(len=*),intent(in) :: path !< the path

    workingdir=path
    call write_log('Set GLIMMER working dir to :'//trim(workingdir))

  end subroutine glimmer_set_path

  !> append path to working dir
  character(200) function process_path(path)

    character(*),intent(in) :: path !< the path to be appended

    character(200) :: alpath

    alpath=adjustl(path)

    if (alpath(1:1)/=dirsep .and. trim(workingdir)/='') then
       process_path=trim(workingdir)//dirsep//alpath
    else
       process_path=alpath
    end if

  end function process_path

  !> returns the next free file unit between 20 and 100
  integer function get_free_unit()

    use glimmer_log


    integer :: unit
    logical :: op

    unit = 20
    do
       inquire(unit,opened=op)
       if (.not.op) exit
       unit=unit+1
       if (unit>=100) then
          call write_log('No file units available',GM_FATAL,__FILE__,__LINE__)
       end if
    end do

    get_free_unit=unit

  end function get_free_unit

end module glimmer_filenames
