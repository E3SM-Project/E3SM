!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_config.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> configuration file parser
!!
!! \author Magnus Hagdorn
!! \date May 2004
!!
!! procedures used to parse configuration files. The file syntax is similar to
!! MS Windows style ini files or files that can be parsed using the Python
!! configuration file parser module. 
!!
!! The file is split up into sections. Each section appears in [] brackets. 
!! Each section can contain a number of key, value pairs. Key, value pairs are
!! separated by : or =.
!! 
!! Strings starting with any of the following characters are ignored 
!! (comments): !, # or ;
!!
!! The sections are stored in a linked list. The key-value pairs of each section
!! are also stored in linked lists. The module provides accessors to query the
!! data structure.
module glimmer_config

  use glimmer_global, only : dp, msg_length

  implicit none

  private :: handle_section, handle_value, InsertSection, InsertValue, dp

  integer, parameter :: namelen=50                 !< the maximum length of key or section
  integer, parameter :: valuelen=400               !< the maximum length of a value
  integer, parameter :: linelen=valuelen+namelen+1 !< the maximum length of a line
  
  !> derived type defining a key-value pair
  type ConfigValue
     character(len=namelen) :: name = ''         !< the key
     character(len=valuelen) :: value            !< the value
     type(ConfigValue), pointer :: next=>NULL()  !< pointer to the next key-value pair
  end type ConfigValue

  !> derived type defining a configuration section 
  type ConfigSection
     character(len=namelen) :: name = ''           !< the section name
     logical :: used = .false.                     !< flag used to check if section is used
     type(ConfigValue), pointer :: values=>NULL()  !< pointer to the first key-value pair
     type(ConfigSection), pointer :: next=>NULL()  !< pointer to the next section
  end type ConfigSection

  !> This type exists so that we can have
  !! arrays of config data, since f90 doesn't
  !! allow arrays of pointers
  type ConfigData
     type(ConfigSection), pointer :: config=>null()
  end type ConfigData

  !> generic interface for the get accessor
  interface GetValue
     module procedure GetValueDouble, GetValueReal, GetValueInt, GetValueChar, GetValueLogical, &
          GetValueDoubleArray, GetValueRealArray, GetValueIntArray, GetValueCharArray
  end interface

  !> generic interface for the set accessor
  interface ConfigSetValue
     module procedure ConfigSetValueData, ConfigSetValueSec
  end interface

  !> generic interface for the combine procedure
  interface ConfigCombine
     module procedure ConfigCombineData, ConfigCombineSec, ConfigCombineDataSec, ConfigCombineSecData
  end interface

contains

  !> read a configuration file
  subroutine ConfigRead(fname,config,fileunit)
    !*FD read configuration file
    use parallel
    use glimmer_log
    implicit none

    character(len=*), intent(in) :: fname    !< the name of the file to be read
    type(ConfigSection), pointer :: config   !< on return this pointer will point to the first section
    integer, optional,intent(in) :: fileunit !< if supplied, open this unit

    ! local variables
    type(ConfigSection), pointer :: this_section
    type(ConfigValue), pointer ::  this_value
    logical there
    integer unit,ios,linenr
    character(len=linelen) :: line
    character(len=msg_length) :: message

    if (main_task) inquire (exist=there,file=fname)
    call broadcast(there)
    if (.not.there) then
       call write_log('Cannot open configuration file '//trim(fname),GM_FATAL)
    end if
    
    unit = 99
    if (present(fileunit)) then
       unit = fileunit
    endif

    if (main_task) open(unit,file=trim(fname),status='old')
    ios=0
    linenr=0
    config=>NULL()
    this_section=>NULL()
    do while(ios == 0)
       if (main_task) read(unit,fmt='(a450)',iostat=ios) line
       call broadcast(line)
       call broadcast(ios)
       line = adjustl(line)
       if (ios /= 0) then
          exit
       end if
       if (.not.(line(1:1) == '!' .or. line(1:1) == '#' .or. line(1:1) == ';' .or. line(1:1) == ' ')) then
          ! handle comments
          if (line(1:1) == '[') then
             ! new section
             call handle_section(linenr,line,this_section)
             this_value=>NULL()
             if (.not.associated(config)) then
                ! this is the first section in config file
                config=>this_section
             end if
          else
             ! handle value
             if (.not.associated(this_section)) then
                call write_log('No section defined yet',GM_ERROR)
                write(message,*) trim(adjustl(fname)), linenr
                call write_log(message,GM_FATAL)
             end if
             call handle_value(linenr,line,this_value)
             if (.not.associated(this_section%values)) then
                this_section%values => this_value
             end if
          end if
       end if
       linenr = linenr + 1
    end do
    if (main_task) close(unit)
    return

  end subroutine ConfigRead

  !> print contents of file
  subroutine PrintConfig(config)
    implicit none
    type(ConfigSection), pointer :: config !< pointer to the first section to be printed

    type(ConfigSection), pointer :: sec
    type(ConfigValue), pointer ::  val
    
    sec=>config
    do while(associated(sec))
       write(*,*) sec%name
       val=>sec%values
       do while(associated(val))
          write(*,*) '  ',trim(val%name),' == ', trim(val%value)
          val=>val%next
       end do
       write(*,*)
       sec=>sec%next
    end do
  end subroutine PrintConfig

  !> serialise config data structure to string
  !! \author Ian Rutt
  subroutine ConfigAsString(config,string)
    use glimmer_global, only: endline
    implicit none
    type(ConfigSection), pointer :: config !< pointer to first section
    character(*),intent(out) :: string     !< on completion this string will hold the conents of the config data structure

    type(ConfigSection), pointer :: sec
    type(ConfigValue), pointer ::  val
    
    string=''

    sec=>config
    do while(associated(sec))
       string=trim(string)//'['//trim(sec%name)//']'//trim(endline)
       val=>sec%values
       do while(associated(val))
          string=trim(string)//trim(val%name)//': '//trim(val%value)//trim(endline)
          val=>val%next
       end do
       sec=>sec%next
    end do
  end subroutine ConfigAsString

  !> Either overwrite a given key-value pair,
  !! or create a new one
  !! \author Ian Rutt
  subroutine ConfigSetValueData(config,secname,valname,value,tag)

    type(ConfigData) :: config               !< 
    character(len=*), intent(in) :: secname  !< name of the section
    character(len=*), intent(in) :: valname  !< name of the key
    character(len=*), intent(in) :: value    !< the value
    character(len=*), intent(in), optional :: tag !< an identifier used to distinguish sections that occur a number of times,e.g. [CF output]       

    call ConfigSetValueSec(config%config,secname,valname,value,tag)

  end subroutine ConfigSetValueData

  !> Either overwrite a given key-value pair,
  !! or create a new one
  !! \author Ian Rutt
  subroutine ConfigSetValueSec(config,secname,valname,value,tag)

    type(ConfigSection), pointer :: config    !< pointer to the first section
    character(len=*), intent(in) :: secname  !< name of the section
    character(len=*), intent(in) :: valname  !< name of the key
    character(len=*), intent(in) :: value    !< the value
    character(len=*), intent(in), optional :: tag !< an identifier used to distinguish sections that occur a number of times,e.g. [CF output]       

    type(ConfigSection), pointer :: found
    type(ConfigSection), pointer :: newsec
    type(ConfigValue), pointer :: val
    type(ConfigValue), pointer :: newval
    type(ConfigValue), pointer :: newtag
    logical :: tagflag

    ! Find or create correct section

    if (.not.associated(config)) allocate(config)

    found=>config
    do
       if (associated(found)) then
          if (present(tag)) then
             tagflag=ConfigSectionHasTag(found,tag)
          else
             tagflag=.true.
          end if
          if ((trim(secname)==trim(found%name)).and.tagflag) then
             exit
          else
             if (associated(found%next)) then
                found=>found%next
             else
                allocate(newsec)
                found%next=>newsec
                found=>found%next
                found%name=trim(secname)
                if (present(tag)) then
                   allocate(newtag)
                   newtag%name='tag'
                   newtag%value=trim(tag)
                   found%values=>newtag
                end if
                exit
             end if
          end if
       else
          exit
       end if
    end do
 
    ! Add or create key-value pair

    if (.not.associated(found%values)) then
       allocate(newval)
       found%values=>newval
       found%values%name=valname
       found%values%value=value
    else
       val=>found%values
       do
          if (trim(valname)==trim(val%name)) then
             val%value=value
             exit
          else
             if (associated(val%next)) then
                val=>val%next
             else
                allocate(newval)
                val%next=>newval
                val%next%name=valname
                val%next%value=value
                exit
             end if
          end if
       end do
    end if

  end subroutine ConfigSetValueSec

  !> Add the contents of config2 to config1, overwriting if necessary
  !! \author Ian Rutt
  subroutine ConfigCombineDataSec(config1,config2)

    type(ConfigData) :: config1
    type(ConfigSection),pointer :: config2

    call ConfigCombineSec(config1%config,config2)

  end subroutine ConfigCombineDataSec

  !> Add the contents of config2 to config1, overwriting if necessary
  !! \author Ian Rutt
  subroutine ConfigCombineSecData(config1,config2)

    type(ConfigSection),pointer :: config1
    type(ConfigData) :: config2

    call ConfigCombineSec(config1,config2%config)

  end subroutine ConfigCombineSecData


  !> Add the contents of config2 to config1, overwriting if necessary
  !! \author Ian Rutt
  subroutine ConfigCombineData(config1,config2)

    type(ConfigData) :: config1
    type(ConfigData) :: config2

    call ConfigCombineSec(config1%config,config2%config)

  end subroutine ConfigCombineData

  !> Add the contents of config2 to config1, overwriting if necessary
  !! \author Ian Rutt
  subroutine ConfigCombineSec(config1,config2)

    type(ConfigSection), pointer :: config1
    type(ConfigSection), pointer :: config2

    type(ConfigSection), pointer :: thissec
    type(ConfigValue),   pointer :: thisval
    character(namelen) :: thisname

    character(150) :: tag

    thissec=>config2
    do
       if (associated(thissec)) then
          thisval=>thissec%values
          thisname=trim(thissec%name)
          do
             if (associated(thisval)) then
                if (ConfigSectionHasValue(thissec,'tag',tag)) then
                   call ConfigSetValue(config1,thisname,trim(thisval%name),trim(thisval%value),tag=tag)
                else
                   call ConfigSetValue(config1,thisname,trim(thisval%name),trim(thisval%value))
                end if
                thisval=>thisval%next
             else
                exit
             end if
          end do
          thissec=>thissec%next
       else
          exit
       end if
    end do

  end subroutine ConfigCombineSec

  !> check if section has specified tag
  !! \author Ian Rutt
  !!
  !! a tag is jus a special key value pair
  logical function ConfigSectionHasTag(section,tag)
    
    type(ConfigSection), pointer :: section !< pointer to section
    character(len=*),intent(in) :: tag  !< the name of the tag
    character(200) :: testtag

    ConfigSectionHasTag=.false.
    if (ConfigSectionHasValue(section,'tag',testtag)) then
       if (trim(tag)==trim(testtag)) then
          ConfigSectionHasTag=.true.
       end if
    end if

  end function ConfigSectionhasTag

  !> check if section has a particular key-value pair
  !! \author Ian Rutt
  logical function ConfigSectionHasValue(section,valname,val)

    type(ConfigSection), pointer :: section  !< pointer to the section to be checked
    type(ConfigValue), pointer :: thisval
    character(len=*), intent(in) :: valname  !< the name of the key
    character(len=*), intent(inout) :: val   !< the value

    ConfigSectionHasValue=.false.
    val=''

    if (.not.associated(section)) return

    thisval=>section%values
    do
       if (.not.associated(thisval)) exit
       if (trim(valname)==trim(thisval%name)) then
          val=trim(thisval%value)
          ConfigSectionHasValue=.true.
          exit
       else
          thisval=>thisval%next
       end if
    end do
 
  end function ConfigSectionHasValue

  !> find a return section
  !! \author Magnus Hagdorn
  subroutine GetSection(config,found,name)
    implicit none
    type(ConfigSection), pointer :: config  !< pointer to the first section
    type(ConfigSection), pointer :: found  
    character(len=*),intent(in) :: name     !< the name of the section to be found

    found=>config
    do while(associated(found))
       if (name == trim(found%name)) then
          found%used = .true.
          return
       end if
       found=>found%next
    end do
  end subroutine GetSection

  !> traverse linked list and check that all sections have been used
  subroutine CheckSections(config)
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: config
    
    ! local variables
    type(ConfigSection), pointer :: cf

    cf=>config
    do while(associated(cf))
       if (.not.cf%used) then
          call write_log('Unused section: '//trim(cf%name),GM_WARNING)
       end if
       cf=>cf%next
    end do
  end subroutine CheckSections

  !> get double array value
  subroutine GetValueDoubleArray(section,name,val,numval)
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section     !< the section from which the value is loaded
    character(len=*),intent(in) :: name         !< the name of the key
    real(kind=dp), pointer, dimension(:) :: val !< on exit this will hold the values
    integer,intent(in), optional :: numval      !< maximum number of values to be read

    ! local variables
    character(len=valuelen) :: value,tmp
    real(kind=dp), dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value == '') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tmp=value(1:ind-1)
       read(tmp,*,err=10)tempval(i)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do
    if (i >= 1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return
    
10  call write_log('Array error in config file - check syntax',GM_FATAL)

  end subroutine GetValueDoubleArray

  !> get real array value
  subroutine GetValueRealArray(section,name,val,numval)
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section !< the section from which the value is loaded
    character(len=*),intent(in) :: name     !< the name of the key
    real, pointer, dimension(:) :: val      !< on exit this will hold the values
    integer,intent(in), optional :: numval  !< maximum number of values to be read

    ! local variables
    character(len=valuelen) :: value,tmp
    real, dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value == '') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tmp=value(1:ind-1)
       read(tmp,*,err=10)tempval(i)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do

    if (i >= 1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return
    
10  call write_log('Array error in config file - check syntax',GM_FATAL)

  end subroutine GetValueRealArray

  !> get integer value array
  subroutine GetValueIntArray(section,name,val,numval)
    !*FD get integer array value
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section !< the section from which the value is loaded
    character(len=*),intent(in) :: name     !< the name of the key
    integer, pointer, dimension(:) :: val   !< on exit this will hold the value
    integer,intent(in), optional :: numval  !< maximum number of values to be read

    ! local variables
    character(len=valuelen) :: value,tmp
    integer, dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value == '') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tmp=value(1:ind-1)
       read(tmp,*,err=10)tempval(i)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do

    if (i >= 1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return
    
10  call write_log('Array error in config file - check syntax',GM_FATAL)

  end subroutine GetValueIntArray

  !> get character array value
  subroutine GetValueCharArray(section,name,val,numval)
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section         !< the section from which the value is loaded
    character(len=*),intent(in) :: name             !< the name of the key
    character(len=80), pointer, dimension(:) :: val !< on exit this will hold the values
    integer,intent(in), optional :: numval          !< maximum number of values to be read

    ! local variables
    character(len=valuelen) :: value
    character(80), dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value == '') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tempval(i)=value(1:ind-1)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do

    if (i >= 1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return
    
10  call write_log('Array error in config file - check syntax',GM_FATAL)

  end subroutine GetValueCharArray

  !> get real value
  subroutine GetValueReal(section,name,val)
    implicit none
    type(ConfigSection), pointer :: section !< the section from which the value is loaded
    character(len=*),intent(in) :: name     !< the name of the key
    real, intent(inout) :: val              !< the value

    ! local variables
    character(len=valuelen) :: value
    real temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueReal

  !> get double value
  subroutine GetValueDouble(section,name,val)
    implicit none
    type(ConfigSection), pointer :: section !< the section from which the value is loaded
    character(len=*),intent(in) :: name     !< the name of the key
    real(kind=dp), intent(inout) :: val     !< the value

    ! local variables
    character(len=valuelen) :: value

    real(kind=dp) :: temp

    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueDouble

  !> get integer value
  subroutine GetValueInt(section,name,val)
    implicit none
    type(ConfigSection), pointer :: section !< the section from which the value is loaded
    character(len=*),intent(in) :: name     !< the name of the key
    integer, intent(inout) :: val           !< the value

    ! local variables
    character(len=valuelen) :: value
    integer temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueInt

  !> get character value
  subroutine GetValueChar(section,name,val)
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section !< the section from which the value is loaded
    character(len=*),intent(in) :: name     !< the name of the key
    character(len=*), intent(inout) :: val  !< the value
    
    type(ConfigValue), pointer :: value

    value=>section%values
    do while(associated(value))
       if (name == trim(value%name)) then
          val = value%value
          if ((len_trim(val) + 1) >= len(val)) then 
            ! Assume that if we get within one space of the variable length (excluding spaces) then we may be truncating the intended value.
            call write_log('The value of config option   ' // trim(name) // '   is too long for the variable.' ,GM_FATAL)
          endif
          return
       end if
       value=>value%next
    end do
  end subroutine GetValueChar

  !> get logical value
  subroutine GetValueLogical(section,name,val)
    implicit none
    type(ConfigSection), pointer :: section !< the section from which the value is loaded
    character(len=*),intent(in) :: name     !< the name of the key
    logical, intent(inout) :: val           !< the value

    ! local variables
    character(len=valuelen) :: value
    integer itemp
    logical ltemp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) itemp
    if (ios==0) then
       val = itemp == 1
    end if
    read(value,*,iostat=ios) ltemp
    if (ios==0) then
       val = ltemp
    end if
  end subroutine GetValueLogical

  !==================================================================================
  ! private procedures
  !==================================================================================

  !> handle line in file containing a section
  subroutine handle_section(linenr,line,section)
    use glimmer_log
    implicit none
    integer, intent(in) :: linenr            !< the line number
    character(len=*), intent(in) :: line     !< buffer containing the line
    type(ConfigSection), pointer :: section  !< pointer to place where new section should be inserted

    ! local variables
    integer i
    character(len=msg_length) :: message

    do i=1,linelen
       if (line(i:i) == ']') then
          exit
       end if
    end do
    if (line(i:i) /= ']') then
       write(message,*) 'Cannot find end of section ',linenr
       call write_log(message,GM_FATAL)
    end if

    call InsertSection(trim(adjustl(line(2:i-1))),section)
  end subroutine handle_section
  
  !> handle line in file containing a key-value pair
  subroutine handle_value(linenr,line,value)
    use glimmer_log
    implicit none
    integer, intent(in) :: linenr        !< the line number
    character(len=*), intent(in) :: line !< buffer containing the line
    type(ConfigValue), pointer :: value  !< pointer to value linked list where value should be added

    ! local variables
    integer i
    character(len=msg_length) :: message
    do i=1,linelen
       if (line(i:i) == '=' .or. line(i:i) == ':') then
          exit
       end if
    end do
    if (.not.(line(i:i) == '=' .or. line(i:i) == ':')) then
       write(message,*) 'Cannot find = or : ',linenr
       call write_log(message,GM_FATAL)
    end if

    call InsertValue(trim(adjustl(line(:i-1))), trim(adjustl(line(i+1:))),value)
  end subroutine handle_value

  !> add a new section
  subroutine InsertSection(name,section)
    !*FD add a new section
    implicit none
    character(len=*), intent(in) :: name    !< name of new section
    type(ConfigSection), pointer :: section !< on entry the element of linked list after which the new element is inserted, on exit: the new element
    type(ConfigSection), pointer :: new_sec

    allocate(new_sec)
    new_sec%name = name

    if (associated(section)) then
       if (associated(section%next)) then
          new_sec%next => section%next
       end if
       section%next=>new_sec
    end if
    section=>new_sec
  end subroutine InsertSection

  !> insert a key-value pair
  subroutine InsertValue(name,val,value)
    use glimmer_log
    implicit none
    character(len=*), intent(in) :: name !< the key
    character(len=*), intent(in) :: val  !< the value
    type(ConfigValue), pointer :: value  !< on entry the element after which the new value should be added, on exit pointer the new element
    type(ConfigValue), pointer :: new_value

    allocate(new_value)
    
    ! Assume that if we get within one space of the variable length (excluding spaces) then we may be truncating the intended value.
    if ((len_trim(val) + 1) >= len(new_value%value)) then 
       call write_log('The value of config option   ' // trim(name) // '   is too long to be read fully.' ,GM_FATAL)  
    endif

    new_value%name = name
    new_value%value = val

    if(associated(value)) then
       if (associated(value%next)) then
          new_value%next => value%next
       end if
       value%next => new_value
    end if
    value=>new_value
  end subroutine InsertValue
end module glimmer_config
