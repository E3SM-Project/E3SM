module e4d_setup

  use e4d_vars

  implicit none
  integer :: ios
  real :: xorig,yorig,zorig

contains
  !________________________________________________________________
  subroutine setup_e4d
    implicit none
    integer :: sz
  
    
    if (my_rank>0) then
       write(*,*) my_rank,"here1!!!!!!!!!!!!!!!!!!!!!!!!!!!1"
       call slave_setup
       return
    end if

    log_file = 'e4d_'//trim(adjustl(pflotran_group_prefix))//'.log'   
    open(13,file=trim(log_file),action='write',status='replace')
    close(13)
!#if 0
    call read_input 

    if (ios .eq. -1) then
       call send_command(0)
       return
    end if
    open(13,file=trim(log_file),action='write',status='old',position='append')
    write(13,*) "Read Input Complete___________________________________"
    close(13)

    call send_info

    call setup_forward
    open(13,file=trim(log_file),action='write',status='old',position='append')
    write(13,*) "Setup Forward Complete___________________________________"
    close(13)
    if (ios .eq. -1) then
       call send_command(0)
       return
    end if
  
!#endif
!    call send_info_geh

    call send_command(0)
    open(13,file=trim(log_file),action='write',status='old',position='append')
    write(13,*) "Setup Phase Complete___________________________________"
    close(13)
  end subroutine setup_e4d
  !________________________________________________________________

  !________________________________________________________________
  subroutine read_input
    implicit none
    
    integer :: nchar,junk,i,check,j,nchr,npre,flg,flnum
    integer :: a,b,m,n
    real :: ex,ey,ez,eflg,vi,wd,tsig
    real, dimension(4) :: etmp 
    logical :: exst

    !get the name of the mesh file
    call elog(0,0,flg)
    if (flg .ne. 0) then
       ios = -1
       return
    end if

    

    open(10,file='e4d.inp',status='old',action='read')
    read(10,*,IOSTAT=ios) mode
    call elog(40,ios,flg) 
    if (flg .eq. -1) then
       ios=-1
       return
    end if
    if(mode == 33) tank_flag = .true.

    read(10,*,IOSTAT=ios) mshfile
    call elog(1,ios,flg); 
    if (flg .eq. -1) then
       ios=-1
       return
    end if

    read(10,*,IOSTAT=ios) efile
    call elog(2,ios,flg); 
    if (flg .eq. -1) then
       ios=-1
       return
    end if

    read(10,*,IOSTAT=ios) sigfile
    call elog(3,ios,flg); 
    if (flg .eq. -1) then
       ios=-1
       return
    end if

    read(10,*,IOSTAT=ios) fmn_file
    call elog(4,ios,flg); 
    if (flg .eq. -1) then
       ios=-1
       return
    end if

    read(10,*,IOSTAT=ios) list_file
    call elog(23,ios,flg)
    if (flg .eq. -1) then
       ios=-1
       return
    end if
    close(10)

    !get the number of character is in the prefix
    nchar = 0
    do i=1,40
       if (mshfile(i:i) == '.') then
          nchar = i
          exit
       end if
    end do
    
    !check to see if the files exist

    !!Allocate/read the electrode positions and survey conf
     open(10,file=efile,status='old',action='read')   
     read(10,*,IOSTAT=ios) ne
     call elog(5,ios,flg)
     if (ios .ne. 0) then
        ios=-1
        return
     end if
       
       
     allocate(e_pos(ne,4))
     do i=1,ne
        read(10,*,IOSTAT=ios) junk,etmp
        if (ios .ne. 0) then
           call elog(6,i,flg)
           ios=-1
           return
        end if
        e_pos(junk,1:4)=etmp
     end do
     
     !get the meshfile prefix
     nchr=len_trim(mshfile)
     do i=1,nchr
        if (mshfile(i:i)=='.') then
           npre=i-1;
           exit
        end if
     end do

     !!translate electrodes
     call elog(7,npre,flg)
     if (flg.eq.-1) then
        ios=-1
        return
     end if

     open(21,file=mshfile(1:npre)//".trn",status="old",action="read")
     read(21,*,IOSTAT=ios) xorig,yorig,zorig
     if (ios .ne. 0) then
        call elog(8,ios,flg)
        ios=-1
        return
     end if
     close(21)
     e_pos(:,1) = e_pos(:,1)-xorig
     e_pos(:,2) = e_pos(:,2)-yorig
     e_pos(:,3) = e_pos(:,3)-zorig

     !!Read in the survey
     read(10,*,IOSTAT=ios) nm
     call elog(9,ios,flg)
     if (ios .ne. 0) then
        ios = -1
        return
     end if

     allocate(s_conf(nm,4))
     do i=1,nm
        read(10,*,IOSTAT=ios) junk,s_conf(i,1:4)
        if (ios.ne.0) then
           call elog(10,i,flg)
           ios=-1
           return
        end if
     end do
     close(10) 
       
       !!Read the conductivity file
       open(10,file=trim(sigfile),status='old',action='read')
       read(10,*,IOSTAT=ios) j !FF,sw_sig,gw_sig
       call elog(11,ios,j)
       if (ios .ne. 0) then
          ios=-1
          close(10)
          return
       end if

       allocate(base_sigma(j))
       do i=1,j
          read(10,*,IOSTAT=ios) base_sigma(i)
          if (ios .ne. 0) then
             call elog(12,i,flg)
             ios=-1
             close(10)
             return
          end if
       end do
       close(10)

       !make sure pf_mesh.txt exists
       call elog(37,ios,flg)
       if(ios.ne.0) then
          ios=-1
          close(1)
          return
       end if

       !Read in the FMN file
       open(10,file=trim(fmn_file),status='old',action='read')
       read(10,*,IOSTAT=ios) j 
       call elog(38,ios,j)
       if (ios .ne. 0) then
          ios=-1
          close(10)
          return
       end if
       if(allocated(FMN)) deallocate(FMN)
       allocate(FMN(j,4))
    
       do i=1,j
           read(10,*,IOSTAT=ios) FMN(i,1:4)
           if (ios .ne. 0) then
              call elog(39,ios,i)
              ios=-1
              close(10)
              return
          end if
       end do
       close(10)

       !!read the list file and all files listed and check for errors
       open(10,file=trim(list_file),status='old',action='read')
       read(10,*,IOSTAT=ios) ntime !,FF,sw_sig,gw_sig
       call elog(24,ios,flg)
       if (flg==-1) then
          ios=-1
          close(10)
          return
       end if
       flnum=0
       do j=1,ntime
          flnum=flnum+1
          read(10,*,IOSTAT=ios) e4d_time,csrv_file,ccond_file
          call elog(25,ios,flnum)
          if (ios .ne. 0) then
             ios=-1
             close(10)
             return
          end if
          call elog(26,ios,flnum)
          if (flnum .ne. 0) then
             ios=-1
             close(10)
             return
          end if
       end do
          
       
  end subroutine read_input
  !________________________________________________________________

  !________________________________________________________________
  subroutine send_info
    implicit none

    call send_command(7)
    call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    call MPI_BCAST(s_conf, 4*nm,MPI_INTEGER,0,E4D_COMM,ierr)

  end subroutine send_info
  !________________________________________________________________

  !____________________________________________________________________
  subroutine setup_forward
    implicit none
    
    integer :: npre,i,nchr,dim,bflag,ns,itmp,ios,flg
    real :: jnk,jnk1 
    integer :: status(MPI_STATUS_SIZE)
    
    !command slaves to do the setup
    call send_command(2)
   
    !read and send the nodes
    !!Determine the prefix of the tetgen output files
    nchr=len_trim(mshfile)
    do i=1,nchr
       if (mshfile(i:i)=='.') then
          npre=i+1;
          exit
       end if
    end do

    !!OPEN AND READ THE NODE POSITIONS 
    !!nodes(x,y,z) and nbounds
    call elog(15,npre,flg)
    if (flg==-1) then
       ios=-1
       return
    end if

    open(10,file=mshfile(1:npre)//".node",status="old",action="read")
    read(10,*,IOSTAT=ios) nnodes,dim,jnk,bflag
    if (ios .ne. 0) then
       call elog(16,ios,flg)
       ios=-1
       return
    end if
  
    allocate(nodes(nnodes,3),nbounds(nnodes))
   
    do i=1,nnodes
       read(10,*,IOSTAT=ios) jnk,nodes(i,1:3),jnk1,nbounds(i)
       if (ios .ne. 0) then
          call elog(17,i,flg)
          ios = -1
          return
       end if
       if (nbounds(i)>2) nbounds(i)=7
    end do
    close(10)

    !!Send the nodes to the slaves
    call MPI_BCAST(nnodes, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, E4D_COMM, ierr )
    call MPI_BCAST(nbounds, nnodes, MPI_INTEGER , 0, E4D_COMM, ierr )
   
    !!Read in the elements and send each slave it's assigned elements
    call elog(18,npre,flg)
    if (flg==-1) then
       ios=-1
       return
    end if

    open(10,file=mshfile(1:npre)//".ele",status="old",action="read")
    read(10,*,IOSTAT=ios) nelem,dim,jnk
    if (ios .ne. 0) then
       call elog(19,ios,flg)
       ios=-1
       return
    end if

    allocate(elements(nelem,4),zones(nelem))
    do i=1,nelem
       read(10,*,IOSTAT=ios) jnk,elements(i,1:4),zones(i)
       if (ios .ne. 0) then
          call elog(20,i,flg)
          ios = -1
          return
       end if
    end do
    close(10)
 
    call MPI_BCAST(nelem, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,E4D_COMM,ierr)
 
    !send the assignments
    call send_dists

    !build the mesh interpolation matrix
    call get_mesh_interp
 
    !get electrode nodes
    call get_electrode_nodes 
 
    call elog(21,ios,flg)
    ios=0
  end subroutine setup_forward
  !____________________________________________________________________
 
  !____________________________________________________________________
  subroutine get_mesh_interp
    implicit none
    integer :: i,ierr,st,en
    integer :: status(MPI_STATUS_SIZE)
    integer, dimension(n_rank-1) :: tnmap
    
    nmap=0
    
    do i=1,n_rank-1
       call MPI_RECV(tnmap(i),1,MPI_INTEGER,i,i,E4D_COMM,status,ierr)      
    end do
   
    nmap=sum(tnmap)

    open(13,file=trim(log_file),status='old',action='write',position='append')
    write(13,*) "Total number of mapping weights: ",nmap
    do i=1,n_rank-1
       write(13,*) "Rank ",i,": ",tnmap(i)
    end do
    close(13)
    
  
    allocate(map_inds(nmap,2),map(nmap))
    
    do i=1,n_rank-1
       if(i==1) then
          st=1
          en=tnmap(1)
          call MPI_RECV(map_inds(st:en,1),tnmap(i),MPI_INTEGER,i,i,E4D_COMM,status,ierr)
          call MPI_RECV(map_inds(st:en,2),tnmap(i),MPI_INTEGER,i,i,E4D_COMM,status,ierr)
          call MPI_RECV(map(st:en),tnmap(i),MPI_REAL,i,i,E4D_COMM,status,ierr)
       else
          st=sum(tnmap(1:i-1))+1
          en=sum(tnmap(1:i))
          call MPI_RECV(map_inds(st:en,1),tnmap(i),MPI_INTEGER,i,i,E4D_COMM,status,ierr)
          call MPI_RECV(map_inds(st:en,2),tnmap(i),MPI_INTEGER,i,i,E4D_COMM,status,ierr)
          call MPI_RECV(map(st:en),tnmap(i),MPI_REAL,i,i,E4D_COMM,status,ierr)
        
       end if
    end do
  

  end subroutine get_mesh_interp
  !____________________________________________________________________
 !____________________________________________________________________
  subroutine send_dists
    !read inputs and distribute the run info to the slave
    implicit none
    integer :: neven,nextra,ce,i,j,k

    
    !!divide the electrodes up among the slave processes
    !!note ne is the number of electrodes defined in read_inp
    tne = ne  
    neven = ne/(n_rank-1)
    nextra = ne - neven*(n_rank-1)  
    if (allocated(eind)) deallocate(eind)
    allocate(eind(n_rank-1,2))
    
    !!Build eind
    if (neven > 0) then
       ce = 1
       do i=1,n_rank-1-nextra
          eind(i,1) = ce
          eind(i,2) = ce+neven-1
          ce = ce+neven
       end do
       
       if (nextra > 0) then
          do i=n_rank-nextra,n_rank-1    
             eind(i,1) = ce
             eind(i,2) = ce+neven
             ce=ce+neven+1
          end do
       end if
    else
       !!in this case there are more processors than electrodes
       !!THIS SHOULD BE AVOIDED
       do i=1,n_rank-1
          eind(i,1) = i
          eind(i,2) = i
       end do
    end if

    !!Build jaco assignments
    if(allocated(jind)) deallocate(jind)
    allocate(jind(n_rank-1,2))
    jind = 0
    if(nelem < (n_rank-1)) then
       do i=1,nelem
          jind(i,1:2)=i
       end do
    else
       neven = nelem/(n_rank-1)
       nextra = nelem - neven*(n_rank-1)
       
       jind=0
       if(neven > 0) then
          ce = 1
          do i=1,n_rank-1-nextra
             jind(i,1) = ce
             jind(i,2) = ce+neven-1
             ce = ce+neven
          end do
          
          if(nextra > 0) then
             do i=n_rank-nextra,n_rank-1    
                jind(i,1) = ce
                jind(i,2) = ce+neven
                ce=ce+neven+1
             end do
          end if
       end if
    end if
    
    !this will be replaced with mesh info from pflotran
    call get_pfmesh
 
    !!send assignments
    call MPI_BCAST(tne, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(e_pos,4*tne ,MPI_REAL , 0, E4D_COMM, ierr )
    call MPI_BCAST(eind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(jind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(pfnx,1,MPI_INTEGER,0,E4D_COMM, ierr)
    call MPI_BCAST(pfny,1,MPI_INTEGER,0,E4D_COMM, ierr)
    call MPI_BCAST(pfnz,1,MPI_INTEGER,0,E4D_COMM, ierr)
    call MPI_BCAST(pfxcb,pfnx+1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(pfycb,pfny+1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(pfzcb,pfnz+1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(mode, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    deallocate(pfxcb)
    deallocate(pfycb)
    deallocate(pfzcb)
  

  end subroutine send_dists
  !______________________________________________________________

 
  !____________________________________________________________________
  subroutine get_pfmesh
    implicit none
    integer :: i
    logical :: chech
    

    open(12,file='pf_mesh.txt',status='old',action='read')
    read(12,*) pfnx
    allocate(pfxcb(pfnx+1))
    read(12,*) pfxcb(1:pfnx+1)
  
    read(12,*) pfny
    allocate(pfycb(pfny+1))
    read(12,*) pfycb(1:pfny+1)

    read(12,*) pfnz
    allocate(pfzcb(pfnz+1))
    read(12,*) pfzcb(1:pfnz+1)

    close(12)
    open(13,file=trim(log_file),action='write',status='old',position='append')
    write(13,*) "________________PF MESH PARAMETERS___________________"
    write(13,*) "NX, NY, NZ: ",pfnx,pfny,pfnz
    write(13,*) "xmin  xmax: ",pfxcb(1),pfxcb(pfnx+1)
    write(13,*) "ymin  ymax: ",pfycb(1),pfycb(pfny+1)
    write(13,*) "zmin  zmax: ",pfzcb(1),pfzcb(pfnz+1)
    write(13,*)
    close(13)
    
  end subroutine get_pfmesh
  !____________________________________________________________________

  !______________________________________________________________
  subroutine get_electrode_nodes
    implicit none
    real, dimension(nnodes) :: dist
    integer :: i
    integer, dimension(1) :: indb
   
    if (.not.allocated(e_nods)) then
       allocate (e_nods(tne))
    end if
  
    do i=1,tne
       dist=sqrt( (e_pos(i,1) - nodes(:,1))**2 + &
            (e_pos(i,2) - nodes(:,2))**2 + &
            (e_pos(i,3) - nodes(:,3))**2)
       
       indb = minloc(dist)
       
       if (my_rank==1) then
          if (dist(indb(1)) > 0.1) then
             !!There should always be a node at all
             !!electrode locations. If not, something is
             !!wrong
             write(*,1001) i,dist(indb)
          end if
       end if
       e_nods(i) = indb(1)

1001   format(" WARNING: ELECTRODE ",I5," HAS BEEN MOVED ",F10.4," TO THE NEAREST NODE")
    end do
    i_zpot(1) = e_nods(tne)

  end subroutine get_electrode_nodes
  !____________________________________________________________________

  !______________________________________________________________
  subroutine send_command(com)
    !!Send a general command to the slaves
    integer :: com
    
    call MPI_BCAST(com,1,MPI_INTEGER,0,E4D_COMM,ierr)
    
  end subroutine send_command
  !________________________________________________________________
  
  !__________________SLAVE ROUTINES________________________________
  !________________________________________________________________
  subroutine slave_setup
    
    implicit none
    integer :: command
    
100 continue
    !Recieve a command from master
    call MPI_BCAST(command,1,MPI_INTEGER,0,E4D_COMM,ierr)
    
    
    !return to main
    if (command == 0) then
       return

    else if (command == 2) then
       call setup_frun
       goto 100

    else if (command == 7) then
       call receive_info
       goto 100
       
    else if (command == 86) then
       call receive_info_geh
       goto 100
       
    else
       goto 100

    end if
  end subroutine slave_setup
  !________________________________________________________________

  !__________________________________________________________________
  subroutine receive_info
    implicit none
    integer :: i
    
    if (allocated(s_conf)) deallocate(s_conf)
    call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    allocate(s_conf(nm,4))
    call MPI_BCAST(s_conf, 4*nm,MPI_INTEGER,0,E4D_COMM,ierr)
    
  end subroutine receive_info
  !__________________________________________________________________
  
  !__________________________________________________________________
  subroutine setup_frun
    implicit none

#include "petsc/finclude/petscsys.h"
    
    integer :: status(MPI_STATUS_SIZE)
    integer :: neven, nextra, ce, i,itmp
    integer :: beg,end
    logical :: eqf = .false.
    real, dimension(:), allocatable :: dist
    integer, dimension(1) :: tmp
    real :: mx,my,mz
    
  
    !receive nodes from master
    call MPI_BCAST(nnodes, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
  

    allocate(nodes(nnodes,3),nbounds(nnodes))
    call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, E4D_COMM, ierr )
    call MPI_BCAST(nbounds, nnodes, MPI_INTEGER , 0, E4D_COMM, ierr )

    !receive elements from master
    call MPI_BCAST(nelem, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    allocate(elements(nelem,4))
    call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,E4D_COMM,ierr)

    !receive the assignments
    call receive_dists
    
    !build the mesh interpolation matrix for my elements
    call send_my_mesh_interp

    !build A_map and delA.......................................
    call build_delA

    !!Initialize the Petsc A matrix
    call MatCreateSeqAIJ(PETSC_COMM_SELF,nnodes,nnodes,d_nz,d_nnz,A, &
                         perr);CHKERRQ(perr)
    call MatSetFromOptions(A,perr);CHKERRQ(perr)
 
    !Get the electrode ownership indexes 
    call get_electrode_nodes
    call MatGetType(A,tp,perr);CHKERRQ(perr)
    
    !Set up the source and solution vectors
    call VecCreate(PETSC_COMM_SELF,X,perr);CHKERRQ(perr)
    call VecSetSizes(X,nnodes,PETSC_DECIDE,perr);CHKERRQ(perr)
    call VecSetFromOptions(X,perr);CHKERRQ(perr)
 
    !allocate and setup the pole solution vector
    call VecCreate(PETSC_COMM_SELF,psol,perr);CHKERRQ(perr)
    call VecSetSizes(psol,nnodes,PETSC_DECIDE,perr);CHKERRQ(perr)
    call VecSetFromOptions(psol,perr);CHKERRQ(perr)
    allocate(poles(nnodes,my_ne))
    poles = 0

    call VecCreate(PETSC_COMM_SELF,B,perr);CHKERRQ(perr)
    call VecSetSizes(B,nnodes,PETSC_DECIDE,perr);CHKERRQ(perr)
    call VecSetFromOptions(B,perr);CHKERRQ(perr)

    !A_map  and S_map given the coupling info so we can deallocate the 
    !nodes and elements
    deallocate(nodes)
    deallocate(elements)
    allocate(sigma(nelem))
    
  end subroutine setup_frun
  !__________________________________________________________________
  
  !__________________________________________________________________
  subroutine receive_dists
    implicit none
    
    integer :: i
   
    if (allocated(eind)) deallocate(eind)
    if (allocated(e_pos)) deallocate(e_pos)
    if (allocated(jind)) deallocate(jind)
    allocate(eind(n_rank-1,2))
    allocate(jind(n_rank-1,2))
    
    call MPI_BCAST(tne, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    allocate(e_pos(tne,4))
    call MPI_BCAST(e_pos,4*tne,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(eind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(jind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(pfnx,1,MPI_INTEGER,0,E4D_COMM, ierr)
    call MPI_BCAST(pfny,1,MPI_INTEGER,0,E4D_COMM, ierr)
    call MPI_BCAST(pfnz,1,MPI_INTEGER,0,E4D_COMM, ierr)
    allocate(pfxcb(pfnx+1),pfycb(pfny+1),pfzcb(pfnz+1))
    call MPI_BCAST(pfxcb,pfnx+1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(pfycb,pfny+1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(pfzcb,pfnz+1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(mode, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    if(mode==33) tank_flag = .true.
    my_ne = eind(my_rank,2)-eind(my_rank,1)+1
   
  end subroutine receive_dists
  !__________________________________________________________________
  !____________________________________________________________________
  subroutine send_my_mesh_interp
    implicit none
    integer :: i,j,my_nelem,cnt,n_tets,indx,indy,indz,k,cct,cpos
    real :: xmn,xmx,ymn,ymx,zmn,zmx,imn,imx
    logical, dimension(:), allocatable :: flags
    real, dimension(:), allocatable :: midx,midy,midz,w
    integer, dimension(:), allocatable :: rw,v
    real, dimension(9,3) :: pts
    real, dimension(4,3) :: mfaces
    real, dimension(3) :: vc
    real, dimension(8) :: wi,C
    integer, dimension(8) :: vi
    real, dimension(72) :: twi
    integer, dimension(72) :: tvi

    
    !!get the pf mesh midpoints
    xmn=1e15
    ymn=xmn
    zmn=xmn
    xmx=-1e15
    ymx=xmx
    zmx=xmx
   

    allocate(midx(pfnx),midy(pfny),midz(pfnz))
    do i=1,pfnx
       midx(i) = 0.5*(pfxcb(i+1)+pfxcb(i))
       if(midx(i)<xmn) xmn=midx(i)
       if(midx(i)>xmx) xmx=midx(i)
      
    end do
   
    do i=1,pfny
       midy(i) = 0.5*(pfycb(i+1)+pfycb(i))
       if(midy(i)<ymn) ymn=midy(i)
       if(midy(i)>ymx) ymx=midy(i)
    
    end do

    do i=1,pfnz
       midz(i) = 0.5*(pfzcb(i+1)+pfzcb(i))
       if(midz(i)<zmn) zmn=midz(i)
       if(midz(i)>zmx) zmx=midz(i)

    end do

    nodes(:,1)=nodes(:,1)+xorig
    nodes(:,2)=nodes(:,2)+yorig
    nodes(:,3)=nodes(:,3)+zorig
    
    my_nelem = jind(my_rank,2)-jind(my_rank,1)+1
    allocate(flags(my_nelem))

    flags=.true.
    cnt=0
    n_tets=0
   

    !!find which elements we need to map
    do i=jind(my_rank,1),jind(my_rank,2)
       cnt=cnt+1
       imx=maxval(nodes(elements(i,1:4),1))
       imn=minval(nodes(elements(i,1:4),1))
       if((imn<xmn .or. imx>xmx) .and. (pfnx > 1)) then
          flags(cnt)=.false.

          goto 100
       end if
     
       imx=maxval(nodes(elements(i,1:4),2))
       imn=minval(nodes(elements(i,1:4),2))
       if((imn<ymn .or. imx>ymx) .and. (pfny > 1)) then
          flags(cnt)=.false.
          goto 100
       end if
  

       imx=maxval(nodes(elements(i,1:4),3))
       imn=minval(nodes(elements(i,1:4),3)) 
    
       if((imn<zmn .or. imx>zmx) .and. (pfnz>1)) then
          flags(cnt)=.false.
          goto 100
       end if

       n_tets=n_tets+1
100    continue
    end do
 
    allocate(rw(72*n_tets),v(72*n_tets),w(72*n_tets))
    cnt=0
    cct=0
    cpos=0
    do i=jind(my_rank,1),jind(my_rank,2)
       cct=cct+1
       if(flags(cct)) then
          !get the nine points
          cnt=cnt+1
          do j=1,3
             pts(1,j)=.25*sum(nodes(elements(i,1:4),j)) !tet midpoint
          end do
        
          do j=1,3
             mfaces(1,j) = sum(nodes(elements(i,[1,2,3]),j))/3
             mfaces(2,j) = sum(nodes(elements(i,[1,2,4]),j))/3
             mfaces(3,j) = sum(nodes(elements(i,[1,3,4]),j))/3
             mfaces(4,j) = sum(nodes(elements(i,[2,3,4]),j))/3       
          end do
          do j=1,4
             vc=mfaces(j,:)-pts(1,:)
             pts(j+1,:) = pts(1,:) + 0.5*vc
             vc=nodes(elements(i,j),:)-pts(1,:)
             pts(j+5,:) = pts(1,:) + 0.5*vc
          end do
          

          !get the weighting function for each point
          do j=1,9
             if(pfnx .eq. 1) then
                indx = 1
             else
                do k=1,pfnx
                   if(midx(k)>pts(j,1)) then
                      indx=k-1
                      exit
                   end if
                end do
             end if

             if(pfny.eq.1) then
                indy = 1
             else
                do k=1,pfny
                   if(midy(k)>pts(j,2)) then
                      indy=k-1
                      exit
                   end if
                end do
             end if
             
             if(pfnz .eq.1) then
                indz = 1
             else
                do k=1,pfnz
                   if(midz(k)>pts(j,3)) then
                      indz=k-1
                      exit
                   end if
                end do
             end if

             
             if(pfny .eq.1) then
                !for 2D problem in x,z
                vi(5) = indx + (indy-1)*pfnx + (indz-1)*pfnx*pfny
                vi(6) = vi(5)
                vi(8) = vi(5)+1
                vi(7) = vi(8);
                vi(1) = vi(5)+pfnx
                vi(2) = vi(1);
                vi(3) = vi(1)+1;
                vi(4) = vi(3);
                
             else
                vi(5) = indx + (indy-1)*pfnx + (indz-1)*pfnx*pfny
                vi(8) = vi(5)+1
                vi(6) = vi(5)+pfnx
                vi(7) = vi(6)+1
                vi(1:4) = vi(5:8) + pfnx*pfny
             end if

              if(pfnx .eq.1) then
                 C(1)=0
              else
                 C(1)=(pts(j,1)-midx(indx))/(midx(indx+1)-midx(indx));
              end if
              C(2)=C(1);
              C(3)=C(1);
              C(4)=C(1);

              if(pfny .eq.1) then
                 C(5)=0
              else
                 C(5)=(pts(j,2)-midy(indy))/(midy(indy+1)-midy(indy));
              end if
              C(6)=C(5);
              
              if(pfnz .eq.1) then
                 C(7)=0
              else
                 C(7)=(pts(j,3)-midz(indz))/(midz(indz+1)-midz(indz));
              end if

              wi(1)=(1-C(1))*(1-C(6))*C(7);
              wi(2)=(1-C(2))*C(6)*C(7);
              wi(3)=C(2)*C(6)*C(7);
              wi(4)=C(1)*C(7)*(1-C(6));
              wi(5)=(1-C(3))*(1-C(5))*(1-C(7));
              wi(6)=C(5)*(1-C(4))*(1-C(7));
              wi(7)=C(4)*C(5)*(1-C(7));
              wi(8)=C(3)*(1-C(5))*(1-C(7));
           
              !write(*,*) vi
              !write(*,*) wi

              tvi(8*(j-1)+1:8*j)=vi
              twi(8*(j-1)+1:8*j)=wi
              
           end do
           
           !consolidate the weights
           do j=1,71
              do k=j+1,72
                 if(tvi(k)==tvi(j)) then
                    twi(j)=twi(j)+twi(k)
                    twi(k)=0
                    tvi(k)=0;
                 end if
              end do
           end do
           twi=twi/9
          
           do j=1,72
              if(tvi(j) .ne. 0) then
                 cpos=cpos+1;
                 rw(cpos)=i
                 v(cpos)=tvi(j)
                 w(cpos)=twi(j)
              end if
           end do
        end if
     
    end do
  
    call MPI_SEND(cpos,1,MPI_INTEGER,0,my_rank,E4D_COMM, ierr)
    call MPI_SEND(rw(1:cpos),cpos,MPI_INTEGER,0,my_rank,E4D_COMM,ierr)
    call MPI_SEND(v(1:cpos),cpos,MPI_INTEGER,0,my_rank,E4D_COMM,ierr)
    call MPI_SEND(w(1:cpos),cpos,MPI_REAL,0,my_rank,E4D_COMM,ierr)
    
    deallocate(flags,midx,midy,midz,rw,v,w)
    nodes(:,1)=nodes(:,1)-xorig
    nodes(:,2)=nodes(:,2)-yorig
    nodes(:,3)=nodes(:,3)-zorig
    
  end subroutine send_my_mesh_interp
  !____________________________________________________________________
!!$  !____________________________________________________________________
!!$  subroutine send_my_mesh_interp
!!$    implicit none
!!$    integer :: i,j,my_nelem,cnt,n_tets,indx,indy,indz,k,cct,cpos
!!$    real :: xmn,xmx,ymn,ymx,zmn,zmx,imn,imx
!!$    logical, dimension(:), allocatable :: flags
!!$    real, dimension(:), allocatable :: midx,midy,midz,w
!!$    integer, dimension(:), allocatable :: rw,v
!!$    real, dimension(9,3) :: pts
!!$    real, dimension(4,3) :: mfaces
!!$    real, dimension(3) :: vc
!!$    real, dimension(8) :: wi,C
!!$    integer, dimension(8) :: vi
!!$    real, dimension(72) :: twi
!!$    integer, dimension(72) :: tvi
!!$
!!$    
!!$    !!get the pf mesh midpoints
!!$    xmn=1e15
!!$    ymn=xmn
!!$    zmn=xmn
!!$    xmx=-1e15
!!$    ymx=xmx
!!$    zmx=xmx
!!$   
!!$
!!$    allocate(midx(pfnx),midy(pfny),midz(pfnz))
!!$    do i=1,pfnx
!!$       midx(i) = 0.5*(pfxcb(i+1)+pfxcb(i))
!!$       if(midx(i)<xmn) xmn=midx(i)
!!$       if(midx(i)>xmx) xmx=midx(i)
!!$      
!!$    end do
!!$   
!!$    do i=1,pfny
!!$       midy(i) = 0.5*(pfycb(i+1)+pfycb(i))
!!$       if(midy(i)<ymn) ymn=midy(i)
!!$       if(midy(i)>ymx) ymx=midy(i)
!!$    
!!$    end do
!!$
!!$    do i=1,pfnz
!!$       midz(i) = 0.5*(pfzcb(i+1)+pfzcb(i))
!!$       if(midz(i)<zmn) zmn=midz(i)
!!$       if(midz(i)>zmx) zmx=midz(i)
!!$
!!$    end do
!!$
!!$    nodes(:,1)=nodes(:,1)+xorig
!!$    nodes(:,2)=nodes(:,2)+yorig
!!$    nodes(:,3)=nodes(:,3)+zorig
!!$    
!!$    my_nelem = jind(my_rank,2)-jind(my_rank,1)+1
!!$    allocate(flags(my_nelem))
!!$
!!$    flags=.true.
!!$    cnt=0
!!$    n_tets=0
!!$   
!!$
!!$    !!find which elements we need to map
!!$    do i=jind(my_rank,1),jind(my_rank,2)
!!$       cnt=cnt+1
!!$       imx=maxval(nodes(elements(i,1:4),1))
!!$       imn=minval(nodes(elements(i,1:4),1))
!!$       if((imx<pfxcb(1) .or. imn>pfxcb(pfnx+1)) .and. (pfnx > 1)) then
!!$          flags(cnt)=.false.
!!$
!!$          goto 100
!!$       end if
!!$     
!!$       imx=maxval(nodes(elements(i,1:4),2))
!!$       imn=minval(nodes(elements(i,1:4),2))
!!$       if((imx<pfycb(1) .or. imn>pfycb(pfny+1)) .and. (pfny > 1)) then
!!$          flags(cnt)=.false.
!!$          goto 100
!!$       end if
!!$  
!!$
!!$       imx=maxval(nodes(elements(i,1:4),3))
!!$       imn=minval(nodes(elements(i,1:4),3))     
!!$       if((imx<pfzcb(1) .or. imn>pfzcb(pfnz+1)) .and. (pfnz>1)) then
!!$          flags(cnt)=.false.
!!$          goto 100
!!$       end if
!!$
!!$       n_tets=n_tets+1
!!$100    continue
!!$    end do
!!$
!!$    allocate(rw(72*n_tets),v(72*n_tets),w(72*n_tets))
!!$    cnt=0
!!$    cct=0
!!$    cpos=0
!!$    do i=jind(my_rank,1),jind(my_rank,2)
!!$       cct=cct+1
!!$       if(flags(cct)) then
!!$          !get the nine points
!!$          cnt=cnt+1
!!$          do j=1,3
!!$             pts(1,j)=.25*sum(nodes(elements(i,1:4),j)) !tet midpoint
!!$          end do
!!$        
!!$          do j=1,3
!!$             mfaces(1,j) = sum(nodes(elements(i,[1,2,3]),j))/3
!!$             mfaces(2,j) = sum(nodes(elements(i,[1,2,4]),j))/3
!!$             mfaces(3,j) = sum(nodes(elements(i,[1,3,4]),j))/3
!!$             mfaces(4,j) = sum(nodes(elements(i,[2,3,4]),j))/3       
!!$          end do
!!$          do j=1,4
!!$             vc=mfaces(j,:)-pts(1,:)
!!$             pts(j+1,:) = pts(1,:) + 0.5*vc
!!$             vc=nodes(elements(i,j),:)-pts(1,:)
!!$             pts(j+5,:) = pts(1,:) + 0.5*vc
!!$          end do
!!$          
!!$
!!$          !get the weighting function for each point
!!$          do j=1,9
!!$             indx = pfnx-1
!!$             if(pfnx .eq. 1) then
!!$                indx = 1
!!$             else
!!$                do k=1,pfnx
!!$                   if(midx(k)>pts(j,1)) then
!!$                      indx=k-1
!!$                      exit
!!$                   end if
!!$                end do
!!$                if(indx .eq. 0) indx=1
!!$                !if(indx .eq. -1) indx=pfnx
!!$             end if
!!$
!!$             indy = pfny-1
!!$             if(pfny.eq.1) then
!!$                indy = 1
!!$             else
!!$                do k=1,pfny
!!$                   if(midy(k)>pts(j,2)) then
!!$                      indy=k-1
!!$                      exit
!!$                   end if
!!$                end do
!!$                if(indy .eq. 0)  indy=1
!!$                !if(indy .eq. -1) indy=pfny
!!$             end if
!!$             
!!$             indz=pfnz-1
!!$             if(pfnz .eq.1) then
!!$                indz = 1
!!$             else
!!$                do k=1,pfnz
!!$                   if(midz(k)>pts(j,3)) then
!!$                      indz=k-1
!!$                      exit
!!$                   end if
!!$                end do
!!$                if(indz .eq. 0) indz=1
!!$                !if(indz .eq.-1) indz=pfnz 
!!$             end if
!!$
!!$             
!!$             if(pfny .eq.1) then
!!$                !for 2D problem in x,z
!!$                vi(5) = indx + (indy-1)*pfnx + (indz-1)*pfnx*pfny
!!$                vi(6) = vi(5)
!!$                vi(8) = vi(5)+1
!!$                vi(7) = vi(8);
!!$                vi(1) = vi(5)+pfnx
!!$                vi(2) = vi(1);
!!$                vi(3) = vi(1)+1;
!!$                vi(4) = vi(3);
!!$                
!!$             else
!!$                vi(5) = indx + (indy-1)*pfnx + (indz-1)*pfnx*pfny
!!$                vi(8) = vi(5)+1
!!$                vi(6) = vi(5)+pfnx
!!$                vi(7) = vi(6)+1
!!$                vi(1:4) = vi(5:8) + pfnx*pfny
!!$             end if
!!$
!!$              if(pfnx .eq.1) then
!!$                 C(1)=0
!!$              else
!!$                 C(1)=(pts(j,1)-midx(indx))/(midx(indx+1)-midx(indx));
!!$              end if
!!$              C(2)=C(1);
!!$              C(3)=C(1);
!!$              C(4)=C(1);
!!$
!!$              if(pfny .eq.1) then
!!$                 C(5)=0
!!$              else
!!$                 C(5)=(pts(j,2)-midy(indy))/(midy(indy+1)-midy(indy));
!!$              end if
!!$              C(6)=C(5);
!!$              
!!$              if(pfnz .eq.1) then
!!$                 C(7)=0
!!$              else
!!$                 C(7)=(pts(j,3)-midz(indz))/(midz(indz+1)-midz(indz));
!!$              end if
!!$
!!$              wi(1)=(1-C(1))*(1-C(6))*C(7);
!!$              wi(2)=(1-C(2))*C(6)*C(7);
!!$              wi(3)=C(2)*C(6)*C(7);
!!$              wi(4)=C(1)*C(7)*(1-C(6));
!!$              wi(5)=(1-C(3))*(1-C(5))*(1-C(7));
!!$              wi(6)=C(5)*(1-C(4))*(1-C(7));
!!$              wi(7)=C(4)*C(5)*(1-C(7));
!!$              wi(8)=C(3)*(1-C(5))*(1-C(7));
!!$           
!!$              !write(*,*) vi
!!$              !write(*,*) wi
!!$
!!$              tvi(8*(j-1)+1:8*j)=vi
!!$              twi(8*(j-1)+1:8*j)=wi
!!$              
!!$           end do
!!$           
!!$           !consolidate the weights
!!$           do j=1,71
!!$              do k=j+1,72
!!$                 if(tvi(k)==tvi(j)) then
!!$                    twi(j)=twi(j)+twi(k)
!!$                    twi(k)=0
!!$                    tvi(k)=0;
!!$                 end if
!!$              end do
!!$           end do
!!$           twi=twi/9
!!$          
!!$           do j=1,72
!!$              if(tvi(j) .ne. 0) then
!!$                 cpos=cpos+1;
!!$                 rw(cpos)=i
!!$                 v(cpos)=tvi(j)
!!$                 w(cpos)=twi(j)
!!$              end if
!!$           end do
!!$        end if
!!$     
!!$    end do
!!$  
!!$    call MPI_SEND(cpos,1,MPI_INTEGER,0,my_rank,E4D_COMM, ierr)
!!$    call MPI_SEND(rw(1:cpos),cpos,MPI_INTEGER,0,my_rank,E4D_COMM,ierr)
!!$    call MPI_SEND(v(1:cpos),cpos,MPI_INTEGER,0,my_rank,E4D_COMM,ierr)
!!$    call MPI_SEND(w(1:cpos),cpos,MPI_REAL,0,my_rank,E4D_COMM,ierr)
!!$    
!!$    deallocate(flags,midx,midy,midz,rw,v,w)
!!$    nodes(:,1)=nodes(:,1)-xorig
!!$    nodes(:,2)=nodes(:,2)-yorig
!!$    nodes(:,3)=nodes(:,3)-zorig
!!$    
!!$  end subroutine send_my_mesh_interp
!!$  !____________________________________________________________________
  !____________________________________________________________________
  subroutine build_delA

    use e4d_mat_inv_module

    implicit none
    
   
    integer :: nnds,i,j,k,l,mcount,rn,cn,ncl,jj
    integer :: lrow
    real, dimension(4,4) :: atild,atildi,temp_a
    integer, dimension(4,4) :: indx
    integer, dimension(:,:),allocatable :: A_tmp
    real :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,evol
    real :: A_kij,lA_kij
    real, dimension(:), allocatable :: t1vec,t2vec
    logical :: ilow,iup
    
    !allocate the petsc matrix preallocation vectors
    allocate(d_nnz(nnodes))
    d_nnz=0
    
    allocate(rows(10*nelem),cols(10*nelem))
    allocate(A_map(10*nelem),S_map(10*nelem))
    allocate(delA(10*nelem))
    
    !!COUNT THE NUMBER OF NODE PAIRS (I.E. COUPLING MATRIX ELEMENTS) AND
    !!REALLOCATE THE COUPLING MATRIX STORAGE VECTORS
    !!These pairs will also accomodate source dependent boundaries
    rows(1)=elements(1,1)
    cols(1)=elements(1,1)
    nvals=1
    mcount = 0
    A_map = 0
    nnds=maxval(elements(:,:))
    
    !!Allocate a temporary matrix to aid in building the 
    !!mapping vector
    allocate(A_tmp(nnds,80))
    A_tmp=0
    A_tmp(rows(1),1) = 1
    A_tmp(rows(1),2) = nvals
    A_tmp(rows(1),3) = cols(1) 
    !!loop over the elements
     
    do k=1,nelem    
       
       !!loop over each node in this element
       do i=1,4
         
          !!rn is a row index of the coupling matrix
          rn=elements(k,i)    
          
          !!loop over each node in this element
          do j=1,4
             
             !!cn is a column index of the coupling matrix
             cn=elements(k,j)
        
          

             !!check to see if this pair (rn,cn) is in the lower
             !!triangle. If not, go to the next pair
             if (cn <= rn) then
                mcount=mcount+1
                !!loop over each pair found thus far (nvals pairs)
                !!and determine if the pair (rn,cn) is already 
                !!represeneted
                if (A_tmp(rn,1) == 0) then
                   nvals = nvals + 1
                   A_tmp(rn,1) = 1
                   A_tmp(rn,2) = nvals
                   A_tmp(rn,3) = cn
                   rows(nvals) = rn
                   cols(nvals) = cn
                   A_map(mcount) = nvals
                           
                else
                   
                   do l=1,A_tmp(rn,1)
                      if (cn == A_tmp(rn,2*l+1)) then
                         A_map(mcount) = A_tmp(rn,2*l)
                         goto 11
                      end if
                   end do
                   
                   !!if were here then no column index was found
                   !!for this row, so we'll add one
                   nvals = nvals+1
                   ncl = A_tmp(rn,1) + 1
                   A_tmp(rn,1) = ncl
                   A_tmp(rn,2*ncl) = nvals
                   A_tmp(rn,2*ncl+1) = cn
                   rows(nvals) = rn
                   cols(nvals) = cn
                   A_map(mcount) = nvals
                
11                 continue
     
                end if

             end if

          end do
          
       end do
       
    end do
    
    !build the petsc allocation vector
    do i=1,nnodes
       !upper part
       d_nnz(i)=A_tmp(i,1)-1;
       do j=1,A_tmp(i,1)
          cn=A_tmp(i,2*j+1)
          d_nnz(cn)=d_nnz(cn)+1
       end do
    end do   
    deallocate(A_tmp)
   
    !!now reallocate for the correct number of matrix values (i.e. pairs)
    allocate(t1vec(nvals),t2vec(nvals))
    t1vec(1:nvals)=rows(1:nvals)
    t2vec(1:nvals)=cols(1:nvals)
    deallocate(rows,cols)
    allocate(rows(nvals),cols(nvals))
    rows=t1vec
    cols=t2vec
    deallocate(t1vec,t2vec)
    allocate(trows(nvals),tcols(nvals))
    
    !!LOOP OVER THE ELEMENTS AND BUILD delA
    mcount = 0
    delA=0;
  
    do k=1,nelem
   
       !!get the 4 nodes for this element 
       x1 = nodes(elements(k,1),1)
       y1 = nodes(elements(k,1),2)
       z1 = nodes(elements(k,1),3)
       
       x2 = nodes(elements(k,2),1)
       y2 = nodes(elements(k,2),2)
       z2 = nodes(elements(k,2),3)
       
       x3 = nodes(elements(k,3),1)
       y3 = nodes(elements(k,3),2)
       z3 = nodes(elements(k,3),3) 
       
       x4 = nodes(elements(k,4),1)
       y4 = nodes(elements(k,4),2)
       z4 = nodes(elements(k,4),3)
       
       !!get the volume of this element
       call get_vol (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,evol)  
   
       if (evol==0) goto 12
       !!build the linear shape function matrix for this element
       atild(1,:) = (/1,1,1,1/)
       atild(2,:) = (/x1,x2,x3,x4/)
       atild(3,:) = (/y1,y2,y3,y4/)
       atild(4,:) = (/z1,z2,z3,z4/)
       temp_a=atild
            
       !!invert temp_a to get atildi
       !!note MIGS routine is located in mat_inv.f      
       call MIGS(temp_a,4,atildi,indx)  
      
12     continue
       !!use atildi to build the coupling coefficients
       !!for this element      
       
       !!loop over each node in this element
       do i=1,4
          !!rn is a row index of the coupling matrix
          rn=elements(k,i)
          
          !!loop over each node in this element
          do j=1,4
             
             !!cn is a column index of the coupling matrix
             cn=elements(k,j)
             
             !!check to see if this pair (rn,cn) is in the lower
             !!triangle. If not, go to the next pair
             if (cn <= rn) then
                mcount=mcount+1
                !!Compute the coupling coefficient for this element
                !!and node pair
                A_kij=0
               
                do l=2,4
                   A_kij=A_kij+(atildi(i,l)*atildi(j,l))
                end do
          
                !!delA(mcount)=lA_kij*evol;
                delA(mcount)=dble(A_kij*evol)
                S_map(mcount)=k
             end if
             
21           continue
          end do
          
       end do
       
    end do  
    
    
  end subroutine build_delA
  !____________________________________________________________________
  
  !____________________________________________________________________
  subroutine get_vol( x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, volume )

    implicit none
    
    real a(4,4)
    real r4_det
    real volume
    real x1
    real x2
    real x3
    real x4
    real y1
    real y2
    real y3
    real y4
    real z1
    real z2
    real z3
    real z4
    
    a(1:4,1) = (/ x1, x2, x3, x4 /)
    a(1:4,2) = (/ y1, y2, y3, y4 /)
    a(1:4,3) = (/ z1, z2, z3, z4 /)
    a(1:4,4) = (/ 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00 /)
    r4_det = det4(a)
    volume = abs ( r4_det ) / 6.0E+00
    
    return
  end subroutine get_vol
  !____________________________________________________________________
    
  !____________________________________________________________________
  function det4 ( a1 )
    implicit none
    
    real :: a1(4,4)
    real*8 :: a(4,4)
    real :: det4
    integer ::i
    real*8 :: c1,c2,c3,c4
    
    a = dble(a1)
    c1=a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*&
         &(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) 
    
    c2=a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*&
         &(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))
    
    c3=a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*&
         &(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
    
    c4=a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(2,2)*&
         &(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
    
    det4 = real(c1 - c2 + c3 - c4)
    
    return
  end function det4
  !__________________________________________________________________________

  !________________________________________________________________
  subroutine send_info_geh
    implicit none

    call send_command(86)
    nelem = 4*4*4
    nm = nelem
    print *, 'master:', nm
    call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)

  end subroutine send_info_geh
  !________________________________________________________________

  !__________________________________________________________________
  subroutine receive_info_geh
    implicit none
    integer :: i
    
    call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    allocate(sigma(nm))
    nelem = nm
    sigma = 0.d0
    !print *, 'slave:', nm
    
  end subroutine receive_info_geh
  !__________________________________________________________________
  
  !________________________________________________________________
  subroutine destroy_e4d
    implicit none
    if (allocated(e4d_ranks)) then
      deallocate(e4d_ranks)
    endif
    if (allocated(pf_e4d_ranks)) then
      deallocate(pf_e4d_ranks)
    endif
    if (allocated(map_inds)) then
      deallocate(map_inds)
    endif
    if (allocated(s_conf)) then
      deallocate(s_conf)
    endif
    if (allocated(eind)) then
      deallocate(eind)
    endif
    if (allocated(nbounds)) then
      deallocate(nbounds)
    endif
    if (allocated(zones)) then
      deallocate(zones)
    endif
    if (allocated(elements)) then
      deallocate(elements)
    endif
    if (allocated(faces)) then
      deallocate(faces)
    endif
    if (allocated(e_nods)) then
      deallocate(e_nods)
    endif
    if (allocated(rows)) then
      deallocate(rows)
    endif
    if (allocated(cols)) then
      deallocate(cols)
    endif
    if (allocated(trows)) then
      deallocate(trows)
    endif
    if (allocated(tcols)) then
      deallocate(tcols)
    endif
    if (allocated(A_map)) then
      deallocate(A_map)
    endif
    if (allocated(S_map)) then
      deallocate(S_map)
    endif
    if (allocated(my_drows)) then
      deallocate(my_drows)
    endif
    if (allocated(e_pos)) then
      deallocate(e_pos)
    endif
    if (allocated(nodes)) then
      deallocate(nodes)
    endif
    if (allocated(poles)) then
      deallocate(poles)
    endif
    if (allocated(pf_porosity)) then
      deallocate(pf_porosity)
    endif
    if (allocated(pf_tracer)) then
      deallocate(pf_tracer)
    endif
    if (allocated(pf_saturation)) then
      deallocate(pf_saturation)
    endif
    if (allocated(pf_saturation_0)) then
      deallocate(pf_saturation_0)
    endif
    if (allocated(pf_temperature)) then
      deallocate(pf_temperature)
    endif
    if (allocated(sigma)) then
      deallocate(sigma)
    endif
    if (allocated(dpred)) then
      deallocate(dpred)
    endif
    if (allocated(dobs)) then
      deallocate(dobs)
    endif
    if (allocated(sd)) then
      deallocate(sd)
    endif
    if (allocated(my_dvals)) then
      deallocate(my_dvals)
    endif
    if (allocated(map)) then
      deallocate(map)
    endif
    if (allocated(base_sigma)) then
      deallocate(base_sigma)
    endif
    if (allocated(d_nnz)) then
      deallocate(d_nnz)
    endif
    if (allocated(delA)) then
      deallocate(delA)
    endif
  end subroutine destroy_e4d
  !________________________________________________________________
 
end module e4d_setup

