module e4d_run

  use e4d_vars
  integer :: pf_com
  logical :: first_sol = .true.
  logical :: sim_e4d = .false.
  integer :: mcomm
  !real*8 :: pf_time
contains
  
  !_____________________________________________________________________
  subroutine run_e4d
  
    implicit none
    logical :: first_flag = .true.
    
    if (my_rank>0) then
       call slave_run
       call cleanup
       return
    end if

    if (.not. allocated(pf_porosity)) allocate(pf_porosity(pflotran_vec_size))  
    if (.not. allocated(pf_tracer)) allocate(pf_tracer(pflotran_vec_size))
    if (.not. allocated(pf_saturation)) &
      allocate(pf_saturation(pflotran_vec_size))
    if (.not. allocated(pf_saturation_0)) &
      allocate(pf_saturation_0(pflotran_vec_size))
    ! if energy is being modeled, pflotran_temperature_vec_mpi will be non-zero
    if (pflotran_temperature_vec_mpi /= 0 .and. &
        .not. allocated(pf_temperature)) &
      allocate(pf_temperature(pflotran_vec_size))
    if (.not. allocated(sigma)) allocate(sigma(nelem))
   

    call get_mcomm
    call cpu_time(Cbeg)
    call get_pf_porosity       !!get the pflotran porosity
  
    do while (mcomm==1)

     
       call get_pf_time      !!get the pflotran solution time
       call get_pf_sol       !!get the pflotran solution
      
       call elog(36,mcomm,mcomm)

       call check_e4d_sim    !!see if we should do an e4d sim for this time
       
       if (sim_e4d) then

          if (first_flag) then
             pf_saturation_0 = pf_saturation
             !call compute_FF
             first_flag = .false.
          end if
 
          call map_pf_e4d       !!map/transform the solution to the E4D mesh
          call send_sigma       !!send the transformed solution to the slaves           
          call send_command(3)  !!instruct slaves to build A matrix
          call send_command(5)  !!instruct slaves to build KSP solver
          call send_command(6)  !!instruct slaves to solve  
         call get_dpred        !!assemble the simulated data
          call cpu_time(Cbeg)
       end if
     
       call get_mcomm
      
    end do

    call send_command(0)  !!instruct slaves to exit
    call cleanup
    return

     
  end subroutine run_e4d
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine check_e4d_sim
    implicit none
    integer :: tmp,i,ios,j,junk,a,b,m,n
   

    open(13,file=trim(log_file),status='old',action='write',position='append')
    write(13,*) "Min/Max Fluid Cond. Value: ",minval(pf_tracer),maxval(pf_tracer)
    close(13)
    sim_e4d = .false.
    
    open(10,file=trim(list_file),status='old',action='read')
    read(10,*) tmp
    do i=1,ntime
       read(10,*) e4d_time,csrv_file,ccond_file
       if (real(e4d_time) .eq. real(pf_time)) then
          call elog(35,i,tmp)
          close(10)
          
          open(10,file=trim(ccond_file),status='old',action='read')
          read(10,*,IOSTAT=ios) nsig
          if (ios .ne. 0)       call elog(31,ios,tmp)
          if (nsig .ne. nelem)  call elog(32,nsig,nelem)
          if (.not.allocated(sigma)) allocate(sigma(nelem))
          do j=1,nelem
             read(10,*) sigma(j)
          end do
          close(10)

          open(10,file=trim(csrv_file),status='old',action='read')
          read(10,*,IOSTAT=ios) tmp
          call elog(27,ios,tmp)
          do j=1,tmp
             read(10,*) junk
          end do
          read(10,*,IOSTAT=ios) tmp
          call elog(29,ios,tmp)
          if (.not. allocated(dobs)) allocate(dobs(nm))
          if (.not. allocated(sd)) allocate(sd(nm))
          do j=1,nm
             read(10,*,IOSTAT=ios) tmp,a,b,m,n,dobs(j),sd(j)
          end do 
          close(10)
          sim_e4d = .true.
          return
       end if
    end do

    close(10)
    open(13,file=trim(log_file),status='old',action='write',position='append')
    write(13,*) "No E4D survey found for pflotran time: ",pf_time
    close(13)
  end subroutine check_e4d_sim
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_mcomm
    implicit none
!geh    call MPI_BCAST(mcomm,1,MPI_INTEGER,0,PFE4D_COMM,ierr)
    call MPI_BCAST(mcomm,1,MPI_INTEGER,0,PFE4D_MASTER_COMM,ierr)
  end subroutine get_mcomm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_pf_time
    implicit none
    call MPI_BCAST(pf_time,1,MPI_DOUBLE_PRECISION,0,PFE4D_MASTER_COMM, &
                   ierr)
   
  end subroutine get_pf_time
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_pf_porosity
    implicit none
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
    integer ::  status(MPI_STATUS_SIZE)
    PetscReal, pointer :: vec_ptr(:)

 
    ! porosity
    ! we actually hijack the tracer vec to transfer porosity
    call VecScatterBegin(pflotran_scatter,pflotran_tracer_vec_mpi, &
                         pflotran_tracer_vec_seq, &
                         INSERT_VALUES,SCATTER_FORWARD,perr);CHKERRQ(perr)
    call VecScatterEnd(pflotran_scatter,pflotran_tracer_vec_mpi, &
                       pflotran_tracer_vec_seq, &
                       INSERT_VALUES,SCATTER_FORWARD,perr);CHKERRQ(perr)
    call VecGetArrayF90(pflotran_tracer_vec_seq,vec_ptr,perr);CHKERRQ(perr)
    pf_porosity = real(vec_ptr)

    call VecRestoreArrayF90(pflotran_tracer_vec_seq,vec_ptr, &
                            perr);CHKERRQ(perr)

  end subroutine get_pf_porosity
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_pf_sol
    implicit none
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
    integer ::  status(MPI_STATUS_SIZE)
    PetscReal, pointer :: vec_ptr(:)

 
    ! tracer
    call VecScatterBegin(pflotran_scatter,pflotran_tracer_vec_mpi, &
                         pflotran_tracer_vec_seq, &
                         INSERT_VALUES,SCATTER_FORWARD,perr);CHKERRQ(perr)
    call VecScatterEnd(pflotran_scatter,pflotran_tracer_vec_mpi, &
                       pflotran_tracer_vec_seq, &
                       INSERT_VALUES,SCATTER_FORWARD,perr);CHKERRQ(perr)
    call VecGetArrayF90(pflotran_tracer_vec_seq,vec_ptr,perr);CHKERRQ(perr)
    pf_tracer = real(vec_ptr)

    call VecRestoreArrayF90(pflotran_tracer_vec_seq,vec_ptr, &
                            perr);CHKERRQ(perr)

    ! saturation                
    call VecScatterBegin(pflotran_scatter,pflotran_saturation_vec_mpi, &
                         pflotran_saturation_vec_seq, &
                         INSERT_VALUES,SCATTER_FORWARD,perr);CHKERRQ(perr)
    call VecScatterEnd(pflotran_scatter,pflotran_saturation_vec_mpi, &
                       pflotran_saturation_vec_seq, &
                       INSERT_VALUES,SCATTER_FORWARD,perr);CHKERRQ(perr)
    call VecGetArrayF90(pflotran_saturation_vec_seq,vec_ptr,perr);CHKERRQ(perr)
    pf_saturation = real(vec_ptr)
    call VecRestoreArrayF90(pflotran_saturation_vec_seq,vec_ptr, &
                            perr);CHKERRQ(perr)

    ! temperature  (only modeled when energy is simulated)
    if (pflotran_temperature_vec_mpi /= 0) then
      call VecScatterBegin(pflotran_scatter,pflotran_temperature_vec_mpi, &
                           pflotran_temperature_vec_seq, &
                           INSERT_VALUES,SCATTER_FORWARD,perr);CHKERRQ(perr)
      call VecScatterEnd(pflotran_scatter,pflotran_temperature_vec_mpi, &
                         pflotran_temperature_vec_seq, &
                         INSERT_VALUES,SCATTER_FORWARD,perr);CHKERRQ(perr)
      call VecGetArrayF90(pflotran_temperature_vec_seq,vec_ptr, &
                          perr);CHKERRQ(perr)
      pf_temperature = real(vec_ptr)
      call VecRestoreArrayF90(pflotran_temperature_vec_seq,vec_ptr, &
                            perr);CHKERRQ(perr)
    endif
 
  end subroutine get_pf_sol
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine map_pf_e4d
    
    implicit none
    integer :: i,cnt,pi,ei
    character*40 :: filename, word
    real :: a,m,n,tc,delSigb
    logical :: tcorr_flag = .false.   
   
    !set the conductivity values that are mapped to zero
    sigma(map_inds(:,1)) = 0
    
    !set the temperature correction flag
    if(pflotran_temperature_vec_mpi .ne. 0) tcorr_flag = .true.

    do i=1,nmap 
       ei=map_inds(i,1)
       pi=map_inds(i,2)
       a=FMN(ei,1)
       m=FMN(ei,2)
       n=FMN(ei,3)
       delsigb = (pf_tracer(pi)*pf_porosity(pi)**m)*(pf_saturation(pi)**n)/a
    
       if(tcorr_flag) delsigb = delsigb*(1+FMN(ei,4)*(pf_temperature(pi)-25))
       sigma(ei)= sigma(ei) + map(i)*delsigb
     
    end do
  
    do i=1,nelem
       if (sigma(i)<1e-6) sigma(i)=1e-6
    end do
 
    !write(*,*) pf_time
    write(word,'(F25.10)') (pf_time)
    filename = 'sigma_' // &
               trim(adjustl(pflotran_group_prefix)) // &
               '_' // &
               trim(adjustl(word)) // &
               '.txt' 
    !write(*,*) filename
    open(unit=86,file=trim(filename),status='replace',action='write')
    write(86,*) nelem, "1", minval(sigma),maxval(sigma)
    do i = 1, nelem
       write(86,*) sigma(i) 
    enddo
    close(86)

  end subroutine map_pf_e4d
  !____________________________________________________________________


  !____________________________________________________________________
  subroutine send_sigma
    implicit none
    
    call send_command(4)
    call MPI_BCAST(sigma, nelem,MPI_REAL,0,E4D_COMM,ierr)
    
  end subroutine send_sigma
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_dpred
    implicit none
    integer :: opt
    integer :: i,j
    integer :: nadd
    integer :: nbuff
    integer, dimension(nm*2) :: ibuff
    real, dimension(nm) :: rbuff
    integer ::  status(MPI_STATUS_SIZE)
    character*40 :: filename, word 


    !instruct slave to assemble and send the prediceted data
    call send_command(8)
    
    !allocate dpred if not already done and zero
    if (.not. allocated(dpred)) then
       allocate(dpred(nm))
    end if
    dpred = 0
    rbuff = 0
    ibuff = 1
    do i=1,n_rank-1 
       call MPI_RECV(nbuff,1,MPI_INTEGER,i,0,E4D_COMM,status,ierr)
       call MPI_RECV(ibuff(1:nbuff),nbuff,MPI_INTEGER,i,0,E4D_COMM,status,ierr)
       call MPI_RECV(rbuff(1:nbuff),nbuff,MPI_REAL,i,0,E4D_COMM,status,ierr)
       
       do j=1,nbuff
          dpred(ibuff(j))=dpred(ibuff(j))+rbuff(j)
       end do
       
    end do

      
    
    write(word,'(F25.10)') (pf_time)
    filename = 'e4d_' // &
               trim(adjustl(pflotran_group_prefix)) // &
               '_' // &
               trim(adjustl(word)) // &
               '.dpd' 
    write(*,*) filename
    open(unit=86,file=trim(filename),status='replace',action='write')
    write(86,*) nm
    do i = 1, nm
       write(86,'(I6,4I5,4E20.8)') i,s_conf(i,1:4),dpred(i),dobs(i),dpred(i)/sd(i),dobs(i)/sd(i)
    enddo
    close(86)
    
    !!output the predicted data to file
    !call output_dpred
  end subroutine get_dpred
  !____________________________________________________________________
  

  !______________________________________________________________
  subroutine send_command(com)
    !!Send a general command to the slaves
    integer :: com
    call MPI_BCAST(com,1,MPI_INTEGER,0,E4D_COMM,ierr)
  end subroutine send_command
  !________________________________________________________________
  



  !___________SLAVE SUBROUTINES________________________________________

  !____________________________________________________________________
  subroutine slave_run
    implicit none
    integer :: command
    
100 continue
    !Recieve a command from master
    call MPI_BCAST(command,1,MPI_INTEGER,0,E4D_COMM,ierr)
    
    !return to main
    if (command == 0) then
       return

    else if (command == 3) then     
       call build_A
       goto 100
         
    else if (command == 4) then
       call receive_sigma
       goto 100

    else if (command == 5) then
       call build_ksp
       goto 100

    else if (command == 6) then
       call forward_run
       goto 100
       
    else if (command == 8) then
       call send_dpred 
       goto 100
       
    else 
       goto 100

    end if

  end subroutine slave_run
  !____________________________________________________________________
  !__________________________________________________________________
  subroutine receive_sigma

    implicit none
    character(len=32) :: filename, word
    integer :: i
    integer, save :: num_calls = 0
    call MPI_BCAST(sigma, nelem,MPI_REAL,0,E4D_COMM,ierr)   

  end subroutine receive_sigma
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine build_A
    implicit none
    
    integer, dimension(nnodes) :: ncolss
    integer, dimension(50) :: colss
    
    integer :: i,lrnl,lrnu,row,col,rbv,cbv,j,ifn
    logical :: ilow,iup
    
    integer :: rw   
    integer :: ncls 
    integer :: cls(50) 
    PetscReal :: val(1)
    
    
    !zero A
    call MatZeroEntries(A,perr);CHKERRQ(perr)
    
    do i=1,10*nelem
       row=rows(A_map(i))
       col=cols(A_map(i))
       rbv=nbounds(row)
       cbv=nbounds(col)
       
       !lower triangle
       if (((rbv>=2 .and. rbv<=6) .or. (cbv>=2 .and. cbv<=6)) .and. .not. tank_flag ) then
          !one or both nodes are on a boundary so set to zero for bc's
          val(1) = 1e-30
       else
          val(1) = sigma(S_map(i))*delA(i)
       end if
       
       prn(1) = row-1
       pcn(1) = col-1
       
       call MatSetValues(A,1,prn,1,pcn,val,ADD_VALUES,perr);CHKERRQ(perr)
       !upper triangle
       if (row .ne. col) then
          call MatSetValues(A,1,pcn,1,prn,val,ADD_VALUES,perr);CHKERRQ(perr)
       end if
       
    end do
    
    !!fill in the diagonal for the zero potential bc's
    do i=1,nnodes
       !I own this node
       
       if((nbounds(i).eq.2) .and. .not. tank_flag) then
          !this node is on a boundary, set val to 1.0 for bc
          prn(1)=i-1;
          pcn(1)=i-1;
          val(1)=1;
          call MatSetValues(A,1,prn,1,pcn,val,ADD_VALUES,perr)
          
       end if
    end do
    
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,perr);CHKERRQ(perr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,perr);CHKERRQ(perr)
    
  end subroutine build_A
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine build_ksp
    implicit none
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscksp.h90"
    real*8 :: rtol = 1e-6
    real*8 :: atol = 1e-35
    real*8 :: dtol = 500
    integer :: maxints = 10000
    !KSPSetTolerances(KSP ksp,double rtol,double atol,double dtol,int maxits);
    
    !Set up the KSP context
    call KSPCreate(PETSC_COMM_SELF,KS,perr);CHKERRQ(perr)
!    call KSPSetOperators(KS,A,A,SAME_PRECONDITIONER,perr)
    call KSPSetOperators(KS,A,A,perr);CHKERRQ(perr)
    call KSPGetPC(KS,P,perr);CHKERRQ(perr)
    !call KSPSetType(KS,KSPGMRES,perr) !use default
    !call KSPGMRESSetRestart(KS,1000,perr);
    !call KSPGetTolerances(KS,rtol,atol,dtol,maxints,perr);CHKERRQ(perr)
    call KSPSetTolerances(KS,rtol,atol,dtol,maxints,perr);CHKERRQ(perr)
    call KSPSetFromOptions(KS,perr);CHKERRQ(perr)
    
  end subroutine build_ksp
  !__________________________________________________________________
 !____________________________________________________________________
  subroutine forward_run
    implicit none
    integer :: i,m,n,niter,j,enum
    integer :: eindx
    real, dimension(2) :: pck
    PetscScalar :: val
    real :: tstart, tend
    
    do i=1,my_ne
       !call cpu_time(tstart)
       
       call VecGetArrayF90(psol,vloc,perr);CHKERRQ(perr)
       vloc(1:nnodes)=dble(poles(:,i))
       call VecRestoreArrayF90(psol,vloc,perr);CHKERRQ(perr)
       enum=eind(my_rank,1)+i-1
          
       val=0.0
       call VecSet(B,val,perr);CHKERRQ(perr)
       
       !if (i_flg) then
       !   call Add_Jpp(i)
       !end if
       
       eindx=e_nods(enum)
       val=1.0
       call VecSetValues(B,1,eindx-1,val,ADD_VALUES,perr);CHKERRQ(perr)

       if (tank_flag) then
          call VecSetValues(B,1,i_zpot-1,-val,ADD_VALUES,perr);CHKERRQ(perr)
       end if
       call VecAssemblyBegin(B,perr);CHKERRQ(perr)
       call VecAssemblyEnd(B,perr);CHKERRQ(perr)
       
       call KSPSolve(KS,B,psol,perr);CHKERRQ(perr)
       !call KSPView(KS,PETSC_VIEWER_STDOUT_SELF,perr);CHKERRQ(perr)

       call VecGetArrayF90(psol,vloc,perr);CHKERRQ(perr)
       poles(:,i)= real(vloc(1:nnodes))
     
       call VecRestoreArrayF90(psol,vloc,perr);CHKERRQ(perr)
      
       call KSPGetIterationNumber(KS,niter,perr);CHKERRQ(perr)

       !write(*,*) my_rank,niter,maxval(poles(:,i)),minval(poles(:,i))
       !call cpu_time(tend)
       !pck(1)=tend-tstart
       !pck(2)=real(niter)
       !call MPI_SEND(pck,2,MPI_REAL,0,1,MPI_COMM_WORLD,perr)
       !write(*,*) "Slave ",my_rank," solved for pole ",eind(my_rank,1)+i-1,'in ',tend-tstart,' seconds and ',niter,' iters'
       
    end do
    
    if (first_sol) then
       call KSPSetInitialGuessNonzero(KS,PETSC_TRUE,perr);CHKERRQ(perr)
       first_sol=.false.
    end if
    
    call KSPDestroy(KS,perr);CHKERRQ(perr)
    
    
  end subroutine forward_run
!_________________________________________________________________________________________________


   !__________________________________________________________________
    subroutine send_dpred
      implicit none
      integer :: flg,i
 
      call assemble_data
      call MPI_SEND(nmy_drows,1,MPI_INTEGER,0,0,E4D_COMM, ierr)
      call MPI_SEND(my_drows,nmy_drows,MPI_INTEGER,0,0,E4D_COMM,ierr)
      call MPI_SEND(my_dvals,nmy_drows,MPI_REAL,0,0,E4D_COMM,ierr)
       
    end subroutine send_dpred
    !__________________________________________________________________
     
    !___________________________________________________________
    subroutine assemble_data
      implicit none
      integer :: i,a,b,m,n,e1,e2,indx,p
      
      
      e1=eind(my_rank,1)
      e2=eind(my_rank,2) 
      if (.not.allocated(my_drows)) then
         nmy_drows=0
         do i=1,nm
            a=s_conf(i,1)
            b=s_conf(i,2)
            if ((a>=e1 .and. a<=e2) .or. (b>=e1.and. b<=e2)) then
               nmy_drows=nmy_drows+1
            end if
         end do
         allocate(my_drows(nmy_drows),my_dvals(nmy_drows))
      end if
      
      indx=0
      my_dvals=0
      
     
      do i=1,nm
         a = s_conf(i,1)
         b = s_conf(i,2)
         m = s_conf(i,3)
         n = s_conf(i,4)
         
         if ((a>=e1 .and. a<=e2) .or. (b>=e1.and. b<=e2)) then
            indx=indx+1
            my_drows(indx)=i
            
            do p=e1,e2
               if (p==a) then
                  if (m.ne.0) my_dvals(indx)=my_dvals(indx) + real(poles(e_nods(m),a-e1+1))
                  if (n.ne.0) my_dvals(indx)=my_dvals(indx) - real(poles(e_nods(n),a-e1+1))
               end if
               if (p==b) then
                  if (m.ne.0) my_dvals(indx)=my_dvals(indx) - real(poles(e_nods(m),b-e1+1))
                  if (n.ne.0) my_dvals(indx)=my_dvals(indx) + real(poles(e_nods(n),b-e1+1))
               end if
            end do
         end if
      end do
     
  end subroutine assemble_data
  !___________________________________________________________  
  !___________________________________________________________
  subroutine cleanup
  
        if(allocated(e4d_ranks)) deallocate(e4d_ranks)
        if(allocated(pf_e4d_ranks)) deallocate(pf_e4d_ranks)
        if(allocated(map_inds)) deallocate(map_inds)
        if(allocated(s_conf)) deallocate(s_conf)
        if(allocated(eind)) deallocate(eind)
        if(allocated(jind)) deallocate(jind)
        if(allocated(nbounds)) deallocate(nbounds)
        if(allocated(zones)) deallocate(zones)
        if(allocated(elements)) deallocate(elements)
        if(allocated(faces)) deallocate(faces)
        if(allocated(e_nods)) deallocate(e_nods)
        if(allocated(rows)) deallocate(rows)
        if(allocated(cols)) deallocate(cols)
        if(allocated(trows)) deallocate(trows)
        if(allocated(tcols)) deallocate(tcols)
        if(allocated(A_map)) deallocate(A_map)
        if(allocated(S_map)) deallocate(S_map)
        if(allocated(my_drows)) deallocate(my_drows)
        if(allocated(e_pos)) deallocate(e_pos)
        if(allocated(nodes)) deallocate(nodes)
        if(allocated(poles)) deallocate(poles)
        if(allocated(pf_tracer)) deallocate(pf_tracer)
        if(allocated(pf_saturation)) deallocate(pf_saturation)
        if(allocated(pf_saturation_0)) deallocate(pf_saturation_0)
        if(allocated(sigma)) deallocate(sigma)
        if(allocated(dpred)) deallocate(dpred)
        if(allocated(dobs)) deallocate(dobs)
        if(allocated(sd)) deallocate(sd)
        if(allocated(my_dvals)) deallocate(my_dvals)
        if(allocated(map)) deallocate(map)
        if(allocated(base_sigma)) deallocate(base_sigma)
        if(allocated(ffac)) deallocate(ffac)
        if(allocated(pfxcb)) deallocate(pfxcb)
        if(allocated(pfycb)) deallocate(pfycb)
        if(allocated(pfzcb)) deallocate(pfzcb)
        if(allocated(d_nnz)) deallocate(d_nnz)
        if(allocated(delA)) deallocate(delA)
  end subroutine cleanup
  !___________________________________________________________
end module e4d_run
