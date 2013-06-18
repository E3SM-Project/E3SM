! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_todynamics
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_constants, only: R_d,R_v,degrad

 implicit none
 private
 public:: physics_addtend


!>\brief interface between the physics parameterizations and the non-hydrostatic dynamical core.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_todynamics:
!> ---------------------------------------
!> physics_addtend: add and mass-weigh tendencies before being added to dynamics tendencies.
!> tend_toEdges   : interpolate wind-tendencies from centers to edges of grid-cells.


 contains
 
!==================================================================================================
 subroutine physics_addtend(mesh, state, diag, tend, tend_physics, mass, mass_edge, rk_step)
!==================================================================================================

!input variables:
!----------------
 type(mesh_type),intent(in):: mesh
 type(state_type),intent(in):: state
 type(diag_type),intent(in):: diag
 type(tend_physics_type),intent(inout):: tend_physics
 integer, intent(in):: rk_step
 real(kind=RKIND),dimension(:,:),intent(in):: mass
 real(kind=RKIND),dimension(:,:),intent(in):: mass_edge

!inout variables:
!----------------
 type(tend_type),intent(inout):: tend

!local variables:
!----------------

 type(block_type),pointer :: block

 integer:: i,iCell,k,n,nCells,nCellsSolve,nEdges,nEdgesSolve,nVertLevels

 real(kind=RKIND),dimension(:,:),pointer:: theta_m,qv
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rqvblten,rqcblten, &
                                           rqiblten,rublten,rvblten
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten, &
                                           rqrcuten,rqicuten,rqscuten, &
                                           rucuten,rvcuten
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw                                    
 
 real(kind=RKIND),dimension(:,:),pointer  :: tend_theta,tend_u
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_scalars

 real(kind=RKIND):: tem
 real(kind=RKIND),dimension(:,:),allocatable:: rublten_Edge,rucuten_Edge

!ldf (2011-12-16):
 real(kind=RKIND),dimension(:,:),allocatable:: theta,tend_th
!ldf end.

!==================================================================================================
 block => mesh % block
 
 nCells      = mesh % nCells
 nEdges      = mesh % nEdges
 nCellsSolve = mesh % nCellsSolve
 nEdgesSolve = mesh % nEdgesSolve
 nVertLevels = mesh % nVertLevels

!theta   => diag % theta % array
 theta_m => state % theta_m % array
 qv      => state % scalars % array(state%index_qv,:,:)

 rublten    => tend_physics % rublten    % array
 rvblten    => tend_physics % rvblten    % array
 rthblten   => tend_physics % rthblten   % array
 rqvblten   => tend_physics % rqvblten   % array
 rqcblten   => tend_physics % rqcblten   % array
 rqiblten   => tend_physics % rqiblten   % array

 rucuten    => tend_physics % rucuten    % array
 rvcuten    => tend_physics % rvcuten    % array
 rthcuten   => tend_physics % rthcuten   % array
 rqvcuten   => tend_physics % rqvcuten   % array
 rqccuten   => tend_physics % rqccuten   % array
 rqrcuten   => tend_physics % rqrcuten   % array
 rqicuten   => tend_physics % rqicuten   % array
 rqscuten   => tend_physics % rqscuten   % array

 rthratenlw => tend_physics % rthratenlw % array
 rthratensw => tend_physics % rthratensw % array

 tend_u       => tend % u % array
 tend_theta   => tend % theta_m % array
 tend_scalars => tend % scalars % array

!initialize the tendency for the potential temperature and all scalars due to PBL, convection,
!and longwave and shortwave radiation:
 allocate(theta(nVertLevels,nCellsSolve)  )
 allocate(tend_th(nVertLevels,nCellsSolve))
 tend_th      = 0.
 tend_scalars = 0.

!add coupled tendencies due to PBL processes:
 if(config_pbl_scheme .ne. 'off') then
    allocate(rublten_Edge(nVertLevels,nEdges))
    rublten_Edge(:,:) = 0.
    call tend_toEdges(mesh,rublten,rvblten,rublten_Edge)
    do i = 1, nEdgesSolve
    do k  = 1, nVertLevels
       tend_u(k,i)=tend_u(k,i)+rublten_Edge(k,i)*mass_edge(k,i)
    enddo
    enddo
    deallocate(rublten_Edge)

    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_th(k,i) = tend_th(k,i)+rthblten(k,i)*mass(k,i)
       tend_scalars(tend%index_qv,k,i)=tend_scalars(tend%index_qv,k,i)+rqvblten(k,i)*mass(k,i)
       tend_scalars(tend%index_qc,k,i)=tend_scalars(tend%index_qc,k,i)+rqcblten(k,i)*mass(k,i)
       tend_scalars(tend%index_qi,k,i)=tend_scalars(tend%index_qi,k,i)+rqiblten(k,i)*mass(k,i)
    enddo
    enddo
 endif

!add coupled tendencies due to convection:
 if(config_conv_deep_scheme .ne. 'off') then

    if(config_conv_deep_scheme .eq. 'tiedtke') then
       allocate(rucuten_Edge(nVertLevels,nEdges))
       rucuten_Edge(:,:) = 0.
       call tend_toEdges(mesh,rucuten,rvcuten,rucuten_Edge)
       do i = 1, nEdgesSolve
       do k  = 1, nVertLevels
          tend_u(k,i)=tend_u(k,i)+rucuten_Edge(k,i)*mass_edge(k,i)
       enddo
       enddo
       deallocate(rucuten_Edge)
    endif

    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_th(k,i)=tend_th(k,i)+rthcuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qv,k,i)=tend_scalars(tend%index_qv,k,i)+rqvcuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qc,k,i)=tend_scalars(tend%index_qc,k,i)+rqccuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qr,k,i)=tend_scalars(tend%index_qr,k,i)+rqrcuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qi,k,i)=tend_scalars(tend%index_qi,k,i)+rqicuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qs,k,i)=tend_scalars(tend%index_qs,k,i)+rqscuten(k,i)*mass(k,i)    
    enddo
    enddo
 endif

!add coupled tendencies due to longwave radiation:
 if(config_radt_lw_scheme .ne. 'off') then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_th(k,i)=tend_th(k,i)+rthratenlw(k,i)*mass(k,i)
    enddo
    enddo
 endif
    
!add coupled tendencies due to shortwave radiation:
 if(config_radt_sw_scheme .ne. 'off') then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_th(k,i)=tend_th(k,i)+rthratensw(k,i)*mass(k,i)
    enddo
    enddo
 endif

!if non-hydrostatic core, convert the tendency for the potential temperature to a
!tendency for the modified potential temperature:
 do i = 1, nCellsSolve
 do k = 1, nVertLevels
    theta(k,i) = theta_m(k,i) / (1. + R_v/R_d * qv(k,i))
    tend_th(k,i) = (1. + R_v/R_d * qv(k,i)) * tend_th(k,i) &
                    + R_v/R_d * theta(k,i) * tend_scalars(tend%index_qv,k,i)
    tend_theta(k,i) = tend_theta(k,i) + tend_th(k,i)
 enddo
 enddo
 deallocate(theta)
 deallocate(tend_th)

 if(rk_step .eq. 3) then
    write(0,*)
    write(0,*) '--- enter subroutine physics_addtend:'
    write(0,*) 'max rthblten   = ',maxval(rthblten(:,1:nCellsSolve))
    write(0,*) 'min rthblten   = ',minval(rthblten(:,1:nCellsSolve))
    write(0,*) 'max rthcuten   = ',maxval(rthcuten(:,1:nCellsSolve))
    write(0,*) 'min rthcuten   = ',minval(rthcuten(:,1:nCellsSolve))
    write(0,*) 'max rthratenlw = ',maxval(rthratenlw(:,1:nCellsSolve))
    write(0,*) 'min rthratenlw = ',minval(rthratenlw(:,1:nCellsSolve))
    write(0,*) 'max rthratensw = ',maxval(rthratensw(:,1:nCellsSolve))
    write(0,*) 'min rthratensw = ',minval(rthratensw(:,1:nCellsSolve))
    write(0,*) '--- end subroutine physics_addtend'
    write(0,*)
 endif

!formats:
 201 format(2i6,10(1x,e15.8))
 202 format(3i6,10(1x,e15.8))

 end subroutine  physics_addtend

!==================================================================================================
 subroutine tend_toEdges(mesh,Ux_tend,Uy_tend,U_tend)
!==================================================================================================

!input arguments:
!----------------
 type(mesh_type),intent(in):: mesh
 real(kind=RKIND),intent(in),dimension(:,:):: Ux_tend,Uy_tend 

!output arguments:
!-----------------
 real(kind=RKIND),intent(out),dimension(:,:):: U_tend

!local variables:
!-----------------
 type(block_type),pointer :: block
 type (field2DReal), pointer :: tempField
 type (field2DReal), target :: tempFieldTarget
 integer:: iCell,iEdge,k,j,nCells,nCellsSolve,nVertLevels
 integer,dimension(:),pointer  :: nEdgesOnCell
 integer,dimension(:,:),pointer:: edgesOnCell

 real(kind=RKIND),dimension(:,:),pointer:: east,north,edge_normal
 real(kind=RKIND),dimension(:,:),allocatable,target:: Ux_tend_halo,Uy_tend_halo
 
!--------------------------------------------------------------------------------------------------

 block => mesh % block
 
 nCells       = mesh % nCells
 nCellsSolve  = mesh % nCellsSolve
 nVertLevels  = mesh % nVertLevels

 east  => mesh % east  % array
 north => mesh % north % array
 edgesOnCell  => mesh % edgesOnCell  % array
 nEdgesOnCell => mesh % nEdgesOnCell % array
 edge_normal  => mesh % edgeNormalVectors % array

 allocate(Ux_tend_halo(nVertLevels,nCells+1))
 allocate(Uy_tend_halo(nVertLevels,nCells+1))

 Ux_tend_halo(:,:) = 0.
 Uy_tend_halo(:,:) = 0.
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       Ux_tend_halo(k,iCell) = Ux_tend(k,iCell)
       Uy_tend_halo(k,iCell) = Uy_tend(k,iCell)
    enddo
 enddo

 tempField => tempFieldTarget 
 tempField % block => block
 tempField % dimSizes(1) = nVertLevels
 tempField % dimSizes(2) = nCellsSolve
 tempField % sendList => block % parinfo % cellsToSend
 tempField % recvList => block % parinfo % cellsToRecv
 tempField % copyList => block % parinfo % cellsToCopy
 tempField % prev => null()
 tempField % next => null()

 tempField % array => Ux_tend_halo
 call mpas_dmpar_exch_halo_field(tempField)
 
 tempField % array => Uy_tend_halo
 call mpas_dmpar_exch_halo_field(tempField)

 U_tend(:,:) = 0.0
 do iCell = 1, nCells
 do j = 1, nEdgesOnCell(iCell)
    iEdge = edgesOnCell(j,iCell)
    do k = 1, nVertLevels
       U_tend(k,iEdge) = U_tend(k,iEdge) & 
                + 0.5 * Ux_tend_halo(k,iCell) * (edge_normal(1,iEdge) * east(1,iCell)   &
                                              +  edge_normal(2,iEdge) * east(2,iCell)   &
                                              +  edge_normal(3,iEdge) * east(3,iCell))  &
                + 0.5 * Uy_tend_halo(k,iCell) * (edge_normal(1,iEdge) * north(1,iCell)  &
                                              +  edge_normal(2,iEdge) * north(2,iCell)  &
                                              +  edge_normal(3,iEdge) * north(3,iCell))
    enddo
 enddo
 enddo

 deallocate(Ux_tend_halo)
 deallocate(Uy_tend_halo)
 
 end subroutine tend_toEdges

!==================================================================================================
 end module mpas_atmphys_todynamics
!==================================================================================================
