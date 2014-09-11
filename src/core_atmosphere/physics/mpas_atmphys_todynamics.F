! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_todynamics
 use mpas_kind_types
 use mpas_grid_types
 use mpas_dmpar

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
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-22.
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-05-15.


 contains

 
!==================================================================================================
 subroutine physics_addtend(block, mesh, state, diag, tend, tend_physics, configs, rk_step)
!==================================================================================================

!input variables:
!----------------
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(inout):: tend_physics
 type(mpas_pool_type),intent(in):: configs
 integer, intent(in):: rk_step

!inout variables:
!----------------
 type(mpas_pool_type),intent(inout):: tend

!local variables:
!----------------


 integer:: i,iCell,k,n
 integer,pointer:: index_qv, index_qc, index_qr, index_qi, index_qs, index_qg
 integer,pointer:: nCells,nCellsSolve,nEdges,nEdgesSolve,nVertLevels

 real(kind=RKIND),dimension(:,:),pointer:: mass          ! time level 2 rho_zz
 real(kind=RKIND),dimension(:,:),pointer:: mass_edge     ! diag rho_edge
 real(kind=RKIND),dimension(:,:),pointer:: theta_m,qv    ! time level 1
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rqvblten,rqcblten, &
                                           rqiblten,rublten,rvblten
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten, &
                                           rqrcuten,rqicuten,rqscuten, &
                                           rucuten,rvcuten
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw                                    
 
 real(kind=RKIND),dimension(:,:),pointer:: tend_theta,tend_u
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_scalars
 real(kind=RKIND),dimension(:,:),pointer:: tend_qv,tend_qc,tend_qr,tend_qi,tend_qs,tend_qg

 real(kind=RKIND):: tem
 real(kind=RKIND),dimension(:,:),allocatable:: rublten_Edge,rucuten_Edge

 character(len=StrKIND), pointer :: config_pbl_scheme, config_conv_deep_scheme, config_radt_lw_scheme, config_radt_sw_scheme

!ldf (2011-12-16):
 real(kind=RKIND),dimension(:,:),allocatable:: theta,tend_th
!ldf end.

!==================================================================================================
 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
 call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
 call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

 call mpas_pool_get_config(configs, 'config_pbl_scheme', config_pbl_scheme)
 call mpas_pool_get_config(configs, 'config_conv_deep_scheme', config_conv_deep_scheme)
 call mpas_pool_get_config(configs, 'config_radt_lw_scheme', config_radt_lw_scheme)
 call mpas_pool_get_config(configs, 'config_radt_sw_scheme', config_radt_sw_scheme)

 call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
 call mpas_pool_get_array(state, 'scalars', scalars, 1)
 call mpas_pool_get_array(state, 'rho_zz', mass, 2)
 call mpas_pool_get_array(diag , 'rho_edge', mass_edge)

 call mpas_pool_get_dimension(state, 'index_qv', index_qv)
 call mpas_pool_get_dimension(state, 'index_qc', index_qc)
 call mpas_pool_get_dimension(state, 'index_qr', index_qr)
 call mpas_pool_get_dimension(state, 'index_qi', index_qi)
 call mpas_pool_get_dimension(state, 'index_qs', index_qs)
 call mpas_pool_get_dimension(state, 'index_qg', index_qg)
 qv => scalars(index_qv,:,:)

 call mpas_pool_get_array(tend_physics, 'rublten', rublten)
 call mpas_pool_get_array(tend_physics, 'rvblten', rvblten)
 call mpas_pool_get_array(tend_physics, 'rthblten', rthblten)
 call mpas_pool_get_array(tend_physics, 'rqvblten', rqvblten)
 call mpas_pool_get_array(tend_physics, 'rqcblten', rqcblten)
 call mpas_pool_get_array(tend_physics, 'rqiblten', rqiblten)

 call mpas_pool_get_array(tend_physics, 'rucuten', rucuten)
 call mpas_pool_get_array(tend_physics, 'rvcuten', rvcuten)
 call mpas_pool_get_array(tend_physics, 'rthcuten', rthcuten)
 call mpas_pool_get_array(tend_physics, 'rqvcuten', rqvcuten)
 call mpas_pool_get_array(tend_physics, 'rqccuten', rqccuten)
 call mpas_pool_get_array(tend_physics, 'rqrcuten', rqrcuten)
 call mpas_pool_get_array(tend_physics, 'rqicuten', rqicuten)
 call mpas_pool_get_array(tend_physics, 'rqscuten', rqscuten)

 call mpas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
 call mpas_pool_get_array(tend_physics, 'rthratensw', rthratensw)

 call mpas_pool_get_array(tend, 'u', tend_u)
 call mpas_pool_get_array(tend, 'theta_m', tend_theta)
 call mpas_pool_get_array(tend, 'scalars_tend', tend_scalars)
 tend_qv => tend_scalars(index_qv,:,:)
 tend_qc => tend_scalars(index_qc,:,:)
 tend_qr => tend_scalars(index_qr,:,:)
 tend_qi => tend_scalars(index_qi,:,:)
 tend_qs => tend_scalars(index_qs,:,:)
 tend_qg => tend_scalars(index_qg,:,:)

!initialize the tendency for the potential temperature and all scalars due to PBL, convection,
!and longwave and shortwave radiation:
 allocate(theta(nVertLevels,nCellsSolve)  )
 allocate(tend_th(nVertLevels,nCellsSolve))
 tend_th = 0._RKIND
 tend_qv = 0._RKIND
 tend_qc = 0._RKIND
 tend_qr = 0._RKIND
 tend_qi = 0._RKIND
 tend_qs = 0._RKIND
 tend_qg = 0._RKIND

!add coupled tendencies due to PBL processes:
 if(config_pbl_scheme .ne. 'off') then
    allocate(rublten_Edge(nVertLevels,nEdges))
    rublten_Edge(:,:) = 0.
    call tend_toEdges(block,mesh,rublten,rvblten,rublten_Edge)
    do i = 1, nEdgesSolve
    do k  = 1, nVertLevels
       tend_u(k,i)=tend_u(k,i)+rublten_Edge(k,i)*mass_edge(k,i)
    enddo
    enddo
    deallocate(rublten_Edge)

    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_th(k,i) = tend_th(k,i) + rthblten(k,i)*mass(k,i)
       tend_qv(k,i) = tend_qv(k,i) + rqvblten(k,i)*mass(k,i)
       tend_qc(k,i) = tend_qc(k,i) + rqcblten(k,i)*mass(k,i)
       tend_qi(k,i) = tend_qi(k,i) + rqiblten(k,i)*mass(k,i)
    enddo
    enddo
 endif

!add coupled tendencies due to convection:
 if(config_conv_deep_scheme .ne. 'off') then

    if(config_conv_deep_scheme .eq. 'tiedtke') then
       allocate(rucuten_Edge(nVertLevels,nEdges))
       rucuten_Edge(:,:) = 0.
       call tend_toEdges(block,mesh,rucuten,rvcuten,rucuten_Edge)
       do i = 1, nEdgesSolve
       do k  = 1, nVertLevels
          tend_u(k,i)=tend_u(k,i)+rucuten_Edge(k,i)*mass_edge(k,i)
       enddo
       enddo
       deallocate(rucuten_Edge)
    endif

    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_th(k,i) = tend_th(k,i) + rthcuten(k,i)*mass(k,i)
       tend_qv(k,i) = tend_qv(k,i) + rqvcuten(k,i)*mass(k,i)
       tend_qc(k,i) = tend_qc(k,i) + rqccuten(k,i)*mass(k,i)
       tend_qr(k,i) = tend_qr(k,i) + rqrcuten(k,i)*mass(k,i)
       tend_qi(k,i) = tend_qi(k,i) + rqicuten(k,i)*mass(k,i)
       tend_qs(k,i) = tend_qs(k,i) + rqscuten(k,i)*mass(k,i)    
    enddo
    enddo
 endif

!add coupled tendencies due to longwave radiation:
 if(config_radt_lw_scheme .ne. 'off') then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_th(k,i) = tend_th(k,i) + rthratenlw(k,i)*mass(k,i)
    enddo
    enddo
 endif
    
!add coupled tendencies due to shortwave radiation:
 if(config_radt_sw_scheme .ne. 'off') then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_th(k,i) = tend_th(k,i) + rthratensw(k,i)*mass(k,i)
    enddo
    enddo
 endif

!if non-hydrostatic core, convert the tendency for the potential temperature to a
!tendency for the modified potential temperature:
 do i = 1, nCellsSolve
 do k = 1, nVertLevels
    theta(k,i) = theta_m(k,i) / (1. + R_v/R_d * qv(k,i))
    tend_th(k,i) = (1. + R_v/R_d * qv(k,i)) * tend_th(k,i) &
                    + R_v/R_d * theta(k,i) * tend_qv(k,i)
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
 subroutine tend_toEdges(block,mesh,Ux_tend,Uy_tend,U_tend)
!==================================================================================================

!input arguments:
!----------------
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: mesh
 real(kind=RKIND),intent(in),dimension(:,:):: Ux_tend,Uy_tend 

!output arguments:
!-----------------
 real(kind=RKIND),intent(out),dimension(:,:):: U_tend

!local variables:
!-----------------
 type (field2DReal), pointer :: tempField
 type (field2DReal), target :: tempFieldTarget
 integer:: iCell,iEdge,k,j
 integer,pointer:: nCells,nCellsSolve,nVertLevels
 integer,dimension(:),pointer  :: nEdgesOnCell
 integer,dimension(:,:),pointer:: edgesOnCell

 real(kind=RKIND),dimension(:,:),pointer:: east,north,edge_normal
 real(kind=RKIND),dimension(:,:),allocatable,target:: Ux_tend_halo,Uy_tend_halo
 
!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

 call mpas_pool_get_array(mesh, 'east', east)
 call mpas_pool_get_array(mesh, 'north', north)
 call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'edgeNormalVectors', edge_normal)

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
