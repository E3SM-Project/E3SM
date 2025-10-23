!
MODULE MOSART_lake_geometry_mod
! Description: calculate lake geometric relationships
! Revised by Hongyi Li, July 2024
!-----------------------------------------------------------------------
    
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use shr_sys_mod   , only : shr_sys_abort
    use RunoffMod     , only : Tctl, TUnit, TRunoff, THeat, TUnit_lake_r, TLake_r, TUnit_lake_t, TLake_t, TPara, rtmCTL
    use RtmVar        , only : iulog, ngeom, nlayers, rstraflag, lakeflag

    implicit none
    private

    real(r8), parameter :: TINYVALUE = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits

! !PUBLIC MEMBER FUNCTIONS:

    public lakegeom_init_r
    public lakegeom_update_r
    public toplayer_merge_r
    public layer_merge_r
    public layer_split_r
	public lake_h2v

    contains  

!----------------------------------------------------------------------------
    subroutine lakegeom_init_r
       
      ! Calculate reservoir layer average area (km2) 
      use shr_sys_mod   , only : shr_sys_flush
      use RunoffMod     , only : rtmCTL
      use RtmVar         , only : iulog, ngeom, nlayers
       
      implicit none
      real(r8) :: M_W,M_L,gm_j,d_res,dd_in
      real(r8),dimension(ngeom+1) :: d_zi0,v_zti0,a_di0
      real(r8) ::dd_zz(ngeom),a_dd(ngeom+1),a_zi(ngeom+1),C_aa(ngeom+1), ar_f = 1.0e6                 ! Factor to convert area to m^2
      real(r8) :: pi      = 3.1415926_r8      ! pi
      real(r8) :: ddz_top    = 3._r8  !2._r8        !0.60_r8                                                   ! Surface layer depth (m)
      real(r8) :: ddz_min    = 5._r8        !minimum layer depth allowed to ensure numerical stability (m)
      real(r8) :: dv,da,dz
      real(r8) :: d_v(nlayers)                           ! Reservoir volume change at layer (m^3)
      real(r8) :: rho_z(nlayers)                           ! Reservoir layer density (kg/m^3), taken constant for furture revision 
      real(r8) :: delta_z                                ! depth change to calculate corresponding area/volume(m)
      integer :: i,j,k,iunit                            ! indices
      integer :: tmp_d_ns                               ! temporary number of layers
      character(len=*),parameter :: subname = '(lakegeom_init_r)'        
                                 
      !**************************************************************************************************************************************************
         
      do iunit = rtmCTL%begr,rtmCTL%endr    
             
          if (TLake_r%d_ns(iunit) >= 1 .and. TUnit_lake_r%lake_flg(iunit) >=1) then ! 
               !if (TUnit_lake_r%h_lake(iunit) <= 2.0_r8) TUnit_lake_r%h_lake(iunit) = 2._r8
               TLake_r%d_lake(iunit) = TUnit_lake_r%h_lake(iunit)         
               d_res = TLake_r%d_lake(iunit)        
               
               ! Uniform subsurface layer depth for initialization    and limit maximum layer thickness 
               dd_in = d_res/ngeom !bottom layers evenly descritized
                                 
               ! Calculate reservoir geometry to establish depth-area-volume relationship            
                         
               ! Calculate depth area    
               do j = 1, ngeom!    
                   a_dd(j) = lake_h2a(dd_in*(ngeom-j+1), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
               end do
                         
               a_dd(ngeom+1) = a_dd(ngeom)*0.001_r8            !Bottom area given non-zero value
                         
               ! Reverse indexing so that bottom is 1 and top is d_n+1 and convert to m2
               do j = 1,ngeom+1
                    k =ngeom+2-j  
                    a_di0(k) = a_dd(j)
                    TLake_r%a_di(iunit,k)  = a_di0(k)     !Area corrected for error for optimal geometry
               end do    
                     
               ! Calculate layer depth,area,and volume from bottom 
               d_zi0(1) = 0._r8
               TLake_r%d_zi(iunit,1)  = d_zi0(1)
               do j = 2, ngeom+1    
                    d_zi0(j) = d_zi0(j-1) + dd_in
                    TLake_r%d_zi(iunit,j)  = d_zi0(j) 
               end do
                             
               ! Calculate layer average area,and total volume from bottom
               !v_zti0(1) = 0.001_r8 * lake_h2v(0.5_r8*(d_zi0(1)+d_zi0(2)), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_biunit))    
               !TLake_r%v_zti(iunit,1) = v_zti0(1)     !lower Volume
               !if (TLake_r%v_zti(iunit,1)==0._r8) TLake_r%v_zti(iunit,1)=1._r8
               v_zti0(1) = 0._r8
               TLake_r%v_zti(iunit,1) = 0._r8
               do j = 2, ngeom+1     
                    !a_zi(j) = 0.5_r8*(TLake_r%a_di(iunit,j)+TLake_r%a_di(iunit,j-1)) !Area converted to m^2
                    v_zti0(j) = lake_h2v(d_zi0(j), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
                    TLake_r%v_zti(iunit,j) = v_zti0(j)   
                    !if (TLake_r%v_zti(iunit,j)==0._r8) TLake_r%v_zti(iunit,j)=1._r8
                    !if (TLake_r%a_di(iunit,j)<0._r8 .or. TLake_r%v_zti(iunit,j)<0._r8) write(iulog,*)subname,'geom negative',iunit
               end do
               
    
               ! treating shallow or high-latitude lakes as single-layer ones (latter can be relaxed when adding a lake-ice module)
               !if(TLake_r%d_ns(iunit) == 1 .or. TUnit_lake_r%h_lake(iunit) <=10._r8 .or. abs(rtmCTL%latc(iunit)) >=40._r8) then
               if(abs(rtmCTL%latc(iunit)) >=40._r8) then
                   TUnit_lake_r%one_layer(iunit) = 1
               end if
               
               ! adjust total number of layers to avoid too small initial layer depth
               if(TLake_r%d_ns(iunit) >= 2) then
                   TLake_r%ddz_local(iunit) = (TLake_r%d_lake(iunit) - ddz_top) / (TLake_r%d_ns(iunit) - 1)
                   if (TLake_r%ddz_local(iunit) - ddz_min < -TINYVALUE) then
                       tmp_d_ns = floor((TLake_r%d_lake(iunit) - ddz_top)/ddz_min)
                       if(tmp_d_ns <= 1) then
                           TLake_r%d_ns(iunit) = 1
                       else
                           TLake_r%d_ns(iunit) = tmp_d_ns
                       end if
                   end if
                   
                   if(TUnit_lake_r%one_layer(iunit) == 1) then
                       TLake_r%d_ns(iunit) = 1
                   end if
               end if
    
               !Initial reservoir storage 
               do j = TLake_r%d_ns(iunit),1,-1
                    if (j == TLake_r%d_ns(iunit) .and. TLake_r%d_ns(iunit) == 1) then
                          TLake_r%dd_z(iunit,j) = TLake_r%d_lake(iunit)
                    elseif (j == TLake_r%d_ns(iunit) .and. TLake_r%d_ns(iunit) > 1) then !top layer depth kept constant
                          TLake_r%dd_z(iunit,j) = ddz_top        !0.6_r8
                    elseif ((TLake_r%d_ns(iunit)>1 .and.j < TLake_r%d_ns(iunit)).and.(TLake_r%d_lake(iunit) - TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)))>0._r8) then
                          TLake_r%dd_z(iunit,j) = (TLake_r%d_lake(iunit) - TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))) / (TLake_r%d_ns(iunit) - 1) !bottom layers evenly descritized
                    endif
               end do    
                 
               if (TLake_r%d_ns(iunit)>1 .and. TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1)<ddz_top)then !layer thickness too small
                    TLake_r%d_ns(iunit)=int((TLake_r%d_lake(iunit)/ddz_top)+1) 
                    do j=1,nlayers!TLake_r%d_ns(iunit)
                          TLake_r%dd_z(iunit,j) = 0._r8
                    end do
                         
                    ! Reinitialize layer thickness
                    do j = TLake_r%d_ns(iunit),1,-1
                          if (j == TLake_r%d_ns(iunit)) then
                              TLake_r%dd_z(iunit,j) = ddz_top      !top layer depth
                          else
                              TLake_r%dd_z(iunit,j) = (TLake_r%d_lake(iunit) - TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)))/(TLake_r%d_ns(iunit) - 1) !bottom layers evenly descritized
                          end if
                    end do    
               end if
             
               if (TLake_r%d_ns(iunit)>1) then
                    TLake_r%ddz_local(iunit) = TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1)
               else
                    TLake_r%ddz_local(iunit) = TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
               end if
                     
               ! Calculate layer depth (minimum at bottom)    
               TLake_r%d_z(iunit,1)=0._r8
               do j = 2, nlayers+1  
                    if (j<=TLake_r%d_ns(iunit)+1) then
                          TLake_r%d_z(iunit,j) = TLake_r%d_z(iunit,j-1) + TLake_r%dd_z(iunit,j-1)
                    else 
                          TLake_r%d_z(iunit,j) = 0._r8
                    end if
               end do    
               ! Assign layer area and volume based on depth-area-volume relationship
               TLake_r%a_d(iunit,1)=(TLake_r%a_di(iunit,1))
               TLake_r%v_zt(iunit,1)=(TLake_r%v_zti(iunit,1))
               do i=2,nlayers+1
                    if (i<=TLake_r%d_ns(iunit)+1) then
                          do j=2,ngeom+1
                               TLake_r%a_d(iunit,i)  = lake_h2a(TLake_r%d_z(iunit,i), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
                               TLake_r%v_zt(iunit,i) = lake_h2v(TLake_r%d_z(iunit,i), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
                          end do
                    else
                          TLake_r%a_d(iunit,i) = 0._r8
                          TLake_r%v_zt(iunit,i) = 0._r8
                    end if
               end do
               ! Calculate layer volume(m^3)
               do j = 1, nlayers    
                    if (j<=TLake_r%d_ns(iunit)) then
                          TLake_r%d_v(iunit,j) = ((TLake_r%v_zt(iunit,j+1) - TLake_r%v_zt(iunit,j)))    
                          TLake_r%dd_z(iunit,j) = TLake_r%d_z(iunit,j+1) - TLake_r%d_z(iunit,j)
                          if (TLake_r%d_v(iunit,j)<0.0_r8) then 
                              write(iulog,*) subname,'Layer volume negative: Check geometry data for lake',iunit
                          end if
                    else
                          TLake_r%d_v(iunit,j) = 0._r8
                          TLake_r%dd_z(iunit,j) = 0._r8
                    end if
               end do        
                     
               !     Intitialize layer temperature and total storage                    
               do j=1,TLake_r%d_ns(iunit)                    
                    rho_z(j) = 1000._r8*( 1._r8 - 1.9549e-05*(abs(TLake_r%temp_lake(iunit,j)-277._r8))**1.68_r8) 
                    TLake_r%v_zn(iunit,j) = (TLake_r%d_v(iunit,j))
                    TLake_r%m_zo(iunit,j) = (TLake_r%d_v(iunit,j)*rho_z(j)) 
                    TLake_r%m_zn(iunit,j) = (TLake_r%d_v(iunit,j)*rho_z(j))
               end do
               
               TUnit_lake_r%h_min(iunit) = TLake_r%d_zi(iunit,ngeom+1) - TUnit%rdepth(iunit)
               if(TUnit_lake_r%h_min(iunit) < TLake_r%d_zi(iunit,1)) then
                   TUnit_lake_r%h_min(iunit) = TLake_r%d_zi(iunit,1)
               end if
               TUnit_lake_r%v_min(iunit) = lake_h2v(TUnit_lake_r%h_min(iunit), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
               
               TLake_r%V_str(iunit) = TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) +1)          
          end if
      end do    
    end subroutine lakegeom_init_r

    subroutine lakegeom_update_r(iunit)
!*******************************************************************************************************
!     Calculate new layer thickness, depth, area/volume after net advection
!*******************************************************************************************************    
        use shr_sys_mod , only : shr_sys_flush
        use RtmVar         , only : iulog, ngeom, nlayers
        
        implicit none
        integer, intent(in) :: iunit              
        !real(r8),dimension(nlayers), intent(in) :: dv_nt    ! layer net inflow/outflow (m3/s)
        real(r8) :: num,dem,delta_z,delta_a    !
        integer :: i,j,k,mm                            ! indices
        
        TLake_r%d_z(iunit,1)=0._r8
        TLake_r%a_d(iunit,1)=TLake_r%a_di(iunit,1)
        TLake_r%v_zt(iunit,1)=TLake_r%v_zti(iunit,1)                                
        do i=1,nlayers             ! check layers for available volume to satisfy net outflow
            if (i<=TLake_r%d_ns(iunit))  then
                if(-TLake_r%dv_nt(iunit,i) > TLake_r%d_v(iunit,i))then !current layer collapses, hence remaining volume taken from next upper/lower layer
                    write(iulog,*) 'lake-r over depleting', iunit, TLake_r%d_ns(iunit), i, TLake_r%d_v(iunit,i), TLake_r%dv_nt(iunit,i) 
                    call shr_sys_abort('lake-r over depleting')
                    TLake_r%v_zn(iunit,i)=0._r8
                    if (i<TLake_r%d_ns(iunit)-1 .and. TLake_r%d_ns(iunit)>1) then
                        TLake_r%v_zn(iunit,i+1)=((TLake_r%v_zn(iunit,i+1) - (-TLake_r%dv_nt(iunit,i)-TLake_r%d_v(iunit,i))))            
                        if (TLake_r%v_zn(iunit,i+1)<0._r8) then 
                           write(iulog,*) 'lake-r negative storage', iunit, TLake_r%d_ns(iunit), i, TLake_r%v_zn(iunit,i+1) 
                           call shr_sys_abort('lake-r negative storage')
                        end if                        
                    elseif (i==TLake_r%d_ns(iunit)-1 .and. TLake_r%d_ns(iunit)>2) then    !avoid top layer collapse, hence remaining mass taken from next lower layer
                        TLake_r%v_zn(iunit,i-1)=((TLake_r%v_zn(iunit,i-1) - (-TLake_r%dv_nt(iunit,i)-TLake_r%d_v(iunit,i))))
                        if (TLake_r%v_zn(iunit,i-1)<0._r8) then 
                           write(iulog,*) 'lake-r negative storage', iunit, TLake_r%d_ns(iunit), i, TLake_r%v_zn(iunit,i-1) 
                           call shr_sys_abort('lake-r negative storage')
                        end if                        
                        mm=i-1
                        do k=mm,TLake_r%d_ns(iunit)    
                            TLake_r%v_zt(iunit,k+1)=(TLake_r%v_zt(iunit,k)+(TLake_r%v_zn(iunit,k)))
                            do j=2,ngeom+1
                                if (TLake_r%v_zt(iunit,k+1)>TLake_r%v_zti(iunit,j-1).and.TLake_r%v_zt(iunit,k+1)<=TLake_r%v_zti(iunit,j))then
                                    num = (TLake_r%v_zt(iunit,k+1)-TLake_r%v_zti(iunit,j-1))
                                    dem = (TLake_r%v_zti(iunit,j)-TLake_r%v_zti(iunit,j-1))
                                    delta_z  = (TLake_r%d_zi(iunit,j)-TLake_r%d_zi(iunit,j-1))*num/dem
                                    TLake_r%d_z(iunit,k+1) = TLake_r%d_zi(iunit,j-1) + delta_z 
                                    TLake_r%a_d(iunit,k+1) = lake_h2a(TLake_r%d_z(iunit,k+1), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
                                elseif (TLake_r%v_zt(iunit,k+1)>TLake_r%v_zti(iunit,ngeom+1))then ! inflow causes maximum storage
                                    num = (TLake_r%v_zt(iunit,k+1)-TLake_r%v_zti(iunit,ngeom))
                                    dem = (TLake_r%v_zti(iunit,ngeom+1)-TLake_r%v_zti(iunit,ngeom))
                                    delta_z  = (TLake_r%d_zi(iunit,ngeom+1)-TLake_r%d_zi(iunit,ngeom))*num/dem
                                    !delta_a  = (TLake_r%a_di(iunit,TLake_r%d_ns(iunit)+1) -TLake_r%a_di(iunit,TLake_r%d_ns(iunit)))*num/dem
                                    TLake_r%d_z(iunit,k+1) = TLake_r%d_zi(iunit,ngeom) + delta_z
                                    !TLake_r%a_d(iunit,k+1) = ((TLake_r%a_di(iunit,ngeom) + max((delta_a),0.0_r8)))
                                    TLake_r%a_d(iunit,k+1) = lake_h2a(TLake_r%d_z(iunit,k+1), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
                                end if
                            end do
                            TLake_r%d_z(iunit,k)=TLake_r%d_z(iunit,k+1)
                            TLake_r%dd_z(iunit,k)=TLake_r%d_z(iunit,k+1)-TLake_r%d_z(iunit,k)
                        end do        
                    elseif ((i==TLake_r%d_ns(iunit)-1 .or. i==TLake_r%d_ns(iunit)) .and. TLake_r%d_ns(iunit)<=2) then    ! Top layer collapses, skip outflow                            
                        if (TLake_r%v_zn(iunit,i)<0._r8) then 
                            !write(unit=1009,fmt="(i4, i4, (e14.6))") 1, i, TLake_r%v_zn(iunit,i)
                           write(iulog,*) 'lake-r negative storage', iunit, TLake_r%d_ns(iunit), i, TLake_r%v_zn(iunit,i) 
                           call shr_sys_abort('lake-r negative storage')
                        end if                        
                    end if
                    TLake_r%v_zt(iunit,i+1)=(TLake_r%v_zt(iunit,i)+(TLake_r%v_zn(iunit,i)))
                    do j=2,ngeom+1
                        if (TLake_r%v_zt(iunit,i+1)>TLake_r%v_zti(iunit,j-1).and.TLake_r%v_zt(iunit,i+1)<=TLake_r%v_zti(iunit,j))then
                            num = (TLake_r%v_zt(iunit,i+1)-TLake_r%v_zti(iunit,j-1))
                            dem = (TLake_r%v_zti(iunit,j)-TLake_r%v_zti(iunit,j-1))
                            delta_z  = (TLake_r%d_zi(iunit,j)-TLake_r%d_zi(iunit,j-1))*num/dem
                            !delta_a  = (TLake_r%a_di(iunit,j)-TLake_r%a_di(iunit,j-1))*num/dem
                            TLake_r%d_z(iunit,i+1) = TLake_r%d_zi(iunit,j-1) + delta_z 
                            !TLake_r%a_d(iunit,i+1) = ((TLake_r%a_di(iunit,j-1) + max((delta_a),0.0_r8)))
                            TLake_r%a_d(iunit,i+1) = lake_h2a(TLake_r%d_z(iunit,i+1), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
                        elseif (TLake_r%v_zt(iunit,i+1)>TLake_r%v_zti(iunit,ngeom+1))then! inflow causes maximum storage
                            num = (TLake_r%v_zt(iunit,i+1)-TLake_r%v_zti(iunit,ngeom))
                            dem = (TLake_r%v_zti(iunit,ngeom+1)-TLake_r%v_zti(iunit,ngeom))
                            delta_z  = (TLake_r%d_zi(iunit,ngeom+1)-TLake_r%d_zi(iunit,ngeom))*num/dem
                            !delta_a  = (TLake_r%a_di(iunit,TLake_r%d_ns(iunit)+1)*1._r8 -TLake_r%a_di(iunit,TLake_r%d_ns(iunit)))*num/dem
                            TLake_r%d_z(iunit,i+1) = TLake_r%d_zi(iunit,ngeom) + delta_z
                            !TLake_r%a_d(iunit,i+1) = ((TLake_r%a_di(iunit,ngeom) + max((delta_a),0.0_r8)))
                            TLake_r%a_d(iunit,i+1) = lake_h2a(TLake_r%d_z(iunit,i+1), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
                        end if
                    end do
                    TLake_r%d_z(iunit,i)  = TLake_r%d_z(iunit,i+1)
                    TLake_r%dd_z(iunit,i) = TLake_r%d_z(iunit,i+1)-TLake_r%d_z(iunit,i)
                else !enough volume, layers don't collapses
                    TLake_r%v_zt(iunit,i+1)=(TLake_r%v_zt(iunit,i) + TLake_r%v_zn(iunit,i))!TLake_r%dv_nt(iunit,i)
                    if (TLake_r%v_zt(iunit,i+1) <= 1.0_r8*TLake_r%v_zti(iunit,ngeom+1)) then ! Allowable maximum storage should include only 5% (freeboard provision)
                        do j=2,ngeom+1
                            if (TLake_r%v_zt(iunit,i+1)>TLake_r%v_zti(iunit,j-1).and.TLake_r%v_zt(iunit,i+1)<=TLake_r%v_zti(iunit,j))then
                                num = (TLake_r%v_zt(iunit,i+1)-TLake_r%v_zti(iunit,j-1))
                                dem = (TLake_r%v_zti(iunit,j)-TLake_r%v_zti(iunit,j-1))
                                delta_z  = (TLake_r%d_zi(iunit,j)-TLake_r%d_zi(iunit,j-1))*num/dem
                                !delta_a  = (TLake_r%a_di(iunit,j)-TLake_r%a_di(iunit,j-1))*num/dem
                                TLake_r%d_z(iunit,i+1) = TLake_r%d_zi(iunit,j-1) + delta_z 
                                !TLake_r%a_d(iunit,i+1) = ((TLake_r%a_di(iunit,j-1) + max((delta_a),0.0_r8)))
                                TLake_r%a_d(iunit,i+1) = lake_h2a(TLake_r%d_z(iunit,i+1), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
                            elseif (TLake_r%v_zt(iunit,i+1)>TLake_r%v_zti(iunit,ngeom+1))then ! inflow causes maximum storage
                                num = (TLake_r%v_zt(iunit,i+1)-TLake_r%v_zti(iunit,ngeom))
                                dem = (TLake_r%v_zti(iunit,ngeom+1)-TLake_r%v_zti(iunit,ngeom))
                                delta_z  = (TLake_r%d_zi(iunit,ngeom+1)-TLake_r%d_zi(iunit,ngeom))*num/dem
                                !delta_a  = (TLake_r%a_di(iunit,TLake_r%d_ns(iunit)+1) -TLake_r%a_di(iunit,TLake_r%d_ns(iunit)))*num/dem
                                TLake_r%d_z(iunit,i+1) = TLake_r%d_zi(iunit,ngeom) + delta_z
                                !TLake_r%a_d(iunit,i+1) = ((TLake_r%a_di(iunit,ngeom) + max((delta_a),0.0_r8)))
                                TLake_r%a_d(iunit,i+1) = lake_h2a(TLake_r%d_z(iunit,i+1), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
                            end if
                        end do
                    else
                        !TLake_r%d_z(iunit,i+1)  = TLake_r%d_zi(iunit,ngeom+1) + (TLake_r%v_zt(iunit,i+1) - TLake_r%v_zti(iunit,ngeom+1))/TLake_r%a_di(iunit,ngeom+1)
                        !TLake_r%a_d(iunit,i+1)  = TLake_r%a_di(iunit,ngeom+1)
                        !!TLake_r%v_zt(iunit,i+1) = 1.05_r8*TLake_r%v_zti(iunit,ngeom+1)
                        write(iulog,*) 'lake-r storage capacity exceedance not handled', iunit, TLake_r%v_zt(iunit,i+1), TLake_r%v_zti(iunit,ngeom+1) 
                        TLake_r%d_z(iunit,i+1)  = TLake_r%d_zi(iunit,ngeom+1) + (TLake_r%v_zt(iunit,i+1) - TLake_r%v_zti(iunit,ngeom+1))/TLake_r%a_di(iunit,ngeom+1)
                        TLake_r%a_d(iunit,i+1)  = TLake_r%a_di(iunit,ngeom+1)
                        call shr_sys_abort('lake-r storage capacity not handled')
                    end if
                end if 
                TLake_r%dd_z(iunit,i) = TLake_r%d_z(iunit,i+1)-TLake_r%d_z(iunit,i)
            elseif (i==TLake_r%d_ns(iunit)+1)then
                TLake_r%d_z(iunit,i)  = TLake_r%d_z(iunit,i)
                TLake_r%a_d(iunit,i)  = TLake_r%a_d(iunit,i)
                TLake_r%v_zt(iunit,i) = TLake_r%v_zt(iunit,i)
                TLake_r%dd_z(iunit,i) = 0._r8
            else
                TLake_r%d_z(iunit,i)  = 0._r8
                TLake_r%a_d(iunit,i)  = 0._r8
                TLake_r%v_zt(iunit,i) = 0._r8
                TLake_r%dd_z(iunit,i) = 0._r8
            end if
        end do
        
        !     Recalculate layer thickness and volume    
        do j = 1,nlayers 
            if (j<=TLake_r%d_ns(iunit))  then   
                TLake_r%d_v(iunit,j)     = TLake_r%v_zt(iunit,j+1) - TLake_r%v_zt(iunit,j)
                TLake_r%v_zn(iunit,j)    = TLake_r%d_v(iunit,j)
            else
                TLake_r%d_v(iunit,j)     = 0._r8
                TLake_r%v_zn(iunit,j)    = 0._r8
            end if
        end do
        
        !    Recalculare reservoir depth
        TLake_r%d_lake(iunit)=0._r8 
        do i=1,TLake_r%d_ns(iunit)
            TLake_r%d_lake(iunit) = TLake_r%d_lake(iunit)+TLake_r%dd_z(iunit,i)
        end do    
                 
    end subroutine lakegeom_update_r
 
    subroutine layer_merge_r(iunit,i)
!*******************************************************************************************************
!     Merge non-top layer if the thickness is less than threshold ddz_min 
!*******************************************************************************************************   
        use shr_sys_mod , only : shr_sys_flush
        use RtmVar         , only : iulog, ngeom, nlayers
        
        implicit none
        integer, intent(in) :: iunit,i              
        real(r8) :: ta,tb,dd_za,dd_zb,d_va,d_vb
        integer :: j,k,m,n                            ! indices                
        real(r8) :: myTINYVALUE ! 
        !real(r8) :: v_sum1, v_sum2, d_sum ! total lake volume [m^3]
        
        myTINYVALUE = 1e-6_r8
        
        if(i<TLake_r%d_ns(iunit)-1 .or. (i==1 .and. TLake_r%d_ns(iunit)==2)) then
            m=i+1
        else 
            m=i-1 !Avoid merging to top layer
        end if
        
        ta         = TLake_r%temp_lake(iunit,i)
        tb         = TLake_r%temp_lake(iunit,m)
        dd_za    = TLake_r%dd_z(iunit,i)
        dd_zb    = TLake_r%dd_z(iunit,m)
        d_va    = TLake_r%d_v(iunit,i)
        d_vb    = TLake_r%d_v(iunit,m)
        
    !    Merge layers
        if(d_va+d_vb > 0._r8)then
            TLake_r%temp_lake(iunit,i)    = (d_va*ta+d_vb*tb)/(d_va+d_vb)
        else
        end if
    
    !     Adjust new layer volume, thickness, mass and inflow/outflow 
        TLake_r%dd_z(iunit,i)    = dd_za+dd_zb
        TLake_r%d_v(iunit,i)    = d_va+d_vb
       
        if(i==1 .and. TLake_r%d_ns(iunit)==2) then
            TLake_r%temp_lake(iunit,i)    = TLake_r%temp_lake(iunit,i)
            TLake_r%d_z(iunit,i+1)        = TLake_r%d_z(iunit,i+2)
            TLake_r%v_zt(iunit,i+1)        = TLake_r%v_zt(iunit,i+2)
            TLake_r%a_d(iunit,i+1)        = TLake_r%a_d(iunit,i+2)                    
            do j=2,nlayers-1
                TLake_r%temp_lake(iunit,j)    = 0._r8
                TLake_r%dd_z(iunit,j)        = 0._r8
                TLake_r%d_v(iunit,j)        = 0._r8
                TLake_r%d_z(iunit,j+1)        = 0._r8
                TLake_r%v_zt(iunit,j+1)        = 0._r8
                TLake_r%a_d(iunit,j+1)        = 0._r8                
            end do
        else
        !    Re-number layers before collapse
            do j=i,nlayers-1
                if (j==i .and. i<(TLake_r%d_ns(iunit)-1)) then !     Identify lower (small and to be merged) and upper (larger) layer 
                    TLake_r%temp_lake(iunit,j)    = TLake_r%temp_lake(iunit,i)
                    TLake_r%dd_z(iunit,j)        = TLake_r%dd_z(iunit,i)
                    TLake_r%d_v(iunit,j)        = TLake_r%d_v(iunit,i)
                    TLake_r%d_z(iunit,j+1)        = TLake_r%d_z(iunit,i+2)
                    TLake_r%v_zt(iunit,j+1)        = TLake_r%v_zt(iunit,i+2)
                    TLake_r%a_d(iunit,j+1)        = TLake_r%a_d(iunit,i+2)        
                elseif (j<TLake_r%d_ns(iunit) .and. i<(TLake_r%d_ns(iunit)-1)) then
                    TLake_r%temp_lake(iunit,j)    = TLake_r%temp_lake(iunit,j+1)
                    TLake_r%dd_z(iunit,j)        = TLake_r%dd_z(iunit,j+1)
                    TLake_r%d_v(iunit,j)        = TLake_r%d_v(iunit,j+1)
                    TLake_r%d_z(iunit,j+1)        = TLake_r%d_z(iunit,j+2)
                    TLake_r%v_zt(iunit,j+1)        = TLake_r%v_zt(iunit,j+2)
                    TLake_r%a_d(iunit,j+1)        = TLake_r%a_d(iunit,j+2)                                    
                elseif (j==i .and. i==(TLake_r%d_ns(iunit)-1)) then
                    TLake_r%temp_lake(iunit,j-1)    = TLake_r%temp_lake(iunit,i)
                    TLake_r%dd_z(iunit,j-1)            = TLake_r%dd_z(iunit,i)
                    TLake_r%d_v(iunit,j-1)            = TLake_r%d_v(iunit,i)
                    TLake_r%d_z(iunit,j)            = TLake_r%d_z(iunit,j+1)
                    TLake_r%v_zt(iunit,j)            = TLake_r%v_zt(iunit,j+1)
                    TLake_r%a_d(iunit,j)            = TLake_r%a_d(iunit,j+1)                
                    TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit)-1)=TLake_r%temp_lake(iunit,TLake_r%d_ns(iunit))
                    TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1)    = TLake_r%dd_z(iunit,TLake_r%d_ns(iunit))
                    TLake_r%d_v(iunit,TLake_r%d_ns(iunit)-1)    = TLake_r%d_v(iunit,TLake_r%d_ns(iunit))
                    TLake_r%d_z(iunit,TLake_r%d_ns(iunit))    = TLake_r%d_z(iunit,TLake_r%d_ns(iunit)+1)
                    TLake_r%v_zt(iunit,TLake_r%d_ns(iunit))    = TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1)
                    TLake_r%a_d(iunit,TLake_r%d_ns(iunit))    = TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)                
                elseif (j>=TLake_r%d_ns(iunit)) then
                    TLake_r%temp_lake(iunit,j)    = 0._r8
                    TLake_r%dd_z(iunit,j)        = 0._r8
                    TLake_r%d_v(iunit,j)        = 0._r8
                    TLake_r%d_z(iunit,j+1)        = 0._r8
                    TLake_r%v_zt(iunit,j+1)        = 0._r8
                    TLake_r%a_d(iunit,j+1)        = 0._r8                
                end if
            end do
        end if
        TLake_r%d_ns(iunit)    = TLake_r%d_ns(iunit)-1        
        if (TLake_r%d_ns(iunit)<1) then
           write(iulog,*) 'Error in lake layer merging ! ', iunit, TLake_r%d_ns(iunit)
           call shr_sys_abort('mosart-lake: no layers in the lake')           
        end if
        if (TLake_r%d_ns(iunit)>=nlayers) then
           write(iulog,*) 'Error in lake layer splitting ! ', iunit, TLake_r%d_ns(iunit)
           call shr_sys_abort('mosart-lake: too many layers in the lake')           
        end if

    end subroutine layer_merge_r

    subroutine toplayer_merge_r(iunit)
!*******************************************************************************************************
!     Merge top layer if the thickness is less than threshold ddz_top_min 
!*******************************************************************************************************   
        use shr_sys_mod , only : shr_sys_flush
        use RtmVar         , only : iulog, ngeom, nlayers
        
        implicit none
        integer, intent(in) :: iunit              
        real(r8) :: ta,tb,dd_za,dd_zb,d_va,d_vb
        integer :: i,j,k,m,n                            ! indices                
        real(r8) :: myTINYVALUE ! 
        real(r8) :: v_sum1, v_sum2, d_sum1, d_sum2 ! total lake volume [m^3]
        
        v_sum1 = sum(TLake_r%d_v(iunit,:))
        d_sum1 = sum(TLake_r%dd_z(iunit,:))
        myTINYVALUE = 1e-6_r8
               
        i = TLake_r%d_ns(iunit)-1
        m = TLake_r%d_ns(iunit)
        ta         = TLake_r%temp_lake(iunit,i)
        tb         = TLake_r%temp_lake(iunit,m)
        dd_za    = TLake_r%dd_z(iunit,i)
        dd_zb    = TLake_r%dd_z(iunit,m)
        d_va    = TLake_r%d_v(iunit,i)
        d_vb    = TLake_r%d_v(iunit,m)
        
    !    Merge layers
        if(d_va+d_vb > 0._r8)then
            TLake_r%temp_lake(iunit,i)    = (d_va*ta+d_vb*tb)/(d_va+d_vb)
        else
            return
        end if
    
    !     Adjust new layer volume, thickness, mass and inflow/outflow 
        TLake_r%dd_z(iunit,i)    = dd_za+dd_zb
        TLake_r%d_v(iunit,i)    = d_va+d_vb
        !if (TLake_r%d_ns(iunit) > 2) then
        !    TLake_r%d_z(iunit,TLake_r%d_ns(iunit)-1)    = TLake_r%d_z(iunit,TLake_r%d_ns(iunit))
        !    TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)-1)    = TLake_r%v_zt(iunit,TLake_r%d_ns(iunit))
        !    TLake_r%a_d(iunit,TLake_r%d_ns(iunit)-1)    = TLake_r%a_d(iunit,TLake_r%d_ns(iunit))
        !end if        
        TLake_r%d_z(iunit,TLake_r%d_ns(iunit))    = TLake_r%d_z(iunit,TLake_r%d_ns(iunit)+1)
        TLake_r%v_zt(iunit,TLake_r%d_ns(iunit))    = TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1)
        TLake_r%a_d(iunit,TLake_r%d_ns(iunit))    = TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)                
        do j=TLake_r%d_ns(iunit),nlayers-1
            TLake_r%temp_lake(iunit,j)    = 0._r8
            TLake_r%dd_z(iunit,j)        = 0._r8
            TLake_r%d_v(iunit,j)        = 0._r8
            TLake_r%d_z(iunit,j+1)        = 0._r8
            TLake_r%v_zt(iunit,j+1)        = 0._r8
            TLake_r%a_d(iunit,j+1)        = 0._r8                
        end do
        
        TLake_r%d_ns(iunit)    = TLake_r%d_ns(iunit)-1

        if (TLake_r%d_ns(iunit)<1) then
           write(iulog,*) 'Error in lake top layer merging ! ', iunit, TLake_r%d_ns(iunit)
           call shr_sys_abort('mosart-lake: no layers in the lake')           
        end if
        if (TLake_r%d_ns(iunit)>=nlayers) then
           write(iulog,*) 'Error in lake layer splitting ! ', iunit, TLake_r%d_ns(iunit)
           call shr_sys_abort('mosart-lake: too many layers in the lake')           
        end if

    end subroutine toplayer_merge_r

    subroutine layer_split_r(iunit,i)
!*******************************************************************************************************
!     Split layer if the thickness is more than threshold ddz_max
!*******************************************************************************************************    
        use shr_sys_mod , only : shr_sys_flush
        use RtmVar         , only : iulog, ngeom, nlayers
        
        implicit none
        integer, intent(in) :: iunit,i             
        real(r8) :: tab,dd_zab,d_vab,delta_z    
        integer :: j,k,m,n                            ! indices
        real(r8) :: myTINYVALUE ! 
        real(r8) :: v_sum1, d_sum1, v_sum2, d_sum2 ! total lake volume [m^3]
        
        myTINYVALUE = 1e-6_r8
               
        ! Calculate layer geometric properties to be splitted
        dd_zab     = TLake_r%dd_z(iunit,i)
        d_vab    = TLake_r%d_v(iunit,i)
        tab        = TLake_r%temp_lake(iunit,i)
        !    Re-number layers before dividing layer 
        TLake_r%d_z(iunit,TLake_r%d_ns(iunit)+2)    = TLake_r%d_z(iunit,TLake_r%d_ns(iunit)+1)
        TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+2)   = TLake_r%v_zt(iunit,TLake_r%d_ns(iunit)+1)
        TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+2)    = TLake_r%a_d(iunit,TLake_r%d_ns(iunit)+1)
        m=i+1
        do j=m,TLake_r%d_ns(iunit)
            k=TLake_r%d_ns(iunit)-j+i+1
            TLake_r%temp_lake(iunit,k+1) = TLake_r%temp_lake(iunit,k)
            TLake_r%d_z(iunit,k+1)       = TLake_r%d_z(iunit,k)
            TLake_r%a_d(iunit,k+1)       = TLake_r%a_d(iunit,k)
            TLake_r%v_zt(iunit,k+1)      = TLake_r%v_zt(iunit,k)
            TLake_r%d_v(iunit,k+1)       = TLake_r%d_v(iunit,k)
            TLake_r%dd_z(iunit,k+1)      = TLake_r%dd_z(iunit,k)
        end do
    !    Divide layer in half and calculate corresponding properties
        TLake_r%dd_z(iunit,i+1)  = 0.5_r8 * dd_zab !ddz_upper
        TLake_r%dd_z(iunit,i)    = 0.5_r8 * dd_zab !dd_zab - ddz_upper
        TLake_r%d_z(iunit,i+1)   = TLake_r%d_z(iunit,i)+TLake_r%dd_z(iunit,i)
		TLake_r%a_d(iunit,i+1)   = lake_h2a(TLake_r%d_z(iunit,i+1), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
		TLake_r%v_zt(iunit,i+1)  = lake_h2v(TLake_r%d_z(iunit,i+1), TUnit_lake_r%para_a(iunit), TUnit_lake_r%para_b(iunit))
        TLake_r%d_v(iunit,i+1)   = TLake_r%v_zt(iunit,i+2)-TLake_r%v_zt(iunit,i+1)
		TLake_r%d_v(iunit,i)     = TLake_r%v_zt(iunit,i+1)-TLake_r%v_zt(iunit,i)
        TLake_r%temp_lake(iunit,i)      = tab
        TLake_r%temp_lake(iunit,i+1)    = tab
        TLake_r%d_ns(iunit)      = TLake_r%d_ns(iunit)+1    

        if (TLake_r%d_ns(iunit)>=nlayers) then
           write(iulog,*) 'Error in lake layer splitting ! ', iunit, TLake_r%d_ns(iunit)
           call shr_sys_abort('mosart-lake: too many layers in the lake')           
        end if
       
    end subroutine layer_split_r
      
  function lake_h2a(h,para_a,para_b) result(area_)
  ! calculate lake surface evaporation
      implicit none
      real(r8), intent(in) :: h      ! lake depth (m)
      real(r8), intent(in) :: para_a ! parameter a for the h~A relationship (power-law)
      real(r8), intent(in) :: para_b ! parameter b for the h~A relationship (power-law)
      real(r8)             :: area_  ! lake area corresponding to h (m2)

      area_ = para_a * h**para_b
	  
	  !area_ = area_ * 1e6 ! km2 --> m2

  end function lake_h2a

  function lake_h2v(h,para_a,para_b) result(vol_)
  ! calculate lake surface evaporation
      implicit none
      real(r8), intent(in) :: h      ! lake depth (s)
      real(r8), intent(in) :: para_a ! parameter a for the h~A relationship (power-law)
      real(r8), intent(in) :: para_b ! parameter b for the h~A relationship (power-law)
      real(r8)             :: vol_   ! lake volume

      real(r8) :: para_c ! parameter a for the h~V relationship
      real(r8) :: para_d ! parameter a for the h~V relationship
      
      para_c = para_a / (1._r8 + para_b)
      para_d = 1._r8 + para_b
      vol_ = para_c * h**para_d
	  
	  !vol_ = vol_ * 1e9 ! km3-m3

  end function lake_h2v
  
end MODULE MOSART_lake_geometry_mod