
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin

 oceanmixed = "oceanmixed_ice.nc"

 f1 = addfile(oceanmixed,"r")
 f3 = addfile("/fis/cgd/cseg/csm/inputdata/ocn/docn7/domain.ocn.gx1v3.051111.nc","r")

 time = (/14., 46., 74., 105., 135., 166., 196., 227., 258., 288., 319., 349./)
 maskr = f1->REGION_MASK
 delete(maskr@coordinates)
 area = f3->area
 delete(area@coordinates)
 tarea = f1->TAREA
 tlon = f1->TLONG
 tlat = f1->TLAT
 xc = flt2dble(tlon)
 yc = flt2dble(tlat)
 dims = dimsizes(xc)
 nj = dims(0)
 ni = dims(1)
 ntime = 12

; Use the annual mean mixed layer depth
 hbltin = f1->HBLT
 Tin = f1->TEMP(:,0,:,:)
 Sin = f1->SALT(:,0,:,:)
 Uin = f2->avXc2i_i_So_u
 Vin = f2->avXc2i_i_So_v
 dhdxin = f2->avXc2i_i_So_dhdx
 dhdyin = f2->avXc2i_i_So_dhdy

; Need to weight the monthly means

 daysinmo = (/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)
 xnp = daysinmo
 xnm = daysinmo
 xnm(1:11) = daysinmo(1:11)+daysinmo(0:10)
 xnm(0) = daysinmo(0)+daysinmo(11)
 xnp(0:10) = daysinmo(0:10)+daysinmo(1:11)
 xnp(11) = daysinmo(11)+daysinmo(0)
 aa = 2.*daysinmo / xnm
 cc = 2.*daysinmo / xnp
 a = aa / 8.
 c = cc / 8.
 b = 1. - a - c
 M = (/(/b(0),c(0),0,0,0,0,0,0,0,0,0,a(0)/), \
       (/a(1),b(1),c(1),0,0,0,0,0,0,0,0,0/), \
       (/0,a(2),b(2),c(2),0,0,0,0,0,0,0,0/), \
       (/0,0,a(3),b(3),c(3),0,0,0,0,0,0,0/), \
       (/0,0,0,a(4),b(4),c(4),0,0,0,0,0,0/), \
       (/0,0,0,0,a(5),b(5),c(5),0,0,0,0,0/), \
       (/0,0,0,0,0,a(6),b(6),c(6),0,0,0,0/), \
       (/0,0,0,0,0,0,a(7),b(7),c(7),0,0,0/), \
       (/0,0,0,0,0,0,0,a(8),b(8),c(8),0,0/), \
       (/0,0,0,0,0,0,0,0,a(9),b(9),c(9),0/), \
       (/0,0,0,0,0,0,0,0,0,a(10),b(10),c(10)/), \
       (/c(11),0,0,0,0,0,0,0,0,0,a(11),b(11)/)/)
 invM = inverse_matrix(M)

 shf = f1->SHF
 qflux = f1->QFLUX
 melth_f = f1->MELTH_F
 resid_t = f1->RESID_T
;tfw_t = f1->TWF_T

 rcp_sw = 1026.*3996.
 surf = shf+qflux
 T1 = Tin
 T1(0:10,:,:) = Tin(1:11,:,:)
 T1(11,:,:) = Tin(0,:,:)
 T2 = Tin
 T2(0,:,:) = Tin(11,:,:)
 T2(1:11,:,:) = Tin(0:10,:,:)
 dT = T1-T2
 release = rcp_sw*dT*hblttmp / (86400.*365./6.)

 ocnheat = surf-release
 maskt = new((/nj,ni/),typeof(ocnheat))
 maskt = 1
 maskt = mask(maskt,ismissing(ocnheat(0,:,:)),False)
 err = new(12,typeof(ocnheat))
 do n=0,ntime-1
    tmp = ndtooned(ocnheat(n,:,:))
    tmp(ind(ismissing(tmp))) = 0.
    err(n) = tmp # ndtooned(tarea) / sum(tarea * maskt)
 end do
 print(err)
 glob = avg(err)
 print(glob)
 ocnheat = ocnheat - glob

 T = new(dimsizes(Tin),typeof(Tin))
 S = new(dimsizes(Sin),typeof(Sin))
 U = new(dimsizes(Uin),typeof(Uin))
 V = new(dimsizes(Vin),typeof(Vin))
 dhdx = new(dimsizes(dhdxin),typeof(dhdxin))
 dhdy = new(dimsizes(dhdyin),typeof(dhdyin))
 hblt = new(dimsizes(hbltin),typeof(hbltin))
 qdp = new(dimsizes(shf),typeof(shf))

 T = 0.
 S = 0.
 U = 0.
 V = 0.
 dhdx = 0.
 dhdy = 0.
 hblt = 0.
 qdp = 0.

 do j=0,ntime-1
 do i=0,ntime-1
    T(j,:,:) = T(j,:,:) + invM(j,i)*Tin(i,:,:)
    S(j,:,:) = S(j,:,:) + invM(j,i)*Sin(i,:,:)
    U(j,:,:) = U(j,:,:) + invM(j,i)*Uin(i,:,:)
    V(j,:,:) = V(j,:,:) + invM(j,i)*Vin(i,:,:)
    dhdx(j,:,:) = dhdx(j,:,:) + invM(j,i)*dhdxin(i,:,:)
    dhdy(j,:,:) = dhdy(j,:,:) + invM(j,i)*dhdyin(i,:,:)
    hblt(j,:,:) = hblt(j,:,:) + invM(j,i)*hblttmp(i,:,:)
    qdp(j,:,:) = qdp(j,:,:) + invM(j,i)*ocnheat(i,:,:)
 end do
 end do

 time@units = "days since 0001-01-01 00:00:00"
 time@long_name = "observation time"
 time@calendar = "noleap"
 
 area@units = "area"
 area@long_name = "area of grid cell in radians squared"
 
 maskr@long_name = "domain maskr"
 maskr@units = "unitless"

 xc@long_name = "longitude of grid cell center"
 xc@units = "degrees east"
 
 yc@long_name = "latitude of grid cell center"
 yc@units = "degrees north"
 
 S@long_name = "salinity"
 S@units = "ppt"
 
 T@long_name = "temperature"
 T@units = "degC"
 
 U@long_name = "u ocean current"
 U@units = "m/s"
 
 V@long_name = "v ocean current"
 V@units = "m/s"
 
 dhdx@long_name = "ocean surface slope: zonal"
 dhdx@units = "m/m"
 
 dhdy@long_name = "ocean surface slope: meridional"
 dhdy@units = "m/m"
 
 hblt@long_name = "boundary layer depth"
 hblt@units = "m"
 
 qdp@long_name = "ocean heat flux convergence"
 qdp@units = "W/m^2"
 
 fout = addfile("oceanmixed_ice.nc","c")
 setfileoption(fout,"DefineMode",True)
 
 fileAtt = True
 fileAtt@title = "Monthly averaged ocean forcing from POP output"
 fileAtt@conventions = "CCSM2.0 data model domain description"
 fileAtt@source = "pop_frc.ncl"
 fileAtt@description = "Input data for CSIM4 mixed layer model from b30.004"
 fileAtt@note1 = "fields computed from 20-yr monthly means from pop"
 fileAtt@note2 = "all fields interpolated to T-grid"
 fileAtt@note3 = "qdp is computed from depth summed ocean column"
 fileAtt@author = "D. Bailey"
 fileAtt@calendar = "standard"
 fileAtt@comment = "This data is on the displaced pole grid gx1v3"
 fileAtt@creation_date = systemfunc("date")
 fileattdef(fout,fileAtt)
 
 dimNames = (/"time","nj","ni"/)
 dimSizes = (/ntime,nj,ni/)
 dimUnlim = (/False,False,False/)
 filedimdef(fout,dimNames,dimSizes,dimUnlim)
 
 filevardef(fout,"area",typeof(area),(/"nj","ni"/))
 filevarattdef(fout,"area",area)
 filevardef(fout,"mask",typeof(maskr),(/"nj","ni"/))
 filevarattdef(fout,"mask",maskr)
 filevardef(fout,"xc",typeof(xc),(/"nj","ni"/))
 filevarattdef(fout,"xc",xc)
 filevardef(fout,"yc",typeof(yc),(/"nj","ni"/))
 filevarattdef(fout,"yc",yc)

 filevardef(fout,"time",typeof(time),"time")
 filevarattdef(fout,"time",time)

 filevardef(fout,"S",typeof(S),dimNames)
 filevarattdef(fout,"S",S)
 filevardef(fout,"T",typeof(T),dimNames)
 filevarattdef(fout,"T",T)
 filevardef(fout,"U",typeof(U),dimNames)
 filevarattdef(fout,"U",U)
 filevardef(fout,"V",typeof(V),dimNames)
 filevarattdef(fout,"V",V)
 filevardef(fout,"dhdx",typeof(dhdx),dimNames)
 filevarattdef(fout,"dhdx",dhdx)
 filevardef(fout,"dhdy",typeof(dhdy),dimNames)
 filevarattdef(fout,"dhdy",dhdy)
 filevardef(fout,"hblt",typeof(hblt),dimNames)
 filevarattdef(fout,"hblt",hblt)
 filevardef(fout,"qdp",typeof(qdp),dimNames)
 filevarattdef(fout,"qdp",qdp)

 fout->area = (/area/)
 fout->mask = (/maskr/)
 fout->xc = (/xc/)
 fout->yc = (/yc/)
 fout->time = (/time/)
 fout->S = (/S/)
 fout->T = (/T/)
 fout->U = (/U/)
 fout->V = (/V/)
 fout->dhdx = (/dhdx/)
 fout->dhdy = (/dhdy/)
 fout->hblt = (/hblt/)
 fout->qdp = (/qdp/)
 
end
