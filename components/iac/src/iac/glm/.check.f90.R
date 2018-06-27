program test

integer,parameter                :: nx=360,ny=180,nh=6
real,dimension(nx,ny)            :: c,p,s,v,w,i,summer
character(len=90)                :: f
character(len=4)                 :: fyear


fyear='1730'

summer=0.0

f='gcrop.' // fyear // '.txt'
open(14,file=f,status='old')
do ih=1,nh
  read(14,*)
end do
do iy=1,ny
  read(14,*) c(:,iy)
end do
close(14)


f='gpast.' // fyear // '.txt'
open(14,file=f,status='old')
do ih=1,nh
  read(14,*)
end do
do iy=1,ny
  read(14,*) p(:,iy)
end do
close(14)

f='gothr.' // fyear // '.txt'
open(14,file=f,status='old')
do ih=1,nh
  read(14,*)
end do
do iy=1,ny
  read(14,*) v(:,iy)
end do
close(14)

f='gsecd.' // fyear // '.txt'
open(14,file=f,status='old')
do ih=1,nh
  read(14,*)
end do
do iy=1,ny
  read(14,*) s(:,iy)
end do
close(14)

f='/raid1/gh/matt/glm/processed/land_a/1deg_grids/gwatr.1700.txt'
open(14,file=f,status='old')
do ih=1,nh
  read(14,*)
end do
do iy=1,ny
  read(14,*) w(:,iy)
end do
close(14)

f='/raid1/gh/matt/glm/processed/land_a/1deg_grids/gicew.1700.txt'
open(14,file=f,status='old')
do ih=1,nh
  read(14,*)
end do
do iy=1,ny
  read(14,*) i(:,iy)
end do
close(14)


summer=c+p+v+s+i+w

do ix=1,nx
do iy=1,ny
  if(summer(ix,iy).gt.1.01.or.summer(ix,iy).lt..99) then
     print *, summer(ix,iy), c(ix,iy), p(ix,iy), s(ix,iy), v(ix,iy), w(ix,iy), i(ix,iy)
  end if
end do
end do
stop
end
