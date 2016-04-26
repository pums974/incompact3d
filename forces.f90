module forces_m
  use derivevitesse_m
implicit none
contains
!***********************************************************************
!
subroutine aerof(ux,uy,uz,ppm,epsi,gx,gy,gz,&
                 tx,ty,tz,fx,fy,fz,di1,di2,&
                 uxm1,uym1,uzm1,uxm2,uym2,uzm2)
!
!***********************************************************************
USE param
USE IBM
USE aeroforce
USE variables
      
implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,tx,ty,tz,fx,fy,fz,di1,di2,epsi
real(8),dimension(nx,ny,nz) :: uxm1,uym1,uzm1,uxm2,uym2,uzm2
real(8),dimension(nxm,nym,nzm) :: ppm
!
if (itime==(idebut)) then
   uxm2=ux
   uym2=uy
   uzm2=uz
   return
endif
if (itime==(idebut+1)) then
   uxm1=ux
   uym1=uy
   uzm1=uz
   return
endif
 
if (nz.gt.1) then
   call rmomentum(ux,uy,uz,epsi,uxm1,uym1,uzm1,uxm2,uym2,uzm2)
   call sforce(ux,uy,uz,ppm,gx,gy,gz,tx,ty,tz,fx,fy,fz,di1,di2)
   xforce = 2.*(f2x-xmom)
   yforce = 2.*(f2y-ymom)
   zforce = 2.*(f2z-zmom)
   write(*,101) (itime-1)*dt,xforce,yforce,zforce
   write(82,100) (itime-1)*dt,xforce,yforce,zforce
else
   call rmomentum2d(ux,uy,uz,epsi,uxm1,uym1,uzm1,uxm2,uym2,uzm2)
   call sforce2d(ux,uy,uz,ppm,gx,gy,gz,tx,ty,tz,fx,fy,fz,di1,di2)
   print *,'F2x,xmom=',f2x,xmom
   print *,'F2y,ymom=',f2y,ymom
   xforce = 2.*(f2x-xmom)
   yforce = 2.*(f2y-ymom)
   write(*,102) (itime-1)*dt,xforce,yforce
   write(82,100) (itime-1)*dt,xforce,yforce
endif

uxm2=uxm1
uym2=uym1
uzm2=uzm1
uxm1=ux
uym1=uy
uzm1=uz

  100 format(4F16.8)
  101 format('Time, Drag, Lift, Yaw = ',4F16.8)
  102 format('Time, Drag, Lift = ',3F16.8)
return
end subroutine aerof
!
!***********************************************************************
!
subroutine rmomentum(ux,uy,uz,epsi,uxm1,uym1,uzm1,uxm2,uym2,uzm2)
!
!***********************************************************************
USE param
USE IBM
USE aeroforce
USE variables
!
implicit none
!
real(8),dimension(nx,ny,nz) :: ux,uy,uz,epsi
real(8),dimension(nx,ny,nz) :: uxm1,uym1,uzm1,uxm2,uym2,uzm2
real(8) :: dudt1,dudt2,dudt3,dudt4,dudt5,dudt6,dudt7,dudt8
real(8) :: uxm,uym,uzm
real(8) :: sumx,sumy,sumz,duxdt,duydt,duzdt
real(8) :: fxa,fya,fza,fxb,fyb,fzb,fxc,fyc,fzc,fxd,fyd,fzd,fxe,fye,fze,fxf,fyf,fzf
integer :: i,j,k

sumx=0.0; sumy=0.0; sumz=0.0
fxa=0.0; fya=0.0; fza=0.0
fxb=0.0; fyb=0.0; fzb=0.0
fxc=0.0; fyc=0.0; fzc=0.0
fxd=0.0; fyd=0.0; fzd=0.0
fxe=0.0; fye=0.0; fze=0.0
fxf=0.0; fyf=0.0; fzf=0.0

do k=kld,kud-1
do j=jld,jud-1
do i=ild,iud-1
      dudt1 = 1.5*ux(i  ,j  ,k  )-2.0*uxm1(i  ,j  ,k  )+0.5*uxm2(i  ,j  ,k  )
      dudt2 = 1.5*ux(i+1,j  ,k  )-2.0*uxm1(i+1,j  ,k  )+0.5*uxm2(i+1,j  ,k  )
      dudt3 = 1.5*ux(i+1,j+1,k  )-2.0*uxm1(i+1,j+1,k  )+0.5*uxm2(i+1,j+1,k  )
      dudt4 = 1.5*ux(i  ,j+1,k  )-2.0*uxm1(i  ,j+1,k  )+0.5*uxm2(i  ,j+1,k  )
      dudt5 = 1.5*ux(i  ,j  ,k+1)-2.0*uxm1(i  ,j  ,k+1)+0.5*uxm2(i  ,j  ,k+1)
      dudt6 = 1.5*ux(i+1,j  ,k+1)-2.0*uxm1(i+1,j  ,k+1)+0.5*uxm2(i+1,j  ,k+1)
      dudt7 = 1.5*ux(i+1,j+1,k+1)-2.0*uxm1(i+1,j+1,k+1)+0.5*uxm2(i+1,j+1,k+1)
      dudt8 = 1.5*ux(i  ,j+1,k+1)-2.0*uxm1(i  ,j+1,k+1)+0.5*uxm2(i  ,j+1,k+1)
      duxdt = (dudt1+dudt2+dudt3+dudt4+dudt5+dudt6+dudt7+dudt8)/8.0/dt
!
      dudt1 = 1.5*uy(i  ,j  ,k  )-2.0*uym1(i  ,j  ,k  )+0.5*uym2(i  ,j  ,k  )
      dudt2 = 1.5*uy(i+1,j  ,k  )-2.0*uym1(i+1,j  ,k  )+0.5*uym2(i+1,j  ,k  )
      dudt3 = 1.5*uy(i+1,j+1,k  )-2.0*uym1(i+1,j+1,k  )+0.5*uym2(i+1,j+1,k  )
      dudt4 = 1.5*uy(i  ,j+1,k  )-2.0*uym1(i  ,j+1,k  )+0.5*uym2(i  ,j+1,k  )
      dudt5 = 1.5*uy(i  ,j  ,k+1)-2.0*uym1(i  ,j  ,k+1)+0.5*uym2(i  ,j  ,k+1)
      dudt6 = 1.5*uy(i+1,j  ,k+1)-2.0*uym1(i+1,j  ,k+1)+0.5*uym2(i+1,j  ,k+1)
      dudt7 = 1.5*uy(i+1,j+1,k+1)-2.0*uym1(i+1,j+1,k+1)+0.5*uym2(i+1,j+1,k+1)
      dudt8 = 1.5*uy(i  ,j+1,k+1)-2.0*uym1(i  ,j+1,k+1)+0.5*uym2(i  ,j+1,k+1)
      duydt = (dudt1+dudt2+dudt3+dudt4+dudt5+dudt6+dudt7+dudt8)/8.0/dt
!
      dudt1 = 1.5*uz(i  ,j  ,k  )-2.0*uzm1(i  ,j  ,k  )+0.5*uzm2(i  ,j  ,k  )
      dudt2 = 1.5*uz(i+1,j  ,k  )-2.0*uzm1(i+1,j  ,k  )+0.5*uzm2(i+1,j  ,k  )
      dudt3 = 1.5*uz(i+1,j+1,k  )-2.0*uzm1(i+1,j+1,k  )+0.5*uzm2(i+1,j+1,k  )
      dudt4 = 1.5*uz(i  ,j+1,k  )-2.0*uzm1(i  ,j+1,k  )+0.5*uzm2(i  ,j+1,k  )
      dudt5 = 1.5*uz(i  ,j  ,k+1)-2.0*uzm1(i  ,j  ,k+1)+0.5*uzm2(i  ,j  ,k+1)
      dudt6 = 1.5*uz(i+1,j  ,k+1)-2.0*uzm1(i+1,j  ,k+1)+0.5*uzm2(i+1,j  ,k+1)
      dudt7 = 1.5*uz(i+1,j+1,k+1)-2.0*uzm1(i+1,j+1,k+1)+0.5*uzm2(i+1,j+1,k+1)
      dudt8 = 1.5*uz(i  ,j+1,k+1)-2.0*uzm1(i  ,j+1,k+1)+0.5*uzm2(i  ,j+1,k+1)
      duzdt = (dudt1+dudt2+dudt3+dudt4+dudt5+dudt6+dudt7+dudt8)/8.0/dt
!
      sumx=sumx+duxdt*dx*(yp(j+1)-yp(j))*dz*(1.-epsi(i,j,k))
      sumy=sumy+duydt*dx*(yp(j+1)-yp(j))*dz*(1.-epsi(i,j,k))
      sumz=sumz+duzdt*dx*(yp(j+1)-yp(j))*dz*(1.-epsi(i,j,k))
enddo
enddo
enddo

do k=kld,kud-1
do j=jld,jud-1
   uxm = (ux(ild,j,k)+ux(ild,j+1,k)+ux(ild,j,k+1)+ux(ild,j+1,k+1))/4.0
   uym = (uy(ild,j,k)+uy(ild,j+1,k)+uy(ild,j,k+1)+uy(ild,j+1,k+1))/4.0
   uzm = (uz(ild,j,k)+uz(ild,j+1,k)+uz(ild,j,k+1)+uz(ild,j+1,k+1))/4.0
   fxa= fxa-uxm*uxm*(yp(j+1)-yp(j))*dz
   fya= fya-uxm*uym*(yp(j+1)-yp(j))*dz
   fza= fza-uxm*uzm*(yp(j+1)-yp(j))*dz
   uxm = (ux(iud,j,k)+ux(iud,j+1,k)+ux(iud,j,k+1)+ux(iud,j+1,k+1))/4.0
   uym = (uy(iud,j,k)+uy(iud,j+1,k)+uy(iud,j,k+1)+uy(iud,j+1,k+1))/4.0
   uzm = (uz(iud,j,k)+uz(iud,j+1,k)+uz(iud,j,k+1)+uz(iud,j+1,k+1))/4.0
   fxb= fxb+uxm*uxm*(yp(j+1)-yp(j))*dz
   fyb= fyb+uxm*uym*(yp(j+1)-yp(j))*dz
   fzb= fzb+uxm*uzm*(yp(j+1)-yp(j))*dz
enddo
enddo
do k=kld,kud-1
do i=ild,iud-1
   uxm = (ux(i,jld,k)+ux(i+1,jld,k)+ux(i,jld,k+1)+ux(i+1,jld,k+1))/4.0
   uym = (uy(i,jld,k)+uy(i+1,jld,k)+uy(i,jld,k+1)+uy(i+1,jld,k+1))/4.0
   uzm = (uz(i,jld,k)+uz(i+1,jld,k)+uz(i,jld,k+1)+uz(i+1,jld,k+1))/4.0
   fxc= fxc-uym*uxm*dx*dz
   fyc= fyc-uym*uym*dx*dz
   fzc= fzc-uym*uzm*dx*dz
   uxm = (ux(i,jud,k)+ux(i+1,jud,k)+ux(i,jud,k+1)+ux(i+1,jud,k+1))/4.0
   uym = (uy(i,jud,k)+uy(i+1,jud,k)+uy(i,jud,k+1)+uy(i+1,jud,k+1))/4.0
   uzm = (uz(i,jud,k)+uz(i+1,jud,k)+uz(i,jud,k+1)+uz(i+1,jud,k+1))/4.0
   fxd= fxd+uym*uxm*dx*dz
   fyd= fyd+uym*uym*dx*dz
   fzd= fzd+uym*uzm*dx*dz
enddo
enddo
do j=jld,jud-1
do i=ild,iud-1
   uxm = (ux(i,j,kld)+ux(i+1,j,kld)+ux(i,j+1,kld)+ux(i+1,j+1,kld))/4.0
   uym = (uy(i,j,kld)+uy(i+1,j,kld)+uy(i,j+1,kld)+uy(i+1,j+1,kld))/4.0
   uzm = (uz(i,j,kld)+uz(i+1,j,kld)+uz(i,j+1,kld)+uz(i+1,j+1,kld))/4.0
   fxe= fxe-uzm*uxm*dx*(yp(j+1)-yp(j))
   fye= fye-uzm*uym*dx*(yp(j+1)-yp(j))
   fze= fze-uzm*uzm*dx*(yp(j+1)-yp(j))
   uxm = (ux(i,j,kud)+ux(i+1,j,kud)+ux(i,j+1,kud)+ux(i+1,j+1,kud))/4.0
   uym = (uy(i,j,kud)+uy(i+1,j,kud)+uy(i,j+1,kud)+uy(i+1,j+1,kud))/4.0
   uzm = (uz(i,j,kud)+uz(i+1,j,kud)+uz(i,j+1,kud)+uz(i+1,j+1,kud))/4.0
   fxf= fxf+uzm*uxm*dx*(yp(j+1)-yp(j))
   fyf= fyf+uzm*uym*dx*(yp(j+1)-yp(j))
   fzf= fzf+uzm*uzm*dx*(yp(j+1)-yp(j))
enddo
enddo

xmom=sumx+fxa+fxb+fxc+fxd+fxe+fxf
ymom=sumy+fya+fyb+fyc+fyd+fye+fyf
zmom=sumz+fza+fzb+fzc+fzd+fze+fzf

return
end subroutine rmomentum
!
!***********************************************************************
!
subroutine sforce(ux,uy,uz,ppm,gx,gy,gz,tx,ty,tz,fx,fy,fz,di1,di2)
!
!***********************************************************************
!
USE param
USE IBM
USE aeroforce
USE variables
!
implicit none
!
real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,tx,ty,tz,fx,fy,fz,di1,di2
real(8),dimension(nxm,nym,nzm) :: ppm
real(8) :: foxa,foxb,foxc,foxd,foxe,foxf,foxpr
real(8) :: foya,foyb,foyc,foyd,foye,foyf,foypr
real(8) :: foza,fozb,fozc,fozd,foze,fozf,fozpr
real(8) :: dudxm,dudym,dudzm
real(8) :: dvdxm,dvdym,dvdzm
real(8) :: dwdxm,dwdym,dwdzm
real(8) :: pud,pld
integer :: i,j,k
! 
print *,'dx,dz=',dx,dz
foxa=0.0; foya=0.0; foza=0.0
foxb=0.0; foyb=0.0; fozb=0.0
foxc=0.0; foyc=0.0; fozc=0.0
foxd=0.0; foyd=0.0; fozd=0.0
foxe=0.0; foye=0.0; foze=0.0
foxf=0.0; foyf=0.0; fozf=0.0
foxpr=0.0; foypr=0.0; fozpr=0.0
!
call derx (tx,ux,di1    ,sx,ffxp,fsxp,fwxp    ,nx,ny,nz,1)
call derx (gx,uy,di1    ,sx,ffxp,fsxp,fwxp    ,nx,ny,nz,1)
call derx (fx,uz,di1    ,sx,ffxp,fsxp,fwxp    ,nx,ny,nz,1)
call dery (ty,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
call dery (gy,uy,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
call dery (fy,uz,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
call derz (tz,ux,di1    ,sz,ffzp,fszp,fwzp    ,nx,ny,nz,1)
call derz (gz,uy,di1    ,sz,ffzp,fszp,fwzp    ,nx,ny,nz,1)
call derz (fz,uz,di1    ,sz,ffzp,fszp,fwzp    ,nx,ny,nz,1)
!
do k=kld,kud-1
do j=jld,jud-1
   dudxm = (tx(ild,j,k)+tx(ild,j+1,k)+tx(ild,j,k+1)+tx(ild,j+1,k+1))/4.0
   dudym = (ty(ild,j,k)+ty(ild,j+1,k)+ty(ild,j,k+1)+ty(ild,j+1,k+1))/4.0
   dvdxm = (gx(ild,j,k)+gx(ild,j+1,k)+gx(ild,j,k+1)+gx(ild,j+1,k+1))/4.0
   dudzm = (tz(ild,j,k)+tz(ild,j+1,k)+tz(ild,j,k+1)+tz(ild,j+1,k+1))/4.0
   dwdxm = (fx(ild,j,k)+fx(ild,j+1,k)+fx(ild,j,k+1)+fx(ild,j+1,k+1))/4.0
   foxa = foxa-xnu*(dudxm+dudxm)*(yp(j+1)-yp(j))*dz
   foya = foya-xnu*(dvdxm+dudym)*(yp(j+1)-yp(j))*dz
   foza = foza-xnu*(dwdxm+dudzm)*(yp(j+1)-yp(j))*dz
   dudxm = (tx(iud,j,k)+tx(iud,j+1,k)+tx(iud,j,k+1)+tx(iud,j+1,k+1))/4.0
   dudym = (ty(iud,j,k)+ty(iud,j+1,k)+ty(iud,j,k+1)+ty(iud,j+1,k+1))/4.0
   dvdxm = (gx(iud,j,k)+gx(iud,j+1,k)+gx(iud,j,k+1)+gx(iud,j+1,k+1))/4.0
   dudzm = (tz(iud,j,k)+tz(iud,j+1,k)+tz(iud,j,k+1)+tz(iud,j+1,k+1))/4.0
   dwdxm = (fx(iud,j,k)+fx(iud,j+1,k)+fx(iud,j,k+1)+fx(iud,j+1,k+1))/4.0
   foxc = foxc+xnu*(dudxm+dudxm)*(yp(j+1)-yp(j))*dz
   foyc = foyc+xnu*(dvdxm+dudym)*(yp(j+1)-yp(j))*dz
   fozc = fozc+xnu*(dwdxm+dudzm)*(yp(j+1)-yp(j))*dz
   pld   = ppm(ild,j,k)/gdt(3)
   pud   = ppm(iud,j,k)/gdt(3)
   foxpr = foxpr+(pld-pud)*(yp(j+1)-yp(j))*dz
enddo
enddo
!
do k=kld,kud-1
do i=ild,iud-1
   dudym = (ty(i,jud,k)+ty(i+1,jud,k)+ty(i,jud,k+1)+ty(i+1,jud,k+1))/4.0
   dvdxm = (gx(i,jud,k)+gx(i+1,jud,k)+gx(i,jud,k+1)+gx(i+1,jud,k+1))/4.0
   dvdym = (gy(i,jud,k)+gy(i+1,jud,k)+gy(i,jud,k+1)+gy(i+1,jud,k+1))/4.0
   dwdym = (fy(i,jud,k)+fy(i+1,jud,k)+fy(i,jud,k+1)+fy(i+1,jud,k+1))/4.0
   dvdzm = (gz(i,jud,k)+gz(i+1,jud,k)+gz(i,jud,k+1)+gz(i+1,jud,k+1))/4.0
   foxb = foxb+xnu*(dudym+dvdxm)*dx*dz
   foyb = foyb+xnu*(dvdym+dvdym)*dx*dz
   fozb = fozb+xnu*(dwdym+dvdzm)*dx*dz
   dudym = (ty(i,jld,k)+ty(i+1,jld,k)+ty(i,jld,k+1)+ty(i+1,jld,k+1))/4.0
   dvdxm = (gx(i,jld,k)+gx(i+1,jld,k)+gx(i,jld,k+1)+gx(i+1,jld,k+1))/4.0
   dvdym = (gy(i,jld,k)+gy(i+1,jld,k)+gy(i,jld,k+1)+gy(i+1,jld,k+1))/4.0
   dwdym = (fy(i,jld,k)+fy(i+1,jld,k)+fy(i,jld,k+1)+fy(i+1,jld,k+1))/4.0
   dvdzm = (gz(i,jld,k)+gz(i+1,jld,k)+gz(i,jld,k+1)+gz(i+1,jld,k+1))/4.0
   foxd = foxd-xnu*(dudym+dvdxm)*dx*dz
   foyd = foyd-xnu*(dvdym+dvdym)*dx*dz
   fozd = fozd-xnu*(dwdym+dvdzm)*dx*dz
   pld  = ppm(i,jld,k)/gdt(3)
   pud  = ppm(i,jud,k)/gdt(3)
   foypr = foypr+(pld-pud)*dx*dz
enddo
enddo
!
do j=jld,jud-1
do i=ild,iud-1
   dudzm = (tz(i,j,kud)+tz(i+1,j,kud)+tz(i,j+1,kud)+tz(i+1,j+1,kud))/4.0
   dvdzm = (gz(i,j,kud)+gz(i+1,j,kud)+gz(i,j+1,kud)+gz(i+1,j+1,kud))/4.0
   dwdxm = (fx(i,j,kud)+fx(i+1,j,kud)+fx(i,j+1,kud)+fx(i+1,j+1,kud))/4.0
   dwdym = (fy(i,j,kud)+fy(i+1,j,kud)+fy(i,j+1,kud)+fy(i+1,j+1,kud))/4.0
   dwdzm = (fz(i,j,kud)+fz(i+1,j,kud)+fz(i,j+1,kud)+fz(i+1,j+1,kud))/4.0
   foxe = foxe+xnu*(dudzm+dwdxm)*dx*(yp(j+1)-yp(j))
   foye = foye+xnu*(dvdzm+dwdym)*dx*(yp(j+1)-yp(j))
   foze = foze+xnu*(dwdzm+dwdzm)*dx*(yp(j+1)-yp(j))
   dudzm = (tz(i,j,kld)+tz(i+1,j,kld)+tz(i,j+1,kld)+tz(i+1,j+1,kld))/4.0
   dvdzm = (gz(i,j,kld)+gz(i+1,j,kld)+gz(i,j+1,kld)+gz(i+1,j+1,kld))/4.0
   dwdxm = (fx(i,j,kld)+fx(i+1,j,kld)+fx(i,j+1,kld)+fx(i+1,j+1,kld))/4.0
   dwdym = (fy(i,j,kld)+fy(i+1,j,kld)+fy(i,j+1,kld)+fy(i+1,j+1,kld))/4.0
   dwdzm = (fz(i,j,kld)+fz(i+1,j,kld)+fz(i,j+1,kld)+fz(i+1,j+1,kld))/4.0
   foxf = foxf-xnu*(dudzm+dwdxm)*dx*(yp(j+1)-yp(j))
   foyf = foyf-xnu*(dvdzm+dwdym)*dx*(yp(j+1)-yp(j))
   fozf = fozf-xnu*(dwdzm+dwdzm)*dx*(yp(j+1)-yp(j))
   pld   = ppm(i,j,kld)/gdt(3)
   pud   = ppm(i,j,kud)/gdt(3)
   fozpr = fozpr+(pld-pud)*dx*(yp(j+1)-yp(j))
enddo
enddo
!
f2x=foxa+foxb+foxc+foxd+foxe+foxf+foxpr
f2y=foya+foyb+foyc+foyd+foye+foyf+foypr
f2z=foza+fozb+fozc+fozd+foze+fozf+fozpr
!
return
end subroutine sforce
!
!***********************************************************************
!
subroutine rmomentum2d(ux,uy,uz,epsi,uxm1,uym1,uzm1,uxm2,uym2,uzm2)
!
!***********************************************************************
USE param
USE IBM
USE aeroforce
USE variables
!
implicit none
!
real(8),dimension(nx,ny,nz) :: ux,uy,uz,epsi
real(8),dimension(nx,ny,nz) :: uxm1,uym1,uzm1,uxm2,uym2,uzm2
real(8) :: dudt1,dudt2,dudt3,dudt4
real(8) :: uxm,uym
real(8) :: sumx,sumy,duxdt,duydt
real(8) :: fxa,fya,fxb,fyb,fxc,fyc,fxd,fyd
integer :: i,j,k

sumx=0.0
sumy=0.0
fxa=0.0
fya=0.0
fxc=0.0
fyc=0.0
fxb=0.0
fyb=0.0
fxd=0.0
fyd=0.0

do j=jld,jud-1
do i=ild,iud-1
      dudt1 = 1.5*ux(i  ,j  ,1)-2.0*uxm1(i  ,j  ,1)+0.5*uxm2(i  ,j  ,1)
      dudt2 = 1.5*ux(i+1,j  ,1)-2.0*uxm1(i+1,j  ,1)+0.5*uxm2(i+1,j  ,1)
      dudt3 = 1.5*ux(i+1,j+1,1)-2.0*uxm1(i+1,j+1,1)+0.5*uxm2(i+1,j+1,1)
      dudt4 = 1.5*ux(i  ,j+1,1)-2.0*uxm1(i  ,j+1,1)+0.5*uxm2(i  ,j+1,1)
      duxdt = (dudt1+dudt2+dudt3+dudt4)/4.0/dt
!
      dudt1 = 1.5*uy(i  ,j  ,1)-2.0*uym1(i  ,j  ,1)+0.5*uym2(i  ,j  ,1)
      dudt2 = 1.5*uy(i+1,j  ,1)-2.0*uym1(i+1,j  ,1)+0.5*uym2(i+1,j  ,1)
      dudt3 = 1.5*uy(i+1,j+1,1)-2.0*uym1(i+1,j+1,1)+0.5*uym2(i+1,j+1,1)
      dudt4 = 1.5*uy(i  ,j+1,1)-2.0*uym1(i  ,j+1,1)+0.5*uym2(i  ,j+1,1)
      duydt = (dudt1+dudt2+dudt3+dudt4)/4.0/dt
!
      sumx=sumx+duxdt*dx*(yp(j+1)-yp(j))*(1.-epsi(i,j,k))
      sumy=sumy+duydt*dx*(yp(j+1)-yp(j))*(1.-epsi(i,j,k))
enddo
enddo

do j=jld,jud-1
   uxm = (ux(ild,j,1)+ux(ild,j+1,1))/2.0
   uym = (uy(ild,j,1)+uy(ild,j+1,1))/2.0
   fxa= fxa-uxm*uxm*(yp(j+1)-yp(j))
   fya= fya-uxm*uym*(yp(j+1)-yp(j))
   uxm = (ux(iud,j,1)+ux(iud,j+1,1))/2.0
   uym = (uy(iud,j,1)+uy(iud,j+1,1))/2.0
   fxc= fxc+uxm*uxm*(yp(j+1)-yp(j))
   fyc= fyc+uxm*uym*(yp(j+1)-yp(j))
enddo
do i=ild,iud-1
   uxm = (ux(i,jud,1)+ux(i+1,jud,1))/2.0
   uym = (uy(i,jud,1)+uy(i+1,jud,1))/2.0
   fxb= fxb+uxm*uym*dx
   fyb= fyb+uym*uym*dx
   uxm = (ux(i,jld,1)+ux(i+1,jld,1))/2.0
   uym = (uy(i,jld,1)+uy(i+1,jld,1))/2.0
   fxd= fxd-uxm*uym*dx
   fyd= fyd-uym*uym*dx
enddo

xmom=sumx+fxa+fxb+fxc+fxd
ymom=sumy+fya+fyb+fyc+fyd

return
end subroutine rmomentum2d
!
!***********************************************************************
!
subroutine sforce2d(ux,uy,uz,ppm,gx,gy,gz,tx,ty,tz,fx,fy,fz,di1,di2)
!
!***********************************************************************
!
USE param
USE IBM
USE aeroforce
USE variables
!
implicit none
!
real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,tx,ty,tz,fx,fy,fz,di1,di2
real(8),dimension(nxm,nym,nzm) :: ppm
real(8) :: foxa,foxb,foxc,foxd,foxpr
real(8) :: foya,foyb,foyc,foyd,foypr
real(8) :: dudxm,dudym
real(8) :: dvdxm,dvdym
real(8) :: pab,pcd,pad,pbc
integer :: i,j,k
! 
foxa=0.0
foya=0.0
foxb=0.0
foyb=0.0
foxc=0.0
foyc=0.0
foxd=0.0
foyd=0.0
foxpr=0.0
foypr=0.0
!
call derx (tx,ux,di1    ,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
call dery (ty,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
call derx (gx,uy,di1    ,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
call dery (gy,uy,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
!
do j=jld,jud-1
   dudxm = (tx(ild,j,1)+tx(ild,j+1,1))/2.0
   dudym = (ty(ild,j,1)+ty(ild,j+1,1))/2.0
   dvdxm = (gx(ild,j,1)+gx(ild,j+1,1))/2.0
   foxa = foxa-2.0*xnu*dudxm*(yp(j+1)-yp(j))
   foya = foya-xnu*(dvdxm+dudym)*(yp(j+1)-yp(j))
   dudxm = (tx(iud,j,1)+tx(iud,j+1,1))/2.0
   dudym = (ty(iud,j,1)+ty(iud,j+1,1))/2.0
   dvdxm = (gx(iud,j,1)+gx(iud,j+1,1))/2.0
   foxc = foxc+2.0*xnu*dudxm*(yp(j+1)-yp(j))
   foyc = foyc+xnu*(dvdxm+dudym)*(yp(j+1)-yp(j))
   pab   = ppm(ild,j,1)/gdt(3)
   pcd   = ppm(iud,j,1)/gdt(3)
   foxpr = foxpr+(pab-pcd)*(yp(j+1)-yp(j))
enddo
!
do i=ild,iud-1
   dudym = (ty(i,jud,1)+ty(i+1,jud,1))/2.0
   dvdxm = (gx(i,jud,1)+gx(i+1,jud,1))/2.0
   dvdym = (gy(i,jud,1)+gy(i+1,jud,1))/2.0
   foxb = foxb+xnu*(dudym+dvdxm)*dx
   foyb = foyb+2.0*xnu*dvdym*dx
   dudym = (ty(i,jld,1)+ty(i+1,jld,1))/2.0
   dvdxm = (gx(i,jld,1)+gx(i+1,jld,1))/2.0
   dvdym = (gy(i,jld,1)+gy(i+1,jld,1))/2.0
   foxd = foxd-xnu*(dudym+dvdxm)*dx
   foyd = foyd-2.0*xnu*dvdym*dx
   pad   = ppm(i,jld,1)/gdt(3)
   pbc   = ppm(i,jud,1)/gdt(3)
   foypr = foypr+(pad-pbc)*dx
enddo
!
f2x=foxa+foxb+foxc+foxd+foxpr
f2y=foya+foyb+foyc+foyd+foypr
!
return
end subroutine sforce2d
!
!***********************************************************************
!
subroutine aerofd(ux,uy,uz,pp,epsi,px,py,pz,&
                  tx,ty,tz,di1,di2)
!
!***********************************************************************

USE param
USE IBM
USE aeroforce
USE variables
      
implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,px,py,pz,tx,ty,tz,di1,di2,epsi
real(8),dimension(nxm,nym,nzm) :: pp
integer :: i,j,k,ip
real(8) :: x2,y2,ds,dp,dudx,dudy,dvdx,dvdy,xforcemp,yforcemp,pp1p
!
call derx (tx,ux,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
call dery (ty,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
call derx (px,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
call dery (py,uy,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
!
xforce=0.
yforce=0.
xforcemp=0.
yforcemp=0.
do i=1,lskin
   x2=vertex(i,1)
   y2=vertex(i,2)
   call interpole(tx,epsi,yp,x2,y2,dudx,nx,ny,nz)
   call interpole(ty,epsi,yp,x2,y2,dudy,nx,ny,nz)
   call interpole(px,epsi,yp,x2,y2,dvdx,nx,ny,nz)
   call interpole(py,epsi,yp,x2,y2,dvdy,nx,ny,nz)
   call interpolep(pp,epsi,ypi,x2,y2,dp,nx,ny,nz,nxm,nym,nzm)
   dp=dp/gdt(3)-pp(1,1,1)
   if (i.ne.lskin) ip=i+1
   if (i.eq.lskin) ip=1
   ds=sqrt((body(ip,1)-body(i,1))**2+(body(ip,2)-body(i,2))**2)
   cdcl(i,1)=ds*((-dp+2.*xnu*dudx)*vertex(i,3)+xnu*(dudy+dvdx)*vertex(i,4))
   cdcl(i,2)=ds*((-dp+2.*xnu*dvdy)*vertex(i,4)+xnu*(dudy+dvdx)*vertex(i,3))
   xforcemp=xforcemp+ds*((2.*xnu*dudx)*vertex(i,3)+xnu*(dudy+dvdx)*vertex(i,4))
   yforcemp=yforcemp+ds*((2.*xnu*dvdy)*vertex(i,4)+xnu*(dudy+dvdx)*vertex(i,3))
enddo

open(220,file='force-inst.dat',form='formatted')
do i=1,lskin
   write(220,*) vertex(i,1),vertex(i,2),cdcl(i,1),cdcl(i,2)
   xforce=xforce+cdcl(i,1)
   yforce=yforce+cdcl(i,2)
enddo
close(220)
write(*,102) xforce,yforce
write(217,101) (itime-1)*dt,xforce,yforce,xforcemp,yforcemp
 101 format(7F18.9)
 102 format('Drag=',F8.5,';    Lift=',F8.5)

return
end subroutine aerofd
!********************************************************************
!
subroutine interpole(u,epsi,yp,x2,y2,us,nx,ny,nz)
! 
!********************************************************************

USE param
USE IBM

implicit none

integer :: nx,ny,nz
real(8),dimension(nx,ny,nz) :: u,epsi
real(8),dimension(ny) :: yp
real(8) :: x2,y2,us
real(8) :: ubr,ubl,utr,utl,ust,usb
real(8) :: xs1,xs2,ys1,ys2
integer :: i2,j

i2   =int(x2/dx)+1
do j=1,ny-1
   if (y2.gt.yp(j).and.y2.lt.yp(j+1)) j2=j
enddo
xs1=abs(x2-(i2-2)*dx)
xs2=abs(x2-(i2+1)*dx)
ys1=abs(y2-yp(j2-1))
ys2=abs(y2-yp(j2+2))
ubr=u(i2+1,j2  ,1)*(1.-epsi(i2+1,j2  ,1))
utr=u(i2+1,j2+1,1)*(1.-epsi(i2+1,j2+1,1))
ubl=u(i2  ,j2  ,1)*(1.-epsi(i2  ,j2  ,1))
utl=u(i2  ,j2+1,1)*(1.-epsi(i2  ,j2+1,1))
!ubr=u(i2+2,j2-1,1)
!utr=u(i2+2,j2+2,1)
!ubl=u(i2-1,j2-1,1)
!utl=u(i2-1,j2+2,1)

ust=(utl*xs2+utr*xs1)/(3.*dx)
usb=(ubl*xs2+ubr*xs1)/(3.*dx)
us=(usb*ys2+ust*ys1)/(yp(j2+2)-yp(j2-1))

return
end subroutine interpole
!********************************************************************
!
subroutine interpolep(u,epsi,yp,x2,y2,us,nx,ny,nz,nxm,nym,nzm)
! 
!********************************************************************

USE param
USE IBM

implicit none
integer :: nxm,nym,nzm,nx,ny,nz
real(8),dimension(nxm,nym,nzm) :: u
real(8),dimension(nx,ny,nz) :: epsi
real(8),dimension(ny) :: yp
real(8) :: x2,y2,us
real(8) :: ubr,ubl,utr,utl,ust,usb
real(8) :: xs1,xs2,ys1,ys2
integer :: i2,j

i2   =int(x2/dx)+1
do j=1,ny-1
   if (y2.gt.yp(j).and.y2.lt.yp(j+1)) j2=j
enddo
xs1=abs(x2-((i2-2)*dx))
xs2=abs(x2-((i2+1)*dx))
ys1=abs(y2-yp(j2-1))
ys2=abs(y2-yp(j2+2))
ubr=u(i2+1,j2  ,1)*(1.-epsi(i2+1,j2  ,1))
utr=u(i2+1,j2+1,1)*(1.-epsi(i2+1,j2+1,1))
ubl=u(i2  ,j2  ,1)*(1.-epsi(i2  ,j2  ,1))
utl=u(i2  ,j2+1,1)*(1.-epsi(i2  ,j2+1,1))
!ubr=u(i2+2,j2-1,1)
!utr=u(i2+2,j2+2,1)
!ubl=u(i2-1,j2-1,1)
!utl=u(i2-1,j2+2,1)

ust=(utl*xs2+utr*xs1)/(3.*dx)
usb=(ubl*xs2+ubr*xs1)/(3.*dx)
us =(usb*ys2+ust*ys1)/(yp(j2+2)-yp(j2-1))

return
end subroutine interpolep
end module forces_m
