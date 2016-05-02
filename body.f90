module body_m
  use navier_m
implicit none
contains
!*******************************************************************
!
subroutine solid_body(ux,uy,uz,epsi,ppm,sy1,sy2,sy3,sy4,sy5,sy6,sy7,sy8,sy9,di1,di2)
!
!*******************************************************************

USE param
USE IBM 
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,epsi
real(8),dimension(nx,ny,nz) :: sy1,sy2,sy3
real(8),dimension(nx,ny,nz) :: sy4,sy5,sy6
real(8),dimension(nx,ny,nz) :: sy7,sy8,sy9
real(8),dimension(nxm,nym,nz) :: ppm
real(8),dimension(mx,my,mz) :: di2,di1 ! TODO

   call gradpression(ppm,sy1,sy2,sy3,sy4,sy5,sy6,sy7,sy8,sy9,di1,di2)
   call corgp_IBM(ux,uy,uz,sy1,sy2,sy3,1)
   call forcage_original(ux,uy,uz,epsi)
   call corgp_IBM(ux,uy,uz,sy1,sy2,sy3,2)

return
end subroutine solid_body

!********************************************************************
!
subroutine corgp_IBM (ux,uy,uz,px,py,pz,nlock)
! 
!********************************************************************

USE param
USE IBM
USE variables

implicit none

integer :: ijk,nlock,nxyz,i,j,k
real(8),dimension(nx,ny,nz) :: ux,uy,uz,px,py,pz

nxyz=nx*ny*nz

if (nlock.eq.1) then
   if (nz.gt.1) then
      do ijk=1,nxyz
         uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1) 
         uz(ijk,1,1)=-pz(ijk,1,1)+uz(ijk,1,1) 
         ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
      enddo
   else
!      do ijk=1,nxyz ! TODO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         uy(i,j,k)=-py(i,j,k)+uy(i,j,k) 
         ux(i,j,k)=-px(i,j,k)+ux(i,j,k)
      enddo
      enddo
      enddo      
   endif
endif
if (nlock.eq.2) then
   if (nz.gt.1) then
      do ijk=1,nxyz
         uy(ijk,1,1)=py(ijk,1,1)+uy(ijk,1,1) 
         uz(ijk,1,1)=pz(ijk,1,1)+uz(ijk,1,1) 
         ux(ijk,1,1)=px(ijk,1,1)+ux(ijk,1,1)
      enddo
   else
!      do ijk=1,nxyz ! TODO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         uy(i,j,k)=py(i,j,k)+uy(i,j,k) 
         ux(i,j,k)=px(i,j,k)+ux(i,j,k)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine corgp_IBM

!*******************************************************************
!
subroutine forcage1(ux,uy,epsi)
! 
!********************************************************************

USE param
USE IBM
USE variables

  implicit none

 real(8),dimension(nx,ny,nz) :: ux,uy,epsi
 real(8) :: r,x,y
 integer :: i,j,k

   epsi(:,:,:)=0.
   do j=1,ny
   do i=1,nx
      y=(j-1)*dy
      x=(i-1)*dx
      if(sqrt((x-cex)**2+(y-cey)**2).le.ra)then
      ux(i,j,1)=0.
      uy(i,j,1)=0.
      epsi(i,j,1)=1.
      endif
   enddo
   enddo
 
   return
 end subroutine forcage1

!*******************************************************************
!
subroutine forcage_original(ux,uy,uz,epsi)
!
!*******************************************************************

USE param
USE IBM 
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,epsi
real(8),dimension(nz) :: uxs,uys,uzs,e
integer :: j, i, k, np, i2, iii, jjj,ngxsym,ngysym,ndxsym,ndysym
real(8) :: r,xm,ym,zero,teta,rsym,xsym,ysym,uuu,ttt
real(8) , external :: epmach

pi=acos(-1.D0)
zero=0.D0
   do j=1,ny 
   do i=1,nx
      xm=(i-1)*dx 
      if (istret.ne.0) ym=yp(j)
      if (istret.eq.0) ym=(j-1)*dy 
!
      r=sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey)) 
      if ((r-ra).lt.zero) then
         if (nz.gt.1) then
            do k=1,nz
               ux(i,j,k)=0.D0
               uy(i,j,k)=0.D0 
               uz(i,j,k)=0.D0 
            enddo
         else
            k=1 
            ux(i,j,k)=0.D0
            uy(i,j,k)=0.D0
         endif
      endif
   enddo
   enddo


if (nz.gt.1) then
   np=0 
   epsi=0.
   do j=1,ny 
   do i=1,nx
      xm=(i-1)*dx 
      if (istret.ne.0) ym=yp(j)
      if (istret.eq.0) ym=(j-1)*dy 

      r=sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey)) 
      if ((r-ra).lt.zero) then
          if (abs(xm-cex).lt.(zero)) then 
                  xm=cex 
          endif
          if (abs(ym-cey).lt.(zero)) then 
                  ym=cey
          endif
          if ((xm-cex).eq.0.) then
               if ((ym-cey).eq.0.) then
                    teta=0.
               else
                    teta=(ym-cey)/(abs(ym-cey))*0.5*pi
               endif
          endif
          if ((xm-cex).gt.0.) then
               teta=atan((ym-cey)/(xm-cex))
          endif
          if ((xm-cex).lt.0.) then
               teta=atan((ym-cey)/(xm-cex))+pi
          endif
          rsym=ra+(ra-r)
          xsym=rsym*cos(teta)+cex
          ysym=rsym*sin(teta)+cey         
          iii=1
          do while (((iii-1)*dx) < xsym)
               iii=iii+1
          enddo
          iii=iii-1
          jjj=1
          do while (yp(jjj) < ysym)
              jjj=jjj+1
          enddo
          jjj=jjj-1                  
          ngxsym=iii
          ngysym=jjj
          ndxsym=ngxsym+1
          ndysym=ngysym+1          
!
          ttt=(xsym-((ngxsym-1)*dx))/dx
          uuu=(ysym-yp(ngysym))/(yp(ndysym)-yp(ngysym))
!
          do k=1,nz                    
!
!*******************************************************************************
! Modele miroir u_t=-u_t;u_n=0     **********************************************
! *******************************************************************************

               uxs(k)=((1-ttt)*(1-uuu)*ux(ngxsym,ngysym,k)+&
                         &ttt*(1-uuu)*ux(ndxsym,ngysym,k)+&
                         &ttt*uuu*ux(ndxsym,ndysym,k)+& 
                         &(1-ttt)*uuu*ux(ngxsym,ndysym,k))
!       
               uys(k)=((1-ttt)*(1-uuu)*uy(ngxsym,ngysym,k)+& 
                         &ttt*(1-uuu)*uy(ndxsym,ngysym,k)+&
                         &ttt*uuu*uy(ndxsym,ndysym,k)+&
                         &(1-ttt)*uuu*uy(ngxsym,ndysym,k))
!
               uzs(k)=((1-ttt)*(1-uuu)*uz(ngxsym,ngysym,k)+& 
                          &ttt*(1-uuu)*uz(ndxsym,ngysym,k)+&
                          &ttt*uuu*uz(ndxsym,ndysym,k)+&
                          &(1-ttt)*uuu*uz(ngxsym,ndysym,k)) 
!
               e(k)=sin(2*pi*r**2)
               ux(i,j,k)=e(k)*(-uxs(k)*sin(teta)**2+&
                               &uys(k)*cos(teta)*sin(teta))
               uy(i,j,k)=e(k)*(uxs(k)*sin(teta)*cos(teta)-&
                               &uys(k)*cos(teta)**2)
               uz(i,j,k)=-e(k)*uzs(k)
!                                   
! *******************************************************************************         
! ******************************************************************************* 
! Modele miroir u_t= -u_t ; u_n=-u_n ********************************************
! *******************************************************************************
!                     ux(i,j,k)=-e(k)*uxs
!                     uy(i,j,k)=-e(k)*uys
!                     uz(i,j,k)=-e(k)*uzs
! 
               np=np+1 
               epsi(i,j,k)=1. 
          enddo
      endif
   enddo
   enddo
else
   np=0
   epsi=0.
   if (itime == 1) print *,"forcage cex,cey,ra :",cex,'',cey,'',ra
   do j=1,ny
   do i=1,nx
      xm=(i-1)*dx
      if (istret.ne.0) ym=yp(j)
      if (istret.eq.0) ym=(j-1)*dy 
!
      r=sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey))
      if ((r-ra).lt.zero) then
          if (abs(xm-cex).lt.(zero)) then 
                  xm=cex 
          endif
          if (abs(ym-cey).lt.(zero)) then 
                  ym=cey
          endif
          if ((xm-cex).eq.0.) then
               if ((ym-cey).eq.0.) then
                    teta=0.
               else
                    teta=(ym-cey)/(abs(ym-cey))*0.5*pi
               endif
          endif
          if ((xm-cex).gt.0.) then
               teta=atan((ym-cey)/(xm-cex))
          endif
          if ((xm-cex).lt.0.) then
               teta=atan((ym-cey)/(xm-cex))+pi
          endif
          rsym=ra+(ra-r)
          xsym=rsym*cos(teta)+cex
          ysym=rsym*sin(teta)+cey
!
          iii=1
          do while (((iii-1)*dx) < xsym)
               iii=iii+1
          enddo
          iii=iii-1
          jjj=1
          do while (yp(jjj) < ysym)
              jjj=jjj+1
          enddo
          jjj=jjj-1                  
          ngxsym=iii
          ngysym=jjj
          ndxsym=ngxsym+1
          ndysym=ngysym+1
!
          ttt=(xsym-((ngxsym-1)*dx))/dx
          uuu=(ysym-yp(ngysym))/(yp(ndysym)-yp(ngysym))
!
          k=1!
          uxs(k)=((1-ttt)*(1-uuu)*ux(ngxsym,ngysym,k)+&
                    &ttt*(1-uuu)*ux(ndxsym,ngysym,k)+&
                    &ttt*uuu*ux(ndxsym,ndysym,k)+& 
                    &(1-ttt)*uuu*ux(ngxsym,ndysym,k))
       
          uys(k)=((1-ttt)*(1-uuu)*uy(ngxsym,ngysym,k)+& 
                    &ttt*(1-uuu)*uy(ndxsym,ngysym,k)+&
                    &ttt*uuu*uy(ndxsym,ndysym,k)+&
                    &(1-ttt)*uuu*uy(ngxsym,ndysym,k))
!
!*******************************************************************************
! Modele miroir u_t=-u_t;u_n=0     **********************************************
! *******************************************************************************
          e(k)=sin(2*pi*r**2/(2.*ra)**2)
          ux(i,j,k)=e(k)*(-uxs(k)*sin(teta)**2+&
                   &uys(k)*cos(teta)*sin(teta))
          uy(i,j,k)=e(k)*(uxs(k)*sin(teta)*cos(teta)-&
                   &uys(k)*cos(teta)**2)
!                                   
! *******************************************************************************         
! ******************************************************************************* 
! Modele miroir u_t= -u_t ; u_n=-u_n ********************************************
! *******************************************************************************
!                     ux(i,j,k)=-e(k)*uxs
!                     uy(i,j,k)=-e(k)*uys
!                     uz(i,j,k)=-e(k)*uzs
          epsi(i,j,k)=1. 
       endif
   enddo
   enddo
endif
print *,"np =",np
!

return  
end subroutine forcage_original

!*******************************************************************
!
subroutine forcage(ux,uy,uz,epsi,lme)
!
!*******************************************************************

USE param
USE IBM 
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,epsi
integer,dimension(nx,ny,nz) :: lme

utrans=0.
vtrans=0.
dAoAdt=0.
AoA=0.

call shape()
call movebody()          
call vertices()
call epsilon(ux,uy,uz,epsi,lme)
call imagebi(ux,uy,uz,epsi,lme)

return
end subroutine forcage

!********************************************************************
!
subroutine shape()
! 
!********************************************************************

USE IBM

implicit none

real(8) :: dalpha,a,b
integer :: i

ra=0.5
a=ra
b=0.05 !! ELLIPSOID
b=ra   !! CYLINDER

dalpha=2.*acos(-1.)/float(lskin)

body(1,1)=ra
body(1,2)=0.
do i=2,lskin
   body(i,1)=cos((i-1)*dalpha)*ra
   body(i,2)=sin((i-1)*dalpha)*ra
enddo
open(20,file='body.dat',form='formatted')
do i=1,lskin
   write(20,*) body(i,1),body(i,2)
enddo

return
end subroutine shape

!********************************************************************
!
subroutine movebody()
! 
!********************************************************************

USE IBM
USE param

implicit none

real(8) :: x1,y1,theta,d
integer :: i

cx=0.
cy=0.
do i=1,lskin
   x1=body(i,1)
   y1=body(i,2)
   body(i,1)=(x1-cx)*cos(AoA)-(y1-cy)*sin(AoA)+cex
   body(i,2)=(x1-cx)*sin(AoA)+(y1-cy)*cos(AoA)+cey
   if ((x1-cx).ne.0) theta=atan((y1-cy)/(x1-cx))
   d=sqrt((x1-cx)**2+(y1-cy)**2)
   if ((x1-cx).lt.0) theta=theta+pi
   urot(i,1)=-d*dAoAdt*sin(AoA+theta)+utrans
   urot(i,2)=+d*dAoAdt*cos(AoA+theta)+vtrans
enddo

return
end subroutine movebody

!********************************************************************
!
subroutine vertices()
! 
!********************************************************************

USE IBM

implicit none

real(8) :: x1,x2,y1,y2,u1,u2,v1,v2,theta
integer :: i,ip
real,parameter :: pi=3.141592654

do i=1,lskin
   ip=i+1
   if (i.eq.lskin) ip=1
   u1=urot(i ,1)
   v1=urot(i ,2)
   u2=urot(ip,1)
   v2=urot(ip,2)
   x1=body(i,1)
   y1=body(i,2)
   x2=body(ip,1)
   y2=body(ip,2)
   vertex(i,1)=(x1+x2)/2.
   vertex(i,2)=(y1+y2)/2.
   if (x1.ne.x2) then 
      theta=atan((y2-y1)/(x2-x1))
      if (x2-x1.lt.0) theta=theta+pi
   else
      if (y1.lt.y2) theta= pi/2.
      if (y1.gt.y2) theta=-pi/2.
   endif
   vertex(i,3)=-sin(theta)
   vertex(i,4)=+cos(theta)
   vertex(i,5)=(u1+u2)/2.
   vertex(i,6)=(v1+v2)/2.
enddo
!

return
end subroutine vertices

!********************************************************************
!
   subroutine epsilon(ux,uy,uz,epsi,lme)
! 
!********************************************************************

USE IBM
USE param
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: epsi,ux,uy,uz
integer,dimension(nx,ny,nz) :: lme
real(8) :: x1,x2,y1,y2,d,dmin,n1,n2,p1,p2,xsum
integer :: i,j,k,l,lmin,nepsi,nghost

nepsi=0
do k=1,nz
do j=1,ny
do i=1,nx
   x1=(i-1)*dx
   y1=yp(j)
   dmin =10000.
   do l=1,lskin
      x2=vertex(l,1)
      y2=vertex(l,2)
      d =sqrt((x2-x1 )**2+(y2-y1 )**2)
      if (d.lt.dmin) then
         dmin=d
         lmin=l
      endif
   enddo
   lme(i,j,k)=lmin
   n1=vertex(lmin,3)
   n2=vertex(lmin,4)
   p1=x1-vertex(lmin,1)
   p2=y1-vertex(lmin,2)
   dmin=(p1*n1+p2*n2+0.000001)/abs(p1*n1+p2*n2+0.000001)
   epsi(i,j,k)=max(0.,dmin)
   ux(i,j,k)=epsi(i,j,k)*vertex(lmin,5)+(1.-epsi(i,j,k))*ux(i,j,k)
   uy(i,j,k)=epsi(i,j,k)*vertex(lmin,6)+(1.-epsi(i,j,k))*uy(i,j,k)
   uz(i,j,k)=0.
   nepsi=nepsi+int(epsi(i,j,k))
enddo
enddo
enddo

return
end subroutine epsilon

!********************************************************************
!
subroutine imagebi(ux,uy,uz,epsi,lme)
! 
!********************************************************************

USE IBM
USE param
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: epsi,ux,uy,uz
integer,dimension(nx,ny,nz) :: lme
real(8),dimension(4) :: x,y,u,v
real(8) :: xim,yim,ui,vi
real(8) :: x2,y2,d,dmin,n1,n2,xsum,dmax,x1,y1
integer :: i,j,k,l,lmin,nepsi,nclearcells,lminp1,m
real(8) :: px1,px2,py1,py2,a,up1,up2,p1,p2,bix,biy,xl1,yl1,bx,by
integer :: iu1,ju1,nghost

if (mod(itime,100).eq.0) then
open(20,file='object.dat')
do l=1,lskin
   write(20,*) body(l,1),body(l,2)
enddo
write(20,*) body(1,1),body(1,2)
close(20)
open(204,file='gc.dat')
open(205,file='bi.dat')
open(206,file='ip.dat')
open(207,file='np.dat')
endif
do j=2,ny-1
do i=2,nx-1
   xsum=(epsi(i-1,j-1,1)+epsi(i-1,j,1)+epsi(i-1,j+1,1)+epsi(i,j+1,1)+&
        &epsi(i+1,j+1,1)+epsi(i+1,j,1)+epsi(i+1,j-1,1)+epsi(i,j-1,1))*epsi(i,j,1)
   if (xsum.gt.0.and.xsum.lt.8) then
      nghost=nghost+1
      x1=(i-1)*dx
      y1=yp(j)
      call normal(lme,x1,y1,bx,by,i,j)

      p1=2*bx-x1
      p2=2*by-y1
      iu1=int(p1/dx+1)
      do jj=1,ny-1
         if (p2.gt.yp(jj).and.p2.lt.yp(jj+1)) ju1=jj
      enddo
! Velocity at the four corners close to the image points
! If one corner inside, then velocity at the closest bi
      xl1=(iu1-1)*dx
      yl1=yp(ju1)
      call normal(lme,xl1,yl1,bix,biy,iu1,ju1)
      x(1)=(1.-epsi(iu1  ,ju1  ,1))*xl1+epsi(iu1  ,ju1  ,1)*bix
      y(1)=(1.-epsi(iu1  ,ju1  ,1))*yl1+epsi(iu1  ,ju1  ,1)*biy
!
      xl1=(iu1-1)*dx
      yl1=yp(ju1+1)
      call normal(lme,xl1,yl1,bix,biy,iu1,ju1+1)
      x(2)=(1.-epsi(iu1  ,ju1+1,1))*xl1+epsi(iu1  ,ju1+1,1)*bix
      y(2)=(1.-epsi(iu1  ,ju1+1,1))*yl1+epsi(iu1  ,ju1+1,1)*biy
!
      xl1=(iu1)*dx
      yl1=yp(ju1+1)
      call normal(lme,xl1,yl1,bix,biy,iu1+1,ju1+1)
      x(3)=(1.-epsi(iu1+1,ju1+1,1))*xl1+epsi(iu1+1,ju1+1,1)*bix
      y(3)=(1.-epsi(iu1+1,ju1+1,1))*yl1+epsi(iu1+1,ju1+1,1)*biy
!
      xl1=(iu1)*dx
      yl1=yp(ju1)
      call normal(lme,xl1,yl1,bix,biy,iu1+1,ju1)
      x(4)=(1.-epsi(iu1+1,ju1  ,1))*xl1+epsi(iu1+1,ju1  ,1)*bix
      y(4)=(1.-epsi(iu1+1,ju1  ,1))*yl1+epsi(iu1+1,ju1  ,1)*biy
!
      u(1)=ux(iu1  ,ju1  ,1)
      u(2)=ux(iu1  ,ju1+1,1)
      u(3)=ux(iu1+1,ju1+1,1)
      u(4)=ux(iu1+1,ju1  ,1)
!
      v(1)=uy(iu1  ,ju1  ,1)
      v(2)=uy(iu1  ,ju1+1,1)
      v(3)=uy(iu1+1,ju1+1,1)
      v(4)=uy(iu1+1,ju1  ,1)
!
      if (mod(itime,100).eq.0) then
      write(207,1102) (x(m),m=1,4),(y(m),m=1,4),(u(m),m=1,4),(v(m),m=1,4)
      endif
      xim=p1
      yim=p2
      call bilinear(u,v,x,y,ui,vi,xim,yim)
      ux(i,j,1)=2*vertex(lme(i,j,1),5)-ui
      uy(i,j,1)=2*vertex(lme(i,j,1),6)-vi

      if (mod(itime,100).eq.0) then
      lmin=lme(i,j,1) ! TODO
      write(204,1101) x1,y1,ux(i,j,1),uy(i,j,1)
      write(205,1101) bx,by,vertex(lmin,5),vertex(lmin,6)
      write(206,1101) p1,p2,ui,vi
      endif
   endif
enddo
enddo
if (mod(itime,100).eq.0) then
close(204)
close(205)
close(206)
close(207)
endif
1101 format(6F18.6)
1102 format(16F18.6)
return
end subroutine imagebi
!********************************************************************
!
subroutine normal(lme,x1,y1,bix,biy,i,j)
! 
!********************************************************************

USE IBM
USE param
USE variables

implicit none

integer,dimension(nx,ny,nz) :: lme
real(8) :: px1,py1,px2,py2,a,bix,biy,x1,y1
integer :: i,j,lmin,lminp1

lmin=lme(i,j,1)
lminp1=lmin+1
if (lmin.eq.lskin) lminp1=1
px1=body(lmin,1)
py1=body(lmin,2)
px2=body(lminp1,1)
py2=body(lminp1,2)
a=((x1-px1)*(px2-px1)+(y1-py1)*(py2-py1))/((px2-px1)**2+(py2-py1)**2)
bix=px1+a*(px2-px1)
biy=py1+a*(py2-py1)

return
end subroutine normal

!********************************************************************
!
subroutine checkimage(ux,uy,uz,epsi,lme)
! 
!********************************************************************

USE IBM
USE param
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: epsi,ux,uy,uz
integer,dimension(nx,ny,nz) :: lme
real(8),dimension(4) :: x,y,u,v
real(8) :: xim,yim,ui,vi
real(8) :: x2,y2,d,dmin,n1,n2,xsum,dmax,x1,y1
integer :: i,j,k,l,lmin,nepsi,nclearcells,lminp1,m
real(8) :: px1,px2,py1,py2,a,up1,up2,p1,p2,bix,biy,xl1,yl1,bx,by
integer :: iu1,ju1,nghost

if (mod(itime,100).eq.0) then
open(205,file='data/bi-check.dat')
open(207,file='data/np-check.dat')
endif
do j=2,ny-1
do i=2,nx-1
   xsum=(epsi(i-1,j-1,1)+epsi(i-1,j,1)+epsi(i-1,j+1,1)+epsi(i,j+1,1)+&
        &epsi(i+1,j+1,1)+epsi(i+1,j,1)+epsi(i+1,j-1,1)+epsi(i,j-1,1))*epsi(i,j,1)
   if (xsum.gt.0.and.xsum.lt.8) then
      nghost=nghost+1
      x1=(i-1)*dx
      y1=yp(j)
      call normal(lme,x1,y1,bx,by,i,j)

      p1=bx
      p2=by
      iu1=int(p1/dx+1)
      do jj=1,ny-1
         if (p2.gt.yp(jj).and.p2.lt.yp(jj+1)) ju1=jj
      enddo
!
      x(1)=(iu1-1)*dx
      x(2)=(iu1-1)*dx
      x(3)=(iu1  )*dx
      x(4)=(iu1  )*dx
!
      y(1)=yp(ju1)
      y(2)=yp(ju1+1)
      y(3)=yp(ju1+1)
      y(4)=yp(ju1)
!
      u(1)=ux(iu1  ,ju1  ,1)
      u(2)=ux(iu1  ,ju1+1,1)
      u(3)=ux(iu1+1,ju1+1,1)
      u(4)=ux(iu1+1,ju1  ,1)
!
      v(1)=uy(iu1  ,ju1  ,1)
      v(2)=uy(iu1  ,ju1+1,1)
      v(3)=uy(iu1+1,ju1+1,1)
      v(4)=uy(iu1+1,ju1  ,1)
!
      if (mod(itime,100).eq.0) then
      write(207,1102) (x(m),m=1,4),(y(m),m=1,4),(u(m),m=1,4),(v(m),m=1,4)
      endif
      xim=bx
      yim=by
      call bilinear(u,v,x,y,ui,vi,xim,yim)

      if (mod(itime,100).eq.0) then
         write(205,1101) bx,by,ui,vi
      endif
   endif
enddo
enddo
if (mod(itime,100).eq.0) then
close(204)
close(205)
close(206)
endif
1101 format(6F18.6)
1102 format(16F18.6)
return
end subroutine checkimage

subroutine bilinear(u,v,x,y,u1,v1,x1,y1)
!
implicit none
!
real(8),dimension(4) :: x,y,u,v
real(8),dimension(4,4) :: a
real(8) :: u1,v1,x1,y1,d
integer,dimension(4) :: indx
integer :: m

!               
do m=1,4
   a(m,1) = x(m)*y(m)
   a(m,2) = x(m)
   a(m,3) = y(m)
   a(m,4) = 1.0
enddo

call ludcmp(a,4,4,indx,d)

call lubksb(a,4,4,indx,u)
u1 = u(1)*x1*y1 + u(2)*x1 + u(3)*y1 + u(4)
call lubksb(a,4,4,indx,v)
v1 = v(1)*x1*y1 + v(2)*x1 + v(3)*y1 + v(4)
      
return
end subroutine bilinear

SUBROUTINE ludcmp(a,n,np,indx,d)
integer :: np,n
real(8),dimension(np,np) :: a
real(8),dimension(30) :: vv
integer,dimension(n) :: indx
real(8) :: d,TINY,aamax,dum,sum
integer :: i,imax,j,k
TINY=1.0e-20

d=1.
do i=1,n
   aamax=0.
   do j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
   enddo
   if (aamax.eq.0.) then 
      stop 'singular matrix in ludcmp'
   endif
   vv(i)=1./aamax
enddo

do j=1,n
   do i=1,j-1
      sum=a(i,j)
      do k=1,i-1
         sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
   enddo
   aamax=0.
   do i=j,n
      sum=a(i,j)
      do k=1,j-1
         sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      dum=vv(i)*abs(sum)
      if (dum.ge.aamax) then
         imax=i
         aamax=dum
      endif
   enddo
   if (j.ne.imax)then
      do k=1,n
         dum=a(imax,k)
         a(imax,k)=a(j,k)
         a(j,k)=dum
      enddo
      d=-d
      vv(imax)=vv(j)
   endif
   indx(j)=imax
   if (a(j,j).eq.0.)a(j,j)=TINY
   if (j.ne.n)then
       dum=1./a(j,j)
       do i=j+1,n
          a(i,j)=a(i,j)*dum
       enddo
   endif
enddo

return
end subroutine ludcmp

SUBROUTINE lubksb(a,n,np,indx,b)
!
implicit none
!
INTEGER n,np
real(8),dimension(np,np) :: a
real(8),dimension(n) :: b
integer,dimension(n) :: indx
INTEGER i,ii,j,ll
real(8) :: sum
!
ii=0
do i=1,n
   ll=indx(i)
   sum=b(ll)
   b(ll)=b(i)
   if (ii.ne.0)then
      do j=ii,i-1
         sum=sum-a(i,j)*b(j)
      enddo
   else if (sum.ne.0.) then
      ii=i
   endif
   b(i)=sum
enddo
!
do i=n,1,-1
   sum=b(i)
   do j=i+1,n
      sum=sum-a(i,j)*b(j)
   enddo
   b(i)=sum/a(i,i)
enddo
!
return
end subroutine lubksb

end module body_m
