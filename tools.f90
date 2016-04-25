!******************************************************************
!
subroutine numcar (num,car)
!
!******************************************************************!

character(len=3) car

if (num.ge.100) then
   write (car,1) num
1  format (i3)
else 
   if (num.ge.10) then
      write (car,2) num
2     format ('0',i2)
   else
      write (car,3) num
3     format ('00',i1)
   endif
endif

return
end subroutine numcar

!****************************************************
!
real(8) function rand2 (idum) 
!
! ****************************************************

implicit none

integer :: idum 
integer,parameter :: im1 = 2147483563 
integer,parameter :: im2 = 2147483399 
integer,parameter :: imm1 = im1 - 1 
integer,parameter :: ia1 = 40014 
integer,parameter :: ia2 = 40692 
integer,parameter :: iq1 = 53668 
integer,parameter :: iq2 = 52774 
integer,parameter :: ir1 = 12211 
integer,parameter :: ir2 = 3791 
integer,parameter :: ntab = 32 
integer,parameter :: ndiv = 1 + imm1/ntab 
real(8),parameter :: am = 1./im1 
real(8),parameter :: eps = 1.2E-14 
real(8),parameter :: rnmx = 1. - eps 
integer :: idum2, j, k 
integer,dimension(ntab) :: iv 
integer :: iy 
save idum2, iv, iy 

data idum2/123456789/ 
data iv/ntab*0/ 
data iy/0/ 

if (idum <= 0) then                        !Initialize. 
   idum=max((-idum),1) 
   idum2=idum 
   do j=ntab + 8,1,-1 
      k=idum/iq1 
      idum=ia1*(idum - k*iq1) - k*ir1 
      if (idum < 0) idum=idum + im1 
      if (j > ntab) cycle  
      iv(j)=idum 
   enddo
endif

k=idum/iq1 
idum=ia1*(idum - k*iq1) - k*ir1 
if (idum < 0) idum=idum + im1 

k=idum2/iq2 
idum2=ia2*(idum2 - k*iq2) - k*ir2 
if (idum2 < 0) idum2=idum2 + im2 

j=1 + iy/ndiv 
iy=iv(j) - idum2 
iv(j)=idum 
if (iy < 1) iy=iy + imm1 
!===========================
rand2=min(am*iy,rnmx) 
!===========================

return  
end function rand2

!************************************************************************************
!
subroutine stretching(yp,ypi,ppy,ppyi,pp2y,pp4y,pp2yi,pp4yi,ny)
!
!************************************************************************************

USE param

implicit none

integer :: ny,j
real(8),dimension(ny) :: ppy,pp2y,pp4y
real(8),dimension(ny) :: ppyi,pp2yi,pp4yi
real(8),dimension(ny) :: yp,ypi,yeta,yetai
real(8) :: yinf,den,xnum,xcx,den1,den3,den4,xnum1,cst,yl2y

yl2y=yly*yly
yinf=-yly/2.
den=2.*beta*yinf
xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
alpha=abs(xnum/den)

xcx=1./beta/alpha

if (alpha.ne.0.) then
   if (istret.eq.1) yp(1)=0.
   if (istret.eq.2) yp(1)=0.
   if (istret.eq.1) yeta(1)=0.
   if (istret.eq.2) yeta(1)=-0.5	
   if (istret.eq.3) yp(1)=0.
   if (istret.eq.3) yeta(1)=-0.5
   do j=2,ny
      if (istret==1) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yeta(j)=((j-1.)*(1./2./ny)-0.5)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=((j-1.)*(1./2./(ny-1.))-0.5)
      endif
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yeta(j))+1.
      xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst-yinf
         if (yeta(j).eq.0.5) yp(j)=0.-yinf
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst+yly
         if (yeta(j).eq.0.5) yp(j)=0.+yly
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yeta(j).lt.0.5) yp(j)=(xnum1-cst+yly)*2.
         if (yeta(j).eq.0.5) yp(j)=(0.+yly)*2.
         if (yeta(j).gt.0.5) yp(j)=(xnum1+cst+yly)*2.
      endif
   enddo
endif
if (alpha.eq.0.) then
   yp(1)=-1.e10
   do j=2,ny
      yeta(j)=(j-1.)*(1./ny)
      yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
   enddo
endif
if (alpha.ne.0.) then
   do j=1,ny
      if (istret==1) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./2./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./2./(ny-1.))-0.5
      endif
      
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yetai(j))+1.
      xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst-yinf
         if (yetai(j).eq.0.5) ypi(j)=0.-yinf
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst+yly
         if (yetai(j).eq.0.5) ypi(j)=0.+yly
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yetai(j).lt.0.5) ypi(j)=(xnum1-cst+yly)*2.
         if (yetai(j).eq.0.5) ypi(j)=(0.+yly)*2.
         if (yetai(j).gt.0.5) ypi(j)=(xnum1+cst+yly)*2.
      endif
   enddo
endif
if (alpha.eq.0.) then
   ypi(1)=-1.e10
   do j=2,ny
      yetai(j)=(j-1.)*(1./ny)
      ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
   enddo
endif
!print *, 'min',yp(ny)-yp(ny-1),'=?',dx

do j=1,ny
   ppy(j)=yly*(alpha/pi+(1./pi/beta)*sin(pi*yeta(j))* &
        sin(pi*yeta(j)))
   pp2y(j)=ppy(j)*ppy(j)
   pp4y(j)=(-2./beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
enddo
do j=1,ny
   ppyi(j)=yly*(alpha/pi+(1./pi/beta)*sin(pi*yetai(j))* &
        sin(pi*yetai(j)))
   pp2yi(j)=ppyi(j)*ppyi(j)
   pp4yi(j)=(-2./beta*cos(pi*yetai(j))*sin(pi*yetai(j)))
enddo



     
return
end subroutine stretching

!****************************************************************************
!
subroutine pre_correc(ux,uy,uz)
!
!****************************************************************************

USE param
USE IBM
USE variables

implicit none

integer :: i,j,k
real(8) :: utin,utout
real(8),dimension(nx,ny,nz) :: ux,uy,uz
real(8) :: x

if ((nclx==2).or.(ncly==2).or.(nclz==2)) then
   do j=1,ny
   do i=1,nx
      dpdxz1(i,j)=dpdxz1(i,j)*gdt(itr)
      dpdxzn(i,j)=dpdxzn(i,j)*gdt(itr)
      dpdyz1(i,j)=dpdyz1(i,j)*gdt(itr)
      dpdyzn(i,j)=dpdyzn(i,j)*gdt(itr)
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      dpdxy1(i,k)=dpdxy1(i,k)*gdt(itr)
      dpdxyn(i,k)=dpdxyn(i,k)*gdt(itr)
      dpdzy1(i,k)=dpdzy1(i,k)*gdt(itr)
      dpdzyn(i,k)=dpdzyn(i,k)*gdt(itr)
   enddo
   enddo
   do k=1,nz
   do j=1,ny
      dpdyx1(j,k)=dpdyx1(j,k)*gdt(itr)
      dpdyxn(j,k)=dpdyxn(j,k)*gdt(itr)
      dpdzx1(j,k)=dpdzx1(j,k)*gdt(itr)
      dpdzxn(j,k)=dpdzxn(j,k)*gdt(itr)
   enddo
   enddo
   if (nz.eq.1) then
      do k=1,nz
      do i=1,nx
         dpdzy1(i,k)=0.
         dpdzyn(i,k)=0.
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         dpdzx1(j,k)=0.
         dpdzxn(j,k)=0.
      enddo
      enddo
   endif
endif
!
utin=0.
utout=0.
do k=1,nz
do j=1,ny-1
   utin =utin +(yp(j+1)-yp(j))*0.5*(bxx1(j,k)+bxx1(j+1,k))
   utout=utout+(yp(j+1)-yp(j))*0.5*(bxxn(j,k)+bxxn(j+1,k))
enddo
enddo
!
utin=utin/utout
!
do k=1,nz
do j=1,ny
   bxxn(j,k)=bxxn(j,k)*utin
enddo
enddo
!

if (nclx.eq.2) then
   do k=1,nz
   do j=1,ny
      ux(1 ,j,k)=bxx1(j,k)
      uy(1 ,j,k)=bxy1(j,k)+dpdyx1(j,k)
      uz(1 ,j,k)=bxz1(j,k)+dpdzx1(j,k)
      ux(nx,j,k)=bxxn(j,k)
      uy(nx,j,k)=bxyn(j,k)+dpdyxn(j,k)
      uz(nx,j,k)=bxzn(j,k)+dpdzxn(j,k)
   enddo
   enddo
endif
if (ncly.eq.2) then
   do k=1,nz
   do i=1,nx
      ux(i,1 ,k)=byx1(i,k)+dpdxy1(i,k)
      uy(i,1 ,k)=byy1(i,k)
      uz(i,1 ,k)=byz1(i,k)+dpdzy1(i,k)
      ux(i,ny,k)=byxn(i,k)+dpdxyn(i,k)
      uy(i,ny,k)=byyn(i,k)
      uz(i,ny,k)=byzn(i,k)+dpdzyn(i,k)
   enddo
   enddo
endif
if (nclz.eq.2) then
   do j=1,ny
   do i=1,nx
      ux(i,j,1 )=bzx1(i,j)+dpdxz1(i,j)
      uy(i,j,1 )=bzy1(i,j)+dpdyz1(i,j)
      uz(i,j,1 )=bzz1(i,j)
      ux(i,j,nz)=bzxn(i,j)+dpdxzn(i,j)
      uy(i,j,nz)=bzyn(i,j)+dpdyzn(i,j)
      uz(i,j,nz)=bzzn(i,j)
   enddo
   enddo
endif

if (iecoule.eq.7) then
   do k=1,nz
   do j=1,ny
      ux(1,j,k)=0.
      uy(1,j,k)=0.
      ux(nx,j,k)=0.
      uy(nx,j,k)=0.             
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ux(i,1,k)=0.
      uy(i,1,k)=0. 
      ux(i,ny,k)=0.!PPmodif 16.*(x*x*x*x-2.*x*x*x+x*x)
      uy(i,ny,k)=0.
   enddo
   enddo
endif
if (iecoule.eq.5) then
   do k=1,nz
   do i=1,nx
      ux(i,1,k)=0.+dpdxy1(i,k)
      uy(i,1,k)=0.
      ux(i,ny,k)=0.+dpdxyn(i,k)
      uy(i,ny,k)=0.
      uz(i,1,k)=0.+dpdzy1(i,k)
      uz(i,ny,k)=0.+dpdzyn(i,k)
   enddo
   enddo
endif

return
end subroutine pre_correc

!*****************************************************************
!
subroutine inversion5(a,ja,jb,sr,a1,b1,mx,ny,mz,e)
!
!*****************************************************************

implicit none

integer :: mx,ny,mz,i,j,k,m,mi,jc,nxy
real(8),dimension(mx,mz,ny/2,5) :: a
real(8),dimension(mx,mz,ny/2) :: e
integer,dimension(2) :: ja,jb
real(8),dimension(mx,mz) :: sr
real(8),dimension(mx,mz) :: a1,b1

nxy=mx*mz

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,ny/2-2
do i=1,2
   mi=m+i
   do k=1,mz
   do j=1,mx 
      if (a(j,k,m,3).ne.0.) sr(j,k)=a(j,k,mi,3-i)/a(j,k,m,3)
      e(j,k,mi)=e(j,k,mi)-sr(j,k)*e(j,k,m)
   enddo
   enddo
   do jc=ja(i),jb(i)
   do k=1,mz
   do j=1,mx
      a(j,k,mi,jc)=a(j,k,mi,jc)-sr(j,k)*a(j,k,m,jc+i)
   enddo
   enddo
   enddo
enddo
enddo
do k=1,mz
do j=1,mx
   if (abs(a(j,k,ny/2-1,3)).gt.1.e-20) then
      sr(j,k)=a(j,k,ny/2,2)/a(j,k,ny/2-1,3)
   else
      sr(j,k)=0.
   endif
   b1(j,k)=a(j,k,ny/2,3)-sr(j,k)*a(j,k,ny/2-1,4)
   if (abs(b1(j,k)).gt.1.e-20) then
      a1(j,k)=sr(j,k)/b1(j,k)
      e(j,k,ny/2)=e(j,k,ny/2)/b1(j,k)-a1(j,k)*e(j,k,ny/2-1)
   else
      a1(j,k)=0.
      e(j,k,ny/2)=0.
   endif
   if (abs(a(j,k,ny/2-1,3)).gt.1.e-20) then
      b1(j,k)=1./a(j,k,ny/2-1,3)
   else
      b1(j,k)=0.
   endif
   a1(j,k)=a(j,k,ny/2-1,4)*b1(j,k)
   e(j,k,ny/2-1)=e(j,k,ny/2-1)*b1(j,k)-a1(j,k)*e(j,k,ny/2)
enddo
enddo
do i=ny/2-2,1,-1
do k=1,mz
do j=1,mx
   if (abs(a(j,k,i,3)).gt.1.e-20) then
      sr(j,k)=1./a(j,k,i,3)
   else
      sr(j,k)=0.
   endif
   a1(j,k)=a(j,k,i,4)*sr(j,k)
   b1(j,k)=a(j,k,i,5)*sr(j,k)
   e(j,k,i)=e(j,k,i)*sr(j,k)-a1(j,k)*e(j,k,i+1)-b1(j,k)*e(j,k,i+2)
enddo
enddo
enddo

return
end subroutine inversion5

!*****************************************************************
!
subroutine inversion5_2d(a,ja,jb,sr,a1,b1,mx,ny,e)
!
!*****************************************************************

implicit none

integer :: mx,ny,i,j,m,mi,jc,nxy
real(8),dimension(mx,ny/2,5) :: a
real(8),dimension(mx,ny/2) :: e
integer,dimension(2) :: ja,jb
real(8),dimension(mx) :: sr,a1,b1

nxy=mx

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,ny/2-2
do i=1,2
   mi=m+i
   do j=1,mx 
      if (a(j,m,3).ne.0.) sr(j)=a(j,mi,3-i)/a(j,m,3)
      e(j,mi)=e(j,mi)-sr(j)*e(j,m)
   enddo
   do jc=ja(i),jb(i)
   do j=1,mx
      a(j,mi,jc)=a(j,mi,jc)-sr(j)*a(j,m,jc+i)
   enddo
   enddo
enddo
enddo
do j=1,mx
   if (abs(a(j,ny/2-1,3)).gt.1.e-20) then
      sr(j)=a(j,ny/2,2)/a(j,ny/2-1,3)
   else
      sr(j)=0.
   endif
   b1(j)=a(j,ny/2,3)-sr(j)*a(j,ny/2-1,4)
   if (abs(b1(j)).gt.1.e-20) then
      a1(j)=sr(j)/b1(j)
      e(j,ny/2)=e(j,ny/2)/b1(j)-a1(j)*e(j,ny/2-1)
   else
      a1(j)=0.
      e(j,ny/2)=0.
   endif
   if (abs(a(j,ny/2-1,3)).gt.1.e-20) then
      b1(j)=1./a(j,ny/2-1,3)
   else
      b1(j)=0.
   endif
   a1(j)=a(j,ny/2-1,4)*b1(j)
   e(j,ny/2-1)=e(j,ny/2-1)*b1(j)-a1(j)*e(j,ny/2)
enddo
do i=ny/2-2,1,-1
do j=1,mx
   if (abs(a(j,i,3)).gt.1.e-20) then
      sr(j)=1./a(j,i,3)
   else
      sr(j)=0.
   endif
   a1(j)=a(j,i,4)*sr(j)
   b1(j)=a(j,i,5)*sr(j)
   e(j,i)=e(j,i)*sr(j)-a1(j)*e(j,i+1)-b1(j)*e(j,i+2)
enddo
enddo

return
end subroutine inversion5_2d

!*****************************************************************
!
subroutine inversion5_i3(a,ja,jb,sr,a1,b1,mx,ny,mz,e)
!
!*****************************************************************

implicit none

integer :: mx,ny,mz,i,j,k,m,mi,jc,nxy
real(8),dimension(mx,mz,ny,5) :: a
real(8),dimension(mx,mz,ny) :: e
integer,dimension(2) :: ja,jb
real(8),dimension(mx,mz) :: sr
real(8),dimension(mx,mz) :: a1,b1

nxy=mx*mz

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,ny-2
do i=1,2
   mi=m+i
   do k=1,mz
   do j=1,mx 
      if (a(j,k,m,3).ne.0.) sr(j,k)=a(j,k,mi,3-i)/a(j,k,m,3)
      e(j,k,mi)=e(j,k,mi)-sr(j,k)*e(j,k,m)
   enddo
   enddo
   do jc=ja(i),jb(i)
   do k=1,mz
   do j=1,mx
      a(j,k,mi,jc)=a(j,k,mi,jc)-sr(j,k)*a(j,k,m,jc+i)
   enddo
   enddo
   enddo
enddo
enddo
do k=1,mz
do j=1,mx
   if (abs(a(j,k,ny-1,3)).gt.1.e-20) then
      sr(j,k)=a(j,k,ny,2)/a(j,k,ny-1,3)
   else
      sr(j,k)=0.
   endif
   b1(j,k)=a(j,k,ny,3)-sr(j,k)*a(j,k,ny-1,4)
   if (abs(b1(j,k)).gt.1.e-20) then
      a1(j,k)=sr(j,k)/b1(j,k)
      e(j,k,ny)=e(j,k,ny)/b1(j,k)-a1(j,k)*e(j,k,ny-1)
   else
      a1(j,k)=0.
      e(j,k,ny)=0.
   endif
   if (abs(a(j,k,ny-1,3)).gt.1.e-20) then
      b1(j,k)=1./a(j,k,ny-1,3)
   else
      b1(j,k)=0.
   endif
   a1(j,k)=a(j,k,ny-1,4)*b1(j,k)
   e(j,k,ny-1)=e(j,k,ny-1)*b1(j,k)-a1(j,k)*e(j,k,ny)
enddo
enddo
do i=ny-2,1,-1
do k=1,mz
do j=1,mx
   if (abs(a(j,k,i,3)).gt.1.e-20) then
      sr(j,k)=1./a(j,k,i,3)
   else
      sr(j,k)=0.
   endif
   a1(j,k)=a(j,k,i,4)*sr(j,k)
   b1(j,k)=a(j,k,i,5)*sr(j,k)
   e(j,k,i)=e(j,k,i)*sr(j,k)-a1(j,k)*e(j,k,i+1)-b1(j,k)*e(j,k,i+2)
enddo
enddo
enddo

return
end subroutine inversion5_i3

!*****************************************************************
!
subroutine inversion5_2d_i3(a,ja,jb,sr,a1,b1,mx,ny,e)
!
!*****************************************************************

implicit none

integer :: mx,ny,i,j,m,mi,jc,nxy
real(8),dimension(mx,ny,5) :: a
real(8),dimension(mx,ny) :: e
integer,dimension(2) :: ja,jb
real(8),dimension(mx) :: sr,a1,b1

nxy=mx
do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,ny-2
do i=1,2
   mi=m+i
   do j=1,mx 
      if (a(j,m,3).ne.0.) sr(j)=a(j,mi,3-i)/a(j,m,3)
      e(j,mi)=e(j,mi)-sr(j)*e(j,m)
   enddo
   do jc=ja(i),jb(i)
   do j=1,mx
         a(j,mi,jc)=a(j,mi,jc)-sr(j)*a(j,m,jc+i)
   enddo
   enddo
enddo
enddo
do j=1,mx
   if (abs(a(j,ny-1,3)).gt.1.e-9) then
      sr(j)=a(j,ny,2)/a(j,ny-1,3)
   else
      sr(j)=0.
   endif
   b1(j)=a(j,ny,3)-sr(j)*a(j,ny-1,4)
   if (abs(b1(j)).gt.1.e-9) then
      a1(j)=sr(j)/b1(j)
      e(j,ny)=e(j,ny)/b1(j)-a1(j)*e(j,ny-1)
   else
      a1(j)=0.
      e(j,ny)=0.
   endif
   if (abs(a(j,ny-1,3)).gt.1.e-9) then
      b1(j)=1./a(j,ny-1,3)
   else
      b1(j)=0.
   endif
   a1(j)=a(j,ny-1,4)*b1(j)
   e(j,ny-1)=e(j,ny-1)*b1(j)-a1(j)*e(j,ny)
enddo
do i=ny-2,1,-1
do j=1,mx
   if (abs(a(j,i,3)).gt.1.e-9) then
      sr(j)=1./a(j,i,3)
   else
      sr(j)=0.
   endif
   a1(j)=a(j,i,4)*sr(j)
   b1(j)=a(j,i,5)*sr(j)
   e(j,i)=e(j,i)*sr(j)-a1(j)*e(j,i+1)-b1(j)*e(j,i+2)
enddo
enddo

return
end subroutine inversion5_2d_i3

!**************************************************************************
!
subroutine matrice_raffinement(a,a2,a3,d1,d)
!
!**************************************************************************

USE param
USE derivX
USE derivY
USE derivZ
USE variables

implicit none

integer :: i,j,k
real(8),dimension(mx) :: transx
real(8),dimension(my) :: transy
real(8),dimension(mz) :: transz
real(8),dimension(mx,mz,ny/2,5) :: a,a2
real(8),dimension(mx,mz,ny,5) :: a3
real(8),dimension(mx,mz,ny) :: d1,d
real(8) :: ytt,xtt,ztt,yt,xt,zt,yt1,xt1,zt1,xa0,xa1 

do k=1,mz
   ztt=(biciz6*2.*cos(ezs(k)*3.*dz/2.)+ciciz6*2.*cos(ezs(k)*5.*dz/2.))
   zt=(aiciz6*2.*cos(ezs(k)*dz/2.))
   zt1=(1.+2.*ailcaiz6*cos(ezs(k)*dz))
   transz(k)=(ztt+zt)/zt1
enddo
do j=1,my
   ytt=(biciy6*2.*cos(eys(j)*3.*dy/2.)+ciciy6*2.*cos(eys(j)*5.*dy/2.))
   yt=(aiciy6*2.*cos(eys(j)*dy/2.))
   yt1=(1.+2.*ailcaiy6*cos(eys(j)*dy))
   transy(j)=(ytt+yt)/yt1
enddo
do i=1,mx
   xtt=(bicix6*2.*cos(exs(i)*3.*dx/2.)+cicix6*2.*cos(exs(i)*5.*dx/2.))
   xt=(aicix6*2.*cos(exs(i)*dx/2.))
   xt1=(1.+2.*ailcaix6*cos(exs(i)*dx))
   transx(i)=(xtt+xt)/xt1
enddo

if ((istret==1).or.(istret==2)) then
   xa0=alpha/pi+1./2./beta/pi
   if (istret==1) xa1=1./4./beta/pi
   if (istret==2) xa1=-1./4./beta/pi   
   
   do k=1,mz
   do j=1,ny/2
   do i=1,mx
      d(i,k,j)=yky(2*j-1)*transx(i)*transz(k)
      d1(i,k,j)=yky(2*j)*transx(i)*transz(k)
   enddo
   enddo
   enddo
!**************************DIAGONALE PRINCIPALE DE LA MATRICE************************************************
   do k=1,mz
   do j=2,ny/2-1
   do i=1,mx
      a(i,k,j,3)=-(xk2(i)*transy(2*j-1)*transy(2*j-1)*transz(k)*transz(k))&
           -(zk2(k)*transy(2*j-1)*transy(2*j-1)*transx(i)*transx(i))&
           -d(i,k,j)*d(i,k,j)*xa0*xa0-xa1*xa1*(d(i,k,j)*d(i,k,j-1)+d(i,k,j)*d(i,k,j+1))
      a2(i,k,j,3)=-(xk2(i)*transy(2*j)*transy(2*j)*transz(k)*transz(k))&
           -(zk2(k)*transy(2*j)*transy(2*j)*transx(i)*transx(i))&
           -d1(i,k,j)*d1(i,k,j)*xa0*xa0-xa1*xa1*(d1(i,k,j)*d1(i,k,j-1)+d1(i,k,j)*d1(i,k,j+1))
   enddo
   enddo
   enddo
   do k=1,mz
   do i=1,mx
      a(i,k,1,3)=-(xk2(i)*transy(1)*transy(1)*transz(k)*transz(k))&
           -(zk2(k)*transy(1)*transy(1)*transx(i)*transx(i))&
           -d(i,k,1)*d(i,k,1)*xa0*xa0-xa1*xa1*(d(i,k,1)*d(i,k,2))
      a(i,k,ny/2,3)=-(xk2(i)*transy(ny-2)*transy(ny-2)*transz(k)*transz(k))&
           -(zk2(k)*transy(ny-2)*transy(ny-2)*transx(i)*transx(i))&
           -d(i,k,ny/2)*d(i,k,ny/2)*xa0*xa0-xa1*xa1*(d(i,k,ny/2)*d(i,k,ny/2-1))
      a2(i,k,1,3)=-(xk2(i)*transy(2)*transy(2)*transz(k)*transz(k))&
           -(zk2(k)*transy(2)*transy(2)*transx(i)*transx(i))&
           -d1(i,k,1)*d1(i,k,1)*(xa0-xa1)*(xa0+xa1)-xa1*xa1*(d1(i,k,1)*d1(i,k,2))
      a2(i,k,ny/2,3)=-(xk2(i)*transy(ny-1)*transy(ny-1)*transz(k)*transz(k))&
           -(zk2(k)*transy(ny-1)*transy(ny-1)*transx(i)*transx(i))&
           -d1(i,k,ny/2)*d1(i,k,ny/2)*(xa0+xa1)*(xa0+xa1)-xa1*xa1*(d1(i,k,ny/2)*d1(i,k,ny/2-1))
   enddo
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP+1 DE LA MATRICE******************************************************
   do k=1,mz
   do i=1,mx
   do j=2,ny/2-1
      a(i,k,j,4)=xa0*xa1*(d(i,k,j)*d(i,k,j+1)+d(i,k,j+1)*d(i,k,j+1))
      a2(i,k,j,4)=xa0*xa1*(d1(i,k,j)*d1(i,k,j+1)+d1(i,k,j+1)*d1(i,k,j+1))
   enddo
   a(i,k,1,4)=(xa0*xa1*(d(i,k,1)*d(i,k,2)+d(i,k,2)*d(i,k,2)))
   a2(i,k,1,4)=(xa0-xa1)*xa1*(d1(i,k,1)*d1(i,k,2))&
        +xa0*xa1*(d1(i,k,2)*d1(i,k,2))
   a2(i,k,ny/2-1,4)=xa0*xa1*(d1(i,k,ny/2-1)*d1(i,k,ny/2))&
        +(xa0+xa1)*xa1*(d1(i,k,ny/2)*d1(i,k,ny/2))
   a2(i,k,ny/2,4)=0.
   enddo
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP+2 DE LA MATRICE******************************************************
   do k=1,mz
   do i=1,mx
   do j=1,ny/2-2
      a(i,k,j,5)=-d(i,k,j+1)*d(i,k,j+2)*xa1*xa1
      a2(i,k,j,5)=-d1(i,k,j+1)*d1(i,k,j+2)*xa1*xa1
   enddo
   a(i,k,1,5)=a(i,k,1,5)*2.
   a(i,k,ny/2-1,5)=0.
   a(i,k,ny/2,5)=0.
   a2(i,k,ny/2-1,5)=0.
   a2(i,k,ny/2,5)=0. 
   enddo
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP-1 DE LA MATRICE******************************************************
   do k=1,mz
   do i=1,mx
   do j=2,ny/2
      a(i,k,j,2)=xa0*xa1*(d(i,k,j)*d(i,k,j-1)+d(i,k,j-1)*d(i,k,j-1))
      a2(i,k,j,2)=xa0*xa1*(d1(i,k,j)*d1(i,k,j-1)+d1(i,k,j-1)*d1(i,k,j-1))
   enddo
   a(i,k,1,2)=0.
   a2(i,k,1,2)=0.
   a2(i,k,2,2)=xa0*xa1*(d1(i,k,2)*d1(i,k,1))&
        +(xa0+xa1)*xa1*(d1(i,k,1)*d1(i,k,1))
   a2(i,k,ny/2,2)=(xa0+xa1)*xa1*(d1(i,k,ny/2)*d1(i,k,ny/2-1))&
        +xa0*xa1*(d1(i,k,ny/2-1)*d1(i,k,ny/2-1))
   enddo
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP-2 DE LA MATRICE******************************************************
   do k=1,mz
   do i=1,mx
   do j=3,ny/2
      a(i,k,j,1)=-d(i,k,j-1)*d(i,k,j-2)*xa1*xa1
      a2(i,k,j,1)=-d1(i,k,j-1)*d1(i,k,j-2)*xa1*xa1
   enddo
   a(i,k,1,1)=0.
   a(i,k,2,1)=0.
   a2(i,k,1,1)=0.
   a2(i,k,2,1)=0.
   enddo
   enddo
!*************************************************************************************************************     
   do k=1,mz
   do i=1,mx
      if ((xkx(i)==0.).and.(zkz(k)==0)) then
         a(i,k,1,3)=1.
         a(i,k,1,4)=0.
         a(i,k,1,5)=0.
      endif
   enddo
   enddo
!*************************************************************************************************************
endif

if (istret==3) then
   xa0=alpha/pi+1./2./beta/pi
   xa1=-1./4./beta/pi 
   do k=1,mz
   do j=1,ny
   do i=1,mx
      d(i,k,j)=yky(j)*transx(i)*transz(k)
   enddo
   enddo
   enddo
!**************************DIAGONALE PRINCIPALE DE LA MATRICE************************************************
   do k=1,mz
   do j=2,ny-1
   do i=1,mx
      a3(i,k,j,3)=-(xk2(i)*transy(j)*transy(j)*transz(k)*transz(k))&
           -(zk2(k)*transy(j)*transy(j)*transx(i)*transx(i))&
           -d(i,k,j)*d(i,k,j)*xa0*xa0-xa1*xa1*(d(i,k,j)*d(i,k,j-1)+d(i,k,j)*d(i,k,j+1))
   enddo
   enddo
   enddo
   do k=1,mz
   do i=1,mx
      a3(i,k,1,3)=-(xk2(i)*transy(1)*transy(1)*transz(k)*transz(k))&
           -(zk2(k)*transy(1)*transy(1)*transx(i)*transx(i))&
           -d(i,k,1)*d(i,k,1)*xa0*xa0-xa1*xa1*(d(i,k,1)*d(i,k,2))
      a3(i,k,ny,3)=-(xk2(i)*transy(ny)*transy(ny)*transz(k)*transz(k))&
           -(zk2(k)*transy(ny)*transy(ny)*transx(i)*transx(i))&
           -d(i,k,ny)*d(i,k,ny)*xa0*xa0-xa1*xa1*(d(i,k,ny)*d(i,k,ny-1))
   enddo
   enddo    
!*************************************************************************************************************     
!**************************DIAGONALE SUP+1 DE LA MATRICE******************************************************
   do k=1,mz
   do i=1,mx
   do j=2,ny-1
      a3(i,k,j,4)=xa0*xa1*(d(i,k,j)*d(i,k,j+1)+d(i,k,j+1)*d(i,k,j+1))
   enddo
   a3(i,k,1,4)=(xa0*xa1*(d(i,k,1)*d(i,k,2)+d(i,k,2)*d(i,k,2)))
   enddo
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP+2 DE LA MATRICE******************************************************
   do k=1,mz
   do i=1,mx
   do j=1,ny-2
      a3(i,k,j,5)=-d(i,k,j+1)*d(i,k,j+2)*xa1*xa1
   enddo
   a3(i,k,1,5)=a(i,k,1,5)*2.
   a3(i,k,ny-1,5)=0.
   a3(i,k,ny,5)=0.
   enddo
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP-1 DE LA MATRICE******************************************************
   do k=1,mz
   do i=1,mx
   do j=2,ny
      a3(i,k,j,2)=xa0*xa1*(d(i,k,j)*d(i,k,j-1)+d(i,k,j-1)*d(i,k,j-1))
   enddo
   a3(i,k,1,2)=0.
   enddo
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP-2 DE LA MATRICE******************************************************
   do k=1,mz
   do i=1,mx
   do j=3,ny
      a3(i,k,j,1)=-d(i,k,j-1)*d(i,k,j-2)*xa1*xa1
   enddo
   a3(i,k,1,1)=0.
   a3(i,k,2,1)=0.
   enddo
   enddo  
!*************************************************************************************************************     
   do k=1,mz
   do i=1,mx
      if ((xkx(i)==0.).and.(zkz(k)==0)) then
         a3(i,k,1,3)=1.
         a3(i,k,1,4)=0.
         a3(i,k,1,5)=0.
      endif
   enddo
   enddo
!*************************************************************************************************************
endif

return
end subroutine matrice_raffinement

!**************************************************************************
!
subroutine matrice_raffinement_2d(a,a2,a3,d1,d)
!
!**************************************************************************

USE param
USE derivX
USE derivY
USE derivZ
USE variables

implicit none

integer :: i,j
real(8),dimension(mx,ny/2,5) :: a,a2
real(8),dimension(mx,ny,5) :: a3
real(8),dimension(mx,ny) :: d1,d
real(8),dimension(my) :: transy
real(8),dimension(mx) :: transx
real(8) :: ytt,xtt,yt,xt,yt1,xt1,xa0,xa1 
      
do j=1,my
   ytt=(biciy6*2.*cos(eys(j)*3.*dy/2.)+ciciy6*2.*cos(eys(j)*5.*dy/2.))
   yt=(aiciy6*2.*cos(eys(j)*dy/2.))
   yt1=(1.+2.*ailcaiy6*cos(eys(j)*dy))
   transy(j)=(ytt+yt)/yt1
enddo
do i=1,mx
   xtt=(bicix6*2.*cos(exs(i)*3.*dx/2.)+cicix6*2.*cos(exs(i)*5.*dx/2.))
   xt=(aicix6*2.*cos(exs(i)*dx/2.))
   xt1=(1.+2.*ailcaix6*cos(exs(i)*dx))
   transx(i)=(xtt+xt)/xt1
enddo

if ((istret.eq.1).or.(istret.eq.2)) then
   xa0=alpha/pi+1./2./beta/pi
   if (istret==1) xa1=1./4./beta/pi
   if (istret==2) xa1=-1./4./beta/pi
   do j=1,ny/2
   do i=1,mx
      d(i,j)=yky(2*j-1)*transx(i)
      d1(i,j)=yky(2*j)*transx(i)
   enddo
   enddo    
!**************************DIAGONALE PRINCIPALE DE LA MATRICE************************************************
   do j=2,ny/2-1
   do i=1,mx
      a(i,j,3)=-(xk2(i)*transy(2*j-1)*transy(2*j-1))&
           -d(i,j)*d(i,j)*xa0*xa0-xa1*xa1*(d(i,j)*d(i,j-1)+d(i,j)*d(i,j+1))
      a2(i,j,3)=-(xk2(i)*transy(2*j)*transy(2*j))&
           -d1(i,j)*d1(i,j)*xa0*xa0-xa1*xa1*(d1(i,j)*d1(i,j-1)+d1(i,j)*d1(i,j+1))
   enddo
   enddo
   do i=1,mx
      a(i,1,3)=-(xk2(i)*transy(1)*transy(1))&
           -d(i,1)*d(i,1)*xa0*xa0-xa1*xa1*(d(i,1)*d(i,2))
      a(i,ny/2,3)=-(xk2(i)*transy(ny-2)*transy(ny-2))&
           -d(i,ny/2)*d(i,ny/2)*xa0*xa0-xa1*xa1*(d(i,ny/2)*d(i,ny/2-1))
      a2(i,1,3)=-(xk2(i)*transy(2)*transy(2))&
           -d1(i,1)*d1(i,1)*(xa0-xa1)*(xa0+xa1)-xa1*xa1*(d1(i,1)*d1(i,2))
      a2(i,ny/2,3)=-(xk2(i)*transy(ny-1)*transy(ny-1))&
           -d1(i,ny/2)*d1(i,ny/2)*(xa0+xa1)*(xa0+xa1)-xa1*xa1*(d1(i,ny/2)*d1(i,ny/2-1))
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP+1 DE LA MATRICE******************************************************
   do i=1,mx
   do j=2,ny/2-1
      a(i,j,4)=xa0*xa1*(d(i,j)*d(i,j+1)+d(i,j+1)*d(i,j+1))
      a2(i,j,4)=xa0*xa1*(d1(i,j)*d1(i,j+1)+d1(i,j+1)*d1(i,j+1))
   enddo
   a(i,1,4)=(xa0*xa1*(d(i,1)*d(i,2)+d(i,2)*d(i,2)))
   a2(i,1,4)=(xa0-xa1)*xa1*(d1(i,1)*d1(i,2))&
        +xa0*xa1*(d1(i,2)*d1(i,2))
   a2(i,ny/2-1,4)=xa0*xa1*(d1(i,ny/2-1)*d1(i,ny/2))&
        +(xa0+xa1)*xa1*(d1(i,ny/2)*d1(i,ny/2))
   a2(i,ny/2,4)=0.
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP+2 DE LA MATRICE******************************************************
   do i=1,mx
   do j=1,ny/2-2
      a(i,j,5)=-d(i,j+1)*d(i,j+2)*xa1*xa1
      a2(i,j,5)=-d1(i,j+1)*d1(i,j+2)*xa1*xa1
   enddo
   a(i,1,5)=a(i,1,5)*2.
   a(i,ny/2-1,5)=0.
   a(i,ny/2,5)=0.
   a2(i,ny/2-1,5)=0.
   a2(i,ny/2,5)=0. 
9999   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP-1 DE LA MATRICE******************************************************
   do i=1,mx
   do j=2,ny/2
      a(i,j,2)=xa0*xa1*(d(i,j)*d(i,j-1)+d(i,j-1)*d(i,j-1))
      a2(i,j,2)=xa0*xa1*(d1(i,j)*d1(i,j-1)+d1(i,j-1)*d1(i,j-1))
   enddo
   a(i,1,2)=0.
   a2(i,1,2)=0.
   a2(i,2,2)=xa0*xa1*(d1(i,2)*d1(i,1))&
        +(xa0+xa1)*xa1*(d1(i,1)*d1(i,1))
   a2(i,ny/2,2)=(xa0+xa1)*xa1*(d1(i,ny/2)*d1(i,ny/2-1))&
        +xa0*xa1*(d1(i,ny/2-1)*d1(i,ny/2-1))
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP-2 DE LA MATRICE******************************************************
   do i=1,mx
   do j=3,ny/2
      a(i,j,1)=-d(i,j-1)*d(i,j-2)*xa1*xa1
      a2(i,j,1)=-d1(i,j-1)*d1(i,j-2)*xa1*xa1
   enddo
   a(i,1,1)=0.
   a(i,2,1)=0.
   a2(i,1,1)=0.
   a2(i,2,1)=0.
   enddo
!*************************************************************************************************************  
   do i=1,mx
      if ((xkx(i)==0.)) then
         a(i,1,3)=1.
         a(i,1,4)=0.
         a(i,1,5)=0.
      endif
   enddo
!*************************************************************************************************************
endif
if (istret.eq.3) then
   xa0=alpha/pi+1./2./beta/pi
   xa1=-1./4./beta/pi 
   do j=1,ny
   do i=1,mx
      d(i,j)=yky(j)*transx(i)
   enddo
   enddo
!**************************DIAGONALE PRINCIPALE DE LA MATRICE************************************************
   do j=2,ny-1
   do i=1,mx
      a3(i,j,3)=-(xk2(i)*transy(j)*transy(j))&
           -d(i,j)*d(i,j)*xa0*xa0-xa1*xa1*(d(i,j)*d(i,j-1)+d(i,j)*d(i,j+1))
   enddo
   enddo
   do i=1,mx
      a3(i,1,3)=-(xk2(i)*transy(1)*transy(1))&
           -d(i,1)*d(i,1)*xa0*xa0-xa1*xa1*(d(i,1)*d(i,2))
      a3(i,ny,3)=-(xk2(i)*transy(ny)*transy(ny))&
           -d(i,ny)*d(i,ny)*xa0*xa0-xa1*xa1*(d(i,ny)*d(i,ny-1))
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP+1 DE LA MATRICE******************************************************
   do i=1,mx
   do j=2,ny-1
      a3(i,j,4)=xa0*xa1*(d(i,j)*d(i,j+1)+d(i,j+1)*d(i,j+1))
   enddo
   a3(i,1,4)=(xa0*xa1*(d(i,1)*d(i,2)+d(i,2)*d(i,2)))
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP+2 DE LA MATRICE******************************************************
   do i=1,mx
   do j=1,ny-2
      a3(i,j,5)=-d(i,j+1)*d(i,j+2)*xa1*xa1
   enddo
   a3(i,1,5)=a(i,1,5)*2.
   a3(i,ny-1,5)=0.
   a3(i,ny,5)=0.
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP-1 DE LA MATRICE******************************************************
   do i=1,mx
   do j=2,ny
      a3(i,j,2)=xa0*xa1*(d(i,j)*d(i,j-1)+d(i,j-1)*d(i,j-1))
   enddo
   a3(i,1,2)=0.
   enddo
!*************************************************************************************************************     
!**************************DIAGONALE SUP-2 DE LA MATRICE******************************************************
   do i=1,mx
   do j=3,ny
      a3(i,j,1)=-d(i,j-1)*d(i,j-2)*xa1*xa1
   enddo
   a3(i,1,1)=0.
   a3(i,2,1)=0.
   enddo
!*************************************************************************************************************     
   do i=1,mx
      if ((xkx(i)==0.)) then
         a3(i,1,3)=1.
         a3(i,1,4)=0.
         a3(i,1,5)=0.
      endif
   enddo
!******************************************************************
endif

return
end subroutine matrice_raffinement_2d

!******************************************************************
!
subroutine numcar2(num,car) 
!
!******************************************************************

implicit none

integer :: num 
character :: car*4 

if (num >= 1000) then 
   write (car,4) num 
4  format(i4) 
else 
   if (num >= 100) then 
      write (car,1) num 
1     format('0',i3) 
   else 
      if (num >= 10) then 
         write (car,2) num 
2        format('00',i2) 
      else 
         write (car,3) num 
3        format('000',i1) 
      endif
   endif
endif

return  
end subroutine numcar2

!*******************************************************************
!
subroutine minmax(u,car)
!
!*******************************************************************
USE param
USE variables

implicit none

integer :: i,j,k,n,nimax,njmax
real(8) :: umax,umin
real(8),dimension(nx,ny,nz) :: u
character*45 :: printminmax
character :: car*2

umax=-10000.
umin= 10000. 

do k=1,nz
do j=1,ny
do i=1,nx
   if (u(i,j,k).gt.umax) umax=u(i,j,k)
   if (u(i,j,k).lt.umin) umin=u(i,j,k)      
enddo
enddo
enddo

 901 format(A,'(min,max)=(',F10.5,';',F10.5,')')
 write(printminmax,901) car,umin,umax
 write(*,*) printminmax
return
end subroutine minmax
