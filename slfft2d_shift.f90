module slfft2d_shift_m
implicit none
contains
!*****************************************************************************
!
subroutine slfft2d_shift(ppm,tb6,tb1,tb2,tb3,tb4,tb5,work,table,nx,nxm,ny,nym,&
  nz,mx,my,mz,nwork,ntrigsX,ntrigsY,isign)
!
!*****************************************************************************

USE param 

implicit none

integer :: mxyz,mxy,i,j,nx,mx,nxm,ny,my,nym,nz,mz,mxxyyz,nwork,iplan
real(8),dimension(nxm,nym,nz) :: ppm
real(8),dimension(nxm,nym) :: tbf
real(8),dimension(nwork) :: work
real(8), dimension(100+2*(nxm+nym+1)) :: table
real(8),dimension(mx,my) :: tb1,tb2,tb3,tb4,tb5,tb6,tb7,tb8 
integer :: ntrigsX,ntrigsY,isign
real(8) :: xx3,xx4
real(8) :: xa,xb,xx1,xx2,xa1,xb1
real(8) :: xx5,xx6,xx7,xx8,yy1,yy2,yy3,yy4,yy5,yy6,yy7,yy8

mxyz=mx*my*mz
mxy=mx*my
mxxyyz=nxm*nym*nz

tb1(:,:)=0. ;  tb2(:,:)=0. ;  tb3(:,:)=0. ;  tb4(:,:)=0. 
if (((nclx.eq.2).or.(nclx.eq.1)).and.((ncly.eq.1).or.(ncly.eq.2))) then
   if (isign==1) then
      do j=1,nym
      do i=1,nxm/2
         tb1(i,j)=ppm(2*(i-1)+1,j,1)
      enddo
      do i=nxm/2+1,nxm
         tb1(i,j)=ppm(2*nxm-2*i+2,j,1)
      enddo
      enddo
      do i=1,nxm
      do j=1,nym/2
         tb2(i,j)=tb1(i,2*(j-1)+1)
      enddo
      do j=nym/2+1,nym
         tb2(i,j)=tb1(i,2*nym-2*j+2)
      enddo
      enddo
      if (ifft==1) then
         call slfft2d(tb2,nxm,nym,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,-1)
      endif
      if (ifft==2) then
         do j=1,nym
         do i=1,nxm
            tbf(i,j)=tb2(i,j)
         enddo
         enddo
         tb2(:,:)=0.
         tb2=tb2/nxm/nym
      endif
      do j=1,nym
      do i=1,nxm/2
         tb3(i,j)=tb2(2*i-1,j)
      enddo
      tb3(nxm/2+1,j)=tb2(nxm+1,j)
      do i=2,nxm/2
         tb3(i+nxm/2,j)=tb3(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,nym
      do i=1,nxm/2
         tb4(i,j)=tb2(2*i,j)
      enddo
      tb4(nxm/2+1,j)=tb2(nx+1,j)
      do i=2,nxm/2
         tb4(i+nxm/2,j)=-tb4(nxm/2-i+2,j)
      enddo
      enddo
      tb1(:,:)=0.; tb2(:,:)=0.
      do j=1,ny
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xx1=tb3(i,j)*xa/2.+tb4(i,j)*xb/2.
         xx2=tb3(i,nym-j+2)*xa/2.-tb4(i,nym-j+2)*xb/2.           
         tb1(i,j)=xx1+xx2      
      enddo
      do i=2,nxm/2
         tb1(i+nxm/2,j)=tb1(nxm/2-i+2,j)
      enddo
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xa1=cos((nym-j+1)/2.*pi/(nym))
         xb1=sin((nym-j+1)/2.*pi/(nym))
         xx1=-tb3(i,j)*xb/2.+tb4(i,j)*xa/2.
         xx2=tb3(i,nym-j+2)*xb/2.+tb4(i,nym-j+2)*xa/2.           
         tb2(i,j)=xx1+xx2 
      enddo
      do i=2,nxm/2
         tb2(i+nxm/2,j)=-tb2(nxm/2-i+2,j)
      enddo
      enddo
      do i=1,nxm
         tb1(i,1)=tb1(i,1)*2.
         tb2(i,1)=tb2(i,1)*2.
      enddo
      tb3(:,:)=0. ; tb4(:,:)=0.
      do j=1,ny
      do i=1,nxm
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         xx1=(tb1(i,j)+tb1(nxm-i+2,j))*xa/2
         xx2=(tb2(i,j)-tb2(nxm-i+2,j))*xb/2.
         tb3(i,j)=xx1+xx2
      enddo
      enddo
      do j=1,ny
         tb3(1,j)=2.*tb3(1,j)
      enddo
      tb3(mx,:)=0.
      tb3(mx-1,:)=0.
      tb3(:,my)=0.
      tb3(:,my-1)=0.
      tb6(:,:)=0. ; tb6(:,:)=tb3(:,:)
   endif
   if (isign==-1) then
      do j=1,ny
      do i=1,nxm/2+1
         tb1(i,j)=tb6(i,j)
         tb3(i,j)=tb6(nxm-i+2,j)
      enddo
      enddo
      tb5(:,:)=0. ; tb6(:,:)=0.; tb4(:,:)=0. ; tb2(:,:)=0.
      do j=1,ny
      do i=1,nxm/2+1
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         xb1=sin((nxm-i+1)/2.*pi/(nxm))
         xa1=cos((nxm-i+1)/2.*pi/(nxm))
         tb5(i,j)=tb1(i,j)*xa+tb3(i,j)*xa1
         tb6(i,j)=tb1(i,j)*xb-tb3(i,j)*xb1
      enddo
      do i=2,nxm/2
         tb5(i+nxm/2,j)=tb5(nxm/2-i+2,j)
         tb6(i+nxm/2,j)=-tb6(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,ny
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xa1=cos((nym-j+1)/2.*pi/(nym))
         xb1=sin((nym-j+1)/2.*pi/(nym))
         xx1=tb5(i,j)*xa-tb6(i,j)*xa1
         xx2=tb5(i,nym-j+2)*xb+tb6(i,nym-j+2)*xa 
         xx3=-tb5(i,nym-j+2)*xb1+tb6(i,nym-j+2)*xb
         xx4=tb5(i,j)*xb+tb6(i,j)*xa
         tb2(2*i-1,j)=xx1+xx2
         tb2(2*i,j)=xx3+xx4
      enddo
      enddo
      do i=1,mx
         tb2(i,ny)=0.
      enddo
      do i=1,mx
         tb2(i,ny)=tb2(i,(ny+1)/2)
      enddo
      if (ifft==1) then
         call slfft2d(tb2,nxm,nym,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,1)
      endif
      if (ifft==2) then
         do j=1,nym
         do i=1,nxm
            tb2(i,j)=tbf(i,j)
         enddo
         enddo
      endif
      tb1(:,:)=0. ; tb3(:,:)=0.
      do j=1,ny
      do i=1,nxm/2
         tb1(2*i-1,j)=tb2(i,j)
         tb1(2*i,j)=tb2(nxm-i+1,j)
      enddo
      enddo
      do j=1,nym/2   
      do i=1,nxm
         ppm(i,2*j-1,1)=tb1(i,j)
         ppm(i,2*j,1)=tb1(i,nym-j+1)
      enddo
      enddo
   endif
endif
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
if (((ncly.eq.2).or.(ncly.eq.1)).and.(nclx.eq.0)) then
   if (isign==1) then
      do j=1,nym/2
      do i=1,nxm
         tb1(i,j)=ppm(i,2*(j-1)+1,1)
      enddo
      enddo
      do j=nym/2+1,nym
      do i=1,nxm
         tb1(i,j)=ppm(i,2*nym-2*j+2,1)
      enddo
      enddo
     if (ifft==1) then
         call slfft2d(tb1,nxm,nym,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,-1)
      endif
      if (ifft==2) then
         do j=1,nym
         do i=1,nxm
            tbf(i,j)=tb1(i,j)
         enddo
         enddo
         tb1(:,:)=0.
         tb1=tb1/nxm/nym
      endif
      do j=1,nym
      do i=1,nxm/2
         tb2(i,j)=tb1(2*i-1,j)
      enddo
      tb2(nxm/2+1,j)=tb1(nxm+1,j)
      do i=2,nxm/2
         tb2(i+nxm/2,j)=tb2(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,nym
      do i=1,nxm/2
         tb3(i,j)=tb1(2*i,j)
      enddo
      tb3(nxm/2+1,j)=tb1(mx,j)
      do i=2,nxm/2
         tb3(i+nxm/2,j)=-tb3(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,my
         tb2(mx,j)=tb2(nxm/2+1,j)
         tb3(mx,j)=-tb3(nxm/2+1,j)
      enddo
      tb1(:,:)=0.
      do j=1,ny
      do i=1,nxm/2
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xx1=tb2(i,j)*xa/2.+tb3(i,j)*xb/2.
         xx2=tb2(i,nym-j+2)*xa/2.-tb3(i,nym-j+2)*xb/2.           
         tb4(i,j)=xx1+xx2      
      enddo
      do i=2,nxm/2+2
         tb4(i+nxm/2,j)=tb4(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,ny
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xx1=tb2(nxm/2+1,j)*xa/2.-tb3(nxm/2+1,nym-j+2)*xb/2.
         xx2=tb2(nxm/2+1,nym-j+2)*xa/2.+tb3(mx,nym-j+2)*xb/2.           
         tb4(nxm/2+1,j)=xx1+xx2
      enddo
      do j=1,ny
      do i=2,nxm/2+2
         tb4(i+nxm/2,j)=tb4(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,ny
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xa1=cos((nym-j+1)/2.*pi/(nym))
         xb1=sin((nym-j+1)/2.*pi/(nym))
         xx1=-tb2(i,j)*xb/2.+tb3(i,j)*xa/2.
         xx2=tb2(i,nym-j+2)*xb/2.+tb3(i,nym-j+2)*xa/2.           
         tb1(i,j)=xx1+xx2 
      enddo
      do i=2,nxm/2+2
         tb1(i+nxm/2,j)=-tb1(nxm/2-i+2,j)
      enddo
      enddo
      do i=1,mx
         tb4(i,1)=tb4(i,1)*2.
         tb1(i,1)=tb1(i,1)*2.
      enddo
      do j=1,ny
      do i=1,nxm/2+2
         xa=cos((i-1)*pi/(nxm))
         xb=sin((i-1)*pi/(nxm))
         xx1=tb4(i,j)*xa/2.+tb4(nxm-i+2,j)*xa/2.
         xx2=-tb1(nxm-i+2,j)*xb/2.+tb1(i,j)*xb/2. 
         xx3=tb4(i,j)*xb/2.+tb4(nxm-i+2,j)*xb/2.
         xx4=tb1(nxm-i+2,j)*xa/2.+tb1(nxm-i+2,j)*xa/2. 
         tb6(2*i-1,j)=xx1+xx2
         tb6(2*i,j)=-(xx3+xx4)
      enddo
      enddo
   endif
   if (isign==-1) then
      tb6(:,ny)=0.
      do j=1,my
         tb6(nxm+1,j)=-tb6(mx,j)
      enddo
      do j=2,my
         tb6(nxm/2+1,j)=tb6(nxm/2+1,j)/2.
         tb6(nxm/2+2,j)=tb6(nxm/2+2,j)/2.
      enddo
      do j=1,nym
      do i=1,nxm/2+1,2
         tb1(i,j)=tb6(i,j)
         tb1(i+1,j)=-tb6(i+1,j)
         tb2(i,j)=tb6(nxm-i+2,j)
         tb2(i+1,j)=-tb6(nxm-i+1,j)
      enddo
      enddo
      do j=1,my
         tb2(nxm/2+2,j)=0.
      enddo
      tb5(:,:)=0. ; tb6(:,:)=0.
      do j=1,ny+1
      do i=1,mx-1,2
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         tb3((i+1)/2,j)=tb1(i+1,j)*xb+tb2(nxm-i+1,j)*xb
         tb4((i+1)/2,j)=-tb1(i+1,j)*xa-tb2(nxm-i+1,j)*xa
         tb5((i+1)/2,j)=tb1(i,j)*xa+tb2(nxm-i+2,j)*xa
         tb6((i+1)/2,j)=tb1(i,j)*xb+tb2(nxm-i+2,j)*xb
      enddo
      enddo
      do j=1,ny
      do i=2,nxm/2
         tb3(i+nxm/2,j)=tb3(nxm/2-i+2,j)
         tb5(i+nxm/2,j)=tb5(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,ny
      do i=2,nxm/2
         tb4(i+nxm/2,j)=-tb4(nxm/2-i+2,j)
         tb6(i+nxm/2,j)=-tb6(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,my
         tb5(nxm/2+1,j)=tb6(nxm/2+1,j)
         tb6(nxm/2+1,j)=0.
      enddo
      tb1(:,:)=0.
      do j=1,ny+1
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xa1=cos((nym-j+1)/2.*pi/(nym))
         xb1=sin((nym-j+1)/2.*pi/(nym))
         xx1=(tb5(i,j)+tb3(i,j))*xa-(tb6(i,j)+tb4(i,j))*xa1
         xx2=(tb5(i,nym-j+2)+tb3(i,nym-j+2))*xb+(tb6(i,nym-j+2)+tb4(i,nym-j+2))*xa 
         xx3=-(tb5(i,nym-j+2)+tb3(i,nym-j+2))*xb1+(tb6(i,nym-j+2)+tb4(i,nym-j+2))*xb
         xx4=(tb5(i,j)+tb3(i,j))*xb+(tb6(i,j)+tb4(i,j))*xa
         tb1(2*i-1,j)=xx1+xx2
         tb1(2*i,j)=xx4+xx3
      enddo
      enddo
      do i=1,mx
         tb1(i,ny)=0.
      enddo
      tb1(nxm/2+1,1)=tb1(nxm/2+1,1)/2.
      tb1(nxm/2+2,1)=tb1(nxm/2+2,1)/2.
      tb1(mx,1)=0.
      if (ifft==1) then
         call slfft2d(tb1,nxm,nym,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,1)
      endif
      if (ifft==2) then
         do j=1,nym
         do i=1,nxm
            tb1(i,j)=tbf(i,j)
         enddo
         enddo
      endif
      do j=1,nym/2
      do i=1,nxm
         ppm(i,2*j-1,1)=tb1(i,j)
      enddo
      enddo
      do j=1,nym/2
      do i=1,nxm
         ppm(i,2*j,1)=tb1(i,nym-j+1) 
      enddo
      enddo
   endif
endif
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
if (((nclx.eq.2).or.(nclx.eq.1)).and.(ncly.eq.0)) then
   if (isign==1) then
      do j=1,nym
      do i=1,nxm/2
         tb1(i,j)=ppm(2*(i-1)+1,j,1)
      enddo
      do i=nxm/2+1,nxm
         tb1(i,j)=ppm(2*nxm-2*i+2,j,1)
      enddo
      enddo
      if (ifft==1) then
         call slfft2d(tb1,nxm,nym,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,-1)
      endif
      if (ifft==2) then
         do j=1,nym
         do i=1,nxm
            tbf(i,j)=tb1(i,j)
         enddo
         enddo
         tb1(:,:)=0.
         tb1=tb1/nxm/nym
      endif
      do j=1,nym
      do i=1,nxm/2
         tb2(i,j)=tb1(2*i-1,j)
         tb3(i,j)=tb1(2*i,j)
      enddo
      tb2(nxm/2+1,j)=tb1(nxm+1,j)
      tb3(nxm/2+1,j)=tb1(nx+1,j)
      do i=2,nxm/2
         tb2(i+nxm/2,j)=tb2(nxm/2-i+2,j)
         tb3(i+nxm/2,j)=-tb3(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,nym
      do i=1,nx+1
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         tb4(i,j)=(tb2(i,j)+tb2(nxm-i+2,j))*&
             xa/2.+(tb3(i,j)-tb3(nxm-i+2,j))*xb/2.
      enddo
      enddo
      do j=1,nym
         xb=sin((nxm/2)/2.*pi/(nxm))
         xa=cos((nxm/2)/2.*pi/(nxm))
         tb4(nx+1,j)=-((tb3(nxm/2+1,j)*xa/2.)+(tb3(nxm/2+1,j)*xb/2.))
      enddo
      do j=1,ny
         tb4(1,j)=2.*tb4(1,j)
         tb4(nx,j)=2.*tb4(nx,j)
      enddo
      do i=1,nx+1
         tb4(i,1)=2.*tb4(i,1)
      enddo
      tb1(:,:)=0.
      do i=1,mx
         tb1(i,nym/2+1)=tb4(i,nym/2+1)
      enddo
      tb6(:,:)=0.
      do j=1,nym/2
      do i=1,nxm/2
         xa=cos((j-1)*pi/(nym))
         xb=sin((j-1)*pi/(nym))
         xx1=tb4(i,j)*xa/2.+tb4(i,nym-j+2)*xa/2.
         xx2=-tb4(nxm-i+2,j)*xb/2.+tb4(nxm-i+2,nym-j+2)*xb/2. 
         xx3=tb4(i,j)*xb/2.+tb4(i,nym-j+2)*xb/2.
         xx4=tb4(nxm-i+2,j)*xa/2.-tb4(nxm-i+2,nym-j+2)*xa/2. 
         tb6(i,2*j-1)=xx1+xx2
         tb6(i,2*j)=-(xx3+xx4)
      enddo
      do i=nxm/2+1,nx+1
         xa=cos((nym-j+1)*pi/(nym))
         xb=sin((nym-j+1)*pi/(nym))
         xx1=tb4(i,j)*xa/2.+tb4(i,nym-j+2)*xa/2.
         xx2=-tb4(nxm-i+2,j)*xb/2.+tb4(nxm-i+2,nym-j+2)*xb/2.
         xx3=tb4(i,j)*xb/2.+tb4(i,nym-j+2)*xb/2.
         xx4=tb4(nxm-i+2,j)*xa/2.-tb4(nxm-i+2,nym-j+2)*xa/2.  
         tb6(i,2*j-1)=-(xx1+xx2)
         tb6(i,2*j)=-(xx3+xx4) 
      enddo
      xa=cos((nym-j+1)*pi/(nym))
      xb=sin((nym-j+1)*pi/(nym))
      xx1=tb4(nxm/2+1,j)*xa/2.+tb4(nxm/2+1,nym-j+2)*xa/2.
      xx2=-tb4(nx+1,j)*xb/2.+tb4(nx+1,nym-j+2)*xb/2.
      xx3=tb4(nxm/2+1,j)*xb/2.+tb4(nxm/2+1,nym-j+2)*xb/2.
      xx4=tb4(nx+1,j)*xa/2.-tb4(nx+1,nym-j+2)*xa/2.
      tb6(nxm/2+1,2*j-1)=-(xx1-xx2)
      tb6(nxm/2+1,2*j)=-(xx3-xx4)
      enddo
      do i=1,nx
         tb6(i,ny+1)=2.*tb1(i,nym/2+1)
      enddo
      tb6(:,my)=0.
      tb6(mx,:)=0.
      tb6(mx-1,:)=0.
   endif
   if (isign==-1) then
      do j=1,ny+1
      do i=1,nxm/2+1
         tb1(i,j)=tb6(i,j)
         tb2(i,j)=-tb6(nxm-i+2,j)
      enddo
     enddo
     tb5(:,:)=0.; tb6(:,:)=0 
     do j=1,ny+1,2
     do i=1,nxm/2+1
        xb=sin((i-1)/2.*pi/(nxm))
        xa=cos((i-1)/2.*pi/(nxm))
        xb1=sin((nxm-i+1)/2.*pi/(nxm))
        xa1=cos((nxm-i+1)/2.*pi/(nxm))
        tb3(i,j)=-tb1(i,j+1)*xb-tb2(i,j+1)*xa
        tb4(i,j)=-tb1(i,j+1)*xa+tb2(i,j+1)*xb
        tb5(i,j)=tb1(i,j)*xa-tb2(i,j)*xb
        tb6(i,j)=tb1(i,j)*xb+tb2(i,j)*xa
     enddo
     enddo
     do j=1,ny
     do i=2,nxm/2
        tb3(i+nxm/2,j)=tb3(nxm/2-i+2,j)
        tb5(i+nxm/2,j)=tb5(nxm/2-i+2,j)
     enddo
     enddo
     do j=1,ny
     do i=2,nxm/2
        tb4(i+nxm/2,j)=-tb4(nxm/2-i+2,j)
        tb6(i+nxm/2,j)=-tb6(nxm/2-i+2,j)
     enddo
     enddo
     tb1(:,:)=0.;tb2(:,:)=0.; tb7(:,:)=0.;tb8(:,:)=0.
     do j=1,ny+1
     do i=1,nxm/2+1
        xa=cos((j-1)/2.*pi/(nym))
        xb=sin((j-1)/2.*pi/(nym))
        xa1=cos((nym-j+1)/2.*pi/(nym))
        xb1=sin((nym-j+1)/2.*pi/(nym))
        xx1=tb3(i,j)*xb-tb4(i,j)*xa
        xx2=tb3(i,nym-j+2)*xa+tb4(i,nym-j+2)*xb 
        xx3=-tb3(i,nym-j+2)*xa1+tb4(i,nym-j+2)*xb1
        xx4=tb3(i,j)*xa+tb4(i,j)*xb
        tb1(2*i-1,j)=xx4
        tb1(2*i,j)=xx1
        tb2(2*i-1,j)=xx3
        tb2(2*i,j)=xx2
     enddo
     enddo
     do i=1,mx
        tb1(i,1)=0.
     enddo
     do j=1,ny+1
     do i=1,nxm/2+1
        xa=cos((j-1)/2.*pi/(nym))
        xb=sin((j-1)/2.*pi/(nym))
        xa1=cos((nym-j+1)/2.*pi/(nym))
        xb1=sin((nym-j+1)/2.*pi/(nym))
        xx1=tb5(i,j)*xa-tb6(i,j)*xa1
        xx2=tb5(i,nym-j+2)*xb+tb6(i,nym-j+2)*xa 
        xx3=-tb5(i,nym-j+2)*xb1+tb6(i,nym-j+2)*xb
        xx4=tb5(i,j)*xb+tb6(i,j)*xa
        tb7(2*i-1,j)=xx1
        tb7(2*i,j)=xx4
        tb8(2*i-1,j)=xx2
        tb8(2*i,j)=xx3
     enddo
     enddo
     do i=1,mx
        tb7(i,ny)=tb7(i,(ny+1)/2)
        tb7(i,ny)=tb7(i,(ny+1)/2)
     enddo
     tb6(:,:)=0.
     do j=1,ny-1,2
     do i=1,mx
        tb6(i,(j+1)/2)=tb7(i,j)+tb1(i,j)         
     enddo
     enddo
     do j=1,ny-1,2
     do i=1,mx
        tb6(i,ny/2+(j+1)/2)=tb8(i,j)+tb2(i,j)
     enddo
     enddo
     tb5(:,:)=0.
     do i=1,mx-1,2
        tb5(i,1)=-(tb8(i+1,1)+tb2(i+1,1))
        tb5(i+1,1)=(tb8(i,1)+tb2(i,1))
     enddo
     do i=1,mx
        tb6(i,nym/2+1)=tb5(i,1)/2.
     enddo
     if (ifft==1) then
         call slfft2d(tb6,nxm,nym,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,1)
      endif
      if (ifft==2) then
         do j=1,nym
         do i=1,nxm
            tb6(i,j)=tbf(i,j)
         enddo
         enddo
      endif
     do j=1,ny
     do i=1,nxm/2
        ppm(2*i-1,j,1)=tb6(i,j)
        ppm(2*i,j,1)=tb6(nxm-i+1,j) 
     enddo
     enddo
  endif
endif
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
if (((ncly.eq.0)).and.(nclx.eq.0)) then
   if (isign==1) then    
      do j=1,nym
      do i=1,nxm
         tb1(i,j)=ppm(i,j,1)
      enddo
      enddo  
      if (ifft==1) then
         call slfft2d(tb1,nxm,nym,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,-1)
      endif
      if (ifft==2) then
         do j=1,nym
         do i=1,nxm
            tbf(i,j)=ppm(i,j,1)
         enddo
         enddo
         tb1=tb1/nxm/nym
      endif
      do j=1,nym
      do i=1,nxm/2+1
         tb2(i,j)=tb1(2*i-1,j)
         tb3(i,j)=tb1(2*i,j)
      enddo
      do i=2,nxm/2
         tb2(i+nxm/2,j)=tb2(nxm/2-i+2,j)
         tb3(i+nxm/2,j)=-tb3(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,my
         tb2(mx,j)=tb2(nxm/2+1,j)
         tb3(mx,j)=-tb3(nxm/2+1,j)
      enddo
      do j=1,nym/2+1
      do i=1,nxm/2+1
         xa=cos((j-1)*pi/(nym))
         xb=sin((j-1)*pi/(nym))
         xx1=tb2(i,j)*xa/2.+tb3(i,j)*xb/2.
         xx2=tb2(i,nym-j+2)*xa/2.-tb3(i,nym-j+2)*xb/2.    
         xx3=-tb2(i,j)*xb/2.+tb3(i,j)*xa/2.
         xx4=-(tb2(i,nym-j+2)*xb/2.+tb3(i,nym-j+2)*xa/2.)          
         tb4(i,2*j-1)=xx1+xx2
         tb4(i,2*j)=xx3+xx4
      enddo
      enddo
      do j=1,nym/2
         xa=cos((j-1)*pi/(nym))
         xb=sin((j-1)*pi/(nym))
         xa1=cos((nym-j+1)*pi/(nym))
         xb1=sin((nym-j+1)*pi/(nym))
         xx1=tb2(nxm/2+1,j)*xa/2.-tb3(nxm/2+1,nym-j+2)*xb/2.
         xx2=tb2(nxm/2+1,nym-j+2)*xa/2.+tb3(mx,nym-j+2)*xb/2.
         xx3=tb2(nxm/2+1,j)*xb/2.+tb3(nxm/2+1,nym-j+2)*xa/2.
         xx4=-tb2(nxm/2+1,nym-j+2)*xb/2.+tb3(mx,nym-j+2)*xa/2.             
         tb4(nxm/2+1,2*j-1)=xx1+xx2
         tb4(nxm/2+1,2*j)=-xx3+xx4
      enddo
      do j=1,my
      do i=2,nxm/2+2
         tb4(i+nxm/2,j)=tb4(nxm/2-i+2,j)
      enddo
      enddo
      tb4(nxm/2+1,:)=tb4(nxm/2+1,:)*0.5
      tb5(:,:)=0.
      do j=1,nym/2+1
      do i=1,nxm/2
         xa=cos((j-1)*pi/(nym))
         xb=sin((j-1)*pi/(nym))
         xa1=cos((nym-j+1)*pi/(nym))
         xb1=sin((nym-j+1)*pi/(nym))
         xx1=-tb2(i,j)*xb/2.+tb3(i,j)*xa/2.
         xx2=tb2(i,nym-j+2)*xb/2.+tb3(i,nym-j+2)*xa/2.     
         xx3=tb2(i,j)*xa/2.+tb3(i,j)*xb/2.
         xx4=-tb2(i,nym-j+2)*xa/2.+tb3(i,nym-j+2)*xb/2.       
         tb5(i,2*j-1)=xx1+xx2
         tb5(i,2*j)=-(xx3+xx4)
      enddo
      enddo
      do j=1,my
      do i=2,nxm/2+2
         tb5(i+nxm/2,j)=-tb5(nxm/2-i+2,j)
      enddo
      enddo
      do i=1,mx
         tb4(i,1)=tb4(i,1)*2.
         tb5(i,1)=tb5(i,1)*2.
      enddo
      tb6(:,:)=0.       
      do j=1,ny+2
      do i=1,nxm/2+2
         xa=cos((i-1)*pi/(nxm))
         xb=sin((i-1)*pi/(nxm))
         xx1=tb4(i,j)*xa/2.+tb4(nxm-i+2,j)*xa/2.
         xx2=-tb5(nxm-i+2,j)*xb/2.+tb5(i,j)*xb/2. 
         xx3=tb4(i,j)*xb/2.+tb4(nxm-i+2,j)*xb/2.
         xx4=tb5(nxm-i+2,j)*xa/2.+tb5(nxm-i+2,j)*xa/2. 
         tb6(2*i-1,j)=xx1+xx2
         tb6(2*i,j)=-(xx3+xx4)
      enddo
      enddo
      do j=1,my
         tb6(1,j)=tb6(1,j)
      enddo
      tb6(:,2)=0.
      tb6(mx,1)=tb6(mx,1)*0.5
   endif
   if (isign==-1) then
      do j=1,my
         tb6(nxm+1,j)=-tb6(mx,j)
      enddo
      do j=2,my
         tb6(nxm/2+1,j)=tb6(nxm/2+1,j)/2.
         tb6(nxm/2+2,j)=tb6(nxm/2+2,j)/2.
      enddo
      do j=1,my
      do i=1,nxm/2+1,2
         tb1(i,j)=tb6(i,j)
         tb1(i+1,j)=-tb6(i+1,j)
         tb2(i,j)=tb6(nxm-i+2,j)
         tb2(i+1,j)=-tb6(nxm-i+1,j)
      enddo
      enddo
      do j=1,my
         tb2(nxm/2+2,j)=0.
      enddo
      tb5(:,:)=0.;tb6(:,:)=0.;tb7(:,:)=0.;tb8(:,:)=0.
      do j=1,my
      do i=1,mx-1,2
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         tb3((i+1)/2,j)=tb1(i+1,j)*xb+tb2(nxm-i+1,j)*xb
         tb4((i+1)/2,j)=-tb1(i+1,j)*xa-tb2(nxm-i+1,j)*xa
         tb5((i+1)/2,j)=tb1(i,j)*xa+tb2(nxm-i+2,j)*xa
         tb6((i+1)/2,j)=tb1(i,j)*xb+tb2(nxm-i+2,j)*xb
      enddo
      enddo
      do j=1,my
      do i=2,nxm/2
         tb3(i+nxm/2,j)=tb3(nxm/2-i+2,j)
         tb5(i+nxm/2,j)=tb5(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,my
      do i=2,nxm/2
         tb4(i+nxm/2,j)=-tb4(nxm/2-i+2,j)
         tb6(i+nxm/2,j)=-tb6(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,my
         tb5(nxm/2+1,j)=tb6(nxm/2+1,j)
         tb6(nxm/2+1,j)=0.
      enddo
      tb1(:,:)=0.
      do j=1,ny+1,2
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/nym)
         xb=sin((j-1)/2.*pi/nym)
         xx1=(tb5(i,j)+tb5(nxm-i+2,j))*xa/2.-(tb6(i,j)-tb6(nxm-i+2,j))*xb/2.
         xx2=(tb5(i,j)+tb5(nxm-i+2,j))*xa/2.+(tb6(i,j)-tb6(nxm-i+2,j))*xb/2.
         xx3=(tb5(i,j+1)+tb5(nxm-i+2,j+1))*xa/2.-(tb6(i,j+1)-tb6(nxm-i+2,j+1))*xb/2.
         xx4=-((tb5(i,j+1)+tb5(nxm-i+2,j+1))*xa/2.+(tb6(i,j+1)-tb6(nxm-i+2,j+1))*xb/2.)
         xx5=-((tb5(i,j)+tb5(nxm-i+2,j))*xb/2.-(tb6(i,j)-tb6(nxm-i+2,j))*xa/2.)
         xx6=(tb5(i,j)+tb5(nxm-i+2,j))*xb/2.+(tb6(i,j)-tb6(nxm-i+2,j))*xa/2.
         xx7=-((tb5(i,j+1)+tb5(nxm-i+2,j+1))*xb/2.-(tb6(i,j+1)-tb6(nxm-i+2,j+1))*xa/2.)
         xx8=-((tb5(i,j+1)+tb5(nxm-i+2,j+1))*xb/2.+(tb6(i,j+1)-tb6(nxm-i+2,j+1))*xa/2.)
         yy1=(tb3(i,j)+tb3(nxm-i+2,j))*xa/2.-(tb4(i,j)-tb4(nxm-i+2,j))*xb/2.
         yy2=(tb3(i,j)+tb3(nxm-i+2,j))*xa/2.+(tb4(i,j)-tb4(nxm-i+2,j))*xb/2.
         yy3=-((tb3(i,j+1)+tb3(nxm-i+2,j+1))*xa/2.+(tb4(i,j+1)-tb4(nxm-i+2,j+1))*xb/2.)
         yy4=(tb3(i,j)+tb3(nxm-i+2,j))*xb/2.+(tb4(i,j)-tb4(nxm-i+2,j))*xa/2.
         yy5=-((tb3(i,j+1)+tb3(nxm-i+2,j+1))*xb/2.-(tb4(i,j+1)-tb4(nxm-i+2,j+1))*xa/2.)
         yy6=-((tb3(i,j+1)+tb3(nxm-i+2,j+1))*xb/2.+(tb4(i,j+1)-tb4(nxm-i+2,j+1))*xa/2.)
         yy7=(tb3(i,j+1)+tb3(nxm-i+2,j+1))*xa/2.-(tb4(i,j+1)-tb4(nxm-i+2,j+1))*xb/2.
         yy8=-((tb3(i,j)+tb3(nxm-i+2,j))*xb/2.-(tb4(i,j)-tb4(nxm-i+2,j))*xa/2.)
         tb1(2*i-1,(j+1)/2)=xx1+xx8+yy1+yy6
         tb1(2*i,(j+1)/2)=xx6+xx3+yy4+yy7
         tb1(2*i-1,nym-(j+1)/2+2)=xx2+xx7+yy2+yy5
         tb1(2*i,nym-(j+1)/2+2)=xx4+xx5+yy3+yy8
      enddo
      enddo
      tb1(nxm/2+1,1)=tb1(nxm/2+1,1)/2.
      tb1(nxm/2+2,1)=tb1(nxm/2+2,1)/2.
      tb1(mx,1)=0.
      tb1(nx+1,1)=tb1(nx+1,1)*2.
      tb1(:,ny+1)=0.
      tb1(1:2,:)=tb1(1:2,:)*2.
      tb1(nx+1:nx+2,:)=tb1(nx+1:nx+2,:)*2.
      if (ifft==1) then
         call slfft2d(tb1,nxm,nym,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,1)
      endif
      if (ifft==2) then
         do j=1,nym
         do i=1,nxm
            tb1(i,j)=tbf(i,j)
         enddo
         enddo
      endif
      do i=1,nxm
      do j=1,nym
         ppm(i,j,1)=tb1(i,j)
      enddo
      enddo
   endif
endif

return
end subroutine slfft2d_shift
end module slfft2d_shift_m