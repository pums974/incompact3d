!*****************************************************************************
!
subroutine slfft3d_shift(ppm,tab1,tr,us,ps,wk1,wk2,work,nx,nxm,ny,nym,nz,nzm,&
     mx,my,mz,nwork,table,isign)
!
!*****************************************************************************
 
USE param 

implicit none

integer :: mxyz,mxy,i,j,k,nx,mx,nxm,ny,my,nym,nzm,nz,mz,nwork,plan
real(8),dimension(nxm,nym,nzm) :: ppm,ppm1
real(8),dimension(mx,my,mz) :: tab1,tr,us,ps,wk1,wk2
real(8), dimension(nwork) :: work
real(8), dimension(100+2*(nxm+nym+nzm)) :: table
integer :: isign,nxyz,nxy,nyz,mxyz1,nvect1,myz
real(8) :: xx3,xx4,x,y,z
real(8) :: xa,xb,xx1,xx2,xa1,xb1
real(8) :: yy1,yy2,yy3,yy4,yy5,yy6,yy7,yy8,xx5,xx6,xx7,xx8,xtoto

mxyz=mx*my*mz
mxyz1=nxm*nym*nzm
nvect1=nym*nzm
mxy=mx*my
myz=my*mz
nyz=ny*nz
nxy=nxm*(nz+1)
nxyz=nx*(nx+1)*(ny/2)


if (((ncly.eq.2).or.(ncly.eq.1)).and.(nclx.eq.0).and.(nclz.eq.0)) then
   if (isign==1) then
      tr(:,:,:)=0.;us(:,:,:)=0.;ps(:,:,:)=0.
      wk1(:,:,:)=0.
      do k=1,nzm
      do j=1,nym/2
      do i=1,nxm
        wk1(i,j,k)=ppm(i,2*(j-1)+1,k)
      enddo
      enddo
      do j=nym/2+1,nym
      do i=1,nxm
         wk1(i,j,k)=ppm(i,2*nym-2*j+2,k)
      enddo
      enddo
      enddo

      call scfft3d(0,nxm,nym,nzm,1./nxm/nym/nzm,wk1,mx,my,wk1,mx/2,my,table,work,0)
      call scfft3d(1,nxm,nym,nzm,1./nxm/nym/nzm,wk1,mx,my,wk1,mx/2,my,table,work,0)

      do k=1,nzm
      do j=1,nym
      do i=1,nxm/2
         tr(i,j,k)=wk1(2*i-1,j,k)
         us(i,j,k)=-wk1(2*i,j,k)
      enddo
      tr(nxm/2+1,j,k)=wk1(nx+1,j,k)
      tr(nxm/2+2,j,k)=wk1(nx+2,j,k)
      us(nxm/2+1,j,k)=-wk1(nx+1,j,k)
      us(nxm/2+2,j,k)=-wk1(nx+2,j,k)
      do i=3,nxm/2
         tr(i+nxm/2,j,k)=tr(nxm/2-i+2,j,k)
         us(i+nxm/2,j,k)=-us(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      wk1(:,:,:)=0. ; wk2(:,:,:)=0.
!**************transformation en y***************************
      do k=1,nz     
      do j=1,ny
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xx1=tr(i,j,k)*xa/2.+us(i,j,k)*xb/2.
         xx2=tr(i,nym-j+2,k)*xa/2.-us(i,nym-j+2,k)*xb/2.           
         wk1(i,j,k)=xx1+xx2      
      enddo
      do i=2,nxm/2
         wk1(i+nxm/2,j,k)=wk1(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xx1=tr(nxm/2+1,j,k)*xa/2.+tr(nxm/2+1,nym-j+2,k)*xa/2.
         xx2=-tr(nxm/2+2,j,k)*xb/2.+tr(nxm/2+2,nym-j+2,k)*xb/2.           
         wk1(nxm/2+1,j,k)=xx1+xx2 
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xx1=-tr(i,j,k)*xb/2.+us(i,j,k)*xa/2.
         xx2=tr(i,nym-j+2,k)*xb/2.+us(i,nym-j+2,k)*xa/2.           
         wk2(i,j,k)=xx1+xx2 
      enddo
      do i=2,nxm/2
         wk2(i+nxm/2,j,k)=-wk2(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=2,ny
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xx1=tr(nxm/2+1,j,k)*xa/2.-tr(nxm/2+1,nym-j+2,k)*xa/2.
         xx2=-tr(nxm/2+2,j,k)*xb/2.-tr(nxm/2+2,nym-j+2,k)*xb/2.
         wk2(nxm/2+1,nym-j+2,k)=xx1+xx2 
      enddo
      enddo
      do k=1,nz
         xb=cos((j-1)/2.*pi/(nym))
         xa=sin((j-1)/2.*pi/(nym))
         xx1=-tr(nxm/2+2,1,k)*xa/2.-tr(nxm/2+1,nym-1+2,k)*xa/2.
         xx2=tr(nxm/2+1,nym-1+2,k)*xb/2.+tr(nxm/2+2,1,k)*xb/2.           
         wk2(nxm/2+1,1,k)=xx1+xx2 
      enddo
      do k=1,mz
      do i=1,mx
         wk1(i,1,k)=wk1(i,1,k)*2.
         wk2(i,1,k)=wk2(i,1,k)*2.
         wk1(i,ny,k)=wk1(i,ny,k)*2.
         wk2(i,ny,k)=wk2(i,ny,k)*2.      
      enddo
      enddo
!********************fin de la tranformation en y***********************************
      do i=1,mxyz
         tr(i,1,1)=0.
         us(i,1,1)=0.
         ps(i,1,1)=0.
      enddo
!********************transformation en x********************************************
      do k=1,nz
      do j=1,ny
      do i=1,nxm/2
         xb=sin((i-1)*pi/(nxm))
         xa=cos((i-1)*pi/(nxm))
         xx1=(wk1(i,j,k)+wk1(nxm-i+2,j,k))*xa/2
         xx2=(wk2(i,j,k)-wk2(nxm-i+2,j,k))*xb/2.
         xx3=(wk1(i,j,k)+wk1(nxm-i+2,j,k))*xb/2
         xx4=(-wk2(i,j,k)+wk2(nxm-i+2,j,k))*xa/2.
         ps(2*i-1,j,k)=xx1+xx2
         ps(2*i,j,k)=-(xx3+xx4)
      enddo
      enddo
      enddo
      do k=1,mz
         ps(nx+1,1,k)=wk1(nxm/2+1,1,k)
         ps(nx+2,1,k)=-wk1(nxm/2+1,ny,k)
      do j=2,my
         ps(nx+1,j,k)=wk1(nxm/2+1,j,k)
         ps(nx+2,j,k)=wk2(nxm/2+1,j,k)
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         ps(1,j,k)=2.*ps(1,j,k)
         ps(2,j,k)=2.*ps(2,j,k)
      enddo
      enddo
      do k=1,mz
      do i=1,mx
         ps(i,ny,k)=0.
      enddo
      enddo
!*******************fin de la transformation en x*********************************
      tab1(:,:,:)=0.
!*******************tranformation en z********************************************
      do j=1,my
      do i=1,mx
         tr(i,j,nzm/2+1)=ps(i,j,nzm/2+1)
      enddo
      enddo
      do k=1,nzm/2+1
      do j=1,nym
      do i=1,nx+1,2
         xa=cos((k-1)*pi/(nzm))
         xb=sin((k-1)*pi/(nzm))
         xx1=ps(i,j,k)*xa/2.+ps(i+1,j,k)*xb/2.
         xx2=ps(i,j,nzm-k+2)*xa/2.-ps(i+1,j,nzm-k+2)*xb/2.
         xx3=ps(i,j,k)*xb/2.-ps(i+1,j,k)*xa/2.
         xx4=ps(i,j,nzm-k+2)*xb/2.+ps(i+1,j,nzm-k+2)*xa/2. 
         tab1(i,j,2*k-1)=xx1+xx2
         tab1(i,j,2*k)=-(xx3+xx4)
         tab1(i+1,j,2*k-1)=xx4-xx3
         tab1(i+1,j,2*k)=xx2-xx1
      enddo
      enddo
      enddo
      do k=1,nzm/2+1
      do j=1,ny
         xa=cos((k-1)*pi/(nzm))
         xb=sin((k-1)*pi/(nzm))
         tab1(nx+2,j,2*k-1)=ps(nx+1,j,k)*xa+ps(nx+2,j,k)*xb
         tab1(nx+2,j,2*k)=ps(nx+1,j,nzm-k+2)*xb+ps(nx+2,j,nzm-k+2)*xa
      enddo
      enddo
      do i=1,mxyz
         wk1(i,1,1)=0.
         wk2(i,1,1)=0.
         ps(i,1,1)=0.
      enddo
      do i=1,mxyz
         ps(i,1,1)=tab1(i,1,1)
      enddo
      ps(:,:,2)=0.
      ps(:,2:ny,1)=2.*ps(:,2:ny,1)
      ps(2:nx,1,1)=2.*ps(2:nx,1,1)
      ps(nx+1,:,:)=0.
      ps(nxm/2+2,1,mz)=ps(nxm/2+2,1,mz)*2.
      ps(1,1,1)=2.*ps(1,1,1)
!*******************fin de la transformation en z********************************* 

!      stop
      tab1(:,:,:)=0. ; tab1(:,:,:)=ps(:,:,:)
   endif
   if (isign==-1) then
      wk1(:,:,:)=0. ; wk1(:,:,:)=tab1(:,:,:)
      wk1(mx,1,1)=2.*wk1(mx,1,1)
      wk1(mx,:,1)=wk1(mx,:,1)/2.
      wk1(nxm/2+1,:,mz)=wk1(nxm/2+1,:,mz)/2.
      wk1(nxm/2+2,:,mz)=wk1(nxm/2+2,:,mz)/2.      
      wk1(nxm/2+1,1,mz)=wk1(nxm/2+1,1,mz)*2.
      do k=1,nz+1,2
      do j=1,my
         wk1(nxm+1,j,k+1)=-wk1(mx,j,k+1)
         wk1(nxm+1,j,k)=wk1(mx,j,k)
      enddo
      enddo
      wk1(nxm/2+1,1,mz)=wk1(nxm/2+1,1,mz)/2.
      wk1(nxm/2+2,1,mz)=wk1(nxm/2+2,1,mz)/2.
      do i=1,mxyz
         tr(i,1,1)=0.
         us(i,1,1)=0.
         wk2(i,1,1)=0.
         ps(i,1,1)=0.
      enddo
      us=0. ; tr=0.
      do k=1,mz
      do j=1,ny
      do i=1,nxm/2+1,2
         tr(i,j,k)=wk1(i,j,k)
         tr(i+1,j,k)=-wk1(i+1,j,k)
         us(i,j,k)=wk1(nxm-i+2,j,k)
         us(i+1,j,k)=-wk1(nxm-i+1,j,k)
      enddo
      enddo
      enddo
      do k=1,mz
      do j=1,my
         us(nxm/2+2,j,k)=0.
      enddo
      enddo
      wk1(:,:,:)=0. ; tab1(:,:,:)=0.
      do k=1,mz
      do j=1,ny+1
      do i=1,mx-1,2
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         wk1((i+1)/2,j,k)=tr(i+1,j,k)*xb+us(nxm-i+1,j,k)*xb
         wk2((i+1)/2,j,k)=-tr(i+1,j,k)*xa-us(nxm-i+1,j,k)*xa
         ps((i+1)/2,j,k)=tr(i,j,k)*xa+us(nxm-i+2,j,k)*xa
         tab1((i+1)/2,j,k)=tr(i,j,k)*xb+us(nxm-i+2,j,k)*xb
      enddo
      enddo
      enddo
      do k=1,nz+1
      do j=1,ny
      do i=2,nxm/2
         wk1(i+nxm/2,j,k)=wk1(nxm/2-i+2,j,k)
         ps(i+nxm/2,j,k)=ps(nxm/2-i+2,j,k)
      enddo
      enddo
      do j=1,ny
      do i=2,nxm/2
         wk2(i+nxm/2,j,k)=-wk2(nxm/2-i+2,j,k)
         tab1(i+nxm/2,j,k)=-tab1(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      do k=1,mz
      do j=1,my
         ps(nxm/2+1,j,k)=tab1(nxm/2+1,j,k)
         tab1(nxm/2+1,j,k)=0.
      enddo
      enddo
      tr(:,:,:)=0.
      do k=1,mz
      do j=1,ny+1
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xa1=cos((nym-j+1)/2.*pi/(nym))
         xb1=sin((nym-j+1)/2.*pi/(nym))
         xx1=(ps(i,j,k)+wk1(i,j,k))*xa-(tab1(i,j,k)+wk2(i,j,k))*xa1
         xx2=(ps(i,nym-j+2,k)+wk1(i,nym-j+2,k))*xb+(tab1(i,nym-j+2,k)+wk2(i,nym-j+2,k))*xa 
         xx3=-(ps(i,nym-j+2,k)+wk1(i,nym-j+2,k))*xb1+(tab1(i,nym-j+2,k)+wk2(i,nym-j+2,k))*xb
         xx4=(ps(i,j,k)+wk1(i,j,k))*xb+(tab1(i,j,k)+wk2(i,j,k))*xa
         tr(2*i-1,j,k)=xx1+xx2
         tr(2*i,j,k)=-(xx4+xx3)
      enddo
      enddo
      enddo
      do k=1,nz+1
      do j=1,ny
      do i=1,mx
        tr(i,ny,k)=0.
      enddo
      tr(nxm/2+1,j,k)=tr(nxm/2+1,j,k)/2.
      tr(nxm/2+2,j,k)=tr(nxm/2+2,j,k)/2.
      tr(mx,1,k)=0.
      enddo
      enddo
      tr(:,my,:)=0.
      ps(:,:,:)=0. ; ps(:,:,:)=tr(:,:,:)
      ps(nz+1:mz,:,nz+1)=0.
!************TRANSFO EN Z*************************************************
      tr(:,:,:)=0. 
      do k=1,nz+1,2
      do j=1,ny
      do i=1,nx+1,2
         xa=cos((k-1)/2.*pi/(nzm))
         xb=sin((k-1)/2.*pi/(nzm))
         xx1=(ps(i,j,k)*xa+ps(i+1,j,k)*xb)+(-ps(i,j,k+1)*xb+ps(i+1,j,k+1)*xa)
         xx2=(-ps(i,j,k)*xb+ps(i+1,j,k)*xa)-(ps(i,j,k+1)*xa+ps(i+1,j,k+1)*xb) 
         xx3=(ps(i,j,k)*xa-ps(i+1,j,k)*xb)-(ps(i,j,k+1)*xb+ps(i+1,j,k+1)*xa) 
         xx4=(ps(i,j,k)*xb+ps(i+1,j,k)*xa)+(ps(i,j,k+1)*xa-ps(i+1,j,k+1)*xb) 
         tr(i,j,(k+1)/2)=xx1
         tr(i+1,j,(k+1)/2)=xx2
         tr(i,j,nzm-(k+1)/2+2)=xx3
         tr(i+1,j,nzm-(k+1)/2+2)=xx4
      enddo
      enddo
      enddo
      tr(:,:,nz+1)=0.
      ps(:,:,:)=0. ; ps(:,:,:)=tr(:,:,:)
      ps(:,ny,nzm/2+1)=0.
   
      call csfft3d(-1,nxm,nym,nzm,1.,ps,mx/2,my,ps,mx,my,table,work,0)
      
      do k=1,nzm
      do j=1,nym/2   
      do i=1,nxm
         ppm(i,2*j-1,k)=ps(i,j,k)
      enddo
      enddo
      enddo
      do k=1,nzm
      do j=1,nym/2
      do i=1,nxm
         ppm(i,2*j,k)=ps(i,nym-j+1,k)
      enddo
      enddo
      enddo 
   endif
endif

if (((nclx.eq.2).or.(nclx.eq.1)).and.(ncly.eq.0).and.(nclz.eq.0)) then
   if (isign==1) then
      tr(:,:,:)=0. ; us(:,:,:)=0. ; ps(:,:,:)=0.
      wk1(:,:,:)=0. ; wk2(:,:,:)=0.       
!**********construction en x************************
      do k=1,nzm
      do j=1,nym
      do i=1,nxm/2
         wk1(i,j,k)=ppm(2*(i-1)+1,j,k)
      enddo
      do i=nxm/2+1,nxm
         wk1(i,j,k)=ppm(2*nxm-2*i+2,j,k)
      enddo
      enddo
      enddo

      call scfft3d(0,nxm,nym,nzm,1./nxm/nym/nzm,wk1,mx,ny+2,wk1,nxm/2+1,ny+2,table,work,0)
      call scfft3d(1,nxm,nym,nzm,1./nxm/nym/nzm,wk1,mx,ny+2,wk1,nxm/2+1,ny+2,table,work,0)!

      do k=1,nzm
      do j=1,nym
      do i=1,nxm/2
         tr(i,j,k)=wk1(2*i-1,j,k)
         us(i,j,k)=-wk1(2*i,j,k)
      enddo
      tr(nxm/2+1,j,k)=wk1(nxm+1,j,k)
      us(nxm/2+1,j,k)=-wk1(nx+1,j,k)
      do i=2,nxm/2
         tr(i+nxm/2,j,k)=tr(nxm/2-i+2,j,k)
      enddo
      do i=2,nxm/2
         us(i+nxm/2,j,k)=-us(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      do k=1,nzm
      do j=1,nym
      do i=1,nxm/2
         tr(i,j,k)=wk1(2*i-1,j,k)
         us(i,j,k)=-wk1(2*i,j,k)
      enddo
      tr(nxm/2+1,j,k)=wk1(nxm+1,j,k)
      us(nxm/2+1,j,k)=-wk1(nx+1,j,k)
      do i=2,nxm/2
         tr(i+nxm/2,j,k)=tr(nxm/2-i+2,j,k)
      enddo
      do i=2,nxm/2
         us(i+nxm/2,j,k)=-us(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      wk1(:,:,:)=0.
!**************transformation en y***************************
      do k=1,nz 
      do j=1,ny/2+1
      do i=1,nxm/2+1
         xa=cos((j-1)*pi/(nym))
         xb=sin((j-1)*pi/(nym))
         xx1=tr(i,j,k)*xa/2.+us(i,j,k)*xb/2.
         xx2=tr(i,nym-j+2,k)*xa/2.-us(i,nym-j+2,k)*xb/2.      
         xx3=tr(i,j,k)*xb/2-us(i,j,k)*xa/2.
         xx4=tr(i,nym-j+2,k)*xb/2.+us(i,nym-j+2,k)*xa/2.
         wk1(i,2*j-1,k)=xx1+xx2      
         wk1(i,2*j,k)=-(xx3+xx4)
      enddo
      do i=1,nxm/2+1
         xa=cos((j-1)*pi/(nym))
         xb=sin((j-1)*pi/(nym))
         xa1=cos((nym-j+1)*pi/(nym))
         xb1=sin((nym-j+1)*pi/(nym))
         xx1=-tr(i,j,k)*xb/2.+us(i,j,k)*xa/2.
         xx2=tr(i,nym-j+2,k)*xb/2.+us(i,nym-j+2,k)*xa/2.           
         xx3=tr(i,j,k)*xa/2+us(i,j,k)*xb/2.
         xx4=-tr(i,nym-j+2,k)*xa/2.+us(i,nym-j+2,k)*xb/2.
         wk2(i,2*j-1,k)=xx1+xx2
         wk2(i,2*j,k)=-(xx3+xx4)
      enddo
      enddo
      do j=1,ny 
      do i=2,nxm/2
         wk1(i+nxm/2,j,k)=wk1(nxm/2-i+2,j,k)
      enddo
      enddo
      do j=1,ny+1
      do i=2,nxm/2
         wk2(i+nxm/2,j,k)=-wk2(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         wk1(i,1,k)=wk1(i,1,k)*2.
         wk2(i,1,k)=wk2(i,1,k)*2.
         wk1(i,ny,k)=wk1(i,ny,k)*2.
         wk2(i,ny,k)=wk2(i,ny,k)*2.
      enddo
      enddo
      do i=1,mxyz
         tr(i,1,1)=0.
         us(i,1,1)=0.
         ps(i,1,1)=0.
      enddo
!********************transformation en x********************************************
      do k=1,nz
      do j=1,ny+2
      do i=1,nx+1
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         xx1=(wk1(i,j,k)+wk1(nxm-i+2,j,k))*xa/2
         xx2=(wk2(i,j,k)-wk2(nxm-i+2,j,k))*xb/2.
         ps(i,j,k)=xx1+xx2
      enddo
      xb=sin((nxm/4.)*pi/(nxm))
      xa=cos((nxm/4.)*pi/(nxm))
      ps(nx+1,j,k)=-((wk2(nxm/2+1,j,k)*xa/2.)+(wk2(nxm/2+1,j,k)*xb/2.))
      ps(1,j,k)=2.*ps(1,j,k)
      ps(nx,j,k)=2.*ps(nx,j,k)
      enddo
      enddo
      ps(:,ny,:)=ps(:,ny,:)/2.
!*******************fin de la transformation en x*********************************
      wk1(:,:,:)=0. ; wk2(:,:,:)=0.
      do k=1,mz
      do i=1,mx
         ps(i,ny+1,k)=-ps(i,ny+2,k)*2.
      enddo
      enddo
      ps(1,ny+1,1)=ps(1,ny+1,1)/2.
      do k=1,mz
      do i=1,mx
         ps(i,my,k)=0.
      enddo
      enddo
!*******************tranformation en z********************************************
      do j=1,my
      do i=1,mx
         tr(i,j,nzm/2+1)=ps(i,j,nzm/2+1)
      enddo
      enddo
      wk1(:,:,:)=0.
      do k=1,nzm/2
      do j=1,nym+1
      do i=1,nxm/2
         xa=cos((k-1)*pi/(nzm))
         xb=sin((k-1)*pi/(nzm))
         xx1=ps(i,j,k)*xa/2.+ps(i,j,nzm-k+2)*xa/2.
         xx2=-ps(nxm-i+2,j,k)*xb/2.+ps(nxm-i+2,j,nzm-k+2)*xb/2. 
         xx3=ps(i,j,k)*xb/2.+ps(i,j,nzm-k+2)*xb/2.
         xx4=ps(nxm-i+2,j,k)*xa/2.-ps(nxm-i+2,j,nzm-k+2)*xa/2. 
         wk1(i,j,2*k-1)=xx1+xx2
         wk1(i,j,2*k)=-(xx3+xx4)
      enddo
      do i=nxm/2+1,nx+1
         xa=cos((nzm-k+1)*pi/(nzm))
         xb=sin((nzm-k+1)*pi/(nzm))
         xx1=ps(i,j,k)*xa/2.+ps(i,j,nzm-k+2)*xa/2.
         xx2=-ps(nxm-i+2,j,k)*xb/2.+ps(nxm-i+2,j,nzm-k+2)*xb/2.
         xx3=ps(i,j,k)*xb/2.+ps(i,j,nzm-k+2)*xb/2.
         xx4=ps(nxm-i+2,j,k)*xa/2.-ps(nxm-i+2,j,nzm-k+2)*xa/2.  
         wk1(i,j,2*k-1)=-(xx1+xx2)
         wk1(i,j,2*k)=-(xx3+xx4) 
      enddo
      xa=cos((nzm-k+1)*pi/(nzm))
      xb=sin((nzm-k+1)*pi/(nzm))
      xx1=ps(nxm/2+1,j,k)*xa/2.+ps(nxm/2+1,j,nzm-k+2)*xa/2.
      xx2=-ps(nx+1,j,k)*xb/2.+ps(nx+1,j,nzm-k+2)*xb/2.
      xx3=ps(nxm/2+1,j,k)*xb/2.+ps(nxm/2+1,j,nzm-k+2)*xb/2.
      xx4=ps(nx+1,j,k)*xa/2.-ps(nx+1,j,nzm-k+2)*xa/2.
      wk1(nxm/2+1,j,2*k-1)=-(xx1-xx2)
      wk1(nxm/2+1,j,2*k)=-(xx3-xx4)    
      enddo
      enddo
      wk1(:,ny+1,1)=ps(:,ny+1,1)
      ps(:,:,:)=0. ; ps(:,:,:)=wk1(:,:,:)
      ps(:,:,nz+1)=tr(:,:,nzm/2+1)
      ps(:,:,2)=0.
      ps(:,2:ny,1)=2.*ps(:,2:ny,1)
      ps(2:nx,1,1)=2.*ps(2:nx,1,1)
      ps(nx+1,:,:)=0.
      ps(:,2,:)=0.
      ps(1,ny+1,2:mz)= ps(1,ny+1,2:mz)/2.
      ps(nxm/2+1,ny+1,2:mz)=ps(nxm/2+1,ny+1,2:mz)*0.5
      ps(mx,:,:)=0.; ps(mx-1,:,:)=0.; ps(:,my,:)=0.; ps(:,:,mz)=0.
      tab1(:,:,:)=0. ; tab1(:,:,:)=ps(:,:,:)

   endif
   if (isign==-1) then
      wk1(:,:,:)=0. ; wk1(:,:,:)=tab1(:,:,:)
      do i=1,mxyz
         tr(i,1,1)=0.
         us(i,1,1)=0.
         wk2(i,1,1)=0.
         ps(i,1,1)=0.
      enddo
      do k=1,nz+1
      do j=1,ny+1
      do i=1,nxm/2+1
         tr(i,j,k)=wk1(i,j,k)
         us(i,j,k)=-wk1(nxm-i+2,j,k)
      enddo
      enddo
      enddo
      do i=1,mxyz
         wk1(i,1,1)=0.
         tab1(i,1,1)=0.
      enddo
!*****************transformation en x**************************
      do k=1,nz+1,2
      do j=1,ny+1
      do i=1,nxm/2+1
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         xb1=sin((nxm-i+1)/2.*pi/(nxm))
         xa1=cos((nxm-i+1)/2.*pi/(nxm))
         wk1(i,j,k)=-tr(i,j,k+1)*xb-us(i,j,k+1)*xa
         wk2(i,j,k)=-(tr(i,j,k+1)*xa-us(i,j,k+1)*xb)
         ps(i,j,k)=tr(i,j,k)*xa-us(i,j,k)*xb
         tab1(i,j,k)=tr(i,j,k)*xb+us(i,j,k)*xa
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny+1
      do i=2,nxm/2
         wk1(i+nxm/2,j,k)=wk1(nxm/2-i+2,j,k)
         ps(i+nxm/2,j,k)=wk1(nxm/2-i+2,j,k)
         wk2(i+nxm/2,j,k)=-wk2(nxm/2-i+2,j,k)
         tab1(i+nxm/2,j,k)=-wk2(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      wk1(:,:,nz)=0. ; wk2(:,:,nz)=0. ; ps(:,:,nz)=0. ; tab1(:,:,nz)=0.
      tr(:,:,:)=0.; tr(:,:,:)=tab1(:,:,:)
      us(:,:,:)=0. ; tab1(:,:,:)=0.
      xtoto=ps(nxm/2+1,1,1)
!***************debut transformation en z****************************
      do k=1,nz+1
      do j=1,ny+1
      do i=1,nxm/2+1
         xa=cos((k-1)/2.*pi/(nzm))
         xb=sin((k-1)/2.*pi/(nzm))
         xa1=cos((nzm-k+1)/2.*pi/(nzm))
         xb1=sin((nzm-k+1)/2.*pi/(nzm))
         xx1=wk1(i,j,k)*xb-wk2(i,j,k)*xa
         xx2=wk1(i,j,nzm-k+2)*xa+wk2(i,j,nzm-k+2)*xb 
         xx3=-wk1(i,j,nzm-k+2)*xa1+wk2(i,j,nzm-k+2)*xb1
         xx4=wk1(i,j,k)*xa+wk2(i,j,k)*xb
         xx5=ps(i,j,k)*xa-tr(i,j,k)*xa1
         xx6=ps(i,j,nzm-k+2)*xb+tr(i,j,nzm-k+2)*xa 
         xx7=-ps(i,j,nzm-k+2)*xb1+tr(i,j,nzm-k+2)*xb
         xx8=ps(i,j,k)*xb+tr(i,j,k)*xa
         tab1(2*i-1,j,k)=xx8+xx1
         tab1(2*i,j,k)=xx5+xx4
         us(2*i-1,j,k)=xx7+xx2
         us(2*i,j,k)=xx6+xx3
      enddo
      enddo
      enddo
      wk1(:,:,:)=0. ; wk2(:,:,:)=0.
      do k=1,nz,2
      do j=1,my
      do i=1,mx
         wk1(i,j,(k+1)/2)=tab1(i,j,k)
         wk1(i,j,nz/2+(k+1)/2)=us(i,j,k)
      enddo
      enddo
      enddo
      ps(:,:,:)=0.
      wk1(:,2,:)=0.
      wk1(:,my,:)=-wk1(:,ny+1,:)
      wk1(:,ny+1,:)=0.
      do k=1,nz
      do j=1,ny+1,2
      do i=1,nx+1,2
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xx1=-(wk1(i,j,k)*xb-wk1(i+1,j,k)*xa)-(wk1(i,j+1,k)*xa+wk1(i+1,j+1,k)*xb)
         xx2=-(wk1(i+1,j,k)*xb+wk1(i,j,k)*xa)-(wk1(i+1,j+1,k)*xa-wk1(i,j+1,k)*xb)
         xx3=(wk1(i,j,k)*xb+wk1(i+1,j,k)*xa)+(wk1(i,j+1,k)*xa-wk1(i+1,j+1,k)*xb)
         xx4=(wk1(i+1,j,k)*xb-wk1(i,j,k)*xa)+(wk1(i+1,j+1,k)*xa+wk1(i,j+1,k)*xb)
         ps(i,nym-(j+1)/2+2,k)=xx3
         ps(i+1,nym-(j+1)/2+2,k)=xx4
         ps(i,(j+1)/2,k)=xx1
         ps(i+1,(j+1)/2,k)=xx2
      enddo
      enddo
      enddo
!****************************************
      tr(:,:,:)=0.
      tr(:,:,:)=ps(:,:,:)
      ps(1,:,nzm/2+1)=tr(2,:,nzm/2+1)
      ps(2,:,nzm/2+1)=-tr(1,:,nzm/2+1)
      do k=1,mz
      do i=nxm/2+1,nxm,2
         ps(i,nym/2+1,k)=tr(i+1,nym/2+1,k)
         ps(i+1,nym/2+1,k)=-tr(i,nym/2+1,k)
      enddo
      enddo
      ps(nx,nym/2+1,1)=-ps(nx,nym/2+1,1)
      do j=1,my
      do i=3,nx+1,2
         ps(i,j,nzm/2+1)=tr(i+1,j,nzm/2+1)
         ps(i+1,j,nzm/2+1)=-tr(i,j,nzm/2+1) 
      enddo
      enddo
      tr(:,:,:)=ps(:,:,:)
      do i=nxm/2+1,nx,2
         ps(i,nym/2+1,:)=-tr(i+1,nym/2+1,:)
         ps(i+1,nym/2+1,:)=tr(i,nym/2+1,:)
      enddo
      ps(nx,nym/2+1,:)=tr(nx,nym/2+1,:)
      ps(nx+1,nym/2+1,:)=tr(nx+1,nym/2+1,:)
      do i=nxm/2+1,nx-1,2
         ps(i,nym/2+1,nzm/2+1)=tr(i,nym/2+1,nzm/2+1)
         ps(i+1,nym/2+1,nzm/2+1)=tr(i+1,nym/2+1,nzm/2+1)
      enddo
      ps(nx,nym/2+1,nzm/2+1)=ps(nx+1,nym/2+1,nzm/2+1)
      ps(nx+1,nym/2+1,nzm/2+1)=0.
      ps(nx,1,1)=xtoto 
      ps(:,1:ny,nz+1)=0.
      ps(:,ny+1,1:nz)=0.
      ps(2,1,nzm/2+1)=0.
      ps(nx,nym/2+1,1)=-ps(nx,nym/2+1,1)/2.
      ps(nx+1,nym/2+1,1)=0.

      call csfft3d(-1,nxm,nym,nzm,1.,ps,nxm/2+1,my,ps,mx,my,table,work,0)!

      do k=1,nzm
      do j=1,nym
      do i=1,nxm/2
         ppm(2*i-1,j,k)=ps(i,j,k)
         ppm(2*i,j,k)=ps(nxm-i+1,j,k)  
      enddo
      enddo
      enddo
   endif
endif
!********************************************************************************************************
!********************************************************************************************************
!********************************************************************************************************
!********************************************************************************************************
!********************************************************************************************************
if (((nclx.eq.2).or.(nclx.eq.1)).and.((ncly.eq.1).or.(ncly.eq.2)).and.((nclz.eq.1).or.(nclz.eq.2))) then
   if (isign==1) then
      tr(:,:,:)=0. ; us(:,:,:)=0. ; ps(:,:,:)=0.
      wk1(:,:,:)=0. ; wk2(:,:,:)=0. 
      do k=1,nzm
      do j=1,nym
      do i=1,nxm/2
         tr(i,j,k)=ppm(2*(i-1)+1,j,k)
      enddo
      do i=nxm/2+1,nxm
         tr(i,j,k)=ppm(2*nxm-2*i+2,j,k)
      enddo
      enddo
      enddo
      do k=1,nzm
      do i=1,nxm
      do j=1,nym/2
        us(i,j,k)=tr(i,2*(j-1)+1,k)
      enddo
      do j=nym/2+1,nym
         us(i,j,k)=tr(i,2*nym-2*j+2,k)
      enddo
      enddo
      enddo
      do j=1,nym
      do i=1,nxm
      do k=1,nzm/2
        ps(i,j,k)=us(i,j,2*(k-1)+1)
      enddo
      do k=nzm/2+1,nzm
         ps(i,j,k)=us(i,j,2*nzm-2*k+2)
      enddo
      enddo
      enddo

      call scfft3d(0,nxm,nym,nzm,1./nxm/nym/nzm,ps,mx,my,ps,mx/2,my,table,work,0)
      call scfft3d(1,nxm,nym,nzm,1./nxm/nym/nzm,ps,mx,my,ps,mx/2,my,table,work,0)
  
      do k=1,nzm
      do j=1,nym
      do i=1,nxm/2
         tr(i,j,k)=ps(2*i-1,j,k)
         us(i,j,k)=-ps(2*i,j,k)
      enddo
      tr(nxm/2+1,j,k)=ps(nxm+1,j,k)
      us(nxm/2+1,j,k)=-ps(nx+1,j,k)
      do i=2,nxm/2
         tr(i+nxm/2,j,k)=tr(nxm/2-i+2,j,k)
         us(i+nxm/2,j,k)=-us(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo   
!**************transformation en z***************************
      do k=1,nz
      do j=1,ny
      do i=1,nxm/2+1
         xa=cos((k-1)/2.*pi/(nzm))
         xb=sin((k-1)/2.*pi/(nzm))
         xx1=tr(i,j,k)*xa/2.+us(i,j,k)*xb/2.
         xx2=tr(i,j,nzm-k+2)*xa/2.-us(i,j,nzm-k+2)*xb/2.           
         xx3=-tr(i,j,k)*xb/2.+us(i,j,k)*xa/2.
         xx4=tr(i,j,nzm-k+2)*xb/2.+us(i,j,nzm-k+2)*xa/2.           
         wk2(i,j,k)=(xx3+xx4) 
         wk1(i,j,k)=xx1+xx2      
      enddo
      do i=2,nxm/2
         wk1(i+nxm/2,j,k)=wk1(nxm/2-i+2,j,k)
         wk2(i+nxm/2,j,k)=-wk2(nxm/2-i+2,j,k) 
      enddo
      enddo
      enddo    
      do j=1,nym
      do i=1,nxm
         wk1(i,j,1)=wk1(i,j,1)*2.
         wk2(i,j,1)=wk2(i,j,1)*2.
      enddo
      enddo
      tr(:,:,:)=0. ; us(:,:,:)=0.
!*************************transformation en y***************************************
      do k=1,nz
      do j=1,ny 
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xa1=cos((nym-j+1)/2.*pi/(nym))
         xb1=sin((nym-j+1)/2.*pi/(nym))
         xx1=wk1(i,j,k)*xa/2.+wk2(i,j,k)*xb/2.
         xx2=wk1(i,nym-j+2,k)*xa/2.-wk2(i,nym-j+2,k)*xb/2.
         xx3=-wk1(i,j,k)*xb/2.+wk2(i,j,k)*xa/2.
         xx4=wk1(i,nym-j+2,k)*xb/2.+wk2(i,nym-j+2,k)*xa/2.           
         us(i,j,k)=xx3+xx4            
         tr(i,j,k)=xx1+xx2      
      enddo
      do i=2,nxm/2
         tr(i+nxm/2,j,k)=tr(nxm/2-i+2,j,k)
         us(i+nxm/2,j,k)=-us(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      do k=1,nzm
      do i=1,nxm
         tr(i,1,k)=tr(i,1,k)*2.
         us(i,1,k)=us(i,1,k)*2.
      enddo
      enddo
!********************fin de la tranformation en y***********************************
      ps(:,:,:)=0.
      do k=1,nz
      do j=1,ny
      do i=1,nx
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         xx1=(tr(i,j,k)+tr(nxm-i+2,j,k))*xa/2
         xx2=(us(i,j,k)-us(nxm-i+2,j,k))*xb/2.
         ps(i,j,k)=xx1+xx2
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         ps(1,j,k)=2.*ps(1,j,k)
         ps(nxm+1,j,k)=0.
      enddo
      enddo
      ps(mx,:,:)=0.; ps(mx,:,:)=0.; ps(:,my,:)=0.; ps(:,my-1,:)=0.; ps(:,:,mz)=0.; ps(:,:,mz-1)=0.
      tab1(:,:,:)=0. ; tab1(:,:,:)=ps(:,:,:)
   endif
   if (isign==-1) then
      us(:,:,:)=tab1(:,:,:)
      wk1(:,:,:)=0. ; wk2(:,:,:)=0.
      do k=1,nz
      do j=1,ny
      do i=1,nxm/2+1
         wk1(i,j,k)=us(i,j,k)
         wk2(i,j,k)=us(nxm-i+2,j,k)
      enddo
      enddo
      enddo
      tr(:,:,:)=0. ; ps(:,:,:)=0.
!*****************transformation en x**************************
     do k=1,nz
     do j=1,ny
     do i=1,nxm/2+1
        xb=sin((i-1)/2.*pi/(nxm))
        xa=cos((i-1)/2.*pi/(nxm))
        xb1=sin((nxm-i+1)/2.*pi/(nxm))
        xa1=cos((nxm-i+1)/2.*pi/(nxm))
       tr(i,j,k)=wk1(i,j,k)*xa+wk2(i,j,k)*xa1
        ps(i,j,k)=(wk1(i,j,k)*xb-wk2(i,j,k)*xb1)
     enddo
     do i=2,nxm/2
        tr(i+nxm/2,j,k)=tr(nxm/2-i+2,j,k)
        ps(i+nxm/2,j,k)=-ps(nxm/2-i+2,j,k)
     enddo
     enddo
     enddo
     do j=1,my
     do i=1,mx
        tr(i,j,nz)=0.
        ps(i,j,nz)=0.
     enddo
     enddo
     do k=1,mz
     do i=1,mx
       tr(i,ny,k)=0.
       ps(i,ny,k)=0.
     enddo
     enddo
!***************fin de la transformation en x************************
     wk1(:,:,:)=0. ; wk2(:,:,:)=0.
!***************debut transformation en y****************************
     do k=1,nz
     do j=1,ny
     do i=1,nxm/2+1
        xa=cos((j-1)/2.*pi/(nym))
        xb=sin((j-1)/2.*pi/(nym))
        xa1=cos((nym-j+1)/2.*pi/(nym))
        xb1=sin((nym-j+1)/2.*pi/(nym))
        xx1=tr(i,j,k)*xa-ps(i,j,k)*xa1
        xx2=tr(i,nym-j+2,k)*xb+ps(i,nym-j+2,k)*xa 
        xx3=-tr(i,nym-j+2,k)*xb1+ps(i,nym-j+2,k)*xb
        xx4=tr(i,j,k)*xb+ps(i,j,k)*xa
        wk1(i,j,k)=xx1+xx2
        wk2(i,j,k)=(xx3+xx4)
     enddo
     do i=2,nxm/2
        wk1(i+nxm/2,j,k)=wk1(nxm/2-i+2,j,k)
        wk2(i+nxm/2,j,k)=-wk2(nxm/2-i+2,j,k)
     enddo
     enddo
     enddo    
     do k=1,mz
     do i=1,mx
       wk1(i,ny,k)=0.
       wk2(i,ny,k)=0.
     enddo
     enddo
!**********************mode kc/2*******************
     ps(:,:,:)=0.
     do j=1,my
        ps(nx,j,1)=wk1(nxm/2+1,j,1)
     enddo
!**************fin transformation en y*******************************
     tr(:,:,:)=0.
!***************debut transformation en z****************************
     do k=1,nz
     do j=1,ny
     do i=1,nxm/2+1
        xa=cos((k-1)/2.*pi/(nzm))
        xb=sin((k-1)/2.*pi/(nzm))
        xa1=cos((nzm-k+1)/2.*pi/(nzm))
        xb1=sin((nzm-k+1)/2.*pi/(nzm))
        xx1=wk1(i,j,k)*xa-wk2(i,j,k)*xa1
        xx2=wk1(i,j,nzm-k+2)*xb+wk2(i,j,nzm-k+2)*xa 
        xx3=-wk1(i,j,nzm-k+2)*xb1+wk2(i,j,nzm-k+2)*xb
        xx4=wk1(i,j,k)*xb+wk2(i,j,k)*xa
        tr(2*i-1,j,k)=xx1+xx2
        tr(2*i,j,k)=-(xx3+xx4)
     enddo
     enddo
     enddo  
     do j=1,my
     do i=1,mx
       tr(i,j,nz)=0.
     enddo
     enddo
     do j=1,my
     do i=1,mx
        tr(i,j,nz)=tr(i,j,(nz+1)/2)
     enddo
     enddo
!***************fin transformation en z******************************* 
     do j=1,my
        tr(nx,j,1)=ps(nx,j,1)
     enddo

     call csfft3d(-1,nxm,nym,nzm,1.,tr,mx/2,my,tr,mx,my,table,work,0)

     wk1(:,:,:)=0. ; wk2(:,:,:)=0.
     do k=1,nz
     do j=1,ny
     do i=1,nxm/2
        wk1(2*i-1,j,k)=tr(i,j,k)
        wk1(2*i,j,k)=tr(nxm-i+1,j,k) 
     enddo
     enddo
     enddo
     do k=1,nzm
     do j=1,nym/2   
     do i=1,nxm
        wk2(i,2*j-1,k)=wk1(i,j,k)
        wk2(i,2*j,k)=wk1(i,nym-j+1,k)
     enddo
     enddo
     enddo
     do k=1,nzm/2
     do j=1,nym   
     do i=1,nxm
        ppm(i,j,2*k-1)=wk2(i,j,k)
        ppm(i,j,2*k)=wk2(i,j,nzm-k+1)
     enddo
     enddo
     enddo
  endif
endif 
!********************************************************************************************************
!********************************************************************************************************
!********************************************************************************************************
!********************************************************************************************************
!********************************************************************************************************
if (((nclx.eq.2).or.(nclx.eq.1)).and.((ncly.eq.1).or.(ncly.eq.2)).and.(nclz.eq.0)) then
   if (isign==1) then
      tr(:,:,:)=0. ; us(:,:,:)=0. ; ps(:,:,:)=0.
      wk1(:,:,:)=0. ; wk2(:,:,:)=0.  
      do k=1,nzm
      do j=1,nym
      do i=1,nxm/2
         tr(i,j,k)=ppm(2*(i-1)+1,j,k)
      enddo
      do i=nxm/2+1,nxm
         tr(i,j,k)=ppm(2*nxm-2*i+2,j,k)
      enddo
      enddo
      enddo
      do k=1,nzm
      do j=1,nym/2
      do i=1,nxm
         ps(i,j,k)=tr(i,2*(j-1)+1,k)
      enddo
      enddo 
      do j=nym/2+1,nym
      do i=1,nxm
         ps(i,j,k)=tr(i,2*nym-2*j+2,k)
      enddo
      enddo
      enddo
   
      call scfft3d(0,nxm,nym,nzm,1./nxm/nym/nzm,ps,mx,my,ps,mx/2,my,table,work,0)
      call scfft3d(1,nxm,nym,nzm,1./nxm/nym/nzm,ps,mx,my,ps,mx/2,my,table,work,0)

!***************************************************

     do k=1,mz
     do j=1,my
     do i=1,nxm/2+1
        tr(i,j,k)=ps(2*i-1,j,k)
        us(i,j,k)=-ps(2*i,j,k)
     enddo
     do i=2,nxm/2
        tr(i+nxm/2,j,k)=tr(nxm/2-i+2,j,k)
        us(i+nxm/2,j,k)=-us(nxm/2-i+2,j,k)
     enddo
     enddo
     enddo
     do k=1,nzm
     do j=1,nym
        tr(nxm/2+1,j,k)=ps(nxm+1,j,k)
        us(nxm/2+1,j,k)=-ps(nx+1,j,k)
     enddo
     enddo

!**************transformation en y***************************
     do k=1,nz
     do j=1,ny    
     do i=1,nxm/2+1
        xa=cos((j-1)/2.*pi/(nym))
        xb=sin((j-1)/2.*pi/(nym))
        xa1=cos((nym-j+1)/2.*pi/(nym))
        xb1=sin((nym-j+1)/2.*pi/(nym))
        xx1=tr(i,j,k)*xa/2.+us(i,j,k)*xb/2.
        xx2=tr(i,nym-j+2,k)*xa/2.-us(i,nym-j+2,k)*xb/2.
        xx3=-tr(i,j,k)*xb/2.+us(i,j,k)*xa/2.
        xx4=tr(i,nym-j+2,k)*xb/2.+us(i,nym-j+2,k)*xa/2.           
        wk1(i,j,k)=xx1+xx2  
        wk2(i,j,k)=xx3+xx4 
     enddo
     do i=2,nxm/2 
        wk1(i+nxm/2,j,k)=wk1(nxm/2-i+2,j,k)
        wk2(i+nxm/2,j,k)=-wk2(nxm/2-i+2,j,k)
     enddo
     enddo
     enddo
     do k=1,mz
     do i=1,mx
        wk1(i,1,k)=wk1(i,1,k)*2.
        wk2(i,1,k)=wk2(i,1,k)*2.
        wk1(i,ny,k)=wk1(i,ny,k)*2.
        wk2(i,ny,k)=wk2(i,ny,k)*2. 
     enddo
     enddo
     tr(:,:,:)=0.;us(:,:,:)=0.;ps(:,:,:)=0.
!******************transformation en x************************
     do k=1,mz
     do j=1,my
     do i=1,nx+1
        xb=sin((i-1)/2.*pi/(nxm))
        xa=cos((i-1)/2.*pi/(nxm))
        xx1=(wk1(i,j,k)+wk1(nxm-i+2,j,k))*xa/2
        xx2=(wk2(i,j,k)-wk2(nxm-i+2,j,k))*xb/2.
        ps(i,j,k)=xx1+xx2
     enddo
     xb=sin((nxm/4.)*pi/(nxm))
     xa=cos((nxm/4.)*pi/(nxm))
     ps(nx+1,j,k)=-((wk2(nxm/2+1,j,k)*xa/2.)+(wk2(nxm/2+1,j,k)*xb/2.))
     ps(1,j,k)=2.*ps(1,j,k)
     ps(nx,j,k)=2.*ps(nx,j,k)
     enddo
     enddo
     do k=1,mz
     do i=1,mx
        ps(i,ny,k)=0.
     enddo
     enddo
     wk1(:,:,:)=0. ; wk2(:,:,:)=0. ; 
!***************tranformation en z******************************
     do i=1,mxy
        tr(i,1,nzm/2+1)=ps(i,1,nzm/2+1)
     enddo
     do k=1,nzm/2
     do j=1,nym
     do i=1,nxm/2
        xa=cos((k-1)*pi/(nzm))
        xb=sin((k-1)*pi/(nzm))
        xx1=ps(i,j,k)*xa/2.+ps(i,j,nzm-k+2)*xa/2.
        xx2=-ps(nxm-i+2,j,k)*xb/2.+ps(nxm-i+2,j,nzm-k+2)*xb/2. 
        xx3=ps(i,j,k)*xb/2.+ps(i,j,nzm-k+2)*xb/2.
        xx4=ps(nxm-i+2,j,k)*xa/2.-ps(nxm-i+2,j,nzm-k+2)*xa/2. 
        us(i,j,2*k-1)=xx1+xx2
        us(i,j,2*k)=-(xx3+xx4)
     enddo
     do i=nxm/2+1,nx+1
        xa=cos((nzm-k+1)*pi/(nzm))
        xb=sin((nzm-k+1)*pi/(nzm))
        xx1=ps(i,j,k)*xa/2.+ps(i,j,nzm-k+2)*xa/2.
        xx2=-ps(nxm-i+2,j,k)*xb/2.+ps(nxm-i+2,j,nzm-k+2)*xb/2.
        xx3=ps(i,j,k)*xb/2.+ps(i,j,nzm-k+2)*xb/2.
        xx4=ps(nxm-i+2,j,k)*xa/2.-ps(nxm-i+2,j,nzm-k+2)*xa/2.  
        us(i,j,2*k-1)=-(xx1+xx2)
        us(i,j,2*k)=-(xx3+xx4) 
     enddo
     xa=cos((nzm-k+1)*pi/(nzm))
     xb=sin((nzm-k+1)*pi/(nzm))
     xx1=ps(nxm/2+1,j,k)*xa/2.+ps(nxm/2+1,j,nzm-k+2)*xa/2.
     xx2=-ps(nx+1,j,k)*xb/2.+ps(nx+1,j,nzm-k+2)*xb/2.
     xx3=ps(nxm/2+1,j,k)*xb/2.+ps(nxm/2+1,j,nzm-k+2)*xb/2.
     xx4=ps(nx+1,j,k)*xa/2.-ps(nx+1,j,nzm-k+2)*xa/2.
     us(nxm/2+1,j,2*k-1)=-(xx1-xx2)
     us(nxm/2+1,j,2*k)=-(xx3-xx4)
     enddo
     enddo
     do i=1,mxyz
        ps(i,1,1)=0.
        ps(i,1,1)=us(i,1,1)
     enddo
     do i=1,mxy
        ps(i,1,nz+1)=tr(i,1,nzm/2+1)
        ps(i,1,2)=0.
     enddo
     ps(:,2:ny,1)=2.*ps(:,2:ny,1)
     ps(2:nx,1,1)=2.*ps(2:nx,1,1)
     do k=1,mz
     do j=1,my
        ps(nx+1,j,k)=0.
     enddo
     enddo
     ps(mx,:,:)=0.; ps(mx-1,:,:)=0.; ps(:,my,:)=0.; ps(:,my-1,:)=0.; ps(:,:,mz)=0.
     tab1(:,:,:)=0.;tab1(:,:,:)=ps(:,:,:)
  endif
  if (isign==-1) then
     ps(:,:,:)=0. ; ps(:,:,:)=tab1(:,:,:)
     tr(:,:,:)=0. ; us(:,:,:)=0. ; wk1(:,:,:)=0.
     do k=1,nz+1
     do j=1,ny
     do i=1,nxm/2+1
        tr(i,j,k)=ps(i,j,k)
        us(i,j,k)=-ps(nxm-i+2,j,k)
     enddo
     enddo
     enddo
     do i=1,mxyz
        ps(i,1,1)=0.
        tab1(i,1,1)=0.
     enddo
!*****************transformation en x**************************
     do k=1,nz+1,2
     do j=1,ny         
     do i=1,nxm/2+1
        xb=sin((i-1)/2.*pi/(nxm))
        xa=cos((i-1)/2.*pi/(nxm))
        xb1=sin((nxm-i+1)/2.*pi/(nxm))
        xa1=cos((nxm-i+1)/2.*pi/(nxm))
        wk1(i,j,k)=-tr(i,j,k+1)*xb-us(i,j,k+1)*xa
        wk2(i,j,k)=-(tr(i,j,k+1)*xa-us(i,j,k+1)*xb)
        ps(i,j,k)=tr(i,j,k)*xa-us(i,j,k)*xb
        tab1(i,j,k)=tr(i,j,k)*xb+us(i,j,k)*xa
     enddo
     enddo
     enddo
     do i=2,nxm/2
     do j=1,myz
        wk1(i+nxm/2,j,1)=wk1(nxm/2-i+2,j,1)
        ps(i+nxm/2,j,1)=wk1(nxm/2-i+2,j,1)
        wk2(i+nxm/2,j,1)=-wk2(nxm/2-i+2,j,1)
        tab1(i+nxm/2,j,1)=-wk2(nxm/2-i+2,j,1)
     enddo
     enddo
!****************************************************************************
     do i=1,mxy
        wk1(i,1,nz)=0.
        wk2(i,1,nz)=0
        ps(i,1,nz)=0.
        tab1(i,1,nz)=0.
     enddo
     do k=1,mz
     do i=1,mx
        wk1(i,ny,k)=0. 
        wk2(i,ny,k)=0. 
        ps(i,ny,k)=0. 
        tab1(i,ny,k)=0.
     enddo
     enddo
     do i=1,mxyz
        tr(i,1,1)=0.
        us(i,1,1)=0.
     enddo
!***************************************************************************
     do i=1,mxy
        tr(i,1,nz+1)=wk1(i,1,nz+1)+ps(i,1,nz+1)
        us(i,1,nz+1)=wk2(i,1,nz+1)+tab1(i,1,nz+1)
     enddo
!********************************************************************
     tr(:,:,:)=0.; tr(:,:,:)=tab1(:,:,:) 
     us(:,:,:)=0. ; tab1(:,:,:)=0.
     xtoto=ps(nxm/2+1,1,1)
!***************debut transformation en z****************************
     do k=1,nz+1
     do j=1,ny
     do i=1,nxm/2+1
        xa=cos((k-1)/2.*pi/(nzm))
        xb=sin((k-1)/2.*pi/(nzm))
        xa1=cos((nzm-k+1)/2.*pi/(nzm))
        xb1=sin((nzm-k+1)/2.*pi/(nzm))
        xx1=wk1(i,j,k)*xb-wk2(i,j,k)*xa
        xx2=wk1(i,j,nzm-k+2)*xa+wk2(i,j,nzm-k+2)*xb 
        xx3=-wk1(i,j,nzm-k+2)*xa1+wk2(i,j,nzm-k+2)*xb1
        xx4=wk1(i,j,k)*xa+wk2(i,j,k)*xb
        xx5=ps(i,j,k)*xa-tr(i,j,k)*xa1
        xx6=ps(i,j,nzm-k+2)*xb+tr(i,j,nzm-k+2)*xa 
        xx7=-ps(i,j,nzm-k+2)*xb1+tr(i,j,nzm-k+2)*xb
        xx8=ps(i,j,k)*xb+tr(i,j,k)*xa
        tab1(2*i-1,j,k)=xx8+xx1
        tab1(2*i,j,k)=xx5+xx4
        us(2*i-1,j,k)=xx7+xx2
        us(2*i,j,k)=xx6+xx3
     enddo
     enddo
     enddo
     wk1(:,:,:)=0. ; wk2(:,:,:)=0.      
     do j=2,my
     do i=1,mx
        tab1(i,j,nz)=tab1(i,j,(nz+1)/2)
     enddo
     enddo
     do k=1,nz,2
     do j=1,my
     do i=1,mx
        wk1(i,j,(k+1)/2)=tab1(i,j,k)
        wk1(i,j,nz/2+(k+1)/2)=us(i,j,k)
     enddo
     enddo
     enddo
     do k=1,mz
     do j=nym/2+1,my
     do i=nxm/2+1,mx
        wk1(i,j,k)=-wk1(i,j,k)
     enddo
     enddo
     enddo
     do k=1,nz
     do j=1,nym/2+1
     do i=1,nxm/2,2
        xa=cos((j-1)/2.*pi/(nym))
        xb=sin((j-1)/2.*pi/(nym))
        xx1=-(wk1(i,j,k)*xb-wk1(i+1,j,k)*xa)+(wk1(i+1,nym-j+2,k)*xb+wk1(i,nym-j+2,k)*xa)
        xx2=-(wk1(i+1,j,k)*xb+wk1(i,j,k)*xa)-(wk1(i,nym-j+2,k)*xb-wk1(i+1,nym-j+2,k)*xa)
        xx3=(wk1(i,j,k)*xb+wk1(i+1,j,k)*xa)+(wk1(i+1,nym-j+2,k)*xb-wk1(i,nym-j+2,k)*xa)
        xx4=(wk1(i+1,j,k)*xb-wk1(i,j,k)*xa)-(wk1(i,nym-j+2,k)*xb+wk1(i+1,nym-j+2,k)*xa)
        ps(i,nym-j+2,k)=xx3
        ps(i+1,nym-j+2,k)=xx4
        ps(i,j,k)=xx1
        ps(i+1,j,k)=xx2
     enddo
     do i=nxm/2+1,nx+1,2
        xa=cos((j-1)/2.*pi/(nym))
        xb=sin((j-1)/2.*pi/(nym))
        xx1=-(wk1(i,j,k)*xb-wk1(i+1,j,k)*xa)-(wk1(i+1,nym-j+2,k)*xb+wk1(i,nym-j+2,k)*xa)
        xx2=-(wk1(i+1,j,k)*xb+wk1(i,j,k)*xa)+(wk1(i,nym-j+2,k)*xb-wk1(i+1,nym-j+2,k)*xa)
        xx3=(wk1(i,j,k)*xb+wk1(i+1,j,k)*xa)-(wk1(i+1,nym-j+2,k)*xb-wk1(i,nym-j+2,k)*xa)
        xx4=(wk1(i+1,j,k)*xb-wk1(i,j,k)*xa)+(wk1(i,nym-j+2,k)*xb+wk1(i+1,nym-j+2,k)*xa)
        ps(i,nym-j+2,k)=xx3
        ps(i+1,nym-j+2,k)=xx4
        ps(i,j,k)=xx1
        ps(i+1,j,k)=xx2
     enddo
     enddo
     enddo
     do i=1,mxyz
        wk1(i,1,1)=0.
     enddo
     do i=1,mxyz
        wk1(i,1,1)=ps(i,1,1)
     enddo
     ps(1,:,nzm/2+1)=wk1(2,:,nzm/2+1)
     ps(2,:,nzm/2+1)=-wk1(1,:,nzm/2+1)
     do k=1,mz
     do i=nxm/2+1,nx,2
        ps(i,nym/2+1,k)=wk1(i+1,nym/2+1,k)
        ps(i+1,nym/2+1,k)=-wk1(i,nym/2+1,k)
     enddo
     enddo
     do j=1,my
     do i=3,nx,2
        ps(i,j,nzm/2+1)=wk1(i+1,j,nzm/2+1)
        ps(i+1,j,nzm/2+1)=-wk1(i,j,nzm/2+1) 
     enddo
     enddo
     do i=1,mxyz
        wk1(i,1,1)=ps(i,1,1)
     enddo
     do i=nxm/2+1,nx-1,2
        ps(i,nym/2+1,nzm/2+1)=wk1(i+1,nym/2+1,nzm/2+1)
        ps(i+1,nym/2+1,nzm/2+1)=-wk1(i,nym/2+1,nzm/2+1)
     enddo
     ps(nx,nym/2+1,nzm/2+1)=ps(nx+1,nym/2+1,nzm/2+1)
     ps(nx+1,nym/2+1,nzm/2+1)=0.
     ps(nx,1,1)=xtoto
     do i=1,mxy
        ps(i,1,nz+1)=0.
     enddo
     ps(2,1,nzm/2+1)=0.

     call csfft3d(-1,nxm,nym,nzm,1.,ps,mx/2,my,ps,mx,my,table,work,0)
   
     do i=1,mxyz
        wk1(i,1,1)=0.
     enddo
      
     do k=1,nz
     do j=1,ny
     do i=1,nxm/2
        wk1(2*i-1,j,k)=ps(i,j,k)
        wk1(2*i,j,k)=ps(nxm-i+1,j,k) 
     enddo
     enddo
     enddo    
     do k=1,nzm
     do j=1,nym/2   
     do i=1,nxm
        ppm(i,2*j-1,k)=wk1(i,j,k)
        ppm(i,2*j,k)=wk1(i,nym-j+1,k)
     enddo
     enddo
     enddo
  endif
endif
!********************************************************************************************************
!********************************************************************************************************
!********************************************************************************************************
!********************************************************************************************************
!********************************************************************************************************
if ((ncly.eq.0).and.(nclx.eq.0).and.(nclz.eq.0)) then
   if (isign==1) then
      tr(:,:,:)=0.;us(:,:,:)=0.;ps(:,:,:)=0.;wk1(:,:,:)=0.;wk2(:,:,:)=0.
      
      do k=1,nzm
      do j=1,nym
      do i=1,nxm
         tr(i,j,k)=ppm(i,j,k)
      enddo
      enddo
      enddo
      
      call scfft3d(0,nxm,nym,nzm,1./nxm/nym/nzm,tr,mx,my,tr,mx/2,my,table,work,0)
      call scfft3d(1,nxm,nym,nzm,1./nxm/nym/nzm,tr,mx,my,tr,mx/2,my,table,work,0)

      do k=1,mz
      do j=1,nym
      do i=1,nxm/2+1
         us(i,j,k)=tr(2*i-1,j,k)
         ps(i,j,k)=-tr(2*i,j,k)
      enddo
      do i=2,nxm/2
         us(i+nxm/2,j,k)=us(nxm/2-i+2,j,k)
         ps(i+nxm/2,j,k)=ps(nxm/2-i+2,j,k)
      enddo
      enddo
      do j=1,my
         us(mx,j,k)=us(nxm/2+1,j,k)
         ps(mx,j,k)=-ps(nxm/2+1,j,k)
      enddo
      enddo

      do k=1,mz
      do j=1,nym/2+1
      do i=1,nxm/2+1
         xa=cos((j-1)*pi/(nym))
         xb=sin((j-1)*pi/(nym))
         xx1=us(i,j,k)*xa/2.+ps(i,j,k)*xb/2.
         xx2=us(i,nym-j+2,k)*xa/2.-ps(i,nym-j+2,k)*xb/2.    
         xx3=-us(i,j,k)*xb/2.+ps(i,j,k)*xa/2.
         xx4=-(us(i,nym-j+2,k)*xb/2.+ps(i,nym-j+2,k)*xa/2.)          
         wk1(i,2*j-1,k)=xx1+xx2
         wk1(i,2*j,k)=xx3+xx4
      enddo
      enddo
      do j=1,nym/2
         xa=cos((j-1)*pi/(nym))
         xb=sin((j-1)*pi/(nym))
         xx1=us(nxm/2+1,j,k)*xa/2.-ps(nxm/2+1,nym-j+2,k)*xb/2.
         xx2=us(nxm/2+1,nym-j+2,k)*xa/2.+ps(mx,nym-j+2,k)*xb/2.
         xx3=us(nxm/2+1,j,k)*xb/2.+ps(nxm/2+1,nym-j+2,k)*xa/2.
         xx4=-us(nxm/2+1,nym-j+2,k)*xb/2.+ps(mx,nym-j+2,k)*xa/2.             
         wk1(nxm/2+1,2*j-1,k)=xx1+xx2
         wk1(nxm/2+1,2*j,k)=-xx3+xx4
      enddo
      do j=1,my
      do i=2,nxm/2+2
         wk1(i+nxm/2,j,k)=wk1(nxm/2-i+2,j,k)
      enddo
      enddo
      wk1(nxm/2+1,1,k)= wk1(nxm/2+1,1,k)*0.5
      enddo
      do k=1,mz
      do j=1,nym/2+1
      do i=1,nxm/2+1
         xa=cos((j-1)*pi/(nym))
         xb=sin((j-1)*pi/(nym))
         xa1=cos((nym-j+1)*pi/(nym))
         xb1=sin((nym-j+1)*pi/(nym))
         xx1=-us(i,j,k)*xb/2.+ps(i,j,k)*xa/2.
         xx2=us(i,nym-j+2,k)*xb/2.+ps(i,nym-j+2,k)*xa/2.     
         xx3=us(i,j,k)*xa/2.+ps(i,j,k)*xb/2.
         xx4=-us(i,nym-j+2,k)*xa/2.+ps(i,nym-j+2,k)*xb/2.       
         wk2(i,2*j-1,k)=xx1+xx2
         wk2(i,2*j,k)=-(xx3+xx4)
      enddo
      enddo
      do j=1,my
      do i=2,nxm/2+2
         wk2(i+nxm/2,j,k)=-wk2(nxm/2-i+2,j,k)
      enddo
      enddo
      enddo
      wk2(nxm/2+1,2,:)=0.                          !!!!
      wk2(nxm/2+1,1,:)=wk2(nxm/2+1,1,:)*0.5       !!!!

      do k=1,mz
      do i=1,mx
         wk1(i,1,k)=wk1(i,1,k)*2.
         wk2(i,1,k)=wk2(i,1,k)*2.
      enddo
      enddo

    us(:,:,:)=0.       
      do k=1,mz
      do j=1,my
      do i=1,nxm/2+2
         xa=cos((i-1)*pi/(nxm))
         xb=sin((i-1)*pi/(nxm))
         xx1=wk1(i,j,k)*xa/2.+wk1(nxm-i+2,j,k)*xa/2.
         xx2=-wk2(nxm-i+2,j,k)*xb/2+wk2(i,j,k)*xb/2. 
         xx3=wk1(i,j,k)*xb/2.+wk1(nxm-i+2,j,k)*xb/2.
         xx4=wk2(nxm-i+2,j,k)*xa/2.+wk2(nxm-i+2,j,k)*xa/2. 
         us(2*i-1,j,k)=xx1+xx2
         us(2*i,j,k)=-(xx3+xx4)
      enddo
      enddo
      enddo
      do k=1,mz
      do j=1,my
         xb=sin(pi/2.)
         xx2=-wk2(nxm/2+1,j,k)*xb/2-wk2(nxm/2+1,j,k)*xb/2. 
         xx3=wk1(nxm/2+1,j,k)*xb/2.+wk1(nxm/2+1,j,k)*xb/2.
         us(mx-1,j,k)=-xx2
         us(mx,j,k)=-xx3
      enddo
      enddo
    
    do k=1,mz          
      do j=1,my
         us(1,j,k)=us(1,j,k)
      enddo
      us(:,2,k)=0.
      us(mx-1,:,k)=us(mx-1,:,k)*0.5 !##############
      us(mx,:,k)=us(mx,:,k)*0.5     !##############
      enddo


      ps(:,:,:)=0.;ps(:,:,:)=us(:,:,:)
      do j=1,my
      do i=1,mx
         tr(i,j,nzm/2+1)=ps(i,j,nzm/2+1)
      enddo
      enddo
      tab1(:,:,:)=0.
      do j=1,my
      do k=1,nzm/2      
      do i=1,mx,2
         xa=cos((k-1)*pi/(nzm))
         xb=sin((k-1)*pi/(nzm))
         xx1=ps(i,j,k)*xa/2.+ps(i,j,nzm-k+2)*xa/2.
         xx2=-ps(i+1,j,k)*xb/2.+ps(i+1,j,nzm-k+2)*xb/2. 
         xx3=ps(i,j,k)*xb/2.+ps(i,j,nzm-k+2)*xb/2.
         xx4=-ps(i+1,j,k)*xa/2.+ps(i+1,j,nzm-k+2)*xa/2.
         xx5=ps(i,j,k)*xb/2.-ps(i,j,nzm-k+2)*xb/2.
         xx6=-ps(i+1,j,k)*xa/2.-ps(i+1,j,nzm-k+2)*xa/2. 
         xx7=ps(i,j,k)*xa/2.-ps(i,j,nzm-k+2)*xa/2.
         xx8=ps(i+1,j,k)*xb/2.+ps(i+1,j,nzm-k+2)*xb/2.
         tab1(i,j,2*k-1)=xx1-xx2
         tab1(i,j,2*k)=-(xx3+xx4)
         tab1(i+1,j,2*k-1)=-(xx5+xx6)
         tab1(i+1,j,2*k)=-(xx7+xx8)
      enddo
      enddo
      enddo

      tab1(:,ny+1,1)=ps(:,ny+1,1)
      ps(:,:,:)=0. ; ps(:,:,:)=tab1(:,:,:)
      ps(:,:,nz+1)=tr(:,:,nzm/2+1)
      ps(:,:,2)=0.
      ps(:,2:my,1)=2.*ps(:,2:my,1)
      ps(2:mx,1,1)=2.*ps(2:mx,1,1)
      ps(nx+1,:,:)=0.
      ps(:,2,:)=0.
      ps(1,ny+1,2:mz)= ps(1,ny+1,2:mz)/2.
      ps(:,:,mz)=-ps(:,:,mz-1)
      ps(:,:,mz-1)=0.      
      tab1(:,:,:)=0. ; tab1(:,:,:)=ps(:,:,:)
   endif
   if  (isign==-1) then
      tr(:,:,:)=0.;tr(:,:,:)=tab1(:,:,:)
      do k=1,mz
      do j=1,my
         tr(nxm+1,j,k)=-tr(mx,j,k)
      enddo
      do j=2,my
         tr(nxm/2+1,j,k)=tr(nxm/2+1,j,k)/2.
         tr(nxm/2+2,j,k)=tr(nxm/2+2,j,k)/2.
      enddo
      enddo
      wk1(:,:,:)=0. ;wk2(:,:,:)=0.
      do k=1,mz
      do j=1,my
      do i=1,nxm/2+1,2
         wk1(i,j,k)=tr(i,j,k)
         wk1(i+1,j,k)=-tr(i+1,j,k)
         wk2(i,j,k)=tr(nxm-i+2,j,k)
         wk2(i+1,j,k)=-tr(nxm-i+1,j,k)
      enddo
      enddo
      do j=1,my
         wk2(nxm/2+2,j,k)=0.
      enddo
      enddo
      us(:,:,:)=0.;ps(:,:,:)=0.;tab1(:,:,:)=0.;tr(:,:,:)=0.
      do k=1,mz
      do j=1,my
      do i=1,mx-1,2
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         us((i+1)/2,j,k)=wk1(i+1,j,k)*xb+wk2(nxm-i+1,j,k)*xb
         ps((i+1)/2,j,k)=-wk1(i+1,j,k)*xa-wk2(nxm-i+1,j,k)*xa
         tab1((i+1)/2,j,k)=wk1(i,j,k)*xa+wk2(nxm-i+2,j,k)*xa
         tr((i+1)/2,j,k)=wk1(i,j,k)*xb+wk2(nxm-i+2,j,k)*xb
      enddo
      enddo
      enddo
      do k=1,mz
      do j=1,my
      do i=2,nxm/2
         us(i+nxm/2,j,k)=us(nxm/2-i+2,j,k)
         tab1(i+nxm/2,j,k)=tab1(nxm/2-i+2,j,k)
         ps(i+nxm/2,j,k)=-ps(nxm/2-i+2,j,k)
         tr(i+nxm/2,j,k)=-tr(nxm/2-i+2,j,k)
      enddo
      enddo
      wk1(:,:,:)=0.
      do j=1,my
         wk1(nxm/2+1,j,k)=tr(nxm/2+1,j,k)
         tab1(nxm/2+1,j,k)=wk1(nxm/2+1,j,k)
         tr(nxm/2+1,j,k)=0.
      enddo
      enddo
      wk2(:,:,:)=0.
      do k=1,mz
      do j=1,ny+1,2
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/nym)
         xb=sin((j-1)/2.*pi/nym)
         xx1=(tab1(i,j,k)+tab1(nxm-i+2,j,k))*xa/2.-(tr(i,j,k)-tr(nxm-i+2,j,k))*xb/2.
         xx2=(tab1(i,j,k)+tab1(nxm-i+2,j,k))*xa/2.+(tr(i,j,k)-tr(nxm-i+2,j,k))*xb/2.
         xx3=(tab1(i,j+1,k)+tab1(nxm-i+2,j+1,k))*xa/2.-(tr(i,j+1,k)-tr(nxm-i+2,j+1,k))*xb/2.
         xx4=-((tab1(i,j+1,k)+tab1(nxm-i+2,j+1,k))*xa/2.+(tr(i,j+1,k)-tr(nxm-i+2,j+1,k))*xb/2.)
         xx5=-((tab1(i,j,k)+tab1(nxm-i+2,j,k))*xb/2.-(tr(i,j,k)-tr(nxm-i+2,j,k))*xa/2.)
         xx6=(tab1(i,j,k)+tab1(nxm-i+2,j,k))*xb/2.+(tr(i,j,k)-tr(nxm-i+2,j,k))*xa/2.
         xx7=-((tab1(i,j+1,k)+tab1(nxm-i+2,j+1,k))*xb/2.-(tr(i,j+1,k)-tr(nxm-i+2,j+1,k))*xa/2.)
         xx8=-((tab1(i,j+1,k)+tab1(nxm-i+2,j+1,k))*xb/2.+(tr(i,j+1,k)-tr(nxm-i+2,j+1,k))*xa/2.)
         yy1=(us(i,j,k)+us(nxm-i+2,j,k))*xa/2.-(ps(i,j,k)-ps(nxm-i+2,j,k))*xb/2.
         yy2=(us(i,j,k)+us(nxm-i+2,j,k))*xa/2.+(ps(i,j,k)-ps(nxm-i+2,j,k))*xb/2.
         yy3=-((us(i,j+1,k)+us(nxm-i+2,j+1,k))*xa/2.+(ps(i,j+1,k)-ps(nxm-i+2,j+1,k))*xb/2.)
         yy4=(us(i,j,k)+us(nxm-i+2,j,k))*xb/2.+(ps(i,j,k)-ps(nxm-i+2,j,k))*xa/2.
         yy5=-((us(i,j+1,k)+us(nxm-i+2,j+1,k))*xb/2.-(ps(i,j+1,k)-ps(nxm-i+2,j+1,k))*xa/2.)
         yy6=-((us(i,j+1,k)+us(nxm-i+2,j+1,k))*xb/2.+(ps(i,j+1,k)-ps(nxm-i+2,j+1,k))*xa/2.)
         yy7=(us(i,j+1,k)+us(nxm-i+2,j+1,k))*xa/2.-(ps(i,j+1,k)-ps(nxm-i+2,j+1,k))*xb/2.
         yy8=-((us(i,j,k)+us(nxm-i+2,j,k))*xb/2.-(ps(i,j,k)-ps(nxm-i+2,j,k))*xa/2.)
         wk2(2*i-1,(j+1)/2,k)=xx1+xx8+yy1+yy6
         wk2(2*i,(j+1)/2,k)=-xx6-xx3-yy4-yy7
         wk2(2*i-1,nym-(j+1)/2+2,k)=xx2+xx7+yy2+yy5
         wk2(2*i,nym-(j+1)/2+2,k)=-xx4-xx5-yy3-yy8
      enddo
      enddo
      enddo
      do k=1,mz
         wk2(nxm/2+1,1,k)=wk2(nxm/2+1,1,k)/2.
         wk2(nxm/2+2,1,k)=wk2(nxm/2+2,1,k)/2.
         wk2(mx,1,k)=0.
         wk2(nx+1,1,k)=wk2(nx+1,1,k)*2.
         wk2(:,ny+1,k)=0.
         wk2(1:2,:,k)=wk2(1:2,:,k)*2.
         wk2(nx+1:nx+2,:,k)=wk2(nx+1:nx+2,:,k)*2.
      enddo
!************TRANSFO EN Z*************************************************
      do k=1,nz+1,2
      do j=1,ny
      do i=1,nx+1,2
         xa=cos((k-1)/2.*pi/(nzm))
         xb=sin((k-1)/2.*pi/(nzm))
         xx1=(wk2(i,j,k)*xa+wk2(i+1,j,k)*xb)+(-wk2(i,j,k+1)*xb+wk2(i+1,j,k+1)*xa)
         xx2=(-wk2(i,j,k)*xb+wk2(i+1,j,k)*xa)-(wk2(i,j,k+1)*xa+wk2(i+1,j,k+1)*xb) 
         xx3=(wk2(i,j,k)*xa-wk2(i+1,j,k)*xb)-(wk2(i,j,k+1)*xb+wk2(i+1,j,k+1)*xa) 
         xx4=(wk2(i,j,k)*xb+wk2(i+1,j,k)*xa)+(wk2(i,j,k+1)*xa-wk2(i+1,j,k+1)*xb) 
         tr(i,j,(k+1)/2)=xx1
         tr(i+1,j,(k+1)/2)=xx2
         tr(i,j,nzm-(k+1)/2+2)=xx3
         tr(i+1,j,nzm-(k+1)/2+2)=xx4
      enddo
      enddo
      enddo
      tr(:,:,nz+1)=0.
      call csfft3d(-1,nxm,nym,nzm,1.,tr,mx/2,my,tr,mx,my,table,work,0)
      do k=1,nzm
      do j=1,nym
      do i=1,nxm
         ppm(i,j,k)=tr(i,j,k)
      enddo
      enddo
      enddo
   endif
endif
!********************************************************************************************************
return
end subroutine slfft3d_shift
