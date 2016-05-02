module poisson_m
  use slfft3d_shift_m
  use tools_m
  use slfft2d_shift_m
implicit none
contains
!****************************************************************
subroutine poisson(ppm,sy7,sy8,sy9,sy10,sy11,sy12,&
                   work,table,a,a3, e,c,a2,sr,a1,b1,d1,d)
!
!****************************************************************

USE param 
USE variables

implicit none

real(8),dimension(nxm,nym,nzm) :: ppm
real(8),dimension(mx,my,mz) :: sy7,sy8,sy9,sy10,sy11,sy12
real(8), dimension(nwork) :: work
real(8), dimension(100+2*(nxm+nym+nzm)) :: table
real(8),dimension(mx,mz,ny/2,5) :: a,a2
real(8),dimension(mx,mz,ny) :: d1,d,e,c
real(8),dimension(mx,mz,ny,5) :: a3
integer,dimension(2) :: ja,jb
real(8),dimension(mx,mz) :: sr
real(8),dimension(mx,mz) :: a1,b1

if (nz.gt.1) then
   if (istret.ne.0) call matrice_raffinement(a,a2,a3,d1,d)
   call poisson_3d(ppm,sy7,sy8,sy9,sy10,sy11,sy12,&
                   work,table,a,a3,e,c,a2,sr,a1,b1)
else    
   if (istret.ne.0) call matrice_raffinement_2d(a,a2,a3,d1,d)
   call  poisson_2d(ppm,sy7,sy8,sy9,sy10,sy11,sy12,&
                   work,table,a,a3,e,c,a2,sr)
endif

return
end subroutine poisson

!****************************************************************
!
subroutine poisson_3d(ppm,sy7,sy8,sy9,sy10,sy11,sy12,&
                      work,table,a,a3, e,c,a2,sr,a1,b1)
!
!****************************************************************

USE param 
USE derivX
USE derivY
USE derivZ
USE variables

implicit none

integer :: i,j,k
real(8),dimension(nxm,nym,nzm) :: ppm
real(8),dimension(mx,my,mz) :: sy7,sy8,sy9,sy10,sy11,sy12
real(8), dimension(nwork) :: work
real(8), dimension(100+2*(nxm+nym+nzm)) :: table
real(8) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2,xyzk
real(8) :: xtt1,ytt1,ztt1,zt1,zt2,xl2,x,y,z,r1
real(8),dimension(mx,mz,ny/2,5) :: a,a2
real(8),dimension(mx,mz,ny) :: e,c
real(8),dimension(mx,mz,ny,5) :: a3
integer,dimension(2) :: ja,jb
real(8),dimension(mx,mz) :: sr
real(8),dimension(mx,mz) :: a1,b1
real(8),external :: rand2
 
call slfft3d_shift(ppm,sy12,sy7,sy8,sy9,sy10,sy11,work,nx,nxm,ny,nym,nz,nzm,mx,my,mz,nwork,table,1)

if (istret.ne.0) then
   if (istret.ne.3) then
      e(:,:,:)=0. ; c(:,:,:)=0.
      do k=1,mz
      do i=1,mx
      do j=1,ny/2
         e(i,k,j)=sy12(i,2*j-1,k)
         c(i,k,j)=sy12(i,2*j,k)
      enddo
      enddo
      enddo
 
      call inversion5(a,ja,jb,sr,a1,b1,mx,ny,mz,e)       
      call inversion5(a2,ja,jb,sr,a1,b1,mx,ny,mz,c)
 
      sy12(:,:,:)=0.
      do k=1,mz
      do j=1,ny-1,2
      do i=1,mx
         sy12(i,j,k)=e(i,k,(j+1)/2)
      enddo
      enddo
      do j=2,ny,2
      do i=1,mx
         sy12(i,j,k)=c(i,k,j/2)
      enddo
      enddo
      enddo
!      do k=1,mz
!      do i=1,mx
!         if ((xkx(i)==0).and.(zkz(k)==0)) then
!            sy12(i,1,1)=0.
!            sy12(i,ny,1)=0.
!         endif
!      enddo
!      enddo
      sy12(mx,:,mz)=0.


   endif
   if (istret.eq.3) then
      e(:,:,:)=0. ;
      do k=1,mz
      do i=1,mx
      do j=1,ny
         e(i,k,j)=sy12(i,j,k)
      enddo
      enddo
      enddo
      call inversion5_i3(a3,ja,jb,sr,a1,b1,mx,ny,mz,e)   
      sy12(:,:,:)=0.
      do k=1,mz
      do j=1,ny
      do i=1,mx
         sy12(i,j,k)=e(i,k,j)
      enddo
      enddo
      enddo
   endif
else
   do k=1,mz
   do j=1,my
   do i=1,mx
      xtt=(bicix6*2.*cos(exs(i)*3.*dx/2.)+cicix6*2.*cos(exs(i)*5.*dx/2.))
      ytt=(biciy6*2.*cos(eys(j)*3.*dy/2.)+ciciy6*2.*cos(eys(j)*5.*dy/2.))
      ztt=(biciz6*2.*cos(ezs(k)*3.*dz/2.)+ciciz6*2.*cos(ezs(k)*5.*dz/2.))
      xtt1=(aicix6*2.*cos(exs(i)*dx/2.))
      ytt1=(aiciy6*2.*cos(eys(j)*dy/2.))
      ztt1=(aiciz6*2.*cos(ezs(k)*dz/2.))
      xt1=(1.+2.*ailcaix6*cos(exs(i)*dx))
      yt1=(1.+2.*ailcaiy6*cos(eys(j)*dy))
      zt1=(1.+2.*ailcaiz6*cos(ezs(k)*dz))
      xt2=xk2(i)*((((ytt1+ytt)/yt1)*((ztt1+ztt)/zt1))**2)
      yt2=yk2(j)*((((xtt1+xtt)/xt1)*((ztt1+ztt)/zt1))**2)
      zt2=zk2(k)*((((xtt1+xtt)/xt1)*((ytt1+ytt)/yt1))**2)
      xyzk=xt2+yt2+zt2
      if (xyzk.lt.1.e-8) then
         sy12(i,j,k)=0.
      else
         sy12(i,j,k)=sy12(i,j,k)/(-xyzk)
      endif
   enddo
   enddo
   enddo



endif

  

call slfft3d_shift(ppm,sy12,sy7,sy8,sy9,sy10,sy11,work,nx,nxm,ny,nym,nz,nzm,mx,my,mz,nwork,table,-1)
return
end subroutine poisson_3d

!********************************************************************
subroutine poisson_2d(ppm,tb6,tb1,tb2,tb3,tb4,tb5,&
                      work,table,a_2d,a3_2d,e,c,a2_2d,sr)
!
!****************************************************************

USE param
USE derivX
USE derivY
USE derivZ
USE variables

implicit none

integer :: i,j
real(8),dimension(nxm,nym,nzm) :: ppm
real(8),dimension(nwork) :: work
real(8),dimension(100+2*(nxm+nym+nzm)) :: table ! TODO
real(8),dimension(mx,my) :: tb1,tb2,tb3,tb4,tb5,tb6 
real(8),dimension(mx,ny/2,5) :: a_2d,a2_2d
real(8),dimension(mx,ny,5) :: a3_2d
real(8),dimension(mx) :: b1,a1
real(8),dimension(mx,ny) :: e,c
real(8),dimension(mx) :: sr
integer,dimension(2) :: ja,jb
real(8) :: ytt,xtt,yt,xt,yt1,xt1,xyzk,xt2,yt2

call slfft2d_shift(ppm,tb6,tb1,tb2,tb3,tb4,tb5,work,table,nx,nxm,ny,nym,nz,mx,my,mz,nwork,ntrigsX,ntrigsY,1)
if (istret.ne.0) then
   if ((istret==1).or.(istret==2)) then
      e(:,:)=0. ; c(:,:)=0.        
      do i=1,mx
      do j=1,ny/2
         e(i,j)=tb6(i,2*j-1)
         c(i,j)=tb6(i,2*j)
      enddo
      enddo
      call inversion5_2d(a_2d,ja,jb,sr,a1,b1,mx,ny,e)  
      call inversion5_2d(a2_2d,ja,jb,sr,a1,b1,mx,ny,c)
      tb5(:,:)=0.
      do j=1,ny-1,2
      do i=1,mx
         tb5(i,j)=e(i,(j+1)/2)
      enddo
      enddo
      do j=2,ny,2
      do i=1,mx
         tb5(i,j)=c(i,j/2)
      enddo
      enddo
      do i=1,mx
         if ((xkx(i)==0)) then
            tb5(i,1)=0.
            tb5(i,ny)=0.
         endif
      enddo
      tb5(mx,:)=0.
   endif
   if (istret==3) then
      e(:,:)=0. 
      do i=1,mx
      do j=1,ny
         e(i,j)=tb6(i,j)
      enddo
      enddo
      call inversion5_2d_i3(a3_2d,ja,jb,sr,a1,b1,mx,ny,e)  
      tb5(:,:)=0.
      do j=1,ny
      do i=1,mx
         tb5(i,j)=e(i,j)
      enddo
      enddo
       do i=1,mx
         if ((xkx(i)==0)) then
            tb5(i,1)=0.
            tb5(i,ny)=0.
         endif
      enddo
   endif
else      
   tb5(:,:)=0.
   do j=1,my
   do i=1,mx
      ytt=(bicix6*2.*cos(exs(i)*3.*dx/2.)+cicix6*2.*cos(exs(i)*5.*dx/2.))
      xtt=(biciy6*2.*cos(eys(j)*3.*dy/2.)+ciciy6*2.*cos(eys(j)*5.*dy/2.))
      yt=yk2(j)*((aicix6*2.*cos(exs(i)*dx/2.)+ytt)**2)
      xt=xk2(i)*((aiciy6*2.*cos(eys(j)*dy/2.)+xtt)**2)
      yt1=(1.+2.*ailcaix6*cos(exs(i)*dx))**2
      xt1=(1.+2.*ailcaiy6*cos(eys(j)*dy))**2
      yt2=yt/yt1
      xt2=xt/xt1
      xyzk=xt2+yt2 
      if (xyzk.eq.0) then
         tb6(i,j)=0.
      else
         tb5(i,j)=tb6(i,j)/(-xyzk)
      endif
   enddo
   enddo
endif

tb6(:,:)=0. ; tb6(:,:)=tb5(:,:) 

call slfft2d_shift(ppm,tb6,tb1,tb2,tb3,tb4,tb5,work,table,nx,nxm,ny,nym,nz,mx,my,mz,nwork,ntrigsX,ntrigsY,-1)
 
return
end subroutine poisson_2d
end module poisson_m
