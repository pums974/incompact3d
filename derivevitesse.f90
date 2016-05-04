module derivevitesse_m
implicit none
contains
!********************************************************************
!
subroutine derx(tx,ux,rx,sx,ffx,fsx,fwx,nx,ny,nz,npaire) 
!
!********************************************************************

USE param
USE derivX

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(8),dimension(nx,ny,nz) :: tx,ux,rx 
real(8),dimension(ny,nz):: sx
real(8),dimension(nx):: ffx,fsx,fwx 

if (nclx==0) then 
   do k=1,nz 
   do j=1,ny 
      tx(1,j,k)=afix*(ux(2,j,k)-ux(nx,j,k))& 
           +bfix*(ux(3,j,k)-ux(nx-1,j,k)) 
      rx(1,j,k)=-1. 
      tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
           +bfix*(ux(4,j,k)-ux(nx,j,k)) 
      rx(2,j,k)=0. 
   do i=3,nx-2 
      tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
           +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
      rx(i,j,k)=0. 
   enddo
      tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
           +bfix*(ux(1,j,k)-ux(nx-3,j,k)) 
      rx(nx-1,j,k)=0. 
      tx(nx,j,k)=afix*(ux(1,j,k)-ux(nx-1,j,k))&
           +bfix*(ux(2,j,k)-ux(nx-2,j,k)) 
      rx(nx,j,k)=alfaix           
   do i=2, nx 
      tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*fsx(i) 
   enddo
      tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      rx(nx,j,k)=rx(nx,j,k)*fwx(nx) 
   do i=nx-1,1,-1 
      tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      rx(i,j,k)=(rx(i,j,k)-ffx(i)*rx(i+1,j,k))*fwx(i) 
   enddo
      sx(j,k)=(tx(1,j,k)-alfaix*tx(nx,j,k))&
           /(1.+rx(1,j,k)-alfaix*rx(nx,j,k)) 
   do i=1,nx 
      tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k) 
   enddo
   enddo
   enddo 
endif

if (nclx==1) then 
   if (npaire==1) then 
      do k=1,nz 
      do j=1,ny 
         tx(1,j,k)=0. 
         tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
              +bfix*(ux(4,j,k)-ux(2,j,k)) 
      do i=3,nx-2 
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
              +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
      enddo
         tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
              +bfix*(ux(nx-1,j,k)-ux(nx-3,j,k)) 
         tx(nx,j,k)=0. 
      do i=2,nx 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      enddo
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      do i=nx-1,1,-1 
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then 
      do k=1,nz 
      do j=1,ny 
         tx(1,j,k)=afix*(ux(2,j,k)+ux(2,j,k))&
              +bfix*(ux(3,j,k)+ux(3,j,k)) 
         tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
              +bfix*(ux(4,j,k)+ux(2,j,k)) 
      do i=3,nx-2 
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
              +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
      enddo
         tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
              +bfix*((-ux(nx-1,j,k))-ux(nx-3,j,k)) 
         tx(nx,j,k)=afix*((-ux(nx-1,j,k))-ux(nx-1,j,k))&
              +bfix*((-ux(nx-2,j,k))-ux(nx-2,j,k)) 
      do i=2,nx 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      enddo
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      do i=nx-1,1,-1 
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      enddo
      enddo
      enddo
   endif
endif

if (nclx==2) then 
   do k=1,nz 
   do j=1,ny 
      tx(1,j,k)=af1x*ux(1,j,k)+bf1x*ux(2,j,k)+cf1x*ux(3,j,k) 
      tx(2,j,k)=af2x*(ux(3,j,k)-ux(1,j,k)) 
   do i=3,nx-2 
      tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
           +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
   enddo
      tx(nx-1,j,k)=afmx*(ux(nx,j,k)-ux(nx-2,j,k)) 
      tx(nx,j,k)=(-afnx*ux(nx,j,k))-bfnx*ux(nx-1,j,k)-cfnx*ux(nx-2,j,k) 
   do i=2,nx 
      tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
   enddo
      tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
   do i=nx-1,1,-1 
      tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
   enddo
   enddo 
   enddo 
endif

return
end subroutine derx 

!********************************************************************
!
subroutine dery(ty,uy,ry,di,sy,ffy,fsy,fwy,ppy,nx,ny,nz,npaire) 
!
!********************************************************************
  
USE param
USE derivY

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(8),dimension(nx,ny,nz) :: ty,uy 
real(8),dimension(nx,nz,ny) :: ry,di
real(8),dimension(nx,nz)  :: sy
real(8),dimension(ny) :: ffy,fsy,fwy,ppy

do k=1,nz 
do j=1,ny 
do i=1,nx 
   di(i,k,j)=uy(i,j,k) 
enddo
enddo 
enddo 

call dery1 (ry,di,ty,sy,ffy,fsy,fwy,nx,ny,nz,npaire) 

if (istret.ne.0) then   
do k=1,nz 
do j=1,ny 
do i=1,nx 
   ty(i,j,k)=ry(i,k,j)*ppy(j) 
enddo
enddo
enddo
endif
if (istret.eq.0) then   
do k=1,nz 
do j=1,ny 
do i=1,nx 
   ty(i,j,k)=ry(i,k,j)
enddo
enddo
enddo
endif

return  
end subroutine dery

!********************************************************************
!
subroutine dery1(ty,uy,ry,sy,ffy,fsy,fwy,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivY 

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(8),dimension(nx,nz,ny) :: ty,uy,ry 
real(8),dimension(nx,nz) :: sy 
real(8),dimension(ny) :: ffy,fsy,fwy 

if (ncly==0) then 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,1)=afjy*(uy(i,k,2)-uy(i,k,ny))&
           +bfjy*(uy(i,k,3)-uy(i,k,ny-1)) 
      ry(i,k,1)=-1. 
      ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
           +bfjy*(uy(i,k,4)-uy(i,k,ny)) 
      ry(i,k,2)=0. 
   enddo
   enddo 
   do j=3,ny-2 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
           +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
      ry(i,k,j)=0. 
   enddo
   enddo 
   enddo
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny-1)=afjy*(uy(i,k,ny)-uy(i,k,ny-2))&
           +bfjy*(uy(i,k,1)-uy(i,k,ny-3)) 
      ry(i,k,ny-1)=0. 
      ty(i,k,ny)=afjy*(uy(i,k,1)-uy(i,k,ny-1))&
           +bfjy*(uy(i,k,2)-uy(i,k,ny-2)) 
      ry(i,k,ny)=alfajy 
   enddo
   enddo 
   do j=2,ny 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
      ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*fsy(j) 
   enddo
   enddo 
   enddo
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
      ry(i,k,ny)=ry(i,k,ny)*fwy(ny) 
   enddo
   enddo 
   do j=ny-1,1,-1 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
      ry(i,k,j)=(ry(i,k,j)-ffy(j)*ry(i,k,j+1))*fwy(j) 
   enddo
   enddo 
   enddo
   do k=1,nz 
   do i=1,nx 
      sy(i,k)=(ty(i,k,1)-alfajy*ty(i,k,ny))&
           /(1.+ry(i,k,1)-alfajy*ry(i,k,ny)) 
   enddo
   enddo 
   do j=1,ny 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j) 
   enddo
   enddo
   enddo
endif

if (ncly==1) then 
   if (npaire==1) then 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,1)=0. 
         ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
              +bfjy*(uy(i,k,4)-uy(i,k,2)) 
      enddo
      enddo 
      do j=3,ny-2 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
              +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
      enddo
      enddo 
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny-1)=afjy*(uy(i,k,ny)-uy(i,k,ny-2))&
              +bfjy*(uy(i,k,ny-1)-uy(i,k,ny-3)) 
         ty(i,k,ny)=0. 
      enddo
      enddo 
      do j=2,ny 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
      enddo
      enddo
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
      enddo 
      enddo 
      do j=ny-1,1,-1 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
      enddo 
      enddo 
      enddo 
   endif
   if (npaire==0) then 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,1)=afjy*(uy(i,k,2)+uy(i,k,2))&
              +bfjy*(uy(i,k,3)+uy(i,k,3)) 
         ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
              +bfjy*(uy(i,k,4)+uy(i,k,2)) 
      enddo
      enddo 
      do j=3,ny-2 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
              +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
      enddo
      enddo 
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny-1)=afjy*(uy(i,k,ny)-uy(i,k,ny-2))&
              +bfjy*((-uy(i,k,ny-1))-uy(i,k,ny-3)) 
         ty(i,k,ny)=afjy*((-uy(i,k,ny-1))-uy(i,k,ny-1))&
              +bfjy*((-uy(i,k,ny-2))-uy(i,k,ny-2)) 
      enddo
      enddo 
      do j=2,ny 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
      enddo 
      enddo 
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
      enddo
      enddo 
      do j=ny-1,1,-1 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
      enddo
      enddo 
      enddo 
   endif
endif

if (ncly==2) then 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,1)=af1y*uy(i,k,1)+bf1y*uy(i,k,2)+cf1y*uy(i,k,3) 
      ty(i,k,2)=af2y*(uy(i,k,3)-uy(i,k,1)) 
   enddo
   enddo
   do j=3,ny-2 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
           +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
   enddo
   enddo 
   enddo
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny-1)=afmy*(uy(i,k,ny)-uy(i,k,ny-2)) 
      ty(i,k,ny)=(-afny*uy(i,k,ny))-bfny*uy(i,k,ny-1)-cfny*uy(i,k,ny-2) 
   enddo
   enddo 
   do j=2,ny 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
   enddo 
   enddo 
   enddo 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
   enddo 
   enddo 
   do j=ny-1,1,-1 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
   enddo 
   enddo 
   enddo
endif

return
end subroutine dery1

!********************************************************************
!
subroutine derz(tz,uz,rz,sz,ffz,fsz,fwz,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(8),dimension(nx,ny,nz) :: tz,uz,rz
real(8),dimension(nx,ny) :: sz 
real(8),dimension(nz) :: ffz,fsz,fwz

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=afkz*(uz(i,j,2)-uz(i,j,nz  ))&
           +bfkz*(uz(i,j,3)-uz(i,j,nz-1))
      rz(i,j,1)=-1.
      tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1 ))&
           +bfkz*(uz(i,j,4)-uz(i,j,nz))
      rz(i,j,2)=0.
   enddo
   enddo
   do k=3,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
               +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)=afkz*(uz(i,j,nz)-uz(i,j,nz-2))&
           +bfkz*(uz(i,j,1 )-uz(i,j,nz-3))
      rz(i,j,nz-1)=0.
      tz(i,j,nz  )=afkz*(uz(i,j,1)-uz(i,j,nz-1))&
           +bfkz*(uz(i,j,2)-uz(i,j,nz-2))
      rz(i,j,nz  )=alfakz
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*fsz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
      rz(i,j,nz)=rz(i,j,nz)*fwz(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
      rz(i,j,k)=(rz(i,j,k)-ffz(k)*rz(i,j,k+1))*fwz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(   tz(i,j,1)-alfakz*tz(i,j,nz))/&
           (1.+rz(i,j,1)-alfakz*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
   enddo
   enddo
   enddo
endif

if (nclz==1) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=0.
         tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1))&
              +bfkz*(uz(i,j,4)-uz(i,j,2))
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
              +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=afkz*(uz(i,j,nz  )-uz(i,j,nz-2))&
              +bfkz*(uz(i,j,nz-1)-uz(i,j,nz-3))
         tz(i,j,nz  )=0.
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=afkz*(uz(i,j,2)+uz(i,j,2))&
              +bfkz*(uz(i,j,3)+uz(i,j,3))
         tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1))&
              +bfkz*(uz(i,j,4)+uz(i,j,2))
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
              +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=afkz*( uz(i,j,nz  )-uz(i,j,nz-2))&
              +bfkz*(-uz(i,j,nz-1)-uz(i,j,nz-3))
         tz(i,j,nz  )=afkz*(-uz(i,j,nz-1)-uz(i,j,nz-1))&
              +bfkz*(-uz(i,j,nz-2)-uz(i,j,nz-2))
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
      enddo
      enddo
      enddo
   endif
endif

if (nclz==2) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=af1z*uz(i,j,1)+bf1z*uz(i,j,2)&
           +cf1z*uz(i,j,3)
      tz(i,j,2)=af2z*(uz(i,j,3)-uz(i,j,1))
   enddo
   enddo
   do k=3,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
           +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)= afmz*(uz(i,j,nz)-uz(i,j,nz-2))
      tz(i,j,nz  )=-afnz*uz(i,j,nz)-bfnz*uz(i,j,nz-1)&
           -cfnz*uz(i,j,nz-2)
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
   enddo
   enddo
   enddo
endif

return
end subroutine derz

!********************************************************************
!
subroutine derxx(tx,ux,rx,sx,sfx,ssx,swx,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(8),dimension(nx,ny,nz) :: tx,ux,rx
real(8),dimension(ny,nz) :: sx
real(8), dimension(nx):: sfx,ssx,swx 

   if (nclx==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=asix*(ux(2,j,k)-ux(1   ,j,k)&
                        -ux(1,j,k)+ux(nx  ,j,k))&
                  +bsix*(ux(3,j,k)-ux(1   ,j,k)&
                        -ux(1,j,k)+ux(nx-1,j,k))&
                  +csix*(ux(4,j,k)-ux(1   ,j,k)&
                        -ux(1,j,k)+ux(nx-2,j,k))
         rx(1,j,k)=-1.
         tx(2,j,k)=asix*(ux(3,j,k)-ux(2   ,j,k)&
                        -ux(2,j,k)+ux(1   ,j,k))&
                  +bsix*(ux(4,j,k)-ux(2   ,j,k)&
                        -ux(2,j,k)+ux(nx  ,j,k))&
                  +csix*(ux(5,j,k)-ux(2   ,j,k)&
                        -ux(2,j,k)+ux(nx-1,j,k))
         rx(2,j,k)=0.
         tx(3,j,k)=asix*(ux(4,j,k)-ux(3 ,j,k)&
                        -ux(3,j,k)+ux(2 ,j,k))&
                  +bsix*(ux(5,j,k)-ux(3 ,j,k)&
                        -ux(3,j,k)+ux(1 ,j,k))&
                  +csix*(ux(6,j,k)-ux(3 ,j,k)&
                        -ux(3,j,k)+ux(nx,j,k))
         rx(3,j,k)=0.
      do i=4,nx-3
         tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-1,j,k))&
                  +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-2,j,k))&
                  +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-3,j,k))
         rx(i,j,k)=0.
      enddo
         tx(nx-2,j,k)=asix*(ux(nx-1,j,k)-ux(nx-2,j,k)&
                           -ux(nx-2,j,k)+ux(nx-3,j,k))&
                     +bsix*(ux(nx  ,j,k)-ux(nx-2,j,k)&
                           -ux(nx-2,j,k)+ux(nx-4,j,k))&
                     +csix*(ux(1   ,j,k)-ux(nx-2,j,k)&
                           -ux(nx-2,j,k)+ux(nx-5,j,k))
         rx(nx-2,j,k)=0.
         tx(nx-1,j,k)=asix*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                           -ux(nx-1,j,k)+ux(nx-2,j,k))&
                     +bsix*(ux(1   ,j,k)-ux(nx-1,j,k)&
                           -ux(nx-1,j,k)+ux(nx-3,j,k))&
                     +csix*(ux(2   ,j,k)-ux(nx-1,j,k)&
                           -ux(nx-1,j,k)+ux(nx-4,j,k))
         rx(nx-1,j,k)=0.
         tx(nx  ,j,k)=asix*(ux(1 ,j,k)-ux(nx  ,j,k)&
                           -ux(nx,j,k)+ux(nx-1,j,k))&
                     +bsix*(ux(2 ,j,k)-ux(nx  ,j,k)&
                           -ux(nx,j,k)+ux(nx-2,j,k))&
                     +csix*(ux(3 ,j,k)-ux(nx  ,j,k)&
                           -ux(nx,j,k)+ux(nx-3,j,k))
         rx(nx  ,j,k)=alsaix
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*ssx(i)
      enddo
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         rx(nx,j,k)=rx(nx,j,k)*swx(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         rx(i,j,k)=(rx(i,j,k)-sfx(i)*rx(i+1,j,k))*swx(i)
      enddo
         sx(j,k)=(   tx(1,j,k)-alsaix*tx(nx,j,k))/&
                 (1.+rx(1,j,k)-alsaix*rx(nx,j,k))
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
      enddo
      enddo
   endif

   if (nclx==1) then
      if (npaire==1) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=asix*(ux(2,j,k)-ux(1,j,k)&
                           -ux(1,j,k)+ux(2,j,k))&
                     +bsix*(ux(3,j,k)-ux(1,j,k)&
                           -ux(1,j,k)+ux(3,j,k))&
                     +csix*(ux(4,j,k)-ux(1,j,k)&
                           -ux(1,j,k)+ux(4,j,k))
            tx(2,j,k)=asix*(ux(3,j,k)-ux(2,j,k)&
                           -ux(2,j,k)+ux(1,j,k))&
                     +bsix*(ux(4,j,k)-ux(2,j,k)&
                           -ux(2,j,k)+ux(2,j,k))&
                     +csix*(ux(5,j,k)-ux(2,j,k)&
                           -ux(2,j,k)+ux(3,j,k))
            tx(3,j,k)=asix*(ux(4,j,k)-ux(3,j,k)&
                           -ux(3,j,k)+ux(2,j,k))&
                     +bsix*(ux(5,j,k)-ux(3,j,k)&
                           -ux(3,j,k)+ux(1,j,k))&
                     +csix*(ux(6,j,k)-ux(3,j,k)&
                           -ux(3,j,k)+ux(2,j,k))
         do i=4,nx-3
            tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-1,j,k))&
                     +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-2,j,k))&
                     +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-3,j,k))
         enddo
            tx(nx-2,j,k)=asix*(ux(nx-1,j,k)-ux(nx-2,j,k)&
                              -ux(nx-2,j,k)+ux(nx-3,j,k))&
                        +bsix*(ux(nx  ,j,k)-ux(nx-2,j,k)&
                              -ux(nx-2,j,k)+ux(nx-4,j,k))&
                        +csix*(ux(nx-1,j,k)-ux(nx-2,j,k)&
                              -ux(nx-2,j,k)+ux(nx-5,j,k))
            tx(nx-1,j,k)=asix*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                              -ux(nx-1,j,k)+ux(nx-2,j,k))&
                        +bsix*(ux(nx-1,j,k)-ux(nx-1,j,k)&
                              -ux(nx-1,j,k)+ux(nx-3,j,k))&
                        +csix*(ux(nx-2,j,k)-ux(nx-1,j,k)&
                              -ux(nx-1,j,k)+ux(nx-4,j,k))
            tx(nx  ,j,k)=asix*(ux(nx-1,j,k)-ux(nx  ,j,k)&
                              -ux(nx  ,j,k)+ux(nx-1,j,k))&
                        +bsix*(ux(nx-2,j,k)-ux(nx  ,j,k)&
                              -ux(nx  ,j,k)+ux(nx-2,j,k))&
                        +csix*(ux(nx-3,j,k)-ux(nx  ,j,k)&
                              -ux(nx  ,j,k)+ux(nx-3,j,k))
         do i=2,nx
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
            tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         do i=nx-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=0.
            tx(2,j,k)=asix*(ux(3,j,k)-ux(2,j,k)&
                           -ux(2,j,k)+ux(1,j,k))&
                     +bsix*(ux(4,j,k)-ux(2,j,k)&
                           -ux(2,j,k)-ux(2,j,k))&
                     +csix*(ux(5,j,k)-ux(2,j,k)&
                           -ux(2,j,k)-ux(3,j,k))
            tx(3,j,k)=asix*(ux(4,j,k)-ux(3,j,k)&
                           -ux(3,j,k)+ux(2,j,k))&
                     +bsix*(ux(5,j,k)-ux(3,j,k)&
                           -ux(3,j,k)+ux(1,j,k))&
                     +csix*(ux(6,j,k)-ux(3,j,k)&
                           -ux(3,j,k)-ux(2,j,k))
         do i=4,nx-3
            tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-1,j,k))&
                     +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-2,j,k))&
                     +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-3,j,k))
         enddo
            tx(nx-2,j,k)=asix*( ux(nx-1,j,k)-ux(nx-2,j,k)&
                               -ux(nx-2,j,k)+ux(nx-3,j,k))&
                        +bsix*( ux(nx  ,j,k)-ux(nx-2,j,k)&
                               -ux(nx-2,j,k)+ux(nx-4,j,k))&
                        +csix*(-ux(nx-1,j,k)-ux(nx-2,j,k)&
                               -ux(nx-2,j,k)+ux(nx-5,j,k))
            tx(nx-1,j,k)=asix*( ux(nx  ,j,k)-ux(nx-1,j,k)&
                               -ux(nx-1,j,k)+ux(nx-2,j,k))&
                        +bsix*(-ux(nx-1,j,k)-ux(nx-1,j,k)&
                               -ux(nx-1,j,k)+ux(nx-3,j,k))&
                        +csix*(-ux(nx-2,j,k)-ux(nx-1,j,k)&
                               -ux(nx-1,j,k)+ux(nx-4,j,k))
            tx(nx  ,j,k)=0.
         do i=2,nx
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
            tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         do i=nx-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
         enddo
         enddo
      endif
   endif

   if (nclx==2) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=as1x*ux(1,j,k)+bs1x*ux(2,j,k)&
                  +cs1x*ux(3,j,k)+ds1x*ux(4,j,k)
         tx(2,j,k)=as2x*(ux(3,j,k)-ux(2,j,k)&
                        -ux(2,j,k)+ux(1,j,k))
         tx(3,j,k)=as3x*(ux(4,j,k)-ux(3,j,k)&
                        -ux(3,j,k)+ux(2,j,k))&
                  +bs3x*(ux(5,j,k)-ux(3,j,k)&
                        -ux(3,j,k)+ux(1,j,k))
      do i=4,nx-3
         tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-1,j,k))&
                  +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-2,j,k))&
                  +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-3,j,k))
      enddo
         tx(nx-2,j,k)=astx*(ux(nx-1,j,k)-ux(nx-2,j,k)&
                           -ux(nx-2,j,k)+ux(nx-3,j,k))&
                     +bstx*(ux(nx  ,j,k)-ux(nx-2,j,k)&
                           -ux(nx-2,j,k)+ux(nx-4,j,k))
         tx(nx-1,j,k)=asmx*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                           -ux(nx-1,j,k)+ux(nx-2,j,k))
         tx(nx  ,j,k)=asnx*ux(nx  ,j,k)+bsnx*ux(nx-1,j,k)&
                     +csnx*ux(nx-2,j,k)+dsnx*ux(nx-3,j,k)
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
      enddo
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
      enddo
      enddo
      enddo
   endif

return  
end subroutine derxx

!********************************************************************
!
subroutine deryy(ty,uy,ry,di,sy,sfy,ssy,swy,nx,ny,nz,npaire) 
!
!********************************************************************

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(8),dimension(nx,ny,nz) :: ty,uy
real(8),dimension(nx,nz,ny) :: ry,di
real(8),dimension(nx,nz) :: sy
real(8),dimension(ny) :: sfy,ssy,swy

do k=1,nz 
do j=1,ny 
do i=1,nx 
   di(i,k,j)=uy(i,j,k) 
enddo
enddo
enddo

call deryy1 (ry,di,ty,sy,sfy,ssy,swy,nx,ny,nz,npaire) 

do k=1,nz 
do j=1,ny 
do i=1,nx 
   ty(i,j,k)=ry(i,k,j) 
enddo
enddo
enddo

return  
end subroutine deryy

!********************************************************************
!
subroutine deryy1(ty,uy,ry,sy,sfy,ssy,swy,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivY 

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(8),dimension(nx,nz,ny) :: ty,uy,ry
real(8),dimension(nx,nz) :: sy
real(8),dimension(ny) :: sfy,ssy,swy

   if (ncly==0) then
      do k=1,nz
      do i=1,nx
         ty(i,k,1)=asjy*(uy(i,k,2)-uy(i,k,1   )&
                        -uy(i,k,1)+uy(i,k,ny  ))&
                  +bsjy*(uy(i,k,3)-uy(i,k,1   )&
                        -uy(i,k,1)+uy(i,k,ny-1))&
                  +csjy*(uy(i,k,4)-uy(i,k,1   )&
                        -uy(i,k,1)+uy(i,k,ny-2))
         ry(i,k,1)=-1.
         ty(i,k,2)=asjy*(uy(i,k,3)-uy(i,k,2 )&
                        -uy(i,k,2)+uy(i,k,1 ))&
                  +bsjy*(uy(i,k,4)-uy(i,k,2 )&
                        -uy(i,k,2)+uy(i,k,ny))&
                  +csjy*(uy(i,k,5)-uy(i,k,2 )&
                        -uy(i,k,2)+uy(i,k,ny-1))
         ry(i,k,2)=0.
         ty(i,k,3)=asjy*(uy(i,k,4)-uy(i,k,3 )&
                        -uy(i,k,3)+uy(i,k,2 ))&
                  +bsjy*(uy(i,k,5)-uy(i,k,3 )&
                        -uy(i,k,3)+uy(i,k,1 ))&
                  +csjy*(uy(i,k,6)-uy(i,k,3 )&
                        -uy(i,k,3)+uy(i,k,ny))
         ry(i,k,3)=0.
      enddo
      enddo
      do j=4,ny-3
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                        -uy(i,k,j  )+uy(i,k,j-1))&
                  +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                        -uy(i,k,j  )+uy(i,k,j-2))&
                  +csjy*(uy(i,k,j+3)-uy(i,k,j  )&
                        -uy(i,k,j  )+uy(i,k,j-3))
         ry(i,k,j)=0.
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny-2)=asjy*(uy(i,k,ny-1)-uy(i,k,ny-2)&
                           -uy(i,k,ny-2)+uy(i,k,ny-3))&
                     +bsjy*(uy(i,k,ny  )-uy(i,k,ny-2)&
                           -uy(i,k,ny-2)+uy(i,k,ny-4))&
                     +csjy*(uy(i,k,1   )-uy(i,k,ny-2)&
                           -uy(i,k,ny-2)+uy(i,k,ny-5))
         ry(i,k,ny-2)=0.
         ty(i,k,ny-1)=asjy*(uy(i,k,ny  )-uy(i,k,ny-1)&
                           -uy(i,k,ny-1)+uy(i,k,ny-2))&
                     +bsjy*(uy(i,k,1   )-uy(i,k,ny-1)&
                           -uy(i,k,ny-1)+uy(i,k,ny-3))&
                     +csjy*(uy(i,k,2   )-uy(i,k,ny-1)&
                           -uy(i,k,ny-1)+uy(i,k,ny-4))
         ry(i,k,ny-1)=0.
         ty(i,k,ny  )=asjy*(uy(i,k,1 )-uy(i,k,ny  )&
                           -uy(i,k,ny)+uy(i,k,ny-1))&
                     +bsjy*(uy(i,k,2 )-uy(i,k,ny  )&
                           -uy(i,k,ny)+uy(i,k,ny-2))&
                     +csjy*(uy(i,k,3 )-uy(i,k,ny  )&
                           -uy(i,k,ny)+uy(i,k,ny-3))
         ry(i,k,ny  )=alsajy
      enddo
      enddo
      do j=2,ny
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
         ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*ssy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny)=ty(i,k,ny)*swy(ny)
         ry(i,k,ny)=ry(i,k,ny)*swy(ny)
      enddo
      enddo
      do j=ny-1,1,-1
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
         ry(i,k,j)=(ry(i,k,j)-sfy(j)*ry(i,k,j+1))*swy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         sy(i,k)=(   ty(i,k,1)-alsajy*ty(i,k,ny))/&
                 (1.+ry(i,k,1)-alsajy*ry(i,k,ny))
      enddo
      enddo
      do j=1,ny
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j)
      enddo
      enddo
      enddo
   endif
!
   if (ncly==1) then
      if (npaire==1) then
         do k=1,nz
         do i=1,nx
            ty(i,k,1)=asjy*(uy(i,k,2)-uy(i,k,1)&
                           -uy(i,k,1)+uy(i,k,2))&
                     +bsjy*(uy(i,k,3)-uy(i,k,1)&
                           -uy(i,k,1)+uy(i,k,3))&
                     +csjy*(uy(i,k,4)-uy(i,k,1)&
                           -uy(i,k,1)+uy(i,k,4))
            ty(i,k,2)=asjy*(uy(i,k,3)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,1))&
                     +bsjy*(uy(i,k,4)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,2))&
                     +csjy*(uy(i,k,5)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,3))
            ty(i,k,3)=asjy*(uy(i,k,4)-uy(i,k,3)&
                           -uy(i,k,3)+uy(i,k,2))&
                     +bsjy*(uy(i,k,5)-uy(i,k,3)&
                           -uy(i,k,3)+uy(i,k,1))&
                     +csjy*(uy(i,k,6)-uy(i,k,3)&
                           -uy(i,k,3)+uy(i,k,2))
         enddo
         enddo
         do j=4,ny-3
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-1))&
                     +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-2))&
                     +csjy*(uy(i,k,j+3)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-3))
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny-2)=asjy*(uy(i,k,ny-1)-uy(i,k,ny-2)&
                              -uy(i,k,ny-2)+uy(i,k,ny-3))&
                        +bsjy*(uy(i,k,ny  )-uy(i,k,ny-2)&
                              -uy(i,k,ny-2)+uy(i,k,ny-4))&
                        +csjy*(uy(i,k,ny-1)-uy(i,k,ny-2)&
                              -uy(i,k,ny-2)+uy(i,k,ny-5))
            ty(i,k,ny-1)=asjy*(uy(i,k,ny  )-uy(i,k,ny-1)&
                              -uy(i,k,ny-1)+uy(i,k,ny-2))&
                        +bsjy*(uy(i,k,ny-1)-uy(i,k,ny-1)&
                              -uy(i,k,ny-1)+uy(i,k,ny-3))&
                        +csjy*(uy(i,k,ny-2)-uy(i,k,ny-1)&
                              -uy(i,k,ny-1)+uy(i,k,ny-4))
            ty(i,k,ny  )=asjy*(uy(i,k,ny-1)-uy(i,k,ny  )&
                              -uy(i,k,ny  )+uy(i,k,ny-1))&
                        +bsjy*(uy(i,k,ny-2)-uy(i,k,ny  )&
                              -uy(i,k,ny  )+uy(i,k,ny-2))&
                        +csjy*(uy(i,k,ny-3)-uy(i,k,ny  )&
                              -uy(i,k,ny  )+uy(i,k,ny-3))
         enddo
         enddo
         do j=2,ny
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny)=ty(i,k,ny)*swy(ny)
         enddo
         enddo
         do j=ny-1,1,-1
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do k=1,nz
         do i=1,nx
            ty(i,k,1)=0.
            ty(i,k,2)=asjy*(uy(i,k,3)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,1))&
                     +bsjy*(uy(i,k,4)-uy(i,k,2)&
                           -uy(i,k,2)-uy(i,k,2))&
                     +csjy*(uy(i,k,5)-uy(i,k,2)&
                           -uy(i,k,2)-uy(i,k,3))
            ty(i,k,3)=asjy*(uy(i,k,4)-uy(i,k,3)&
                           -uy(i,k,3)+uy(i,k,2))&
                     +bsjy*(uy(i,k,5)-uy(i,k,3)&
                           -uy(i,k,3)+uy(i,k,1))&
                     +csjy*(uy(i,k,6)-uy(i,k,3)&
                           -uy(i,k,3)-uy(i,k,2))
         enddo
         enddo
         do j=4,ny-3
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-1))&
                     +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-2))&
                     +csjy*(uy(i,k,j+3)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-3))
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny-2)=asjy*( uy(i,k,ny-1)-uy(i,k,ny-2)&
                               -uy(i,k,ny-2)+uy(i,k,ny-3))&
                        +bsjy*( uy(i,k,ny  )-uy(i,k,ny-2)&
                               -uy(i,k,ny-2)+uy(i,k,ny-4))&
                        +csjy*(-uy(i,k,ny-1)-uy(i,k,ny-2)&
                               -uy(i,k,ny-2)+uy(i,k,ny-5))
            ty(i,k,ny-1)=asjy*( uy(i,k,ny  )-uy(i,k,ny-1)&
                               -uy(i,k,ny-1)+uy(i,k,ny-2))&
                        +bsjy*(-uy(i,k,ny-1)-uy(i,k,ny-1)&
                               -uy(i,k,ny-1)+uy(i,k,ny-3))&
                        +csjy*(-uy(i,k,ny-2)-uy(i,k,ny-1)&
                               -uy(i,k,ny-1)+uy(i,k,ny-4))
            ty(i,k,ny  )=0.
         enddo
         enddo
         do j=2,ny
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny)=ty(i,k,ny)*swy(ny)
         enddo
         enddo
         do j=ny-1,1,-1
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
         enddo
         enddo
         enddo
      endif
   endif
!
   if (ncly==2) then
      do k=1,nz
      do i=1,nx
         ty(i,k,1)=as1y*uy(i,k,1)+bs1y*uy(i,k,2)&
                  +cs1y*uy(i,k,3)+ds1y*uy(i,k,4)
         ty(i,k,2)=as2y*(uy(i,k,3)-uy(i,k,2)&
                        -uy(i,k,2)+uy(i,k,1))
         ty(i,k,3)=as3y*(uy(i,k,4)-uy(i,k,3)&
                        -uy(i,k,3)+uy(i,k,2))&
                  +bs3y*(uy(i,k,5)-uy(i,k,3)&
                        -uy(i,k,3)+uy(i,k,1))
      enddo
      enddo
      do j=4,ny-3
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                        -uy(i,k,j  )+uy(i,k,j-1))&
                  +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                        -uy(i,k,j  )+uy(i,k,j-2))&
                  +csjy*(uy(i,k,j+3)-uy(i,k,j  )&
                        -uy(i,k,j  )+uy(i,k,j-3))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny-2)=asty*(uy(i,k,ny-1)-uy(i,k,ny-2)&
                           -uy(i,k,ny-2)+uy(i,k,ny-3))&
                     +bsty*(uy(i,k,ny  )-uy(i,k,ny-2)&
                           -uy(i,k,ny-2)+uy(i,k,ny-4))
         ty(i,k,ny-1)=asmy*(uy(i,k,ny  )-uy(i,k,ny-1)&
                           -uy(i,k,ny-1)+uy(i,k,ny-2))
         ty(i,k,ny  )=asny*uy(i,k,ny  )+bsny*uy(i,k,ny-1)&
                     +csny*uy(i,k,ny-2)+dsny*uy(i,k,ny-3)
      enddo
      enddo
      do j=2,ny
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny)=ty(i,k,ny)*swy(ny)
      enddo
      enddo
      do j=ny-1,1,-1
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
      enddo
      enddo
      enddo
   endif



return  
end subroutine deryy1

!********************************************************************
!
subroutine derzz(tz,uz,rz,sz,sfz,ssz,swz,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(8),dimension(nx,ny,nz) :: tz,uz,rz
real(8),dimension(nx,ny) :: sz 
real(8),dimension(nz) :: sfz,ssz,swz

   if (nclz==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=askz*(uz(i,j,2)-uz(i,j,1   )&
                        -uz(i,j,1)+uz(i,j,nz  ))&
                  +bskz*(uz(i,j,3)-uz(i,j,1   )&
                        -uz(i,j,1)+uz(i,j,nz-1))&
                  +cskz*(uz(i,j,4)-uz(i,j,1   )&
                        -uz(i,j,1)+uz(i,j,nz-2))
         rz(i,j,1)=-1.
         tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2 )&
                        -uz(i,j,2)+uz(i,j,1 ))&
                  +bskz*(uz(i,j,4)-uz(i,j,2 )&
                        -uz(i,j,2)+uz(i,j,nz))&
                  +cskz*(uz(i,j,5)-uz(i,j,2 )&
                        -uz(i,j,2)+uz(i,j,nz-1))
         rz(i,j,2)=0.
         tz(i,j,3)=askz*(uz(i,j,4)-uz(i,j,3 )&
                        -uz(i,j,3)+uz(i,j,2 ))&
                  +bskz*(uz(i,j,5)-uz(i,j,3 )&
                        -uz(i,j,3)+uz(i,j,1 ))&
                  +cskz*(uz(i,j,6)-uz(i,j,3 )&
                        -uz(i,j,3)+uz(i,j,nz))
         rz(i,j,3)=0.
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-1))&
                  +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-2))&
                  +cskz*(uz(i,j,k+3)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-3))
          rz(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-2)=askz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
                           -uz(i,j,nz-2)+uz(i,j,nz-3))&
                     +bskz*(uz(i,j,nz  )-uz(i,j,nz-2)&
                           -uz(i,j,nz-2)+uz(i,j,nz-4))&
                     +cskz*(uz(i,j,1   )-uz(i,j,nz-2)&
                           -uz(i,j,nz-2)+uz(i,j,nz-5))
         rz(i,j,nz-2)=0.
         tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-2))&
                     +bskz*(uz(i,j,1   )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-3))&
                     +cskz*(uz(i,j,2   )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-4))
         rz(i,j,nz-1)=0.
         tz(i,j,nz  )=askz*(uz(i,j,1 )-uz(i,j,nz  )&
                           -uz(i,j,nz)+uz(i,j,nz-1))&
                     +bskz*(uz(i,j,2 )-uz(i,j,nz  )&
                           -uz(i,j,nz)+uz(i,j,nz-2))&
                     +cskz*(uz(i,j,3 )-uz(i,j,nz  )&
                           -uz(i,j,nz)+uz(i,j,nz-3))
         rz(i,j,nz  )=alsakz
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
         rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
         rz(i,j,nz)=rz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
         rz(i,j,k)=(rz(i,j,k)-sfz(k)*rz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         sz(i,j)=(   tz(i,j,1)-alsakz*tz(i,j,nz))/&
                 (1.+rz(i,j,1)-alsakz*rz(i,j,nz))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
      enddo
      enddo
      enddo
   endif

   if (nclz==1) then
      if (npaire==1) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=askz*(uz(i,j,2)-uz(i,j,1)&
                           -uz(i,j,1)+uz(i,j,2))&
                     +bskz*(uz(i,j,3)-uz(i,j,1)&
                           -uz(i,j,1)+uz(i,j,3))&
                     +cskz*(uz(i,j,4)-uz(i,j,1)&
                           -uz(i,j,1)+uz(i,j,4))
            tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,1))&
                     +bskz*(uz(i,j,4)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,2))&
                     +cskz*(uz(i,j,5)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,3))
            tz(i,j,3)=askz*(uz(i,j,4)-uz(i,j,3)&
                           -uz(i,j,3)+uz(i,j,2))&
                     +bskz*(uz(i,j,5)-uz(i,j,3)&
                           -uz(i,j,3)+uz(i,j,1))&
                     +cskz*(uz(i,j,6)-uz(i,j,3)&
                           -uz(i,j,3)+uz(i,j,2))
         enddo
         enddo
         do k=4,nz-3
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-1))&
                     +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-2))&
                     +cskz*(uz(i,j,k+3)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-3))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-2)=askz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
                              -uz(i,j,nz-2)+uz(i,j,nz-3))&
                        +bskz*(uz(i,j,nz  )-uz(i,j,nz-2)&
                              -uz(i,j,nz-2)+uz(i,j,nz-4))&
                        +cskz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
                              -uz(i,j,nz-2)+uz(i,j,nz-5))
            tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-2))&
                        +bskz*(uz(i,j,nz-1)-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-3))&
                        +cskz*(uz(i,j,nz-2)-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-4))
            tz(i,j,nz  )=askz*(uz(i,j,nz-1)-uz(i,j,nz  )&
                              -uz(i,j,nz  )+uz(i,j,nz-1))&
                        +bskz*(uz(i,j,nz-2)-uz(i,j,nz  )&
                              -uz(i,j,nz  )+uz(i,j,nz-2))&
                        +cskz*(uz(i,j,nz-3)-uz(i,j,nz  )&
                              -uz(i,j,nz  )+uz(i,j,nz-3))
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*swz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=0.
            tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,1))&
                     +bskz*(uz(i,j,4)-uz(i,j,2)&
                           -uz(i,j,2)-uz(i,j,2))&
                     +cskz*(uz(i,j,5)-uz(i,j,2)&
                           -uz(i,j,2)-uz(i,j,3))
            tz(i,j,3)=askz*(uz(i,j,4)-uz(i,j,3)&
                           -uz(i,j,3)+uz(i,j,2))&
                     +bskz*(uz(i,j,5)-uz(i,j,3)&
                           -uz(i,j,3)+uz(i,j,1))&
                     +cskz*(uz(i,j,6)-uz(i,j,3)&
                           -uz(i,j,3)-uz(i,j,2))
         enddo
         enddo
         do k=4,nz-3
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-1))&
                     +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-2))&
                     +cskz*(uz(i,j,k+3)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-3))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-2)=askz*( uz(i,j,nz-1)-uz(i,j,nz-2)&
                               -uz(i,j,nz-2)+uz(i,j,nz-3))&
                        +bskz*( uz(i,j,nz  )-uz(i,j,nz-2)&
                               -uz(i,j,nz-2)+uz(i,j,nz-4))&
                        +cskz*(-uz(i,j,nz-1)-uz(i,j,nz-2)&
                               -uz(i,j,nz-2)+uz(i,j,nz-5))
            tz(i,j,nz-1)=askz*( uz(i,j,nz  )-uz(i,j,nz-1)&
                               -uz(i,j,nz-1)+uz(i,j,nz-2))&
                        +bskz*(-uz(i,j,nz-1)-uz(i,j,nz-1)&
                               -uz(i,j,nz-1)+uz(i,j,nz-3))&
                        +cskz*(-uz(i,j,nz-2)-uz(i,j,nz-1)&
                               -uz(i,j,nz-1)+uz(i,j,nz-4))
            tz(i,j,nz  )=0.
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*swz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
         enddo
         enddo
         enddo
      endif
   endif

   if (nclz==2) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=as1z*uz(i,j,1)+bs1z*uz(i,j,2)&
                  +cs1z*uz(i,j,3)+ds1z*uz(i,j,4)
         tz(i,j,2)=as2z*(uz(i,j,3)-uz(i,j,2)&
                        -uz(i,j,2)+uz(i,j,1))
         tz(i,j,3)=as3z*(uz(i,j,4)-uz(i,j,3)&
                        -uz(i,j,3)+uz(i,j,2))&
                  +bs3z*(uz(i,j,5)-uz(i,j,3)&
                        -uz(i,j,3)+uz(i,j,1))
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-1))&
                  +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-2))&
                  +cskz*(uz(i,j,k+3)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-3))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-2)=astz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
                           -uz(i,j,nz-2)+uz(i,j,nz-3))&
                     +bstz*(uz(i,j,nz  )-uz(i,j,nz-2)&
                           -uz(i,j,nz-2)+uz(i,j,nz-4))
         tz(i,j,nz-1)=asmz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-2))
         tz(i,j,nz  )=asnz*uz(i,j,nz  )+bsnz*uz(i,j,nz-1)&
                     +csnz*uz(i,j,nz-2)+dsnz*uz(i,j,nz-3)
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
   endif



return  
end subroutine derzz

!********************************************************************
!
subroutine decx6(tx,ux,rx,sx,cfx6,csx6,cwx6,nx,nxm,ny,nz,npaire)  
!
!********************************************************************

USE param 
USE derivX 

implicit none

integer :: nx,nxm,ny,nz,npaire
real(8),dimension(nxm,ny,nz) :: tx
real(8),dimension(nx,ny,nz) :: ux,rx
real(8),dimension(ny,nz) :: sx
real(8),dimension(nxm) :: cfx6,csx6,cwx6
integer :: i,j,k,nyz

nyz=ny*nz

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=acix6*(ux(2,j,k)-ux(1  ,j,k))&
           +bcix6*(ux(3,j,k)-ux(nx,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=acix6*(ux(3,j,k)-ux(2 ,j,k))&
           +bcix6*(ux(4,j,k)-ux(1,j,k))
      rx(2,j,k)=0.
   do i=3,nx-2
      tx(i,j,k)=acix6*(ux(i+1,j,k)-ux(i,j,k))&
           +bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
      rx(i,j,k)=0.
   enddo
      tx(nx-1,j,k)=acix6*(ux(nx,j,k)-ux(nx-1,j,k))&
           +bcix6*(ux(1 ,j,k)-ux(nx-2,j,k))
      rx(nx-1,j,k)=0.
      tx(nx  ,j,k)=acix6*(ux(1,j,k)-ux(nx,j,k))&
           +bcix6*(ux(2,j,k)-ux(nx-1,j,k))
      rx(nx  ,j,k)=alcaix6
   do i=2,nx
      tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
      rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*csx6(i)
   enddo
      tx(nx,j,k)=tx(nx,j,k)*cwx6(nx)
      rx(nx,j,k)=rx(nx,j,k)*cwx6(nx)
   do i=nx-1,1,-1
      tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
      rx(i,j,k)=(rx(i,j,k)-cfx6(i)*rx(i+1,j,k))*cwx6(i)
   enddo
      sx(j,k)=(tx(1,j,k)-alcaix6*tx(nx,j,k))/&
           (1.+rx(1,j,k)-alcaix6*rx(nx,j,k))
   do i=1,nx
      tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
   enddo
   enddo
   enddo
endif

if ((nclx==1).or.(nclx==2)) then
   if (npaire==1) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=acix6*(ux(2,j,k)-ux(1,j,k))&
              +bcix6*(ux(3,j,k)-ux(2,j,k))             
         tx(2,j,k)=acix6*(ux(3,j,k)-ux(2,j,k))&
              +bcix6*(ux(4,j,k)-ux(1,j,k)) 
         do i=3,nxm-2
            tx(i,j,k)=acix6*(ux(i+1,j,k)-ux(i,j,k))&
                 +bcix6*(ux(i+2,j,k)-ux(i-1,j,k))              
         enddo
         tx(nxm-1,j,k)=acix6*(ux(nxm,j,k)-ux(nxm-1,j,k))&
              +bcix6*(ux(nx,j,k)-ux(nxm-2,j,k))         
         tx(nxm,j,k)=acix6*(ux(nx,j,k)-ux(nxm,j,k))&
              +bcix6*(ux(nxm,j,k)-ux(nxm-1,j,k))
         do i=2,nxm
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
         enddo
         tx(nxm,j,k)=tx(nxm,j,k)*cwx6(nxm)
         do i=nxm-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
         enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=acix6*(ux(2,j,k)-ux(1,j,k))&
              +bcix6*(ux(3,j,k)-2.*ux(1,j,k)+ux(2,j,k)) 
         tx(2,j,k)=acix6*(ux(3,j,k)-ux(2,j,k))&
              +bcix6*(ux(4,j,k)-ux(1,j,k))
         do i=3,nxm-2
            tx(i,j,k)=acix6*(ux(i+1,j,k)-ux(i,j,k))&
                 +bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
         enddo
         tx(nxm-1,j,k)=acix6*(ux(nxm,j,k)-ux(nxm-1,j,k))&
              +bcix6*(ux(nx,j,k)-ux(nxm-2,j,k)) 
         tx(nxm,j,k)=acix6*(ux(nx,j,k)-ux(nxm,j,k))&
              +bcix6*(2.*ux(nx,j,k)-ux(nxm,j,k)-ux(nxm-1,j,k)) 
         do i=2,nxm
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
         enddo
         tx(nxm,j,k)=tx(nxm,j,k)*cwx6(nxm)
         do i=nxm-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
         enddo
      enddo
      enddo
   endif
endif

return     
end subroutine decx6

!********************************************************************
!
subroutine inter6(tx,ux,rx,sx,cifx6,cisx6,ciwx6,nx,nxm,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivX 

implicit none

integer :: nx,nxm,ny,nz,npaire,i,j,nyz,k
real(8),dimension(nxm,ny,nz) :: tx
real(8),dimension(nx,ny,nz) :: ux,rx
real(8),dimension(ny,nz) :: sx
real(8),dimension(nxm) :: cifx6,cisx6,ciwx6

nyz=ny*nz

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=aicix6*(ux(2,j,k)+ux(1  ,j,k))&
           +bicix6*(ux(3,j,k)+ux(nx,j,k))&
           +cicix6*(ux(4,j,k)+ux(nx-1,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=aicix6*(ux(3,j,k)+ux(2 ,j,k))&
           +bicix6*(ux(4,j,k)+ux(1,j,k))&
           +cicix6*(ux(5,j,k)+ux(nx,j,k))
      rx(2,j,k)=0.
   do i=3,nx-3   
      tx(i,j,k)=aicix6*(ux(i+1,j,k)+ux(i,j,k))&
           +bicix6*(ux(i+2,j,k)+ux(i-1,j,k))&
           +cicix6*(ux(i+3,j,k)+ux(i-2,j,k))
      rx(i,j,k)=0.
   enddo
      tx(nx-2,j,k)=aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k))&
           +bicix6*(ux(nx ,j,k)+ux(nx-3,j,k))&
           +cicix6*(ux(1,j,k)+ux(nx-4,j,k))
      rx(nx-2,j,k)=0.
      tx(nx-1,j,k)=aicix6*(ux(nx,j,k)+ux(nx-1,j,k))&
           +bicix6*(ux(1 ,j,k)+ux(nx-2,j,k))&
           +cicix6*(ux(2,j,k)+ux(nx-3,j,k))
      rx(nx-1,j,k)=0.
      tx(nx  ,j,k)=aicix6*(ux(1,j,k)+ux(nx,j,k))&
           +bicix6*(ux(2,j,k)+ux(nx-1,j,k))&
           +cicix6*(ux(3,j,k)+ux(nx-2,j,k))
      rx(nx  ,j,k)=ailcaix6
   do i=2,nx
      tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*cisx6(i)
      rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*cisx6(i)
   enddo
      tx(nx,j,k)=tx(nx,j,k)*ciwx6(nx)
      rx(nx,j,k)=rx(nx,j,k)*ciwx6(nx)
   do i=nx-1,1,-1
      tx(i,j,k)=(tx(i,j,k)-cifx6(i)*tx(i+1,j,k))*ciwx6(i)
      rx(i,j,k)=(rx(i,j,k)-cifx6(i)*rx(i+1,j,k))*ciwx6(i)
   enddo
      sx(j,k)=(tx(1,j,k)-ailcaix6*tx(nx,j,k))/&
           (1.+rx(1,j,k)-ailcaix6*rx(nx,j,k))
   do i=1,nx
      tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
   enddo
   enddo
   enddo
endif

if ((nclx==1).or.(nclx==2)) then
   if (npaire==1) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=aicix6*(ux(2,j,k)+ux(1,j,k))&
              +bicix6*(ux(3,j,k)+ux(2,j,k))&
              +cicix6*(ux(4,j,k)+ux(3,j,k))             
         tx(2,j,k)=aicix6*(ux(3,j,k)+ux(2,j,k))&
              +bicix6*(ux(4,j,k)+ux(1,j,k))&
              +cicix6*(ux(5,j,k)+ux(2,j,k))               
         do i=3,nxm-2
         tx(i,j,k)=aicix6*(ux(i+1,j,k)+ux(i,j,k))&
              +bicix6*(ux(i+2,j,k)+ux(i-1,j,k))&
              +cicix6*(ux(i+3,j,k)+ux(i-2,j,k))
      enddo
         tx(nxm-1,j,k)=aicix6*(ux(nxm,j,k)+ux(nxm-1,j,k))&
              +bicix6*(ux(nx,j,k)+ux(nxm-2,j,k))&
              +cicix6*(ux(nxm,j,k)+ux(nxm-3,j,k))
         tx(nxm,j,k)=aicix6*(ux(nx,j,k)+ux(nxm,j,k))&
              +bicix6*(ux(nxm,j,k)+ux(nxm-1,j,k))&
              +cicix6*(ux(nxm-1,j,k)+ux(nxm-2,j,k))            
      do i=2,nxm
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*cisx6(i)
      enddo
         tx(nxm,j,k)=tx(nxm,j,k)*ciwx6(nxm)
      do i=nxm-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-cifx6(i)*tx(i+1,j,k))*ciwx6(i)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine inter6 

!********************************************************************
!
subroutine deci6(tx,ux,rx,sx,cfi6,csi6,cwi6,cfx6,csx6,cwx6,nxm,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivX 

implicit none

integer :: nx,nxm,ny,nz,npaire
real(8),dimension(nx,ny,nz) :: tx
real(8),dimension(nxm,ny,nz) :: ux,rx
real(8),dimension(ny,nz) :: sx
real(8),dimension(nx) :: cfi6,csi6,cwi6
real(8),dimension(nx) :: cfx6,csx6,cwx6
integer :: i,j,k

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=acix6*(ux(1,j,k)-ux(nx  ,j,k))&
           +bcix6*(ux(2,j,k)-ux(nx-1,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=acix6*(ux(2,j,k)-ux(1 ,j,k))&
           +bcix6*(ux(3,j,k)-ux(nx,j,k))
      rx(2,j,k)=0.
   do i=3,nx-2
      tx(i,j,k)=acix6*(ux(i,j,k)-ux(i-1,j,k))&
           +bcix6*(ux(i+1,j,k)-ux(i-2,j,k))
      rx(i,j,k)=0.
   enddo
      tx(nx-1,j,k)=acix6*(ux(nx-1,j,k)-ux(nx-2,j,k))&
           +bcix6*(ux(nx ,j,k)-ux(nx-3,j,k))
      rx(nx-1,j,k)=0.
      tx(nx  ,j,k)=acix6*(ux(nx,j,k)-ux(nx-1,j,k))&
           +bcix6*(ux(1,j,k)-ux(nx-2,j,k))
      rx(nx  ,j,k)=alcaix6
   do i=2,nx
      tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
      rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*csx6(i)
   enddo
      tx(nx,j,k)=tx(nx,j,k)*cwx6(nx)
      rx(nx,j,k)=rx(nx,j,k)*cwx6(nx)
   do i=nx-1,1,-1
      tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
      rx(i,j,k)=(rx(i,j,k)-cfx6(i)*rx(i+1,j,k))*cwx6(i)
   enddo
      sx(j,k)=(tx(1,j,k)-alcaix6*tx(nx,j,k))/&
           (1.+rx(1,j,k)-alcaix6*rx(nx,j,k))
   do i=1,nx
      tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
   enddo
   enddo
   enddo
endif

if ((nclx==1).or.(nclx==2)) then
   if (npaire==1) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=0.  
         tx(2,j,k)=acix6*(ux(2,j,k)-ux(1,j,k))&
              +bcix6*(ux(3,j,k)-ux(1,j,k))  
      do i=3,nx-2
         tx(i,j,k)=acix6*(ux(i,j,k)-ux(i-1,j,k))&
              +bcix6*(ux(i+1,j,k)-ux(i-2,j,k))
      enddo
         tx(nx-1,j,k)=acix6*(ux(nx-1,j,k)-ux(nx-2,j,k))&
              +bcix6*(ux(nx-1,j,k)-ux(nx-3,j,k))
         tx(nx,j,k)=0.
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csi6(i)
      enddo
         tx(nx,j,k)=tx(nx,j,k)*cwi6(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-cfi6(i)*tx(i+1,j,k))*cwi6(i)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine deci6 

!********************************************************************
!
subroutine interi6(tx,ux,rx,sx,cifi6,cisi6,ciwi6,cifx6,cisx6,ciwx6,nxm,nx,ny,nz,npaire) 
!
!********************************************************************
  
USE param 
USE derivX 

implicit none

integer :: nx,nxm,ny,nz,npaire
real(8),dimension(nx,ny,nz) :: tx,rx
real(8),dimension(nxm,ny,nz) :: ux
real(8),dimension(ny,nz) :: sx
real(8),dimension(nx) :: cifi6,cisi6,ciwi6
real(8),dimension(nx) :: cifx6,cisx6,ciwx6
integer :: i,j,k

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=aicix6*(ux(1,j,k)+ux(nx  ,j,k))&
           +bicix6*(ux(2,j,k)+ux(nx-1,j,k))&
           +cicix6*(ux(3,j,k)+ux(nx-2,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=aicix6*(ux(2,j,k)+ux(1 ,j,k))&
           +bicix6*(ux(3,j,k)+ux(nx,j,k))&
           +cicix6*(ux(4,j,k)+ux(nx-1,j,k))
      rx(2,j,k)=0.
      tx(3,j,k)=aicix6*(ux(3,j,k)+ux(2 ,j,k))&
           +bicix6*(ux(4,j,k)+ux(1,j,k))&
           +cicix6*(ux(5,j,k)+ux(nx,j,k))
      rx(3,j,k)=0.
   do i=4,nx-2
      tx(i,j,k)=aicix6*(ux(i,j,k)+ux(i-1,j,k))&
           +bicix6*(ux(i+1,j,k)+ux(i-2,j,k))&
           +cicix6*(ux(i+2,j,k)+ux(i-3,j,k))
      rx(i,j,k)=0.
   enddo
      tx(nx-1,j,k)=aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k))&
           +bicix6*(ux(nx ,j,k)+ux(nx-3,j,k))&
           +cicix6*(ux(1,j,k)+ux(nx-4,j,k))
      rx(nx-1,j,k)=0.
      tx(nx  ,j,k)=aicix6*(ux(nx,j,k)+ux(nx-1,j,k))&
           +bicix6*(ux(1,j,k)+ux(nx-2,j,k))&
           +cicix6*(ux(2,j,k)+ux(nx-3,j,k))
      rx(nx  ,j,k)=ailcaix6
   do i=2,nx
      tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*cisx6(i)
      rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*cisx6(i)
   enddo
      tx(nx,j,k)=tx(nx,j,k)*ciwx6(nx)
      rx(nx,j,k)=rx(nx,j,k)*ciwx6(nx)
   do i=nx-1,1,-1
      tx(i,j,k)=(tx(i,j,k)-cifx6(i)*tx(i+1,j,k))*ciwx6(i)
      rx(i,j,k)=(rx(i,j,k)-cifx6(i)*rx(i+1,j,k))*ciwx6(i)
   enddo
      sx(j,k)=(tx(1,j,k)-ailcaix6*tx(nx,j,k))/&
           (1.+rx(1,j,k)-ailcaix6*rx(nx,j,k))
   do i=1,nx
      tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
   enddo
   enddo
   enddo
endif

if ((nclx==1).or.(nclx==2)) then
   if (npaire==1) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=aicix6*(ux(1,j,k)+ux(1,j,k))&
              +bicix6*(ux(2,j,k)+ux(2,j,k))&
              +cicix6*(ux(3,j,k)+ux(3,j,k))     
         tx(2,j,k)=aicix6*(ux(2,j,k)+ux(1,j,k))&
              +bicix6*(ux(3,j,k)+ux(1,j,k))&
              +cicix6*(ux(4,j,k)+ux(2,j,k))
         tx(3,j,k)=aicix6*(ux(3,j,k)+ux(2,j,k))&
              +bicix6*(ux(4,j,k)+ux(1,j,k))&
              +cicix6*(ux(5,j,k)+ux(1,j,k))
      do i=4,nx-3
         tx(i,j,k)=aicix6*(ux(i,j,k)+ux(i-1,j,k))&
              +bicix6*(ux(i+1,j,k)+ux(i-2,j,k))&
              +cicix6*(ux(i+2,j,k)+ux(i-3,j,k))
      enddo
         tx(nx-2,j,k)=aicix6*(ux(nx-2,j,k)+ux(nx-3,j,k))&
              +bicix6*(ux(nx-1,j,k)+ux(nx-4,j,k))&
              +cicix6*(ux(nx-1,j,k)+ux(nx-5,j,k))
         tx(nx-1,j,k)=aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k))&
              +bicix6*(ux(nx-1,j,k)+ux(nx-3,j,k))&
              +cicix6*(ux(nx-2,j,k)+ux(nx-4,j,k))
         tx(nx,j,k)=aicix6*(ux(nx-1,j,k)+ux(nx-1,j,k))&
              +bicix6*(ux(nx-2,j,k)+ux(nx-2,j,k))&
              +cicix6*(ux(nx-3,j,k)+ux(nx-3,j,k))
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*cisi6(i)
      enddo
         tx(nx,j,k)=tx(nx,j,k)*ciwi6(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-cifi6(i)*tx(i+1,j,k))*ciwi6(i)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine interi6

!********************************************************************
!
subroutine intery6(ty,uy,ry,di,sy,cify6,cisy6,ciwy6,nx,ny,nym,nz,npaire) 
!
!********************************************************************

implicit none

integer :: nx,ny,nym,nz,npaire
real(8),dimension(nx,nym,nz) :: ty
real(8),dimension(nx,ny,nz) :: uy
real(8),dimension(nx,nz,nym) :: ry
real(8),dimension(nx,nz,ny) :: di
real(8),dimension(nx,nz) :: sy
real(8),dimension(nym) :: cify6,cisy6,ciwy6
integer :: i,j,k

do k=1,nz
do j=1,ny
do i=1,nx
   di(i,k,j)=uy(i,j,k)
enddo
enddo
enddo
         
call intery16(ry,di,ty,sy,cify6,cisy6,ciwy6,nx,ny,nym,nz,npaire)

do k=1,nz
do j=1,nym
do i=1,nx
   ty(i,j,k)=ry(i,k,j)
enddo
enddo
enddo

return
end subroutine intery6

!********************************************************************
!
subroutine intery16(ty,uy,ry,sy,cify6,cisy6,ciwy6,nx,ny,nym,nz,npaire) 
!
!********************************************************************
  
USE param 
USE derivY 

implicit none

integer :: nx,ny,nym,nz,npaire
real(8),dimension(nx,nz,nym) :: ty
real(8),dimension(nx,nz,ny) :: uy,ry
real(8),dimension(nx,nz) :: sy
real(8),dimension(nym) :: cify6,cisy6,ciwy6
integer :: i,j,k

if (ncly==0) then
   do k=1,nz 
   do i=1,nx 
      ty(i,k,1)=aiciy6*(uy(i,k,2)+uy(i,k,1))&
           +biciy6*(uy(i,k,3)+uy(i,k,ny))&
           +ciciy6*(uy(i,k,4)+uy(i,k,ny-1)) 
      ry(i,k,1)=-1.
      ty(i,k,2)=aiciy6*(uy(i,k,3)+uy(i,k,2))&
           +biciy6*(uy(i,k,4)+uy(i,k,1))&
           +ciciy6*(uy(i,k,5)+uy(i,k,ny))  
      ry(i,k,2)=0. 
   enddo
   enddo
   do j=3,ny-3
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=aiciy6*(uy(i,k,j+1)+uy(i,k,j))&
           +biciy6*(uy(i,k,j+2)+uy(i,k,j-1))&
           +ciciy6*(uy(i,k,j+3)+uy(i,k,j-2))  
      ry(i,k,j)=0. 
   enddo
   enddo
   enddo
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny-2)=aiciy6*(uy(i,k,ny-1)+uy(i,k,ny-2))&
           +biciy6*(uy(i,k,ny)+uy(i,k,ny-3))&
           +ciciy6*(uy(i,k,1)+uy(i,k,ny-4)) 
      ry(i,k,ny-2)=0. 
      ty(i,k,ny-1)=aiciy6*(uy(i,k,ny)+uy(i,k,ny-1))&
           +biciy6*(uy(i,k,1)+uy(i,k,ny-2))&
           +ciciy6*(uy(i,k,2)+uy(i,k,ny-3)) 
      ry(i,k,ny-1)=0. 
      ty(i,k,ny)=aiciy6*(uy(i,k,1)+uy(i,k,ny))&
           +biciy6*(uy(i,k,2)+uy(i,k,ny-1))&
           +ciciy6*(uy(i,k,3)+uy(i,k,ny-2))  
      ry(i,k,ny)=ailcaiy6
   enddo
   enddo
   do j=2,ny 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*cisy6(j) 
      ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*cisy6(j) 
   enddo
   enddo
   enddo
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny)=ty(i,k,ny)*ciwy6(ny) 
      ry(i,k,ny)=ry(i,k,ny)*ciwy6(ny) 
   enddo
   enddo
   do j=ny-1,1,-1 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=(ty(i,k,j)-cify6(j)*ty(i,k,j+1))*ciwy6(j) 
      ry(i,k,j)=(ry(i,k,j)-cify6(j)*ry(i,k,j+1))*ciwy6(j) 
   enddo
   enddo
   enddo
   do k=1,nz 
   do i=1,nx 
      sy(i,k)=(ty(i,k,1)-ailcaiy6*ty(i,k,ny))&
           /(1.+ry(i,k,1)-ailcaiy6*ry(i,k,ny)) 
   enddo 
   enddo
   do j=1,ny 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j) 
   enddo
   enddo
   enddo
endif
if ((ncly==1).or.(ncly==2)) then
   if (npaire==1) then
      do k=1,nz
      do i=1,nx
         ty(i,k,1)=aiciy6*(uy(i,k,2)+uy(i,k,1))&
              +biciy6*(uy(i,k,3)+uy(i,k,2))&
              +ciciy6*(uy(i,k,4)+uy(i,k,3))  
         ty(i,k,2)=aiciy6*(uy(i,k,3)+uy(i,k,2))&
              +biciy6*(uy(i,k,4)+uy(i,k,1))&
              +ciciy6*(uy(i,k,5)+uy(i,k,2))               
      enddo
      enddo
      do j=3,nym-3
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=aiciy6*(uy(i,k,j+1)+uy(i,k,j))&
              +biciy6*(uy(i,k,j+2)+uy(i,k,j-1))&
              +ciciy6*(uy(i,k,j+3)+uy(i,k,j-2))                
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,nym-2)=aiciy6*(uy(i,k,nym-1)+uy(i,k,nym-2))&
              +biciy6*(uy(i,k,nym)+uy(i,k,nym-3))&
              +ciciy6*(uy(i,k,ny)+uy(i,k,nym-4))  
         ty(i,k,nym-1)=aiciy6*(uy(i,k,nym)+uy(i,k,nym-1))&
              +biciy6*(uy(i,k,ny)+uy(i,k,nym-2))&
              +ciciy6*(uy(i,k,nym)+uy(i,k,nym-3))  
         ty(i,k,nym)=aiciy6*(uy(i,k,ny)+uy(i,k,nym))&
              +biciy6*(uy(i,k,nym)+uy(i,k,nym-1))&
              +ciciy6*(uy(i,k,nym-1)+uy(i,k,nym-2))              
      enddo
      enddo
      do j=2,nym
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*cisy6(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,nym)=ty(i,k,nym)*ciwy6(nym)
      enddo
      enddo
      do j=nym-1,1,-1
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=(ty(i,k,j)-cify6(j)*ty(i,k,j+1))*ciwy6(j)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine intery16

!********************************************************************
!
subroutine decy6(ty,uy,ry,di,sy,cfy6,csy6,cwy6,ppyi,nx,ny,nym,nz,npaire) 
!
!********************************************************************
  
USE param
USE derivY

implicit none

integer :: nx,ny,nym,nz,npaire
real(8),dimension(nx,nym,nz) :: ty
real(8),dimension(nx,ny,nz) :: uy
real(8),dimension(nx,nz,ny) :: di
real(8),dimension(nx,nz,nym) :: ry
real(8),dimension(nx,nz) :: sy
real(8),dimension(nym) :: cfy6,csy6,cwy6,ppyi
integer :: i,j,k

do k=1,nz
do j=1,ny
do i=1,nx
   di(i,k,j)=uy(i,j,k)
enddo
enddo
enddo
         
call decy16(ry,di,ty,sy,cfy6,csy6,cwy6,nx,ny,nym,nz,npaire)

if (istret.ne.0) then   
do k=1,nz 
do j=1,nym 
do i=1,nx 
   ty(i,j,k)=ry(i,k,j)*ppyi(j) 
enddo
enddo
enddo
endif
if (istret.eq.0) then   
do k=1,nz 
do j=1,nym
do i=1,nx 
   ty(i,j,k)=ry(i,k,j)
enddo
enddo
enddo
endif

return
end subroutine decy6 

!********************************************************************
!
subroutine decy16(ty,uy,ry,sy,cfy6,csy6,cwy6,nx,ny,nym,nz,npaire) 
!
!********************************************************************
  
USE param
USE derivY 

implicit none

integer :: nx,ny,nym,nz,npaire
real(8),dimension(nx,nz,nym) :: ty
real(8),dimension(nx,nz,ny) :: uy
real(8),dimension(nx,nz,ny) :: ry
real(8),dimension(nx,nz) :: sy
real(8),dimension(nym) :: cfy6,csy6,cwy6
integer :: i,j,k

if (ncly==0) then
   do k=1,nz 
   do i=1,nx 
      ty(i,k,1)=aciy6*(uy(i,k,2)-uy(i,k,1))&
           +bciy6*(uy(i,k,3)-uy(i,k,ny)) 
      ry(i,k,1)=-1.
      ty(i,k,2)=aciy6*(uy(i,k,3)-uy(i,k,2))&
           +bciy6*(uy(i,k,4)-uy(i,k,1)) 
      ry(i,k,2)=0. 
   enddo
   enddo 
   do j=3,ny-2 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=aciy6*(uy(i,k,j+1)-uy(i,k,j))&
           +bciy6*(uy(i,k,j+2)-uy(i,k,j-1)) 
      ry(i,k,j)=0.
   enddo
   enddo
   enddo
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny-1)=aciy6*(uy(i,k,ny)-uy(i,k,ny-1))&
           +bciy6*(uy(i,k,1)-uy(i,k,ny-2)) 
      ry(i,k,ny-1)=0. 
      ty(i,k,ny)=aciy6*(uy(i,k,1)-uy(i,k,ny))&
           +bciy6*(uy(i,k,2)-uy(i,k,ny-1)) 
      ry(i,k,ny)=alcaiy6
   enddo
   enddo
   do j=2,ny 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*csy6(j) 
      ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*csy6(j) 
   enddo
   enddo
   enddo
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny)=ty(i,k,ny)*cwy6(ny) 
      ry(i,k,ny)=ry(i,k,ny)*cwy6(ny) 
   enddo
   enddo
   do j=ny-1,1,-1 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=(ty(i,k,j)-cfy6(j)*ty(i,k,j+1))*cwy6(j) 
      ry(i,k,j)=(ry(i,k,j)-cfy6(j)*ry(i,k,j+1))*cwy6(j) 
   enddo
   enddo
   enddo
   do k=1,nz 
   do i=1,nx 
      sy(i,k)=(ty(i,k,1)-alcaiy6*ty(i,k,ny))&
           /(1.+ry(i,k,1)-alcaiy6*ry(i,k,ny)) 
   enddo
   enddo
   do j=1,ny 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j) 
   enddo
   enddo
   enddo
endif
if ((ncly==1).or.(ncly==2)) then
   if (npaire==0) then
      do k=1,nz
      do i=1,nx
         ty(i,k,1)=aciy6*(uy(i,k,2)-uy(i,k,1))&
              +bciy6*(uy(i,k,3)-2.*uy(i,k,1)+uy(i,k,2)) 
         ty(i,k,2)=aciy6*(uy(i,k,3)-uy(i,k,2))&
              +bciy6*(uy(i,k,4)-uy(i,k,1))
      enddo
      enddo
      do j=3,nym-2
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=aciy6*(uy(i,k,j+1)-uy(i,k,j))&
              +bciy6*(uy(i,k,j+2)-uy(i,k,j-1))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,nym-1)=aciy6*(uy(i,k,nym)-uy(i,k,nym-1))&
              +bciy6*(uy(i,k,ny)-uy(i,k,nym-2)) 
         ty(i,k,nym)=aciy6*(uy(i,k,ny)-uy(i,k,nym))&
              +bciy6*(2.*uy(i,k,ny)-uy(i,k,nym)-uy(i,k,nym-1)) 
      enddo
      enddo
      do j=2,nym
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*csy6(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,nym)=ty(i,k,nym)*cwy6(nym)
      enddo
      enddo
      do j=nym-1,1,-1
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=(ty(i,k,j)-cfy6(j)*ty(i,k,j+1))*cwy6(j)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine decy16

!********************************************************************
!
subroutine interiy6(ty,uy,ry,di,sy,cifi6y,cisi6y,ciwi6y,cify6,cisy6,ciwy6,nx,nym,ny,nz,npaire) 
!
!********************************************************************
 
implicit none

integer :: nx,ny,nym,nz,npaire
real(8),dimension(nx,ny,nz) :: ty
real(8),dimension(nx,nym,nz) :: uy
real(8),dimension(nx,nz,ny) :: ry
real(8),dimension(nx,nz,nym) :: di
real(8),dimension(nx,nz) :: sy
real(8),dimension(ny) :: cifi6y,cisi6y,ciwi6y
real(8),dimension(ny) :: cify6,cisy6,ciwy6
integer :: i,j,k

do k=1,nz
do j=1,nym
do i=1,nx
   di(i,k,j)=uy(i,j,k)
enddo
enddo
enddo
 
call interiy16(ry,di,ty,sy,cifi6y,cisi6y,ciwi6y,&
     cify6,cisy6,ciwy6,nx,nym,ny,nz,npaire)

do k=1,nz
do j=1,ny
do i=1,nx
   ty(i,j,k)=ry(i,k,j)
enddo
enddo
enddo

return 
end subroutine interiy6

!********************************************************************
!
subroutine interiy16(ty,uy,ry,sy,cifi6y,cisi6y,ciwi6y,cify6,cisy6,ciwy6,nx,nym,ny,nz,npaire) 
!
!********************************************************************
  
USE param
USE derivY

implicit none

integer :: nx,ny,nym,nz,npaire
real(8),dimension(nx,nz,ny) :: ty
real(8),dimension(nx,nz,nym) :: uy
real(8),dimension(nx,nz,ny) :: ry
real(8),dimension(nx,nz) :: sy
real(8),dimension(ny) :: cifi6y,cisi6y,ciwi6y
real(8),dimension(ny) :: cify6,cisy6,ciwy6
integer :: i,j,k
    
if (ncly==0) then
   do k=1,nz
   do i=1,nx
      ty(i,k,1)=aiciy6*(uy(i,k,1)+uy(i,k,ny))&
           +biciy6*(uy(i,k,2)+uy(i,k,ny-1))&
           +ciciy6*(uy(i,k,3)+uy(i,k,ny-2))
      ry(i,k,1)=-1.
      ty(i,k,2)=aiciy6*(uy(i,k,2)+uy(i,k,1))&
           +biciy6*(uy(i,k,3)+uy(i,k,ny))&
           +ciciy6*(uy(i,k,4)+uy(i,k,ny-1))
      ry(i,k,2)=0.
      ty(i,k,3)=aiciy6*(uy(i,k,3)+uy(i,k,2))&
           +biciy6*(uy(i,k,4)+uy(i,k,1))&
           +ciciy6*(uy(i,k,5)+uy(i,k,ny))
      ry(i,k,3)=0.
   enddo
   enddo
   do j=4,ny-2
   do k=1,nz
   do i=1,nx
      ty(i,k,j)=aiciy6*(uy(i,k,j)+uy(i,k,j-1))&
           +biciy6*(uy(i,k,j+1)+uy(i,k,j-2))&
           +ciciy6*(uy(i,k,j+2)+uy(i,k,j-3))
      ry(i,k,j)=0.
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,k,ny-1)=aiciy6*(uy(i,k,ny-1)+uy(i,k,ny-2))&
           +biciy6*(uy(i,k,ny)+uy(i,k,ny-3))&
           +ciciy6*(uy(i,k,1)+uy(i,k,ny-4))
      ry(i,k,ny-1)=0.
      ty(i,k,ny)=aiciy6*(uy(i,k,ny)+uy(i,k,ny-1))&
           +biciy6*(uy(i,k,1)+uy(i,k,ny-2))&
           +ciciy6*(uy(i,k,2)+uy(i,k,ny-3))
      ry(i,k,ny)=ailcaiy6
   enddo
   enddo
   do j=2,ny
   do k=1,nz
   do i=1,nx
      ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*cisy6(j)
      ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*cisy6(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,k,ny)=ty(i,k,ny)*ciwy6(ny)
      ry(i,k,ny)=ry(i,k,ny)*ciwy6(ny)
   enddo
   enddo
   do j=ny-1,1,-1
   do k=1,nz
   do i=1,nx
      ty(i,k,j)=(ty(i,k,j)-cify6(j)*ty(i,k,j+1))*ciwy6(j)
      ry(i,k,j)=(ry(i,k,j)-cify6(j)*ry(i,k,j+1))*ciwy6(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      sy(i,k)=(ty(i,k,1)-ailcaiy6*ty(i,k,ny))/&
           (1.+ry(i,k,1)-ailcaiy6*ry(i,k,ny))
   enddo
   enddo
   do j=1,ny
   do k=1,nz
   do i=1,nx
      ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j)
   enddo
   enddo
   enddo   
endif
if ((ncly==1).or.(ncly==2)) then
   if (npaire==1) then
      do k=1,nz
      do i=1,nx
         ty(i,k,1)=aiciy6*(uy(i,k,1)+uy(i,k,1))&
              +biciy6*(uy(i,k,2)+uy(i,k,2))&
              +ciciy6*(uy(i,k,3)+uy(i,k,3))
         ty(i,k,2)=aiciy6*(uy(i,k,2)+uy(i,k,1))&
              +biciy6*(uy(i,k,3)+uy(i,k,1))&
              +ciciy6*(uy(i,k,4)+uy(i,k,2))
         ty(i,k,3)=aiciy6*(uy(i,k,3)+uy(i,k,2))&
              +biciy6*(uy(i,k,4)+uy(i,k,1))&
              +ciciy6*(uy(i,k,5)+uy(i,k,1))
      enddo
      enddo
      do j=4,ny-3
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=aiciy6*(uy(i,k,j)+uy(i,k,j-1))&
              +biciy6*(uy(i,k,j+1)+uy(i,k,j-2))&
              +ciciy6*(uy(i,k,j+2)+uy(i,k,j-3))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny-2)=aiciy6*(uy(i,k,ny-2)+uy(i,k,ny-3))&
              +biciy6*(uy(i,k,ny-1)+uy(i,k,ny-4))&
              +ciciy6*(uy(i,k,ny-1)+uy(i,k,ny-5))
         ty(i,k,ny-1)=aiciy6*(uy(i,k,ny-1)+uy(i,k,ny-2))&
              +biciy6*(uy(i,k,ny-1)+uy(i,k,ny-3))&
              +ciciy6*(uy(i,k,ny-2)+uy(i,k,ny-4))
         ty(i,k,ny)=aiciy6*(uy(i,k,ny-1)+uy(i,k,ny-1))&
              +biciy6*(uy(i,k,ny-2)+uy(i,k,ny-2))&
              +ciciy6*(uy(i,k,ny-3)+uy(i,k,ny-3))
      enddo
      enddo
      do j=2,ny
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*cisi6y(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny)=ty(i,k,ny)*ciwi6y(ny)
      enddo
      enddo
      do j=ny-1,1,-1
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=(ty(i,k,j)-cifi6y(j)*ty(i,k,j+1))*ciwi6y(j)
      enddo
      enddo
      enddo
   endif
endif

return  
end subroutine interiy16 

!********************************************************************
!
subroutine deciy6(ty,uy,ry,di,sy,cfi6y,csi6y,cwi6y,cfy6,csy6,cwy6,ppy,nx,nym,ny,nz,npaire) 
!
!********************************************************************
  
USE param
USE derivY

implicit none

integer :: nx,ny,nym,nz,npaire
real(8),dimension(nx,ny,nz) :: ty
real(8),dimension(nx,nym,nz) :: uy
real(8),dimension(nx,nz,ny) :: ry
real(8),dimension(nx,nz,nym) :: di
real(8),dimension(nx,nz) ::sy
real(8),dimension(ny) :: cfi6y,csi6y,cwi6y
real(8),dimension(ny) :: cfy6,csy6,cwy6,ppy
integer :: i,j,k

do k=1,nz
do j=1,nym
do i=1,nx
   di(i,k,j)=uy(i,j,k)
enddo
enddo
enddo
 
call deciy16(ry,di,ty,sy,cfi6y,csi6y,cwi6y,&
     cfy6,csy6,cwy6,nx,nym,ny,nz,npaire)

if (istret.ne.0) then   
do k=1,nz 
do j=1,ny
do i=1,nx 
   ty(i,j,k)=ry(i,k,j)*ppy(j) 
enddo
enddo
enddo
endif
if (istret.eq.0) then   
do k=1,nz 
do j=1,ny
do i=1,nx 
   ty(i,j,k)=ry(i,k,j)
enddo
enddo
enddo
endif

return
end subroutine deciy6

!********************************************************************
!
subroutine deciy16(ty,uy,ry,sy,cfi6y,csi6y,cwi6y,cfy6,csy6,cwy6,nx,nym,ny,nz,npaire) 
!
!********************************************************************
  
USE param
USE derivY 

implicit none

integer :: nx,ny,nym,nz,npaire
real(8),dimension(nx,nz,ny) :: ty
real(8),dimension(nx,nz,nym) :: uy
real(8),dimension(nx,nz,ny) :: ry
real(8),dimension(nx,nz) :: sy
real(8),dimension(ny) :: cfi6y,csi6y,cwi6y
real(8),dimension(nym) :: cfy6,csy6,cwy6
integer :: i,j,k

if (ncly==0) then
   do k=1,nz 
   do i=1,nx 
      ty(i,k,1)=aciy6*(uy(i,k,1)-uy(i,k,ny))&
           +bciy6*(uy(i,k,2)-uy(i,k,ny-1)) 
      ry(i,k,1)=-1.
      ty(i,k,2)=aciy6*(uy(i,k,2)-uy(i,k,1))&
           +bciy6*(uy(i,k,3)-uy(i,k,ny)) 
      ry(i,k,2)=0.
   enddo
   enddo 
   do j=3,ny-2 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=aciy6*(uy(i,k,j)-uy(i,k,j-1))&
           +bciy6*(uy(i,k,j+1)-uy(i,k,j-2)) 
      ry(i,k,j)=0. 
   enddo
   enddo
   enddo
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny-1)=aciy6*(uy(i,k,ny-1)-uy(i,k,ny-2))&
           +bciy6*(uy(i,k,ny)-uy(i,k,ny-3)) 
      ry(i,k,ny-1)=0. 
      ty(i,k,ny)=aciy6*(uy(i,k,ny)-uy(i,k,ny-1))&
           +bciy6*(uy(i,k,1)-uy(i,k,ny-2)) 
      ry(i,k,ny)=alcaiy6
   enddo
   enddo
   do j=2,ny 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*csy6(j) 
      ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*csy6(j) 
   enddo
   enddo
   enddo
   do k=1,nz 
   do i=1,nx 
      ty(i,k,ny)=ty(i,k,ny)*cwy6(ny) 
      ry(i,k,ny)=ry(i,k,ny)*cwy6(ny) 
   enddo
   enddo
   do j=ny-1,1,-1 
   do k=1,nz 
   do i=1,nx 
      ty(i,k,j)=(ty(i,k,j)-cfy6(j)*ty(i,k,j+1))*cwy6(j) 
      ry(i,k,j)=(ry(i,k,j)-cfy6(j)*ry(i,k,j+1))*cwy6(j) 
   enddo
   enddo
   enddo
   do k=1,nz 
   do i=1,nx 
      sy(i,k)=(ty(i,k,1)-alcaiy6*ty(i,k,ny))&
           /(1.+ry(i,k,1)-alcaiy6*ry(i,k,ny)) 
   enddo
   enddo
   do j=1,ny 
   do k=1,nz 
   do i=1,nx
      ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j)
   enddo
   enddo
   enddo
endif  
if ((ncly==1).or.(ncly==2)) then
   if (npaire==1) then
      do k=1,nz
      do i=1,nx
         ty(i,k,1)=0.            
         ty(i,k,2)=aciy6*(uy(i,k,2)-uy(i,k,1))&
              +bciy6*(uy(i,k,3)-uy(i,k,1))             
      enddo
      enddo
      do j=3,ny-2
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=aciy6*(uy(i,k,j)-uy(i,k,j-1))&
              +bciy6*(uy(i,k,j+1)-uy(i,k,j-2))              
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny-1)=aciy6*(uy(i,k,ny-1)-uy(i,k,ny-2))&
              +bciy6*(uy(i,k,ny-1)-uy(i,k,ny-3))         
         ty(i,k,ny)=0.          
      enddo
      enddo
      do j=2,ny
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*csi6y(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny)=ty(i,k,ny)*cwi6y(ny)
      enddo
      enddo
      do j=ny-1,1,-1
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=(ty(i,k,j)-cfi6y(j)*ty(i,k,j+1))*cwi6y(j)
      enddo
      enddo
      enddo
   endif
endif

return  
end subroutine deciy16

!********************************************************************
!
subroutine decz6(tz,uz,rz,sz,cfz6,csz6,cwz6,nx,ny,nz,nzm,npaire) 
!
!********************************************************************
  
USE param
USE derivZ 

implicit none

integer :: nx,ny,nz,nzm,npaire
real(8),dimension(nx,ny,nz) :: tz
real(8),dimension(nx,ny,nz) :: uz
real(8),dimension(nx,ny,nz) :: rz
real(8),dimension(nx,ny) :: sz
real(8),dimension(nz) :: cfz6,csz6,cwz6
integer :: i,j,k

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=aciz6*(uz(i,j,2)-uz(i,j,1))&
           +bciz6*(uz(i,j,3)-uz(i,j,nz))
      rz(i,j,1)=-1.
      tz(i,j,2)=aciz6*(uz(i,j,3)-uz(i,j,2))&
           +bciz6*(uz(i,j,4)-uz(i,j,1))
      rz(i,j,2)=0.
   enddo
   enddo
   do k=3,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=aciz6*(uz(i,j,k+1)-uz(i,j,k))&
           +bciz6*(uz(i,j,k+2)-uz(i,j,k-1))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)=aciz6*(uz(i,j,nz)-uz(i,j,nz-1))&
           +bciz6*(uz(i,j,1)-uz(i,j,nz-2))
      rz(i,j,nz-1)=0.
      tz(i,j,nz)=aciz6*(uz(i,j,1)-uz(i,j,nz))&
           +bciz6*(uz(i,j,2)-uz(i,j,nz-1))
      rz(i  ,j,nz)=alcaiz6
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*csz6(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*csz6(k)
   enddo
   enddo
   enddo
   do i=1,nx
   do j=1,ny
      tz(i,j,nz)=tz(i,j,nz)*cwz6(nz)
      rz(i,j,nz)=rz(i,j,nz)*cwz6(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-cfz6(k)*tz(i,j,k+1))*cwz6(k)
      rz(i,j,k)=(rz(i,j,k)-cfz6(k)*rz(i,j,k+1))*cwz6(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(tz(i,j,1)-alcaiz6*tz(i,j,nz))/&
           (1.+rz(i,j,1)-alcaiz6*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
   enddo
   enddo
   enddo
endif
if ((nclz==1).or.(nclz==2)) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=aciz6*(uz(i,j,2)-uz(i,j,1))&
              +bciz6*(uz(i,j,3)-uz(i,j,2))             
         tz(i,j,2)=aciz6*(uz(i,j,3)-uz(i,j,2))&
              +bciz6*(uz(i,j,4)-uz(i,j,1))             
      enddo
      enddo
      do k=3,nzm-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=aciz6*(uz(i,j,k+1)-uz(i,j,k))&
              +bciz6*(uz(i,j,k+2)-uz(i,j,k-1))              
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm-1)=aciz6*(uz(i,j,nzm)-uz(i,j,nzm-1))&
              +bciz6*(uz(nz,j,k)-uz(nzm-2,j,k))         
         tz(i,j,nzm)=aciz6*(uz(i,j,nz)-uz(i,j,nzm))&
              +bciz6*(uz(i,j,nzm)-uz(i,j,nzm-1))            
      enddo
      enddo
      do k=2,nzm
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*csz6(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm)=tz(i,j,nzm)*cwz6(nzm)
      enddo
      enddo
      do k=nzm-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-cfz6(k)*tz(i,j,k+1))*cwz6(k)
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=aciz6*(uz(i,j,2)-uz(i,j,1))&
              +bciz6*(uz(i,j,3)-2.*uz(i,j,1)+uz(i,j,2)) 
         tz(i,j,2)=aciz6*(uz(i,j,3)-uz(i,j,2))&
              +bciz6*(uz(i,j,4)-uz(i,j,1))
      enddo
      enddo
      do k=3,nzm-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=aciz6*(uz(i,j,k+1)-uz(i,j,k))&
              +bciz6*(uz(i,j,k+2)-uz(i,j,k-1))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm-1)=aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2))&
              +bciz6*(uz(i,j,nz)-uz(i,j,nz-3)) 
         tz(i,j,nzm)=aciz6*(uz(i,j,nz)-uz(i,j,nz-1))&
              +bciz6*(2.*uz(i,j,nz)-uz(i,j,nz-1)-uz(i,j,nz-2)) 
      enddo
      enddo
      do k=2,nzm
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*csz6(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm)=tz(i,j,nzm)*cwz6(nzm)
      enddo
      enddo
      do k=nzm-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-cfz6(k)*tz(i,j,k+1))*cwz6(k)
      enddo
      enddo
      enddo
   endif
endif

return     
end subroutine decz6

!********************************************************************
!
subroutine interz6(tz,uz,rz,sz,cifz6,cisz6,ciwz6,nx,ny,nz,nzm,npaire) 
!
!********************************************************************
  
USE param 
USE derivZ 

implicit none

integer :: nx,ny,nz,nzm,npaire
real(8),dimension(nx,ny,nz) :: tz
real(8),dimension(nx,ny,nz) :: uz,rz
real(8),dimension(nx,ny) :: sz
real(8),dimension(nz) :: cifz6,cisz6,ciwz6
integer :: i,j,k

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=aiciz6*(uz(i,j,2)+uz(i,j,1))&
           +biciz6*(uz(i,j,3)+uz(i,j,nz))&
           +ciciz6*(uz(i,j,4)+uz(i,j,nz-1))
      rz(i,j,1)=-1.
      tz(i,j,2)=aiciz6*(uz(i,j,3)+uz(i,j,2))&
           +biciz6*(uz(i,j,4)+uz(i,j,1))&
           +ciciz6*(uz(i,j,5)+uz(i,j,nz))
      rz(i,j,2)=0.
   enddo
   enddo
   do k=3,nz-3
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=aiciz6*(uz(i,j,k+1)+uz(i,j,k))&
           +biciz6*(uz(i,j,k+2)+uz(i,j,k-1))&
           +ciciz6*(uz(i,j,k+3)+uz(i,j,k-2))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-2)=aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2))&
           +biciz6*(uz(i,j,nz)+uz(i,j,nz-3))&
           +ciciz6*(uz(i,j,1)+uz(i,j,nz-4))
      rz(i,j,nz-2)=0.
      tz(i,j,nz-1)=aiciz6*(uz(i,j,nz)+uz(i,j,nz-1))&
           +biciz6*(uz(i,j,1)+uz(i,j,nz-2))&
           +ciciz6*(uz(i,j,2)+uz(i,j,nz-3))
      rz(i,j,nz-1)=0.
      tz(i,j,nz)=aiciz6*(uz(i,j,1)+uz(i,j,nz))&
           +biciz6*(uz(i,j,2)+uz(i,j,nz-1))&
           +ciciz6*(uz(i,j,3)+uz(i,j,nz-2))
      rz(i  ,j,nz)=ailcaiz6
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*cisz6(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*cisz6(k)
   enddo
   enddo
   enddo
   do i=1,nx
   do j=1,ny
      tz(i,j,nz)=tz(i,j,nz)*ciwz6(nz)
      rz(i,j,nz)=rz(i,j,nz)*ciwz6(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-cifz6(k)*tz(i,j,k+1))*ciwz6(k)
      rz(i,j,k)=(rz(i,j,k)-cifz6(k)*rz(i,j,k+1))*ciwz6(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(tz(i,j,1)-ailcaiz6*tz(i,j,nz))/&
           (1.+rz(i,j,1)-ailcaiz6*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
   enddo
   enddo
   enddo
endif
if ((nclz==1).or.(nclz==2)) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=aiciz6*(uz(i,j,2)+uz(i,j,1))&
              +biciz6*(uz(i,j,3)+uz(i,j,2))&
              +ciciz6*(uz(i,j,4)+uz(i,j,3))             
         tz(i,j,2)=aiciz6*(uz(i,j,3)+uz(i,j,2))&
              +biciz6*(uz(i,j,4)+uz(i,j,1))&
              +ciciz6*(uz(i,j,5)+uz(i,j,2))             
      enddo
      enddo
      do k=3,nzm-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=aiciz6*(uz(i,j,k+1)+uz(i,j,k))&
              +biciz6*(uz(i,j,k+2)+uz(i,j,k-1))&
              +ciciz6*(uz(i,j,k+3)+uz(i,j,k-2))              
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm-1)=aiciz6*(uz(i,j,nzm)+uz(i,j,nzm-1))&
              +biciz6*(uz(i,j,nz)+uz(i,j,nzm-2))&
              +ciciz6*(uz(i,j,nzm)+uz(i,j,nzm-3))
         tz(i,j,nzm)=aiciz6*(uz(i,j,nz)+uz(i,j,nzm))&
              +biciz6*(uz(i,j,nzm)+uz(i,j,nzm-1))&
              +ciciz6*(uz(i,j,nzm-1)+uz(i,j,nzm-2))            
      enddo
      enddo
      do k=2,nzm
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*cisz6(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm)=tz(i,j,nzm)*ciwz6(nzm)
      enddo
      enddo
      do k=nzm-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-cifz6(k)*tz(i,j,k+1))*ciwz6(k)
      enddo
      enddo
      enddo
   endif
endif

return     
end subroutine interz6 

!********************************************************************
!
subroutine deciz6(tz,uz,rz,sz,cfiz6,csiz6,cwiz6,cfz6,csz6,cwz6,nx,ny,nzm,nz,npaire) 
!
!********************************************************************
  
USE param 
USE derivZ 

implicit none

integer :: nx,nzm,ny,nz,npaire
real(8),dimension(nx,ny,nz) :: tz
real(8),dimension(nx,ny,nz) :: uz,rz
real(8),dimension(nx,ny) :: sz
real(8),dimension(nz) :: cfiz6,csiz6,cwiz6
real(8),dimension(nz) :: cfz6,csz6,cwz6
integer :: i,j,k

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=aciz6*(uz(i,j,1)-uz(i,j,nz))&
           +bciz6*(uz(i,j,2)-uz(i,j,nz-1))
      rz(i,j,1)=-1.
      tz(i,j,2)=aciz6*(uz(i,j,2)-uz(i,j,1))&
           +bciz6*(uz(i,j,3)-uz(i,j,nz))
      rz(i,j,2)=0.
   enddo
   enddo
   do k=3,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=aciz6*(uz(i,j,k)-uz(i,j,k-1))&
           +bciz6*(uz(i,j,k+1)-uz(i,j,k-2))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)=aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2))&
           +bciz6*(uz(i,j,nz)-uz(i,j,nz-3))
      rz(i,j,nz-1)=0.
      tz(i,j,nz)=aciz6*(uz(i,j,nz)-uz(i,j,nz-1))&
           +bciz6*(uz(i,j,1)-uz(i,j,nz-2))
      rz(i,j,nz)=alcaiz6
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*csz6(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*csz6(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*cwz6(nz)
      rz(i,j,nz)=rz(i,j,nz)*cwz6(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-cfz6(k)*tz(i,j,k+1))*cwz6(k)
      rz(i,j,k)=(rz(i,j,k)-cfz6(k)*rz(i,j,k+1))*cwz6(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(tz(i,j,1)-alcaiz6*tz(i,j,nz))/&
           (1.+rz(i,j,1)-alcaiz6*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
   enddo
   enddo
   enddo
endif

if ((nclz==1).or.(nclz==2)) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=0.
         tz(i,j,2)=aciz6*(uz(i,j,2)-uz(i,j,1))&
              +bciz6*(uz(i,j,3)-uz(i,j,1))  
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=aciz6*(uz(i,j,k)-uz(i,j,k-1))&
              +bciz6*(uz(i,j,k+1)-uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2))&
              +bciz6*(uz(i,j,nz-1)-uz(i,j,nz-3))
         tz(i,j,nz)=0.
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*csiz6(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*cwiz6(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-cfiz6(k)*tz(i,j,k+1))*cwiz6(k)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine deciz6 

!********************************************************************
!
subroutine interiz6(tz,uz,rz,sz,cifiz6,cisiz6,ciwiz6,cifz6,cisz6,ciwz6,nx,ny,nzm,nz,npaire) 
!
!********************************************************************
  
USE param
USE derivZ 

implicit none

integer :: nx,ny,nz,nzm,npaire
real(8),dimension(nx,ny,nz) :: tz
real(8),dimension(nx,ny,nz) :: uz
real(8),dimension(nx,ny,nz) :: rz
real(8),dimension(nx,ny) :: sz
real(8),dimension(nz) :: cifiz6,cisiz6,ciwiz6
real(8),dimension(nz) :: cifz6,cisz6,ciwz6
integer :: i,j,k

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=aiciz6*(uz(i,j,1)+uz(i,j,nz))&
           +biciz6*(uz(i,j,2)+uz(i,j,nz-1))&
           +ciciz6*(uz(i,j,3)+uz(i,j,nz-2))
      rz(i,j,1)=-1.
      tz(i,j,2)=aiciz6*(uz(i,j,2)+uz(i,j,1))&
           +biciz6*(uz(i,j,3)+uz(i,j,nz))&
           +ciciz6*(uz(i,j,4)+uz(i,j,nz-1))
      rz(i,j,2)=0.
     tz(i,j,3)=aiciz6*(uz(i,j,3)+uz(i,j,2))&
           +biciz6*(uz(i,j,4)+uz(i,j,1))&
           +ciciz6*(uz(i,j,5)+uz(i,j,nz))
      rz(i,j,3)=0.
   enddo
   enddo
   do k=4,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=aiciz6*(uz(i,j,k)+uz(i,j,k-1))&
           +biciz6*(uz(i,j,k+1)+uz(i,j,k-2))&
           +ciciz6*(uz(i,j,k+2)+uz(i,j,k-3))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)=aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2))&
           +biciz6*(uz(i,j,nz)+uz(i,j,nz-3))&
           +ciciz6*(uz(i,j,1)+uz(i,j,nz-4))
      rz(i,j,nz-1)=0.
      tz(i,j,nz)=aiciz6*(uz(i,j,nz)+uz(i,j,nz-1))&
           +biciz6*(uz(i,j,1)+uz(i,j,nz-2))&
           +ciciz6*(uz(i,j,2)+uz(i,j,nz-3))
      rz(i,j,nz)=ailcaiz6
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*cisz6(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*cisz6(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*ciwz6(nz)
      rz(i,j,nz)=rz(i,j,nz)*ciwz6(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-cifz6(k)*tz(i,j,k+1))*ciwz6(k)
      rz(i,j,k)=(rz(i,j,k)-cifz6(k)*rz(i,j,k+1))*ciwz6(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(tz(i,j,1)-ailcaiz6*tz(i,j,nz))/&
           (1.+rz(i,j,1)-ailcaiz6*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
   enddo
   enddo
   enddo
endif

if ((nclz==1).or.(nclz==2)) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=aiciz6*(uz(i,j,1)+uz(i,j,1))&
              +biciz6*(uz(i,j,2)+uz(i,j,2))&
              +ciciz6*(uz(i,j,3)+uz(i,j,3))     
         tz(i,j,2)=aiciz6*(uz(i,j,2)+uz(i,j,1))&
              +biciz6*(uz(i,j,3)+uz(i,j,1))&
              +ciciz6*(uz(i,j,4)+uz(i,j,2))
         tz(i,j,3)=aiciz6*(uz(i,j,3)+uz(i,j,2))&
              +biciz6*(uz(i,j,4)+uz(i,j,1))&
              +ciciz6*(uz(i,j,5)+uz(i,j,1))
      enddo
      enddo         
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=aiciz6*(uz(i,j,k)+uz(i,j,k-1))&
              +biciz6*(uz(i,j,k+1)+uz(i,j,k-2))&
              +ciciz6*(uz(i,j,k+2)+uz(i,j,k-3))
      enddo
      enddo
      enddo     
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-2)=aiciz6*(uz(i,j,nz-2)+uz(i,j,nz-3))&
              +biciz6*(uz(i,j,nz-1)+uz(i,j,nz-4))&
              +ciciz6*(uz(i,j,nz-1)+uz(i,j,nz-5))
         tz(i,j,nz-1)=aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2))&
              +biciz6*(uz(i,j,nz-1)+uz(i,j,nz-3))&
              +ciciz6*(uz(i,j,nz-2)+uz(i,j,nz-4))
         tz(i,j,nz)=aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-1))&
              +biciz6*(uz(i,j,nz-2)+uz(i,j,nz-2))&
              +ciciz6*(uz(i,j,nz-3)+uz(i,j,nz-3))
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*cisiz6(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*ciwiz6(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-cifiz6(k)*tz(i,j,k+1))*ciwiz6(k)
      enddo
      enddo
      enddo 
   endif
endif

return
end subroutine interiz6
end module derivevitesse_m
