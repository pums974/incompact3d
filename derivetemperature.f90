module derivetemperature_m
implicit none
contains
!********************************************************************
!
subroutine dertx (tx,ux,rx,sx,ffx,fsx,fwx,nx,ny,nz,npaire) 
!
!********************************************************************

USE param
USE convection
USE derivX

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(8),dimension(nx,ny,nz) :: tx,ux,rx 
real(8),dimension(ny,nz):: sx
real(8),dimension(nx):: ffx,fsx,fwx 

if (ncltx==0) then 
   do k=1,nz 
   do j=1,ny 
      tx(1,j,k)=afix*(ux(2,j,k)-ux(nx,j,k))& 
           +bfix*(ux(3,j,k)-ux(nx-1,j,k)) 
      rx(1,j,k)=-1. 
      tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
           +bfix*(ux(4,j,k)-ux(nx,j,k)) 
      rx(2,j,k)=0. 
   enddo
   enddo
   do i=3,nx-2 
   do k=1,nz 
   do j=1,ny 
      tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
           +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
      rx(i,j,k)=0. 
   enddo
   enddo
   enddo
   do k=1,nz 
   do j=1,ny 
      tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
           +bfix*(ux(1,j,k)-ux(nx-3,j,k)) 
      rx(nx-1,j,k)=0. 
      tx(nx,j,k)=afix*(ux(1,j,k)-ux(nx-1,j,k))&
           +bfix*(ux(2,j,k)-ux(nx-2,j,k)) 
      rx(nx,j,k)=alfaix           
   enddo
   enddo 
   do i=2, nx 
   do k=1, nz 
   do j=1, ny 
      tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*fsx(i) 
   enddo
   enddo
   enddo
   do k=1,nz 
   do j=1,ny 
      tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      rx(nx,j,k)=rx(nx,j,k)*fwx(nx) 
   enddo
   enddo 
   do i=nx-1,1,-1 
   do k=1,nz 
   do j=1,ny 
      tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      rx(i,j,k)=(rx(i,j,k)-ffx(i)*rx(i+1,j,k))*fwx(i) 
   enddo
   enddo
   enddo
   do k=1,nz 
   do j=1,ny 
      sx(j,k)=(tx(1,j,k)-alfaix*tx(nx,j,k))&
           /(1.+rx(1,j,k)-alfaix*rx(nx,j,k)) 
   enddo
   enddo
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k) 
   enddo
   enddo
   enddo 
endif

if (ncltx==1) then 
   if (npaire==1) then !Neumann
      do k=1,nz 
      do j=1,ny 
         tx(1,j,k)=0. 
         tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
              +bfix*(ux(4,j,k)-ux(2,j,k)) 
      enddo
      enddo
      do i=3,nx-2 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
              +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
      enddo
      enddo
      enddo
      do k=1,nz 
      do j=1,ny 
         tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
              +bfix*(ux(nx-1,j,k)-ux(nx-3,j,k)) 
         tx(nx,j,k)=0. 
      enddo
      enddo 
      do i=2,nx 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      enddo
      enddo
      enddo
      do k=1,nz 
      do j=1,ny 
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      enddo
      enddo
      do i=nx-1,1,-1 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then !Dirichlet
      do k=1,nz 
      do j=1,ny 
         tx(1,j,k)=afix*(ux(2,j,k)+ux(2,j,k))&
              +bfix*(ux(3,j,k)+ux(3,j,k)) 
         tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
              +bfix*(ux(4,j,k)+ux(2,j,k)) 
      enddo
      enddo
      do i=3,nx-2 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
              +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
      enddo
      enddo
      enddo
      do k=1,nz 
      do j=1,ny 
         tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
              +bfix*((-ux(nx-1,j,k))-ux(nx-3,j,k)) 
         tx(nx,j,k)=afix*((-ux(nx-1,j,k))-ux(nx-1,j,k))&
              +bfix*((-ux(nx-2,j,k))-ux(nx-2,j,k)) 
      enddo
      enddo
      do i=2,nx 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      enddo
      enddo
      enddo
      do k=1,nz 
      do j=1,ny 
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      enddo
      enddo
      do i=nx-1,1,-1 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      enddo
      enddo
      enddo
   endif
endif

if (ncltx==2) then 
   do k=1,nz 
   do j=1,ny 
      tx(1,j,k)=af1x*ux(1,j,k)+bf1x*ux(2,j,k)+cf1x*ux(3,j,k) 
      tx(2,j,k)=af2x*(ux(3,j,k)-ux(1,j,k)) 
   enddo
   enddo 
   do i=3,nx-2 
   do k=1,nz 
   do j=1,ny 
      tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
           +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
   enddo
   enddo 
   enddo 
   do k=1,nz 
   do j=1,ny 
      tx(nx-1,j,k)=afmx*(ux(nx,j,k)-ux(nx-2,j,k)) 
      tx(nx,j,k)=(-afnx*ux(nx,j,k))-bfnx*ux(nx-1,j,k)-cfnx*ux(nx-2,j,k) 
   enddo
   enddo
   do i=2,nx 
   do k=1,nz 
   do j=1,ny 
      tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
   enddo
   enddo
   enddo
   do k=1,nz 
   do j=1,ny 
      tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
   enddo
   enddo
   do i=nx-1,1,-1 
   do k=1,nz 
   do j=1,ny 
      tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
   enddo
   enddo 
   enddo 
endif

return
end subroutine dertx 

!********************************************************************
!
subroutine derty(ty,uy,ry,di,sy,ffy,fsy,fwy,ppy,nx,ny,nz,npaire) 
!
!********************************************************************
  
USE param
USE convection
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

call derty1 (ry,di,ty,sy,ffy,fsy,fwy,nx,ny,nz,npaire) 

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
end subroutine derty

!********************************************************************
!
subroutine derty1(ty,uy,ry,sy,ffy,fsy,fwy,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE convection
USE derivY 

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(8),dimension(nx,nz,ny) :: ty,uy,ry 
real(8),dimension(nx,nz) :: sy 
real(8),dimension(ny) :: ffy,fsy,fwy 

if (nclty==0) then 
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

if (nclty==1) then 
   if (npaire==1) then 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,1)=0. 
         ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
              +bfjy*(uy(i,k,4)-uy(i,k,2)) 
      enddo
      enddo 
      do k=1,nz 
      do i=1,nx 
      do j=3,ny-2 
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
      do k=1,nz 
      do i=1,nx 
      do j=3,ny-2 
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

if (nclty==2) then 
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
end subroutine derty1

!********************************************************************
!
subroutine dertz(tz,uz,rz,sz,ffz,fsz,fwz,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE convection
USE derivZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(8),dimension(nx,ny,nz) :: tz,uz,rz
real(8),dimension(nx,ny) :: sz 
real(8),dimension(nz) :: ffz,fsz,fwz

if (ncltz==0) then
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

if (ncltz==1) then
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

if (ncltz==2) then
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
end subroutine dertz

!********************************************************************
!
subroutine dertxx(tx,ux,rx,sx,sfx,ssx,swx,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE convection
USE derivX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(8),dimension(nx,ny,nz) :: tx,ux,rx
real(8),dimension(ny,nz) :: sx
real(8), dimension(nx):: sfx,ssx,swx 

   if (ncltx==0) then
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
      enddo
      enddo
      do i=4,nx-3
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-1,j,k))&
                  +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-2,j,k))&
                  +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-3,j,k))
         rx(i,j,k)=0.
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
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
      enddo
      enddo
      do i=2,nx
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*ssx(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         rx(nx,j,k)=rx(nx,j,k)*swx(nx)
      enddo
      enddo
      do i=nx-1,1,-1
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         rx(i,j,k)=(rx(i,j,k)-sfx(i)*rx(i+1,j,k))*swx(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         sx(j,k)=(   tx(1,j,k)-alsaix*tx(nx,j,k))/&
                 (1.+rx(1,j,k)-alsaix*rx(nx,j,k))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
      enddo
      enddo
   endif
!
   if (ncltx==1) then
      if (npaire==1) then !Neumann
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
         enddo
         enddo
         do i=4,nx-3
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-1,j,k))&
                     +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-2,j,k))&
                     +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-3,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
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
         enddo
         enddo
         do i=2,nx
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         enddo
         enddo
         do i=nx-1,1,-1
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then !Dirichlet
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
         enddo
         enddo
         do i=4,nx-3
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-1,j,k))&
                     +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-2,j,k))&
                     +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-3,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
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
         enddo
         enddo
         do i=2,nx
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         enddo
         enddo
         do i=nx-1,1,-1
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
         enddo
         enddo
      endif
   endif
!
   if (ncltx==2) then
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
      enddo
      enddo
!
      do i=4,nx-3
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-1,j,k))&
                  +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-2,j,k))&
                  +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-3,j,k))
!         tx(i,j,k)=asix*(ux(i+1,j,k)+ux(i-1,j,k))&
!                  +bsix*(ux(i+2,j,k)+ux(i-2,j,k))&
!                  +csix*(ux(i+3,j,k)+ux(i-3,j,k))&
!                  -(asix*2.+bsix*2.+csix*2.)*ux(i,j,k)
      enddo
      enddo
      enddo
!
      do k=1,nz
      do j=1,ny
         tx(nx-2,j,k)=astx*(ux(nx-1,j,k)-ux(nx-2,j,k)&
                           -ux(nx-2,j,k)+ux(nx-3,j,k))&
                     +bstx*(ux(nx  ,j,k)-ux(nx-2,j,k)&
                           -ux(nx-2,j,k)+ux(nx-4,j,k))
         tx(nx-1,j,k)=asmx*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                           -ux(nx-1,j,k)+ux(nx-2,j,k))
         tx(nx  ,j,k)=asnx*ux(nx  ,j,k)+bsnx*ux(nx-1,j,k)&
                     +csnx*ux(nx-2,j,k)+dsnx*ux(nx-3,j,k)
      enddo
      enddo
      do i=2,nx
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
      enddo
      enddo
      do i=nx-1,1,-1
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
      enddo
      enddo
      enddo
   endif

return  
end subroutine dertxx

!********************************************************************
!
subroutine dertyy(ty,uy,ry,di,sy,sfy,ssy,swy,nx,ny,nz,npaire) 
!
!********************************************************************
!
USE param 
USE convection
USE derivY 

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

call dertyy1 (ry,di,ty,sy,sfy,ssy,swy,nx,ny,nz,npaire) 

do k=1,nz 
do j=1,ny 
do i=1,nx 
   ty(i,j,k)=ry(i,k,j) 
enddo
enddo
enddo

return  
end subroutine dertyy

!********************************************************************
!
subroutine dertyy1(ty,uy,ry,sy,sfy,ssy,swy,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE convection
USE derivY 

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(8),dimension(nx,nz,ny) :: ty,uy,ry
real(8),dimension(nx,nz) :: sy
real(8),dimension(ny) :: sfy,ssy,swy

   if (nclty==0) then
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
   if (nclty==1) then
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
         do k=1,nz
         do i=1,nx
         do j=4,ny-3
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
         do k=1,nz
         do i=1,nx
         do j=4,ny-3
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
   if (nclty==2) then
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
      do k=1,nz
      do i=1,nx
      do j=4,ny-3
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
end subroutine dertyy1

!********************************************************************
!
subroutine dertzz(tz,uz,rz,sz,sfz,ssz,swz,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE convection
USE derivZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(8),dimension(nx,ny,nz) :: tz,uz,rz
real(8),dimension(nx,ny) :: sz 
real(8),dimension(nz) :: sfz,ssz,swz

   if (ncltz==0) then
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
!
   if (ncltz==1) then
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
!
   if (ncltz==2) then
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
end subroutine dertzz
end module derivetemperature_m