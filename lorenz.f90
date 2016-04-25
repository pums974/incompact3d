module lorenz_m
implicit none
contains
!***********************************************************
!
subroutine lorenz(ux,uy,uz,fx,fy,fz)
!
!************************************************************
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,fx,fy,fz
integer :: i,j,k

do k=1,nz
do j=1,ny
do i=1,nx
   ux(i,j,k)=ux(i,j,k)+fx(i,j,k)
   uy(i,j,k)=uy(i,j,k)+fy(i,j,k)
   uz(i,j,k)=uz(i,j,k)+fz(i,j,k)
enddo
enddo
enddo

return
end subroutine lorenz

!***********************************************************
!
subroutine initmagnet(fx,fy,fz)
!
!************************************************************

USE magnets
USE param
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: fx,fy,fz
integer :: i

open(201,file="input.txt")
read(201,*) nMagnets
read(201,901) jX,jY,jZ
read(201,*) nBlocks
 
do i=1,nMagnets
   read(201,901) xMag(i),yMag(i),zMag(i),angle(i),Br(i),a(i),b(i),thickness(i)
   xMag(i)=(xMag(i)*0.001/lref)+0.5*xlx
   yMag(i)=yMag(i)*0.001/lref
   zMag(i)=(zMag(i)*0.001/lref)+0.5*zlz
   a(i)=a(i)*0.001/lref
   b(i)=b(i)*0.001/lref
   thickness(i)=thickness(i)*0.001/lref
enddo
call storeForcing(fx,fy,fz)
 
if (nBlocks>1) then
   read(201,*) duration
   read(201,*) nextTime
   duration=duration*uref/lref
   nextTime=nextTime*uref/lref
   do i=1,nMagnets
      read(201,901) xMag(i),yMag(i),zMag(i),angle(i),Br(i),a(i),b(i),thickness(i)
   end do
   call storeForcing(fx,fy,fz)
endif
901 FORMAT(8(F8.3,1x))

open(20,file='forcage.dat',form='unformatted')
read(20) fx,fy,fz
close(20)
 
return
end subroutine initmagnet

!***********************************************************
!
subroutine magnet(ux,uy,uz,fx,fy,fz)
!
!************************************************************
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,fx,fy,fz
integer :: i,j,k

do k=1,nz
do j=1,ny
do i=1,nx
   ux(i,j,k)=ux(i,j,k)+fx(i,j,k)
   uy(i,j,k)=uy(i,j,k)+fy(i,j,k)
   uz(i,j,k)=uz(i,j,k)+fz(i,j,k)
enddo
enddo
enddo

return
end subroutine magnet

!************************************************************
!
subroutine magFieldSurface(xBs,yBs,zBs,Br,a,b,Bs)
!
!************************************************************
implicit none
real(8), intent(in) :: xBs,yBs,zBs,Br,a,b  !dimensions in mm
real(8), dimension(3), intent(out) :: Bs
real(8), parameter :: Pi=3.14159265

real(8), dimension(2) :: S,T
real(8), dimension(2,2) :: R
real(8) :: sumX,sumY,sumZ
integer :: i,j

do i=1,2
S(i)=zBs-(b*0.5)*(-1)**(i-1)
T(i)=xBs-(a*0.5)*(-1)**(i-1)
end do

do j=1,2
  do i=1,2
     R(i,j)=sqrt(S(i)**2+T(j)**2+yBs**2)
  end do
end do

sumX=0.0
sumY=0.0
sumZ=0.0

do j=1,2
  do i=1,2
     sumX=sumX+log(R(i,j)-S(i))*(-1)**(i+j)
     sumY=sumY+atan(S(i)*T(j)/(yBs*R(i,j)))*(-1)**(i+j)
     sumZ=sumZ+log(R(i,j)-T(j))*(-1)**(i+j)
  end do
end do

Bs(1)=(Br/(2*Pi))*sumX
Bs(2)=(Br/(2*Pi))*sumY
Bs(3)=(Br/(2*Pi))*sumZ

return
end subroutine magFieldSurface

!***********************************************************
!
subroutine magFieldMagnet(xBs,yBs,zBs,Br,a,b,thickness,Bs)
!
!***********************************************************
implicit none
real(8), intent(in) :: xBs,yBs,zBs,Br,a,b,thickness
real(8), dimension(3), intent(out) :: Bs

real(8), dimension(3) :: BsUpperSurface,BsLowerSurface

call magFieldSurface(xBs,yBs,zBs,0.5*Br,a,b,BsUpperSurface)
call magFieldSurface(xBs,yBs+thickness,zBs,0.5*Br,a,b,BsLowerSurface)

Bs(1)=BsUpperSurface(1)-BsLowerSurface(1)
Bs(2)=BsUpperSurface(2)-BsLowerSurface(2)
Bs(3)=BsUpperSurface(3)-BsLowerSurface(3)

return
end subroutine magFieldMagnet



!***********************************************************
!
subroutine changeFrame(xGlobal,zGlobal,angle,xLocal,zLocal) 
!
!***********************************************************
implicit none
real(8), intent(in) :: xGlobal,zGlobal,angle
real(8), intent(out) :: xLocal,zLocal
real(8), parameter :: Pi=3.14159265
real(8) :: angleRad

angleRad=angle*Pi*(1.0/180)

xLocal=xGlobal*cos(angleRad)-zGlobal*sin(angleRad)
zLocal=xGlobal*sin(angleRad)+zGlobal*cos(angleRad)

return
end subroutine changeFrame

!***********************************************************
!
subroutine forcing(xfl,yfl,zfl,f)
!
!***********************************************************

USE param
USE variables
USE magnets
implicit none
real(8), intent(in) :: xfl,yfl,zfl
real(8), dimension(3) :: f

real(8), dimension(3) :: BsLocal,BsGlobal,sumBs
real(8) :: xLocal,zLocal
integer :: i

sumBs(1)=0.0
sumBs(2)=0.0
sumBs(3)=0.0

do i=1,nMagnets
   call changeFrame(xfl-xMag(i),zfl-zMag(i),angle(i),xLocal,zLocal)
   call magFieldMagnet(xLocal,yfl-yMag(i),zLocal,Br(i),a(i),b(i),thickness(i),BsLocal)
   call changeFrame(BsLocal(1),BsLocal(3),-angle(i),BsGlobal(1),BsGlobal(3))
   BsGlobal(2)=BsLocal(2)
   sumBs(1)=sumBs(1)+BsGlobal(1)
   sumBs(2)=sumBs(2)+BsGlobal(2)
   sumBs(3)=sumBs(3)+BsGlobal(3)  
end do

f(1)=jY*sumBs(3)-jZ*sumBs(2)
f(2)=-jX*sumBs(3)+jZ*sumBs(1)
f(3)=jX*sumBs(2)-jY*sumBs(1)

return
end subroutine forcing

!***********************************************************
!
subroutine storeForcing(fx,fy,fz)
!
!***********************************************************

USE variables
USE param
USE magnets

implicit none

real(8), dimension(nx,ny,nz) :: fx,fy,fz
real(8), dimension(3) :: f
integer :: i,j,k

do k=1,nz
do j=1,ny/2
do i=1,nx
   call forcing((i-1)*dx,yp(j),(k-1)*dz,f)
   fx(i,j,k)=f(1)*lref/(rho*uref**2)
   fy(i,j,k)=f(2)*lref/(rho*uref**2)
   fz(i,j,k)=f(3)*lref/(rho*uref**2)
enddo
enddo
enddo
!
do k=1,nz
do j=1,ny/2
do i=1,nx
   fx(i,ny-j+1,k)=fx(i,j,k)
   fy(i,ny-j+1,k)=-fy(i,j,k)
   fz(i,ny-j+1,k)=fz(i,j,k)
enddo
enddo
enddo
!
do k=1,nz
do i=1,nx
   fx(i,ny/2+1,k)=fx(i,ny/2,k)
   fy(i,ny/2+1,k)=0.0
   fz(i,ny/2+1,k)=fz(i,ny/2,k)
enddo
enddo

return
end subroutine storeForcing

!***********************************************************
!
subroutine interpolateForcing(time,time1,time2,fx1,fy1,fz1,fx2,fy2,fz2,fx,fy,fz)
!
!***********************************************************
USE variables
implicit none
real(8) :: time,time1,time2
real(8), dimension(nx,ny,nz) :: fx1,fy1,fz1,fx2,fy2,fz2
real(8), dimension(nx,ny,nz) :: fx,fy,fz
real(8) :: a,b
integer :: i,j,k

a=(time-time1)/(time2-time1)
b=(1.-a)

do k=1,nz
do j=1,ny
do i=1,nx
   fx(i,j,k)=b*fx1(i,j,k)+a*fx2(i,j,k)
   fy(i,j,k)=b*fy1(i,j,k)+a*fy2(i,j,k)
   fz(i,j,k)=b*fz1(i,j,k)+a*fz2(i,j,k)
enddo
enddo
enddo

return
end subroutine interpolateForcing

!***********************************************************
!
subroutine forcingWithInterpolation(timeStep,fx1,fy1,fz1,fx2,fy2,fz2,fx,fy,fz)
!
!***********************************************************
USE param
USE variables
USE magnets

implicit none

real(8) :: timeStep
real(8), dimension(nx,ny,nz) :: fx1,fy1,fz1,fx2,fy2,fz2
real(8), dimension(nx,ny,nz) :: fx,fy,fz
integer :: i,j,k

!check if we are in last block
if (mem1==nBlocks) then
   fx=fx1
   fy=fy1
   fz=fz1
   return
else
   if(t<nextTime-0.5*timeStep) then
     if(timeIni+duration+0.5*timeStep<t) then
       call interpolateForcing(t,timeIni+duration,nextTime,fx1,fy1,fz1,fx2,fy2,fz2,fx,fy,fz)
       return
     else
       fx=fx1
       fy=fy1
       fz=fz1
       return
     endif
   else
     mem1=mem1+1
     fx1=fx2
     fy1=fy2
     fz1=fz2
     if(mem1<nBlocks) then
       timeIni=nextTime
       read(201,*) duration
       read(201,*) nextTime

       if(nextTime<duration) write(21,*) "warning : negative interval of interpolation at block ",mem1

       duration=duration*uref/lref
       nextTime=nextTime*uref/lref

       do i=1,nMagnets
          read(201,FMT='(8(F8.3,1x))') xMag(i),yMag(i),zMag(i),angle(i),Br(i),a(i),b(i),thickness(i)
       enddo
       call storeForcing(fx2,fy2,fz2)
       xMag(i)=(xMag(i)*0.001/lref)+0.5*xlx
       yMag(i)=yMag(i)*0.001/lref
       zMag(i)=(zMag(i)*0.001/lref)+0.5*zlz
       a(i)=a(i)*0.001/lref
       b(i)=b(i)*0.001/lref
       thickness(i)=thickness(i)*0.001/lref
     endif
     fx=fx1
     fy=fy1
     fz=fz1
     return
   endif
endif

end subroutine forcingWithInterpolation
end module lorenz_m
