module convection_m
  use derivetemperature_m
  use navier_m
implicit none
contains
!*******************************************************************
!
subroutine force_gravite (sy8,sy9,temp)
!
!*******************************************************************

USE param
USE IBM 
USE variables
USE convection

implicit none

real(8),dimension(nx,ny,nz) :: sy8,sy9,temp
real(8) :: xm,ym,xc,yc,r
integer :: j, i, k

if (nz.gt.1) then !gravite suivant z
   do k=1,nz
   do j=1,ny 
   do i=1,nx
      sy8(i,j,k)=1./T0*gravite*(temp(i,j,k)-T0)+sy8(i,j,k) 
   enddo
   enddo
   enddo
else
   do j=1,ny !gravite suivant y
   do i=1,nx
      sy8(i,j,1)=1./T0*gravite*(temp(i,j,1)-T0)+sy8(i,j,1)
   enddo
   enddo
endif
return
end subroutine force_gravite

!
!*******************************************************************
!
subroutine energie (ux,uy,uz,temp,gtemp)
!
!*******************************************************************

USE param
USE IBM 
USE variables
USE convection

implicit none

integer :: nxyz1,ijk,i,j,k
real(8),dimension(nx,ny,nz) :: sy1,sy2,sy3,sy4,sy5,sy6,sy7,sy8,sy9,ux,uy,uz,di1,di2
real(8),dimension(nx,ny,nz) :: temp,temps,gtemp,som

real(8) :: cfl,cfl1,xnb,xm,ym,umoy,umax,tmoy,normu,normt,tmax,Reynolds,r,xc,yc
real(8) :: reynolds_cyl
nxyz1=nx*ny*nz
! 
   call limittemp (temp)
! Udelta T

if (nz.gt.1) then
   call dertx (sy4,temp,di1,    sx,ffxp,fsxp,fwxp,    nx,ny,nz,1)
   call derty (sy5,temp,sy1,di1,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
   call dertz (sy6,temp,sy2,    sz,ffzp,fszp,fwzp,    nx,ny,nz,1)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      temps(i,j,k)=-(sy4(i,j,k)*ux(i,j,k)+sy5(i,j,k)*uy(i,j,k)+sy6(i,j,k)*uz(i,j,k))
      sy1(i,j,k)=temp(i,j,k)
   enddo
   enddo
   enddo
else
   call dertx (sy4,temp,di1,    sx,ffxp,fsxp,fwxp,    nx,ny,nz,1)
   call derty (sy5,temp,sy1,di1,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
   do j=1,ny
   do i=1,nx
      temps(i,j,1)=-(sy4(i,j,1)*ux(i,j,1)+sy5(i,j,1)*uy(i,j,1))
   enddo
   enddo
endif
!diff      
diffu=xnu/Pr
if (nz.gt.1) then
   call dertxx (sy2,temp,di1,sx,sfxp,ssxp,swxp,nx,ny,nz,1)
   if (istret.ne.0) then 
      call dertyy (sy3,temp,di1,sy5,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
      call derty (sy4,temp,di1,sy5,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         sy3(i,j,k)=sy3(i,j,k)*pp2y(j)-pp4y(j)*sy4(i,j,k)
      enddo
      enddo
      enddo
   else
      call dertyy (sy3,temp,di1,sy5,sy,sfyp,ssyp,swyp,nx,ny,nz,1) 
   endif
   call dertzz (sy4,temp,di1,sz,sfzp,sszp,swzp,nx,ny,nz,1)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      temps(i,j,k)=diffu*(sy2(i,j,k)+sy3(i,j,k)+sy4(i,j,k))+temps(i,j,k)
   enddo
   enddo
   enddo
else
   call dertxx (sy2,temp,di1,sx,sfxp,ssxp,swxp,nx,ny,nz,1)
   if (istret.ne.0) then 
      call dertyy (sy3,temp,di1,sy5,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
      call derty (sy4,temp,di1,sy5,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
      do j=1,ny
      do i=1,nx
         sy3(i,j,1)=sy3(i,j,1)*pp2y(j)-pp4y(j)*sy4(i,j,1)
      enddo
      enddo
   else
      call dertyy (sy3,temp,di1,sy5,sy,sfyp,ssyp,swyp,nx,ny,nz,1) 
   endif
   do j=1,ny
   do i=1,nx
      temps(i,j,1)=diffu*(sy2(i,j,1)+sy3(i,j,1)+sy4(i,j,1))+temps(i,j,1)
   enddo
   enddo
endif
!int temp
call intt_temperature (temp,gtemp,temps)
!FORCAGE
if (nz.gt.1) then
   do k=1,nz
   do j=1,ny 
   do i=1,nx
      xm=(i-1)*dx 
      ym=yp(j)
      if(sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey)).lt.ra)then
         temp(i,j,k)=T1
      endif
   enddo
   enddo
   enddo
else
   do j=1,ny
   do i=1,nx
      xm=(i-1)*dx
      ym=yp(j)
      if(sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey)).lt.ra)then
         temp(i,j,1)=T1
      endif
   enddo
   enddo 
endif
! Calcul du Reynolds
  umoy=0. 
  tmoy=0.
  umax=-10000000. 
  tmax=umax
  xnb=0.
if (nz.gt.1) then
   do k=1,nz
   do j=1,ny 
   do i=1,nx
      xm=(i-1)*dx
      ym=yp(j)
    if(sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey)).gt.ra)then
        normu=sqrt(ux(i,j,k)*ux(i,j,k)+uy(i,j,k)*uy(i,j,k)+uz(i,j,k)*uz(i,j,k))
        umoy=umoy+normu
        normt= sqrt(temp(i,j,k)*temp(i,j,k))
        tmoy=tmoy+normt          
        xnb=xnb+1.
        if (umax.lt.normu) umax=normu
        if (tmax.lt.normt) tmax=normt
    endif
   enddo
   enddo
  enddo
else
   do j=1,ny
   do i=1,nx
      xm=(i-1)*dx
      ym=yp(j)
      if(sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey)).gt.ra)then
          normu=sqrt(ux(i,j,1)*ux(i,j,1)+uy(i,j,1)*uy(i,j,1))
          umoy=umoy+normu
          normt= sqrt(temp(i,j,1)*temp(i,j,1))
          tmoy=tmoy+normt          
          xnb=xnb+1.
          if (umax.lt.normu) umax=normu
          if (tmax.lt.normt) tmax=normt
      endif
   enddo
   enddo
endif
  umoy=umoy/xnb
  tmoy=tmoy/xnb
  reynolds_cyl=umax*0.5/xnu
  cfl=dx*dx*reynolds_cyl/4.
! TODO
  if (umax>0.) then
     cfl1=2./(umax*umax*reynolds_cyl)
  else
      cfl1=cfl
  endif
  if (cfl.gt.cfl1) then
      cfl=cfl1
  endif
! TODO
  print *,"Umoy :",umoy
  print *,"Umax :",umax
  print *,"Tmoy :",tmoy
  print *,"Tmax :",tmax
  print *,"Reynolds basé sur D :", reynolds_cyl
  print *,"Reynolds basé sur H :", umax*3./xnu
  print *,"Rayleigh :", Rayleigh
  print *,"Gr :", Gr
  print *,"Cfl :", cfl 
  !do while (dt.ge.cfl) 
  !     dt=dt/2.
  !     print *,"dt modifié:",dt
  !enddo
return
end subroutine energie
!*******************************************************************
!
subroutine limittemp (temp)
!
!*******************************************************************

USE param
USE IBM 
USE variables
USE convection

implicit none

real(8),dimension(nx,ny,nz) :: temp
integer :: j, i, k

if (iecoule.eq.7) then ! condition isotherme 
   if (nz.gt.1) then
      do k=1,nz
      do j=1,ny 
         temp(1,j,k)=T0
         temp(nx,j,k)=T0 
      enddo
      enddo
!
      do k=1,nz
      do i=1,nx 
         temp(i,1,k)=T0
         temp(i,ny,k)=T0 
      enddo
      enddo
!
     do j=1,ny
     do i=1,nx 
         temp(i,j,1)=T0
         temp(i,j,nz)=T0 
      enddo
      enddo  
   else
      do j=1,ny
         temp(1,j,1)=T0
         temp(nx,j,1)=T0 
      enddo
      do i=1,nx
         temp(i,1,1)=T0
         temp(i,ny,1)=T0 
      enddo
   endif
endif
if (iecoule.eq.1) then ! condition entree/sortie
stop "currently in progress"
   if (nz.gt.1) then
      do k=1,nz
      do j=1,ny 
         temp(1,j,k)=T0
      enddo
      enddo
!
      do k=1,nz
      do i=1,nx 
         temp(i,1,k)=0.
      enddo
      enddo
!
     do j=1,ny
     do i=1,nx 
         temp(i,j,1)=0.
      enddo
      enddo  
   else
      do j=1,ny
         temp(1,j,1)=T0
      enddo
      do i=1,nx
         temp(i,1,1)=0.
      enddo
   endif
!
 endif
return

end subroutine limittemp
!
!*******************************************************************
!
end module convection_m
