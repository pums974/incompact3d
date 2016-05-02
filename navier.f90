module navier_m
  use derivevitesse_m
  use scalar_m
  use tools_m
implicit none
contains
!********************************************************************
!
subroutine initial(ux,uy,uz,gx,gy,gz,fx,fy,fz,phi,ppm,temp,gtemp)
!
!********************************************************************

USE param
USE IBM
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,fx,fy,fz
real(8),dimension(nx,ny,nz) :: temp,gtemp,visct,condt
real(8),dimension(nx,ny,nz,nphi) :: phi
real(8),dimension(nxm,nym,nzm) :: ppm
real(8) :: y,r,um,r1,r2,r3
integer :: k,j,i,l
integer :: l1,l2,nxyz
character*80 nfichier

nxyz=nx*ny*nz

if (ilit==0) then
   ux=0.;uy=0.;uz=0.
   call ecoule(ux,uy,uz)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      ux(i,j,k)=ux(i,j,k)+bxx1(j,k)
      uy(i,j,k)=uy(i,j,k)+bxy1(j,k)
      uz(i,j,k)=uz(i,j,k)+bxz1(j,k)
      gx(i,j,k)=ux(i,j,k)
      gy(i,j,k)=uy(i,j,k)
      gz(i,j,k)=uz(i,j,k)
   enddo
   enddo
   enddo
   if (iscalaire==1) call initscalar(phi)
   if (iconvect==1) call inittemp(temp)
endif

if (ilit.eq.1) then

open(11,file='restart',form='unformatted',status='unknown')
if (iscalaire==0) then
  if (iconvect == 0) then
   if (nz.gt.1) then
      read(11) ux,uy,uz,ppm,gx,gy,gz,dpdyx1,dpdyxn,&
              dpdzx1,dpdzxn,dpdxy1,dpdxyn,dpdzy1,dpdzyn,dpdxz1,dpdxzn,dpdyz1,dpdyzn
   else
      read(11) ux,uy,ppm,gx,gy,dpdyx1,dpdyxn,dpdxy1,dpdxyn
   endif
  else
   if (nz.gt.1) then
      read(11) ux,uy,uz,ppm,gx,gy,gz,temp,gtemp,dpdyx1,dpdyxn,&
               dpdzx1,dpdzxn,dpdxy1,dpdxyn,dpdzy1,dpdzyn,dpdxz1,dpdxzn,dpdyz1,dpdyzn
   else
      read(11) ux,uy,ppm,gx,gy,temp,gtemp,dpdyx1,dpdyxn,dpdxy1,dpdxyn
   endif
  endif
else
  if (iconvect == 0) then
   if (nz.gt.1) then
      read(11) ux,uy,uz,phi,ppm,gx,gy,gz,dpdyx1,dpdyxn,dpdzx1,dpdzxn,&
               dpdxy1,dpdxyn,dpdzy1,dpdzyn,dpdxz1,dpdxzn,dpdyz1,dpdyzn
   else
      read(11) ux,uy,phi,ppm,gx,gy,dpdyx1,dpdyxn,dpdxy1,dpdxyn
   endif
  else
   if (nz.gt.1) then
      read(11) ux,uy,uz,phi,ppm,gx,gy,gz,temp,gtemp,dpdyx1,dpdyxn,&
               dpdzx1,dpdzxn,dpdxy1,dpdxyn,dpdzy1,dpdzyn,dpdxz1,dpdxzn,dpdyz1,dpdyzn
   else
      read(11) ux,uy,phi,ppm,gx,gy,temp,gtemp,dpdyx1,dpdyxn,dpdxy1,dpdxyn
   endif
  endif
endif
close(11)

endif

return
end subroutine initial

!********************************************************************
!                
subroutine  intt (ux,uy,uz,gx,gy,gz,hx,hy,hz)
! 
!********************************************************************

USE param
USE variables

implicit none

integer :: ijk,nxyz,i,j,k
real(8),dimension(nx,ny,nz) :: ux,uy,uz,hx,hy,hz,gx,gy,gz

nxyz=nx*ny*nz

if ((nschema.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
     (nschema.eq.2.and.itr.eq.1)) then
   if (nz.gt.1) then
      do ijk=1,nxyz
         ux(ijk,1,1)=gdt(itr)*hx(ijk,1,1)+ux(ijk,1,1)
         uy(ijk,1,1)=gdt(itr)*hy(ijk,1,1)+uy(ijk,1,1) 
         uz(ijk,1,1)=gdt(itr)*hz(ijk,1,1)+uz(ijk,1,1)
         gx(ijk,1,1)=hx(ijk,1,1)
         gy(ijk,1,1)=hy(ijk,1,1)
         gz(ijk,1,1)=hz(ijk,1,1)            
      enddo
   else
!      do ijk=1,nxyz ! TODO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         ux(i,j,k)=gdt(itr)*hx(i,j,k)+ux(i,j,k)
         uy(i,j,k)=gdt(itr)*hy(i,j,k)+uy(i,j,k)   
         gx(i,j,k)=hx(i,j,k)
         gy(i,j,k)=hy(i,j,k)
      enddo
      enddo
      enddo
   endif
else
   if (nz.gt.1) then
      do ijk=1,nxyz
         ux(ijk,1,1)=adt(itr)*hx(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+ux(ijk,1,1)
         uy(ijk,1,1)=adt(itr)*hy(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+uy(ijk,1,1)   
         uz(ijk,1,1)=adt(itr)*hz(ijk,1,1)+bdt(itr)*gz(ijk,1,1)+uz(ijk,1,1)
         gx(ijk,1,1)=hx(ijk,1,1)
         gy(ijk,1,1)=hy(ijk,1,1)
         gz(ijk,1,1)=hz(ijk,1,1)            
      enddo
   else
!      do ijk=1,nxyz ! TODO
      do k=1,nz
      do j=1,ny
      do i=1,nx
         ux(i,j,k)=adt(itr)*hx(i,j,k)+bdt(itr)*gx(i,j,k)+ux(i,j,k)
         uy(i,j,k)=adt(itr)*hy(i,j,k)+bdt(itr)*gy(i,j,k)+uy(i,j,k)   
         gx(i,j,k)=hx(i,j,k)
         gy(i,j,k)=hy(i,j,k)
      enddo
      enddo
      enddo
   endif
endif

if (nschema.eq.3) then 
   if (nz.gt.1) then
      if (adt(itr)==0.) then
         do ijk=1,nxyz
            gx(ijk,1,1)=dt*hx(ijk,1,1)
            gy(ijk,1,1)=dt*hy(ijk,1,1)
            gz(ijk,1,1)=dt*hz(ijk,1,1)
         enddo
      else
         do ijk=1,nxyz
            gx(ijk,1,1)=adt(itr)*gx(ijk,1,1)+dt*hx(ijk,1,1)
            gy(ijk,1,1)=adt(itr)*gy(ijk,1,1)+dt*hy(ijk,1,1)
            gz(ijk,1,1)=adt(itr)*gz(ijk,1,1)+dt*hz(ijk,1,1)
         enddo
      endif
      do ijk=1,nxyz
         ux(ijk,1,1)=ux(ijk,1,1)+bdt(itr)*gx(ijk,1,1)
         uy(ijk,1,1)=uy(ijk,1,1)+bdt(itr)*gy(ijk,1,1)
         uz(ijk,1,1)=uz(ijk,1,1)+bdt(itr)*gz(ijk,1,1)
      enddo
   else
      if (adt(itr)==0.) then
         do ijk=1,nxyz
            gx(ijk,1,1)=dt*hx(ijk,1,1)
            gy(ijk,1,1)=dt*hy(ijk,1,1)
         enddo
      else
         do ijk=1,nxyz
            gx(ijk,1,1)=adt(itr)*gx(ijk,1,1)+dt*hx(ijk,1,1)
            gy(ijk,1,1)=adt(itr)*gy(ijk,1,1)+dt*hy(ijk,1,1)
         enddo
      endif
      do ijk=1,nxyz
         ux(ijk,1,1)=ux(ijk,1,1)+bdt(itr)*gx(ijk,1,1)
         uy(ijk,1,1)=uy(ijk,1,1)+bdt(itr)*gy(ijk,1,1)
      enddo
   endif
endif

return
end subroutine intt
!********************************************************************
!                
subroutine  intt_temperature (temp,gtemp,htemp)
! 
!********************************************************************

USE param
USE variables

implicit none

integer :: ijk,nxyz,i,j,k
real(8),dimension(nx,ny,nz) :: temp, gtemp, htemp

nxyz=nx*ny*nz

if ((nschema.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
     (nschema.eq.2.and.itr.eq.1)) then
!   do ijk=1,nxyz ! TODO
    do k=1,nz
    do j=1,ny
    do i=1,nx
      temp(i,j,k)=gdt(itr)*htemp(i,j,k)+temp(i,j,k)
      gtemp(i,j,k)=htemp(i,j,k)            
   enddo
   enddo
   enddo
else
!   do ijk=1,nxyz ! TODO
    do k=1,nz
    do j=1,ny
    do i=1,nx
      temp(i,j,k)=adt(itr)*htemp(i,j,k)+bdt(itr)*gtemp(i,j,k)+temp(i,j,k)
      gtemp(i,j,k)=htemp(i,j,k)
   enddo
   enddo
   enddo
endif

if (nschema.eq.3) then 
      if (adt(itr)==0.) then
         do ijk=1,nxyz
            gtemp(ijk,1,1)=dt*htemp(ijk,1,1)
         enddo
      else
         do ijk=1,nxyz
            gtemp(ijk,1,1)=adt(itr)*gtemp(ijk,1,1)+dt*htemp(ijk,1,1)
         enddo
      endif
      do ijk=1,nxyz
         temp(ijk,1,1)=temp(ijk,1,1)+bdt(itr)*gtemp(ijk,1,1)
      enddo
endif

return
end subroutine intt_temperature
!
!*******************************************************************
!
subroutine inittemp (temp)
!
!*******************************************************************

USE param
USE IBM 
USE variables
USE convection

implicit none

real(8),dimension(nx,ny,nz) :: temp
integer :: j, i, k
if (nz.gt.1) then
   do k=1,nz
   do j=1,ny 
   do i=1,nx
      temp(i,j,k)=T0 
   enddo
   enddo
   enddo
else
   do j=1,ny
   do i=1,nx
      temp(i,j,1)=T0
   enddo
   enddo
endif

return

end subroutine inittemp

!********************************************************************
!
subroutine corgp (ux,uy,uz,px,py,pz)
! 
!********************************************************************

USE param
USE variables
USE convection


implicit none

integer :: ijk,nxyz,i,j,k
real(8),dimension(nx,ny,nz) :: ux,uy,uz,px,py,pz
real(8) :: can,ut3,ut

nxyz=nx*ny*nz

if (nz.gt.1) then
   do ijk=1,nxyz
      uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1) 
      uz(ijk,1,1)=-pz(ijk,1,1)+uz(ijk,1,1) 
      ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
   enddo
else
!   do ijk=1,nxyz ! TODO
    do k=1,nz
    do j=1,ny
    do i=1,nx
      uy(i,j,k)=-py(i,j,k)+uy(i,j,k) 
      ux(i,j,k)=-px(i,j,k)+ux(i,j,k)
   enddo
   enddo
   enddo
endif

if (iecoule.eq.5) then
   ut3=0.
   do k=1,nz
   do i=1,nx
      ut=0.
      do j=1,ny-1
         ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-0.5*(ux(i,j+1,k)-ux(i,j,k)))
      enddo
      ut=ut/yly
      ut3=ut3+ut
   enddo
   enddo
   ut3=ut3/nz/nx

   can=-(2./3.-ut3) ! gradient de pression constant

   do k=1,nz
   do i=1,nx
   do j=2,ny-1
      ux(i,j,k)=-can+ux(i,j,k)
   enddo
   enddo
   enddo
endif
if (iconvect == 1) then
   if (iecoule == 10) then
      if (nz.gt.1) then
         do k=1,nz 
         do i=1,nx 
            uy(i,1,k)=0.
            uy(i,ny,k)=0.
         enddo
         enddo
      else
         do i=1,nx 
            uy(i,1,1)=0.
            uy(i,ny,1)=0.
         enddo
      endif
   endif
endif
return
end subroutine corgp

!*********************************************************
!
subroutine inflow (ux,uy,uz)
!  
!*********************************************************

USE param
USE IBM
USE variables

implicit none

integer  :: k,j,idum
real(8),dimension(nx,ny,nz) :: ux,uy,uz
real(8) :: r1,r2,r3,y,um

idum=-67

call ecoule(ux,uy,uz)

if (ientree.eq.1) then  
   if (nz.gt.1) then
      do k=1,nz
      do j=1,ny
         r1=rand2(idum)-0.5
         r2=rand2(idum)-0.5
         r3=rand2(idum)-0.5
         bxx1(j,k)=bxx1(j,k)+r1*bruit
         bxy1(j,k)=bxy1(j,k)+r2*bruit
         bxz1(j,k)=bxz1(j,k)+r3*bruit
      enddo
      enddo
   else
      do j=1,ny
         y=yp(j)-yly/2.
         um=exp(-0.2*y*y)
         r1=rand2(idum)-0.5
         r2=rand2(idum)-0.5
         bxx1(j,1)=bxx1(j,1)+r1*bruit
         bxy1(j,1)=bxy1(j,1)+r2*bruit
      enddo
   endif
endif

return
end subroutine inflow 

!*********************************************************
!
subroutine outflow (ux,uy,uz)
!
!*********************************************************

USE param
USE variables

implicit none

integer :: j,k,i
real(8),dimension(nx,ny,nz) :: ux,uy,uz
real(8) :: udx,udy,udz,uddx,uddy,uddz,uxmax,uxmin,vphase,cx,coef,uxmax1,uxmin1


byxn(:,:)=0.;byyn(:,:)=0.;byzn(:,:)=0.;bzxn(:,:)=0.;bzyn(:,:)=0.;bzzn(:,:)=0.;
udx=1./dx
udy=1./dy
udz=1./dz
uddx=0.5/dx
uddy=0.5/dy
uddz=0.5/dz
cx=0.5*(u1+u2)*gdt(itr)*udx

uxmax=-1.e10
uxmin=1.e10
do k=1,nz
do j=1,ny
   if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
   if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
enddo
enddo
vphase=0.5*(uxmax+uxmin)
!cx=vphase*gdt(itr)*udx
coef=gdt(itr)*dx/(dx+cx*gdt(itr))

if (iecoule.ne.9) then
   if (nz.gt.1) then
      do k=1,nz
      do j=1,ny
         bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
         bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
         bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
      enddo
      enddo
   else
      do j=1,ny
         bxxn(j,1)=ux(nx,j,1)-cx*(ux(nx,j,1)-ux(nx-1,j,1))
         bxyn(j,1)=uy(nx,j,1)-cx*(uy(nx,j,1)-uy(nx-1,j,1))
      enddo
   endif
else
   do k=2,nz-1
   do i=2,nx-1
      byxn(i,k)=ux(i,ny-1,k)-(uy(i+1,ny-1,k)-uy(i-1,ny-1,k))*udx*0.5
      byzn(i,k)=uz(i,ny-1,k)-(uz(i-1,ny-1,k)-uz(i-1,ny-1,k))*udz*0.5
      byyn(i,k)=0.
   enddo
   enddo
endif

return
end subroutine outflow 

!********************************************************************
!
subroutine gradpression(ppm,sy1,sy2,sy3,sy4,sy5,sy6,sy7,sy8,sy9,di1,di2)
!
!********************************************************************

USE param 
USE variables

implicit none

real(8),dimension(mx,my,mz) :: di2,di1
real(8),dimension(nxm,ny,nzm) :: sy4
real(8),dimension(nxm,ny,nz) :: sy7
real(8),dimension(nx,nym,nzm) :: sy8,sy9
real(8),dimension(nx,ny,nzm) :: sy6
real(8),dimension(nx,nym,nz) :: sy5
real(8),dimension(nx,ny,nz) :: sy1,sy2,sy3
real(8),dimension(nxm,nym,nzm) :: ppm
integer :: i,j,k,mxyz,ijk
real(8):: x,y,z

call interiy6(sy4,ppm,di1,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,nxm,nym,ny,nz,1)
call interi6(sy5,ppm,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,nxm,nx,nym,nz,1)
if (nz.gt.1) then
   call interi6(sy9,ppm,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,nxm,nx,nym,nz,1) 
   call interiy6(sy6,sy9,di1,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,nx,nym,ny,nz,1)
   call interiz6(sy7,sy4,di1,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,nxm,ny,nzm,nz,1)
   call interiz6(sy8,sy5,di1,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,nx,nym,nzm,nz,1)
   call deciz6(sy3,sy6,di1,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,nx,ny,nzm,nz,1)
endif
if (nz.gt.1) then
   call deci6(sy1,sy7,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,nxm,nx,ny,nz,1)
   call deciy6(sy2,sy8,di1,di2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,nx,nym,ny,nz,1)
else
   call deci6(sy1,sy4,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,nxm,nx,ny,nz,1)
   call deciy6(sy2,sy5,di1,di2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,nx,nym,ny,nz,1)
endif

!      byzn(i,k)=0.

do j=1,ny
do i=1,nx
   dpdxz1(i,j)=sy1(i,j,1)/gdt(itr)
   dpdyz1(i,j)=sy2(i,j,1)/gdt(itr)
   dpdxzn(i,j)=sy1(i,j,nz)/gdt(itr)
   dpdyzn(i,j)=sy2(i,j,nz)/gdt(itr)
enddo
enddo

do k=1,nz
do i=1,nx
   dpdxy1(i,k)=sy1(i,1,k)/gdt(itr)
   dpdxyn(i,k)=sy1(i,ny,k)/gdt(itr)
   dpdzy1(i,k)=sy3(i,1,k)/gdt(itr)
   dpdzyn(i,k)=sy3(i,ny,k)/gdt(itr)
enddo
enddo
do k=1,nz
do j=1,ny    
   dpdyxn(j,k)=sy2(nx,j,k)/gdt(itr)
   dpdyx1(j,k)=sy2(1,j,k)/gdt(itr)
   dpdzx1(j,k)=sy3(1,j,k)/gdt(itr)
   dpdzxn(j,k)=sy3(nx,j,k)/gdt(itr)
enddo
enddo

return
end subroutine gradpression

!********************************************************************
!
subroutine divergence (ppm,ux,uy,uz,sy4,sy5,sy6,di1,di2,&
                       sy7,sy8,sy9,sy10,sy11,sy12,&
                       epsi,sy1,sy2,sy3,work,table,nlock)    
! 
!********************************************************************

USE param
USE IBM
USE variables
   
implicit none

integer :: i,j,k,ijk,nxyz,l,imax,kmax
integer :: nlock
real(8),dimension(nwork) :: work,table
real(8),dimension(nx,ny,nz) :: ux,uy,uz,di1,di2,epsi,sy1,sy2,sy3
real(8),dimension(nx,nym,nzm) :: sy10
real(8),dimension(nxm,ny,nzm) ::sy11
real(8),dimension(nx,nym,nz) :: sy7
real(8),dimension(nxm,ny,nz) :: sy8,sy12
real(8),dimension(nxm,nym,nz) :: sy9
real(8),dimension(nxm,nym,nzm) :: ppm,sy4,sy5,sy6!,err
real(8) :: tmax,tmoy,tmax1,tmoy1,xl2,x,y,z
real(8),dimension(mx,my) :: tb11

nxyz=nx*ny*nz

if (ivirtuel.eq.0) epsi(:,:,:)=0.

if (nz.gt.1) then
   do i=1,nxyz
      sy1(i,1,1)=(1.-epsi(i,1,1))*ux(i,1,1)
      sy2(i,1,1)=(1.-epsi(i,1,1))*uy(i,1,1)
      sy3(i,1,1)=(1.-epsi(i,1,1))*uz(i,1,1)
   enddo
else
!   do i=1,nxyz ! TODO
   do k=1,nz
   do j=1,ny
   do i=1,nx
      sy1(i,j,k)=(1.-epsi(i,j,k))*ux(i,j,k)
      sy2(i,j,k)=(1.-epsi(i,j,k))*uy(i,j,k)
!      sy1(i,1,1)=ux(i,1,1)
!      sy2(i,1,1)=uy(i,1,1)
   enddo
   enddo
   enddo
endif



call intery6(sy7,sy1,di1,di2,sy,cifyp6,cisyp6,ciwyp6,nx,ny,nym,nz,1)
call inter6(sy8,sy2,di1,sx,cifxp6,cisxp6,ciwxp6,nx,nxm,ny,nz,1)
if (nz.gt.1) then
   call inter6(sy12,sy3,di1,sx,cifxp6,cisxp6,ciwxp6,nx,nxm,ny,nz,1)
   call intery6(sy9,sy12,di1,di2,sy,cifyp6,cisyp6,ciwyp6,nxm,ny,nym,nz,1)
   call interz6(sy10,sy7,di1,sz,cifzp6,ciszp6,ciwzp6,nx,nym,nz,nzm,1)    
   call interz6(sy11,sy8,di1,sz,cifzp6,ciszp6,ciwzp6,nxm,ny,nz,nzm,1)
   call decz6(sy6,sy9,di1,sz,cfz6,csz6,cwz6,nxm,nym,nz,nzm,0)
endif
if (nz.gt.1) then
   call decx6(sy4,sy10,di1,sx,cfx6,csx6,cwx6,nx,nxm,nym,nz,0)
   call decy6(sy5,sy11,di1,di2,sy,cfy6,csy6,cwy6,ppyi,nxm,ny,nym,nz,0)
else
   call decx6(sy4,sy7,di1,sx,cfx6,csx6,cwx6,nx,nxm,nym,nz,0)
   call decy6(sy5,sy8,di1,di2,sy,cfy6,csy6,cwy6,ppyi,nxm,ny,nym,nz,0) 
endif

if (nz.gt.1) then
   nxyz=nxm*nym*nzm
   do ijk=1,nxyz
      ppm(ijk,1,1)=sy4(ijk,1,1)+sy5(ijk,1,1)+sy6(ijk,1,1)
   enddo
   if (nlock==2) then
      xl2=ppm(1,1,1)
      do k=1,nzm
      do j=1,nym
      do i=1,nxm
         ppm(i,j,k)=ppm(i,j,k)-xl2
      enddo
      enddo
      enddo
   endif
else
   do j=1,nym
   do i=1,nxm
      ppm(i,j,1)=sy4(i,j,1)+sy5(i,j,1)
   enddo
   enddo
   if (nlock==2) then
      xl2=ppm(1,1,1)
      do j=1,nym
      do i=1,nxm
         ppm(i,j,1)=ppm(i,j,1)-xl2
!         err(i,j,1)=ppm(i,j,1)
      enddo
      enddo
   endif
endif

tmax=-1609.
tmoy=0.
do k=1,nzm
do j=1,nym
do i=1,nxm
   if (ppm(i,j,k).gt.tmax) tmax=ppm(i,j,k)     
   tmoy=tmoy+abs(ppm(i,j,k))
enddo
enddo
enddo
tmoy=tmoy/nxm/nym/nzm
if (nlock==2) then
   write(*,1001) tmax,tmoy
 1001 format('Div(u) final Max, Moy= ',2E10.3)
else
   write(*,1002) tmax,tmoy
 1002 format('Div(u*) Max, Moy     = ',2E10.3)
endif

return
end subroutine divergence

!**********************************************************************
!
subroutine ecoule (ux,uy,uz)
!
!**********************************************************************

USE param
USE IBM
USE variables

implicit none

integer  :: i,j,k,jj1,jj2 
real(8),dimension(nx,ny,nz) :: ux,uy,uz
real(8)  :: x,y,ym
real(8) :: r1,r2,r3,r
real(8) :: uh,ud,um,xxk1,xxk2,xv,bruit1

bxx1=0.;bxy1=0.;bxz1=0.
byx1=0.;byy1=0.;byz1=0.
bzx1=0.;bzy1=0.;bzz1=0. 

!IECOULE=1 --> Constant flow field
!IECOULE=2 --> Mixing layer hyperbolic tangent profile
!IECOULE=3 --> Wake flow
!IECOULE=4 --> Mixing layer with splitter plate
!IECOULE=5 --> Channel flow
!IECOULE=6 --> Taylor Green vortices
!IECOULE=7 --> Cavity flow
!IECOULE=8 --> Flat plate Boundary layer
!IECOULE=9 --> Tank 
!IECOULE=10 -> Convection-Neuman ou periodique 
if (iecoule.eq.1) then
   um=0.5*(u1+u2)
   do k=1,nz
   do j=1,ny
      bxx1(j,k)=um
      bxy1(j,k)=0.
      bxz1(j,k)=0.
   enddo
   enddo
endif

if (iecoule.eq.2) then
   uh=0.5*(u1+u2)
   ud=0.5*(u1-u2)
   do k=1,nz
   do j=1,ny
      y=yp(j)-yly/2.
      bxx1(j,k)=uh+ud*tanh(2.*y)
      bxy1(j,k)=0.
      bxz1(j,k)=0.
   enddo
   enddo
endif

if (iecoule.eq.3) then
   do k=1,nz
   do j=1,ny
      y=yp(j)-yly/2.
      bxx1(j,k)=1-exp(-log(2.)*y**2)
      bxy1(j,k)=0.
      bxz1(j,k)=0.
   enddo
   enddo
endif

if (iecoule.eq.4) then
   if (istret.ne.1) then
       print *,'ECOULEMENT IMPOSSIBLE SI ISTRET <> 1!!'
   endif
   i5=int(cex/dx)
   jj=int(cey/dy)+1
   j1=1
   do j=2,jj
      r=abs(yp(j)-cey)
      if (r.gt.ra) j1=j+1
   enddo
   j2=1
   do j=ny,jj,-1
      r=abs(yp(j)-cey)
      if (r.gt.ra) j2=j-1
   enddo

   um=0.5*(u1+u2)
   do k=1,nz
   do j=1,ny
      bxx1(j,k)=um
      bxy1(j,k)=0.
      bxz1(j,k)=0.
   enddo
   enddo
   do j=ny,jj,-1
      r=abs(yp(j)-cey)
      if (r.gt.(0.5+(11./19.))) jj2=j
   enddo
   do k=1,nz
   do j=j2,jj2
      ym=(yp(j)-yp(j2))/(yp(jj2)-yp(j2))
      ym=ym*0.772588621523 
      bxx1(j,k)=2.*ym*u1-5.*ym*ym*ym*ym*u1+6.*ym*ym*ym*ym*ym*u1-2.*ym*ym*ym*ym*ym*ym*u1
      if ((j>j2+1).and.(bxx1(j,k) < bxx1(j-1,k))) bxx1(j,k)=u1
   enddo
   do j=jj2+1,ny
      bxx1(j,k)=u1
   enddo
   do j=j1,j2
      bxx1(j,k)=0.
   enddo
   enddo
   do j=2,jj
      r=abs(yp(j)-cey)
      if (r.gt.(ra+(9.5/19.))) jj1=j+1
   enddo
   do k=1,nz
   do j=j1,jj1,-1
      ym=(yp(j)-yp(j1))/(yp(jj1)-yp(j1))
      ym=ym*0.772588621523
      bxx1(j,k)=2.*ym*u2-5.*ym*ym*ym*ym*u2+6.*ym*ym*ym*ym*ym*u2-2.*ym*ym*ym*ym*ym*ym*u2
      if ((j<j1-1).and.(bxx1(j,k) < bxx1(j+1,k))) bxx1(j,k)=u2
   enddo
   do j=1,jj1-1
      bxx1(j,k)=u2
   enddo
   enddo
endif

if (iecoule.eq.5) then
   bruit1=0.125
   if (ilit.eq.0) then
      do k=1,nz
      do j=1,ny
      do i=1,nx
         x=(i-1)*dx-xlx/2.
         y=yp(j)-yly/2.
         if (istret.eq.1.or.istret.eq.3) then
            print *,'Not possible with this refinement parameter'
            stop
         endif
         ux(i,j,k)=1.-y*y+bruit1*r1
         uy(i,j,k)=bruit1*r2
         uz(i,j,k)=bruit1*r3
      enddo
      enddo
      enddo
   endif
endif

if (iecoule.eq.6) then
   if (nz.gt.1) then 
      print *,'ECOULEMENT UNIQUEMENT 2D'
      stop
   else
   xv=1./100.
   xxk1=twopi/xlx
   xxk2=twopi/yly
   do k=1,nz
   do j=1,ny
      y=yp(j)
      do i=1,nx
         x=(i-1)*dx
         ux(i,j,k)=sin(xxk1*x)*cos(xxk2*y)*exp(-(xxk1*xxk1+xxk2*xxk2)*xv*t)
         uy(i,j,k)=-xxk1/xxk2*sin(xxk2*y)*cos(xxk1*x)*exp(-(xxk1*xxk1+xxk2*xxk2)*xv*t)
         bxx1(j,k)=0.
         bxy1(j,k)=0.
      enddo
   enddo
   enddo
   endif
endif

if (iecoule.eq.7) then
   if (nz.gt.1) then 
   do j=1,ny
   do i=1,nx
      ux(i,j,1)=0.
      uy(i,j,1)=0.
      uz(i,j,1)=0.
      ux(i,j,nz)=0.
      uy(i,j,nz)=0.             
      uz(i,j,nz)=0.             
   enddo
   enddo
   do k=1,nz
   do j=1,ny
      ux(1,j,k)=0.
      uy(1,j,k)=0.
      uz(1,j,k)=0.
      ux(nx,j,k)=0.
      uy(nx,j,k)=0.             
      uz(nx,j,k)=0.             
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ux(i,1,k)=0.
      uy(i,1,k)=0.
      uz(i,1,k)=0.
      ux(i,ny,k)=0.
      uy(i,ny,k)=0.             
      uz(i,ny,k)=0.             
   enddo
   enddo
   else
!      do k=1,nz
!      do j=1,ny
!         y=(j-1)*dy
!         do i=1,nx
!            x=(i-1)*dx
           ! ux(i,j,k)=   8.*(x*x*x*x-2.*x*x*x+x*x)*(4.*y*y*y-2.*y)
           ! uy(i,j,k)=  -8.*(4.*x*x*x-6.*x*x+2.*x)*(y*y*y*y-y*y)
           ! bxx1(j,k)=0.
           ! bxy1(j,k)=0.
           ! bxz1(j,k)=0.
!         enddo
!      enddo
!      enddo
   do j=1,ny
      ux(1,j,1)=0.
      uy(1,j,1)=0.
      ux(nx,j,1)=0.
      uy(nx,j,1)=0.             
   enddo
   do i=1,nx
      ux(i,1,1)=0.
      uy(i,1,1)=0.
      ux(i,ny,1)=0.
      uy(i,ny,1)=0.             
   enddo
   endif
endif

if (iecoule.eq.8) then
   do k=1,nz
   do j=1,ny
      ym=yp(j)*2.
      bxx1(j,k)=2.*ym*u1-5.*ym*ym*ym*ym*u1+6.*ym*ym*ym*ym*ym*u1-2.*ym*ym*ym*ym*ym*ym*u1
      if (ym.gt.1.) bxx1(j,k)=u1
      bxy1(j,k)=0.
      bxz1(j,k)=0.
      byx1(j,k)=1.
      byy1(j,k)=0.
      byz1(j,k)=0.
   enddo
   enddo
endif

if (iecoule.eq.9) then
   bxx1=0.
   bxy1=0.
   bxz1=0.
   bzx1=0.
   bzy1=0.
   bzz1=0.
endif

if (iecoule.eq.10) then
   bxx1=0.
   bxy1=0.
   bxz1=0.
   bzx1=0.
   bzy1=0.
   bzz1=0.
endif

return
end subroutine ecoule

!********************************************************************
!
subroutine cavite (uy,nx,ny,nz)
!
!********************************************************************

USE param

implicit none

integer :: nx,ny,nz
real(8),dimension(nx,ny,nz) :: cav,uy
integer :: nxyz,j,i,k
real(8) :: x,y,xy1,xy2,xy3,xy4,xy5,xy6,xy7,xy8

nxyz=nx*ny*nz

do k=1,nz
do j=1,ny
   y=(j-1.)*dy
   do i=1,nx
      x=(i-1.)*dx
      xy1=x*x*x*x-2.*x*x*x+x*x
      xy2=y*y*y*y-y*y
      xy3=0.2*x*x*x*x*x-0.5*x*x*x*x+(1./3.)*x*x*x
      xy4=-4.*x*x*x*x*x*x+12.*x*x*x*x*x-14.*x*x*x*x+8.*x*x*x-2.*x*x
      xy5=0.5*xy1*xy1
      xy6=-24.*y*y*y*y*y+8.*y*y*y-4.*y
      xy7=4.*x*x*x-6.*x*x+2.*x
      xy8=4.*y*y*y-2.*y
      cav(i,j,k)=-8.*xnu*(24.*xy3+2.*xy7*(12.*y*y-2.)+xy2*(24.*x-12.))&
           -64.*(xy5*xy6-xy2*xy8*xy4)
   enddo
enddo
enddo
do i=1,nxyz
   uy(i,1,1)=-cav(i,1,1)+uy(i,1,1)
enddo

return
end subroutine cavite
end module navier_m
