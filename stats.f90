module stats_m
  use tools_m
  use forces_m
implicit none
contains
!***************************************************************************
!
subroutine stats(ux,uy,uz,gx,gy,gz,ppm,phi,phiss,temp,gtemp,&
                 epsi,tx,ty,tz,fx,fy,fz,di1,di2,px,py,pz,&
                 uxm1,uym1,uzm1,uxm2,uym2,uzm2)
!
!***************************************************************************

USE param
USE variables
USE aeroforce

implicit none

integer  :: i,j,k
real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,epsi
real(8),dimension(nx,ny,nz) :: tx,ty,tz,fx,fy,fz,di1,di2,px,py,pz
real(8),dimension(nx,ny,nz) :: uxm1,uym1,uzm1,uxm2,uym2,uzm2
real(8),dimension(nx,ny,nz) :: temp,gtemp
real(8),dimension(nx,ny,nz,nphi) :: phi,phiss
real(8),dimension(nxm,nym,nzm) :: ppm
integer :: nxyz

nxyz=nx*ny*nz
call minmax(ux,'Ux')
call minmax(uy,'Uy')
if (nz.gt.1) call minmax(uz,'Uz')

if (mod(itime,100).eq.0) then
!   open(20,file='tecplot.dat',form='formatted')
!   write(20,1001)
!   write(20,1002) nx,ny
!   do j=1,ny
!   do i=1,nx
!      write(20,901) (i-1)*dx,yp(j),ux(i,j,1),uy(i,j,1),epsi(i,j,1)
!   enddo
!   enddo
!   close(20)
! 1001 format('VARIABLES = "X","Y","U","V","e"')
! 1002 format('ZONE T="vitesse", I=',I7,' J=',I7,' F=POINT')
! 901  format(10F14.8)

   call save_restart(ux,uy,uz,gx,gy,gz,ppm,phi,phiss,temp,gtemp) 

endif

if (mod(itime,imodulo).eq.0 .and. itime.ge.idebmod) then 
   call snapshots(ux,uy,uz)
   isave=isave+1
endif

if ((ivirtuel.eq.1).and.(iaero.eq.1))then
   call aerof(ux,uy,uz,ppm,epsi,px,py,pz,&
              tx,ty,tz,fx,fy,fz,di1,di2,&
              uxm1,uym1,uzm1,uxm2,uym2,uzm2)
   call aerofd(ux,uy,uz,ppm,epsi,px,py,pz,&
               tx,ty,tz,di1,di2)
endif
return
end subroutine stats

!********************************************************************
!
subroutine snapshots(ux,uy,uz)
!
!********************************************************************

USE param
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz
integer :: longueur,num 
real(8) :: wzmin, wzmax 
character(len=4) suffix
character(len=20) nfichier

num=isave
call numcar (num,suffix)
longueur=index(nchamp,' ')-1
nfichier=nchamp(1:longueur)//suffix
longueur=index(nfichier,' ')-1


open(12,file=nfichier(1:longueur),form='unformatted',status='unknown')
if (nz.gt.1) then
   write(12) ux,uy,uz
else
   write(12) ux,uy
endif
close(12)
!
return
end subroutine snapshots


!********************************************************************
!
subroutine save_restart(ux,uy,uz,gx,gy,gz,ppm,phi,phiss,temp,gtemp)
!
!*******************************************************************

USE param
USE variables

implicit none

integer :: num,longueur
real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,temp,gtemp
real(8),dimension(nx,ny,nz,nphi) :: phi,phiss
real(8),dimension(nxm,nym,nzm) :: ppm
character(len=4) suffix
character(len=20) nfichier


open(11,file='restart',form='unformatted',status='unknown')
if (iscalaire==0) then
  if (iconvect == 0) then
   if (nz.gt.1) then
      write(11) ux,uy,uz,ppm,gx,gy,gz,dpdyx1,dpdyxn,dpdzx1,dpdzxn,dpdxy1,&
                dpdxyn,dpdzy1,dpdzyn,dpdxz1,dpdxzn,dpdyz1,dpdyzn
   else
      write(11) ux,uy,ppm,gx,gy,dpdyx1,dpdyxn,dpdxy1,dpdxyn
   endif
  else
   if (nz.gt.1) then
      write(11) ux,uy,uz,ppm,gx,gy,gz,temp,gtemp,dpdyx1,dpdyxn,dpdzx1,&
                dpdzxn,dpdxy1,dpdxyn,dpdzy1,dpdzyn,dpdxz1,dpdxzn,dpdyz1,dpdyzn
   else
      write(11) ux,uy,ppm,gx,gy,temp,gtemp,dpdyx1,dpdyxn,dpdxy1,dpdxyn
   endif
  endif
else
  if (iconvect == 0) then
   if (nz.gt.1) then
      write(11) ux,uy,uz,phi,ppm,gx,gy,gz,dpdyx1,dpdyxn,dpdzx1,&
                 dpdzxn,dpdxy1,dpdxyn,dpdzy1,dpdzyn,dpdxz1,dpdxzn,dpdyz1,dpdyzn
   else
      write(11) ux,uy,phi,ppm,gx,gy,dpdyx1,dpdyxn,dpdxy1,dpdxyn
   endif
  else
   if (nz.gt.1) then
      write(11) ux,uy,uz,phi,ppm,gx,gy,gz,temp,gtemp,dpdyx1,&
                dpdyxn,dpdzx1,dpdzxn,dpdxy1,dpdxyn,dpdzy1,dpdzyn,dpdxz1,dpdxzn,dpdyz1,dpdyzn
   else
      write(11) ux,uy,phi,ppm,gx,gy,temp,gtemp,dpdyx1,dpdyxn,dpdxy1,dpdxyn
   endif
  endif
endif
close(11)


return
end subroutine save_restart
end module stats_m