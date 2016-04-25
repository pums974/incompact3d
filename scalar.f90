module scalar_m
implicit none
contains
!************************************************************
!
subroutine scalar(ux,uy,uz,phi,phis,phiss,sy1,sy2,sy3,sy4,sy5,sy6,di1)
!
!************************************************************

USE param
USE variables

implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,sy1,sy2,sy3,sy4,sy5,sy6,di1
real(8),dimension(nx,ny,nz,nphi) :: phi,phis,phiss
integer :: ns,i,j,k,nxyz

if (iscalaire.eq.1) then
do ns=1,nphi
   do k=1,nz
   do j=1,ny
   do i=1,nx
      sy1(i,j,k)=ux(i,j,k)*phi(i,j,k,ns)
      sy2(i,j,k)=uy(i,j,k)*phi(i,j,k,ns)
      sy3(i,j,k)=uz(i,j,k)*phi(i,j,k,ns)
   enddo
   enddo
   enddo
   call derx (sy4,sy1,di1,    sx,ffxp,fsxp,fwxp,    nx,ny,nz,1)
   call dery (sy5,sy2,sy1,di1,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
   call derz (sy6,sy3,sy2,    sz,ffzp,fszp,fwzp,    nx,ny,nz,1)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      phis(i,j,k,ns)=-(sy4(i,j,k)+sy5(i,j,k)+sy6(i,j,k))
      sy1(i,j,k)=phi(i,j,k,ns)
   enddo
   enddo
   enddo
   call derxx (sy2,sy1,di1,sx,sfxp,ssxp,swxp,nx,ny,nz,1)
   call deryy (sy3,sy1,di1,sy5,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
   if (istret.ne.0) then 
      call deryy (sy3,sy1,di1,sy5,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
      call dery (sy4,sy1,di1,sy5,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         sy3(i,j,k)=sy3(i,j,k)*pp2y(j)-pp4y(j)*sy4(i,j,k)
      enddo
      enddo
      enddo
   else
      call deryy (sy3,sy1,di1,sy5,sy,sfyp,ssyp,swyp,nx,ny,nz,1) 
   endif
   call derzz (sy4,sy1,di1,sz,sfzp,sszp,swzp,nx,ny,nz,1)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      phis(i,j,k,ns)=xnu/sc*(sy2(i,j,k)+sy3(i,j,k)+sy4(i,j,k))+phis(i,j,k,ns)
   enddo
   enddo
   enddo
   if (itr.eq.1) then
      do k=1,nz
      do j=1,ny
      do i=1,nx
         phi(i,j,k,ns)=gdt(itr)*phis(i,j,k,ns)+phi(i,j,k,ns)
         phiss(i,j,k,ns)=phis(i,j,k,ns)
      enddo
      enddo
      enddo
   else
      do k=1,nz
      do j=1,ny
      do i=1,nx
         phi(i,j,k,ns)=adt(itr)*phis(i,j,k,ns)+bdt(itr)*phiss(i,j,k,ns)+phi(i,j,k,ns)
         phiss(i,j,k,ns)=phis(i,j,k,ns)
      enddo
      enddo
      enddo
   endif
enddo

endif

return
end subroutine scalar

!************************************************************
!
subroutine initscalar(phi)
!
!************************************************************

USE variables

implicit none

real(8),dimension(nx,ny,nz,nphi) :: phi

phi=1.

return
end subroutine initscalar

end module scalar_m