!*****************************************************************************************
!
subroutine SLFFT2D(x,nx,ny,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,isign)
!
!*****************************************************************************************

USE param
  
implicit none

integer, intent(in) :: nx,ny,mx,my,nxm,nym,nwork,ntrigsX,ntrigsY
integer, intent(in) :: isign
real(8),  intent(inout), dimension(mx,my) :: x
real(8), dimension(nwork) :: work,table
real(8), dimension(mx,my) :: x1,x2
real(8), dimension(0:ntrigsX-1) :: trigsX
real(8), dimension(0:ntrigsY-1) :: trigsY
integer, parameter :: nfax = 19
integer, dimension(0:nfax-1) :: ifaxX, ifaxY
integer :: inc, jump
integer :: i, j

!************************FFT2D***(0-0)***********************************
    
   ! Initialisation des factorisations et tableaux trigonométrique
   ! Note : L'appel a jmcftfax calcule les cos et sin(2*pi/n), tandis
   ! que les appels suivants calculent les cos(pi/n)
   ! Il y a un peu de gaspillage, mais c'est plus simple comme ca
   
   call fftfax(nxm,ifaxX,trigsX)
   call cftfax(nym,ifaxY,trigsY)

   if (isign==-1) then
      ! On applique une TF reelle->complexe aux colonnes de ce tableau
      inc = 1
      jump = mx
      call rfftmlt(x,work,trigsX,ifaxX,inc,jump,nx,ny,isign)
      ! On applique une TF complexe->complexe aux lignes de ce tableau
      inc = mx
      jump = 1
      do j=1,my
      do i=1,mx/2
         x1(i,j)=x(2*i-1,j)
         x2(i,j)=x(2*i,j)
      enddo
      enddo
      call cfftmlt(x1,x2,work,trigsY,ifaxY,inc,jump,ny,nx,isign)
      do j=1,my
      do i=1,mx/2
         x(2*i-1,j)=x1(i,j)
         x(2*i,j)=x2(i,j)
      enddo
      enddo
      x(:,:)=x(:,:)/nym
   endif

   if (isign==1) then
      ! On applique une TF complexe->complexe aux lignes de ce tableau
      inc = mx
      jump = 1
      do j=1,my
      do i=1,mx/2
         x1(i,j)=x(2*i-1,j)
         x2(i,j)=x(2*i,j)
      enddo
      enddo
      call cfftmlt(x1,x2,work,trigsY,ifaxY,inc,jump,ny,nx,isign)
      do j=1,my
      do i=1,mx/2+2
         x(2*i-1,j)=x1(i,j)
         x(2*i,j)=x2(i,j)
      enddo
      enddo
      ! On applique une TF reelle->complexe aux colonnes de ce tableau
      inc = 1
      jump = mx
      call rfftmlt(x,work,trigsX,ifaxX,inc,jump,nx,ny,isign)
      x(:,:)=x(:,:)
   endif
!
!*************VERSION QUI MARCHE***********************************    
!        call scfft2d(0,nx,ny,1.,x,mx,x,nxm/2+1,table,work,0)
!     if ( isign== -1) then
!        call scfft2d(isign,nx,ny,1.,x,mx,x,nxm/2+1,table,work,0)
!        x(:,:)=x(:,:)/nxm/nym
!     endif
!     if ( isign==  1) then 
!        call csfft2d(isign,nx,ny,1.,x,nxm/2+1,x,mx,table,work,0)
!        x(:,:)=x(:,:)
!     endif
! 
!******************************************************************
!
return
end subroutine SLFFT2D
