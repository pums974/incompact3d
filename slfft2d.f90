module slfft2d_m
  use fft_m
implicit none
contains
!*****************************************************************************************
!
subroutine SLFFT2D(x,nx,ny,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,isign)
!
!*****************************************************************************************

USE param
  
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

integer, intent(in) :: nx,ny,mx,my,nxm,nym,nwork,ntrigsX,ntrigsY
integer, intent(in) :: isign
real(8),  intent(inout), dimension(mx,my) :: x
real(8),  dimension(mx,my) :: x5
real(8), dimension(nwork) :: work,table
real(8), dimension(mx,my) :: x1,x2
real(8), dimension(nx,ny) :: x3,x4
integer :: inc, jump
integer :: i, j

  complex(C_DOUBLE_COMPLEX) :: fftw_cmp_x(nx),fftw_cmp_y(ny)
  real(C_double) :: fftw_real_x(nx),fftw_real_y(ny)
  type(C_PTR) :: fftw_plan_x1,fftw_plan_y,fftw_plan_x2


!************************FFT2D***(0-0)***********************************
    
      fftw_plan_x1 = fftw_plan_dft_r2c_1d(nx, &
                                    fftw_real_x, fftw_cmp_x, &
                                    FFTW_ESTIMATE)
                                    
      fftw_plan_x2 = fftw_plan_dft_c2r_1d(nx, &
                                    fftw_cmp_x, fftw_real_x, &
                                    FFTW_ESTIMATE)
                                    
     fftw_plan_y = fftw_plan_dft_1d(ny,   &
                                fftw_cmp_y, fftw_cmp_y,&
                                isign, FFTW_ESTIMATE);

   if (isign==-1) then
      ! On applique une TF reelle->complexe aux colonnes de ce tableau
      do j=1,ny
        fftw_real_x=x(1:nx,j)
        call fftw_execute_dft_r2c(fftw_plan_x1, fftw_real_x, fftw_cmp_x)
        x3(:,j)=real(fftw_cmp_x)/nym/2
        x4(:,j)=aimag(fftw_cmp_x)/nym/2
      enddo
      
      ! On applique une TF complexe->complexe aux lignes de ce tableau
      do i=1,nx
        fftw_cmp_y=cmplx(x3(i,1:ny),x4(i,1:ny))
        call fftw_execute_dft(fftw_plan_y, fftw_cmp_y, fftw_cmp_y)
        x3(i,:)=real(fftw_cmp_y)
        x4(i,:)=aimag(fftw_cmp_y)
      enddo        
      
      do j=1,ny
      do i=1,mx/2
         x(2*i-1,j)=x3(i,j)
         x(2*i,j)=x4(i,j)
      enddo
      enddo
      x(:,:)=x(:,:)/nym
   endif
   
   if (isign==1) then
      ! On applique une TF complexe->complexe aux lignes de ce tableau
      do j=1,my
      do i=1,mx/2
         x1(i,j)=x(2*i-1,j)
         x2(i,j)=x(2*i,j)
      enddo
      enddo
      
      do i=1,nx
        fftw_cmp_y=cmplx(x1(i,1:ny),x2(i,1:ny))
        call fftw_execute_dft(fftw_plan_y, fftw_cmp_y, fftw_cmp_y)
        x3(i,:)=real(fftw_cmp_y)
        x4(i,:)=aimag(fftw_cmp_y)
      enddo      

      ! On applique une TF reelle<-complexe aux colonnes de ce tableau
      do j=1,ny
        fftw_cmp_x=cmplx(x3(1:nx,j),x4(1:nx,j))
        call fftw_execute_dft_c2r(fftw_plan_x2, fftw_cmp_x, fftw_real_x)
        x(1:nx,j)=fftw_real_x
      enddo
      
   endif
   
      call fftw_destroy_plan(fftw_plan_x1)
      call fftw_destroy_plan(fftw_plan_x2)
      call fftw_destroy_plan(fftw_plan_y)
   
!
return
end subroutine SLFFT2D
end module slfft2d_m
