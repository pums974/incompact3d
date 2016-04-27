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
real(8), dimension(nwork) :: work,table
real(8),  intent(inout), dimension(mx,my) :: x
real(8), dimension(mx,my) :: x1,x2
integer :: i, j

  complex(C_DOUBLE_COMPLEX) :: fftw_cmp_x(nx,ny),fftw_cmp_y(ny,nx)
  type(C_PTR) :: fftw_plan_x1,fftw_plan_y,fftw_plan_x2


!************************FFT2D***(0-0)***********************************
    
      fftw_plan_x1 = fftw_plan_dft_r2c_1d(nx, &
                                    x(1:nx,1), fftw_cmp_x(1:nx,1), &
                                    FFTW_ESTIMATE)
                                    
      fftw_plan_x2 = fftw_plan_dft_c2r_1d(nx, &
                                    fftw_cmp_x(1:nx,1), x(1:nx,1), &
                                    FFTW_ESTIMATE)
                                    
     fftw_plan_y = fftw_plan_dft_1d(ny,   &
                                fftw_cmp_y(1:ny,1), fftw_cmp_y(1:ny,1),&
                                isign, FFTW_ESTIMATE);

   if (isign==-1) then
      ! On applique une TF reelle->complexe aux colonnes de ce tableau
      do j=1,ny
        call fftw_execute_dft_r2c(fftw_plan_x1, x(1:nx,j), fftw_cmp_x(1:nx,j))
      enddo
      
      !renormalisation
      fftw_cmp_x=fftw_cmp_x/nym/2
      fftw_cmp_y=transpose(fftw_cmp_x)
      
      ! On applique une TF complexe->complexe aux lignes de ce tableau
      do i=1,nx
        call fftw_execute_dft(fftw_plan_y, fftw_cmp_y(1:ny,i), fftw_cmp_y(1:ny,i))
      enddo        
      
      do j=1,ny
      do i=1,mx/2
         x(2*i-1,j)=real(fftw_cmp_y(j,i))/nym
         x(2*i,j)=aimag(fftw_cmp_y(j,i))/nym
      enddo
      enddo
   endif
   
   if (isign==1) then
      ! On applique une TF complexe->complexe aux lignes de ce tableau
      do j=1,ny
      do i=1,mx/2
         fftw_cmp_y(j,i)=cmplx(x(2*i-1,j),x(2*i,j))
      enddo
      enddo
      
      do i=1,nx
        call fftw_execute_dft(fftw_plan_y, fftw_cmp_y(1:ny,i), fftw_cmp_y(1:ny,i))
      enddo      

      fftw_cmp_x=transpose(fftw_cmp_y)

      ! On applique une TF reelle<-complexe aux colonnes de ce tableau
      do j=1,ny
        call fftw_execute_dft_c2r(fftw_plan_x2, fftw_cmp_x(1:nx,j), x(1:nx,j))
      enddo
      
   endif
   
      call fftw_destroy_plan(fftw_plan_x1)
      call fftw_destroy_plan(fftw_plan_x2)
      call fftw_destroy_plan(fftw_plan_y)
   
!
return
end subroutine SLFFT2D
end module slfft2d_m
