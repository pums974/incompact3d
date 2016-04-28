module slfft2d_m
  use fft_m
  use, intrinsic :: iso_c_binding
implicit none
  include 'fftw3.f03'
  private
  type(C_PTR) :: fftw_plan_x1,fftw_plan_x2,fftw_plan_y1,fftw_plan_y2
  complex(C_DOUBLE_COMPLEX), pointer :: fftw_cmp(:,:)
  type(C_PTR) :: ptr
  
  public SLFFT2D, init_SLFFT2D,end_SLFFT2D
contains
!*****************************************************************************************
!
subroutine SLFFT2D(x,nx,ny,mx,my,nxm,nym,nwork,work,table,ntrigsX,ntrigsY,isign)
!
!*****************************************************************************************

USE param
  
  implicit none

integer, intent(in) :: nx,ny,mx,my,nxm,nym,nwork,ntrigsX,ntrigsY
integer, intent(in) :: isign
real(8), dimension(nwork) :: work,table
real(8),  intent(inout), dimension(mx,my) :: x
integer :: i, j


!************************FFT2D***(0-0)***********************************
    
   if (isign==-1) then
   
      ! On applique une TF reelle->complexe aux colonnes de ce tableau
      call fftw_execute_dft_r2c(fftw_plan_x1, x, fftw_cmp)
      
      !renormalisation
      fftw_cmp=fftw_cmp/nym/2
      
      ! On applique une TF complexe->complexe aux lignes de ce tableau
      call fftw_execute_dft(fftw_plan_y1, fftw_cmp, fftw_cmp)
      
      do j=1,ny
      do i=1,mx/2
         x(2*i-1,j)=real(fftw_cmp(i,j))/nym
         x(2*i,j)=aimag(fftw_cmp(i,j))/nym
      enddo
      enddo
   endif
   
   if (isign==1) then
      do j=1,ny
      do i=1,mx/2
         fftw_cmp(i,j)=cmplx(x(2*i-1,j),x(2*i,j))
      enddo
      enddo

      ! On applique une TF complexe->complexe aux lignes de ce tableau      
      call fftw_execute_dft(fftw_plan_y2, fftw_cmp, fftw_cmp)   

      ! On applique une TF complexe->reelle aux colonnes de ce tableau
      call fftw_execute_dft_c2r(fftw_plan_x2, fftw_cmp, x)
      
   endif
   
!
return
end subroutine SLFFT2D

subroutine init_SLFFT2D(nx,ny,mx,my)
implicit none
integer,intent(in) :: nx,ny,mx,my
real(8) :: x(mx,my)
integer :: err

!integer,parameter :: FFTW_FLAG=FFTW_ESTIMATE
!integer,parameter :: FFTW_FLAG=FFTW_MEASURE
integer,parameter :: FFTW_FLAG=FFTW_PATIENT

       ptr = fftw_alloc_complex(int(nx*ny, C_SIZE_T))
       call c_f_pointer(ptr, fftw_cmp, [nx,ny])

     err=fftw_init_threads()
     call fftw_plan_with_nthreads(1)


     fftw_plan_x1 = fftw_plan_many_dft_r2c(1, [nx], ny,   &
                                  x, [nx],                &
                                  1, mx,                  &
                                  fftw_cmp, [nx],         &
                                  1, nx,                  &
                                  FFTW_FLAG)

     fftw_plan_x2 = fftw_plan_many_dft_c2r(1, [nx], ny,   &
                                  fftw_cmp, [nx],         &
                                  1, nx,                  &
                                  x, [nx],                &
                                  1, mx,                  &
                                  FFTW_FLAG)

     fftw_plan_y1 = fftw_plan_many_dft(1, [ny], nx,   &
                                  fftw_cmp, [ny],     &
                                  nx, 1,              &
                                  fftw_cmp, [ny],     &
                                  nx, 1,              &
                                  FFTW_FORWARD, FFTW_FLAG)
                          
     fftw_plan_y2 = fftw_plan_many_dft(1, [ny], nx,   &
                                  fftw_cmp, [ny],     &
                                  nx, 1,              &
                                  fftw_cmp, [ny],     &
                                  nx, 1,              &
                                  FFTW_BACKWARD, FFTW_FLAG)
                                
end subroutine init_SLFFT2D

subroutine end_SLFFT2D
implicit none

      call fftw_destroy_plan(fftw_plan_x1)
      call fftw_destroy_plan(fftw_plan_x2)
      call fftw_destroy_plan(fftw_plan_y1)
      call fftw_destroy_plan(fftw_plan_y2)
      
           call fftw_free(ptr)
      
      call fftw_cleanup_threads()
      
end subroutine end_SLFFT2D

end module slfft2d_m
