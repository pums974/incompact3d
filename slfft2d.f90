module slfft2d_m
  use fft_m
  use, intrinsic :: iso_c_binding
implicit none
  include 'fftw3.f03'
  private
  type(C_PTR) :: fftw_plan_x1,fftw_plan_x2,fftw_plan_y1,fftw_plan_y2
  complex(C_DOUBLE_COMPLEX), pointer :: fftw_cmp_x(:,:),fftw_cmp_y(:,:)
  real(C_DOUBLE), pointer :: fftw_real(:,:)
  type(C_PTR) :: p1,p2,p3
  real(C_DOUBLE), pointer :: ptr_real(:)
  complex(C_DOUBLE_COMPLEX), pointer :: ptr_cmpl1(:),ptr_cmpl2(:)
  
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
   
       fftw_real=x(1:nx,1:ny)
   
      ! On applique une TF reelle->complexe aux colonnes de ce tableau
      do j=1,ny
        ptr_real=>fftw_real(1:nx,j)
        ptr_cmpl1=>fftw_cmp_x(1:nx,j)
        call fftw_execute_dft_r2c(fftw_plan_x1, ptr_real, ptr_cmpl1)
      enddo
      
      !renormalisation
      fftw_cmp_x=fftw_cmp_x/nym/2
      fftw_cmp_y=transpose(fftw_cmp_x)
      
      ! On applique une TF complexe->complexe aux lignes de ce tableau
      do i=1,nx
        ptr_cmpl1=>fftw_cmp_y(1:ny,i)
        ptr_cmpl2=>fftw_cmp_y(1:ny,i)
        call fftw_execute_dft(fftw_plan_y1, ptr_cmpl1,ptr_cmpl2)
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
        ptr_cmpl1=>fftw_cmp_y(1:ny,i)
        ptr_cmpl2=>fftw_cmp_y(1:ny,i)
        call fftw_execute_dft(fftw_plan_y2, ptr_cmpl1,ptr_cmpl2)
      enddo      

      fftw_cmp_x=transpose(fftw_cmp_y)

      ! On applique une TF reelle<-complexe aux colonnes de ce tableau
      do j=1,ny
        ptr_real=>fftw_real(1:nx,j)
        ptr_cmpl1=>fftw_cmp_x(1:nx,j)
        call fftw_execute_dft_c2r(fftw_plan_x2, ptr_cmpl1, ptr_real)
      enddo
      
      x(1:nx,1:ny)=fftw_real
      
   endif
   
!
return
end subroutine SLFFT2D

subroutine init_SLFFT2D(nx,ny)
implicit none
integer,intent(in) :: nx,ny
integer :: err


       p1 = fftw_alloc_real(int(nx * ny, C_SIZE_T))
       p2 = fftw_alloc_complex(int(nx*ny, C_SIZE_T))
       p3 = fftw_alloc_complex(int(ny*nx, C_SIZE_T))
       call c_f_pointer(p1, fftw_real, [nx,ny])
       call c_f_pointer(p2, fftw_cmp_x, [nx,ny])
       call c_f_pointer(p3, fftw_cmp_y, [ny,nx])


        ptr_real=>fftw_real(1:nx,1)
        ptr_cmpl1=>fftw_cmp_x(1:nx,1)
        ptr_cmpl2=>fftw_cmp_y(1:ny,1)

!     err=fftw_init_threads()

!     call fftw_plan_with_nthreads(1)

      fftw_plan_x1 = fftw_plan_dft_r2c_1d(nx, &
                                    ptr_real, ptr_cmpl1, &
                                    FFTW_ESTIMATE)
                                    
      fftw_plan_x2 = fftw_plan_dft_c2r_1d(nx, &
                                    ptr_cmpl1, ptr_real, &
                                    FFTW_ESTIMATE)
                                    
     fftw_plan_y1 = fftw_plan_dft_1d(ny,   &
                                ptr_cmpl2, ptr_cmpl2,&
                                -1, FFTW_ESTIMATE);
     fftw_plan_y2 = fftw_plan_dft_1d(ny,   &
                                ptr_cmpl2, ptr_cmpl2,&
                                +1, FFTW_ESTIMATE);
end subroutine init_SLFFT2D

subroutine end_SLFFT2D
implicit none

      call fftw_destroy_plan(fftw_plan_x1)
      call fftw_destroy_plan(fftw_plan_x2)
      call fftw_destroy_plan(fftw_plan_y1)
      call fftw_destroy_plan(fftw_plan_y2)
      
           call fftw_free(p1)
           call fftw_free(p2)
           call fftw_free(p3)
      
!      call fftw_cleanup_threads()
      
end subroutine end_SLFFT2D

end module slfft2d_m
