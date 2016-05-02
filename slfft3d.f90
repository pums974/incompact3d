module slfft3d_m
  use fft_m
  use, intrinsic :: iso_c_binding
implicit none
  include 'fftw3.f03'
  private
  type(C_PTR) :: fftw_plan_x1,fftw_plan_x2,fftw_plan_y1,fftw_plan_y2
  complex(C_DOUBLE_COMPLEX), pointer :: fftw_cmp(:,:)
  type(C_PTR) :: ptr
  
  public SLFFT3D, init_SLFFT3D,end_SLFFT3D

  logical,parameter :: fftw=.false.
contains

!*****************************************************************************************
!
subroutine SLFFT3D(x,y,nxm,nym,nzm,mxx,mxy,myx,myy,work,table,isign,scale,mx,my,mz,nwork)
!
!*****************************************************************************************

USE param
  
  implicit none

integer, intent(in) :: nxm,nym,nzm,mxx,myx,mxy,myy,mx,my,mz,nwork
integer, intent(in) :: isign
real(8), intent(in) :: scale
real(8),  intent(inout), dimension(mx,my,mz) :: x,y
real(8), dimension(nwork) :: work
real(8), dimension(100+2*(nxm+nym+nzm)) :: table

!************************FFT2D***(0-0)***********************************
    
if(fftw) then ! FFTW
!   if (isign==-1) then
!   
!      ! On applique une TF reelle->complexe aux colonnes de ce tableau
!      call fftw_execute_dft_r2c(fftw_plan_x1, x, fftw_cmp)
!      
!      !renormalisation
!      fftw_cmp=0.5*fftw_cmp/nym
!      
!      ! On applique une TF complexe->complexe aux lignes de ce tableau
!      call fftw_execute_dft(fftw_plan_y1, fftw_cmp, fftw_cmp)
!      
!      do j=1,ny
!      do i=1,mx/2
!         x(2*i-1,j)=real(fftw_cmp(i,j))/nym
!         x(2*i,j)=aimag(fftw_cmp(i,j))/nym
!      enddo
!      enddo
!   endif
!   
!   if (isign==1) then
!      do j=1,ny
!      do i=1,mx/2
!         fftw_cmp(i,j)=cmplx(x(2*i-1,j),x(2*i,j))
!      enddo
!      enddo

!      ! On applique une TF complexe->complexe aux lignes de ce tableau      
!      call fftw_execute_dft(fftw_plan_y2, fftw_cmp, fftw_cmp)   

!      ! On applique une TF complexe->reelle aux colonnes de ce tableau
!      call fftw_execute_dft_c2r(fftw_plan_x2, fftw_cmp, x)
!      
!   endif

else

   if (isign==-1) then
   
      call scfft3d(0,nxm,nym,nzm,scale,x,mxx,mxy,y,myx,myy,table,work,0)
      call scfft3d(1,nxm,nym,nzm,scale,x,mxx,mxy,y,myx,myy,table,work,0)

   endif
   
   if (isign==1) then

      call csfft3d(-1,nxm,nym,nzm,scale,x,mxx,mxy,y,myx,myy,table,work,0)     

   endif


endif
   
!
return
end subroutine SLFFT3D

subroutine init_SLFFT3D(nx,ny,mx,my)
implicit none
integer,intent(in) :: nx,ny,mx,my
real(8) :: x(mx,my)
integer :: err

!integer,parameter :: FFTW_FLAG=FFTW_ESTIMATE
integer,parameter :: FFTW_FLAG=FFTW_MEASURE
!integer,parameter :: FFTW_FLAG=FFTW_PATIENT

if(fftw) then ! FFTW
print*,"Preparing FFTW ..."

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
                                
endif        
end subroutine init_SLFFT3D

subroutine end_SLFFT3D
implicit none

if(fftw) then ! FFTW
      call fftw_destroy_plan(fftw_plan_x1)
      call fftw_destroy_plan(fftw_plan_x2)
      call fftw_destroy_plan(fftw_plan_y1)
      call fftw_destroy_plan(fftw_plan_y2)
      
           call fftw_free(ptr)
      
      call fftw_cleanup_threads()
      
endif
end subroutine end_SLFFT3D

end module slfft3d_m
