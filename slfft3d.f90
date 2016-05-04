module slfft3d_m
  use fft_m
  use, intrinsic :: iso_c_binding
implicit none
  include 'fftw3.f03'
  private
  type(C_PTR) :: fftw_plan1,fftw_plan2
  complex(C_DOUBLE_COMPLEX), pointer :: fftw_cmp(:,:,:)
  type(C_PTR) :: ptr
  
  public SLFFT3D, init_SLFFT3D,end_SLFFT3D

  logical,parameter :: fftw=.true.

!  integer,parameter :: FFTW_FLAG=FFTW_ESTIMATE
  integer,parameter :: FFTW_FLAG=FFTW_MEASURE
  !integer,parameter :: FFTW_FLAG=FFTW_PATIENT
contains

!*****************************************************************************************
!
subroutine SLFFT3D(x,nxm,nym,nzm,mxx,mxy,myx,myy,work,table,isign,scale,mx,my,mz,nwork)
!
!*****************************************************************************************

USE param
  
  implicit none

integer, intent(in) :: nxm,nym,nzm,mxx,myx,mxy,myy,mx,my,mz,nwork
integer, intent(in) :: isign
real(8), intent(in) :: scale
real(8),  intent(inout), dimension(mx,my,mz) :: x
real(8), dimension(nwork) :: work
real(8), dimension(100+2*(nxm+nym+nzm)) :: table

integer :: i,j,k
!************************FFT2D***(0-0)***********************************
    
if(fftw) then ! FFTW

   if (isign==1) then

      call fftw_execute_dft_r2c(fftw_plan1, x, fftw_cmp)

      do k=1,nzm
      do j=1,nym
      do i=1,mx/2
         x(2*i-1,j,k)=real(fftw_cmp(i,j,k))*scale
         x(2*i,j,k)=-aimag(fftw_cmp(i,j,k))*scale
      enddo
      enddo
      enddo   
   
   
   endif

   if (isign==-1) then
   
      do k=1,nzm
      do j=1,nym
      do i=1,mx/2
         fftw_cmp(i,j,k)=cmplx(x(2*i-1,j,k),-x(2*i,j,k))
      enddo
      enddo
      enddo   
      
      call fftw_execute_dft_c2r(fftw_plan2, fftw_cmp, x)

      do k=1,nzm
      do j=1,nym
      do i=1,nxm
         x(i,j,k)=x(i,j,k)*scale
      enddo
      enddo
      enddo   

   endif

else

   if (isign==1) then
   
      call scfft3d(0,nxm,nym,nzm,scale,x,mxx,mxy,x,myx,myy,table,work,0)
      call scfft3d(isign,nxm,nym,nzm,scale,x,mxx,mxy,x,myx,myy,table,work,0)

   endif
   
   if (isign==-1) then

      call csfft3d(isign,nxm,nym,nzm,scale,x,mxx,mxy,x,myx,myy,table,work,0)     

   endif


endif
   
!
return
end subroutine SLFFT3D

subroutine init_SLFFT3D(nxm,nym,nzm,mx,my,mz)
implicit none
integer,intent(in) :: nxm,nym,nzm,mx,my,mz
real(8), dimension(mx,my,mz) :: x
integer :: err

if(fftw) then ! FFTW
print*,"Preparing FFTW ..."

       ptr = fftw_alloc_complex(int(nxm*nym*nzm, C_SIZE_T))
       call c_f_pointer(ptr, fftw_cmp, [mx/2,nym,nzm])

!     err=fftw_init_threads()
!     call fftw_plan_with_nthreads(1)


     fftw_plan1 = fftw_plan_many_dft_r2c(3, [nzm,nym,nxm], 1,   &
                                      x, [mz,my,mx],            &
                                      1, 1,                     &
                                      fftw_cmp, [nzm,nym,mx/2], &
                                      1, 1,                     &
                                      FFTW_FLAG);

     fftw_plan2 = fftw_plan_many_dft_c2r(3, [nzm,nym,nxm], 1,   &
                                      fftw_cmp, [nzm,nym,mx/2], &
                                      1, 1,                     &
                                      x, [mz,my,mx],            &
                                      1, 1,                     &
                                      FFTW_FLAG);
endif        
end subroutine init_SLFFT3D

subroutine end_SLFFT3D
implicit none

if(fftw) then ! FFTW
      call fftw_destroy_plan(fftw_plan1)
      call fftw_destroy_plan(fftw_plan2)
      
           call fftw_free(ptr)
      
!      call fftw_cleanup_threads()
      
endif
end subroutine end_SLFFT3D

end module slfft3d_m
