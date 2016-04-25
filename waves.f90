!***********************************************************
!
subroutine waves ()
!
!***********************************************************
!
USE derivX 
USE derivY 
USE derivZ 
USE param
USE variables
!
implicit none
!
integer :: i,j,k
real(8) :: w,wp 
!
xkx(:)=0. ; xk2(:)=0. ; yky(:)=0. ; yk2(:)=0.
zkz(:)=0. ; zk2(:)=0.
!
if (nclx==0) then
   do i=1,mx,2
      w=twopi*0.5*(i-1)/(mx-2)
      wp=acix6*2.*dx*sin(w/2.)+(bcix6*2.*dx)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaix6*cos(w))
      xkx(i)=(mx-2)*wp/xlx
      xkx(i+1)=xkx(i)
      exs(i)=(mx-2)*w/xlx
      exs(i+1)=(mx-2)*w/xlx
      xk2(i)=xkx(i)*xkx(i)
      xk2(i+1)=xk2(i) 
   enddo
endif
if ((nclx==1).or.(nclx==2)) then
   do i=1,mx
      w=twopi*0.5*(i-1)/(mx-2)
      wp=acix6*2.*dx*sin(w/2.)+(bcix6*2.*dx)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaix6*cos(w))
      xkx(i)=(mx-2)*wp/xlx
      exs(i)=(mx-2)*w/xlx
      xk2(i)=xkx(i)*xkx(i)
   enddo
   xkx(mx-1)=0.
   xk2(mx-1)=0.
   xkx(mx  )=0.
   xk2(mx  )=0.
endif
if (ncly==0) then
   do j=1,my,2
      w=twopi*0.5*(j-1)/(my-2)
      wp=aciy6*2.*dy*sin(w/2.)+(bciy6*2.*dy)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaiy6*cos(w))
      if (istret.ne.0) yky(j)=(my-2)*wp
      if (istret.eq.0) yky(j)=(my-2)*wp/yly
      if (istret.ne.0) yky(j+1)=(my-2)*wp
      if (istret.eq.0) yky(j+1)=(my-2)*wp/yly
      eys(j)=(my-2)*w/yly
      eys(j+1)=(my-2)*w/yly
      yk2(j)=yky(j)*yky(j)
      yk2(j+1)=yk2(j) 
   enddo
endif
if ((ncly==1).or.(ncly==2)) then
   do j=1,my
      w=twopi*0.5*(j-1)/(my-2)
      wp=aciy6*2.*dy*sin(w/2.)+(bciy6*2.*dy)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaiy6*cos(w))
      if (istret.ne.0) yky(j)=(my-2)*wp
      if (istret.eq.0) yky(j)=(my-2)*wp/yly
      eys(j)=(my-2)*w/yly
      yk2(j)=yky(j)*yky(j)
   enddo
   yky(my-1)=0.
   yk2(my-1)=0.
   yky(my  )=0.
   yk2(my  )=0.  
endif
if (mz.gt.1) then
   if (nclz==0) then
      do k=1,mz,2
         w=twopi*0.5*(k-1)/(mz-2)
         wp=aciz6*2.*dz*sin(w/2.)+(bciz6*2.*dz)*sin(3./2.*w)
         wp=wp/(1.+2.*alcaiz6*cos(w))
         zkz(k)=(mz-2)*wp/zlz
         zkz(k+1)=zkz(k)
         ezs(k)=(mz-2)*w/zlz
         ezs(k+1)=(mz-2)*w/zlz
         zk2(k)=zkz(k)*zkz(k)
         zk2(k+1)=zk2(k) 
      enddo
   endif
   if ((nclz==1).or.(nclz==2)) then
      do k=1,mz
         w=twopi*0.5*(k-1)/(mz-2)
         wp=aciz6*2.*dz*sin(w/2.)+(bciz6*2.*dz)*sin(3./2.*w)
         wp=wp/(1.+2.*alcaiz6*cos(w))
         zkz(k)=(mz-2)*wp/zlz
         ezs(k)=(mz-2)*w/zlz
         zk2(k)=zkz(k)*zkz(k)
      enddo
      zkz(mz-1)=0.
      zk2(mz-1)=0.
      zkz(mz  )=0.
      zk2(mz  )=0. 
   endif
endif
!
return
end subroutine waves
