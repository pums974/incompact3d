!********************************************************************
!
PROGRAM incompact3d
!
!********************************************************************
USE param
USE parfiX;USE parfiY;USE parfiZ
USE derivX;USE derivY;USE derivZ
USE aeroforce
USE variables
USE convection
use parametre_m
use filtre_m
use schemas_m
use waves_m
use navier_m
use geom_complex_m
use convdiff_m
use scalar_m
use convection_m
use body_m
use tools_m
use poisson_m
use stats_m
use derivevitesse_m

implicit none

real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,fx,fy,fz
real(8),dimension(nx,ny,nz) :: temp,gtemp,htemp
real(8),dimension(nx,ny,nz,nphi) :: phi,phis,phiss
real(8),dimension(nx,ny,nz) :: sy1,sy2,sy3,sy4,sy5,sy6
real(8),dimension(mx,my,mz) :: sy7,sy8,sy9,sy10,sy11,sy12
real(8),dimension(nx,ny,nz) :: epsi
real(8),dimension(mx,my,mz) :: di2,di1 ! TODO
real(8),dimension(nxm,nym,nzm) :: ppm,ppm1!,err
real(8),dimension(nx,ny,nz) :: uxm1,uym1,uzm1,uxm2,uym2,uzm2
!
real(8),dimension(mx,mz,ny/2,5) :: a,a2
real(8),dimension(mx,mz,ny,5) :: a3 
real(8),dimension(mx,mz,ny) :: d,d1,e,c
real(8),dimension(mx,mz) :: sr
real(8),dimension(mx,mz) :: a1,b1
!
real(8),dimension(nwork) :: work,table
!real(8),dimension(100+2*(nxm+nym+nzm)) :: table
!
real(8),dimension(nxm,nym,nzm) :: epsidec
real(8),dimension(nx,ny,nz) :: px,py,pz
integer                     :: i,j,numvis
integer, parameter          :: debug=0
character(len=4)            :: suffixe
!
real(8),dimension(nx,ny,nz) :: Udamien,Vdamien
real(8)                     :: aaaa
!
   call parametre()
   call filtres()
   call schemas()
   call waves ()  
   call initial (ux,uy,uz,gx,gy,gz,fx,fy,fz,phi,ppm,temp,gtemp)
   
   call init_SLFFT2D(nxm,nym,mx,my)
   if (ivirtuel.eq.2) then !new method
       call gene_epsi(epsi,epsidec)
       call verif_epsi(epsi)
   endif
!
   numvis=1
   do itime=idebut,ifin
      t=(itime-1)*dt
      write(*,1001) itime,t
 1001 format('Time step =',i7,', Time unit =',F9.3)
      call time_measure (0)
      if (itime==idebut+1) then
        ntics=0
        average_t=0.
      endif
      do itr=1,iavance_temps
!        
         call inflow    (ux,uy,uz)
         call outflow   (ux,uy,uz)
!
         call convdiff  (ux,uy,uz,sy1,sy2,sy3,sy4,sy5,sy6,sy7,sy8,sy9,sy10,sy11,sy12,di1,di2)
         call scalar    (ux,uy,uz,phi,phis,phiss,sy1,sy2,sy3,sy4,sy5,sy6,di1)
         if (iconvect==1)  then
             call force_gravite (sy8,sy9,temp)
             call energie (ux,uy,uz,temp,gtemp)
         endif
             call intt (ux,uy,uz,gx,gy,gz,sy7,sy8,sy9)
         
         if (ivirtuel==1) then !old school
             call solid_body(ux,uy,uz,epsi,ppm,sy1,sy2,sy3,sy4,sy5,sy6,sy7,sy8,sy9,di1,di2)
         endif
         call pre_correc(ux,uy,uz)
         call divergence(ppm,ux,uy,uz,sy4,sy5,sy6,di1,di2,&
                         sy7,sy8,sy9,sy10,sy11,sy12,&
                         epsi,sy1,sy2,sy3,work,table,1)

         call poisson   (ppm,sy7,sy8,sy9,sy10,sy11,sy12,&
                         work,table,a,a3,e,c,a2,sr,a1,b1,d1,d)
         if (iconvect==1)  then
             ppm=ppm/rho
         endif
         call gradpression(ppm,sy1,sy2,sy3,sy4,sy5,sy6,sy7,sy8,sy9,di1,di2)
         call corgp       (ux,uy,uz,sy1,sy2,sy3)
          
         call divergence  (ppm1,ux,uy,uz,sy4,sy5,sy6,di1,di2,&
                           sy7,sy8,sy9,sy10,sy11,sy12,&
                           epsi,sy1,sy2,sy3,work,table,2)
      enddo
      call time_measure (1)
      call stats       (ux,uy,uz,gx,gy,gz,ppm,phi,phiss,temp,gtemp,epsi,&
                        sy1,sy2,sy3,sy4,sy5,sy6,di1,di2,sy7,sy8,sy9,&
                        uxm1,uym1,uzm1,uxm2,uym2,uzm2)
!     AFFICHAGE :
!
!        CHAMP DE VORTICITÉ
           if(mod(itime,100).eq.0) then
!            call forcage_flugrange_x(uy)!IBM
            call derx (py,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
!            call forcage_flugrange_y(ux)!IBM
            call dery (px,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1) 
            pz(:,:,:)=py(:,:,:)-px(:,:,:)
            open(10,file='tampon.dat',form='formatted')
            open(20,file='temp.dat',form='formatted')
            do j=1,ny
               do i=1,nx
                  write(10,*)(i-1)*dx,yp(j),ux(i,j,max(1,nz/2))
                  write(20,*) dx*(i-1),yp(j),temp(i,j,max(1,nz/2)) !,uy(i,j,1)
               enddo
               write(10,*)
               write(20,*)
            enddo 
            close(10)
            close(20)
!!
            call system('gnuplot visu2.gnu')
            write(suffixe,'(i4)') numvis+1000
            call system('mv tampon.jpeg wz'//suffixe(1:4)//'.jpg')
            call system('mv temp.jpeg temperature'//suffixe(1:4)//'.jpg')
            numvis=numvis+1
!!
         endif
!stop
!        CHAMP DE PRESSION
!         if(mod(itime,10000).eq.0.and.itime.ge.1) then
!            open(10,file='tampon.dat',form='formatted')
!            do j=1,nym
!               do i=1,nxm
!                  write(10,*) (i-1)*dx+dx/2.,ppm(i,(8+3.*dy)/dy,1)/dt
!                  write(10,*) yp(j)+dy/2.,ppm(6/dx,j,1)/dt
!                  write(10,*) (i-1)*dx+dx/2.,yp(j)+dy/2.,ppm(i,j,1)/dt
!               enddo
!               write(10,*)
!            enddo 
!            close(10)
!!
!            call system('gnuplot visu2.gnu')
!!            write(suffixe,'(i4)') numvis+1000
!!            call system('mv tampon.jpeg image'//suffixe(1:4)//'.jpeg')
!!            numvis=numvis+1
!        endif
!stop
!        CONDITION AUX LIMITES SUR LE GRADIENT DE PRESSION :
!  123    if(mod(itime,1).eq.0.and.itime.ge.1) then
!            open(10,file='tampon.dat',form='formatted')
!            call forcage_flugrange_x(ux,xi,xf,nobjx,ny,nz)!IBM
!            call derxx (du2dx,ux,di1,sx,sfx ,ssx ,swx ,nx,ny,nz,0)
!            call forcage_flugrange_y(ux,yi,yf,nobjy,nx,nz)!IBM
!            call deryy (du2dy,ux,di1,di2,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
!            call forcage_flugrange_x(uy,xi,xf,nobjx,ny,nz)!IBM
!            call derxx (dv2dx,uy,di1,sx,sfxp,ssxp,swxp,nx,ny,nz,1)
!            call forcage_flugrange_y(uy,yi,yf,nobjy,nx,nz)!IBM
!            call deryy (dv2dy,uy,di1,di2,sy,sfy,ssy,swy,nx,ny,nz,0)
!!            gpx(:,:,:)=sy1(:,:,:)/dt-xnu*(du2dx(:,:,:)+du2dy(:,:,:))
!!            gpy(:,:,:)=sy2(:,:,:)/dt-xnu*(dv2dx(:,:,:)+dv2dy(:,:,:))
!            gpx(:,:,:)=sy1(:,:,:)/dt
!            gpy(:,:,:)=sy2(:,:,:)/dt
!            do j=1,ny
!               do i=1,nx
!                  write(10,*)dx*(i-1),yp(j),gpx(i,j,1)
!               enddo
!              write(10,*)
!            enddo 
!            close(10)
!            call system('gnuplot visu2.gnu')
!         endif
!stop
!        CHAMPS DE VITESSE :
!        if(mod(itime,1000).eq.0.and.itime.ge.1) then
!         open(10,file='ux.dat',form='formatted')  
!         open(20,file='temp.dat',form='formatted')
!         do j=1,ny
!            do i=1,nx
!               write(10,*) dx*(i-1),yp(j),ux(i,j,1)
!               write(20,*) dx*(i-1),yp(j),temp(i,j,1)
!            enddo
!            write(10,*)
!            write(20,*)
!         enddo
!         close(10)
!         close(20)
!        endif
!
!         open(20,file='uy_hd1.dat',form='formatted')
!         do j=1,ny
!            y=yp(j)
!            do i=1,nx
!               x=dx*(i-1)
!               write(20,*)x,y,uy(i,j,1)
!            enddo
!         enddo
!         close(20)
!stop
!        CONDITIONS DE DIVERGENCE, CARTE D'ERREURS :
!         if(mod(itime,1).eq.0.and.itime.ge.1) then
!         open(10,file='tampon.dat',form='formatted')
!            do j=1,nym
!               do i=1,nxm
!                  write(10,*) (i-1)*dx+dx/2.,yp(j)+dy/2.,err(i,j,1)/dt
!               enddo
!               write(10,*)
!            enddo 
!            close(10)
!            call system('gnuplot visu.gnu')
!         endif
!stop
!
!        ÉVOLUTION - STATIONNARITÉ :
!         if(itime.eq.10)then
!            ux_old(:,:,:)=ux(:,:,:)
!            uy_old(:,:,:)=uy(:,:,:)
!            er_ux=0.
!            er_uy=0.
!         endif
!         if(mod(itime,10).eq.0.and.itime.gt.10) then
!            er_ux=0.
!            er_uy=0.
!            do j=1,ny
!            do i=1,nx
!               if(esp(i,j,1).eq.0.)then
!                  er_ux=er_ux+(ux(i,j,1)-ux_old(i,j,1))**2
!                  er_uy=er_uy+(uy(i,j,1)-uy_old(i,j,1))**2
!               endif
!            enddo
!            enddo
!            er_ux=sqrt(er_ux/(nx*ny))
!            er_uy=sqrt(er_uy/(nx*ny))
!            open(43,file='evol_ux_r8.dat',form='formatted')
!            open(86,file='evol_uy_r8.dat',form='formatted')
!            write(43,*)t,er_ux
!            write(86,*)t,er_ux
!            ux_old(:,:,:)=ux(:,:,:)
!            uy_old(:,:,:)=uy(:,:,:)
!         endif
!stop
!        VORTICITÉ DAMIEN :
! 23      open(100,file='U.dat',form='formatted')
!         open(200,file='V.dat',form='formatted')
!         if(mod(itime,1).eq.0.and.itime.ge.1) then
!            do j=1,ny
!               do i=1,nx
!                  read(100,*)aaaa,aaaa,Udamien(i,j,1)
!                  read(200,*)aaaa,aaaa,Vdamien(i,j,1)
!               enddo
!               read(100,*)
!               read(200,*)
!            enddo
!!            call forcage_flugrange_x(Vdamien,xi,xf,nobjx,ny,nz)!IBM
!            call derx (py,Vdamien,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
!!            call forcage_flugrange_y(Udamien,yi,yf,nobjy,nx,nz)!IBM
!            call dery (px,Udamien,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1) 
!            pz(:,:,:)=py(:,:,:)-px(:,:,:)
!!
!            open(10,file='tampon.dat',form='formatted')
!            do j=1,ny
!               do i=1,nx
!                  write(10,*) (i-1)*dx,yp(j),pz(i,j,1)
!               enddo
!               write(10,*)
!            enddo 
!            close(10)
!            call system('gnuplot visu2.gnu')
!         endif
!stop
   enddo
   call end_SLFFT2D
!
end PROGRAM incompact3d

