module convdiff_m
  use geom_complex_m
  use derivevitesse_m
  use filtre_m
implicit none
contains
!********************************************************************
!
subroutine convdiff(ux,uy,uz,sy1,sy2,sy3,sy4,sy5,sy6,sy7,sy8,sy9,sy10,sy11,sy12,di1,di2)
! 
!********************************************************************
USE param
USE variables
USE convection
!
implicit none
!
integer :: nxyz1,ijk,i,j,k,ns
real(8),dimension(nx,ny,nz) :: sy1,sy2,sy3,sy4,sy5,sy6,sy7,sy8,sy9,sy10,sy11,sy12,ux,uy,uz,di1,di2
!
nxyz1=nx*ny*nz
 if (ivirtuel==2) then 
  if (iskew==0) then !UROTU!
      if (nz.gt.1) then!3D
         call forcage_flugrange_y(uz)!IBM
         call dery (sy1,uz,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call forcage_flugrange_z(uy)!IBM
         call derz (sy2,uy,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call forcage_flugrange_z(ux)!IBM
         call derz (sy3,ux,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call forcage_flugrange_x(uz)!IBM
         call derx (sy7,uz,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_x(uy)!IBM
         call derx (sy8,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1) 
         call forcage_flugrange_y(ux)!IBM         
         call dery (sy9,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy1(ijk,1,1)=sy1(ijk,1,1)-sy2(ijk,1,1)
            sy3(ijk,1,1)=sy3(ijk,1,1)-sy7(ijk,1,1)
            sy8(ijk,1,1)=sy8(ijk,1,1)-sy9(ijk,1,1)
         enddo
         do ijk=1,nxyz1
            sy2(ijk,1,1)=sy3(ijk,1,1)*uz(ijk,1,1)-sy8(ijk,1,1)*uy(ijk,1,1)
            sy7(ijk,1,1)=sy8(ijk,1,1)*ux(ijk,1,1)-sy1(ijk,1,1)*uz(ijk,1,1)
            sy9(ijk,1,1)=sy1(ijk,1,1)*uy(ijk,1,1)-sy3(ijk,1,1)*ux(ijk,1,1)
         enddo
      else!2D
         call forcage_flugrange_x(uy)!IBM
         call derx (sy8,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(ux)!IBM
         call dery (sy9,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1) 
         do ijk=1,nxyz1
            sy8(ijk,1,1)= sy8(ijk,1,1)-sy9(ijk,1,1)
            sy2(ijk,1,1)=-sy8(ijk,1,1)*uy(ijk,1,1)
            sy7(ijk,1,1)= sy8(ijk,1,1)*ux(ijk,1,1)
         enddo
      endif
   else!skew
      if (nz.gt.1) then!3D
         do ijk=1,nxyz1
            sy1(ijk,1,1)=ux(ijk,1,1)*ux(ijk,1,1)
            sy4(ijk,1,1)=uy(ijk,1,1)*uy(ijk,1,1)
            sy3(ijk,1,1)=uz(ijk,1,1)*uz(ijk,1,1)
            sy5(ijk,1,1)=uy(ijk,1,1)*uz(ijk,1,1)
            sy8(ijk,1,1)=ux(ijk,1,1)*uz(ijk,1,1)
            sy6(ijk,1,1)=ux(ijk,1,1)*uy(ijk,1,1)
         enddo
         call forcage_flugrange_x(sy1)!IBM
         call derx (sy2,sy1,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(sy4)!IBM
         call dery (sy7,sy4,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call forcage_flugrange_z(sy3)!IBM
         call derz (sy9,sy3,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call forcage_flugrange_x(sy6)!IBM
         call derx (sy1,sy6,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(sy6)!IBM
         call dery (sy4,sy6,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call forcage_flugrange_x(sy8)!IBM
         call derx (sy3,sy8,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=sy2(ijk,1,1)+sy4(ijk,1,1)
            sy7(ijk,1,1)=sy7(ijk,1,1)+sy1(ijk,1,1)
            sy9(ijk,1,1)=sy9(ijk,1,1)+sy3(ijk,1,1)
         enddo
         call forcage_flugrange_z(sy8)!IBM
         call derz (sy3,sy8,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call forcage_flugrange_y(sy5)!IBM
         call dery (sy1,sy5,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call forcage_flugrange_z(sy5)!IBM
         call derz (sy4,sy5,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=sy2(ijk,1,1)+sy3(ijk,1,1)
            sy7(ijk,1,1)=sy7(ijk,1,1)+sy4(ijk,1,1)
            sy9(ijk,1,1)=sy9(ijk,1,1)+sy1(ijk,1,1)
         enddo
         call forcage_flugrange_x(ux)!IBM
         call derx (sy1,ux,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(uy)!IBM
         call dery (sy4,uy,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call forcage_flugrange_z(uz)!IBM
         call derz (sy3,uz,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call forcage_flugrange_x(uy)!IBM
         call derx (sy6,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(uz)!IBM
         call dery (sy5,uz,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call forcage_flugrange_z(ux)!IBM
         call derz (sy8,ux,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=0.5*(sy2(ijk,1,1)+ux(ijk,1,1)*sy1(ijk,1,1)+uz(ijk,1,1)*sy8(ijk,1,1))
            sy7(ijk,1,1)=0.5*(sy7(ijk,1,1)+uy(ijk,1,1)*sy4(ijk,1,1)+ux(ijk,1,1)*sy6(ijk,1,1))
            sy9(ijk,1,1)=0.5*(sy9(ijk,1,1)+uz(ijk,1,1)*sy3(ijk,1,1)+uy(ijk,1,1)*sy5(ijk,1,1))
         enddo
!
         call forcage_flugrange_x(uz)!IBM
         call derx (sy3,uz,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(ux)!IBM
         call dery (sy1,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call forcage_flugrange_z(uy)!IBM
         call derz (sy4,uy,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=sy2(ijk,1,1)+0.5*(uy(ijk,1,1)*sy1(ijk,1,1))
            sy7(ijk,1,1)=sy7(ijk,1,1)+0.5*(uz(ijk,1,1)*sy4(ijk,1,1))
            sy9(ijk,1,1)=sy9(ijk,1,1)+0.5*(ux(ijk,1,1)*sy3(ijk,1,1))
         enddo
!
      else!2D
         do ijk=1,nxyz1
            sy1(ijk,1,1)=ux(ijk,1,1)*ux(ijk,1,1)
            sy4(ijk,1,1)=uy(ijk,1,1)*uy(ijk,1,1)
            sy3(ijk,1,1)=ux(ijk,1,1)*uy(ijk,1,1)
         enddo
         call forcage_flugrange_x(sy1)!IBM
         call derx (sy2,sy1,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(sy4)!IBM
         call dery (sy7,sy4,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call forcage_flugrange_x(sy3)!IBM
         call derx (sy1,sy3,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(sy3)!IBM
         call dery (sy4,sy3,di1,di2,sy,ffy,fsy,fwy,ppy,nx,ny,nz,0)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=sy2(ijk,1,1)+sy4(ijk,1,1)
            sy7(ijk,1,1)=sy7(ijk,1,1)+sy1(ijk,1,1)
         enddo
         call forcage_flugrange_x(ux)!IBM
         call derx (sy1,ux,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(uy)!IBM
         call dery (sy4,uy,di1,di2,sy,ffy,fsy,fwy,ppy,nx,ny,nz,0)
         call forcage_flugrange_x(uy)!IBM
         call derx (sy3,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call forcage_flugrange_y(ux)!IBM
         call dery (sy5,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=0.5*(sy2(ijk,1,1)+ux(ijk,1,1)*sy1(ijk,1,1)+uy(ijk,1,1)*sy5(ijk,1,1))
            sy7(ijk,1,1)=0.5*(sy7(ijk,1,1)+uy(ijk,1,1)*sy4(ijk,1,1)+ux(ijk,1,1)*sy3(ijk,1,1))
         enddo
      endif
   endif
   if (ifiltre==1)then!filtrage
      if (nz.gt.1) then!3D
         call filx(sy1,sy2,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,0)
         call filx(sy3,sy7,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,1)
         call filx(sy4,sy9,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,1)
         if (ncly.eq.0) then
            call fily(sy5,sy1,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,1)
            call fily(sy6,sy3,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,0)
            call fily(sy8,sy4,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,1)
         endif
         if ((ncly.eq.1).or.(ncly.eq.2)) then 
            call fily(sy5,sy1,di1,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,filayp,fiz1y,fiz2y,nx,ny,nz,1)
            call fily(sy6,sy3,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,0)
            call fily(sy8,sy4,di1,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,filayp,fiz1y,fiz2y,nx,ny,nz,1)
         endif
         if (nclz.eq.0) then
            call filz(sy2,sy5,di1,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nx,ny,nz,1)
            call filz(sy7,sy6,di1,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nx,ny,nz,1)
            call filz(sy9,sy8,di1,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nx,ny,nz,0)
         endif
         if ((nclz.eq.1).or.(nclz.eq.2)) then
            call filz(sy2,sy5,di1,sz,vz,fiffzp,fifzp,ficzp,fibzp,fibbzp,filazp,fiz1z,fiz2z,nx,ny,nz,1)
            call filz(sy7,sy6,di1,sz,vz,fiffzp,fifzp,ficzp,fibzp,fibbzp,filazp,fiz1z,fiz2z,nx,ny,nz,1)
            call filz(sy9,sy8,di1,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nx,ny,nz,0)
         endif
      else!2D
         call filx(sy1,sy2,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,0)
         call filx(sy3,sy7,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,1)
         if (ncly.eq.0) then   
            call fily(sy2,sy1,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,1)
            call fily(sy7,sy3,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,0)
         endif
         if ((ncly.eq.1).or.(ncly.eq.2)) then
            call fily(sy2,sy1,di1,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,filayp,fiz1y,fiz2y,nx,ny,nz,1)
            call fily(sy7,sy3,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,0)
         endif
      endif
   endif
   sy8(:,:,:)=sy7(:,:,:)
   call forcage_flugrange_x(ux)!IBM
   call derxx (sy1,ux,di1,sx,sfx ,ssx ,swx ,nx,ny,nz,0)
   if (istret.ne.0) then !raffinement
      call forcage_flugrange_y(ux)!IBM
      call deryy (sy4,ux,di1,di2,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
      call forcage_flugrange_y(ux)!IBM
      call dery (sy5,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         sy4(i,j,k)=sy4(i,j,k)*pp2y(j)-pp4y(j)*sy5(i,j,k)
      enddo
      enddo
      enddo
   else! pas de raffinement
      call forcage_flugrange_y(ux)!IBM
      call deryy (sy4,ux,di1,di2,sy,sfyp,ssyp,swyp,nx,ny,nz,1) 
   endif
   if (nz.gt.1) then!3D
      call forcage_flugrange_z(ux)!IBM
      call derzz (sy6,ux,di1,sz,sfzp,sszp,swzp,nx,ny,nz,1)
      do ijk=1,nxyz1
         sy7(ijk,1,1)=xnu*(sy1(ijk,1,1)+sy4(ijk,1,1)+sy6(ijk,1,1))-sy2(ijk,1,1)
      enddo
   else!2D
      do j=1,ny
      do i=1,nx
         sy7(i,j,1)=xnu*(sy1(i,j,1)+sy4(i,j,1))-sy2(i,j,1)
      enddo
      enddo
   endif
!
   call forcage_flugrange_x(uy)!IBM
   call derxx (sy1,uy,di1,sx,sfxp,ssxp,swxp,nx,ny,nz,1)
   if (istret.ne.0) then !raffinement
      call forcage_flugrange_y(uy)!IBM
      call deryy (sy4,uy,di1,di2,sy,sfy,ssy,swy,nx,ny,nz,0)
      call forcage_flugrange_y(uy)!IBM
      call dery (sy5,uy,di1,di2,sy,ffy,fsy,fwy,ppy,nx,ny,nz,0)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         sy4(i,j,k)=sy4(i,j,k)*pp2y(j)-pp4y(j)*sy5(i,j,k)
      enddo
      enddo
      enddo
   else!pas de raffinement
      call forcage_flugrange_y(uy)!IBM
      call deryy (sy4,uy,di1,di2,sy,sfy,ssy,swy,nx,ny,nz,0)
   endif
   if (nz.gt.1) then!3D
      call forcage_flugrange_z(uy)!IBM
      call derzz (sy6,uy,di1,sz,sfzp,sszp,swzp,nx,ny,nz,1)
      do ijk=1,nxyz1
         sy8(ijk,1,1)=xnu*(sy1(ijk,1,1)+sy4(ijk,1,1)+sy6(ijk,1,1))-sy8(ijk,1,1)
      enddo
   else!2D
      do j=1,ny
      do i=1,nx
         sy8(i,j,1)=xnu*(sy1(i,j,1)+sy4(i,j,1))-sy8(i,j,1)
      enddo
      enddo
   endif
!
   if (nz.gt.1) then!3D
      call forcage_flugrange_x(uz)!IBM
      call derxx (sy1,uz,di1,sx,sfxp,ssxp,swxp,nx,ny,nz,1)
      if (istret.ne.0) then!raffinement
         call forcage_flugrange_y(uz)!IBM
         call deryy (sy4,uz,di1,di2,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
         call forcage_flugrange_y(uz)!IBM
         call dery (sy5,uz,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         do k=1,nz
         do j=1,ny
         do i=1,nx
            sy4(i,j,k)=sy4(i,j,k)*pp2y(j)-pp4y(j)*sy5(i,j,k)
         enddo
         enddo
         enddo
      else!pas de raffinement
         call forcage_flugrange_y(uz)!IBM
         call deryy (sy4,uz,di1,di2,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
      endif
      call forcage_flugrange_z(uz)!IBM
      call derzz (sy6,uz,di1,sz,sfz ,ssz ,swz ,nx,ny,nz,0)
      do ijk=1,nxyz1
         sy9(ijk,1,1)=xnu*(sy1(ijk,1,1)+sy4(ijk,1,1)+sy6(ijk,1,1))-sy9(ijk,1,1)
      enddo
   endif
  else
   if (iskew==0) then !UROTU!
      if (nz.gt.1) then!3D
         call dery (sy1,uz,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call derz (sy2,uy,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call derz (sy3,ux,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call derx (sy7,uz,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call derx (sy8,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)          
         call dery (sy9,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy1(ijk,1,1)=sy1(ijk,1,1)-sy2(ijk,1,1)
            sy3(ijk,1,1)=sy3(ijk,1,1)-sy7(ijk,1,1)
            sy8(ijk,1,1)=sy8(ijk,1,1)-sy9(ijk,1,1)
         enddo
         do ijk=1,nxyz1
            sy2(ijk,1,1)=sy3(ijk,1,1)*uz(ijk,1,1)-sy8(ijk,1,1)*uy(ijk,1,1)
            sy7(ijk,1,1)=sy8(ijk,1,1)*ux(ijk,1,1)-sy1(ijk,1,1)*uz(ijk,1,1)
            sy9(ijk,1,1)=sy1(ijk,1,1)*uy(ijk,1,1)-sy3(ijk,1,1)*ux(ijk,1,1)
         enddo
      else!2D
         call derx (sy8,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy9,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1) 
         do ijk=1,nxyz1
            sy8(ijk,1,1)= sy8(ijk,1,1)-sy9(ijk,1,1)
            sy2(ijk,1,1)=-sy8(ijk,1,1)*uy(ijk,1,1)
            sy7(ijk,1,1)= sy8(ijk,1,1)*ux(ijk,1,1)
         enddo
      endif
   else!skew
      if (nz.gt.1) then!3D
         do ijk=1,nxyz1
            sy1(ijk,1,1)=ux(ijk,1,1)*ux(ijk,1,1)
            sy4(ijk,1,1)=uy(ijk,1,1)*uy(ijk,1,1)
            sy3(ijk,1,1)=uz(ijk,1,1)*uz(ijk,1,1)
            sy5(ijk,1,1)=uy(ijk,1,1)*uz(ijk,1,1)
            sy8(ijk,1,1)=ux(ijk,1,1)*uz(ijk,1,1)
            sy6(ijk,1,1)=ux(ijk,1,1)*uy(ijk,1,1)
         enddo
         call derx (sy2,sy1,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy7,sy4,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call derz (sy9,sy3,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call derx (sy1,sy6,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy4,sy6,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call derx (sy3,sy8,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=sy2(ijk,1,1)+sy4(ijk,1,1)
            sy7(ijk,1,1)=sy7(ijk,1,1)+sy1(ijk,1,1)
            sy9(ijk,1,1)=sy9(ijk,1,1)+sy3(ijk,1,1)
         enddo
         call derz (sy3,sy8,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call dery (sy1,sy5,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call derz (sy4,sy5,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=sy2(ijk,1,1)+sy3(ijk,1,1)
            sy7(ijk,1,1)=sy7(ijk,1,1)+sy4(ijk,1,1)
            sy9(ijk,1,1)=sy9(ijk,1,1)+sy1(ijk,1,1)
         enddo
         call derx (sy1,ux,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy4,uy,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call derz (sy3,uz,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         call derx (sy6,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy5,uz,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call derz (sy8,ux,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=0.5*(sy2(ijk,1,1)+ux(ijk,1,1)*sy1(ijk,1,1)+uz(ijk,1,1)*sy8(ijk,1,1))
            sy7(ijk,1,1)=0.5*(sy7(ijk,1,1)+uy(ijk,1,1)*sy4(ijk,1,1)+ux(ijk,1,1)*sy6(ijk,1,1))
            sy9(ijk,1,1)=0.5*(sy9(ijk,1,1)+uz(ijk,1,1)*sy3(ijk,1,1)+uy(ijk,1,1)*sy5(ijk,1,1))
         enddo
!
         call derx (sy3,uz,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy1,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call derz (sy4,uy,di1,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         do ijk=1,nxyz1
            sy2(ijk,1,1)=sy2(ijk,1,1)+0.5*(uy(ijk,1,1)*sy1(ijk,1,1))
            sy7(ijk,1,1)=sy7(ijk,1,1)+0.5*(uz(ijk,1,1)*sy4(ijk,1,1))
            sy9(ijk,1,1)=sy9(ijk,1,1)+0.5*(ux(ijk,1,1)*sy3(ijk,1,1))
         enddo
!
      else!2D
!         do ijk=1,nxyz1 ! TODO
         do k=1,nz
         do j=1,ny
         do i=1,nx
            sy1(i,j,k)=ux(i,j,k)*ux(i,j,k)
            sy4(i,j,k)=uy(i,j,k)*uy(i,j,k)
            sy3(i,j,k)=ux(i,j,k)*uy(i,j,k)
         enddo
         enddo
         enddo
         call derx (sy2,sy1,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy7,sy4,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         call derx (sy1,sy3,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy4,sy3,di1,di2,sy,ffy,fsy,fwy,ppy,nx,ny,nz,0)
!         do ijk=1,nxyz1 ! TODO
         do k=1,nz
         do j=1,ny
         do i=1,nx
            sy2(i,j,k)=sy2(i,j,k)+sy4(i,j,k)
            sy7(i,j,k)=sy7(i,j,k)+sy1(i,j,k)
         enddo
         enddo
         enddo
         call derx (sy1,ux,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy4,uy,di1,di2,sy,ffy,fsy,fwy,ppy,nx,ny,nz,0)
         call derx (sy3,uy,di1,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (sy5,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
!         do ijk=1,nxyz1 ! TODO
         do k=1,nz
         do j=1,ny
         do i=1,nx
            sy2(i,j,k)=0.5*(sy2(i,j,k)+ux(i,j,k)*sy1(i,j,k)+uy(i,j,k)*sy5(i,j,k))
            sy7(i,j,k)=0.5*(sy7(i,j,k)+uy(i,j,k)*sy4(i,j,k)+ux(i,j,k)*sy3(i,j,k))
         enddo
         enddo
         enddo
      endif
   endif
   if (ifiltre==1)then!filtrage
      if (nz.gt.1) then!3D
         call filx(sy1,sy2,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,0)
         call filx(sy3,sy7,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,1)
         call filx(sy4,sy9,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,1)
         if (ncly.eq.0) then
            call fily(sy5,sy1,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,1)
            call fily(sy6,sy3,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,0)
            call fily(sy8,sy4,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,1)
         endif
         if ((ncly.eq.1).or.(ncly.eq.2)) then 
            call fily(sy5,sy1,di1,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,filayp,fiz1y,fiz2y,nx,ny,nz,1)
            call fily(sy6,sy3,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,0)
            call fily(sy8,sy4,di1,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,filayp,fiz1y,fiz2y,nx,ny,nz,1)
         endif
         if (nclz.eq.0) then
            call filz(sy2,sy5,di1,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nx,ny,nz,1)
            call filz(sy7,sy6,di1,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nx,ny,nz,1)
            call filz(sy9,sy8,di1,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nx,ny,nz,0)
         endif
         if ((nclz.eq.1).or.(nclz.eq.2)) then
            call filz(sy2,sy5,di1,sz,vz,fiffzp,fifzp,ficzp,fibzp,fibbzp,filazp,fiz1z,fiz2z,nx,ny,nz,1)
            call filz(sy7,sy6,di1,sz,vz,fiffzp,fifzp,ficzp,fibzp,fibbzp,filazp,fiz1z,fiz2z,nx,ny,nz,1)
            call filz(sy9,sy8,di1,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nx,ny,nz,0)
         endif
      else!2D
         call filx(sy1,sy2,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,0)
         call filx(sy3,sy7,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx,ny,nz,1)
         if (ncly.eq.0) then   
            call fily(sy2,sy1,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,1)
            call fily(sy7,sy3,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,0)
         endif
         if ((ncly.eq.1).or.(ncly.eq.2)) then
            call fily(sy2,sy1,di1,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,filayp,fiz1y,fiz2y,nx,ny,nz,1)
            call fily(sy7,sy3,di1,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,nx,ny,nz,0)
         endif
      endif
   endif
   sy8(:,:,:)=sy7(:,:,:)
   call derxx (sy1,ux,di1,sx,sfx ,ssx ,swx ,nx,ny,nz,0)
   if (istret.ne.0) then !raffinement
      call deryy (sy4,ux,di1,di2,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
      call dery (sy5,ux,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         sy4(i,j,k)=sy4(i,j,k)*pp2y(j)-pp4y(j)*sy5(i,j,k)
      enddo
      enddo
      enddo
   else! pas de raffinement
      call deryy (sy4,ux,di1,di2,sy,sfyp,ssyp,swyp,nx,ny,nz,1) 
   endif
   if (nz.gt.1) then!3D
      call derzz (sy6,ux,di1,sz,sfzp,sszp,swzp,nx,ny,nz,1)
      do ijk=1,nxyz1
         sy7(ijk,1,1)=xnu*(sy1(ijk,1,1)+sy4(ijk,1,1)+sy6(ijk,1,1))-sy2(ijk,1,1)
      enddo
   else!2D
      do j=1,ny
      do i=1,nx
         sy7(i,j,1)=xnu*(sy1(i,j,1)+sy4(i,j,1))-sy2(i,j,1)
      enddo
      enddo
   endif
!
   call derxx (sy1,uy,di1,sx,sfxp,ssxp,swxp,nx,ny,nz,1)
   if (istret.ne.0) then !raffinement
      call deryy (sy4,uy,di1,di2,sy,sfy,ssy,swy,nx,ny,nz,0)
      call dery (sy5,uy,di1,di2,sy,ffy,fsy,fwy,ppy,nx,ny,nz,0)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         sy4(i,j,k)=sy4(i,j,k)*pp2y(j)-pp4y(j)*sy5(i,j,k)
      enddo
      enddo
      enddo
   else!pas de raffinement
      call deryy (sy4,uy,di1,di2,sy,sfy,ssy,swy,nx,ny,nz,0)
   endif
   if (nz.gt.1) then!3D
      call derzz (sy6,uy,di1,sz,sfzp,sszp,swzp,nx,ny,nz,1)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         sy8(ijk,1,1)=xnu*(sy1(ijk,1,1)+sy4(ijk,1,1)+sy6(ijk,1,1))-sy8(ijk,1,1)
      enddo
      enddo
      enddo
   else!2D
      do j=1,ny
      do i=1,nx
         sy8(i,j,1)=xnu*(sy1(i,j,1)+sy4(i,j,1))-sy8(i,j,1)
      enddo
      enddo
   endif
!
   if (nz.gt.1) then!3D
      call derxx (sy1,uz,di1,sx,sfxp,ssxp,swxp,nx,ny,nz,1)
      if (istret.ne.0) then!raffinement
         call deryy (sy4,uz,di1,di2,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
         call dery (sy5,uz,di1,di2,sy,ffyp,fsyp,fwyp,ppy,nx,ny,nz,1)
         do k=1,nz
         do j=1,ny
         do i=1,nx
            sy4(i,j,k)=sy4(i,j,k)*pp2y(j)-pp4y(j)*sy5(i,j,k)
         enddo
         enddo
         enddo
      else!pas de raffinement
         call deryy (sy4,uz,di1,di2,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
      endif
      call derzz (sy6,uz,di1,sz,sfz ,ssz ,swz ,nx,ny,nz,0)
      do ijk=1,nxyz1
         sy9(ijk,1,1)=xnu*(sy1(ijk,1,1)+sy4(ijk,1,1)+sy6(ijk,1,1))-sy9(ijk,1,1)
      enddo
   endif
  endif
return
end subroutine convdiff
end module convdiff_m