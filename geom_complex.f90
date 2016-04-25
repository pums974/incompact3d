!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine gene_epsi(epsi,epsidec)
!
USE param
USE variables
   implicit none
!
   real(8),dimension(nx,ny,nz)    :: epsi
   real(8),dimension(nxm,nym,nzm) :: epsidec
!
   integer                        :: i,j,k
   real                           :: dxfin,dyfin,dzfin
   real(8)                        :: x,y,z
   real(8)                        :: xc,yc,r
   real(8),dimension(nxfin,ny,nz) :: xepsi
   real(8),dimension(nx,nyfin,nz) :: yepsi
   real(8),dimension(nx,ny,nzfin) :: zepsi
   integer                        :: inum,jnum,knum
   integer                        :: nobjxmax,nobjymax,nobjzmax
   integer                        :: iexa
!
!  générer paramètres (2):
   xc=2.25
   yc=0.75
   r=0.25
   iexa=0        !=0 : valeurs approchées,=1 : valeurs exactes
!
!  générer 'epsi' :
   epsi(:,:,:)=0.
   call disque_2d(epsi,nx,ny,dx,dy,r,xc,yc,1.)

!  générer 'xepsi' :
   if(nclx.eq.0)then
      dxfin=xlx/nxfin
   elseif(nclx.eq.1.or.nclx.eq.2)then
      dxfin=xlx/(nxfin-1)
   endif
!
   xepsi(:,:,:)=0.
   call disque_2d(xepsi,nxfin,ny,dxfin,dy,r,xc,yc,1.)

!  générer 'yepsi' :
   if(ncly.eq.0)then
      dyfin=yly/nyfin
   elseif(ncly.eq.1.or.ncly.eq.2)then
      dyfin=yly/(nyfin-1)
   endif
!
   yepsi(:,:,:)=0.
   call disque_2d(yepsi,nx,nyfin,dx,dyfin,r,xc,yc,1.)
!
!  générer 'zepsi' (3D) :
   if(nz.gt.1)then
!
      zepsi(:,:,:)=0.
      if(nclz.eq.0)then
         dzfin=zlz/nzfin
      elseif(nclz.eq.1.or.nclz.eq.2)then
         dzfin=zlz/(nzfin-1)
      endif
!
      do k=1,nzfin
      do j=1,ny
      do i=1,nx
         z=dzfin*(k-1)
         y=yp(j)
         x=dx*(i-1)
         if(sqrt((x-xc)**2+(y-yc)**2).le.r)then
            zepsi(i,j,k)=1.
         endif
      enddo
      enddo
      enddo
!
   endif
!
!  Rechercher le nombre max d'obstacles par tranche 1D :
!  selon x :
   nobjx(:,:)=0
   nobjxmax=0
   do k=1,nz
   do j=1,ny
      inum=0
      if(epsi(1,j,k).eq.1.)then
         inum=1
         nobjx(j,k)=1
      endif
      do i=1,nx-1
         if(epsi(i,j,k).eq.0..and.epsi(i+1,j,k).eq.1.)then
            inum=inum+1
            nobjx(j,k)=nobjx(j,k)+1
         endif
      enddo
      if(inum.gt.nobjxmax)then
         nobjxmax=inum
      endif
   enddo
   enddo
print*,'nobjxmax=',nobjxmax
!
!  selon y :
   nobjy(:,:)=0
   nobjymax=0
   do k=1,nz
   do i=1,nx
      jnum=0
      if(epsi(i,1,k).eq.1.)then
         jnum=1
         nobjy(i,k)=1
      endif
      do j=1,ny-1
         if(epsi(i,j,k).eq.0..and.epsi(i,j+1,k).eq.1.)then
            jnum=jnum+1
            nobjy(i,k)=nobjy(i,k)+1
         endif
      enddo
      if(jnum.gt.nobjymax)then
         nobjymax=jnum
      endif
   enddo
   enddo
print*,'nobjymax=',nobjymax
!
!  selon z :
   nobjz(:,:)=0
   nobjzmax=0
   do j=1,ny
   do i=1,nx
      knum=0
      if(epsi(i,j,1).eq.1.)then
         knum=1
         nobjz(i,j)=1
      endif
      do k=1,nz-1
         if(epsi(i,j,k).eq.0..and.epsi(i,j,k+1).eq.1.)then
            knum=knum+1
            nobjz(i,j)=nobjz(i,j)+1
         endif
      enddo
      if(knum.gt.nobjzmax)then
         nobjzmax=knum
      endif
   enddo
   enddo
print*,'nobjzmax=',nobjzmax
!
!  Valeurs exactes des points d'adhérence (cylindre) :
   if(iexa.eq.1)then
!
!  selon x:
      do k=1,nz
      do j=1,ny
      y=yp(j)
      inum=0
      do i=1,nx-1
         x=sqrt(r*r-(y-yc)*(y-yc))
         if(epsi(i,j,k).eq.0..and.epsi(i+1,j,k).eq.1.)then
            inum=inum+1
            xi(inum,j,k)=xc-x
         elseif(epsi(i,j,k).eq.1..and.epsi(i+1,j,k).eq.0.)then
            xf(inum,j,k)=xc+x
         endif
      enddo
      enddo
      enddo
!
!  selon y:
      do k=1,nz
      do i=1,nx
      x=dx*(i-1)
      jnum=0
      do j=1,ny-1
         y=sqrt(r*r-(x-xc)*(x-xc))
         if(epsi(i,j,k).eq.0..and.epsi(i,j+1,k).eq.1.)then
            jnum=jnum+1
            yi(jnum,i,k)=yc-y
         elseif(epsi(i,j,k).eq.1..and.epsi(i,j+1,k).eq.0.)then
            yf(jnum,i,k)=yc+y
         endif
      enddo
      enddo
      enddo
!
!  Rechercher les positions approchées des points d'adhérence :
   else
!
!  selon x :
      do k=1,nz
      do j=1,ny
         inum=0
         if(xepsi(1,j,k).eq.1.)then
            inum=inum+1
            xi(inum,j,k)=-dx/2.
         endif
         do i=1,nxfin-1
            if(xepsi(i,j,k).eq.0..and.xepsi(i+1,j,k).eq.1.)then
               inum=inum+1
               xi(inum,j,k)=dxfin*(i-1)+dxfin/2.
            elseif(xepsi(i,j,k).eq.1..and.xepsi(i+1,j,k).eq.0.)then
               xf(inum,j,k)=dxfin*(i-1)+dxfin/2.
            endif
         enddo
         if(xepsi(nxfin,j,k).eq.1.)then
            xf(inum,j,k)=xlx+dx/2.
         endif
      enddo
      enddo
!
!  selon y :
      do k=1,nz
      do i=1,nx
         jnum=0
         if(yepsi(i,1,k).eq.1.)then
            jnum=jnum+1
            yi(jnum,i,k)=-dy/2.
         endif
         do j=1,nyfin-1
            if(yepsi(i,j,k).eq.0..and.yepsi(i,j+1,k).eq.1.)then
               jnum=jnum+1
               yi(jnum,i,k)=dyfin*(j-1)+dyfin/2.
            elseif(yepsi(i,j,k).eq.1..and.yepsi(i,j+1,k).eq.0.)then
               yf(jnum,i,k)=dyfin*(j-1)+dyfin/2.
            endif
         enddo
         if(yepsi(i,nyfin,k).eq.1.)then
            yf(jnum,i,k)=yly+dy/2.
         endif
      enddo
      enddo
!
!  selon z :
      do j=1,ny
      do i=1,nx
         knum=0
         if(zepsi(i,j,1).eq.1.)then
            knum=knum+1
            zi(knum,i,j)=-dz/2.
         endif
         do k=1,nzfin-1
            if(zepsi(i,j,k).eq.0..and.zepsi(i,j,k+1).eq.1.)then
               knum=knum+1
               zi(knum,i,j)=dzfin*(k-1)+dzfin/2.
            elseif(zepsi(i,j,k).eq.1..and.zepsi(i,j,k+1).eq.0.)then
               zf(knum,i,j)=dzfin*(k-1)+dzfin/2.
            endif
         enddo
         if(zepsi(i,j,nzfin).eq.1.)then
            zf(knum,i,j)=zlz+dz/2.
         endif
      enddo
      enddo
!
   endif
!
   return
end subroutine gene_epsi
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine verif_epsi(epsi)
!
USE param
USE variables
   implicit none
!
   real(8),dimension(nx,ny,nz) :: epsi
!
   integer                     :: i,j,k
!
!  Repérer les singularités à éviter sur epsilon:
!  balayage sur x :
   do k=1,nz
   do j=1,ny
!
      if(epsi(1,j,k).eq.0..and.&
         epsi(2,j,k).eq.1..or. &
         epsi(1,j,k).eq.0..and.&
         epsi(2,j,k).eq.0..and.&
         epsi(3,j,k).eq.1.)then
         write(*,*)'probleme de singularite sur epsilon (x1)'
         stop
      endif
!
      if(epsi(nx  ,j,k).eq.0..and.&
         epsi(nx-1,j,k).eq.1..or. &
         epsi(nx  ,j,k).eq.0..and.&
         epsi(nx-1,j,k).eq.0..and.&
         epsi(nx-2,j,k).eq.1.)then
         write(*,*)'probleme de singularite sur epsilon (x2)'
         stop
      endif
!
      do i=2,nx-1
         if(epsi(i-1,j,k).eq.1..and.&
            epsi(i  ,j,k).eq.0..and.&
            epsi(i+1,j,k).eq.1.)then
            write(*,*)'probleme de singularite sur epsilon (x3)'
            stop
         endif
      enddo
!
      do i=2,nx-2
         if(epsi(i-1,j,k).eq.1..and.&
            epsi(i  ,j,k).eq.0..and.&
            epsi(i+1,j,k).eq.0..and.&
            epsi(i+2,j,k).eq.1.)then
            write(*,*)'probleme de singularite sur epsilon (x4)'
            stop
         endif
      enddo
!
   enddo
   enddo
!
!  balayage sur y :
   do k=1,nz
   do i=1,nx
!
      if(epsi(i,1,k).eq.0..and.&
         epsi(i,2,k).eq.1..or. &
         epsi(i,1,k).eq.0..and.&
         epsi(i,2,k).eq.0..and.&
         epsi(i,3,k).eq.1.)then
         write(*,*)'probleme de singularite sur epsilon (y1)'
         stop
      endif
!
      if(epsi(i,ny  ,k).eq.0..and.&
         epsi(i,ny-1,k).eq.1..or. &
         epsi(i,ny  ,k).eq.0..and.&
         epsi(i,ny-1,k).eq.0..and.&
         epsi(i,ny-2,k).eq.1.)then
         write(*,*)'probleme de singularite sur epsilon (y2)'
         stop
      endif
!
      do j=2,ny-1
         if(epsi(i,j-1,k).eq.1..and.&
            epsi(i,j  ,k).eq.0..and.&
            epsi(i,j+1,k).eq.1.)then
            write(*,*)'probleme de singularite sur epsilon (y3)'
            print*,dx*(i-1),dy*(j-1)
            stop
         endif
      enddo
!
      do j=2,ny-2
         if(epsi(i,j-1,k).eq.1..and.&
            epsi(i,j  ,k).eq.0..and.&
            epsi(i,j+1,k).eq.0..and.&
            epsi(i,j+2,k).eq.1.)then
            write(*,*)'probleme de singularite sur epsilon (y4)'
            stop
         endif
      enddo
!
   enddo
   enddo
!
!  balayage sur z :
   if(nz.gt.1)then
      do j=1,ny
      do i=1,nx
!
         if(epsi(i,j,1).eq.0..and.&
            epsi(i,j,2).eq.1..or. &
            epsi(i,j,1).eq.0..and.&
            epsi(i,j,2).eq.0..and.&
            epsi(i,j,3).eq.1.)then
            write(*,*)'probleme de singularite sur epsilon (z1)'
            stop
         endif
!
         if(epsi(i,j,nz  ).eq.0..and.&
            epsi(i,j,nz-1).eq.1..or. &
            epsi(i,j,nz  ).eq.0..and.&
            epsi(i,j,nz-1).eq.0..and.&
            epsi(i,j,nz-2).eq.1.)then
            write(*,*)'probleme de singularite sur epsilon (z2)'
            stop
         endif
!
         do k=2,nz-1
            if(epsi(i,j,k-1).eq.1..and.&
               epsi(i,j,k  ).eq.0..and.&
               epsi(i,j,k+1).eq.1.)then
               write(*,*)'probleme de singularite sur epsilon (z3)'
               stop
            endif
         enddo
!
         do k=2,nz-2
            if(epsi(i,j,k-1).eq.1..and.&
               epsi(i,j,k  ).eq.0..and.&
               epsi(i,j,k+1).eq.0..and.&
               epsi(i,j,k+2).eq.1.)then
               write(*,*)'probleme de singularite sur epsilon (z4)'
               stop
            endif
         enddo
!
      enddo
      enddo
   endif  
!
   return
end subroutine verif_epsi
!
!
!---------------------------------------------------------------------------
!*      *      *      *      *      *      *      *      *      *      *    
!*   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *  
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
!************* subroutines de forçage / polynôme de Lagrange ***************
!*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!*   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *  
!*      *      *      *      *      *      *      *      *      *      *    
!---------------------------------------------------------------------------
!
!
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine forcage_flugrange_x(u)
!
USE param
USE variables
   implicit none
!
   real(8),dimension(nx,ny,nz) :: u
!
   integer                     :: i,j,k
   real(8)                     :: x,y,z
   integer                     :: ix          != position du point "zappé"
   integer                     :: ipif,ipol
   integer                     :: npif        != Nombre de Points d'Influence Fluides par frontière
   integer                     :: ipoli,ipolf != positions Initiales et Finales du POLynôme considéré
!
   real(8)                  :: xpol,ypol,dypol !|variables concernant les polynômes
   real(8),dimension(10)    :: xa,ya           !|de Lagrange. A mettre imérativement en 
   integer                  :: ia,na           !|double précision
!
!  parametres :
   npif=2
!
!  double-boucle (y,z) :
   do k=1,nz
   do j=1,ny
      z=dz*(k-1)
      y=dy*(j-1)
      if(nobjx(j,k).ne.0)then
         ia=0
         do i=1,nobjx(j,k)          !boucle sur le nombre d'objets par (j,k)
         !  1ère frontière
            ia=ia+1
            xa(ia)=xi(i,j,k)
            ya(ia)=0.
            if(xi(i,j,k).gt.0.)then !pt.d'inf.flu. (objet immergé)
               ix=xi(i,j,k)/dx+1
               ipoli=ix+1
               do ipif=1,npif
                  ia=ia+1
                  xa(ia)=(ix-1)*dx-ipif*dx
                  ya(ia)=u(ix-ipif,j,k)
!                  xa(ia)=(ix-1)*dx-(ipif-1)*dx !no zap
!                  ya(ia)=u(ix-ipif+1,j,k)      !no zap
               enddo
            else                    !pt.d'inf.flu. (frontière = bordure du domaine)
               ipoli=1
               do ipif=1,npif
                  ia=ia+1
                  xa(ia)=-ipif*dx
                  ya(ia)=0.
               enddo
            endif
         !  2ème frontière
            ia=ia+1
            xa(ia)=xf(i,j,k)
            ya(ia)=0.
            if(xf(i,j,k).lt.xlx)then !objet immergé
               ix=(xf(i,j,k)+dx)/dx+1
               ipolf=ix-1
               do ipif=1,npif
                  ia=ia+1
                  xa(ia)=(ix-1)*dx+ipif*dx
                  ya(ia)=u(ix+ipif,j,k)
!                  xa(ia)=(ix-1)*dx+(ipif-1)*dx !no zap
!                  ya(ia)=u(ix+ipif-1,j,k)      !no zap
               enddo
            else                     !frontière = bordure du domaine
               ipolf=nx
               do ipif=1,npif
                  ia=ia+1
                  xa(ia)=xlx+ipif*dx
                  ya(ia)=0.
               enddo
            endif
         !  cas (très) particulier (ne marche que pour les frontières exactes)
            if(xi(i,j,k).eq.xf(i,j,k))then
               u(ipoli-1,j,k)=0.     !ou u(ipolf+1,j,k)=0.
            endif
         !  calcul du polynôme
            na=ia
            do ipol=ipoli,ipolf
               xpol=dx*(ipol-1)
               call polint(xa,ya,na,xpol,ypol,dypol)
               u(ipol,j,k)=ypol
            enddo
            ia=0
         enddo
      endif
   enddo
   enddo
!
   return
end subroutine forcage_flugrange_x
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine forcage_flugrange_y(u)
!
USE param
USE variables
   implicit none
!
   real(8),dimension(nx,ny,nz) :: u
!
   integer                     :: i,j,k
   real(8)                     :: x,y,z
   integer                     :: jy          != position du point "zappé"
   integer                     :: jpif,jpol
   integer                     :: npif        != Nombre de Points d'Influence Fluides par frontière
   integer                     :: jpoli,jpolf != positions Initiales et Finales du POLynôme considéré
!
   real(8)                     :: xpol,ypol,dypol !|variables concernant les polynômes
   real(8),dimension(10)       :: xa,ya           !|de Lagrange. A mettre imérativement en 
   integer                     :: ia,na           !|double précision
!
!  parametres :
   npif=2
!
!  double-boucle (y,z) :
   do k=1,nz
   do i=1,nx
      z=dz*(k-1)
      x=dx*(i-1)
      if(nobjy(i,k).ne.0)then
         ia=0
         do j=1,nobjy(i,k)          !boucle sur le nombre d'objets par (j,k)
         !  1ère frontière
            ia=ia+1
            xa(ia)=yi(j,i,k)
            ya(ia)=0.
            if(yi(j,i,k).gt.0.)then !pt.d'inf.flu. (objet immergé)
               jy=yi(j,i,k)/dy+1
               jpoli=jy+1
               do jpif=1,npif
                  ia=ia+1
                  xa(ia)=(jy-1)*dy-jpif*dy
                  ya(ia)=u(i,jy-jpif,k)
!                  xa(ia)=(jy-1)*dy-(jpif-1)*dy !no zap
!                  ya(ia)=u(i,jy-jpif+1,k)      !no zap
               enddo
            else                    !pt.d'inf.flu. (frontière = bordure du domaine)
               jpoli=1
               do jpif=1,npif
                  ia=ia+1
                  xa(ia)=-jpif*dy
                  ya(ia)=0.
               enddo
            endif
         !  2ème frontière
            ia=ia+1
            xa(ia)=yf(j,i,k)
            ya(ia)=0.
            if(yf(j,i,k).lt.yly)then !objet immergé
               jy=(yf(j,i,k)+dy)/dy+1
               jpolf=jy-1
               do jpif=1,npif
                  ia=ia+1
                  xa(ia)=(jy-1)*dy+jpif*dy
                  ya(ia)=u(i,jy+jpif,k)
!                  xa(ia)=(jy-1)*dy+(jpif-1)*dy !no zap
!                  ya(ia)=u(i,jy+jpif-1,k)      !no zap
               enddo
            else                     !frontière = bordure du domaine
               jpolf=ny
               do jpif=1,npif
                  ia=ia+1
                  xa(ia)=yly+jpif*dy
                  ya(ia)=0.
               enddo
            endif
         !  cas (très) particulier (ne marche que pour les frontières exactes)
            if(yi(j,i,k).eq.yf(j,i,k))then
               u(i,jpoli-1,k)=0.     !ou u(ipolf+1,j,k)=0.
            endif
         !  calcul du polynôme
            na=ia
            do jpol=jpoli,jpolf
               xpol=dy*(jpol-1)
               call polint(xa,ya,na,xpol,ypol,dypol)
               u(i,jpol,k)=ypol
            enddo
            ia=0
         enddo
      endif
   enddo
   enddo
!
   return
end subroutine forcage_flugrange_y
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine forcage_flugrange_z(u)
!
USE param
USE variables
   implicit none
!
   real(8),dimension(nx,ny,nz) :: u
!
   integer                     :: i,j,k
   real(8)                     :: x,y,z
   integer                     :: kz          != position du point "zappé"
   integer                     :: kpif,kpol
   integer                     :: npif        != Nombre de Points d'Influence Fluides par frontière
   integer                     :: kpoli,kpolf != positions Initiales et Finales du POLynôme considéré
!
   real(8)                     :: xpol,ypol,dypol !|variables concernant les polynômes
   real(8),dimension(10)       :: xa,ya           !|de Lagrange. A mettre imérativement en 
   integer                     :: ia,na           !|double précision
!
!  parametres :
   npif=2
!
!  double-boucle (y,z) :
   do j=1,ny
   do i=1,nx
      y=dy*(j-1)
      x=dx*(i-1)
      if(nobjz(i,j).ne.0)then
         ia=0
         do k=1,nobjz(i,j)          !boucle sur le nombre d'objets par couple (i,j)
         !  1ère frontière
            ia=ia+1
            xa(ia)=zi(k,i,j)
            ya(ia)=0.
            if(zi(k,i,j).gt.0.)then !pt.d'inf.flu. (objet immergé)
               kz=zi(k,i,j)/dz+1
               kpoli=kz+1
               do kpif=1,npif
                  ia=ia+1
                  xa(ia)=(kz-1)*dz-kpif*dz
                  ya(ia)=u(i,j,kz-kpif)
               enddo
            else                    !pt.d'inf.flu. (frontière = bordure du domaine)
               kpoli=1
               do kpif=1,npif
                  ia=ia+1
                  xa(ia)=-kpif*dz
                  ya(ia)=0.
               enddo
            endif
         !  2ème frontière
            ia=ia+1
            xa(ia)=zf(k,i,j)
            ya(ia)=0.
            if(zf(k,i,j).lt.zlz)then !objet immergé
               kz=(zf(k,i,j)+dz)/dz+1
               kpolf=kz-1
               do kpif=1,npif
                  ia=ia+1
                  xa(ia)=(kz-1)*dz+kpif*dz
                  ya(ia)=u(i,j,kz+kpif)
               enddo
            else                     !frontière = bordure du domaine
               kpolf=nz
               do kpif=1,npif
                  ia=ia+1
                  xa(ia)=zlz+kpif*dz
                  ya(ia)=0.
               enddo
            endif
         !  cas (très) particulier (ne marche que pour les frontières exactes)
            if(zi(k,i,j).eq.zf(k,i,j))then
               u(i,j,kpoli-1)=0.     !ou u(ipolf+1,j,k)=0.
            endif
         !  calcul du polynôme
            na=ia
            do kpol=kpoli,kpolf
               xpol=dz*(kpol-1)
               call polint(xa,ya,na,xpol,ypol,dypol)
               u(i,j,kpol)=ypol
            enddo
            ia=0
         enddo
      endif
   enddo
   enddo
!
   return
end subroutine forcage_flugrange_z
!
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine polint(xa,ya,n,x,y,dy)
   implicit none
!
   integer,parameter       :: nmax=30
   integer                 :: n,i,m,ns
   real(8)                 :: dy,x,y,den,dif,dift,ho,hp,w
   real(8),dimension(nmax) :: c,d
   real(8),dimension(n)    :: xa,ya
!
   ns=1
   dif=abs(x-xa(1))
!
   do i=1,n
      dift=abs(x-xa(i))
      if(dift.lt.dif)then
         ns=i
         dif=dift
      endif
      c(i)=ya(i)
      d(i)=ya(i)
   enddo
!
   y=ya(ns)
   ns=ns-1
!
   do m=1,n-1
      do i=1,n-m
         ho=xa(i)-x
         hp=xa(i+m)-x
         w=c(i+1)-d(i)
         den=ho-hp
         if(den.eq.0)read(*,*)
         den=w/den
         d(i)=hp*den
         c(i)=ho*den
      enddo
      if (2*ns.lt.n-m)then
         dy=c(ns+1)
      else
         dy=d(ns)
         ns=ns-1
      endif
      y=y+dy
   enddo
!
   return
end subroutine polint
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!                          géométries complexes
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine disque_2d(epsi,nx,ny,dx,dy,r,xc,yc,remp)
   implicit none
!
   integer                    :: nx,ny
   real(8),dimension(nx,ny,1) :: epsi
   real(8)                    :: dx,dy
   real(8)                    :: xc,yc,r
   real(8)                    :: remp
!
   integer                    :: i,j
   real(8)                    :: x,y
!
   do j=1,ny
      y=dy*(j-1)
      do i=1,nx
         x=dx*(i-1)
         if(sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)).le.r)then
            epsi(i,j,1)=remp
         endif
      enddo
   enddo
!
   return
end subroutine disque_2d
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine rectangle_2d(epsi,nx,ny,dx,dy,xri,xrs,yri,yrs,remp)
   implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer                    :: nx,ny
   real(8),dimension(nx,ny,1) :: epsi
   real(8),dimension(ny)      :: yp
   real(8)                    :: dx,dy
   real(8)                    :: xri,xrs
   real(8)                    :: yri,yrs
   real(8)                    :: remp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer                    :: i,j
   real(8)                    :: x,y
!-----------------------------------------------
!
   do j=1,ny
      y=dy*(j-1)
      do i=1,nx
         x=dx*(i-1)
         if(x.ge.xri.and.&
            x.le.xrs.and.&
            y.ge.yri.and.&
            y.le.yrs)then
            epsi(i,j,1)=remp
         endif
      enddo
   enddo
!
   return
end subroutine rectangle_2d
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine quadri_2d(epsi,nx,ny,dx,dy,x1,y1,x2,y2,x3,y3,x4,y4,remp)
   implicit none
!
   integer                    :: nx,ny
   real(8),dimension(nx,ny,1) :: epsi
   real(8)                    :: dx,dy
   real(8)                    :: x1,y1
   real(8)                    :: x2,y2
   real(8)                    :: x3,y3
   real(8)                    :: x4,y4
   real(8)                    :: remp
!
   integer                    :: i,j
   real(8)                    :: x,y
   real(8)                    :: a1,b1
   real(8)                    :: a2,b2
   real(8)                    :: a3,b3
   real(8)                    :: a4,b4
!
   a1=(y2-y1)/(x2-x1)
   b1=y1-a1*x1
   a2=(y2-y3)/(x2-x3)
   b2=y2-a2*x2
   a3=(y3-y4)/(x3-x4)
   b3=y3-a3*b3
   a4=(y1-y4)/(x1-x4)
   b4=y4-a4*b4
!
   do j=1,ny
      y=dy*(j-1)
      do i=1,nx
         x=dx*(i-1)
         if(y.ge.a1*x+b1.and.&
            y.le.a2*x+b2.and.&
            y.le.a3*x+b3.and.&
            y.ge.a4*x+b4)then
            epsi(i,j,1)=remp
         endif
      enddo
   enddo
!
   return
end subroutine quadri_2d
!
!***************************************************************************
!***************************************************************************
!***************************************************************************
!
subroutine ellipse_2d(epsi,nx,ny,dx,dy,xa,ya,xc,yc,remp)
   implicit none
!
   integer                    :: nx,ny
   real(8),dimension(nx,ny,1) :: epsi
   real(8)                    :: dx,dy
   real(8)                    :: xa,ya !longueurs des axes selon x et y
   real(8)                    :: xc,yc !coordonnées du centre de l'ellipse
   real(8)                    :: remp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer                     :: i,j
   real(8)                     :: x,y
!-----------------------------------------------
!
   do j=1,ny
      y=dy*(j-1)
      do i=1,nx
         x=dx*(i-1)
         if(sqrt((x-xc)*(x-xc)/(xa*xa)+(y-yc)*(y-yc)/(ya*ya)).le.1.)then
            epsi(i,j,1)=remp
         endif
      enddo
   enddo
!
   return
end subroutine ellipse_2d
