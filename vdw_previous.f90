!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! This gets the list of the gauss points within a certain
! distance from each gauss point.
SUBROUTINE make_near(vdw1,mesh0,x0,r_c,nskip,nneigh)
USE data_mesh
USE data_vdw
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(mesh) :: mesh0
TYPE(vdw_data):: vdw1
REAL(8) :: x0(3*(mesh0%numnods+mesh0%nedge))
DIMENSION :: xneigh(12,3), neigh(12)


if (nskip.eq.1) then
   write(*,*) ' We make near approximately '
   do ielem=1,mesh0%numele
     do ig=1,vdw1%ngauss_vdw
       ihh=(ielem-1)*vdw1%ngauss_vdw+ig
       do jg=ig+1,vdw1%ngauss_vdw
         jhh=(ielem-1)*vdw1%ngauss_vdw+jg
         if(vdw1%near(ihh,0)+1.gt.nneigh-1) STOP ' Trouble with nneigh'
         vdw1%near(ihh,0)=vdw1%near(ihh,0)+1
         vdw1%near(ihh,vdw1%near(ihh,0))=jhh
       enddo
     enddo
     do inei=1,12
       jnei=mesh0%connect(ielem)%neigh_elem(inei)
!!$       write(*,*) jnei, ' we do its neighboring elements '
       if ((jnei.ne.0).and.(jnei.gt.ielem)) then
         do ig=1,vdw1%ngauss_vdw
           ihh=(ielem-1)*vdw1%ngauss_vdw+ig
           do jg=1,vdw1%ngauss_vdw
             jhh=(jnei-1)*vdw1%ngauss_vdw+jg
             if(vdw1%near(ihh,0)+1.gt.nneigh-1) STOP ' Trouble with nneigh'
             vdw1%near(ihh,0)=vdw1%near(ihh,0)+1
             vdw1%near(ihh,vdw1%near(ihh,0))=jhh
           enddo
!!$           ipos=vdw1%near((ielem-1)*vdw1%ngauss_vdw+ig,0)+1
!!$           vdw1%near((ielem-1)*vdw1%ngauss_vdw+ig,0)= &
!!$           vdw1%near((ielem-1)*vdw1%ngauss_vdw+ig,0) + vdw1%ngauss_vdw
!!$           vdw1%near((ielem-1)*vdw1%ngauss_vdw+ig,ipos:ipos+vdw1%ngauss_vdw-1)= &
!!$           (jnei-1)*vdw1%ngauss_vdw + [(ijk,ijk=1,vdw1%ngauss_vdw)]
         enddo
       endif
     enddo
   enddo
else
   write(*,*) ' We make near properly '
   call ghost_nodes(mesh0,x0)
   ! Compute the position of all the gauss points
   vdw1%x=0.d0
   do ielem=1,mesh0%numele
      neigh=mesh0%connect(ielem)%neigh_vert
      xneigh(:,1)=x0(3*neigh(:)-2)
      xneigh(:,2)=x0(3*neigh(:)-1)
      xneigh(:,3)=x0(3*neigh(:))
      do ig=1,vdw1%ngauss_vdw
         igg=(ielem-1)*vdw1%ngauss_vdw+ig
         do idim=1,3
            do inod=1,12
               vdw1%x(igg,idim)= vdw1%x(igg,idim) + vdw1%shapef(ig,inod)*xneigh(inod,idim)
            enddo
         enddo
      enddo
   enddo

   ! For each gp, Make a list with the gp that are too close
   ! to be considered "non-bonded". These are excluded later.
   vdw1%near=0
   do i=1,vdw1%ng_tot
      do j=i+1,vdw1%ng_tot
         dist=sqrt( (vdw1%x(i,1)-vdw1%x(j,1))**2 + (vdw1%x(i,2)-vdw1%x(j,2))**2 + &
              (vdw1%x(i,3)-vdw1%x(j,3))**2 )
         if (dist.le.r_c) then
            if(vdw1%near(i,0)+1.gt.nneigh-1) STOP ' Trouble with nneigh'
            vdw1%near(i,0)=vdw1%near(i,0)+1
            vdw1%near(i,vdw1%near(i,0))=j
         endif
      enddo
   enddo
endif

END SUBROUTINE make_near
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE Vvdw(a,A_p,sig,y0,Vs)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8), INTENT(IN) :: a,A_p,sig,y0
REAL(8), INTENT(OUT) :: Vs(2)

a1=sig/a
a6=(a1)**6
a12=a6*a6
y06=y0**6

Vs(1)=0.5d0*y06*a12 - a6
Vs(2)=a1*6.d0/sig*(-y06*a12 + a6)

Vs=Vs*A_p/(sig**6)

END SUBROUTINE Vvdw
