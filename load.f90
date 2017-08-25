!*********************************************************************88
!*********************************************************************88
!*********************************************************************88
! define rotation
SUBROUTINE load_pre(x0,mesh0,BCs,xlength,ylength,& 
                    nrow,ncol,nborder,imesh,angle2)
USE data_mesh
USE data_BC
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
PARAMETER (pi0=3.141592653589793238d0)
TYPE(BC_data):: BCs
TYPE(mesh) :: mesh0
REAL(8) :: x0(3*(mesh0%numnods))
REAL(8) xb(3)

if ((BCs%nCodeLoad.eq.0)) then
  continue
else if ((BCs%nCodeLoad.eq.1).or.(BCs%nCodeLoad.eq.10)) then
  ! Bending
  BCs%value=BCs%value*pi0/180.
  BCs%rotation(1,1:3)=[ cos(BCs%value/BCs%nloadstep),sin(BCs%value/BCs%nloadstep),0.d0]
  BCs%rotation(2,1:3)=[-sin(BCs%value/BCs%nloadstep),cos(BCs%value/BCs%nloadstep),0.d0]
  BCs%rotation(3,1:3)=[0.d0,0.d0,1.d0]
else if (BCs%nCodeLoad.eq.13) then
  ! Bending
  BCs%value=BCs%value*pi0/180.
  BCs%rotation(1,1:3)=[ cos(BCs%value/BCs%nloadstep),sin(BCs%value/BCs%nloadstep),0.d0]
  BCs%rotation(2,1:3)=[-sin(BCs%value/BCs%nloadstep),cos(BCs%value/BCs%nloadstep),0.d0]
  BCs%rotation(3,1:3)=[angle2/2.d0/BCs%nloadstep,0.d0,0.d0]
else if (BCs%nCodeLoad.eq.2) then
  ! Twisting
  BCs%value=BCs%value*pi0/180.
  BCs%rotation(1,1:3)=[1.d0,0.d0,0.d0]
  BCs%rotation(2,1:3)=[0.d0, cos(BCs%value/BCs%nloadstep),sin(BCs%value/BCs%nloadstep)]
  BCs%rotation(3,1:3)=[0.d0,-sin(BCs%value/BCs%nloadstep),cos(BCs%value/BCs%nloadstep)]
else if ((BCs%nCodeLoad.eq.3).or.(BCs%nCodeLoad.eq.4).or.(BCs%nCodeLoad.eq.7)) then
  ! Axially deforming
  ! In this event, the increment vector is stored
  ! in the first row of rotation
  ! Remember that angle is the total percentage of compression
  BCs%rotation(1,1:3)=[0.d0,0.d0,BCs%value/BCs%nloadstep]
else if (BCs%nCodeLoad.eq.5) then
  ! One end is locked and the other has a force...
  ! The increment force is computed in the first row of rotation
  BCs%rotation(1,1:3)=[0.d0,BCs%value/BCs%nloadstep*ylength/nrow,0.d0]
else if (BCs%nCodeLoad.eq.6) then
  ! Squeeze the tube
  BCs%rotation(1,1:3)=[0.d0,-BCs%value*ylength/pi0/200.d0/BCs%nloadstep,0.d0]
else if (BCs%nCodeLoad.eq.666) then
  ! Squeeze the multi-tube
  continue
else if (BCs%nCodeLoad.eq.8) then
  BCs%rotation(1,1:3)=[0.d0,BCs%value/BCs%nloadstep*ylength/nrow,0.d0]
else
  STOP ' Loading option not implemented2'
end if

! Coordinates of the center of the nanotube
BCs%xc=0.
do inod=1,mesh0%numnods
  xb=x0(3*inod-2:3*inod)
  BCs%xc=BCs%xc+xb
end do
BCs%xc=BCs%xc/mesh0%numnods

if ((BCs%nCodeLoad.eq.0)) then
  BCs%value=0.
  BCs%nloadstep=1
  BCs%mnodBC(1:(BCs%nnodBC),1)= [(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]
  BCs%mnodBC(:,2)= 0

  if (imesh.eq.1) then
    BCs%mdofBC(1:(BCs%ndofBC-3))= 3*[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]-2
    BCs%mdofBC(BCs%ndofBC-2:BCs%ndofBC)=[3, 3*(1+(ncol+1)*nrow/4) -1 ,3*(1+(ncol+1)*nrow/2)]
  else
    BCs%mdofBC(1:(BCs%ndofBC))= 3*[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]-2
  endif
else if ((BCs%nCodeLoad.eq.7)) then
  BCs%mnodBC(1:4,1)= [(1+(ncol+1)*nrow/4)+nborder,(1+(ncol+1)*nrow/4*3)+nborder, &
                     (1+(ncol+1)*nrow/4)+ncol-nborder,(1+(ncol+1)*nrow/4*3)+ncol-nborder]
  BCs%mnodBC(1:2,2)= 1
  BCs%mnodBC(3:4,2)= 2

  BCs%mdofBC(1:3)  = [(ijk,ijk=3*BCs%mnodBC(1,1)-2,3*BCs%mnodBC(1,1))]
  BCs%mdofBC(4:6)  = [(ijk,ijk=3*BCs%mnodBC(2,1)-2,3*BCs%mnodBC(2,1))]
  BCs%mdofBC(7:9)  = [(ijk,ijk=3*BCs%mnodBC(3,1)-2,3*BCs%mnodBC(3,1))]
  BCs%mdofBC(10:12)= [(ijk,ijk=3*BCs%mnodBC(4,1)-2,3*BCs%mnodBC(4,1))]

else if ((BCs%nCodeLoad.eq.1).or.(BCs%nCodeLoad.eq.2).or. &
         (BCs%nCodeLoad.eq.3).or.(BCs%nCodeLoad.eq.13)) then

  ! all the edge nodes are checked out

!!$  BCs%mnodBC(1:BCs%nnodBC/2,1)=[(ijk,ijk=1,nrow*(ncol+1)+1,ncol+1)]+nborder
!!$  BCs%mnodBC(1:BCs%nnodBC/2,2)=1
!!$  BCs%mnodBC(BCs%nnodBC/2+1:BCs%nnodBC,1)=[(ijk,ijk=1,nrow*(ncol+1)+1,ncol+1)]+ &
!!$                                          ncol-nborder
!!$  BCs%mnodBC(BCs%nnodBC/2+1:BCs%nnodBC,2)=2
!!$
!!$  BCs%mdofBC(1:BCs%ndofBC/2:3)=3*([(ijk,ijk=1,nrow*(ncol+1)+1,ncol+1)]+nborder)-2
!!$  BCs%mdofBC(2:BCs%ndofBC/2:3)=3*([(ijk,ijk=1,nrow*(ncol+1)+1,ncol+1)]+nborder)-1
!!$  BCs%mdofBC(3:BCs%ndofBC/2:3)=3*([(ijk,ijk=1,nrow*(ncol+1)+1,ncol+1)]+nborder)
!!$
!!$  BCs%mdofBC(BCs%ndofBC/2+1:BCs%ndofBC:3)=3*([(ijk,ijk=1,nrow*(ncol+1)+1,ncol+1)]+ncol-nborder)-2
!!$  BCs%mdofBC(BCs%ndofBC/2+2:BCs%ndofBC:3)=3*([(ijk,ijk=1,nrow*(ncol+1)+1,ncol+1)]+ncol-nborder)-1
!!$  BCs%mdofBC(BCs%ndofBC/2+3:BCs%ndofBC:3)=3*([(ijk,ijk=1,nrow*(ncol+1)+1,ncol+1)]+ncol-nborder)

!!$  write(*,*)BCs%nnodBC
!!$  write(*,*)BCs%mnodBC(1:BCs%nnodBC,1)
!!$  write(*,*)
!!$  write(*,*)BCs%mnodBC(1:BCs%nnodBC,2)
!!$  write(*,*)
!!$  write(*,*)BCs%ndofBC
!!$  write(*,*)BCs%mdofBC(1:BCs%ndofBC)



else if ((BCs%nCodeLoad.eq.10)) then

  BCs%mdofBC(1:3*nrow:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)-2
  BCs%mdofBC(2:3*nrow:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)-1
  BCs%mdofBC(3:3*nrow:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)

  BCs%mdofBC(3*nrow+1:6*nrow:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol-nborder)-2
  BCs%mdofBC(3*nrow+2:6*nrow:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol-nborder)-1
  BCs%mdofBC(3*nrow+3:6*nrow:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol-nborder)

  BCs%mdofBC(6*nrow+1:6*nrow+(ncol-1))=3*[(ijk,ijk=2,ncol)]
  BCs%mdofBC(6*nrow+(ncol-1)+1:6*nrow+2*(ncol-1))=3* &
        [(ijk,ijk=(ncol+1)*(nrow/2)+2,(ncol+1)*(nrow/2)+ncol)]

  BCs%mnodBC(1:nrow,1)=       [(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder
  BCs%mnodBC(1:nrow,2)=1
  BCs%mnodBC(nrow+1:2*nrow,1)=[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]-nborder+ncol
  BCs%mnodBC(nrow+1:2*nrow,2)=2
  BCs%mnodBC(2*nrow+1:2*nrow+(ncol-1),1)=[(ijk,ijk=2,ncol)]
  BCs%mnodBC(2*nrow+(ncol-1)+1:2*nrow+2*(ncol-1),1)= &
                     [(ijk,ijk=(ncol+1)*(nrow/2)+2,(ncol+1)*(nrow/2)+ncol)]
  BCs%mnodBC(2*nrow+1:2*nrow+2*(ncol-1),2)=0

else if ((BCs%nCodeLoad.eq.4)) then

  BCs%mnodBC(1:BCs%nnodBC/2,1)=[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder
  BCs%mnodBC(1:BCs%nnodBC/2,2)=1
  BCs%mnodBC(BCs%nnodBC/2+1:BCs%nnodBC,1)=[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol-nborder
  BCs%mnodBC(BCs%nnodBC/2+1:BCs%nnodBC,2)=2

  BCs%mdofBC(1:(BCs%ndofBC-3)/2)=      3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)-2
  BCs%mdofBC((BCs%ndofBC-3)/2+1:BCs%ndofBC-3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol -nborder)-2  
  BCs%mdofBC(BCs%ndofBC-2:BCs%ndofBC)=[3, 3*(1+(ncol+1)*nrow/4) -1 ,3*(1+(ncol+1)*nrow/2)]
else if ((BCs%nCodeLoad.eq.5)) then

  BCs%mnodBC(1:BCs%nnodBC/2,1)=[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder
  BCs%mnodBC(1:BCs%nnodBC/2,2)=1
  BCs%mnodBC(BCs%nnodBC/2+1:BCs%nnodBC,1)=[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol-nborder
  BCs%mnodBC(BCs%nnodBC/2+1:BCs%nnodBC,2)=2

  BCs%mdofBC(1:BCs%ndofBC:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)-2
  BCs%mdofBC(2:BCs%ndofBC:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)-1
  BCs%mdofBC(3:BCs%ndofBC:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)
else if ((BCs%nCodeLoad.eq.6).and.(imesh.eq.1)) then
  ! Constrained nodes
  BCs%mnodBC(1:(BCs%nnodBC-nrow)/2+1,1)=[(ijk,ijk=1,ncol+1)]
  BCs%mnodBC(  (BCs%nnodBC-nrow)/2+2:BCs%nnodBC-nrow+2,1)= &
        [(ijk,ijk=(ncol+1)*(nrow/2)+1,(ncol+1)*(nrow/2)+1+ncol)]
  BCs%mnodBC(BCs%nnodBC-nrow+3:BCs%nnodBC-nrow+1+nrow/2,1)=[(ijk,ijk=1+(ncol+1),(nrow/2-1)*(ncol+1)+1,ncol+1)]
  BCs%mnodBC(BCs%nnodBC-nrow+2+nrow/2:BCs%nnodBC,1)=[(ijk,ijk=(nrow/2+1)*(ncol+1)+1,(nrow-1)*(ncol+1)+1,ncol+1)]
  ! Codes
  BCs%mnodBC(:,2)=0
  BCs%mnodBC(1:(BCs%nnodBC-nrow)/2+1,2)=1
  BCs%mnodBC(  (BCs%nnodBC-nrow)/2+2:BCs%nnodBC-nrow+2,2)=2


  ! Constrained d.o.f.'s
  BCs%mdofBC(1:(BCs%ndofBC-nrow)/4)=3*[(ijk,ijk=1,ncol+1)]-1
  BCs%mdofBC(  (BCs%ndofBC-nrow)/4+1:(BCs%ndofBC-nrow)/2)=3* &
                  [(ijk,ijk=(ncol+1)*(nrow/2)+1,(ncol+1)*(nrow/2)+1+ncol)]-1
  BCs%mdofBC((BCs%ndofBC-nrow)/2+1:(BCs%ndofBC-nrow)/4*3)=3*[(ijk,ijk=1,ncol+1)]
  BCs%mdofBC((BCs%ndofBC-nrow)/4*3+1:(BCs%ndofBC-nrow))=3* &
                  [(ijk,ijk=(ncol+1)*(nrow/2)+1,(ncol+1)*(nrow/2)+1+ncol)]
  BCs%mdofBC(BCs%ndofBC-nrow+1:BCs%ndofBC)=3*[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]-2

else if ((BCs%nCodeLoad.eq.6).and.(imesh.ne.1)) then
  ! Constrained nodes
  BCs%mnodBC(1:(BCs%nnodBC-nrow)/2+1,1)=[(ijk,ijk=1,ncol+1)]
  BCs%mnodBC(  (BCs%nnodBC-nrow)/2+2:BCs%nnodBC-nrow+2,1)= &
        [(ijk,ijk=(ncol+1)*(nrow/2)+1,(ncol+1)*(nrow/2)+1+ncol)]
  BCs%mnodBC(BCs%nnodBC-nrow+3:BCs%nnodBC-nrow+1+nrow/2,1)=[(ijk,ijk=1+(ncol+1),(nrow/2-1)*(ncol+1)+1,ncol+1)]
  BCs%mnodBC(BCs%nnodBC-nrow+2+nrow/2:BCs%nnodBC,1)=[(ijk,ijk=(nrow/2+1)*(ncol+1)+1,(nrow-1)*(ncol+1)+1,ncol+1)]
  ! Codes
  BCs%mnodBC(:,2)=0

  ! Constrained d.o.f.'s
  BCs%mdofBC(1:(BCs%ndofBC-nrow)/2)=3*[(ijk,ijk=1,ncol+1)]
  BCs%mdofBC(  (BCs%ndofBC-nrow)/2+1:(BCs%ndofBC-nrow))=3* &
                  [(ijk,ijk=(ncol+1)*(nrow/2)+1,(ncol+1)*(nrow/2)+1+ncol)]
  BCs%mdofBC(BCs%ndofBC-nrow+1:BCs%ndofBC)=3*[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]-2
else if ((BCs%nCodeLoad.eq.666)) then

  BCs%mnodBC(1:BCs%nnodBC/2,1)=[(ijk,ijk=1,ncol+1)]
  BCs%mnodBC(BCs%nnodBC/2+1:BCs%nnodBC,1)=[(ijk,ijk=(ncol+1)*(nrow/2)+1,(ncol+1)*(nrow/2)+1+ncol)]
  BCs%mnodBC(1:BCs%nnodBC/2,2)=0
  BCs%mnodBC(BCs%nnodBC/2+1:BCs%nnodBC,2)=0

  BCs%mdofBC(1:(BCs%ndofBC-2)/2)=3*[(ijk,ijk=1,ncol+1)]
  BCs%mdofBC((BCs%ndofBC-2)/2+1:(BCs%ndofBC-2))=3* &
                  [(ijk,ijk=(ncol+1)*(nrow/2)+1,(ncol+1)*(nrow/2)+1+ncol)]

  BCs%mdofBC(BCs%ndofBC-1:BCs%ndofBC)=[1,3*((ncol+1)*(nrow/2)+1)-2]
else if ((BCs%nCodeLoad.eq.8)) then

  ! Constrained nodes
  BCs%mnodBC(1:nrow,1)=[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder
  BCs%mnodBC(1:nrow,2)=1
  BCs%mnodBC(nrow+1:2*nrow,1)=[(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol-nborder
  BCs%mnodBC(nrow+1:2*nrow,2)=2

!  BCs%mdofBC(1:nrow*3:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)-2
!  BCs%mdofBC(2:nrow*3:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)-1
!  BCs%mdofBC(3:nrow*3:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+nborder)
!  BCs%mdofBC(nrow*3+1:nrow*6:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol-nborder)-2
!  BCs%mdofBC(nrow*3+2:nrow*6:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol-nborder)-1
!  BCs%mdofBC(nrow*3+3:nrow*6:3)=3*([(ijk,ijk=1,(nrow-1)*(ncol+1)+1,ncol+1)]+ncol-nborder)


  ! Constrained nodes
  BCs%mnodBC(nrow*2+1:nrow*2+(ncol-1),1)=[(ijk,ijk=2,ncol)]
  BCs%mnodBC(nrow*2+(ncol-1)+1:nrow*2+2*(ncol-1),1)= &
        [(ijk,ijk=(ncol+1)*(nrow/2)+2,(ncol+1)*(nrow/2)+ncol)]
  BCs%mnodBC(nrow*2+2*(ncol-1)+1,1)=(int(ncol/2)+1)+(ncol+1)
  BCs%mnodBC(nrow*2+2*(ncol-1)+2,1)=(int(ncol/2)+1)-(ncol+1) + (ncol+1)*nrow
  ! Codes
  BCs%mnodBC(nrow*2+1:nrow*2+2*(ncol-1)+2,2)=0
  ! Constrained d.o.f.'s
  BCs%mdofBC(1:(ncol-1))=3*[(ijk,ijk=2,ncol)]
  BCs%mdofBC((ncol-1)+1:2*(ncol-1))=3* &
        [(ijk,ijk=(ncol+1)*(nrow/2)+2,(ncol+1)*(nrow/2)+ncol)]
  BCs%mdofBC(2*(ncol-1)+1)=3*(int(ncol/2)+1)-2
  BCs%mdofBC(2*(ncol-1)+2)=3*(int(ncol/2)+1)-1
  BCs%mdofBC(2*(ncol-1)+3)=3*((int(ncol/2)+1)+(ncol+1))-1
  BCs%mdofBC(2*(ncol-1)+4)=3*((int(ncol/2)+1)-(ncol+1) + (ncol+1)*nrow)-1

else
  STOP ' Loading option not implemented3'
end if


!!$j=0
!!$do i=1,3*mesh0%numnods
!!$  if (all(i.ne.BCs%mdofBC)) then
!!$    j=j+1
!!$    BCs%mdofOP(j)=i
!!$  end if
!!$end do

END SUBROUTINE load_pre


