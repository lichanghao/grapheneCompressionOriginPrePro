SUBROUTINE mesh_gen_square(nrow,ncol,xlength,ylength,xx0,mesh0)
USE data_mesh
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(tri) :: void
TYPE(mesh) :: mesh0
REAL(8) ::  xx0(3*(mesh0%numnods))

! An empty triangle
void%vertices=[0,0,0]
void%num_neigh_elem=0
void%num_neigh_vert=0
void%neigh_elem=0
void%neigh_vert=0
void%code_bc = [0,0,0]

!initialize
mesh0%connect=void
mesh0%ntable=0
xx0=0.

! nodal positions
xstep=xlength/ncol
ystep=ylength/nrow

do i=1,nrow+1
  j=1
  inode = (i-1)*(ncol+1)+j
  ix0=inode*3-2
  iy=inode*3-1
  iz=inode*3
  xx0(ix0)=0.
  xx0(iy)=(i-1)*ystep
  xx0(iz)=0.
  do j=2,ncol
    inode = (i-1)*(ncol+1)+j
	ix=inode*3-2
	iy=inode*3-1
	iz=inode*3
	xx0(ix)=(j-1)*xstep
    xx0(iy)=(i-1)*ystep
    xx0(iz)=0.
  end do

  j=ncol+1
  inode = (i-1)*(ncol+1)+j
  ix=inode*3-2
  iy=inode*3-1
  iz=inode*3
  xx0(ix)=xlength
  xx0(iy)=(i-1)*ystep
  xx0(iz)=0.

end do

! Basic connectivity
do i=1,nrow
  do j=1,ncol
    ielem=((i-1)*ncol+j)*2 - 1
	in1=(i-1)*(ncol+1)+j
	in2=in1+1
	in3=in1+ncol+1
	in4=in3+1
    mesh0%connect(ielem)%vertices=[in1,in2,in3]
    mesh0%connect(ielem+1)%vertices=[in2,in4,in3]
  end do
end do

! Tag on boundary edges
nedg=0
do i=1,ncol
  ielem=i*2-1
  mesh0%connect(ielem)%code_bc(1)=1
  nedg=nedg+1
  ielem=2*ncol*(nrow-1)+i*2
  mesh0%connect(ielem)%code_bc(2)=1
  nedg=nedg+1
end do
do i=1,nrow
  ielem=(i-1)*ncol*2+1
  mesh0%connect(ielem)%code_bc(3)=1
  nedg=nedg+1
  ielem=ielem+ncol*2-1
  mesh0%connect(ielem)%code_bc(1)=1
  nedg=nedg+1
end do

mesh0%nedge=nedg


END SUBROUTINE mesh_gen_square

SUBROUTINE meshg_gen_square(nrow,ncol,xlength,ylength,xx0,mesh0)
USE data_mesh
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(tri) :: void
TYPE(mesh) :: mesh0
REAL(8) ::  xx0(3*(mesh0%numnods))


! An empty triangle
void%vertices=[0,0,0]
void%num_neigh_elem=0
void%num_neigh_vert=0
void%neigh_elem=0
void%neigh_vert=0
void%code_bc = [0,0,0]

!initialize
mesh0%connect=void
mesh0%ntable=0
xx0=0.

! nodal positions
xstep=xlength/ncol
ystep=ylength/nrow

do i=1,nrow+1
  j=1
  inode = (i-1)*(ncol+1)+j
  ix0=inode*3-2
  iy=inode*3-1
  iz=inode*3
  xx0(ix0)=0.
  xx0(iy)=(i-1)*ystep
  xx0(iz)=0.
  do j=2,ncol
    inode = (i-1)*(ncol+1)+j
	ix=inode*3-2
	iy=inode*3-1
	iz=inode*3
	xx0(ix)=(j-1)*xstep
    xx0(iy)=(i-1)*ystep
    xx0(iz)=0.
  end do

  j=ncol+1
  inode = (i-1)*(ncol+1)+j
  ix=inode*3-2
  iy=inode*3-1
  iz=inode*3
  xx0(ix)=xlength
  xx0(iy)=(i-1)*ystep
  xx0(iz)=0.

end do

! Basic connectivity
do i=1,nrow
  do j=1,ncol
    ielem=((i-1)*ncol+j)*2 - 1
	in1=(i-1)*(ncol+1)+j
	in2=in1+1
	in3=in1+ncol+1
	in4=in3+1
    mesh0%connect(ielem)%vertices=[in1,in2,in3]
    mesh0%connect(ielem+1)%vertices=[in2,in4,in3]
  end do
end do

! Tag on boundary edges
nedg=0
do i=1,ncol
   ielem=i*2-1
   mesh0%connect(ielem)%code_bc(1)=1
   nedg=nedg+1
   ielem=2*ncol*(nrow-1)+i*2
   mesh0%connect(ielem)%code_bc(2)=1
   nedg=nedg+1
end do
do i=1,nrow
   ielem=(i-1)*ncol*2+1
   mesh0%connect(ielem)%code_bc(3)=1
   nedg=nedg+1
   ielem=ielem+ncol*2-1
   mesh0%connect(ielem)%code_bc(1)=1
   nedg=nedg+1
end do

mesh0%nedge=nedg


END SUBROUTINE meshg_gen_square
