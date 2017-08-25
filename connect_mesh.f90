SUBROUTINE connect_mesh(meshh)
USE data_mesh
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
logical ask, mask, flag
type(tri) :: triang, tri_master, tri_slave
TYPE(mesh) :: meshh
dimension  ivert(3), mask(maxneigh_vert), iaux(3)
integer, dimension(:,:), allocatable :: mtable
dimension iperm(2,3), neigh(3), neigh2(3)
dimension ilist_el(4), ilist_ve(4)
dimension :: iright(6), ileft(6), ihelp(3), ikk(3,2)
data iperm(1,1:3) /2, 3, 1/
data iperm(2,1:3) /3, 1, 2/


ALLOCATE(mtable(meshh%numnods,maxneigh_vert+1),STAT=istat)
if (istat/=0) STOP '**** Not enough memory ****'
write(*,*)maxneigh_vert,meshh%numnods

! Find all nodes connected to a particular node
! Find all the elements connected to particular node
meshh%ntable=0
mtable=0
do inode=1,meshh%numnods
  do ielem=1,meshh%numele
     ivert=meshh%connect(ielem)%vertices
     ask = any(ivert.eq.inode)
	 if (ask) then
       mtable(inode,1)=mtable(inode,1)+1
	   mtable(inode,1+mtable(inode,1))=ielem
	 end if
     if (ask) then
       do i=1,3
          iv=ivert(i)
          mask=meshh%ntable(inode,2:meshh%ntable(inode,1)+1).eq.iv
          ask=any(mask)
          if ((.not.ask).and.(iv/=inode)) then
            meshh%ntable(inode,1)=meshh%ntable(inode,1)+1
		    meshh%ntable(iv,1)=meshh%ntable(iv,1)+1
		    meshh%ntable(inode,1+meshh%ntable(inode,1))=iv
		    meshh%ntable(iv,1+meshh%ntable(iv,1))=inode
          end if
        end do
      end if
  end do
end do


write(*,*)'finish generate tables '

! Fill in gaps in connect

do ielem=1,meshh%numele
  triang=meshh%connect(ielem)                         ! data of ielem stored in triang
  call find_elem_adj(ielem,triang,mtable,neigh,meshh%numnods) 
  ivert=triang%vertices
  iaux=meshh%ntable(ivert(1:3),1)
  if (maxval(iaux).gt.6) then
    write(*,*) 'Error, program not prepared for valences > 6'
  end if

  if (neigh(1)/=0) then
!   Search on first edge to the right
    istep=0
    ilist_el=[3,2,1,12]
    ilist_ve=[3,1,2,5]
    ipivot=ivert(1)             ! we are pivoting around this vertex
    master=ielem                ! master element is for the moment the original one 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)  ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=i_master                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or.(mslave.eq.neigh(3))      ! if active edge of master is boundary or ...
      if (flag) exit                   ! ... reach edge 3, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(1,i_slave)) ! set neighbor vertex data
      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do

!   Search on the first edge to the left
    istep=0
    ilist_el=[4,5,6,0]
    ilist_ve=[6,10,11,0]
    ipivot=ivert(2)             ! we are pivoting around this vertex
    master=triang%neigh_elem(3)  ! master element is for the moment the original one 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)                ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=iperm(2,i_master)                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or.(mslave.eq.neigh(2))        ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... we reach edge 2 of ielem, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(2,i_slave)) ! set neighbor vertex data

      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do

  endif


  if (neigh(2)/=0) then
!   Search on second edge to the right
    istep=0
    ilist_el=[7,6,5,4]
    ilist_ve=[11,10,6,3]
    ipivot=ivert(2)             ! we are pivoting around this vertex
    master=ielem                ! master element is for the moment the original one 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)  ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=i_master                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or. &
  	     (triang%neigh_elem(ilist_el(istep+1))/=0)    ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... the data already filled, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      if ((triang%neigh_vert(ilist_ve(istep))/=0).and. &
	     (triang%neigh_vert(ilist_ve(istep))/=tri_slave%vertices(iperm(1,i_slave)))) &
		 write(*,*) ' Gross error in second right '
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(1,i_slave)) ! set neighbor vertex data
      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do
 
!   Search on the second edge left
    istep=0
    ilist_el=[8,9,10,0]
    ilist_ve=[12,9,5,0]
    ipivot=ivert(3)             ! we are pivoting around this vertex
    master=triang%neigh_elem(7)         ! master element 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)                ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=iperm(2,i_master)                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or.(mslave.eq.neigh(3))        ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... we reach edge 2 of ielem, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      if ((triang%neigh_vert(ilist_ve(istep))/=0).and. &
	     (triang%neigh_vert(ilist_ve(istep))/=tri_slave%vertices(iperm(2,i_slave)))) &
		 write(*,*) ' Gross error in second left '
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(2,i_slave)) ! set neighbor vertex data

      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do
  end if


  if (neigh(3)/=0) then
!   Search on third edge to the right
    istep=0
    ilist_el=[11,10,9,8]
    ilist_ve=[5,9,12,11]
    ipivot=ivert(3)             ! we are pivoting around this vertex
    master=ielem                ! master element is for the moment the original one 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)  ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=i_master                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or. &
  	     (triang%neigh_elem(ilist_el(istep+1))/=0)    ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... the data already filled, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      if ((triang%neigh_vert(ilist_ve(istep))/=0).and. &
	     (triang%neigh_vert(ilist_ve(istep))/=tri_slave%vertices(iperm(1,i_slave)))) &
		 write(*,*) ' Gross error in third right '
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(1,i_slave)) ! set neighbor vertex data
      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do
 
!   Search on the third edge left
    istep=0
    ilist_el=[12,1,2,0]
    ilist_ve=[2,1,3,0]
    ipivot=ivert(1)             ! we are pivoting around this vertex
    master=triang%neigh_elem(11)         ! master element 
    tri_master=meshh%connect(master)
    call glob_loc(ipivot,tri_master,i_master)                ! orient master element wrt ipivot
    do
      call find_elem_adj(master,tri_master,mtable,neigh2,meshh%numnods)  ! find neighbors of master element
      iact=iperm(2,i_master)                                ! identify active edge of master
      mslave=neigh2(iact)                               ! identify salve element of master
      flag=(mslave.eq.0).or.(mslave.eq.neigh(1))        ! if active edge of master is boundary or ...
      if (flag) exit                                    ! ... we reach edge 1 of ielem, finish with edge 1
      istep=istep+1                                     ! we advance one step
      tri_slave=meshh%connect(mslave)
      call glob_loc(ipivot,tri_slave,i_slave)           ! orient slave element wrt ipivot
      triang%neigh_elem(ilist_el(istep))=mslave         ! set neighbor element data
      if ((triang%neigh_vert(ilist_ve(istep))/=0).and. &
	     (triang%neigh_vert(ilist_ve(istep))/=tri_slave%vertices(iperm(2,i_slave)))) &
		 write(*,*) ' Gross error in third left '
      triang%neigh_vert(ilist_ve(istep))=tri_slave%vertices(iperm(2,i_slave)) ! set neighbor vertex data

      master=mslave
      tri_master=tri_slave
      i_master=i_slave
    end do
  end if

  triang%neigh_vert(4)=triang%vertices(1)
  triang%neigh_vert(7)=triang%vertices(2)
  triang%neigh_vert(8)=triang%vertices(3)
  triang%num_neigh_elem=count(triang%neigh_elem/=0)
  triang%num_neigh_vert=count(triang%neigh_vert/=0)

  meshh%connect(ielem)=triang

end do


END SUBROUTINE connect_mesh


SUBROUTINE connect_orig_mesh(mesh0,meshg,ncol,nrow)
USE data_mesh
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
logical ask, mask, flag
type(tri) :: triang,triangg
TYPE(mesh) :: mesh0,meshg
dimension  iperm(2,3)
dimension  ivert(3)
dimension  node((nrow+1)*(ncol+1)),nodeg((nrow+3)*(ncol+3))
dimension  node2g(2,(nrow+3)*(ncol+3))
dimension  nelem(2*nrow*ncol),nelemg(2*(nrow+2)*(ncol+2))
dimension  nelemfg(2*nrow*ncol),nelem2g(2*(nrow+2)*(ncol+2))
data iperm(1,1:3) /2, 3, 1/
data iperm(2,1:3) /3, 1, 2/

numno=(nrow+1)*(ncol+1)
numnog=(nrow+3)*(ncol+3)
numghost=2*nrow+2*ncol+6 ! numghost = numno-numnog-2
numel=2*nrow*ncol

nelem2g(1:2*(nrow+2)*(ncol+2))=0
do irow = 1, nrow
   do icol = 1,2*ncol
      ielem = (irow-1)*ncol*2+icol
      jelem = irow*2*(ncol+2)+2+icol
      nelemfg(ielem)=jelem ! transfer the element from real to aux
      nelem2g(jelem)=ielem
   enddo
enddo

node2g(1:2,1:(nrow+3)*(ncol+3))=0
!node2g(1,:)=0 : no entry
!           =1 : real node
!           =2 : ghost nodes

do irow = 1,nrow+1
   do icol = 1,ncol+1
      inode = (irow-1)*(ncol+1)+icol
      jnode = irow*(ncol+3)+1+icol
      node2g(1,jnode)=1
      node2g(2,jnode)=inode  ! tranfer the node index from aux to real
   enddo
enddo


irow =1
ighost=0
do icol =1,ncol+2
   ighost =ighost+1
   node2g(1,icol+1)=2  
   node2g(2,icol+1)=ighost
enddo

do irow=2,nrow+2
   inode=(irow-1)*(ncol+3)+1
   ighost=ighost+1
   node2g(1,inode)=2
   node2g(2,inode)=ighost

   inode=irow*(ncol+3)
   ighost=ighost+1
   node2g(1,inode)=2
   node2g(2,inode)=ighost

enddo

irow=nrow+3
do icol=1,ncol+2
   inode=(irow-1)*(ncol+3)+icol
   ighost=ighost+1
   node2g(1,inode)=2
   node2g(2,inode)=ighost

enddo


do ielem = 1, mesh0%numele
   triang=mesh0%connect(ielem)
   jelem=nelemfg(ielem)
   triangg=meshg%connect(jelem)

 !  write(*,*)ielem,jelem
 !  write(*,*)triangg%neigh_vert(1:12)
 !  write(*,*)

   if(triangg%num_neigh_vert/=12) then
      write(*,*) 'stop ! the element has less than 12 neighboring nodes'
      stop
   endif

   do jj=1,12
 !     write(*,*)node2g(1,triangg%neigh_vert(jj)),node2g(2,triangg%neigh_vert(jj))
      if(node2g(1,triangg%neigh_vert(jj))==0) write(*,*)'stop!, wrong interpretation'
      if(node2g(1,triangg%neigh_vert(jj))==1) then
         triang%neigh_vert(jj)= node2g(2,triangg%neigh_vert(jj))  !insert regular nodes
 !        write(*,*)triang%neigh_vert(jj)
      elseif(node2g(1,triangg%neigh_vert(jj))==2)then
         triang%neigh_vert(jj)= mesh0%numnods+node2g(2,triangg%neigh_vert(jj)) !insert ghost nodes
 !        write(*,*)triang%neigh_vert(jj)
      endif
   enddo
   
   if(triangg%num_neigh_elem/=12)then
      write(*,*) 'stop ! the element has less than 12 neighboring elements'
      stop
   endif
   do jj=1,12
      triang%neigh_elem(jj)= nelem2g(triangg%neigh_elem(jj))  !insert regular nodes
   enddo
   triang%num_neigh_elem=count(triang%neigh_elem/=0)
   triang%num_neigh_vert=count(triang%neigh_vert/=0)
 !  write(*,*) triang%num_neigh_elem,triang%num_neigh_vert
   mesh0%connect(ielem)=triang
enddo


!!$do ielem=1,mesh0%numele
!!$   triang=mesh0%connect(ielem)
!!$   
!!$   write(*,*)
!!$   write(*,*) ielem,nelemfg(ielem), triang%vertices(1:3)
!!$   write(*,*) triang%num_neigh_elem
!!$   write(*,*) triang%num_neigh_vert
!!$   do jj =1,12
!!$      write(*,*)triang%neigh_elem(jj),triang%neigh_vert(jj)
!!$   enddo
!!$
!!$enddo


! insert the ghost_node table
jghost=1
do ielem=1,2*ncol,2 ! the first row
   if (mesh0%connect(ielem)%code_bc(1)==1) then
      ivert=mesh0%connect(ielem)%vertices
      jghost=jghost+1
      ied =1 
      i1=ied
      i2=iperm(1,ied)
      i3=iperm(2,ied)
      mesh0%nghost_tab(jghost,:)=[ivert(i1),ivert(i2),ivert(i3)]
   else
      write(*,*)'stop!, this is not an exposed edgeA!'
   endif
enddo

jghost=ncol+3
do ielem=1,mesh0%numele,2*ncol ! the first col
   if (mesh0%connect(ielem)%code_bc(3)==1) then
      jghost=jghost+2
      ivert=mesh0%connect(ielem)%vertices
      ied=3
      i1=ied
      i2=iperm(1,ied)
      i3=iperm(2,ied)
      mesh0%nghost_tab(jghost,:)=[ivert(i1),ivert(i2),ivert(i3)]
   else
      write(*,*)'stop!, this is not an exposed edgeB!'
   endif
enddo


jghost=ncol+2
do ielem=2*ncol,mesh0%numele,2*ncol ! the last col
   if (mesh0%connect(ielem)%code_bc(1)==1) then
      jghost=jghost+2
      ivert=mesh0%connect(ielem)%vertices
      ied=1
      i1=ied
      i2=iperm(1,ied)
      i3=iperm(2,ied)
      mesh0%nghost_tab(jghost,:)=[ivert(i1),ivert(i2),ivert(i3)]
   else
      write(*,*)'stop!, this is not an exposed edgeC!'
   endif
enddo


jghost=ncol+2+2*(nrow+1)+1
do ielem=(nrow-1)*2*ncol+2,mesh0%numele,2 ! the last row
   if (mesh0%connect(ielem)%code_bc(2)==1) then
      jghost=jghost+1
      ivert=mesh0%connect(ielem)%vertices
      ied=2
      i1=ied
      i2=iperm(1,ied)
      i3=iperm(2,ied)
      mesh0%nghost_tab(jghost,:)=[ivert(i1),ivert(i2),ivert(i3)]
   else
      write(*,*)ielem, 'stop!, this is not an exposed edgeD!'
   endif
enddo

! special treatment for the corners: 6 corner ghost nodes
jghost=1 ! first corner, first ghost node
ivert=mesh0%connect(1)%vertices
mesh0%nghost_tab(jghost,:)=[ivert(1),ivert(1),ivert(3)]

jghost=ncol+3  ! first corner,second ghost node
ivert=mesh0%connect(1)%vertices
mesh0%nghost_tab(jghost,:)=[ivert(1),ivert(1),ivert(2)]

jghost=ncol+2 ! second corner
mesh0%nghost_tab(jghost,:)=[mesh0%numnods+ncol+1,ncol+1,ncol]

jghost=ncol+2+2*(nrow+1)+1 ! third corner
mesh0%nghost_tab(jghost,:)=[(ncol+1)*nrow+1,mesh0%numnods+jghost-2,(ncol+1)*(nrow-1)+1]

jghost=(ncol+2)+2*(nrow+1) ! last corner
mesh0%nghost_tab(jghost,:)=[mesh0%numnods,mesh0%numnods,mesh0%numnods-1]

jghost=2*(ncol+2)+2*(nrow+1) ! last corner
mesh0%nghost_tab(jghost,:)=[mesh0%numnods,mesh0%numnods,mesh0%numnods-ncol-1]

!!$do jghost=1,2*(ncol+2)+2*(nrow+1)
!!$   write(*,*)jghost,mesh0%nghost_tab(jghost,:)
!!$enddo


END SUBROUTINE CONNECT_ORIG_MESH


!*********************************************************************************
!*********************************************************************************
!*********************************************************************************

! Finds elements that share an edge with a particular element
SUBROUTINE find_elem_adj(ielem,triang,mtable,neigh,numnods)
USE data_tri
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
logical ask
TYPE(tri) :: triang
dimension mtable(numnods,15), iperm(2,3), neigh(3)
data iperm(1,1:3) /2, 3, 1/
data iperm(2,1:3) /3, 1, 2/

neigh=0

do ivertex=1,3
  iv1=triang%vertices(ivertex)
  iv2=triang%vertices(iperm(1,ivertex))
  if (triang%code_bc(ivertex)/=1) then
! not a boundary edge
    do k=1,mtable(iv1,1)
       jelem=mtable(iv1,k+1)
       ask=any(mtable(iv2,2:(mtable(iv2,1)+1)).eq.jelem)
       if (ask.and.(jelem/=ielem)) then
         neigh(ivertex)=jelem
    	 exit
       end if
    end do
    if (neigh(ivertex).eq.0) write(*,*) ' Error in neighbors of',ielem,ivertex
  end if
end do


END SUBROUTINE find_elem_adj

!*********************************************************************************
!*********************************************************************************
!*********************************************************************************

! given a global numbering, find the local numbering in a given element
! 
SUBROUTINE glob_loc(iglob,triang,iloc)
USE data_tri
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(tri) :: triang

iloc=0

if (iglob.eq.(triang%vertices(1))) iloc=1 
if (iglob.eq.(triang%vertices(2))) iloc=2 
if (iglob.eq.(triang%vertices(3))) iloc=3 
if (iloc.eq.0) write(*,*) 'Error in glob_loc'

END SUBROUTINE glob_loc

!*********************************************************************************
!*********************************************************************************
!*********************************************************************************

SUBROUTINE ghost_nodes(meshh,xx)
USE data_mesh
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(mesh) :: meshh
REAL(8) ::  xx(3*(meshh%numnods+meshh%nedge)), xaux(3)
!REAL(8) :: a(3), b(3), c(3), d(3)
!dimension  ivert(3), iperm(2,3)
!data iperm(1,1:3) /2, 3, 1/
!data iperm(2,1:3) /3, 1, 2/


do ijk=1,meshh%nedge
  i1=meshh%nghost_tab(ijk,1)
  i2=meshh%nghost_tab(ijk,2)
  i3=meshh%nghost_tab(ijk,3)
  xaux=xx(3*i1-2:3*i1)+xx(3*i2-2:3*i2)-xx(3*i3-2:3*i3)    ! parallelogram
  xx(3*(meshh%numnods+ijk)-2:3*(meshh%numnods+ijk))=xaux  ! position of ghost
end do


END SUBROUTINE ghost_nodes

!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
