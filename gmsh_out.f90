SUBROUTINE gmsh_out(xx,meshh,name,ntab)
USE data_mesh
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
CHARACTER(11), INTENT(IN) :: name
TYPE(mesh) :: meshh
REAL(8) ::  xx(3*(meshh%numnods))
!!$INTEGER(4), ALLOCATABLE :: ntab(:)
INTEGER(4) :: ntab(meshh%numnods)

open(unit=99,file=name,status='unknown')

!!$write(*,*) ' Allocate ntab'
!!$ALLOCATE(ntab(meshh%numnods),STAT=istat)
!!$if (istat/=0) STOP '**** Not enough memory ****'

! renumber
if (meshh%nelem_ghost.gt.0) then
  ii=0
  do inode=1,meshh%numnods
    if (all(inode.ne.meshh%node_ghost(:))) then
      ii=ii+1
      ntab(inode)=ii
    endif
  end do
else
  do inode=1,meshh%numnods
      ntab(inode)=inode
  end do
endif


! NODES
write(99,'(A)') '$NOD'
write(99,'(I10)') meshh%numnods-meshh%nnode_ghost
do inode=1,meshh%numnods
    if (all(inode.ne.meshh%node_ghost(:))) then
      inoden=ntab(inode)
      ix=inode*3-2
      iy=inode*3-1
      iz=inode*3
      write(99,'(I10,3E15.6)') inoden,xx(ix),xx(iy),xx(iz)
    endif
end do
write(99,'(A)') '$ENDNOD'

! ELEMENTS
write(99,'(A)') '$ELM'
write(99,'(I10)') meshh%numele-meshh%nelem_ghost

ii=0
do ielem=1,meshh%numele
  if (meshh%nelem_ghost.gt.0) then
    if (all(ielem.ne.meshh%elem_ghost(:))) then
      ii=ii+1
      write(99,'(8I10)') ii,2,100,6,3,ntab(meshh%connect(ielem)%vertices(1)), &
                                      ntab(meshh%connect(ielem)%vertices(2)), &
                                      ntab(meshh%connect(ielem)%vertices(3))
    endif
  else
    ii=ii+1
    write(99,'(8I10)') ii,2,100,6,3,meshh%connect(ielem)%vertices
  endif
end do
write(99,'(A)') '$ENDELM'
write(99,'(I10)') meshh%nedge


!!$write(*,*) ' Deallocate ntab'
!!$DEALLOCATE(ntab)
!!$write(*,*) ' Deallocated ntab'
close(99)


END SUBROUTINE gmsh_out
