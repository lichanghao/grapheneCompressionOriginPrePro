MODULE data_tri
  integer, parameter :: maxneigh_elem=12, maxneigh_vert=12
  TYPE tri
    INTEGER, DIMENSION(3) :: vertices
    INTEGER :: num_neigh_elem
    INTEGER :: num_neigh_vert
    INTEGER, DIMENSION(maxneigh_elem) :: neigh_elem
    INTEGER, DIMENSION(maxneigh_vert) :: neigh_vert
    INTEGER, DIMENSION(3) :: code_bc
  END TYPE tri
END MODULE

MODULE data_mesh
  USE data_tri
  TYPE mesh
    INTEGER :: numele
    INTEGER :: numnods
    INTEGER :: nedge
    TYPE(tri), DIMENSION(:), POINTER :: connect
    INTEGER, POINTER :: ntable(:,:) 
    INTEGER, POINTER :: nghost_tab(:,:)
    INTEGER :: nelem_ghost            ! number of ghost elements to remove boundary effect
    INTEGER, POINTER :: elem_ghost(:) ! list of ghost elements
    INTEGER :: nnode_ghost            ! number of ghost elements to remove boundary effect
    INTEGER, POINTER :: node_ghost(:) ! list of ghost nodes
  END TYPE
END MODULE

MODULE data_BC
  TYPE BC_data
    INTEGER :: nloadstep
    INTEGER :: nCodeLoad
    INTEGER :: nnodBC
    INTEGER :: ndofBC
    INTEGER :: ndofOP
    INTEGER, POINTER :: mdofBC(:)
    INTEGER, POINTER :: mnodBC(:,:)
    INTEGER, POINTER :: mdofOP(:)
    REAL(8) :: rotation(3,3), xc(3), value
  END TYPE
END MODULE

MODULE data_tensor22
   TYPE tensor22
     REAL(8) :: val(2,2)
   END TYPE
END MODULE

MODULE data_vector2
   TYPE vector2
     REAL(8) :: val(2)
   END TYPE
END MODULE

MODULE data_vector3
   TYPE vector3
     REAL(8) :: val(3)
   END TYPE
END MODULE

MODULE data_mat
  TYPE material
    INTEGER(4) :: nCode_Pot
    REAL(8) :: E(3,2), A0, s0, A1
    REAL(8) :: Vs(3)
    REAL(8) :: Va(3)
  END TYPE
END MODULE

MODULE data_vdw
  TYPE vdw_data
    INTEGER(4) :: nvdw           ! 1 if vdw are active
    INTEGER(4) :: neval          ! number of times evaluated.
    INTEGER(4) :: meval          ! How often update neighbors
    INTEGER(4) :: flag           ! indicates neighbor list has been updated
    INTEGER(4) :: ngauss_vdw     ! how many gp per element
    INTEGER(4) :: ng_tot         ! total number of gauss points for vdw
    INTEGER(4) :: nx(3)          ! Number of bins in each direction
    REAL(8) :: r_cut, r_bond, sig,  a, y0   ! data of the 12-6 potential, rvdw is the fraction of 
                                 ! sig at which is the cut-off radius, i.e. r_cut=r_vdw*sig
    REAL(8) :: Vcut(2)           ! Potential and derivative at cut-off, for smoothness
    INTEGER, POINTER :: near(:,:)       ! list of "bonded" neighbors of each gp, 
                                        ! the 0th entry of second index is number of neighbors 
    INTEGER, POINTER :: bin(:,:,:,:)    ! for a given bin, 0-# of gp that belong, followed by list
    INTEGER, POINTER :: gbin(:,:)       ! given gp, bin at which it belongs
    INTEGER, POINTER :: neigh(:,:)      ! given gp, 0th is # of valid neighs, followed by list
    REAL(8), POINTER :: shapef(:,:), weight(:)
    REAL(8), POINTER :: x(:,:), W(:)    ! positions and density 
  END TYPE
END MODULE

