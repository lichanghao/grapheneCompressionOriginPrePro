!***************************************************************************
!***************************************************************************
!***************************************************************************
! This Routine computes F0, the inverse of the tangent between T_bar
! and T0 in cartesian in both domains. This coincides with the componets
! of the deformation gradient in cartesian/convected
! It also computes J0, the jacobian to be used to integrate in T0 with 
! respect to the parametric space T_bar (dT0=J0*dT_bar)
SUBROUTINE Def_Grad_Cart_Conv(temp,mesh0,x0,J0,F0)
USE data_mesh
USE data_tensor22
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(mesh) :: mesh0
TYPE(tensor22) :: F0(mesh0%numele), temp(mesh0%numele)
REAL(8) :: x0(3*(mesh0%numnods+mesh0%nedge)), J0(mesh0%numele)

forall (i=1:mesh0%numele)
  temp(i)%val(1,1)=x0(3*mesh0%connect(i)%vertices(2)-2)-x0(3*mesh0%connect(i)%vertices(1)-2)
  temp(i)%val(1,2)=x0(3*mesh0%connect(i)%vertices(3)-2)-x0(3*mesh0%connect(i)%vertices(1)-2)
  temp(i)%val(2,1)=x0(3*mesh0%connect(i)%vertices(2)-1)-x0(3*mesh0%connect(i)%vertices(1)-1)
  temp(i)%val(2,2)=x0(3*mesh0%connect(i)%vertices(3)-1)-x0(3*mesh0%connect(i)%vertices(1)-1)
end forall
J0(:)=(temp(:)%val(1,1)*temp(:)%val(2,2)-temp(:)%val(1,2)*temp(:)%val(2,1))
F0(:)%val(1,1)=1./J0(:)*temp(:)%val(2,2)
F0(:)%val(1,2)=-1./J0(:)*temp(:)%val(1,2)
F0(:)%val(2,1)=-1./J0(:)*temp(:)%val(2,1)
F0(:)%val(2,2)=1./J0(:)*temp(:)%val(1,1)
END SUBROUTINE Def_Grad_Cart_Conv

