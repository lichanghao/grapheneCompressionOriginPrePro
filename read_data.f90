SUBROUTINE read_data(nsheets,nrow,ncol,xlength,ylength, &
           numel,numno,numed,mat1,stretch_ini, &
		   angle,angle2,nloadstep,crit,ngauss,nW_hat,nCodeLoad, &
                   imperfect,fact_imp,nborder,vdw1)
USE data_mat
USE data_vdw
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(material) :: mat1
TYPE(vdw_data):: vdw1
DIMENSION :: crit(2), Vaux(2)
REAL(8), ALLOCATABLE :: xn1(:), xn2(:)
REAL(8) :: xlength(nsheets),ylength(nsheets)
INTEGER(4) :: nrow(nsheets), ncol(nsheets), numel(nsheets), numno(nsheets), numed(nsheets)

ALLOCATE(xn1(nsheets),xn2(nsheets),STAT=istat)
if (istat/=0) STOP '**** Not enough memory ****'

write(*,*) 'how many ?', nsheets
read(33,*)
read(33,*)
read(33,*)
read(33,*) ngauss
read(33,*)
do i=1,nsheets 
  ! for the moment ncol should be the same for all tubes
  read(33,*) nrow(i), ncol(i) ! nrow is in the diameter, ncol in length
enddo
read(33,*)
read(33,*) xxlength,yylength
read(33,*)
read(33,*) stretch_ini
read(33,*)
read(33,*) mat1%A0
read(33,*)
do i=1,nsheets
  read(33,*) xn1(i), xn2(i)
enddo
read(33,*)
read(33,*) nchir
read(33,*) xorient
read(33,*)
read(33,*) mat1%nCode_Pot
read(33,*)
read(33,*) nCodeLoad
read(33,*)
read(33,*) angle
if (nCodeLoad.eq.13) read(33,*) angle2
read(33,*)
read(33,*) nloadstep
read(33,*)
read(33,*) nW_hat
read(33,*)
read(33,*) crit(1)
read(33,*)
read(33,*) crit(2)
read(33,*)
read(33,*) imperfect
read(33,*)
read(33,*) fact_imp
read(33,*)
read(33,*) nborder
read(33,*)
read(33,*) nvdw
if (nvdw.eq.1) then ! Van der Waals things
  vdw1%nvdw=nvdw
  read(33,*)
  read(33,*) vdw1%ngauss_vdw
  read(33,*)
  read(33,*) vdw1%r_cut
  read(33,*)
  read(33,*) vdw1%r_bond
  read(33,*)
  read(33,*) vdw1%a
  read(33,*)
  read(33,*) vdw1%sig
  read(33,*)
  read(33,*) vdw1%y0
  read(33,*)
  read(33,*) vdw1%meval
  vdw1%Vcut=0.d0
  call Vvdw(vdw1%r_cut,vdw1%a,vdw1%sig,vdw1%y0,Vaux)
  vdw1%Vcut=Vaux
endif

pi=dacos(-1.D0)

do i=1,nsheets
  xlength(i)=xxlength
  ylength(i)=yylength
enddo

if (nchir.eq.1) then
  ! Chirality introduced explicitely
  theta=xorient*pi/180.d0
else
  theta=datan(-dsqrt(3.d0)*xn2(1)/(2*xn1(1)+xn2(1)))
end if
mat1%E(1,:)=[dcos(theta), dsin(theta)]
mat1%E(2,:)=[dcos(theta+2.*pi/3.), dsin(theta+2.*pi/3.)]
mat1%E(3,:)=[dcos(theta-2.*pi/3.), dsin(theta-2.*pi/3.)]
mat1%s0=3.d0*dsqrt(3.d0)/2.d0*mat1%A0*mat1%A0

xlength=xlength*(ncol(1)+2*nborder)/ncol(1)

if(nborder == 0)then
   write(*,*)'we do not have boder elements'
endif

do i=1,nsheets
  nrow(i)=nrow(i)+2*nborder
  ncol(i)=ncol(i)+2*nborder
  numel(i)=nrow(i)*ncol(i)*2
  numno(i)=(nrow(i)+1)*(ncol(i)+1)
  numed(i)=2*(nrow(i)+ncol(i)+3)
enddo
END SUBROUTINE read_data



