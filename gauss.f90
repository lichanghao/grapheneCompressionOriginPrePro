SUBROUTINE gauss(ngauss,shapef,weight)
  implicit REAL(8) (a-h,o-z)
  implicit INTEGER*4 (i-n)
  REAL(8) shapef(ngauss,12,6), weight(ngauss)

  if (ngauss==1) then
     call BSpline(shapef(1,:,1),1./3.,1./3.)
     call DBSpline(shapef(1,:,2:3),1./3.,1./3.)
     call DDBSpline(shapef(1,:,4:6),1./3.,1./3.)
     weight(1)=1.d0
  else if (ngauss==2) then

     pos1=1.d0/6.d0
     pos2=2.d0/3.d0

     call BSpline(shapef(1,:,1),pos1,pos2)
     call DBSpline(shapef(1,:,2:3),pos1,pos2)
     call DDBSpline(shapef(1,:,4:6),pos1,pos2)

     call BSpline(shapef(2,:,1),pos2,pos1)
     call DBSpline(shapef(2,:,2:3),pos2,pos1)
     call DDBSpline(shapef(2,:,4:6),pos2,pos1)

     weight=[1.d0/2.d0,1.d0/2.d0]

  else if (ngauss==3) then

     pos1=1.d0/6.d0
     pos2=2.d0/3.d0


     call BSpline(shapef(1,:,1),pos1,pos1)
     call DBSpline(shapef(1,:,2:3),pos1,pos1)
     call DDBSpline(shapef(1,:,4:6),pos1,pos1)

     call BSpline(shapef(2,:,1),pos2,pos1)
     call DBSpline(shapef(2,:,2:3),pos2,pos1)
     call DDBSpline(shapef(2,:,4:6),pos2,pos1)

     call BSpline(shapef(3,:,1),pos1,pos2)
     call DBSpline(shapef(3,:,2:3),pos1,pos2)
     call DDBSpline(shapef(3,:,4:6),pos1,pos2)

     !  call BSpline(shapef(1,:,1),1./2.,0.)
     !  call DBSpline(shapef(1,:,2:3),1./2.,0.)
     !  call DDBSpline(shapef(1,:,4:6),1./2.,0.)

     !  call BSpline(shapef(2,:,1),0.,1./2.)
     !  call DBSpline(shapef(2,:,2:3),0.,1./2.)
     !  call DDBSpline(shapef(2,:,4:6),0.,1./2.)

     !  call BSpline(shapef(3,:,1),1./2.,1./2.)
     !  call DBSpline(shapef(3,:,2:3),1./2.,1./2.)
     !  call DDBSpline(shapef(3,:,4:6),1./2.,1./2.)

     weight=[1.d0/3.d0,1.d0/3.d0,1.d0/3.d0]
  else
     STOP ' Number of Gauss points not implemented'
  end if

END SUBROUTINE gauss

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE gauss_vdw(ngauss,shapef,weight)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8) shapef(ngauss,12), weight(ngauss)
REAL(8), ALLOCATABLE :: position_in(:), position_out(:), weight_in(:), weight_out(:)

if (ngauss==1) then
  call BSpline(shapef(1,:),1./3.,1./3.)
  weight(1)=1.d0
else if (ngauss==2) then

  pos1=1.d0/6.d0
  pos2=2.d0/3.d0

  call BSpline(shapef(1,:),pos1,pos2)
  call BSpline(shapef(2,:),pos2,pos1)

  weight=[1.d0/2.d0,1.d0/2.d0]

else if (ngauss==3) then

  pos1=1.d0/6.d0
  pos2=2.d0/3.d0


  call BSpline(shapef(1,:),pos1,pos1)
  call BSpline(shapef(2,:),pos2,pos1)
  call BSpline(shapef(3,:),pos1,pos2)

  weight=[1.d0/3.d0,1.d0/3.d0,1.d0/3.d0]

else if (ngauss==6) then

!  pos1=0.659027622374092d0
!  pos2=0.231933368553031d0
!  pos3=0.109039009072877d0
  pos1=0.816847572980459d0
  pos2=0.091576213509771d0
  pos11=0.108103018168070d0
  pos22=0.445948490915965d0

  call BSpline(shapef(1,:),pos1,pos2)
  call BSpline(shapef(2,:),pos2,pos1)
  call BSpline(shapef(3,:),pos2,pos2)
  call BSpline(shapef(4,:),pos11,pos22)
  call BSpline(shapef(5,:),pos22,pos11)
  call BSpline(shapef(6,:),pos22,pos22)

  w1=0.109951743655322d0
  w2=0.223381589678011d0
  weight=[w1,w1,w1,w2,w2,w2]

else if (ngauss==9) then

  pos1=0.124949503233232d0
  pos2=0.437525248383384d0
  pos11=0.797112651860071d0
  pos22=0.165409927389841d0
  pos33=0.037477420750088d0

  call BSpline(shapef(1,:),pos1,pos2)
  call BSpline(shapef(2,:),pos2,pos1)
  call BSpline(shapef(3,:),pos2,pos2)
  call BSpline(shapef(4,:),pos11,pos22)
  call BSpline(shapef(5,:),pos11,pos33)
  call BSpline(shapef(6,:),pos22,pos11)
  call BSpline(shapef(7,:),pos22,pos33)
  call BSpline(shapef(8,:),pos33,pos11)
  call BSpline(shapef(9,:),pos33,pos22)

  weight=[ 0.205950504760887d0, 0.205950504760887d0, 0.205950504760887d0, &
           0.063691414286223d0, 0.063691414286223d0, 0.063691414286223d0, &
           0.063691414286223d0, 0.063691414286223d0, 0.063691414286223d0]

else if (ngauss==12) then

  pos1=0.873821971016996d0
  pos2=0.063089014491502d0

  call BSpline(shapef(1,:),pos1,pos2)
  call BSpline(shapef(2,:),pos2,pos1)
  call BSpline(shapef(3,:),pos2,pos2)

  pos11=0.501426509658179d0
  pos22=0.249286745170910d0

  call BSpline(shapef(4,:),pos11,pos22)
  call BSpline(shapef(5,:),pos22,pos11)
  call BSpline(shapef(6,:),pos22,pos22)

  pos111=0.636502499121399d0
  pos222=0.310352451033785d0
  pos333=0.053145049844816d0

  call BSpline(shapef(7,:),pos111,pos222)
  call BSpline(shapef(8,:),pos111,pos333)
  call BSpline(shapef(9,:),pos222,pos111)
  call BSpline(shapef(10,:),pos222,pos333)
  call BSpline(shapef(11,:),pos333,pos111)
  call BSpline(shapef(12,:),pos333,pos222)

  weight(1:3)=0.050844906370207d0
  weight(4:6)=0.116786275726379d0
  weight(7:12)=0.082851075618374d0

else if (ngauss==24) then  ! this is one subdivision from 6 gauss points
  pos1=0.816847572980459d0/2.d0
  pos2=0.091576213509771d0/2.d0
  pos11=0.108103018168070d0/2.d0
  pos22=0.445948490915965d0/2.d0

  call BSpline(shapef(1,:),pos1,pos2)
  call BSpline(shapef(2,:),pos2,pos1)
  call BSpline(shapef(3,:),pos2,pos2)
  call BSpline(shapef(4,:),pos11,pos22)
  call BSpline(shapef(5,:),pos22,pos11)
  call BSpline(shapef(6,:),pos22,pos22)

  call BSpline(shapef( 7,:),pos1+.5d0,pos2)
  call BSpline(shapef( 8,:),pos2+.5d0,pos1)
  call BSpline(shapef( 9,:),pos2+.5d0,pos2)
  call BSpline(shapef(10,:),pos11+.5d0,pos22)
  call BSpline(shapef(11,:),pos22+.5d0,pos11)
  call BSpline(shapef(12,:),pos22+.5d0,pos22)

  call BSpline(shapef(13,:),pos1,pos2+.5d0)
  call BSpline(shapef(14,:),pos2,pos1+.5d0)
  call BSpline(shapef(15,:),pos2,pos2+.5d0)
  call BSpline(shapef(16,:),pos11,pos22+.5d0)
  call BSpline(shapef(17,:),pos22,pos11+.5d0)
  call BSpline(shapef(18,:),pos22,pos22+.5d0)

  call BSpline(shapef(19,:),.5d0-pos1,.5d0-pos2)
  call BSpline(shapef(20,:),.5d0-pos2,.5d0-pos1)
  call BSpline(shapef(21,:),.5d0-pos2,.5d0-pos2)
  call BSpline(shapef(22,:),.5d0-pos11,.5d0-pos22)
  call BSpline(shapef(23,:),.5d0-pos22,.5d0-pos11)
  call BSpline(shapef(24,:),.5d0-pos22,.5d0-pos22)

  w1=0.109951743655322d0/4.d0
  w2=0.223381589678011d0/4.d0
  weight(1:6)=[w1,w1,w1,w2,w2,w2]
  weight(7:12)=[w1,w1,w1,w2,w2,w2]
  weight(13:18)=[w1,w1,w1,w2,w2,w2]
  weight(19:24)=[w1,w1,w1,w2,w2,w2]

else if (ngauss==16) then ! this is two subdivisions from 1 gauss points
  pos1=1.d0/3.d0
  ngauss_in=1
  ALLOCATE(position_in(ngauss_in*2),weight_in(ngauss_in), &
           position_out(ngauss_in*8),weight_out(ngauss_in*4),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'
  position_in=[pos1,pos1]
  weight_in=[1.d0]
  call subdiv_gauss(ngauss_in,position_in,weight_in,ngauss_out,position_out,weight_out)
  DEALLOCATE(position_in,weight_in)
  ALLOCATE(position_in(ngauss_out*8),weight_in(ngauss_out*4),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'
  call subdiv_gauss(ngauss_out,position_out,weight_out,ngauss_in,position_in,weight_in)
  if(ngauss_in.ne.16) STOP ' no chuto el 16'
  weight(:)=weight_in(:)

  do ikk=1,ngauss
    call BSpline(shapef(ikk,:),position_in(2*ikk-1),position_in(2*ikk))  
  enddo
!!$  open(unit=65,file='g_vdw.dat',status='unknown')
!!$     do jhh=1,ngauss
!!$       write(65,*) position_in(2*jhh-1),position_in(2*jhh)
!!$     enddo
!!$     write(65,*)
!!$     write(65,*) weight(:)
!!$  close(65)
else if (ngauss==32) then ! this is two subdivisions from 2 gauss points
  pos1=1.d0/6.d0
  pos2=2.d0/3.d0
  ngauss_in=2
  ALLOCATE(position_in(ngauss_in*2),weight_in(ngauss_in), &
           position_out(ngauss_in*8),weight_out(ngauss_in*4),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'
  position_in=[pos1,pos2,pos2,pos1]
  weight_in=[1.d0/2.d0,1.d0/2.d0]
  call subdiv_gauss(ngauss_in,position_in,weight_in,ngauss_out,position_out,weight_out)
  DEALLOCATE(position_in,weight_in)
  ALLOCATE(position_in(ngauss_out*8),weight_in(ngauss_out*4),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'
  call subdiv_gauss(ngauss_out,position_out,weight_out,ngauss_in,position_in,weight_in)
  if(ngauss_in.ne.32) STOP ' no chuto el 32'
  weight(:)=weight_in(:)

  do ikk=1,ngauss
    call BSpline(shapef(ikk,:),position_in(2*ikk-1),position_in(2*ikk))  
  enddo
else if (ngauss==48) then ! this is two subdivisions from 3 gauss points
  pos1=1.d0/6.d0
  pos2=2.d0/3.d0
  ngauss_in=3
  ALLOCATE(position_in(ngauss_in*2),weight_in(ngauss_in), &
           position_out(ngauss_in*8),weight_out(ngauss_in*4),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'
  position_in=[pos1,pos1,pos2,pos1,pos1,pos2]
  weight_in=[1.d0/3.d0,1.d0/3.d0,1.d0/3.d0]
  call subdiv_gauss(ngauss_in,position_in,weight_in,ngauss_out,position_out,weight_out)
  DEALLOCATE(position_in,weight_in)
  ALLOCATE(position_in(ngauss_out*8),weight_in(ngauss_out*4),STAT=istat)
  if (istat/=0) STOP '**** Not enough memory ****'
  call subdiv_gauss(ngauss_out,position_out,weight_out,ngauss_in,position_in,weight_in)
  if(ngauss_in.ne.48) STOP ' no chuto el 48'
  weight(:)=weight_in(:)
  do ikk=1,ngauss
    call BSpline(shapef(ikk,:),position_in(2*ikk-1),position_in(2*ikk))  
  enddo
else
  STOP ' Number of Gauss points not implemented'
end if

END SUBROUTINE gauss_vdw

!*************************************************************
!*************************************************************
SUBROUTINE subdiv_gauss(ngauss_in,position_in,weight_in,ngauss_out,position_out,weight_out)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8) position_in(ngauss_in*2), weight_in(ngauss_in)
REAL(8) position_out(ngauss_in*2*4), weight_out(ngauss_in*4)

ngauss_out=ngauss_in*4

position_out(1:ngauss_in*2)=position_in(:)/2.d0

position_out(ngauss_in*2+1:ngauss_in*4:2)=position_in(1:ngauss_in*2:2)/2.d0+.5d0
position_out(ngauss_in*2+2:ngauss_in*4:2)=position_in(2:ngauss_in*2:2)/2.d0

position_out(ngauss_in*4+1:ngauss_in*6:2)=position_in(1:ngauss_in*2:2)/2.d0
position_out(ngauss_in*4+2:ngauss_in*6:2)=position_in(2:ngauss_in*2:2)/2.d0+.5d0

position_out(ngauss_in*6+1:ngauss_in*8:2)=.5d0-position_in(1:ngauss_in*2:2)/2.d0
position_out(ngauss_in*6+2:ngauss_in*8:2)=.5d0-position_in(2:ngauss_in*2:2)/2.d0

weight_out(1:ngauss_in)=weight_in(1:ngauss_in)/4.d0
weight_out(ngauss_in+1:2*ngauss_in)=weight_in(1:ngauss_in)/4.d0
weight_out(2*ngauss_in+1:3*ngauss_in)=weight_in(1:ngauss_in)/4.d0
weight_out(3*ngauss_in+1:4*ngauss_in)=weight_in(1:ngauss_in)/4.d0

END SUBROUTINE subdiv_gauss

