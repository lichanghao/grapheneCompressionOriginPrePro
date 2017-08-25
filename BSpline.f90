SUBROUTINE BSpline(N,v,w)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8) :: N(12)
      t1 = 1.D0-v-w
      t2 = t1**2
      t3 = t2**2
      t4 = t2*t1
      t5 = t4*v
      t6 = 2.D0*t5
      t8 = t4*w
      t9 = 2.D0*t8
      t11 = 6.D0*t5
      t13 = t2*v*w
      t14 = 6.D0*t13
      t15 = v**2
      t16 = t2*t15
      t17 = 12.D0*t16
      t19 = t1*t15*w
      t20 = 6.D0*t19
      t21 = t15*v
      t22 = t1*t21
      t23 = 6.D0*t22
      t24 = t21*w
      t25 = 2.D0*t24
      t26 = t15**2
      t30 = w**2
      t31 = t2*t30
      t32 = 24.D0*t31
      t33 = t30*w
      t34 = t1*t33
      t36 = t30**2
      t40 = t1*v*t30
      t41 = 36.D0*t40
      t42 = v*t33
      t43 = 6.D0*t42
      t44 = 24.D0*t16
      t45 = 36.D0*t19
      t46 = t15*t30
      t47 = 12.D0*t46
      t49 = 6.D0*t24
      t50 = 6.D0*t3+24.D0*t8+t32+8.D0*t34+t36+24.D0*t5+60.D0*t13+t41+t43 &
           +t44+t45+t47+8.D0*t22+t49+t26
      t51 = 6.D0*t8
      t52 = 12.D0*t31
      t53 = 6.D0*t34
      t54 = 6.D0*t40
      t55 = 2.D0*t42
      t57 = 2.D0*t22
      t60 = 36.D0*t13
      t63 = 24.D0*t46
      t67 = t3+t51+t52+t53+t36+8.D0*t5+t60+t41+8.D0*t42+t44+60.D0*t19+ &
	        t63+24.D0*t22+24.D0*t24+6.D0*t26
      t74 = t3+8.D0*t8+t32+24.D0*t34+6.D0*t36+t11+t60+60.D0*t40+24.D0*t42 &
	        +t17+t45+t63+t23+8.D0*t24+t26
      t75 = 2.D0*t34
      N(1) = t3+t6
      N(2) = t3+t9
      N(3) = t3+t9+t11+t14+t17+t20+t23+t25+t26
      N(4) = t50
      N(5) = t3+t51+t52+t53+t36+t6+t14+t54+t55
      N(6) = t57+t26
      N(7) = t67
      N(8) = t74
      N(9) = t75+t36
      N(10) = t25+t26
      N(11) = t75+t36+t54+t43+t20+t47+t57+t49+t26
      N(12) = t36+t55

     N=N/12.d0

END SUBROUTINE BSpline


SUBROUTINE DBSpline(DN,v,w)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8) :: DN(12,2)

      t1 = 1.D0-v-w
      t2 = t1**2
      t3 = t2*t1
      t4 = 2.D0*t3
      t5 = t2*v
      t6 = 6.D0*t5
      t8 = 4.D0*t3
      t10 = t2*w
      t11 = 6.D0*t10
      t14 = v**2
      t15 = t1*t14
      t16 = 6.D0*t15
      t17 = t14*v
      t18 = 2.D0*t17
      t20 = 12.D0*t5
      t22 = t1*v*w
      t23 = 12.D0*t22
      t25 = t14*w
      t26 = 6.D0*t25
      t27 = 4.D0*t17
      t29 = 12.D0*t10
      t30 = w**2
      t31 = t1*t30
      t32 = 12.D0*t31
      t33 = t30*w
      t34 = 2.D0*t33
      t35 = 24.D0*t5
      t36 = 48.D0*t22
      t37 = v*t30
      t38 = 12.D0*t37
      t39 = 24.D0*t15
      t40 = 18.D0*t25
      t42 = 24.D0*t10
      t43 = 24.D0*t31
      t44 = 4.D0*t33
      t45 = 18.D0*t37
      t46 = 12.D0*t15
      t47 = 12.D0*t25
      t50 = 6.D0*t37
      t52 = 6.D0*t31
      DN(1,1) = -t4-t6
      DN(1,2) = -t8-t6
      DN(2,1) = -t8-t11
      DN(2,2) = -t4-t11
      DN(3,1) = t4+t6-t16-t18
      DN(3,2) = -t4-t11-t20-t23-18.D0*t15-t26-t27
      DN(4,1) = -t29-t32-t34-t35-t36-t38-t39-t40-t27
      DN(4,2) = -t42-t43-t44-t20-t36-t45-t46-t47-t18
      DN(5,1) = -t4-t29-18.D0*t31-t44-t6-t23-t50
      DN(5,2) = t4+t11-t52-t34
      DN(6,1) = t18+t16
      DN(6,2) = -t18
      DN(7,1) = t8+18.D0*t10+t32+t34+t35+t36+t38+t39+t47
      DN(7,2) = t4+t11-t52-t34+t20-t38+t46-t47
      DN(8,1) = t4+t29+t32+t6-t38-t16-t47-t18
      DN(8,2) = t8+t42+t43+18.D0*t5+t36+t38+t46+t47+t18
      DN(9,1) = -t34
      DN(9,2) = t34+t52
      DN(10,1) = t26+t27
      DN(10,2) = t18
      DN(11,1) = t44+t45+t52+t47+t23+t18+t16
      DN(11,2) = t34+t52+t38+t23+t40+t16+t27
      DN(12,1) = t34
      DN(12,2) = t44+t50


     DN=DN/12.d0


END SUBROUTINE DBSpline

SUBROUTINE DDBSpline(DDN,v,w)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8) :: DDN(12,3)
!1 is v,v
!2 is w,w
!3 is v,w or w,v

      t1 = 1.D0-v-w
      t2 = t1*v
      t3 = 12.D0*t2
      t4 = t1**2
      t6 = 6.D0*t4
      t8 = t1*w
      t10 = 12.D0*t8
      t12 = 24.D0*t2
      t13 = v*w
      t14 = v**2
      t15 = t8+t2+t13+t14
      t16 = 6.D0*t14
      t18 = 24.D0*t8
      t19 = 24.D0*t4
      t20 = 12.D0*t13
      t21 = 12.D0*t14
      t23 = w**2
      t24 = 12.D0*t23
      t26 = 12.D0*t4
      t27 = 6.D0*t23
      t28 = 24.D0*t13
      t30 = t8+t23+t2+t13
      DDN(1,1) = t3
      DDN(1,2) = 12.D0*t4+12.D0*t2
      DDN(1,3) = t6+t3
      DDN(2,1) = 12.D0*t4+12.D0*t8
      DDN(2,2) = t10
      DDN(2,3) = t6+t10
      DDN(3,1) = -t12
      DDN(3,2) = 12.D0*t15
      DDN(3,3) = -t6-t3+t16
      DDN(4,1) = -t18-t19+t20+t21
      DDN(4,2) = -t19+t24-t12+t20
      DDN(4,3) = -t26+t27+t28+t16
      DDN(5,1) = 12.D0*t30
      DDN(5,2) = -t18
      DDN(5,3) = -t6-t10+t27
      DDN(6,1) = t3
      DDN(6,2) = 0.D0
      DDN(6,3) = -t16
      DDN(7,1) = t26+t10-t28-24.D0*t14
      DDN(7,2) = -24.D0*t15
      DDN(7,3) = t6-t10-t27-t28-t21
      DDN(8,1) = -24.D0*t30
      DDN(8,2) = t26-24.D0*t23+t3-t28
      DDN(8,3) = t6-t24-t3-t28-t16
      DDN(9,1) = 0.D0
      DDN(9,2) = t10
      DDN(9,3) = -t27
      DDN(10,1) = 12.D0*t13+12.D0*t14
      DDN(10,2) = 0.D0
      DDN(10,3) = t16
      DDN(11,1) = 12.D0*t30
      DDN(11,2) = 12.D0*t15
      DDN(11,3) = t27+t28+t10+t16+t3
      DDN(12,1) = 0.D0
      DDN(12,2) = 12.D0*t23+12.D0*t13
      DDN(12,3) = t27

     DDN=DDN/12.d0


END SUBROUTINE DDBSpline


