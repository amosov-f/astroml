MODULE OGOROD

! Функции модели Огородникова-Милна 
! Матрица деформации - только члены M^+

IMPLICIT NONE


Type  TVel  ! Скорость движения Солнца в [км/с]
  Sequence
  Real(4) :: X,Y,Z
End Type    

Type  TUVW  ! Скорость движения Солнца
  Sequence
  Real(8) :: U,V,W
End Type    


Type  TRot  ! Угловая скорость
  Sequence
  Real(8) :: X,Y,Z
End Type    

Type  TDef  ! Матрица деформации M^+
  Sequence
  Real(4) :: M12, M13, M23, M11, M22, M33
End Type    


Type  TDef1  ! Матрица деформации M^+ (модифицированная)
  Sequence
  Real(8) :: M12, M13, M23, Ms, Mx, My
End Type    


REAL(8), PARAMETER :: PI = 3.1415926535897932384626433832795

 real(4), parameter, private :: K=4.738
 real(4), parameter :: KK=4.738

character(33), parameter :: UVW = "U  V  W  W1 W2 W3 M13M23M12M11X  "

PRIVATE R, LAM, TAU, P, Q

CONTAINS

! Вычисляет лучевую скорость, умноженную на параллакс в мсд в модели Огородникова-Милна
! модификация (13)
REAL(4) FUNCTION VrPx_1(S,D,l,b,px)
Type(TUVW),  INTENT (IN) :: S ! Компоненты движения солнца
Type(TDef1), INTENT (IN) :: D ! Комопненты матрицы деформации
real(8), INTENT (IN)    :: l,b,px

 VrPx_1=-S.U*px*cosd(l)*cosd(b)-S.V*px*sind(l)*cosd(B)-S.W*px*sind(B) &
     +D.M13*sind(2*b)*cosd(l)+D.M23*sind(2*b)*sind(l)+D.M12*cosd(b)**2*sind(2*l)+ &
     +0.5*D.Ms*cosd(b)**2*cosd(2*l)+D.Mx*(sind(b)**2-1.0/3.0)+D.My
END FUNCTION VrPx_1

! Вычисляет собственное движение по долготе в км/с/кпк в модели Огородникова-Милна
! модификация (11)
REAL(8) FUNCTION kmul_1(S,W,D,l,b,px)
Type(TUVW), INTENT (IN) :: S ! Компоненты движения солнца
Type(TRot), INTENT (IN) :: W ! Комопненты вектора вращения
Type(TDef1),INTENT (IN) :: D ! Комопненты матрицы деформации
real(8), INTENT (IN)    :: l,b,px

 kmul_1=S.U*px*dsind(l)-S.V*px*dcosd(l) &
    -W.X*dsind(b)*dcosd(l)-W.Y*dsind(b)*dsind(l)+W.Z*dcosd(b) &
	-D.M13*dsind(b)*dsind(l)+D.M23*dsind(b)*dcosd(l)+D.M12*dcosd(b)*dcosd(2*l) &
	-0.5*D.Ms*dcosd(b)*dsind(2*l)
END FUNCTION kmul_1

! Вычисляет собственное движение по долготе в км/с/кпк в модели Огородникова-Милна
! модификация
REAL(8) FUNCTION kmub_1(S,W,D,l,b,px)
Type(TUVW), INTENT (IN) :: S ! Компоненты движения солнца
Type(TRot), INTENT (IN) :: W ! Комопненты вектора вращения
Type(TDef1),INTENT (IN) :: D ! Комопненты матрицы деформации
real(8), INTENT (IN)    :: l,b,px

 kmub_1=S.U*px*dcosd(l)*dsind(b)+S.V*px*dsind(l)*dsind(b)-S.W*px*dcosd(b) &
    +W.X*dsind(l)-W.Y*dcosd(l) &
	-0.5*D.M12*dsind(2*b)*dsind(2*l)+D.M13*dcosd(2*b)*dcosd(l)+D.M23*dcosd(2*b)*dsind(l) &
	-0.5*D.Ms*dsind(2*b)*dcosd(l)**2+0.5*D.Mx*dsind(2*b)

END FUNCTION kmub_1


! Вычисляет лучевую скорость, умноженную на параллакс в мсд в модели Огородникова-Милна
REAL(4) FUNCTION Vr(S,D,l,b,px)
Type(TUVW), INTENT (IN) :: S ! Компоненты движения солнца
Type(TDef), INTENT (IN) :: D ! Комопненты матрицы деформации
real(8), INTENT (IN)    :: l,b,px

 Vr=-S.U*cosd(l)*cosd(b)-S.V*sind(l)*cosd(B)-S.W*sind(B) &
     +(D.M13*sind(2*b)*cosd(l)+D.M23*sind(2*b)*sind(l)+D.M12*cosd(b)**2*sind(2*l)+ &
     +D.M11*cosd(b)**2*cosd(l)**2+D.M22*cosd(b)**2*sind(l)**2+D.M33*sind(b)**2)/px

END FUNCTION Vr



! Вычисляет лучевую скорость, умноженную на параллакс в мсд в модели Огородникова-Милна
REAL(4) FUNCTION VrPx(S,D,l,b,px)
Type(TUVW), INTENT (IN) :: S ! Компоненты движения солнца
Type(TDef), INTENT (IN) :: D ! Комопненты матрицы деформации
real(8), INTENT (IN)    :: l,b,px

 VrPx=-S.U*px*cosd(l)*cosd(b)-S.V*px*sind(l)*cosd(B)-S.W*px*sind(B) &
     +D.M13*sind(2*b)*cosd(l)+D.M23*sind(2*b)*sind(l)+D.M12*cosd(b)**2*sind(2*l)+ &
     +D.M11*cosd(b)**2*cosd(l)**2+D.M22*cosd(b)**2*sind(l)**2+D.M33*sind(b)**2

END FUNCTION VrPx


! Вычисляет собственное движение по долготе в км/с/кпк в модели Огородникова-Милна
REAL(4) FUNCTION kmul(S,W,D,l,b,px)
Type(TUVW), INTENT (IN) :: S ! Компоненты движения солнца
Type(TRot), INTENT (IN) :: W ! Комопненты вектора вращения
Type(TDef), INTENT (IN) :: D ! Комопненты матрицы деформации
real(8), INTENT (IN)    :: l,b,px

 kmul=S.U*px*sind(l)-S.V*px*cosd(l) &
    -W.X*sind(b)*cosd(l)-W.Y*sind(b)*sind(l)+W.Z*cosd(b) &
	-D.M13*sind(b)*sind(l)+D.M23*sind(b)*cosd(l)+D.M12*cosd(b)*cosd(2*l) &
	-0.5*D.M11*cosd(b)*sind(2*l)+0.5*D.M22*cosd(b)*sind(2*l)
 
END FUNCTION kmul

! Вычисляет собственное движение по долготе в км/с/кпк в модели Огородникова-Милна
REAL(4) FUNCTION kmub(S,W,D,l,b,px)
Type(TUVW), INTENT (IN) :: S ! Компоненты движения солнца
Type(TRot), INTENT (IN) :: W ! Комопненты вектора вращения
Type(TDef), INTENT (IN) :: D ! Комопненты матрицы деформации
real(8), INTENT (IN)    :: l,b,px

 kmub=S.U*px*cosd(l)*sind(b)+S.V*px*sind(l)*sind(b)-S.W*px*cosd(b) &
    +W.X*sind(l)-W.Y*cosd(l) &
	-0.5*D.M12*sind(2*b)*sind(2*l)+D.M13*cosd(2*b)*cosd(l)+D.M23*cosd(2*b)*sind(l) &
	-0.5*D.M11*sind(2*b)*cosd(l)**2-0.5*D.M22*sind(2*b)*sind(l)**2+0.5*D.M33*sind(2*b)

END FUNCTION kmub


real(4) function solar(kind,j,l,b,px)
integer, intent(in) :: kind ! Vr*px =0, mul, mub, vr
integer, intent(in) :: j    ! U (=1) , V (=2), W (=3)
real(8), intent(in) :: l,b
real(4), intent(in) :: px

solar=0.0

SELECT CASE(kind)
 CASE (0) ! Vr*px
   select case(j)
   case(1) 
    solar=-px*cosd(l)*cosd(b)
   case(2) 
	solar=-px*sind(l)*cosd(b)
   case(3)
    solar=-px*sind(b)
   end select 
 
 CASE (1) ! mul
   select case(j)
   case(1) 
    solar=+px*sind(l)
   case(2) 
	solar=-px*cosd(l)
   case(3)
    solar=0.0
   end select

 CASE (2) ! mub
   select case(j)
   case(1) 
    solar=+px*cosd(l)*sind(b)
   case(2) 
	solar=+px*sind(l)*sind(b)
   case(3)
    solar=-px*cosd(b)
   end select 

 CASE (3) ! Vr
   select case(j)
   case(1) 
    solar=-cosd(l)*cosd(b)
   case(2) 
	solar=-sind(l)*cosd(b)
   case(3)
    solar=-sind(b)
   end select 
   

END SELECT

END FUNCTION SOLAR



 REAL(4) FUNCTION SQVR(i,j,k,l,b,px)
 ! Возвращает возвращает функцию при d M_i,j/d r_k
 integer, intent(in) :: i,j,k
 real(8), intent(in) :: l,b
 real(8), intent(in) :: px

 SQVR=R(k,l,b)*R(i,l,b)*R(j,l,b)
 IF (I/=J) SQVR=2.0*SQVR 

 END FUNCTION SQVR

 REAL(4) FUNCTION SQPMROT_L(i,k,l,b,px)
 ! Возвращает возвращает функцию при d W_i/d r_k
 integer, intent(in) :: i,k
 real(8), intent(in) :: l,b
 real(8), intent(in) :: px

 SQPMROT_L=-R(k,l,b)*LAM(i,l,b)
 IF (i==3) SQPMROT_L=-SQPMROT_L
  
END FUNCTION SQPMROT_L


REAL(4) FUNCTION SQPMROT_B(i,k,l,b,px)
 ! Возвращает возвращает функцию при d W_i/d r_k
 integer, intent(in) :: i,k
 real(8), intent(in) :: l,b
 real(8), intent(in) :: px

 IF (I==3) THEN 
   SQPMROT_B=0.0
 ELSE
  SQPMROT_B=R(k,l,b)*TAU(i,l,b)
  IF (i==2) SQPMROT_B=-SQPMROT_B
 ENDIF
  
END FUNCTION SQPMROT_B


REAL(4) FUNCTION SQPMDEF_L(i,j,k,l,b,px)
 ! Возвращает возвращает функцию при d W_i/d r_k
 integer, intent(in) :: i,j,k
 real(8), intent(in) :: l,b
 real(8), intent(in) :: px

 SQPMDEF_L=R(k,l,b)*P(i,j,l,b)
   
END FUNCTION SQPMDEF_L


REAL(4) FUNCTION SQPMDEF_B(i,j,k,l,b,px)
 ! Возвращает возвращает функцию при d W_i/d r_k
 integer, intent(in) :: i,j,k
 real(8), intent(in) :: l,b
 real(8), intent(in) :: px

  SQPMDEF_B=R(k,l,b)*Q(i,j,l,b)
  
END FUNCTION SQPMDEF_B


 REAL(4) FUNCTION R(k,l,b)
  integer, intent(in) :: k
  real(8), intent(in) :: l,b
  SELECT CASE(k)
   CASE (1) 
    R=cosd(b)*cosd(l)
   CASE (2)
    R=cosd(b)*sind(l)
   CASE (3)
    R=sind(b)
  END SELECT  
 END FUNCTION R

 REAL(4) FUNCTION LAM(k,l,b)
  integer, intent(in) :: k
  real(8), intent(in) :: l,b
  SELECT CASE(k)
   CASE (1) 
    LAM=sind(b)*cosd(l)
   CASE (2)
    LAM=sind(b)*sind(l)
   CASE (3)
    LAM=cosd(b)
  END SELECT  
  END FUNCTION LAM

 REAL(4) FUNCTION TAU(k,l,b)
  integer, intent(in) :: k
  real(8), intent(in) :: l,b
  SELECT CASE(k)
   CASE (1) 
    TAU=sind(l)
   CASE (2)
    TAU=cosd(l)
   CASE (3)
    TAU=0.0
  END SELECT  
  END FUNCTION TAU


 REAL(4) FUNCTION P(i,j,l,b)
  integer, intent(in) :: i,j
  real(8), intent(in) :: l,b
  P=0.0
  SELECT CASE(i)
   CASE (1) 
    SELECT CASE(j)
	 CASE(1)
	  P=-0.5*cosd(b)*sind(2*l)
	 CASE(2)
	  P=cosd(b)*cosd(2*l)
	 CASE(3)
	  P=-sind(b)*sind(l)
	END SELECT
   CASE (2)
    SELECT CASE(j)
	 CASE(2)
	  P=0.5*cosd(b)*sind(2*l)
	 CASE(3)
	  P=sind(b)*cosd(l)
	END SELECT
  END SELECT  
  END FUNCTION P

 REAL(4) FUNCTION Q(i,j,l,b)
  integer, intent(in) :: i,j
  real(8), intent(in) :: l,b
  Q=0.0
  SELECT CASE(i)
   CASE (1) 
    SELECT CASE(j)
	 CASE(1)
	  Q=-0.5*sind(2*b)*(cosd(l)**2)
	 CASE(2)
	  Q=-0.5*sind(2*b)*sind(2*l)
	 CASE(3)
	  Q=cosd(2*b)*cosd(l)
	END SELECT
   CASE (2)
    SELECT CASE(j)
	 CASE(2)
	  Q=-0.5*sind(2*b)*(sind(l)**2)
	 CASE(3)
	  Q=cosd(2*b)*sind(l)
	END SELECT
   CASE (3)
    SELECT CASE(j)
	 CASE(3)
	  Q=0.5*sind(2*b)
	END SELECT
    
  END SELECT  
  END FUNCTION Q



! Вычисляет базисную функцию модели Огородникова-Милна
REAL(4) FUNCTION kmul_base(j,l,b,px)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b,px ! гал. координаты в градусах и mas

SELECT CASE (j)
CASE (1) ! U
 kmul_base =+px*sind(l)
CASE (2) ! V
 kmul_base =-px*cosd(l)
CASE (3) ! W
 kmul_base = 0.0
CASE (4) ! Wx
 kmul_base =-sind(b)*cosd(l)
CASE (5) ! Wy
 kmul_base =-sind(b)*sind(l)
CASE (6) ! Wz
 kmul_base =+cosd(b)
CASE (7) ! M12
 kmul_base =+cosd(b)*cosd(2*l)
CASE (8) ! M13
 kmul_base =-sind(b)*sind(l)
CASE (9) ! M23
 kmul_base =+sind(b)*cosd(l) 
CASE (10) ! M11
 kmul_base = -0.5*cosd(b)*sind(2*l)
CASE (11) ! M22
 kmul_base = +0.5*cosd(b)*sind(2*l)
CASE (12) ! M33
 kmul_base = 0.0
CASE DEFAULT
 kmul_base = 0.0
END SELECT 
END FUNCTION kmul_base

! Вычисляет базисную функцию модели Огородникова-Милна
REAL(4) FUNCTION kmub_base(j,l,b,px)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b,px ! гал. координаты в градусах и mas

SELECT CASE (j)
CASE (1) ! U
 kmub_base = +px*cosd(l)*sind(b)
CASE (2) ! V
 kmub_base = +px*sind(l)*sind(b)
CASE (3) ! W
 kmub_base = -px*cosd(b)
CASE (4) ! Wx
 kmub_base =  sind(l)
CASE (5) ! Wy
 kmub_base = -cosd(l)
CASE (6) ! Wz
 kmub_base = 0.0
CASE (7) ! M12
 kmub_base = -0.5*sind(2*b)*sind(2*l)
CASE (8) ! M13
 kmub_base = +cosd(2*b)*cosd(l)
CASE (9) ! M23
 kmub_base = +cosd(2*b)*sind(l)
CASE (10) ! M11
 kmub_base = -0.5*sind(2*b)*cosd(l)**2
CASE (11) ! M22
 kmub_base = -0.5*sind(2*b)*sind(l)**2
CASE (12) ! M33
 kmub_base = +0.5*sind(2*b)
CASE DEFAULT
 kmub_base = 0.0
END SELECT 
END FUNCTION kmub_base

 
! Вычисляет базисную функцию модели Огородникова-Милна
REAL(4) FUNCTION vrpx_base(j,l,b,px)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b,px ! гал. координаты в градусах и mas

SELECT CASE (j)
CASE (1) ! U
 vrpx_base = -cosd(l)*cosd(b)
CASE (2) ! V
 vrpx_base = -sind(l)*cosd(b)
CASE (3) ! W
 vrpx_base = -sind(B)
CASE (4) ! Wx
 vrpx_base = 0.0
CASE (5) ! Wy
 vrpx_base = 0.0
CASE (6) ! Wz
 vrpx_base = 0.0
CASE (7) ! M12
 vrpx_base = cosd(b)**2*sind(2*l)/px
CASE (8) ! M13
 vrpx_base = sind(2*b)*cosd(l)/px
CASE (9) ! M23
 vrpx_base = sind(2*b)*sind(l)/px
CASE (10) ! M11
 vrpx_base =  cosd(b)**2*cosd(l)**2/px
CASE (11) ! M22
 vrpx_base = cosd(b)**2*sind(l)**2/px
CASE (12) ! M33
 vrpx_base = sind(b)**2/px
END SELECT 
END FUNCTION vrpx_base



! Вычисляет базисную функцию модели Огородникова-Милна
! с 11 параметрами
REAL(8) FUNCTION kmul_base_11(j,l,b,px)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b,px ! гал. координаты в градусах и mas

SELECT CASE (j)
CASE (1) ! U
 kmul_base_11 =+px*dsind(l)
CASE (2) ! V
 kmul_base_11 =-px*dcosd(l)
CASE (3) ! W
 kmul_base_11 = 0.0
CASE (4) ! Wx
 kmul_base_11 =-dsind(b)*dcosd(l)
CASE (5) ! Wy
 kmul_base_11 =-dsind(b)*dsind(l)
CASE (6) ! Wz
 kmul_base_11 =+dcosd(b)
CASE (7) ! M13
 kmul_base_11 =-dsind(b)*dsind(l)
CASE (8) ! M23
 kmul_base_11 =+dsind(b)*dcosd(l) 
CASE (9) ! M12
 kmul_base_11 =+dcosd(b)*dcosd(2*l)
CASE (10) ! M11*
 kmul_base_11 = -0.5*dcosd(b)*dsind(2*l)
CASE (11) ! M33
 kmul_base_11 = 0.0
CASE DEFAULT
 kmul_base_11 = 0.0
END SELECT 
END FUNCTION kmul_base_11

! Вычисляет базисную функцию модели Огородникова-Милна
REAL(8) FUNCTION kmub_base_11(j,l,b,px)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b,px ! гал. координаты в градусах и mas

SELECT CASE (j)
CASE (1) ! U
 kmub_base_11 = +px*dcosd(l)*dsind(b)
CASE (2) ! V
 kmub_base_11 = +px*dsind(l)*dsind(b)
CASE (3) ! W
 kmub_base_11 = -px*dcosd(b)
CASE (4) ! Wx
 kmub_base_11 =  dsind(l)
CASE (5) ! Wy
 kmub_base_11 = -dcosd(l)
CASE (6) ! Wz
 kmub_base_11 = 0.0
CASE (7) ! M13
 kmub_base_11 = +dcosd(2*b)*dcosd(l)
CASE (8) ! M23
 kmub_base_11 = +dcosd(2*b)*dsind(l)
CASE (9) ! M12
 kmub_base_11 = -0.5*dsind(2*b)*dsind(2*l)
CASE (10) ! M11*
 kmub_base_11 = -0.5*dsind(2*b)*dcosd(l)**2
CASE (11) ! M33
 kmub_base_11 = +0.5*dsind(2*b)
CASE DEFAULT
 kmub_base_11 = 0.0
END SELECT 
END FUNCTION kmub_base_11


! Вычисляет модифицированную базисную функцию модели Огородникова-Милна
! с 11 параметрами
REAL(8) FUNCTION kmul_base_11_mod(j,l,b)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b ! гал. координаты в градусах и mas

SELECT CASE (j)
CASE (1) ! U
 kmul_base_11_mod = dsind(l)*dsind(dabs(b))
CASE (2) ! V
 kmul_base_11_mod =-dcosd(l)*dsind(dabs(b))
CASE (3) ! W
 kmul_base_11_mod = 0.0
CASE (4) ! Wx
 kmul_base_11_mod =-dsind(b)*dcosd(l)
CASE (5) ! Wy
 kmul_base_11_mod =-dsind(b)*dsind(l)
CASE (6) ! Wz
 kmul_base_11_mod =+dcosd(b)
CASE (7) ! M13
 kmul_base_11_mod =-dsind(b)*dsind(l)
CASE (8) ! M23
 kmul_base_11_mod =+dsind(b)*dcosd(l) 
CASE (9) ! M12
 kmul_base_11_mod =+dcosd(b)*dcosd(2*l)
CASE (10) ! M11*
 kmul_base_11_mod = -0.5*dcosd(b)*dsind(2*l)
CASE (11) ! M33
 kmul_base_11_mod = 0.0
CASE DEFAULT
 kmul_base_11_mod = 0.0
END SELECT 
END FUNCTION kmul_base_11_mod

! Вычисляет модифицированную базисную функцию модели Огородникова-Милна
REAL(8) FUNCTION kmub_base_11_mod(j,l,b)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b ! гал. координаты в градусах и mas

SELECT CASE (j)
CASE (1) ! U
 kmub_base_11_mod = dcosd(l)*dsind(b)*dsind(dabs(b))
CASE (2) ! V
 kmub_base_11_mod = dsind(l)*dsind(b)*dsind(dabs(b))
CASE (3) ! W
 kmub_base_11_mod =-dcosd(b)*dsind(dabs(b))
CASE (4) ! Wx
 kmub_base_11_mod =  dsind(l)
CASE (5) ! Wy
 kmub_base_11_mod = -dcosd(l)
CASE (6) ! Wz
 kmub_base_11_mod = 0.0
CASE (7) ! M13
 kmub_base_11_mod = +dcosd(2*b)*dcosd(l)
CASE (8) ! M23
 kmub_base_11_mod = +dcosd(2*b)*dsind(l)
CASE (9) ! M12
 kmub_base_11_mod = -0.5*dsind(2*b)*dsind(2*l)
CASE (10) ! M11*
 kmub_base_11_mod = -0.25*dsind(2*b)*dcosd(2*l)**2
CASE (11) ! M33
 kmub_base_11_mod = +0.5*dsind(2*b)
CASE DEFAULT
 kmub_base_11_mod = 0.0
END SELECT 
END FUNCTION kmub_base_11_mod


! ============= 1 ИЮЛЯ 2012 ======================================================

! Вычисляет базисную функцию модели Огородникова-Милна с 11 параметрами 
REAL(8) FUNCTION kmul_base_prim(j,l,b,px)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b ! гал. координаты в градусах 
real(4), INTENT (IN)    :: px !  mas

SELECT CASE (j)
CASE (1) ! U
 kmul_base_prim =+px*dsind(l)
CASE (2) ! V
 kmul_base_prim =-px*dcosd(l)
CASE (3) ! W
 kmul_base_prim = 0.0
CASE (4) ! Wx
 kmul_base_prim =-dsind(b)*dcosd(l)
CASE (5) ! Wy
 kmul_base_prim =-dsind(b)*dsind(l)
CASE (6) ! Wz
 kmul_base_prim =+dcosd(b)
CASE (7) ! M13+
 kmul_base_prim =-dsind(b)*dsind(l)
CASE (8) ! M23+
 kmul_base_prim =+dsind(b)*dcosd(l) 
CASE (9) ! M12+
 kmul_base_prim =+dcosd(b)*dcosd(2*l)
CASE (10) ! M11*
 kmul_base_prim = -0.5*dcosd(b)*dsind(2*l)
CASE (11) ! M33
 kmul_base_prim = 0.0
CASE DEFAULT
 kmul_base_prim = 0.0
END SELECT 
END FUNCTION kmul_base_prim

! Вычисляет базисную функцию модели Огородникова-Милна
REAL(8) FUNCTION kmub_base_prim(j,l,b,px)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b ! гал. координаты в градусах и mas
real(4), INTENT (IN)    :: px !  mas

SELECT CASE (j)
CASE (1) ! U
 kmub_base_prim = +px*dcosd(l)*dsind(b)
CASE (2) ! V
 kmub_base_prim = +px*dsind(l)*dsind(b)
CASE (3) ! W
 kmub_base_prim = -px*dcosd(b)
CASE (4) ! Wx
 kmub_base_prim =  dsind(l)
CASE (5) ! Wy
 kmub_base_prim = -dcosd(l)
CASE (6) ! Wz
 kmub_base_prim = 0.0
CASE (7) ! M13
 kmub_base_prim = +dcosd(2*b)*dcosd(l)
CASE (8) ! M23
 kmub_base_prim = +dcosd(2*b)*dsind(l)
CASE (9) ! M12
 kmub_base_prim = -0.5*dsind(2*b)*dsind(2*l)
CASE (10) ! M11*
 kmub_base_prim = -0.25*dsind(2*b)*dcosd(2*l)
CASE (11) ! X
 kmub_base_prim = +0.5*dsind(2*b)
CASE DEFAULT
 kmub_base_prim = 0.0
END SELECT 
END FUNCTION kmub_base_prim


END MODULE OGOROD