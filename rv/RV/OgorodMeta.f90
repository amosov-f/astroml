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


! ===================== 4 НОЯБРЯ 2013 ================================

! Вычисляет базисную функцию модели Огородникова-Милна
REAL(8) FUNCTION VR_base(j,l,b,R)
integer, intent(IN) :: j ! Номер функции
real(8), INTENT (IN)    :: l,b,R ! гал. координаты в градусах и пк

SELECT CASE (j)
CASE (1) ! U
 vr_base = -dcosd(l)*dcosd(b)
CASE (2) ! V
 vr_base = -dsind(l)*dcosd(b)
CASE (3) ! W
 vr_base = -dsind(B)
CASE (4) ! M12
 vr_base = r*dcosd(b)**2*dsind(2*l)
CASE (5) ! M13
 vr_base = r*dsind(2*b)*dcosd(l)
CASE (6) ! M23
 vr_base = r*dsind(2*b)*dsind(l)
CASE (7) ! M11
 vr_base = r*dcosd(b)**2*dcosd(l)**2
CASE (8) ! M22
 vr_base = r*dcosd(b)**2*dsind(l)**2
CASE (9) ! M33
 vr_base = r*dsind(b)**2
END SELECT 
END FUNCTION vr_base