Program GetOgorod

! Определение параметров модели Огородинкова-Милна по с.д. звезд Hipparcos

Use LSQ
Use Ogorod
Use GALACTIC

IMPLICIT NONE


  CHARACTER(*), PARAMETER :: NAME = 'Sigma=10' 

  integer, parameter :: N = 9 ! количество коэффициентов уравнений
 
  real(8) :: RA,DE,DIST,HVr,HVr_R
  real(8) :: L,B
  INTEGER(4) :: NUM


  real(8),allocatable, dimension (:,:) :: a   ! матрица системы 
  real(8),allocatable, dimension (:) :: y, w  ! правые части и веса
  
  real(8), dimension (N,N) :: r   ! корреляционная матрица
  real(8) :: Cond ! Число обусловленности1
  real(8) :: s0 ! ошибка единицы веса
  real(8) :: MaxCor(N) ! Максимальные коэф. корелляции
  integer :: IndCor(N) ! Номер максимального коэф.
  real(8), dimension (N) :: v,sv ! Ответы и их ошибки
  character(27), parameter :: XZ = "U  V  W  M12M13M23M11M22M33"


  integer :: M ! количество уравнений и удвоенное количество

  integer :: i,j,jf

  ! Считывание из первой строчки числа звезд
  OPEN(UNIT=1,FILE=TRIM(NAME) // '.txt')
  READ(1,"(I10)") M
  PRINT *, 'Number of plates', M
  IF (M==0) STOP "M is zero"

  allocate(a(M,N), Y(M), W(M) )
  W=1.0

 ! Нахождения параметров модели Огородникова-Милна

   i=0 
   DO WHILE (I<M)
    READ(1,"(6X,2F10.4,20X,F10.4,F10.4,F12.4,40X,I8)") RA,DE,DIST,HVr,HVr_R,Num
	IF (NUM==0) CYCLE
    I=I+1

    DIST=1.0_8 

    CALL GALAXY(RA,DE,L,B)

    DO J=1,N
	  IF (j>3) THEN
	   Jf=J+3
      ELSE
	   Jf=J
	  END IF
	  A(i  ,j)  =vrpx_base(Jf,L,B,1.0/Dist)
    END DO

    Y(i)    = HVr_R
   END DO !  WHILE
   
   CLOSE(UNIT=1)
  
   PRINT *,I

  Call LSQM(a,y,w, V,SV, s0, r, Cond)
  PRINT *,'Parameters has been found'



  MaxCor=0.0
  do i=1,N
   do j=1,N
    if (i/=j) then
     if (abs(r(i,j))>abs(MaxCor(i))) then
       MaxCor(i)=r(i,j)!
	   IndCor(i)=j
     endif
    endif 
   enddo
  enddo   
  
  
  OPEN(unit=1,FILE=TRIM(NAME) // '.Vr1')
  DO J=1,N
   WRITE(1,"(A3,2F9.3,F7.2,2X,A3))") XZ(3*j-2:3*j), V(J),SV(J),MaxCor(j),XZ(3*IndCor(j)-2:3*IndCor(j))
  END DO
  
  WRITE(1,"('Cond=',F12.3)") COND
  
  CLOSE(unit=1)

  Deallocate(A, Y, W)  

  
 END Program
