MODULE VSF

! Реализация векторных сферических функций

IMPLICIT NONE

 REAL(8), PARAMETER ::  PI=3.1415926535897932384626433832795_8

CONTAINS

 REAL(8) FUNCTION FACT(N)
 ! Факториал 
 INTEGER(4), INTENT(IN) :: N
 INTEGER(4) ::  I
 FACT=1.0_8
 DO I=1, N
  FACT=FACT*I
 END DO
 END FUNCTION FACT

 REAL(8) FUNCTION FACT2(N1,N2)
 ! Произведение чисел от N1 до N2 
 INTEGER(4), INTENT(IN) :: N1,N2
 INTEGER(4) ::  I
 FACT2=1.0_8
 DO I=N1, N2
  FACT2=FACT2*I
 END DO
 END FUNCTION FACT2


 REAL(8) RECURSIVE FUNCTION PLR(N,K,D)
 ! Присоединенный полином Лежанжра P_nk (Рекурсивная форма)

  INTEGER(4), INTENT(IN) :: N,K  ! Индексы полинома
  REAL(8), INTENT(IN) :: D       ! Склонение или галактическая широта в градусах!
  
   PLR=0.0
   IF (N<K) RETURN 
   
   IF (N==K) THEN ! P_kk
    PLR=FACT2(K+1,2*K)/2**K*(DCOSD(D)**K)
   ELSE IF (N==K+1) THEN
    PLR=FACT2(K+2,2*K+2)/(2**(K+1))*(DCOSD(D)**K)*DSIND(D)
   ELSE
    PLR=DSIND(D)*DFLOAT(2*N-1)/DFLOAT(N-K)*PLR(N-1,K,D)-DFLOAT(N+K-1)/DFLOAT(N-K)*PLR(N-2,K,D)
   END IF
 END FUNCTION PLR

 REAL(8) FUNCTION PL(N,K,D)
  INTEGER(4), INTENT(IN) :: N,K  ! Индексы полинома
  REAL(8), INTENT(IN) :: D       ! Склонение или галактическая широта в градусах!
  REAL(8) X
  REAL(8) Z
  INTEGER(4) F
  INTEGER M

  PL=0.0
  IF (N<K) RETURN 

  X=DSIND(D)
  
  F=(N-K)/2

  IF (F==0) THEN
    Z= X**(N-K) * DSQRT(1-X**2)**K
  ELSE
    Z=X**(N-K)
	DO M=1,F
	 Z=Z+A(M,N,K)*X**(N-K-2*M)
	END DO
	Z=Z*DSQRT(1-X**2)**K
  END IF
 
  PL=NORM(N,K)*Z
  
  CONTAINS
  
   REAL(8) FUNCTION NORM(N,K)
   INTEGER(4), INTENT(IN) :: N,K
   INTEGER(4) I
   REAL(8) :: Up, Dn

   Up=1.0_8
   DO I=N+1,2*N
    Up=Up*I
   END DO

   Dn=1.0_8
   DO I=1,N-K
    Dn=Dn*I
   END DO
  
   Up=Up/Dn

   DO I=1,N
    Up=Up*0.5
   END DO

   NORM=Up

   END FUNCTION NORM

   REAL(8) FUNCTION A(M,N,K)
 
    INTEGER(4), INTENT(IN) :: M,N,K
	INTEGER(4) :: I
    REAL(8) :: Up, Dn
	
	Up=1.0
    DO I=0,2*M-1
	 Up=Up*(N-K-I)
	END DO
	
	Dn=1
    DO I=1,M
	 Dn=Dn*2*I*(2*(N-I)+1)
	END DO
	
	IF ( MOD(M,2)== 0) THEN
	  A=Up/Dn
    ELSE
	  A=-Up/Dn
    END IF
 
   END FUNCTION A
 
 END FUNCTION PL


 INTEGER(4) FUNCTION INDEXJ(N,K,L)
 ! Вычисление индекса J=N**2+2*K+L-1
 INTEGER, INTENT (IN) :: N,K,L
  INDEXJ=N**2+2*K+L-1
 END FUNCTION INDEXJ

 SUBROUTINE INDEXES(J,N,K,L)
 ! Вычисление индексов N,K,L по индексу J (J=N**2+2*K+L-1)
 INTEGER, INTENT (IN) :: J 
 INTEGER, INTENT (OUT) :: N,K,L
  N=SQRT(FLOAT(J))
  K=J-N**2
  IF (MOD(K,2)==0) THEN
   L=1
  ELSE
   L=0
  ENDIF
  K=(K-L+1)/2
 END SUBROUTINE

 REAL(8) FUNCTION FK(N,K,L,A,D)
 ! Ненормированная сферическая функция
 INTEGER, INTENT (IN) :: N,K,L
 REAL(8), INTENT (IN) :: A,D
 IF (K==0) THEN
   FK=PL(N,0,D)
 ELSE
   IF (L==0) THEN
    FK=PL(N,K,D)*DSIND(K*A)
   ELSE 
    FK=PL(N,K,D)*DCOSD(K*A)
   END IF
 END IF
 END FUNCTION FK


 REAL(8) FUNCTION FR(N,K)
 ! Норма сферической функции
 INTEGER, INTENT (IN) :: N,K
 FR=DSQRT((2*N+1)/(4.0*PI))
 IF (K>0) THEN
  FR=FR*DSQRT(2.0/FACT2(N-K+1,N+K))
 END IF 

 END FUNCTION FR

 REAL(8) FUNCTION FV(N,K,L,A,D)
 ! Нормированная сферическая функция
 INTEGER, INTENT (IN) :: N,K,L
 REAL(8), INTENT (IN) :: A,D ! В градусах

 FV=FR(N,K)*FK(N,K,L,A,D)
 END FUNCTION FV

 REAL(8) FUNCTION FVJ(J,A,D)
 ! Нормированная сферическая функция от одного индекса
 INTEGER, INTENT (IN) :: J
 REAL(8), INTENT (IN) :: A,D ! В градусах
 INTEGER N,K,L
 CALL INDEXES(J,N,K,L)
 FVJ=FV(N,K,L,A,D)
 END FUNCTION FVJ

 REAL(8) FUNCTION TA(N,K,L,A,D)
 ! Компонент торроидальной функции по A
 INTEGER, INTENT (IN) :: N,K,L
 REAL(8), INTENT (IN) :: A,D ! В градусах
 IF (K==0) THEN
  TA=PL(N,1,D)
 ELSE
  IF (L==0) THEN
   TA=(-K*DTAND(D)*PL(N,K,D)+PL(N,K+1,D))*DSIND(K*A)
  ELSE
   TA=(-K*DTAND(D)*PL(N,K,D)+PL(N,K+1,D))*DCOSD(K*A)
  END IF
 END IF 
  TA=FR(N,K)/DSQRT(DFLOAT(N*(N+1)))*TA 
 END FUNCTION TA

 REAL(8) FUNCTION TAJ(J,A,D)
  ! Компонент торроидальной функции по A от одного индекса
 INTEGER, INTENT (IN) :: J
 REAL(8), INTENT (IN) :: A,D ! В градусах
 INTEGER N,K,L
 CALL INDEXES(J,N,K,L)
 TAJ=TA(N,K,L,A,D)
 END FUNCTION TAJ

 REAL(8) FUNCTION TD(N,K,L,A,D)
 ! Компонент торроидальной функции по D
 INTEGER, INTENT (IN) :: N,K,L
 REAL(8), INTENT (IN) :: A,D ! В градусах
 IF (K==0) THEN
  TD=0.0
 ELSE
  IF (L==0) THEN
   TD=-K/DCOSD(D)*PL(N,K,D)*DCOSD(K*A)
  ELSE
   TD=+K/DCOSD(D)*PL(N,K,D)*DSIND(K*A)
  END IF
  TD=FR(N,K)/DSQRT(DFLOAT(N*(N+1)))*TD
 END IF 

 END FUNCTION TD

 REAL(8) FUNCTION TDJ(J,A,D)
  ! Компонент торроидальной функции по D от одного индекса
 INTEGER, INTENT (IN) :: J
 REAL(8), INTENT (IN) :: A,D ! В градусах
 INTEGER N,K,L
 CALL INDEXES(J,N,K,L)
 TDJ=TD(N,K,L,A,D)
 END FUNCTION TDJ

 REAL(8) FUNCTION SA(N,K,L,A,D)
 ! Компонент сфероидальной функции по A
 INTEGER, INTENT (IN) :: N,K,L
 REAL(8), INTENT (IN) :: A,D ! В градусах
 IF (K==0) THEN
  SA=0.0
 ELSE
  IF (L==0) THEN
   SA=+K/DCOSD(D)*PL(N,K,D)*DCOSD(K*A)
  ELSE
   SA=-K/DCOSD(D)*PL(N,K,D)*DSIND(K*A)
  END IF
  SA=FR(N,K)/DSQRT(DFLOAT(N*(N+1)))*SA
 END IF 
 END FUNCTION SA

 REAL(8) FUNCTION SAJ(J,A,D)
  ! Компонент сфероидальной функции по A от одного индекса
 INTEGER, INTENT (IN) :: J
 REAL(8), INTENT (IN) :: A,D ! В градусах
 INTEGER N,K,L
 CALL INDEXES(J,N,K,L)
 SAJ=SA(N,K,L,A,D)
 END FUNCTION SAJ

 REAL(8) FUNCTION SD(N,K,L,A,D)
 ! Компонент сфероидальной функции по D
 INTEGER, INTENT (IN) :: N,K,L
 REAL(8), INTENT (IN) :: A,D ! В градусах
 IF (K==0) THEN
  SD=PL(N,1,D)
 ELSE
  IF (L==0) THEN
   SD=(-K*DTAND(D)*PL(N,K,D)+PL(N,K+1,D))*DSIND(K*A)
  ELSE
   SD=(-K*DTAND(D)*PL(N,K,D)+PL(N,K+1,D))*DCOSD(K*A)
  END IF
 END IF 
   SD=FR(N,K)/DSQRT(DFLOAT(N*(N+1)))*SD
 END FUNCTION SD

 REAL(8) FUNCTION SDJ(J,A,D)
  ! Компонент сфероидальной функции по D от одного индекса
 INTEGER, INTENT (IN) :: J
 REAL(8), INTENT (IN) :: A,D ! В градусах
 INTEGER N,K,L
 CALL INDEXES(J,N,K,L)
 SDJ=SD(N,K,L,A,D)
 END FUNCTION SDJ

END MODULE VSF