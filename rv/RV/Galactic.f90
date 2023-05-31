MODULE GALACTIC
USE DFLIB

IMPLICIT NONE



 real(8), parameter :: Leo = 282.85948083 ! °
 real(8), parameter :: L0  =  32.931918056 ! ° 
 real(8), parameter :: si  = 0.88998807641_8 ! sin 62.871748611° 
 real(8), parameter :: ci  = 0.45598379779_8 ! cos 62.871748611° 
 
 real(8), parameter :: Pi  = 3.1415926535897932384626433832795_8


CONTAINS





INTEGER(4) FUNCTION Round(x)
REAL(8), INTENT(IN) :: x

 Round=INT(ANINT(x))
 
END FUNCTION



 Subroutine Galaxy(a,d,l,b)
 real(8),intent(in) :: a,d
 real(8),intent(out) :: l,b
 real(8) :: al,sa,ca,sd,cd 

  al=a-Leo
  sa=dsind(al);  ca=dcosd(al)
  sd=dsind(d);   cd=dcosd(d)
  b=dasind(sd*ci-cd*si*sa)
  l=datan2d(sd*si+cd*ci*sa,cd*ca)+L0
  if (l<0) then 
   l=l+360.0_8
  endif
 end subroutine


 Subroutine GalaxMu(mua,mud,l,b,d, mul,mub)
 real(8),intent(in) :: mua,mud,l,b,d 
 real(8),intent(out) :: mul,mub
 real(8) :: cd,sfi,cfi
   cd=dcosd(d);
   sfi=si*dcosd(l-L0)/cd
   cfi=(dcosd(b)*ci-dsind(b)*si*dsind(l-L0))/cd
   mul= cfi*mua+sfi*mud
   mub=-sfi*mua+cfi*mud
 end subroutine

 Subroutine EquatorialMu(mul,mub,l,b,d, mua,mud)
 real(8),intent(in) :: mul,mub,l,b,d 
 real(8),intent(out) :: mua,mud
 real(8) :: cd,sfi,cfi
   cd=dcosd(d);
   sfi=si*dcosd(l-L0)/cd
   cfi=(dcosd(b)*ci-dsind(b)*si*dsind(l-L0))/cd
   mua= cfi*mul-sfi*mub
   mud= sfi*mul+cfi*mub
 end subroutine



 REAL(4) function Ranorm(s)
 ! s - δθροεπρθÿ
 real(4), intent(in) :: s
 real(4) :: a
 integer, save :: z=37843827
 integer(4) ::   i;
  a=0.0
  DO i=1,12 
   a=a+ran(z)
  END DO
  Ranorm=(a-6.0)*s
 end FUNCTION RANORM

END MODULE