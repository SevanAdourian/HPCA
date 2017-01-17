include	"modules.f90"
include "nrtype.f90"

recursive subroutine kenpsv(n,m, k, Ruac, Tuac, Rdac, Tdac)

  use constants
  implicit none

  integer 		:: n, m
  real			:: k
  real, dimension(2,2)	:: Ruac, Tuac, Rdac, Tdac

  !local variables
  real, dimension(2,2)	:: dummy, TEMP
  real, dimension(2,2)	:: Ruab, Tuab, Rdab, Tdab
  real, dimension(2,2)	:: Rubc, Tubc, Rdbc, Tdbc
  real, dimension(2,2)	:: Rum,  Tum,  Rdm,  Tdm, Emm1

  !	print *, n, m

  if( m < n .or. m < 1 .or. n < 0 ) then
     print *, "function 'kenpsv' called with n=%d and m=%d\n", n, m
     print *, "Error: 'n' should be '>=' than '0', \n"
     print *, "Error: 'm' should be '>=' than '1'  and \n"
     print *, "Error: 'n' should be '<=' than 'm' \n"
     stop
  endif

  if ((m == n) .or. (m == 1 .and. n == 0)) then
     Tdac  = Id2
     Rdac  = Ze2
     Tuac  = Id2
     Ruac  = Ze2
     if(n == 0) then
	!               when called with n=0; only 'Ruac' is 
        call RERPSV(1,dummy,Ruac)
     endif
  else
     call ERTUDPSV(m, k, Emm1,Rum,Tum,Rdm,Tdm)

     Rubc  = Rum
     Tubc  = matmul(Emm1,Tum)
     Rdbc  = matmul(Emm1, matmul(Rdm,Emm1))
     Tdbc  = matmul(Tdm,Emm1)

     call kenpsv(n,m-1,k,Ruab, Tuab, Rdab, Tdab)

     TEMP  = matmul(Ruab,Rdbc)
     Tdac  =      matmul(matmul(Tdbc,inv22(Id2-TEMP)),Tdab)
     Rdac  = Rdab+matmul(matmul(matmul(Tuab,Rdbc),inv22(Id2-TEMP)),Tdab)
     TEMP  = matmul(Rdbc,Ruab)
     Tuac  =      matmul(matmul(Tuab,inv22(Id2-TEMP)),Tubc)
     Ruac  = Rubc+matmul(matmul(matmul(Tdbc,Ruab),inv22(Id2-TEMP)),Tubc)
  endif

contains

  function inv22(a) result (ia)

    use constants
    implicit none
    real, dimension(2,2)	:: a, ia

    !local variables
    real			:: det

    det = a(1,1)*a(2,2)-a(1,2)*a(2,1)

    if( det < DET_MIN) then
       print *, "Error: inv22 called with a singular matrix"
       stop
    endif

    ia(1,1) =  a(2,2)
    ia(1,2) = -a(1,2)
    ia(2,1) = -a(2,1)
    ia(2,2) =  a(1,1)

    ia = ia/det

  end function inv22

end subroutine kenpsv

recursive subroutine kensh(n,m, k, Ruac, Tuac, Rdac, Tdac) 

  use constants
  implicit none

  integer 		:: n, m
  real			:: k
  real			:: Ruac, Tuac, Rdac, Tdac

  !local variables
  real  			:: dummy
  real  			:: Ruab, Tuab, Rdab, Tdab
  real  			:: Rubc, Tubc, Rdbc, Tdbc
  real  			:: Rum,  Tum,  Rdm,  Tdm, Emm1

  !print *, "function 'kensh(n,m)' called with n=%d and m=%d\n", n, m
  if( m < n .or. m < 1 .or. n < 0 ) then
     print *, "function 'kensh' called with n=%d and m=%d\n", n, m
     print *, "Error: 'n' should be '>=' than '0', \n"
     print *, "Error: 'm' should be '>=' than '1'  and \n"
     print *, "Error: 'n' should be '<=' than 'm' \n"
     stop
  endif

  if ((m == n) .or. (m == 1 .and. n == 0)) then
     Tdac  = 1
     Rdac  = 0
     Tuac  = 1
     Ruac  = 0
     if(n == 0) then
	!               when called with n=0; only 'Ruac' is 
        call RERSH(dummy, Ruac)
     endif
  else
     call ERTUDSH(m, k, Emm1,Rum,Tum,Rdm,Tdm)
     Rubc  = Rum
     Tubc  = Emm1*Tum
     Rdbc  = Emm1*Rdm*Emm1
     Tdbc  = Tdm*Emm1
     
     call kensh(n,m-1, k, Ruab, Tuab, Rdab, Tdab)
     Tdac  = Tdbc/(1.-Ruab*Rdbc)*Tdab
     Rdac  = Rdab+Tuab*Rdbc/(1-Ruab*Rdbc)*Tdab
     Tuac  = Tuab/(1.-Rdbc*Ruab)*Tubc
     Ruac  = Rubc+Tdbc*Ruab/(1-Rdbc*Ruab)*Tubc
  endif

end subroutine kensh

FUNCTION bessj0_s(x)
  USE nrtype; USE nrutil, ONLY : poly
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: x
  REAL(SP) :: bessj0_s
  REAL(SP) :: ax,xx,z
  REAL(DP) :: y
  REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
       0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
  REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
       0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
       -0.934945152e-7_dp/)
  REAL(DP), DIMENSION(6) :: r = (/57568490574.0_dp,-13362590354.0_dp,&
       651619640.7_dp,-11214424.18_dp,77392.33017_dp,&
       -184.9052456_dp/)
  REAL(DP), DIMENSION(6) :: s = (/57568490411.0_dp,1029532985.0_dp,&
       9494680.718_dp,59272.64853_dp,267.8532712_dp,1.0_dp/)
  if (abs(x) < 8.0) then
     y=x**2
     bessj0_s=poly(y,r)/poly(y,s)
  else
     ax=abs(x)
     z=8.0_sp/ax
     y=z**2
     xx=ax-0.785398164_sp
     bessj0_s=sqrt(0.636619772_sp/ax)*(cos(xx)*&
          poly(y,p)-z*sin(xx)*poly(y,q))
  end if
END FUNCTION bessj0_s



FUNCTION bessj0_v(x)
  USE nrtype; USE nrutil, ONLY : poly
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: x
  REAL(SP), DIMENSION(size(x)) :: bessj0_v
  REAL(SP), DIMENSION(size(x)) :: ax,xx,z
  REAL(DP), DIMENSION(size(x)) :: y
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
       0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
  REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
       0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
       -0.934945152e-7_dp/)
  REAL(DP), DIMENSION(6) :: r = (/57568490574.0_dp,-13362590354.0_dp,&
       651619640.7_dp,-11214424.18_dp,77392.33017_dp,&
       -184.9052456_dp/)
  REAL(DP), DIMENSION(6) :: s = (/57568490411.0_dp,1029532985.0_dp,&
       9494680.718_dp,59272.64853_dp,267.8532712_dp,1.0_dp/)
  mask = (abs(x) < 8.0)
  where (mask)
     y=x**2
     bessj0_v=poly(y,r,mask)/poly(y,s,mask)
  elsewhere
     ax=abs(x)
     z=8.0_sp/ax
     y=z**2
     xx=ax-0.785398164_sp
     bessj0_v=sqrt(0.636619772_sp/ax)*(cos(xx)*&
          poly(y,p,.not. mask)-z*sin(xx)*poly(y,q,.not. mask))
  end where
END FUNCTION bessj0_v
FUNCTION bessj1_s(x)
  USE nrtype; USE nrutil, ONLY : poly
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: x
  REAL(SP) :: bessj1_s
  REAL(SP) :: ax,xx,z
  REAL(DP) :: y
  REAL(DP), DIMENSION(6) :: r = (/72362614232.0_dp,&
       -7895059235.0_dp,242396853.1_dp,-2972611.439_dp,&
       15704.48260_dp,-30.16036606_dp/)
  REAL(DP), DIMENSION(6) :: s = (/144725228442.0_dp,2300535178.0_dp,&
       18583304.74_dp,99447.43394_dp,376.9991397_dp,1.0_dp/)
  REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
       -0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
  REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
       -0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
       0.105787412e-6_dp/)
  if (abs(x) < 8.0) then
     y=x**2
     bessj1_s=x*(poly(y,r)/poly(y,s))
  else
     ax=abs(x)
     z=8.0_sp/ax
     y=z**2
     xx=ax-2.356194491_sp
     bessj1_s=sqrt(0.636619772_sp/ax)*(cos(xx)*&
          poly(y,p)-z*sin(xx)*poly(y,q))*sign(1.0_sp,x)
  end if
END FUNCTION bessj1_s


FUNCTION bessj1_v(x)
  USE nrtype; USE nrutil, ONLY : poly
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: x
  REAL(SP), DIMENSION(size(x)) :: bessj1_v
  REAL(SP), DIMENSION(size(x)) :: ax,xx,z
  REAL(DP), DIMENSION(size(x)) :: y
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  REAL(DP), DIMENSION(6) :: r = (/72362614232.0_dp,&
       -7895059235.0_dp,242396853.1_dp,-2972611.439_dp,&
       15704.48260_dp,-30.16036606_dp/)
  REAL(DP), DIMENSION(6) :: s = (/144725228442.0_dp,2300535178.0_dp,&
       18583304.74_dp,99447.43394_dp,376.9991397_dp,1.0_dp/)
  REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
       -0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
  REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
       -0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
       0.105787412e-6_dp/)
  mask = (abs(x) < 8.0)
  where (mask)
     y=x**2
     bessj1_v=x*(poly(y,r,mask)/poly(y,s,mask))
  elsewhere
     ax=abs(x)
     z=8.0_sp/ax
     y=z**2
     xx=ax-2.356194491_sp
     bessj1_v=sqrt(0.636619772_sp/ax)*(cos(xx)*&
          poly(y,p,.not. mask)-z*sin(xx)*poly(y,q,.not. mask))*&
          sign(1.0_sp,x)
  end where
END FUNCTION bessj1_v
FUNCTION bessj_s(n,x)
  USE nrtype; USE nrutil, ONLY : assert
  USE nr, ONLY : bessj0,bessj1
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: n
  REAL(SP), INTENT(IN) :: x
  REAL(SP) :: bessj_s
  INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
  INTEGER(I4B) :: j,jsum,m
  REAL(SP) :: ax,bj,bjm,bjp,summ,tox
  call assert(n >= 2, 'bessj_s args')
  ax=abs(x)
  if (ax*ax <= 8.0_sp*tiny(x)) then
     bessj_s=0.0
  else if (ax > real(n,sp)) then
     tox=2.0_sp/ax
     bjm=bessj0(ax)
     bj=bessj1(ax)
     do j=1,n-1
        bjp=j*tox*bj-bjm
        bjm=bj
        bj=bjp
     end do
     bessj_s=bj
  else
     tox=2.0_sp/ax
     m=2*((n+int(sqrt(real(IACC*n,sp))))/2)
     bessj_s=0.0
     jsum=0
     summ=0.0
     bjp=0.0
     bj=1.0
     do j=m,1,-1
        bjm=j*tox*bj-bjp
        bjp=bj
        bj=bjm
        if (exponent(bj) > IEXP) then
           bj=scale(bj,-IEXP)
           bjp=scale(bjp,-IEXP)
           bessj_s=scale(bessj_s,-IEXP)
           summ=scale(summ,-IEXP)
        end if
        if (jsum /= 0) summ=summ+bj
        jsum=1-jsum
        if (j == n) bessj_s=bjp
     end do
     summ=2.0_sp*summ-bj
     bessj_s=bessj_s/summ
  end if
  if (x < 0.0 .and. mod(n,2) == 1) bessj_s=-bessj_s
END FUNCTION bessj_s

FUNCTION bessj_v(n,xx)
  USE nrtype; USE nrutil, ONLY : assert
  USE nr, ONLY : bessj0,bessj1
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: n
  REAL(SP), DIMENSION(:), INTENT(IN) :: xx
  REAL(SP), DIMENSION(size(xx)) :: bessj_v
  INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(xx)/2
  REAL(SP), DIMENSION(size(xx)) :: ax
  LOGICAL(LGT), DIMENSION(size(xx)) :: mask,mask0
  REAL(SP), DIMENSION(:), ALLOCATABLE :: x,bj,bjm,bjp,summ,tox,bessjle
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE :: renorm
  INTEGER(I4B) :: j,jsum,m,npak
  call assert(n >= 2, 'bessj_v args')
  ax=abs(xx)
  mask = (ax <= real(n,sp))
  mask0 = (ax*ax <= 8.0_sp*tiny(xx))
  bessj_v=bessjle_v(n,ax,logical(mask .and. .not.mask0, kind=lgt))
  bessj_v=merge(bessjgt_v(n,ax,.not. mask),bessj_v,.not. mask)
  where (mask0) bessj_v=0.0
  where (xx < 0.0 .and. mod(n,2) == 1) bessj_v=-bessj_v
CONTAINS
  !BL
  !BL
  FUNCTION bessjgt_v(n,xx,mask)
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(:), INTENT(IN) :: xx
    LOGICAL(LGT), DIMENSION(size(xx)), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(xx)) :: bessjgt_v
    npak=count(mask)
    if (npak == 0) RETURN
    allocate(x(npak),bj(npak),bjm(npak),bjp(npak),tox(npak))
    x=pack(xx,mask)
    tox=2.0_sp/x
    bjm=bessj0(x)
    bj=bessj1(x)
    do j=1,n-1
       bjp=j*tox*bj-bjm
       bjm=bj
       bj=bjp
    end do
    bessjgt_v=unpack(bj,mask,0.0_sp)
    deallocate(x,bj,bjm,bjp,tox)
  END FUNCTION bessjgt_v
  !BL
  FUNCTION bessjle_v(n,xx,mask)
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(:), INTENT(IN) :: xx
    LOGICAL(LGT), DIMENSION(size(xx)), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(xx)) :: bessjle_v
    npak=count(mask)
    if (npak == 0) RETURN
    allocate(x(npak),bj(npak),bjm(npak),bjp(npak),summ(npak), &
         bessjle(npak),tox(npak),renorm(npak))
    x=pack(xx,mask)
    tox=2.0_sp/x
    m=2*((n+int(sqrt(real(IACC*n,sp))))/2)
    bessjle=0.0
    jsum=0
    summ=0.0
    bjp=0.0
    bj=1.0
    do j=m,1,-1
       bjm=j*tox*bj-bjp
       bjp=bj
       bj=bjm
       renorm = (exponent(bj)>IEXP)
       bj=merge(scale(bj,-IEXP),bj,renorm)
       bjp=merge(scale(bjp,-IEXP),bjp,renorm)
       bessjle=merge(scale(bessjle,-IEXP),bessjle,renorm)
       summ=merge(scale(summ,-IEXP),summ,renorm)
       if (jsum /= 0) summ=summ+bj
       jsum=1-jsum
       if (j == n) bessjle=bjp
    end do
    summ=2.0_sp*summ-bj
    bessjle=bessjle/summ
    bessjle_v=unpack(bessjle,mask,0.0_sp)
    deallocate(x,bj,bjm,bjp,summ,bessjle,tox,renorm)
  END FUNCTION bessjle_v
END FUNCTION bessj_v

subroutine calc_BCD0E(q, w, v, qE, wE, B, C, D, E, F, nk2, kv2)

  use nrtype
  use nr
  use wavenumber
  use parameters
  use constants

  implicit none
  integer                 :: nk2
  real, dimension(0:2,nk2):: q, w, v
  real, dimension(nk2)    :: qE, wE
  real			  :: B(0:2), C(0:2), D(3:4), E, F
  real, allocatable       :: kv2(:)

  B(0) = sum(w(0,:)*kv2)
  B(1) = 0.
  B(2) = 0.
  C(0) = 0.
  C(1) = sum((q(1,:)-v(1,:))*kv2)/2.
  C(2) = 0.
  D(3) = C(1)
  D(4) = 0.

  E    = sum(wE(:)*kv2)
  F    = 0.

end subroutine calc_BCD0E
!	include 'modules.inc'
subroutine calc_coeffsE(q, w, v, qE, wE, r, B, C, D, E, F, nk3, kv3)

  use nrtype
  use nr
  use wavenumber
  use parameters
  use constants

  implicit none 
  integer                 :: nk3
  real, dimension(0:2,nk3) :: q, w, v
  real, dimension(nk3)	  :: qE, wE
  real			  :: r
  real			  :: B(0:2), C(0:2), D(3:4), E, F

  ! local variables
  real			  :: tmp1, tmp2, tmp3
  real, allocatable	  :: j0(:), j1(:), j2(:)
  real, allocatable       :: kv3(:)

  allocate(j0(nk3), j1(nk3), j2(nk3))
  j0 = bessj0( kv3*r)
  j1 = bessj1( kv3*r)
  j2 = bessj(2,kv3*r)
  B(0) =    sum(   w(0,:)        *j0*kv3)
  C(0) =   -sum(   q(0,:)        *j1*kv3)

  B(1) =    sum(   w(1,:)        *j1*kv3)
  tmp1 =    sum((  q(1,:)+v(1,:))*j1)/r
  tmp2 =    sum(kv3*q(1,:)        *j0)
  tmp3 =    sum(kv3*       v(1,:) *j0)
  C(1) =   -tmp1 + tmp2
  D(3) =    tmp1 - tmp3

  B(2) =    sum(   w(2,:)        *j2*kv3)
  tmp1 =  2*sum((  q(2,:)+v(2,:))*j2)/r
  tmp2 =    sum(kv3*q(2,:)        *j1)
  tmp3 =    sum(kv3*       v(2,:) *j1)
  C(2) =   -tmp1 + tmp2
  D(4) =    tmp1 - tmp3

  E    =    sum(  wE(:)        *j0*kv3)
  F    =   -sum(  qE(:)        *j1*kv3)
  deallocate(j0, j1, j2)

end subroutine calc_coeffsE
!	include "modules.inc"

subroutine ERTUDPSV(jn, k, EM, RU, TU, RD, TD)
  use model
  implicit none
  integer 		:: jn
  real			:: k
  real, dimension(2,2) 	:: EM, RU, TU, RD, TD

  ! local variables
  real			:: mu_1, Delta_1, d1
  real			:: mu_2, Delta_2

  real			:: ex1, Id2(2,2)
  real			:: a, b, c, d, e, f, DR

  real			:: rppp, rpsp, rspp, rssp
  real			:: rppm, rpsm, rspm, rssm
  real			:: tppp, tpsp, tspp, tssp
  real			:: tppm, tpsm, tspm, tssm
  data Id2(1,:)/1., 0./, Id2(2,:)/0., 1./

  mu_1    = mu(jn-1)
  Delta_1 = delta(jn-1)
  d1      = th(jn-1)
  mu_2    = mu(jn)
  Delta_2 = delta(jn)


  ex1 	= exp(-k*d1)
  EM   	= ex1*Id2

  a 	=                 (mu_1/mu_2+Delta_2)/(1+Delta_2)
  b 	= -2*Delta_1*k*d1*(mu_1/mu_2+Delta_2)/(1+Delta_2)
  c 	=         (mu_1*Delta_1/mu_2-Delta_2)/(1+Delta_2)
  d 	=               (mu_1*Delta_1/mu_2+1)/(1+Delta_2)
  e 	=                       (mu_1/mu_2-1)/(1+Delta_2)
  f 	=        2*Delta_1*k*d1*(mu_1/mu_2-1)/(1+Delta_2)

  DR   	=  a*d

  rppp 	= b*e/DR
  rpsp 	= -a*e/DR
  rspp 	= -(d*c-b*f)/DR
  rssp 	= -a*f/DR
  RD(1,:)	= (/rppp, rspp/)
  RD(2,:)	= (/rpsp, rssp/)

  rppm 	= 0.
  rpsm 	= d*e/DR
  rspm 	= a*c/DR
  !	rssm 	= (-a*f-b*e)/DR
  rssm 	= 0.
  RU(1,:)	= (/rppm, rspm/)
  RU(2,:)	= (/rpsm, rssm/)

  tppp 	=  a*(a*d-c*e)/DR
  !	tpsp 	=  e*(b*e+a*f)/DR
  tpsp 	=  0.
  tspp 	= -a*(b*d+c*f)/DR
  tssp 	= ((a*d-c*e)*d+(b*e+a*f)*f)/DR
  TD(1,:) = (/tppp, tspp/)
  TD(2,:) = (/tpsp, tssp/)

  tppm 	= d/DR
  tpsm 	= 0.
  tspm 	= -b/DR
  tssm 	= a/DR
  TU(1,:)	= (/tppm, tspm/)
  TU(2,:)	= (/tpsm, tssm/)
END subroutine ERTUDPSV

subroutine ERTUDSH(jn, k, EM, RU, TU, RD, TD)
  use model
  implicit none
  integer 		:: jn
  real			:: k
  real			:: EM, RU, TU, RD, TD

  ! local variables
  real			:: mu_1, d1
  real			:: mu_2
  real			:: DL

  mu_1    = mu(jn-1)
  d1      = th(jn-1)
  mu_2    = mu(jn)

  DL      =  mu_1+mu_2

  EM      = exp(-k*d1)
  RD      = (mu_1-mu_2)/DL
  TD      = 2*mu_1/DL
  RU      = -RD
  TU      = 2*mu_2/DL

END subroutine ERTUDSH

subroutine RERPSV(jn, REV, R)
  use model
  implicit none
  integer			:: jn
  real, dimension(2,2)	:: REV, R

  ! local variables
  real			:: idelta

  idelta	 = 1./delta(jn)
  REV(1,:) = (/ idelta, -1. /)
  REV(2,:) = (/ idelta,  1. /)
  REV	 = REV*(1.+delta(jn))

  R(1,:)	 = (/ 0.,  -delta(jn) /)
  R(2,:)	 = (/-idelta,	 0.  /)

END subroutine RERPSV

subroutine RERSH(REV, R)
  implicit none
  real			:: REV, R

  REV	 = 2.
  R	 = 1.

END subroutine RERSH
! program f2
!   integer i,j


!   open(unit=7, file='sequen', form='unformatted')
!   do i=1,20
!      read(7, end=10) j
!      write(*,*) j
!   end do
!   close(7)

! 10 continue

!   open(unit=8, file='direct', form='unformatted', recl=4, access='direct')
!   do i=1,20
!      read(8,rec=i,end=20) j
!      write(*,*) j
!   end do
!   close(8)

! 20 continue
! end program f2

