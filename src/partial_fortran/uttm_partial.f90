MODULE allFuncs

IMPLICIT NONE

CONTAINS

FUNCTION cross (arr1, arr2)

REAL(8), DIMENSION(3) :: cross
REAL(8), DIMENSION(3), INTENT(IN) :: arr1, arr2

cross(1) = arr1(2) * arr2(3) - arr1(3) * arr2(2)
cross(2) = arr1(3) * arr2(1) - arr1(1) * arr2(3)
cross(3) = arr1(1) * arr2(2) - arr1(2) * arr2(1)

END FUNCTION

END MODULE allFuncs

PROGRAM uttm_partial

! This program is intended to run as backend to accelerate the python code
! It runs the lagrangian part of the code responsible for computing the 
! thin tube velocity

! Import modules
USE allFuncs

! Declarations
IMPLICIT NONE

! Number of points
INTEGER :: NP
! Others
INTEGER :: i, i_j, flag_overflow, ll, X, Y, Z, i_b, n_b, pp, j, image
REAL(8) :: delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, PI, ln_o_ln_coeff, gamma_tmp1
REAL(8) :: one_o_norm3, norm3, kernel_tanh_coeff1, kernel_tanh_coeff2

! Create arrays to store data
REAL(8), DIMENSION(3) :: vec_auto, left_right, vec_tmp, vec_tmp1, point_c, vec_auto1, vec_auto2, vec_tmp2, left_right_b
REAL(8), DIMENSION(:), ALLOCATABLE :: Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA
REAL(8), DIMENSION(:), ALLOCATABLE :: ux_tmp, uy_tmp, uz_tmp, uux_s, uuy_s, uuz_s, sigma_u

! Denote X, Y, Z as 1, 2, 3
X = 1
Y = 2
Z = 3

! Read input from file
OPEN(1, FILE = '../../fio.dat', STATUS = 'old')
READ(1, *) NP
READ(1, *) delta_ttm_tmp
READ(1, *) sigma_1
READ(1, *) sigma_2
READ(1, *) gamma_tmp
READ(1, *) ds
READ(1, *) n_b
READ(1, *) image

!PRINT *, 'Num points:', NP
!PRINT *, 'delta ttm:', delta_ttm_tmp
!PRINT *, 'sigma 1:', sigma_1
!PRINT *, 'sigma 2:', sigma_2
!PRINT *, 'Gamma:', gamma_tmp
!PRINT *, 'ds:', ds
!PRINT *, 'nImages:', n_b
!PRINT *, 'image on:', image

! Allocate memory
ALLOCATE(Ux(NP+2))
ALLOCATE(Uy(NP+2))
ALLOCATE(Uz(NP+2))
ALLOCATE(Ux_s(NP+2))
ALLOCATE(Uy_s(NP+2))
ALLOCATE(Uz_s(NP+2))
ALLOCATE(SIGMA(NP+2))

ALLOCATE(ux_tmp(NP+2))
ALLOCATE(uy_tmp(NP+2))
ALLOCATE(uz_tmp(NP+2))
ALLOCATE(uux_s(NP+2))
ALLOCATE(uuy_s(NP+2))
ALLOCATE(uuz_s(NP+2))
ALLOCATE(sigma_u(NP+2))

!PRINT *, 'Memory allocated..'

! Read other data points
DO i = 1, NP+2
READ(1, *) Ux(i), Uy(i), Uz(i), Ux_s(i), Uy_s(i), Uz_s(i), SIGMA(i)
!PRINT *, i, Ux(i), Uy(i), Uz(i), Ux_s(i), Uy_s(i), Uz_s(i), SIGMA(i)
END DO

CLOSE(1)

! Define
PI=4.D0*DATAN(1.D0)

! Compute uttm for all points and write output
OPEN(2, FILE = 'fout.dat', STATUS = 'replace')

DO i = 1, NP+1
!DO i = 1, 3

! Something is not getting reset correctly

vec_auto(X) = 0.
vec_auto(Y) = 0.
vec_auto(Z) = 0.

vec_auto1(X) = 0.
vec_auto1(Y) = 0.
vec_auto1(Z) = 0.

vec_auto2(X) = 0.
vec_auto2(Y) = 0.
vec_auto2(Z) = 0.

ux_tmp(:) = 0.
uy_tmp(:) = 0.
uz_tmp(:) = 0.

uux_s(:) = 0.
uuy_s(:) = 0.
uuz_s(:) = 0.

sigma_u(:) = 0.

left_right(X) = 0.
left_right(Y) = 0.
left_right(Z) = 0.

vec_tmp(X) = 0.
vec_tmp(Y) = 0.
vec_tmp(Z) = 0.

vec_tmp1(X) = 0.
vec_tmp1(Y) = 0.
vec_tmp1(Z) = 0.

vec_tmp2(X) = 0.
vec_tmp2(Y) = 0.
vec_tmp2(Z) = 0.

point_c(X) = 0.
point_c(Y) = 0.
point_c(Z) = 0.

point_c(X) = Ux(i)
point_c(Y) = Uy(i)
point_c(Z) = Uz(i)

ln_o_ln_coeff = 0.
gamma_tmp1 = 0.
one_o_norm3 = 0.
norm3 = 0.
kernel_tanh_coeff1 = 0.
kernel_tanh_coeff2 = 0.
left_right_b(X) = 0.
left_right_b(Y) = 0.
left_right_b(Z) = 0.

left_right(X) = 0.
left_right(Y) = 0.
left_right(Z) = 0.

ll = 0

!PRINT *, vec_auto ! Checked
!PRINT *, vec_auto1, vec_auto2 ! Checked

!DO j = 1, NP+2
!PRINT*, ux_tmp(j), uy_tmp(j), uz_tmp(j), uux_s(j), uuy_s(j), uuz_s(j), sigma_u(j)
!END DO ! Checked

!PRINT *, vec_tmp, vec_tmp1, vec_tmp2, point_c ! Checked

!PRINT *, point_c(X), point_c(Y), point_c(Z)
!! Checked

left_right(X) = Ux(NP) - Ux(1)
left_right(Y) = Uy(NP) - Uy(1)
left_right(Z) = Uz(NP) - Uz(1)
!PRINT *, left_right(X), left_right(Y), left_right(Z)
!!PRINT *, Ux(NP+1), Ux(2), Uz(NP+1), Uz(2)
!!PRINT *, ((NP+2)/2)+1
!! Checked
flag_overflow = 0

DO i_j = 1, ((NP+1)/2)

ll = i_j + i - 1

IF (ll > NP) THEN

ll = ll - (NP-1)
flag_overflow = 1

END IF

IF (flag_overflow == 0) THEN

ux_tmp(((NP+1)/2) + i_j - 1) = Ux(ll)
uy_tmp(((NP+1)/2) + i_j - 1) = Uy(ll)
uz_tmp(((NP+1)/2) + i_j - 1) = Uz(ll)

ELSE

ux_tmp(((NP+1)/2) + i_j - 1) = Ux(ll) + left_right(X)
uy_tmp(((NP+1)/2) + i_j - 1) = Uy(ll) + left_right(Y)
uz_tmp(((NP+1)/2) + i_j - 1) = Uz(ll) + left_right(Z)

END IF

uux_s(((NP+1)/2) + i_j - 1) = Ux_s(ll)
uuy_s(((NP+1)/2) + i_j - 1) = Uy_s(ll)
uuz_s(((NP+1)/2) + i_j - 1) = Uz_s(ll)
sigma_u(((NP+1)/2) + i_j - 1) = SIGMA(ll)

!PRINT *, i_j, ux_tmp(i_j), uy_tmp(i_j), uz_tmp(i_j)
!PRINT *, i_j, uux_s(i_j), uuy_s(i_j), uuz_s(i_j)
!PRINT *, i_j, sigma_u(i_j)
!PRINT *, ll
! Checked. Checked again. Seems all right. Checked uux_s as well. Seems alright.
! sigma_u is correct as well.

END DO

flag_overflow = 0

!!PRINT *, (((NP+1)/2)-1)+1

DO i_j = 1, (((NP+1)/2)-1)

ll = i - i_j

IF (ll < 1) THEN

ll = (NP-1) + ll
flag_overflow = 1

END IF

IF (flag_overflow == 0) THEN

ux_tmp(((NP+1)/2) - i_j) = Ux(ll)
uy_tmp(((NP+1)/2) - i_j) = Uy(ll)
uz_tmp(((NP+1)/2) - i_j) = Uz(ll)

ELSE

ux_tmp(((NP+1)/2) - i_j) = Ux(ll) - left_right(X)
uy_tmp(((NP+1)/2) - i_j) = Uy(ll) - left_right(Y)
uz_tmp(((NP+1)/2) - i_j) = Uz(ll) - left_right(Z)

END IF

uux_s(((NP+1)/2) - i_j) = Ux_s(ll)
uuy_s(((NP+1)/2) - i_j) = Uy_s(ll)
uuz_s(((NP+1)/2) - i_j) = Uz_s(ll)
sigma_u(((NP+1)/2) - i_j) = SIGMA(ll)

!PRINT *, i_j, ux_tmp(i_j), uy_tmp(i_j), uz_tmp(i_j)
!PRINT *, i_j, uux_s(i_j), uuy_s(i_j), uuz_s(i_j)
!PRINT *, i_j, sigma_u(i_j)
!PRINT *, ll
! Checked. Checked again. Seems all right. Checked uux_s also. Seems alright.
! sigma_u is correct as well.

END DO

ln_o_ln_coeff = LOG(sigma_1/delta_ttm_tmp)/LOG(sigma_2/sigma_1)
gamma_tmp1 = gamma_tmp / (4 * PI)
left_right(X) = ux_tmp(NP) - ux_tmp(1)
left_right(Y) = uy_tmp(NP) - uy_tmp(1)
left_right(Z) = uz_tmp(NP) - uz_tmp(1)

!PRINT *, left_right(X), left_right(Y), left_right(Z)
!PRINT *, vec_auto1(X), vec_auto1(Y), vec_auto1(Z)
! Checked

! Central part of the domain

!PRINT *, (NP+1)/2

DO i_j = 1, NP

IF (i_j /= ((NP+1)/2)) THEN

vec_tmp(X) = uux_s(i_j)
vec_tmp(Y) = uuy_s(i_j)
vec_tmp(Z) = uuz_s(i_j)

vec_tmp1(X) = ux_tmp(i_j)
vec_tmp1(Y) = uy_tmp(i_j)
vec_tmp1(Z) = uz_tmp(i_j)

vec_tmp1 = point_c - vec_tmp1
one_o_norm3 = SQRT(vec_tmp1(X)**2 + vec_tmp1(Y)**2 + vec_tmp1(Z)**2)
norm3 = one_o_norm3**3
one_o_norm3 = 1/norm3
kernel_tanh_coeff1 = TANH(norm3/(sigma_1**3))
kernel_tanh_coeff2 = TANH(norm3/(sigma_2**3))

vec_tmp1 = one_o_norm3 * vec_tmp1

vec_tmp2 = cross(vec_tmp, vec_tmp1)
vec_auto1 = kernel_tanh_coeff1 * vec_tmp2
vec_auto2 = kernel_tanh_coeff2 * vec_tmp2

vec_tmp2 = vec_auto1 + (vec_auto1 - vec_auto2) * ln_o_ln_coeff

vec_auto(X) = vec_auto(X) + gamma_tmp1 * ds * vec_tmp2(X)
vec_auto(Y) = vec_auto(Y) + gamma_tmp1 * ds * vec_tmp2(Y)
vec_auto(Z) = vec_auto(Z) + gamma_tmp1 * ds * vec_tmp2(Z)

END IF

!!PRINT *, i_j, gamma_tmp1, ds ! This is ok
!!PRINT *, i_j, vec_tmp2(X), vec_tmp2(Y), vec_tmp2(Z)
!!PRINT *, i_j, ln_o_ln_coeff ! This is ok
!!PRINT *, i_j, vec_auto1(X), vec_auto1(Y), vec_auto1(Z)
!!PRINT *, i_j, vec_auto2(X), vec_auto2(Y), vec_auto2(Z)
!!PRINT *, i_j, vec_tmp(X), vec_tmp(Y), vec_tmp(Z) ! This is ok
!!PRINT *, i_j, one_o_norm3 ! This is ok
!!PRINT *, i_j, point_c(X) - vec_tmp1(X), point_c(Y) - vec_tmp1(Y), point_c(Z) - vec_tmp1(Z) ! This is ok
!!PRINT *, i_j, kernel_tanh_coeff1, kernel_tanh_coeff2

!PRINT *, vec_auto(X), vec_auto(Y), vec_auto(Z)
!!PRINT *, vec_auto1(X), vec_auto1(Y), vec_auto1(Z)
!!PRINT *, vec_auto2(X), vec_auto2(Y), vec_auto2(Z)
!!PRINT *, i_j, vec_tmp(X), vec_tmp(Y), vec_tmp(Z)
!!PRINT *, i_j, vec_tmp1(X), vec_tmp1(Y), vec_tmp1(Z)
!!PRINT *, i_j
!!PRINT *, kernel_tanh_coeff1, kernel_tanh_coeff2
!!PRINT *, ln_o_ln_coeff, gamma_tmp, ds
!!PRINT *, i_j, TANH(norm3/(sigma_2**3))
!!PRINT *, i_j, cross(vec_tmp, vec_tmp1)
!!PRINT *, cross(vec_tmp, vec_tmp1)

!! Checked

END DO

DO i_b = 1, n_b - 1

left_right_b = i_b * left_right

! Left image

DO i_j = 1, NP-1

vec_tmp(X) = uux_s(i_j)
vec_tmp(Y) = uuy_s(i_j)
vec_tmp(Z) = uuz_s(i_j)

vec_tmp1(X) = ux_tmp(i_j) - left_right_b(X)
vec_tmp1(Y) = uy_tmp(i_j) - left_right_b(Y)
vec_tmp1(Z) = uz_tmp(i_j) - left_right_b(Z)

vec_tmp1 = point_c - vec_tmp1
one_o_norm3 = SQRT(vec_tmp1(X)**2 + vec_tmp1(Y)**2 + vec_tmp1(Z)**2)
one_o_norm3 = 1/one_o_norm3**3

vec_tmp1 = one_o_norm3 * vec_tmp1
vec_tmp2 = cross(vec_tmp, vec_tmp1)

vec_auto(X) = vec_auto(X) + gamma_tmp1 * ds * vec_tmp2(X)
vec_auto(Y) = vec_auto(Y) + gamma_tmp1 * ds * vec_tmp2(Y)
vec_auto(Z) = vec_auto(Z) + gamma_tmp1 * ds * vec_tmp2(Z)

END DO

!PRINT *, left_right_b
!PRINT *, vec_auto(X), vec_auto(Y), vec_auto(Z)
!PRINT *, vec_tmp(X), vec_tmp(Y), vec_tmp(Z)
!PRINT *, vec_tmp1(X), vec_tmp1(Y), vec_tmp1(Z)
!PRINT *, vec_tmp2(X), vec_tmp2(Y), vec_tmp2(Z)

! Right image

DO i_j = 2, NP

vec_tmp(X) = uux_s(i_j)
vec_tmp(Y) = uuy_s(i_j)
vec_tmp(Z) = uuz_s(i_j)

vec_tmp1(X) = ux_tmp(i_j) + left_right_b(X)
vec_tmp1(Y) = uy_tmp(i_j) + left_right_b(Y)
vec_tmp1(Z) = uz_tmp(i_j) + left_right_b(Z)

vec_tmp1 = point_c - vec_tmp1
one_o_norm3 = SQRT(vec_tmp1(X)**2 + vec_tmp1(Y)**2 + vec_tmp1(Z)**2)
one_o_norm3 = 1/one_o_norm3**3

vec_tmp1 = one_o_norm3 * vec_tmp1
vec_tmp2 = cross(vec_tmp, vec_tmp1)

vec_auto(X) = vec_auto(X) + gamma_tmp1 * ds * vec_tmp2(X)
vec_auto(Y) = vec_auto(Y) + gamma_tmp1 * ds * vec_tmp2(Y)
vec_auto(Z) = vec_auto(Z) + gamma_tmp1 * ds * vec_tmp2(Z)

!PRINT *, vec_auto(X), vec_auto(Y), vec_auto(Z)

END DO

!PRINT *, vec_auto(X), vec_auto(Y), vec_auto(Z)
!PRINT *, left_right_b(X), left_right_b(Y), left_right_b(Z)

END DO

!PRINT *, vec_auto(X), vec_auto(Y), vec_auto(Z)
!! Minor difference is there between python and fortran. Needs to be resolved.
!! Difference has been resolved. Not completely. Seems that this works for one iteration only. 
!! Difference resolved completely


! This portion of the code is when image vortices beyond the wall are necessary

IF (image == 1) THEN

DO i_j = 1, NP

IF (i_j /= ((NP+1)/2)) THEN

vec_tmp(X) = uux_s(i_j)
vec_tmp(Y) = uuy_s(i_j)
vec_tmp(Z) = uuz_s(i_j)

vec_tmp1(X) = ux_tmp(i_j)
vec_tmp1(Y) = -uy_tmp(i_j)
vec_tmp1(Z) = uz_tmp(i_j)

vec_tmp1 = point_c - vec_tmp1
one_o_norm3 = SQRT(vec_tmp1(X)**2 + vec_tmp1(Y)**2 + vec_tmp1(Z)**2)
norm3 = one_o_norm3**3
one_o_norm3 = 1/norm3
!kernel_tanh_coeff1 = TANH(norm3/(sigma_1**3))
!kernel_tanh_coeff2 = TANH(norm3/(sigma_2**3))

vec_tmp1 = one_o_norm3 * vec_tmp1

vec_tmp2 = cross(vec_tmp, vec_tmp1)
!vec_auto1 = kernel_tanh_coeff1 * vec_tmp2
!vec_auto2 = kernel_tanh_coeff2 * vec_tmp2

!vec_tmp2 = vec_auto1 + (vec_auto1 - vec_auto2) * ln_o_ln_coeff

vec_auto(X) = vec_auto(X) + gamma_tmp1 * ds * vec_tmp2(X)
vec_auto(Y) = vec_auto(Y) - gamma_tmp1 * ds * vec_tmp2(Y)
vec_auto(Z) = vec_auto(Z) + gamma_tmp1 * ds * vec_tmp2(Z)

END IF

!!PRINT *, i_j, gamma_tmp1, ds ! This is ok
!!PRINT *, i_j, vec_tmp2(X), vec_tmp2(Y), vec_tmp2(Z)
!!PRINT *, i_j, ln_o_ln_coeff ! This is ok
!!PRINT *, i_j, vec_auto1(X), vec_auto1(Y), vec_auto1(Z)
!!PRINT *, i_j, vec_auto2(X), vec_auto2(Y), vec_auto2(Z)
!!PRINT *, i_j, vec_tmp(X), vec_tmp(Y), vec_tmp(Z) ! This is ok
!!PRINT *, i_j, one_o_norm3 ! This is ok
!!PRINT *, i_j, point_c(X) - vec_tmp1(X), point_c(Y) - vec_tmp1(Y), point_c(Z) - vec_tmp1(Z) ! This is ok
!!PRINT *, i_j, kernel_tanh_coeff1, kernel_tanh_coeff2

!PRINT *, vec_auto(X), vec_auto(Y), vec_auto(Z)
!!PRINT *, vec_auto1(X), vec_auto1(Y), vec_auto1(Z)
!!PRINT *, vec_auto2(X), vec_auto2(Y), vec_auto2(Z)
!!PRINT *, i_j, vec_tmp(X), vec_tmp(Y), vec_tmp(Z)
!!PRINT *, i_j, vec_tmp1(X), vec_tmp1(Y), vec_tmp1(Z)
!!PRINT *, i_j
!!PRINT *, kernel_tanh_coeff1, kernel_tanh_coeff2
!!PRINT *, ln_o_ln_coeff, gamma_tmp, ds
!!PRINT *, i_j, TANH(norm3/(sigma_2**3))
!!PRINT *, i_j, cross(vec_tmp, vec_tmp1)
!!PRINT *, cross(vec_tmp, vec_tmp1)

!! Checked

END DO

DO i_b = 1, n_b - 1

left_right_b = i_b * left_right

! Left image

DO i_j = 1, NP-1

vec_tmp(X) = uux_s(i_j)
vec_tmp(Y) = uuy_s(i_j)
vec_tmp(Z) = uuz_s(i_j)

vec_tmp1(X) = ux_tmp(i_j) - left_right_b(X)
vec_tmp1(Y) = -uy_tmp(i_j) - left_right_b(Y)
vec_tmp1(Z) = uz_tmp(i_j) - left_right_b(Z)

vec_tmp1 = point_c - vec_tmp1
one_o_norm3 = SQRT(vec_tmp1(X)**2 + vec_tmp1(Y)**2 + vec_tmp1(Z)**2)
one_o_norm3 = 1/one_o_norm3**3

vec_tmp1 = one_o_norm3 * vec_tmp1
vec_tmp2 = cross(vec_tmp, vec_tmp1)

vec_auto(X) = vec_auto(X) + gamma_tmp1 * ds * vec_tmp2(X)
vec_auto(Y) = vec_auto(Y) - gamma_tmp1 * ds * vec_tmp2(Y)
vec_auto(Z) = vec_auto(Z) + gamma_tmp1 * ds * vec_tmp2(Z)

END DO

!PRINT *, left_right_b
!PRINT *, vec_auto(X), vec_auto(Y), vec_auto(Z)
!PRINT *, vec_tmp(X), vec_tmp(Y), vec_tmp(Z)
!PRINT *, vec_tmp1(X), vec_tmp1(Y), vec_tmp1(Z)
!PRINT *, vec_tmp2(X), vec_tmp2(Y), vec_tmp2(Z)

! Right image

DO i_j = 2, NP

vec_tmp(X) = uux_s(i_j)
vec_tmp(Y) = uuy_s(i_j)
vec_tmp(Z) = uuz_s(i_j)

vec_tmp1(X) = ux_tmp(i_j) + left_right_b(X)
vec_tmp1(Y) = -uy_tmp(i_j) + left_right_b(Y)
vec_tmp1(Z) = uz_tmp(i_j) + left_right_b(Z)

vec_tmp1 = point_c - vec_tmp1
one_o_norm3 = SQRT(vec_tmp1(X)**2 + vec_tmp1(Y)**2 + vec_tmp1(Z)**2)
one_o_norm3 = 1/one_o_norm3**3

vec_tmp1 = one_o_norm3 * vec_tmp1
vec_tmp2 = cross(vec_tmp, vec_tmp1)

vec_auto(X) = vec_auto(X) + gamma_tmp1 * ds * vec_tmp2(X)
vec_auto(Y) = vec_auto(Y) - gamma_tmp1 * ds * vec_tmp2(Y)
vec_auto(Z) = vec_auto(Z) + gamma_tmp1 * ds * vec_tmp2(Z)

!PRINT *, vec_auto(X), vec_auto(Y), vec_auto(Z)

END DO

!PRINT *, vec_auto(X), vec_auto(Y), vec_auto(Z)
!PRINT *, left_right_b(X), left_right_b(Y), left_right_b(Z)

END DO

END IF

WRITE(2, *) vec_auto(X), vec_auto(Y), vec_auto(Z)

END DO

CLOSE(2)

!PRINT *, 'Done..'

END PROGRAM uttm_partial
