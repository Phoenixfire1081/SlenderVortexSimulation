MODULE allFuncs

IMPLICIT NONE

CONTAINS

FUNCTION cross (arr1, arr2)

REAL(8), DIMENSION(3) :: cross
REAL(8), DIMENSION(3), INTENT(IN) :: arr1, arr2

cross(1) = arr1(2) * arr2(3) - arr1(3) * arr2(2)
cross(2) = arr1(3) * arr2(1) - arr1(1) * arr2(3)
cross(3) = arr1(1) * arr2(2) - arr1(2) * arr2(1)

END FUNCTION cross

FUNCTION uttm (additional, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA, NP, n_b, &
image, delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, ux_tmp, &
uy_tmp, uz_tmp, uux_s, uuy_s, uuz_s, sigma_u)

! Define stuff
INTEGER :: i, i_j, flag_overflow, ll, X, Y, Z, i_b, pp, j
REAL(8) :: one_o_norm3, norm3, kernel_tanh_coeff1, kernel_tanh_coeff2, PI, ln_o_ln_coeff, gamma_tmp1
! Create arrays to store data
REAL(8), DIMENSION(3) :: vec_auto, left_right, vec_tmp, vec_tmp1, point_c, vec_auto1, vec_auto2, vec_tmp2, left_right_b

REAL(8), DIMENSION(NP+2, 3) :: uttm
INTEGER, INTENT(IN) :: NP, n_b, image
REAL(8), INTENT(IN) :: delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds
REAL(8), DIMENSION(NP+2), INTENT(IN) :: Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA
REAL(8), DIMENSION(NP+2, 3), INTENT(IN) :: additional
REAL(8), DIMENSION(NP+2), INTENT(INOUT) :: ux_tmp, uy_tmp, uz_tmp, uux_s, uuy_s, uuz_s, sigma_u

! Define
PI=4.D0*DATAN(1.D0)

! Denote X, Y, Z as 1, 2, 3
X = 1
Y = 2
Z = 3

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

point_c(X) = Ux(i) + additional(i, X)
point_c(Y) = Uy(i) + additional(i, Y)
point_c(Z) = Uz(i) + additional(i, Z)

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
!PRINT *, 
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

vec_tmp1 = one_o_norm3 * vec_tmp1
vec_tmp2 = cross(vec_tmp, vec_tmp1)

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

uttm(i, X) = vec_auto(X)
uttm(i, Y) = vec_auto(Y)
uttm(i, Z) = vec_auto(Z)

END DO

END FUNCTION uttm

FUNCTION CV(delta_bar)

REAL(8) :: CV
REAL(8), INTENT(IN) :: delta_bar

CV = ((1.0 + 0.577 - LOG(2.0))/2.0) - LOG(delta_bar)
!CV = ((0.5) + (0.25) - LOG(delta_bar))
!return ((1/2) + (1/4) - np.log(delta_bar))

END FUNCTION CV

END MODULE allFuncs

PROGRAM uttm_full

! This program is intended to run as backend to accelerate the python code
! It runs the lagrangian part of the code responsible for computing the 
! thin tube velocity

! Import modules
USE allFuncs

! Declarations
IMPLICIT NONE

! Define stuff
INTEGER :: i, i_j, flag_overflow, ll, X, Y, Z, i_b, pp, j

!! Number of points
INTEGER :: NP
!! Others
INTEGER :: n_b, image, flowCondition, interpArraySize, nsteps, step, writeEveryNSteps
REAL(8) :: delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, Ub, maxV, yShear, dt
REAL(8) :: S_0, S_tmp, S, S0_o_S, int_S_o_S0, delta_bar, delta_0_bar_param
REAL(8) :: nu_bar_param, c_v, beta, PI, gamma_param, epsilonVal, delta_ttm
!REAL(8) :: start, finish, elapsed
REAL(8) :: c_ttm, sigma_max, sigma3_1, sigma3_2, one_o_2ds, one_o_ds2, sigma_tmp, one_o_sigma3
REAL(8), DIMENSION(:), ALLOCATABLE :: yPlus, uVelBL, vVelBL, wVelBL
REAL(8), DIMENSION(:), ALLOCATABLE :: Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA, Ux_ss, Uy_ss, Uz_ss
REAL(8), DIMENSION(:), ALLOCATABLE :: ux_tmp, uy_tmp, uz_tmp, uux_s, uuy_s, uuz_s, sigma_u
REAL(8), DIMENSION(:), ALLOCATABLE :: Ux_0, Uy_0, Uz_0
REAL(8), DIMENSION(:, :), ALLOCATABLE :: VxVyVz, vec_local, additional, k1, k2, k3, k4, k5, k6, vec_add
REAL(8), DIMENSION(:, :), ALLOCATABLE :: VxVyVz_0, VxVyVz_m, VxVyVz_m1, VxVyVz_m2
REAL(8), DIMENSION(3) :: left_rightBC, point_m, point_c, point_p, vec_tmp, vec_tmp1, vec_tmp2

! Start timer

!CALL CPU_TIME(start)

! Define
PI=4.D0*DATAN(1.D0)

X = 1
Y = 2
Z = 3

! Read input from file
OPEN(1, FILE = '../../fio.dat', STATUS = 'old')
READ(1, *) writeEveryNSteps
READ(1, *) S_0
READ(1, *) delta_0_bar_param
READ(1, *) nu_bar_param
READ(1, *) gamma_param
READ(1, *) epsilonVal
READ(1, *) c_ttm
READ(1, *) nsteps
READ(1, *) NP
READ(1, *) delta_ttm_tmp
READ(1, *) sigma_1
READ(1, *) sigma_2
READ(1, *) gamma_tmp
READ(1, *) ds
READ(1, *) n_b
READ(1, *) dt
READ(1, *) image
READ(1, *) flowCondition

one_o_2ds = 1 / (ds * 2)
one_o_ds2 = 1 / (ds * ds)

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
ALLOCATE(Ux_ss(NP+2))
ALLOCATE(Uy_ss(NP+2))
ALLOCATE(Uz_ss(NP+2))
ALLOCATE(SIGMA(NP+2))

ALLOCATE(ux_tmp(NP+2))
ALLOCATE(uy_tmp(NP+2))
ALLOCATE(uz_tmp(NP+2))
ALLOCATE(uux_s(NP+2))
ALLOCATE(uuy_s(NP+2))
ALLOCATE(uuz_s(NP+2))
ALLOCATE(sigma_u(NP+2))

ALLOCATE(Ux_0(NP+2))
ALLOCATE(Uy_0(NP+2))
ALLOCATE(Uz_0(NP+2))

ALLOCATE(VxVyVz(NP+2, 3))
ALLOCATE(vec_local(NP+2, 3))
ALLOCATE(additional(NP+2, 3))
ALLOCATE(k1(NP+2, 3))
ALLOCATE(k2(NP+2, 3))
ALLOCATE(k3(NP+2, 3))
ALLOCATE(k4(NP+2, 3))
ALLOCATE(k5(NP+2, 3))
ALLOCATE(k6(NP+2, 3))
ALLOCATE(vec_add(NP+2, 3))
ALLOCATE(VxVyVz_0(NP+2, 3))
ALLOCATE(VxVyVz_m(NP+2, 3))
ALLOCATE(VxVyVz_m1(NP+2, 3))
ALLOCATE(VxVyVz_m2(NP+2, 3))

! Set others to 0 as well (outside time loop)
int_S_o_S0 = 0.
VxVyVz_0(:, :) = 0.
VxVyVz_m(:, :) = 0.
VxVyVz_m1(:, :) = 0.
VxVyVz_m2(:, :) = 0.

OPEN(2, FILE = 'fout.dat', STATUS = 'replace')
CLOSE(2)

! Start main time loop
DO step = 1, nsteps

! Set all additional values to 0 first
additional(:, :) = 0.
VxVyVz(:, :) = 0.
Ux_0(:) = 0.
Uy_0(:) = 0.
Uz_0(:) = 0.

!PRINT *, 'Memory allocated..'

! Impose BC after step 1

IF (step > 1) THEN

!left_rightBC(X) = Ux(NP) - Ux(1)
!left_rightBC(Y) = Uy(NP) - Uy(1)
!left_rightBC(Z) = Uz(NP) - Uz(1)

!Ux(1) = Ux(NP-1) - left_rightBC(X)
!Uy(1) = Uy(NP-1) - left_rightBC(Y)
!Uz(1) = Uz(NP-1) - left_rightBC(Z)

!Ux(NP+1) = Ux(2) + left_rightBC(X)
!Uy(NP+1) = Uy(2) + left_rightBC(Y)
!Uz(NP+1) = Uz(2) + left_rightBC(Z)

!Ux_s(1) = Ux_s(NP-1)
!Uy_s(1) = Uy_s(NP-1)
!Uz_s(1) = Uz_s(NP-1)
!SIGMA(1) = SIGMA(NP-1)

!Ux_s(NP+1) = Ux_s(2)
!Uy_s(NP+1) = Uy_s(2)
!Uz_s(NP+1) = Uz_s(2)
!SIGMA(NP+1) = SIGMA(2)

!Ux_ss(1) = Ux_ss(NP-1)
!Uy_ss(1) = Uy_ss(NP-1)
!Uz_ss(1) = Uz_ss(NP-1)

!!PRINT *, Ux_ss(1), Uy_ss(1), Uz_ss(1), SIGMA(1)

!Ux_ss(NP+1) = Ux_ss(2)
!Uy_ss(NP+1) = Uy_ss(2)
!Uz_ss(NP+1) = Uz_ss(2)

!Ux(NP) = Ux(1) + left_rightBC(X)
!Uy(NP) = Uy(1) + left_rightBC(Y)
!Uz(NP) = Uz(1) + left_rightBC(Z)

!!PRINT *, Ux(NP), Uy(NP), Uz(NP)

!Ux_s(NP) = Ux_s(1)
!Uy_s(NP) = Uy_s(1)
!Uz_s(NP) = Uz_s(1)
!SIGMA(NP) = SIGMA(1)

!Ux_ss(NP) = Ux_ss(1)
!Uy_ss(NP) = Uy_ss(1)
!Uz_ss(NP) = Uz_ss(1)

! Compute beta

beta = gamma_param / (4*PI)

! M1 method of Knio Klein

delta_ttm = epsilonVal
delta_ttm = delta_ttm * EXP(c_ttm - c_v + 1)

delta_ttm_tmp = 0.
delta_ttm_tmp = delta_ttm_tmp + delta_ttm 

sigma_max = MAXVAL(SIGMA)
sigma_max = ds * sigma_max

sigma_1 = 3 * sigma_max
sigma_2 = 2 * sigma_1

sigma3_1 = sigma_1**3
sigma3_2 = sigma_2**3

!PRINT *, 'dttm', epsilonVal
!PRINT *, 'c_v', c_v
!PRINT *, gamma_param

! Compute derivatives with finite difference

DO i = 2, NP+1

point_m(:) = 0.
point_c(:) = 0.
point_p(:) = 0.

point_m(X) = point_m(X) + Ux(i-1)
point_m(Y) = point_m(Y) + Uy(i-1)
point_m(Z) = point_m(Z) + Uz(i-1)

point_c(X) = point_c(X) + Ux(i)
point_c(Y) = point_c(Y) + Uy(i)
point_c(Z) = point_c(Z) + Uz(i)

point_p(X) = point_p(X) + Ux(i+1)
point_p(Y) = point_p(Y) + Uy(i+1)
point_p(Z) = point_p(Z) + Uz(i+1)

vec_tmp = point_p - point_m

Ux_s(i) = vec_tmp(X) * one_o_2ds
Uy_s(i) = vec_tmp(Y) * one_o_2ds
Uz_s(i) = vec_tmp(Z) * one_o_2ds

SIGMA(i) = NORM2(vec_tmp) * one_o_2ds
sigma_tmp = SIGMA(i)

one_o_sigma3 = 1 / (sigma_tmp**3)

vec_tmp1 = point_p + point_m
vec_tmp1 = vec_tmp1 - point_c
vec_tmp1 = vec_tmp1 - point_c

vec_tmp2 = CROSS(vec_tmp, vec_tmp1)

vec_tmp2(X) = vec_tmp2(X) * one_o_2ds * one_o_ds2
Ux_ss(i) = vec_tmp2(X) * one_o_sigma3

vec_tmp2(Y) = vec_tmp2(Y) * one_o_2ds * one_o_ds2
Uy_ss(i) = vec_tmp2(Y) * one_o_sigma3

vec_tmp2(Z) = vec_tmp2(Z) * one_o_2ds * one_o_ds2
Uz_ss(i) = vec_tmp2(Z) * one_o_sigma3

END DO

! Impose BC again

!left_rightBC(X) = Ux(NP) - Ux(1)
!left_rightBC(Y) = Uy(NP) - Uy(1)
!left_rightBC(Z) = Uz(NP) - Uz(1)

!Ux(1) = Ux(NP-1) - left_rightBC(X)
!Uy(1) = Uy(NP-1) - left_rightBC(Y)
!Uz(1) = Uz(NP-1) - left_rightBC(Z)

!Ux(NP+1) = Ux(2) + left_rightBC(X)
!Uy(NP+1) = Uy(2) + left_rightBC(Y)
!Uz(NP+1) = Uz(2) + left_rightBC(Z)

!Ux_s(1) = Ux_s(NP-1)
!Uy_s(1) = Uy_s(NP-1)
!Uz_s(1) = Uz_s(NP-1)
!SIGMA(1) = SIGMA(NP-1)

!Ux_s(NP+1) = Ux_s(2)
!Uy_s(NP+1) = Uy_s(2)
!Uz_s(NP+1) = Uz_s(2)
!SIGMA(NP+1) = SIGMA(2)

!Ux_ss(1) = Ux_ss(NP-1)
!Uy_ss(1) = Uy_ss(NP-1)
!Uz_ss(1) = Uz_ss(NP-1)

!Ux_ss(NP+1) = Ux_ss(2)
!Uy_ss(NP+1) = Uy_ss(2)
!Uz_ss(NP+1) = Uz_ss(2)

!Ux(NP) = Ux(1) + left_rightBC(X)
!Uy(NP) = Uy(1) + left_rightBC(Y)
!Uz(NP) = Uz(1) + left_rightBC(Z)

!Ux_s(NP) = Ux_s(1)
!Uy_s(NP) = Uy_s(1)
!Uz_s(NP) = Uz_s(1)
!SIGMA(NP) = SIGMA(1)

!Ux_ss(NP) = Ux_ss(1)
!Uy_ss(NP) = Uy_ss(1)
!Uz_ss(NP) = Uz_ss(1)

END IF

IF (step == 1 .AND. flowCondition == 1) THEN

READ(1, *) Ub
PRINT *, Ub

END IF

IF (step == 1 .AND. flowCondition == 2) THEN

READ(1, *) maxV
READ(1, *) yShear
PRINT *, maxV
PRINT *, yShear

END IF

IF (step == 1 .AND. flowCondition == 3) THEN

READ(1, *) interpArraySize
PRINT *, interpArraySize

ALLOCATE(yPlus(interpArraySize))
ALLOCATE(uVelBL(interpArraySize))
ALLOCATE(vVelBL(interpArraySize))
ALLOCATE(wVelBL(interpArraySize))

DO i = 1, interpArraySize

READ(1, *) yPlus(i), uVelBL(i), vVelBL(i), wVelBL(i)

END DO

END IF

IF (step == 1) THEN

! Read other data points
DO i = 1, NP+2
READ(1, *) Ux(i), Uy(i), Uz(i), Ux_s(i), Uy_s(i), Uz_s(i), SIGMA(i)
!PRINT *, i, Ux(i), Uy(i), Uz(i), Ux_s(i), Uy_s(i), Uz_s(i), SIGMA(i)
END DO

! Write the first time step
OPEN(2, FILE = 'fout.dat', STATUS = 'old', POSITION = 'append')
DO i = 1, NP+2
WRITE(2, *) Ux(i), Uy(i), Uz(i)
END DO
CLOSE(2)

CLOSE(1)

END IF

! Compute velocity

!DO i = 1, NP+2
!PRINT *, Ux(i), Uy(i), Uz(i)
!!PRINT *, Ux_s(i), Uy_s(i), Uz_s(i)
!!PRINT *, SIGMA(i), 
!!PRINT *, ux_tmp(i), uy_tmp(i), uz_tmp(i), uux_s(i), uuy_s(i), uuz_s(i), sigma_u(i) ! unaffected
!END DO

VxVyVz = uttm(additional, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA, NP, n_b, image, &
delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, ux_tmp, uy_tmp, &
uz_tmp, uux_s, uuy_s, uuz_s, sigma_u)

!DO i = 1, NP+2
!PRINT *, VxVyVz(i, :)
!PRINT *, Ux(i), Uy(i), Uz(i)
!PRINT *, Ux_s(i), Uy_s(i), Uz_s(i)
!PRINT *, SIGMA(i), 
!PRINT *, ux_tmp(i), uy_tmp(i), uz_tmp(i), uux_s(i), uuy_s(i), uuz_s(i), sigma_u(i) ! unaffected
!END DO

!PRINT *, sigma_u

! Compute position based on flow condition

VxVyVz = -VxVyVz

IF (flowCondition == 1) THEN

VxVyVz(:, X) = VxVyVz(:, X) + Ub

END IF

IF (flowCondition == 2) THEN

DO i = 1, NP+2

IF (Uy(i) > yShear) THEN

VxVyVz(i, X) = VxVyVz(i, X) + maxV

ELSE

VxVyVz(i, X) = VxVyVz(i, X) + (Uy(i) * maxV / yShear)

END IF

END DO

END IF

! For the first time step use, RK45

IF (step == 1) THEN

! -------- Compute k1 --------

vec_local = uttm(additional, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA, NP, n_b, image, &
delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, ux_tmp, uy_tmp, &
uz_tmp, uux_s, uuy_s, uuz_s, sigma_u)

! Compute position based on flow condition

IF (flowCondition == 1) THEN

vec_local(:, X) = vec_local(:, X) + Ub

END IF

IF (flowCondition == 2) THEN

DO i = 1, NP+2

IF (Uy(i) > yShear) THEN

vec_local(i, X) = vec_local(i, X) + maxV

ELSE

vec_local(i, X) = vec_local(i, X) + (Uy(i) * maxV / yShear)

END IF

END DO

END IF

k1 = dt * vec_local

! -------- Compute k2 --------

vec_local(:, :) = 0.
additional(:, :) = 0.25*k1(:, :)

vec_local = uttm(additional, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA, NP, n_b, image, &
delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, ux_tmp, uy_tmp, &
uz_tmp, uux_s, uuy_s, uuz_s, sigma_u)

! Compute position based on flow condition

IF (flowCondition == 1) THEN

vec_local(:, X) = vec_local(:, X) + Ub

END IF

IF (flowCondition == 2) THEN

DO i = 1, NP+2

IF (Uy(i) > yShear) THEN

vec_local(i, X) = vec_local(i, X) + maxV

ELSE

vec_local(i, X) = vec_local(i, X) + (Uy(i) * maxV / yShear)

END IF

END DO

END IF

k2 = dt * vec_local

! -------- Compute k3 --------

vec_local(:, :) = 0.
additional(:, :) = 0.

additional = (0.09375*k1) + (0.28125*k2)

vec_local = uttm(additional, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA, NP, n_b, image, &
delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, ux_tmp, uy_tmp, &
uz_tmp, uux_s, uuy_s, uuz_s, sigma_u)

! Compute position based on flow condition

IF (flowCondition == 1) THEN

vec_local(:, X) = vec_local(:, X) + Ub

END IF

IF (flowCondition == 2) THEN

DO i = 1, NP+2

IF (Uy(i) > yShear) THEN

vec_local(i, X) = vec_local(i, X) + maxV

ELSE

vec_local(i, X) = vec_local(i, X) + (Uy(i) * maxV / yShear)

END IF

END DO

END IF

k3 = dt * vec_local

! -------- Compute k4 --------

vec_local(:, :) = 0.
additional(:, :) = 0.

additional = (0.8793809740555303*k1) - (3.277196176604461*k2) + (3.3208921256258535*k3)

vec_local = uttm(additional, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA, NP, n_b, image, &
delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, ux_tmp, uy_tmp, &
uz_tmp, uux_s, uuy_s, uuz_s, sigma_u)

! Compute position based on flow condition

IF (flowCondition == 1) THEN

vec_local(:, X) = vec_local(:, X) + Ub

END IF

IF (flowCondition == 2) THEN

DO i = 1, NP+2

IF (Uy(i) > yShear) THEN

vec_local(i, X) = vec_local(i, X) + maxV

ELSE

vec_local(i, X) = vec_local(i, X) + (Uy(i) * maxV / yShear)

END IF

END DO

END IF

k4 = dt * vec_local

! -------- Compute k5 --------

vec_local(:, :) = 0.
additional(:, :) = 0.

additional = (2.0324074074074074*k1) - (8.0*k2) + (7.173489278752436*k3) - (0.20589668615984405*k4)

vec_local = uttm(additional, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA, NP, n_b, image, &
delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, ux_tmp, uy_tmp, &
uz_tmp, uux_s, uuy_s, uuz_s, sigma_u)

! Compute position based on flow condition

IF (flowCondition == 1) THEN

vec_local(:, X) = vec_local(:, X) + Ub

END IF

IF (flowCondition == 2) THEN

DO i = 1, NP+2

IF (Uy(i) > yShear) THEN

vec_local(i, X) = vec_local(i, X) + maxV

ELSE

vec_local(i, X) = vec_local(i, X) + (Uy(i) * maxV / yShear)

END IF

END DO

END IF

k5 = dt * vec_local

! -------- Compute k6 --------

vec_local(:, :) = 0.
additional(:, :) = 0.

additional = (-0.2962962962962963*k1) + (2.0*k2) - (1.3816764132553607*k3) + (0.4529727095516569*k4) - (0.275*k5)

vec_local = uttm(additional, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, SIGMA, NP, n_b, image, &
delta_ttm_tmp, sigma_1, sigma_2, gamma_tmp, ds, ux_tmp, uy_tmp, &
uz_tmp, uux_s, uuy_s, uuz_s, sigma_u)

! Compute position based on flow condition

IF (flowCondition == 1) THEN

vec_local(:, X) = vec_local(:, X) + Ub

END IF

IF (flowCondition == 2) THEN

DO i = 1, NP+2

IF (Uy(i) > yShear) THEN

vec_local(i, X) = vec_local(i, X) + maxV

ELSE

vec_local(i, X) = vec_local(i, X) + (Uy(i) * maxV / yShear)

END IF

END DO

END IF

k6 = dt * vec_local

! Add all
vec_add = 0.11851851851851852*k1 + 0.5189863547758284*k3 + 0.5061314903420167*k4 - 0.18*k5 + 0.03636363636363636*k6

Ux_0 = Ux + vec_add(:, X)
Uy_0 = Uy + vec_add(:, Y)
Uz_0 = Uz + vec_add(:, Z)

ELSE IF (step == 2) THEN

vec_add = dt * (((1.5) * VxVyVz) + ((-0.5) * VxVyVz_0))

Ux_0 = Ux + vec_add(:, X)
Uy_0 = Uy + vec_add(:, Y)
Uz_0 = Uz + vec_add(:, Z)

!PRINT *, vec_add(:, X)

ELSE IF (step == 3) THEN

vec_add = dt * (((1.9166666666666667) * VxVyVz) + ((-1.3333333333333333) * VxVyVz_0) + ((0.4166666666666667) * VxVyVz_m))

Ux_0 = Ux + vec_add(:, X)
Uy_0 = Uy + vec_add(:, Y)
Uz_0 = Uz + vec_add(:, Z)

ELSE IF (step == 4) THEN

vec_add = dt * (((2.2916666666666665) * VxVyVz) + ((-2.4583333333333335) * VxVyVz_0) + ((1.5416666666666667) &
* VxVyVz_m) + ((-0.375) * VxVyVz_m1))

Ux_0 = Ux + vec_add(:, X)
Uy_0 = Uy + vec_add(:, Y)
Uz_0 = Uz + vec_add(:, Z)

ELSE

vec_add = dt * (((2.640277777777778) * VxVyVz) + ((-3.852777777777778) * VxVyVz_0) + ((3.6333333333333333) * VxVyVz_m) &
+ ((-1.7694444444444444) * VxVyVz_m1) + ((0.3486111111111111) * VxVyVz_m2))

Ux_0 = Ux + vec_add(:, X)
Uy_0 = Uy + vec_add(:, Y)
Uz_0 = Uz + vec_add(:, Z)

END IF ! Other time steps are run with Adam Bashforth scheme

! Update Ux, Uy, Uz

Ux(:) = 0.
Uy(:) = 0.
Uz(:) = 0.

! Update others for Adam Bashforth Scheme

VxVyVz_0 = VxVyVz
VxVyVz_m = VxVyVz_0
VxVyVz_m1 = VxVyVz_m
VxVyVz_m2 = VxVyVz_m1

DO i = 1, NP+2

Ux(i) = Ux(i) + Ux_0(i)
Uy(i) = Uy(i) + Uy_0(i)
Uz(i) = Uz(i) + Uz_0(i)

END DO

IF (MOD(step, writeEveryNSteps) == 0) THEN

! Write data to file
OPEN(2, FILE = 'fout.dat', STATUS = 'old', POSITION = 'append')
DO i = 1, NP+2
WRITE(2, *) Ux(i), Uy(i), Uz(i)
END DO
CLOSE(2)

END IF

! Find new length

S_tmp = 0.

DO i = 1, NP

S_tmp = S_tmp + SIGMA(i)

END DO

! Find new stretching terms

S = S_tmp * ds

S0_o_S = S_0 / S
int_S_o_S0 = int_S_o_S0 + (dt * (1 / S0_o_S))

! Find new core thickness

delta_bar = delta_0_bar_param**2
delta_bar = delta_bar + 4 * nu_bar_param * int_S_o_S0
delta_bar = SQRT(delta_bar * S0_o_S)

! Find new c_v. c_w is assumed 0

c_v = CV(delta_bar)

!PRINT *, 'S_tmp:', S_tmp
!PRINT *, 'S0_o_S:', S0_o_S
!PRINT *, 'int_S_o_S0:', int_S_o_S0
!PRINT *, SIGMA

PRINT *, 'Time step:', step, '/', nsteps
PRINT *, 'Time (actual): ', step * dt 

!CALL CPU_TIME(elapsed)

!PRINT *, (elapsed - start) * (nsteps - step)

END DO ! Main time loop

!! Compute uttm for all points and write output
!OPEN(2, FILE = 'fout.dat', STATUS = 'replace')

!DO i = 1, NP+1
!WRITE(2, *) VxVyVz(i, X), VxVyVz(i, Y), VxVyVz(i, Z)
!PRINT *, Ux(i), Uy(i), Uz(i)
!END DO
!CLOSE(2)

!PRINT *, 'Done..'

! end timer

!CALL CPU_TIME(finish)

!PRINT *, 'Total time', finish-start

END PROGRAM uttm_full
