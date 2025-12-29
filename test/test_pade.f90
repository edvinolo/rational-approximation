module test_pade
use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
use testdrive, only: new_unittest, unittest_type, error_type, check
use kinds
use robust_pade, only: pade_re, pade_cp
implicit none

private

real(sp), parameter :: rel_tol_sp = 1.0e-5_sp
real(dp), parameter :: rel_tol_dp = 1.0e-13_dp

public :: collect_pade_suite

contains

subroutine collect_pade_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
        new_unittest('pade_exp_re_sp test', pade_exp_re_sp_test),&
        new_unittest('pade_exp_re_dp test', pade_exp_re_dp_test),&
        new_unittest('pade_exp_cp_sp test', pade_exp_cp_sp_test),&
        new_unittest('pade_exp_cp_dp test', pade_exp_cp_dp_test),&
        new_unittest('pade_poly_re_sp test', pade_poly_re_sp_test),&
        new_unittest('pade_poly_re_dp test', pade_poly_re_dp_test),&
        new_unittest('pade_poly_cp_sp test', pade_poly_cp_sp_test),&
        new_unittest('pade_poly_cp_dp test', pade_poly_cp_dp_test),&
        new_unittest('noise_re_sp_test', noise_re_sp_test),&
        new_unittest('noise_re_dp_test', noise_re_dp_test),&
        new_unittest('noise_cp_sp_test', noise_cp_sp_test),&
        new_unittest('noise_cp_dp_test', noise_cp_dp_test)&
    ]
end subroutine collect_pade_suite

subroutine pade_exp_re_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp

     integer, parameter :: N_coeff = 6
    integer, parameter :: m = 3
    integer, parameter :: n = 2

    real(wp), parameter :: c_re(N_coeff) = 1.0_wp/[1.0_wp,1.0_wp,2.0_wp,6.0_wp,24.0_wp,120.0_wp]
    type(pade_re(wp)) :: pade

    real(wp), parameter :: p_ref(4) = [1.0_wp,3.0_wp/5.0_wp,3.0_wp/20.0_wp,1.0_wp/60.0_wp]
    real(wp), parameter :: q_ref(3) = [1.0_wp,-2.0_wp/5.0_wp,1.0_wp/20.0_wp]

    real(wp) :: y, y_vec(2)
    real(wp), parameter :: y_ref = 2.71794871794872_wp
    real(wp), parameter :: y_ref_vec(2) = [2.71794871794872_wp,1.99997203177995_wp]

    real(wp) :: rel_error

    call pade%init(m,n,c_re)

    call check(error, m == pade%mu)
    if (allocated(error)) return

    call check(error, n == pade%nu)
    if (allocated(error)) return

    rel_error = maxval(abs(p_ref-pade%p)/abs(pade%p))
    call check(error,rel_error <= 6*rel_tol_sp)
    if (allocated(error)) return

    rel_error = maxval(abs(q_ref-pade%q)/abs(pade%q))
    call check(error,rel_error <= 6*rel_tol_sp)
    if (allocated(error)) return

    y = pade%eval(1.0_wp)
    y_vec = pade%eval([1.0_wp,log(2.0_wp)])

    rel_error = abs(y_ref-y)/abs(y)
    call check(error,rel_error <= rel_tol_sp)
    if (allocated(error)) return

    rel_error = maxval(abs(y_ref_vec-y_vec)/abs(y_vec))
    call check(error,rel_error <= rel_tol_sp)
    if (allocated(error)) return
end subroutine pade_exp_re_sp_test

subroutine pade_exp_cp_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp

    integer, parameter :: N_coeff = 6
    integer, parameter :: m = 3
    integer, parameter :: n = 2

    real(wp), parameter :: c_re(N_coeff) = 1.0_wp/[1.0_wp,1.0_wp,2.0_wp,6.0_wp,24.0_wp,120.0_wp]
    complex(wp), parameter :: c_cp(N_coeff) = cmplx(c_re,y=0,kind=wp)
    type(pade_cp(wp)) :: pade

    complex(wp), parameter :: p_ref(4) = [1.0_wp,3.0_wp/5.0_wp,3.0_wp/20.0_wp,1.0_wp/60.0_wp]
    complex(wp), parameter :: q_ref(3) = [1.0_wp,-2.0_wp/5.0_wp,1.0_wp/20.0_wp]

    complex(wp) :: y, y_vec(2)
    complex(wp), parameter :: y_ref = (2.71794871794872_wp,0.0_wp)
    complex(wp), parameter :: y_ref_vec(2) = [(1.44687791239515_wp,0.790462876669773_wp),(2.38615664845173_wp,-1.30327868852459_wp)]

    real(wp) :: rel_error

    call pade%init(m,n,c_cp)

    call check(error, m == pade%mu)
    if (allocated(error)) return

    call check(error, n == pade%nu)
    if (allocated(error)) return

    rel_error = maxval(abs(p_ref-pade%p)/abs(pade%p))
    call check(error,rel_error <= 6*rel_tol_sp)
    if (allocated(error)) return

    rel_error = maxval(abs(q_ref-pade%q)/abs(pade%q))
    call check(error,rel_error <= 6*rel_tol_sp)
    if (allocated(error)) return

    y = pade%eval((1.0_wp,0.0_wp))
    y_vec = pade%eval([(0.5_wp,0.5_wp),(1.0_wp,-0.5_wp)])

    rel_error = abs(y_ref-y)/abs(y)
    call check(error,rel_error <= rel_tol_sp)
    if (allocated(error)) return

    rel_error = maxval(abs(y_ref_vec-y_vec)/abs(y_vec))
    call check(error,rel_error <= rel_tol_sp)
    if (allocated(error)) return
end subroutine pade_exp_cp_sp_test

subroutine pade_exp_re_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp

     integer, parameter :: N_coeff = 6
    integer, parameter :: m = 3
    integer, parameter :: n = 2

    real(wp), parameter :: c_re(N_coeff) = 1.0_wp/[1.0_wp,1.0_wp,2.0_wp,6.0_wp,24.0_wp,120.0_wp]
    type(pade_re(wp)) :: pade

    real(wp), parameter :: p_ref(4) = [1.0_wp,3.0_wp/5.0_wp,3.0_wp/20.0_wp,1.0_wp/60.0_wp]
    real(wp), parameter :: q_ref(3) = [1.0_wp,-2.0_wp/5.0_wp,1.0_wp/20.0_wp]

    real(wp) :: y, y_vec(2)
    real(wp), parameter :: y_ref = 2.71794871794872_wp
    real(wp), parameter :: y_ref_vec(2) = [2.71794871794872_wp,1.99997203177995_wp]

    real(wp) :: rel_error

    call pade%init(m,n,c_re)

    call check(error, m == pade%mu)
    if (allocated(error)) return

    call check(error, n == pade%nu)
    if (allocated(error)) return

    rel_error = maxval(abs(p_ref-pade%p)/abs(pade%p))
    call check(error,rel_error <= 6*rel_tol_dp)
    if (allocated(error)) return

    rel_error = maxval(abs(q_ref-pade%q)/abs(pade%q))
    call check(error,rel_error <= 6*rel_tol_dp)
    if (allocated(error)) return

    y = pade%eval(1.0_wp)
    y_vec = pade%eval([1.0_wp,log(2.0_wp)])

    rel_error = abs(y_ref-y)/abs(y)
    call check(error,rel_error <= rel_tol_dp)
    if (allocated(error)) return

    rel_error = maxval(abs(y_ref_vec-y_vec)/abs(y_vec))
    call check(error,rel_error <= rel_tol_dp)
    if (allocated(error)) return
end subroutine pade_exp_re_dp_test

subroutine pade_exp_cp_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp

    integer, parameter :: N_coeff = 6
    integer, parameter :: m = 3
    integer, parameter :: n = 2

    real(wp), parameter :: c_re(N_coeff) = 1.0_wp/[1.0_wp,1.0_wp,2.0_wp,6.0_wp,24.0_wp,120.0_wp]
    complex(wp), parameter :: c_cp(N_coeff) = cmplx(c_re,y=0,kind=wp)
    type(pade_cp(wp)) :: pade

    complex(wp), parameter :: p_ref(4) = [1.0_wp,3.0_wp/5.0_wp,3.0_wp/20.0_wp,1.0_wp/60.0_wp]
    complex(wp), parameter :: q_ref(3) = [1.0_wp,-2.0_wp/5.0_wp,1.0_wp/20.0_wp]

    complex(wp) :: y, y_vec(2)
    complex(wp), parameter :: y_ref = (2.71794871794872_wp,0.0_wp)
    complex(wp), parameter :: y_ref_vec(2) = [(1.44687791239515_wp,0.790462876669773_wp),(2.38615664845173_wp,-1.30327868852459_wp)]

    real(wp) :: rel_error

    call pade%init(m,n,c_cp)

    call check(error, m == pade%mu)
    if (allocated(error)) return

    call check(error, n == pade%nu)
    if (allocated(error)) return

    rel_error = maxval(abs(p_ref-pade%p)/abs(pade%p))
    call check(error,rel_error <= 6*rel_tol_dp)
    if (allocated(error)) return

    rel_error = maxval(abs(q_ref-pade%q)/abs(pade%q))
    call check(error,rel_error <= 6*rel_tol_dp)
    if (allocated(error)) return

    y = pade%eval((1.0_wp,0.0_wp))
    y_vec = pade%eval([(0.5_wp,0.5_wp),(1.0_wp,-0.5_wp)])

    rel_error = abs(y_ref-y)/abs(y)
    call check(error,rel_error <= rel_tol_dp)
    if (allocated(error)) return

    rel_error = maxval(abs(y_ref_vec-y_vec)/abs(y_vec))
    call check(error,rel_error <= rel_tol_dp)
    if (allocated(error)) return
end subroutine pade_exp_cp_dp_test


subroutine pade_poly_re_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp

    real(wp) :: c(35) = 0
    integer, parameter :: m = 19
    integer, parameter :: n = 15
    type(pade_re(wp)) :: pade

    integer, parameter :: mu = 11
    integer, parameter :: nu = 0
    real(wp) :: p_ref(mu+1) = 0
    real(wp) :: q_ref(nu+1) = 0
    real(wp) :: rel_error

    c(1) = 1.0_wp
    c(5) = 1.0_wp
    c(9) = 1.0_wp
    c(12) = 1.0_wp

    p_ref(1) = 1.0_wp
    p_ref(5) = 1.0_wp
    p_ref(9) = 1.0_wp
    p_ref(12) = 1.0_wp

    q_ref(1) = 1.0_wp

    call pade%init(m,n,c)

    call check(error, pade%mu == mu)
    if (allocated(error)) return

    call check(error, pade%nu == nu)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%p-p_ref)/abs(pade%p))
    call check(error,rel_error <= 0.1*rel_tol_sp)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%q-q_ref))
    call check(error,rel_error <= 0.1*rel_tol_sp)
    if (allocated(error)) return

end subroutine pade_poly_re_sp_test
subroutine pade_poly_re_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp

    real(wp) :: c(35) = 0
    integer, parameter :: m = 19
    integer, parameter :: n = 15
    type(pade_re(wp)) :: pade

    integer, parameter :: mu = 11
    integer, parameter :: nu = 0
    real(wp) :: p_ref(mu+1) = 0
    real(wp) :: q_ref(nu+1) = 0
    real(wp) :: rel_error

    c(1) = 1.0_wp
    c(5) = 1.0_wp
    c(9) = 1.0_wp
    c(12) = 1.0_wp

    p_ref(1) = 1.0_wp
    p_ref(5) = 1.0_wp
    p_ref(9) = 1.0_wp
    p_ref(12) = 1.0_wp

    q_ref(1) = 1.0_wp

    call pade%init(m,n,c)

    call check(error, pade%mu == mu)
    if (allocated(error)) return

    call check(error, pade%nu == nu)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%p-p_ref)/abs(pade%p))
    call check(error,rel_error <= 0.1*rel_tol_dp)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%q-q_ref))
    call check(error,rel_error <= 0.1*rel_tol_dp)
    if (allocated(error)) return

end subroutine pade_poly_re_dp_test
subroutine pade_poly_cp_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp

    complex(wp) :: c(35) = 0
    integer, parameter :: m = 19
    integer, parameter :: n = 15
    type(pade_cp(wp)) :: pade

    integer, parameter :: mu = 11
    integer, parameter :: nu = 0
    complex(wp) :: p_ref(mu+1) = 0
    complex(wp) :: q_ref(nu+1) = 0
    real(wp) :: rel_error

    c(1) = 1.0_wp
    c(5) = 1.0_wp
    c(9) = 1.0_wp
    c(12) = 1.0_wp

    p_ref(1) = 1.0_wp
    p_ref(5) = 1.0_wp
    p_ref(9) = 1.0_wp
    p_ref(12) = 1.0_wp

    q_ref(1) = 1.0_wp

    call pade%init(m,n,c)

    call check(error, pade%mu == mu)
    if (allocated(error)) return

    call check(error, pade%nu == nu)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%p-p_ref)/abs(pade%p))
    call check(error,rel_error <= 0.1*rel_tol_sp)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%q-q_ref))
    call check(error,rel_error <= 0.1*rel_tol_sp)
    if (allocated(error)) return

end subroutine pade_poly_cp_sp_test
subroutine pade_poly_cp_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp

    complex(wp) :: c(35) = 0
    integer, parameter :: m = 19
    integer, parameter :: n = 15
    type(pade_cp(wp)) :: pade

    integer, parameter :: mu = 11
    integer, parameter :: nu = 0
    complex(wp) :: p_ref(mu+1) = 0
    complex(wp) :: q_ref(nu+1) = 0
    real(wp) :: rel_error

    c(1) = 1.0_wp
    c(5) = 1.0_wp
    c(9) = 1.0_wp
    c(12) = 1.0_wp

    p_ref(1) = 1.0_wp
    p_ref(5) = 1.0_wp
    p_ref(9) = 1.0_wp
    p_ref(12) = 1.0_wp

    q_ref(1) = 1.0_wp

    call pade%init(m,n,c)

    call check(error, pade%mu == mu)
    if (allocated(error)) return

    call check(error, pade%nu == nu)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%p-p_ref)/abs(pade%p))
    call check(error,rel_error <= 0.1*rel_tol_dp)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%q-q_ref))
    call check(error,rel_error <= 0.1*rel_tol_dp)
    if (allocated(error)) return

end subroutine pade_poly_cp_dp_test

subroutine noise_re_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    type(pade_re(sp)) :: pade

    integer, parameter :: m = 8
    integer, parameter :: n = 6
    integer, parameter :: i_max = n+m+1

    real(sp) :: c(i_max) = 0.0_sp
    real(sp), parameter :: tol = 1.0e-5_sp
    real(sp), parameter :: noise_lvl = 1.0e-6_sp
    real(sp) :: noise

    integer, parameter :: mu_ref = 0
    integer, parameter :: nu_ref = 1
    real(sp), parameter :: p_ref = 1.0_sp
    real(sp), parameter :: q_ref(2) = [1.0_sp,-1.0_sp]
    real(sp) :: rel_error

    integer :: i

    call random_seed()
    do i = 1,i_max
        call random_number(noise)
        noise = noise_lvl*2.0_sp*(noise-0.5_sp)
        c(i) = 1.0_sp + noise
    end do

    call pade%init(m,n,c,tol=tol)

    ! Function should be 1/(1-z) with some noise
    call check(error, pade%mu == 0)
    if (allocated(error)) return

    call check(error, pade%nu == 1)
    if (allocated(error)) return

    rel_error = abs(pade%p(1)-p_ref)/abs(p_ref)
    call check(error, rel_error <= tol)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%q-q_ref)/abs(q_ref))
    call check(error, rel_error <= tol)
    if (allocated(error)) return
end subroutine noise_re_sp_test

subroutine noise_re_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    type(pade_re(dp)) :: pade

    integer, parameter :: m = 8
    integer, parameter :: n = 6
    integer, parameter :: i_max = n+m+1

    real(dp) :: c(i_max) = 0.0_dp
    real(dp), parameter :: tol = 1.0e-5_dp
    real(dp), parameter :: noise_lvl = 1.0e-6_dp
    real(dp) :: noise

    integer, parameter :: mu_ref = 0
    integer, parameter :: nu_ref = 1
    real(dp), parameter :: p_ref = 1.0_dp
    real(dp), parameter :: q_ref(2) = [1.0_dp,-1.0_dp]
    real(dp) :: rel_error

    integer :: i

    call random_seed()
    do i = 1,i_max
        call random_number(noise)
        noise = noise_lvl*2.0_dp*(noise-0.5_dp)
        c(i) = 1.0_dp + noise
    end do

    call pade%init(m,n,c,tol=tol)

    ! Function should be 1/(1-z) with some noise
    call check(error, pade%mu == 0)
    if (allocated(error)) return

    call check(error, pade%nu == 1)
    if (allocated(error)) return

    rel_error = abs(pade%p(1)-p_ref)/abs(p_ref)
    call check(error, rel_error <= tol)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%q-q_ref)/abs(q_ref))
    call check(error, rel_error <= tol)
    if (allocated(error)) return
end subroutine noise_re_dp_test

subroutine noise_cp_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    type(pade_cp(sp)) :: pade

    integer, parameter :: m = 8
    integer, parameter :: n = 6
    integer, parameter :: i_max = n+m+1

    complex(sp) :: c(i_max) = 0.0_sp
    real(sp), parameter :: tol = 1.0e-5_sp
    real(sp), parameter :: noise_lvl = 1.0e-6_sp
    real(sp) :: noise

    integer, parameter :: mu_ref = 0
    integer, parameter :: nu_ref = 1
    complex(sp), parameter :: p_ref = 1.0_sp
    complex(sp), parameter :: q_ref(2) = [1.0_sp,-1.0_sp]
    real(sp) :: rel_error

    integer :: i

    call random_seed()
    do i = 1,i_max
        call random_number(noise)
        noise = noise_lvl*2.0_sp*(noise-0.5_sp)
        c(i) = 1.0_sp + noise
        call random_number(noise)
        noise = noise_lvl*2.0_sp*(noise-0.5_sp)
        c(i) = c(i) + noise*(0.0_sp,1.0_sp)
    end do

    call pade%init(m,n,c,tol=tol)

    ! Function should be 1/(1-z) with some noise
    call check(error, pade%mu == 0)
    if (allocated(error)) return

    call check(error, pade%nu == 1)
    if (allocated(error)) return

    rel_error = abs(pade%p(1)-p_ref)/abs(p_ref)
    call check(error, rel_error <= tol)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%q-q_ref)/abs(q_ref))
    call check(error, rel_error <= tol)
    if (allocated(error)) return
end subroutine noise_cp_sp_test

subroutine noise_cp_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    type(pade_cp(dp)) :: pade

    integer, parameter :: m = 8
    integer, parameter :: n = 6
    integer, parameter :: i_max = n+m+1

    complex(dp) :: c(i_max) = 0.0_dp
    real(dp), parameter :: tol = 1.0e-5_dp
    real(dp), parameter :: noise_lvl = 1.0e-6_dp
    real(dp) :: noise

    integer, parameter :: mu_ref = 0
    integer, parameter :: nu_ref = 1
    complex(dp), parameter :: p_ref = 1.0_dp
    complex(dp), parameter :: q_ref(2) = [1.0_dp,-1.0_dp]
    real(dp) :: rel_error

    integer :: i

    call random_seed()
    do i = 1,i_max
        call random_number(noise)
        noise = noise_lvl*2.0_dp*(noise-0.5_dp)
        c(i) = 1.0_dp + noise
        call random_number(noise)
        noise = noise_lvl*2.0_dp*(noise-0.5_dp)
        c(i) = c(i) + noise*(0.0_dp,1.0_dp)
    end do

    call pade%init(m,n,c,tol=tol)

    ! Function should be 1/(1-z) with some noise
    call check(error, pade%mu == 0)
    if (allocated(error)) return

    call check(error, pade%nu == 1)
    if (allocated(error)) return

    rel_error = abs(pade%p(1)-p_ref)/abs(p_ref)
    call check(error, rel_error <= tol)
    if (allocated(error)) return

    rel_error = maxval(abs(pade%q-q_ref)/abs(q_ref))
    call check(error, rel_error <= tol)
    if (allocated(error)) return
end subroutine noise_cp_dp_test

end module test_pade