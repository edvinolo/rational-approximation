module test_polynomial_eval
use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
use testdrive, only: new_unittest, unittest_type, error_type, check
use kinds
use polynomial_eval, only: poly_eval, rat_eval
implicit none

private

real(sp), parameter :: rel_tol_sp = 1.0e-5_sp
real(dp), parameter :: rel_tol_dp = 1.0e-13_dp

public :: collect_polynomial_eval_suite

contains

subroutine collect_polynomial_eval_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
        new_unittest('poly_eval_re_sp test', poly_eval_re_sp_test),&
        new_unittest('poly_eval_re_dp test', poly_eval_re_dp_test),&
        new_unittest('poly_eval_cp_sp test', poly_eval_cp_sp_test),&
        new_unittest('poly_eval_cp_dp test', poly_eval_cp_dp_test),&
        new_unittest('rat_eval_re_sp test', rat_eval_re_sp_test),&
        new_unittest('rat_eval_re_dp test', rat_eval_re_dp_test),&
        new_unittest('rat_eval_cp_sp test', rat_eval_cp_sp_test),&
        new_unittest('rat_eval_cp_dp test', rat_eval_cp_dp_test)&
    ]
end subroutine collect_polynomial_eval_suite

subroutine poly_eval_re_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp

    real(wp) :: c(9) = 0.0_wp
    real(wp), parameter :: x_ref(3) = [-2.0_wp,1.0_wp,3.0_wp]
    real(wp), parameter :: y_ref(3) = [-149.0_wp, 46.0_wp,8706.0_wp]
    real(wp) :: y(3)

    integer :: i
    real(wp) :: rel_err

    c(1) = 75.0_wp
    c(3) = -40.0_wp
    c(6) = 10.0_wp
    c(9) = 1.0_wp

    do i = 1,3
        y(i) = poly_eval(c,x_ref(i))
    end do

    rel_err = maxval(abs(y_ref-y)/abs(y))

    call check(error,rel_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine poly_eval_re_sp_test

subroutine poly_eval_cp_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp

    complex(wp) :: c(9) = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_ref(3) = [(-2.0_wp,2.0_wp),(1.0_wp,1.0_wp),(3.0_wp,-3.0_wp)]
    complex(wp), parameter :: y_ref(3) = [(2075.0_wp,3136.0_wp), (-145.0_wp,-104.0_wp),(-8025.0_wp,115416.0_wp)]
    complex(wp) :: y(3)

    integer :: i
    real(wp) :: rel_err

    c(1) = (75.0_wp,0.0_wp)
    c(3) = (-40.0_wp,90.0_wp)
    c(6) = (10.0_wp,0.0_wp)
    c(9) = (0.0_wp,1.0_wp)

    do i = 1,3
        y(i) = poly_eval(c,x_ref(i))
    end do

    rel_err = maxval(abs(y_ref-y)/abs(y))

    call check(error,rel_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine poly_eval_cp_sp_test

subroutine rat_eval_re_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp

    real(wp) :: p(6) = 0.0_wp
    real(wp) :: q(3) = 0.0_wp
    real(wp), parameter :: x_ref(3) = [-2.0_wp,1.0_wp,3.0_wp]
    real(wp), parameter :: y_ref(3) = [-33.0_wp/5.0_wp, 0.0_wp,121.0_wp/5.0_wp]
    real(wp) :: y(3)

    integer :: i
    real(wp) :: rel_err

    p(1) = -1.0_wp
    p(6) = 1.0_wp
    q(1) = 1.0_wp
    q(3) = 1.0_wp

    do i = 1,3
        y(i) = rat_eval(p,q,x_ref(i))
    end do

    rel_err = maxval(abs(y_ref-y)/abs(y))

    call check(error,rel_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine rat_eval_re_sp_test

subroutine rat_eval_cp_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp

    complex(wp) :: p(6) = (0.0_wp,0.0_wp)
    complex(wp) :: q(3) = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_ref(3) = [(-2.0_wp,2.0_wp),(1.0_wp,1.0_wp),(3.0_wp,-3.0_wp)]
    real(wp), parameter :: y_real(3) = [133.0_wp/5.0_wp,-9.0_wp/5.0_wp,-967.0_wp/15.0_wp]
    real(wp), parameter :: y_imag(3) = [128.0_wp/5.0_wp,4.0_wp/5.0_wp,-324.0_wp/5.0_wp]
    complex(wp), parameter :: y_ref(3) = cmplx(y_real, y = y_imag, kind = wp)
    complex(wp) :: y(3)

    integer :: i
    real(wp) :: rel_err

    p(1) = (0.0_wp,-5.0_wp)
    p(6) = (1.0_wp,0.0_wp)
    q(1) = (0.0_wp,3.0_wp)
    q(3) = (1.0_wp,0.0_wp)

    do i = 1,3
        y(i) = rat_eval(p,q,x_ref(i))
    end do

    rel_err = maxval(abs(y_ref-y)/abs(y))

    call check(error,rel_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine rat_eval_cp_sp_test

subroutine poly_eval_re_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp

    real(wp) :: c(9) = 0.0_wp
    real(wp), parameter :: x_ref(3) = [-2.0_wp,1.0_wp,3.0_wp]
    real(wp), parameter :: y_ref(3) = [-149.0_wp, 46.0_wp,8706.0_wp]
    real(wp) :: y(3)

    integer :: i
    real(wp) :: rel_err

    c(1) = 75.0_wp
    c(3) = -40.0_wp
    c(6) = 10.0_wp
    c(9) = 1.0_wp

    do i = 1,3
        y(i) = poly_eval(c,x_ref(i))
    end do

    rel_err = maxval(abs(y_ref-y)/abs(y))

    call check(error,rel_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine poly_eval_re_dp_test

subroutine poly_eval_cp_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp

    complex(wp) :: c(9) = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_ref(3) = [(-2.0_wp,2.0_wp),(1.0_wp,1.0_wp),(3.0_wp,-3.0_wp)]
    complex(wp), parameter :: y_ref(3) = [(2075.0_wp,3136.0_wp), (-145.0_wp,-104.0_wp),(-8025.0_wp,115416.0_wp)]
    complex(wp) :: y(3)

    integer :: i
    real(wp) :: rel_err

    c(1) = (75.0_wp,0.0_wp)
    c(3) = (-40.0_wp,90.0_wp)
    c(6) = (10.0_wp,0.0_wp)
    c(9) = (0.0_wp,1.0_wp)

    do i = 1,3
        y(i) = poly_eval(c,x_ref(i))
    end do

    rel_err = maxval(abs(y_ref-y)/abs(y))

    call check(error,rel_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine poly_eval_cp_dp_test

subroutine rat_eval_re_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp

    real(wp) :: p(6) = 0.0_wp
    real(wp) :: q(3) = 0.0_wp
    real(wp), parameter :: x_ref(3) = [-2.0_wp,1.0_wp,3.0_wp]
    real(wp), parameter :: y_ref(3) = [-33.0_wp/5.0_wp, 0.0_wp,121.0_wp/5.0_wp]
    real(wp) :: y(3)

    integer :: i
    real(wp) :: rel_err

    p(1) = -1.0_wp
    p(6) = 1.0_wp
    q(1) = 1.0_wp
    q(3) = 1.0_wp

    do i = 1,3
        y(i) = rat_eval(p,q,x_ref(i))
    end do

    rel_err = maxval(abs(y_ref-y)/abs(y))

    call check(error,rel_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine rat_eval_re_dp_test

subroutine rat_eval_cp_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp

    complex(wp) :: p(6) = (0.0_wp,0.0_wp)
    complex(wp) :: q(3) = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_ref(3) = [(-2.0_wp,2.0_wp),(1.0_wp,1.0_wp),(3.0_wp,-3.0_wp)]
    real(wp), parameter :: y_real(3) = [133.0_wp/5.0_wp,-9.0_wp/5.0_wp,-967.0_wp/15.0_wp]
    real(wp), parameter :: y_imag(3) = [128.0_wp/5.0_wp,4.0_wp/5.0_wp,-324.0_wp/5.0_wp]
    complex(wp), parameter :: y_ref(3) = cmplx(y_real, y = y_imag, kind = wp)
    complex(wp) :: y(3)

    integer :: i
    real(wp) :: rel_err

    p(1) = (0.0_wp,-5.0_wp)
    p(6) = (1.0_wp,0.0_wp)
    q(1) = (0.0_wp,3.0_wp)
    q(3) = (1.0_wp,0.0_wp)

    do i = 1,3
        y(i) = rat_eval(p,q,x_ref(i))
    end do

    rel_err = maxval(abs(y_ref-y)/abs(y))

    call check(error,rel_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine rat_eval_cp_dp_test

end module test_polynomial_eval