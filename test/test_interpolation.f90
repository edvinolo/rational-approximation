module test_interpolation
use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
use testdrive, only: new_unittest, unittest_type, error_type, check
use kinds
use rational_interpolation, only: thiele_interp_re, thiele_interp_cp, MTT_interp_re, MTT_interp_cp
implicit none

private

real(sp), parameter :: rel_tol_sp = 1.0e-5_sp
real(dp), parameter :: rel_tol_dp = 1.0e-13_dp

public :: collect_interpolation_suite

contains

subroutine collect_interpolation_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
        new_unittest('thiele_re_sp test', thiele_re_sp_test),&
        new_unittest('thiele_re_dp test', thiele_re_dp_test),&
        new_unittest('thiele_cp_sp test', thiele_cp_sp_test),&
        new_unittest('thiele_cp_dp test', thiele_cp_dp_test),&
        new_unittest('MTT_re_sp test', MTT_re_sp_test),&
        new_unittest('MTT_re_dp test', MTT_re_dp_test),&
        new_unittest('MTT_cp_sp test', MTT_cp_sp_test),&
        new_unittest('MTT_cp_dp test', MTT_cp_dp_test),&
        new_unittest('MTT_t_test_re_sp',MTT_t_test_re_sp),&
        new_unittest('MTT_t_test_re_dp',MTT_t_test_re_dp),&
        new_unittest('MTT_t_test_cp_sp',MTT_t_test_cp_sp),&
        new_unittest('MTT_t_test_cp_dp',MTT_t_test_cp_dp)&
    ]
end subroutine collect_interpolation_suite

subroutine thiele_re_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: N_int = 8
    integer, parameter :: N = 51
    integer, parameter :: wp = sp
    real(wp), parameter :: dx_int = 0.4_wp
    real(wp), parameter :: dx = 0.13_wp
    real(wp), parameter :: x_int_start = 0.0_wp
    real(wp), parameter :: x_start = -0.2_wp

    real(wp), allocatable :: x_int(:)
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: y_int(:)
    real(wp), allocatable :: y(:)
    type(thiele_interp_re(wp)) :: interp

    integer :: i
    real(wp) :: rel_interp_err

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = cos(x_int(i))
    end do

    call interp%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = interp%eval(x(i))
    end do

    y_int = interp%eval(x_int)
    rel_interp_err = maxval(abs(cos(x_int)-y_int)/abs(y_int))

    call check(error,rel_interp_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine thiele_re_sp_test

subroutine thiele_cp_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: N_int = 8
    integer, parameter :: N = 51
    integer, parameter :: wp = sp
    complex(wp), parameter :: dx_int = (0.4_wp,0.1_wp)
    complex(wp), parameter :: dx = (0.13_wp,0.1_wp)
    complex(wp), parameter :: x_int_start = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_start = (-0.2_wp,-0.2_wp)

    complex(wp), allocatable :: x_int(:)
    complex(wp), allocatable :: x(:)
    complex(wp), allocatable :: y_int(:)
    complex(wp), allocatable :: y(:)
    type(thiele_interp_cp(wp)) :: interp

    integer :: i
    real(wp) :: rel_interp_err

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = cos(x_int(i))
    end do

    call interp%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = interp%eval(x(i))
    end do

    y_int = interp%eval(x_int)
    rel_interp_err = maxval(abs(cos(x_int)-y_int)/abs(y_int))

    call check(error,rel_interp_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine thiele_cp_sp_test

subroutine thiele_re_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: N_int = 8
    integer, parameter :: N = 51
    integer, parameter :: wp = dp
    real(wp), parameter :: dx_int = 0.4_wp
    real(wp), parameter :: dx = 0.13_wp
    real(wp), parameter :: x_int_start = 0.0_wp
    real(wp), parameter :: x_start = -0.2_wp

    real(wp), allocatable :: x_int(:)
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: y_int(:)
    real(wp), allocatable :: y(:)
    type(thiele_interp_re(wp)) :: interp

    integer :: i
    real(wp) :: rel_interp_err

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = cos(x_int(i))
    end do

    call interp%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = interp%eval(x(i))
    end do

    y_int = interp%eval(x_int)
    rel_interp_err = maxval(abs(cos(x_int)-y_int)/abs(y_int))

    call check(error,rel_interp_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine thiele_re_dp_test

subroutine thiele_cp_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: N_int = 8
    integer, parameter :: N = 51
    integer, parameter :: wp = dp
    complex(wp), parameter :: dx_int = (0.4_wp,0.1_wp)
    complex(wp), parameter :: dx = (0.13_wp,0.1_wp)
    complex(wp), parameter :: x_int_start = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_start = (-0.2_wp,-0.2_wp)

    complex(wp), allocatable :: x_int(:)
    complex(wp), allocatable :: x(:)
    complex(wp), allocatable :: y_int(:)
    complex(wp), allocatable :: y(:)
    type(thiele_interp_cp(wp)) :: interp

    integer :: i
    real(wp) :: rel_interp_err

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = cos(x_int(i))
    end do

    call interp%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = interp%eval(x(i))
    end do

    y_int = interp%eval(x_int)
    rel_interp_err = maxval(abs(cos(x_int)-y_int)/abs(y_int))

    call check(error,rel_interp_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine thiele_cp_dp_test

subroutine MTT_re_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: N_int = 8
    integer, parameter :: N = 51
    integer, parameter :: wp = sp
    real(wp), parameter :: dx_int = 0.4_wp
    real(wp), parameter :: dx = 0.13_wp
    real(wp), parameter :: x_int_start = 0.0_wp
    real(wp), parameter :: x_start = -0.2_wp

    real(wp), allocatable :: x_int(:)
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: y_int(:)
    real(wp), allocatable :: y(:)
    type(MTT_interp_re(wp)) :: interp

    integer :: i
    real(wp) :: rel_interp_err

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = cos(x_int(i))
    end do

    call interp%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = interp%eval(x(i))
    end do

    y_int = interp%eval(x_int)
    rel_interp_err = maxval(abs(cos(x_int)-y_int)/abs(y_int))

    call check(error,rel_interp_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine MTT_re_sp_test

subroutine MTT_cp_sp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: N_int = 8
    integer, parameter :: N = 51
    integer, parameter :: wp = sp
    complex(wp), parameter :: dx_int = (0.4_wp,0.1_wp)
    complex(wp), parameter :: dx = (0.13_wp,0.1_wp)
    complex(wp), parameter :: x_int_start = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_start = (-0.2_wp,-0.2_wp)

    complex(wp), allocatable :: x_int(:)
    complex(wp), allocatable :: x(:)
    complex(wp), allocatable :: y_int(:)
    complex(wp), allocatable :: y(:)
    type(MTT_interp_cp(wp)) :: interp

    integer :: i
    real(wp) :: rel_interp_err

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = cos(x_int(i))
    end do

    call interp%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = interp%eval(x(i))
    end do

    y_int = interp%eval(x_int)
    rel_interp_err = maxval(abs(cos(x_int)-y_int)/abs(y_int))

    call check(error,rel_interp_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine MTT_cp_sp_test

subroutine MTT_re_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: N_int = 8
    integer, parameter :: N = 51
    integer, parameter :: wp = dp
    real(wp), parameter :: dx_int = 0.4_wp
    real(wp), parameter :: dx = 0.13_wp
    real(wp), parameter :: x_int_start = 0.0_wp
    real(wp), parameter :: x_start = -0.2_wp

    real(wp), allocatable :: x_int(:)
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: y_int(:)
    real(wp), allocatable :: y(:)
    type(MTT_interp_re(wp)) :: interp

    integer :: i
    real(wp) :: rel_interp_err

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = cos(x_int(i))
    end do

    call interp%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = interp%eval(x(i))
    end do

    y_int = interp%eval(x_int)
    rel_interp_err = maxval(abs(cos(x_int)-y_int)/abs(y_int))

    call check(error,rel_interp_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine MTT_re_dp_test

subroutine MTT_cp_dp_test(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: N_int = 8
    integer, parameter :: N = 51
    integer, parameter :: wp = dp
    complex(wp), parameter :: dx_int = (0.4_wp,0.1_wp)
    complex(wp), parameter :: dx = (0.13_wp,0.1_wp)
    complex(wp), parameter :: x_int_start = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_start = (-0.2_wp,-0.2_wp)

    complex(wp), allocatable :: x_int(:)
    complex(wp), allocatable :: x(:)
    complex(wp), allocatable :: y_int(:)
    complex(wp), allocatable :: y(:)
    type(MTT_interp_cp(wp)) :: interp

    integer :: i
    real(wp) :: rel_interp_err

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = cos(x_int(i))
    end do

    call interp%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = interp%eval(x(i))
    end do

    y_int = interp%eval(x_int)
    rel_interp_err = maxval(abs(cos(x_int)-y_int)/abs(y_int))

    call check(error,rel_interp_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine MTT_cp_dp_test


subroutine MTT_t_test_re_sp(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp
    real(wp), parameter :: x_int(5) = [-1.0_wp,0.0_wp,1.0_wp,2.0_wp,3.0_wp]
    real(wp), parameter :: y_int(5) = (x_int**2-1.0_wp)/(x_int+2.0_wp)
    type(MTT_interp_re(wp)) :: interp

    real(sp) :: abs_err
    integer, parameter :: t_expected = 4

    call interp%init(x_int,y_int)

    call check(error, interp%t == t_expected)
    if (allocated(error)) return

    abs_err = maxval(abs(interp%eval(x_int)-y_int))
    call check(error,abs_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine MTT_t_test_re_sp

subroutine MTT_t_test_re_dp(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp
    real(wp), parameter :: x_int(5) = [-1.0_wp,0.0_wp,1.0_wp,2.0_wp,3.0_wp]
    real(wp), parameter :: y_int(5) = (x_int**2-1.0_wp)/(x_int+2.0_wp)
    type(MTT_interp_re(wp)) :: interp

    real(dp) :: abs_err
    integer, parameter :: t_expected = 4

    call interp%init(x_int,y_int)

    call check(error, interp%t == t_expected)
    if (allocated(error)) return

    abs_err = maxval(abs(interp%eval(x_int)-y_int))
    call check(error,abs_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine MTT_t_test_re_dp


subroutine MTT_t_test_cp_sp(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = sp
    complex(wp), parameter :: x_int(5) = [(-1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(1.0_wp,0.0_wp),(2.0_wp,0.0_wp),(3.0_wp,0.0_wp)]
    complex(wp), parameter :: y_int(5) = (x_int**2-(1.0_wp,0.0_wp))/(x_int+(2.0_wp,0.0_wp))
    type(MTT_interp_cp(wp)) :: interp

    real(sp) :: abs_err
    integer, parameter :: t_expected = 4

    call interp%init(x_int,y_int)

    call check(error, interp%t == t_expected)
    if (allocated(error)) return

    abs_err = maxval(abs(interp%eval(x_int)-y_int))
    call check(error,abs_err <= rel_tol_sp)
    if (allocated(error)) return
end subroutine MTT_t_test_cp_sp

subroutine MTT_t_test_cp_dp(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: wp = dp
    complex(wp), parameter :: x_int(5) = [(-1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(1.0_wp,0.0_wp),(2.0_wp,0.0_wp),(3.0_wp,0.0_wp)]
    complex(wp), parameter :: y_int(5) = (x_int**2-(1.0_wp,0.0_wp))/(x_int+(2.0_wp,0.0_wp))
    type(MTT_interp_cp(wp)) :: interp

    real(dp) :: abs_err
    integer, parameter :: t_expected = 4

    call interp%init(x_int,y_int)

    call check(error, interp%t == t_expected)
    if (allocated(error)) return

    abs_err = maxval(abs(interp%eval(x_int)-y_int))
    call check(error,abs_err <= rel_tol_dp)
    if (allocated(error)) return
end subroutine MTT_t_test_cp_dp

end module test_interpolation