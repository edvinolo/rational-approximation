module polynomial_eval
    use kinds
    implicit none

    private

    public :: poly_eval
    public :: rat_eval

    interface poly_eval
        module procedure :: poly_eval_re_sp
        module procedure :: poly_eval_re_dp
        module procedure :: poly_eval_cp_sp
        module procedure :: poly_eval_cp_dp
    end interface poly_eval

    interface rat_eval
        module procedure :: rat_eval_re_sp
        module procedure :: rat_eval_re_dp
        module procedure :: rat_eval_cp_sp
        module procedure :: rat_eval_cp_dp
    end interface rat_eval

contains

    pure function poly_eval_re_sp(c,x) result(res)
        real(sp), intent(in) :: c(:)
        real(sp), intent(in) :: x
        real(sp) :: res

        res = horner_re_sp(c,x)
    end function

    pure function poly_eval_re_dp(c,x) result(res)
        real(dp), intent(in) :: c(:)
        real(dp), intent(in) :: x
        real(dp) :: res

        res = horner_re_dp(c,x)
    end function

    pure function poly_eval_cp_sp(c,x) result(res)
        complex(sp), intent(in) :: c(:)
        complex(sp), intent(in) :: x
        complex(sp) :: res

        res = horner_cp_sp(c,x)
    end function

    pure function poly_eval_cp_dp(c,x) result(res)
        complex(dp), intent(in) :: c(:)
        complex(dp), intent(in) :: x
        complex(dp) :: res

        res = horner_cp_dp(c,x)
    end function


    pure function rat_eval_re_sp(p,q,x) result(res)
        real(sp), intent(in) :: p(:)
        real(sp), intent(in) :: q(:)
        real(sp), intent(in) :: x
        real(sp) :: res

        res = poly_eval_re_sp(p,x)/poly_eval_re_sp(q,x)
    end function

    pure function rat_eval_re_dp(p,q,x) result(res)
        real(dp), intent(in) :: p(:)
        real(dp), intent(in) :: q(:)
        real(dp), intent(in) :: x
        real(dp) :: res

        res = poly_eval_re_dp(p,x)/poly_eval_re_dp(q,x)
    end function

    pure function rat_eval_cp_sp(p,q,x) result(res)
        complex(sp), intent(in) :: p(:)
        complex(sp), intent(in) :: q(:)
        complex(sp), intent(in) :: x
        complex(sp) :: res

        res = poly_eval_cp_sp(p,x)/poly_eval_cp_sp(q,x)
    end function

    pure function rat_eval_cp_dp(p,q,x) result(res)
        complex(dp), intent(in) :: p(:)
        complex(dp), intent(in) :: q(:)
        complex(dp), intent(in) :: x
        complex(dp) :: res

        res = poly_eval_cp_dp(p,x)/poly_eval_cp_dp(q,x)
    end function


    ! Use Horner's method to evaluate a polynomial.
    pure function horner_re_sp(c,x) result(res)
        real(sp), intent(in) :: c(:)
        real(sp), intent(in) :: x
        real(sp) :: res

        integer :: i

        res = c(size(c))

        do i = size(c)-1,1,-1
            res = c(i) + x*res
        end do
    end function

    pure function horner_re_dp(c,x) result(res)
        real(dp), intent(in) :: c(:)
        real(dp), intent(in) :: x
        real(dp) :: res

        integer :: i

        res = c(size(c))

        do i = size(c)-1,1,-1
            res = c(i) + x*res
        end do
    end function

    pure function horner_cp_sp(c,x) result(res)
        complex(sp), intent(in) :: c(:)
        complex(sp), intent(in) :: x
        complex(sp) :: res

        integer :: i

        res = c(size(c))

        do i = size(c)-1,1,-1
            res = c(i) + x*res
        end do
    end function

    pure function horner_cp_dp(c,x) result(res)
        complex(dp), intent(in) :: c(:)
        complex(dp), intent(in) :: x
        complex(dp) :: res

        integer :: i

        res = c(size(c))

        do i = size(c)-1,1,-1
            res = c(i) + x*res
        end do
    end function


end module polynomial_eval