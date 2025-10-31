module rational_interpolation
    use kinds
    implicit none

    integer, parameter :: re = 1
    integer, parameter :: cp = 2

    type, abstract, public :: interp_base(wp,domain)
        integer, kind :: wp
        integer, kind :: domain
        integer :: N
        real(wp), allocatable :: x_r(:)
        real(wp), allocatable :: y_r(:)
        complex(wp), allocatable :: x_z(:)
        complex(wp), allocatable :: y_z(:)
    contains
        generic :: init => init_re_sp, init_cp_sp, init_re_dp, init_cp_dp, init_re_qp, init_cp_qp
        generic :: eval => eval_vec_re_sp, eval_scalar_re_sp,eval_vec_cp_sp, eval_scalar_cp_sp, eval_vec_re_dp,&
            & eval_scalar_re_dp,eval_vec_cp_dp, eval_scalar_cp_dp, eval_vec_re_qp, eval_scalar_re_qp,eval_vec_cp_qp,&
            & eval_scalar_cp_qp

        procedure(re_init_sp), deferred :: init_re_sp
        procedure(cp_init_sp), deferred :: init_cp_sp
        procedure(re_eval_scalar_sp), private, deferred :: eval_scalar_re_sp
        procedure(re_eval_vec_sp), private, deferred :: eval_vec_re_sp
        procedure(cp_eval_scalar_sp), private, deferred :: eval_scalar_cp_sp
        procedure(cp_eval_vec_sp), private, deferred :: eval_vec_cp_sp
        procedure(re_init_dp), deferred :: init_re_dp
        procedure(cp_init_dp), deferred :: init_cp_dp
        procedure(re_eval_scalar_dp), private, deferred :: eval_scalar_re_dp
        procedure(re_eval_vec_dp), private, deferred :: eval_vec_re_dp
        procedure(cp_eval_scalar_dp), private, deferred :: eval_scalar_cp_dp
        procedure(cp_eval_vec_dp), private, deferred :: eval_vec_cp_dp
        procedure(re_init_qp), deferred :: init_re_qp
        procedure(cp_init_qp), deferred :: init_cp_qp
        procedure(re_eval_scalar_qp), private, deferred :: eval_scalar_re_qp
        procedure(re_eval_vec_qp), private, deferred :: eval_vec_re_qp
        procedure(cp_eval_scalar_qp), private, deferred :: eval_scalar_cp_qp
        procedure(cp_eval_vec_qp), private, deferred :: eval_vec_cp_qp
    end type interp_base

    type, extends(interp_base), public :: thiele_interp
        real(wp), allocatable :: rho_r(:)
        complex(wp), allocatable :: rho_z(:)
    contains
        procedure :: init_re_sp => thiele_init_re_sp
        procedure :: init_cp_sp => thiele_init_cp_sp
        procedure :: eval_scalar_re_sp => thiele_eval_scalar_re_sp
        procedure :: eval_vec_re_sp => thiele_eval_vec_re_sp
        procedure :: eval_scalar_cp_sp => thiele_eval_scalar_cp_sp
        procedure :: eval_vec_cp_sp => thiele_eval_vec_cp_sp
        procedure :: init_re_dp => thiele_init_re_dp
        procedure :: init_cp_dp => thiele_init_cp_dp
        procedure :: eval_scalar_re_dp => thiele_eval_scalar_re_dp
        procedure :: eval_vec_re_dp => thiele_eval_vec_re_dp
        procedure :: eval_scalar_cp_dp => thiele_eval_scalar_cp_dp
        procedure :: eval_vec_cp_dp => thiele_eval_vec_cp_dp
        procedure :: init_re_qp => thiele_init_re_qp
        procedure :: init_cp_qp => thiele_init_cp_qp
        procedure :: eval_scalar_re_qp => thiele_eval_scalar_re_qp
        procedure :: eval_vec_re_qp => thiele_eval_vec_re_qp
        procedure :: eval_scalar_cp_qp => thiele_eval_scalar_cp_qp
        procedure :: eval_vec_cp_qp => thiele_eval_vec_cp_qp
    end type thiele_interp

    abstract interface
        subroutine re_init_sp(this,x,y)
            import :: sp
            import :: re
            import :: interp_base
            class(interp_base(sp,re)), intent(out) :: this
            real(sp), intent(in) :: x(:)
            real(sp), intent(in) :: y(:)
        end subroutine re_init_sp

        subroutine cp_init_sp(this,x,y)
            import :: sp
            import :: cp
            import :: interp_base
            class(interp_base(sp,cp)), intent(out) :: this
            complex(sp), intent(in) :: x(:)
            complex(sp), intent(in) :: y(:)
        end subroutine cp_init_sp

        function re_eval_vec_sp(this,x) result(res)
            import :: sp
            import :: re
            import :: interp_base
            class(interp_base(sp,re)), intent(in) :: this
            real(sp), intent(in) :: x(:)
            real(sp), allocatable :: res(:)
        end function re_eval_vec_sp

        function cp_eval_vec_sp(this,x) result(res)
            import :: sp
            import :: cp
            import :: interp_base
            class(interp_base(sp,cp)), intent(in) :: this
            complex(sp), intent(in) :: x(:)
            complex(sp), allocatable :: res(:)
        end function cp_eval_vec_sp

        function re_eval_scalar_sp(this,x) result(res)
            import :: sp
            import :: re
            import :: interp_base
            class(interp_base(sp,re)), intent(in) :: this
            real(sp), intent(in) :: x
            real(sp) :: res
        end function re_eval_scalar_sp

        function cp_eval_scalar_sp(this,x) result(res)
            import :: sp
            import :: cp
            import :: interp_base
            class(interp_base(sp,cp)), intent(in) :: this
            complex(sp), intent(in) :: x
            complex(sp) :: res
        end function cp_eval_scalar_sp

        subroutine re_init_dp(this,x,y)
            import :: dp
            import :: re
            import :: interp_base
            class(interp_base(dp,re)), intent(out) :: this
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: y(:)
        end subroutine re_init_dp

        subroutine cp_init_dp(this,x,y)
            import :: dp
            import :: cp
            import :: interp_base
            class(interp_base(dp,cp)), intent(out) :: this
            complex(dp), intent(in) :: x(:)
            complex(dp), intent(in) :: y(:)
        end subroutine cp_init_dp

        function re_eval_vec_dp(this,x) result(res)
            import :: dp
            import :: re
            import :: interp_base
            class(interp_base(dp,re)), intent(in) :: this
            real(dp), intent(in) :: x(:)
            real(dp), allocatable :: res(:)
        end function re_eval_vec_dp

        function cp_eval_vec_dp(this,x) result(res)
            import :: dp
            import :: cp
            import :: interp_base
            class(interp_base(dp,cp)), intent(in) :: this
            complex(dp), intent(in) :: x(:)
            complex(dp), allocatable :: res(:)
        end function cp_eval_vec_dp

        function re_eval_scalar_dp(this,x) result(res)
            import :: dp
            import :: re
            import :: interp_base
            class(interp_base(dp,re)), intent(in) :: this
            real(dp), intent(in) :: x
            real(dp) :: res
        end function re_eval_scalar_dp

        function cp_eval_scalar_dp(this,x) result(res)
            import :: dp
            import :: cp
            import :: interp_base
            class(interp_base(dp,cp)), intent(in) :: this
            complex(dp), intent(in) :: x
            complex(dp) :: res
        end function cp_eval_scalar_dp

        subroutine re_init_qp(this,x,y)
            import :: qp
            import :: re
            import :: interp_base
            class(interp_base(qp,re)), intent(out) :: this
            real(qp), intent(in) :: x(:)
            real(qp), intent(in) :: y(:)
        end subroutine re_init_qp

        subroutine cp_init_qp(this,x,y)
            import :: qp
            import :: cp
            import :: interp_base
            class(interp_base(qp,cp)), intent(out) :: this
            complex(qp), intent(in) :: x(:)
            complex(qp), intent(in) :: y(:)
        end subroutine cp_init_qp

        function re_eval_vec_qp(this,x) result(res)
            import :: qp
            import :: re
            import :: interp_base
            class(interp_base(qp,re)), intent(in) :: this
            real(qp), intent(in) :: x(:)
            real(qp), allocatable :: res(:)
        end function re_eval_vec_qp

        function cp_eval_vec_qp(this,x) result(res)
            import :: qp
            import :: cp
            import :: interp_base
            class(interp_base(qp,cp)), intent(in) :: this
            complex(qp), intent(in) :: x(:)
            complex(qp), allocatable :: res(:)
        end function cp_eval_vec_qp

        function re_eval_scalar_qp(this,x) result(res)
            import :: qp
            import :: re
            import :: interp_base
            class(interp_base(qp,re)), intent(in) :: this
            real(qp), intent(in) :: x
            real(qp) :: res
        end function re_eval_scalar_qp

        function cp_eval_scalar_qp(this,x) result(res)
            import :: qp
            import :: cp
            import :: interp_base
            class(interp_base(qp,cp)), intent(in) :: this
            complex(qp), intent(in) :: x
            complex(qp) :: res
        end function cp_eval_scalar_qp

    end interface

contains
    subroutine thiele_init_re_sp(this,x,y)
        class(thiele_interp(sp,re)), intent(out) :: this
        real(sp), intent(in) :: x(:)
        real(sp), intent(in) :: y(:)

        real(sp), allocatable :: rho_temp(:,:)
        integer :: i

        this%x_r = x
        this%y_r = y
        this%N = size(x)

        allocate(rho_temp(this%N,this%N))
        rho_temp = 0
        rho_temp(:,1) = this%y_r
        rho_temp(:this%N-1,2) = (this%x_r(:this%N-1) - this%x_r(2:))/(rho_temp(:this%N-1,1) - rho_temp(2:,1))

        do i = 3,this%N
            rho_temp(:this%N-i+1,i) = (this%x_r(:this%N-i+1) - this%x_r(i:))&
                                        /(rho_temp(:this%N-i+1,i-1) - rho_temp(2:this%N-i+2,i-1)) &
                                        + rho_temp(2:this%N-i+2,i-2)
        end do

        this%rho_r = rho_temp(1,:)
    end subroutine thiele_init_re_sp

    function thiele_eval_scalar_re_sp(this,x) result(res)
        class(thiele_interp(sp,re)), intent(in) :: this
        real(sp), intent(in) :: x
        real(sp) :: res

        integer :: i
        res = 0

        do i = this%N,3,-1
            res = (x - this%x_r(i-1))/(this%rho_r(i) - this%rho_r(i-2) + res)
        end do

        res = this%rho_r(1) + (x - this%x_r(1))/(this%rho_r(2) + res)

    end function thiele_eval_scalar_re_sp

    function thiele_eval_vec_re_sp(this,x) result(res)
        class(thiele_interp(sp,re)), intent(in) :: this
        real(sp), intent(in) :: x(:)
        real(sp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function thiele_eval_vec_re_sp

subroutine thiele_init_cp_sp(this,x,y)
        class(thiele_interp(sp,cp)), intent(out) :: this
        complex(sp), intent(in) :: x(:)
        complex(sp), intent(in) :: y(:)

        complex(sp), allocatable :: rho_temp(:,:)
        integer :: i

        this%x_r = x
        this%y_r = y
        this%N = size(x)

        allocate(rho_temp(this%N,this%N))
        rho_temp = 0
        rho_temp(:,1) = this%y_z
        rho_temp(:this%N-1,2) = (this%x_z(:this%N-1) - this%x_z(2:))/(rho_temp(:this%N-1,1) - rho_temp(2:,1))

        do i = 3,this%N
            rho_temp(:this%N-i+1,i) = (this%x_z(:this%N-i+1) - this%x_z(i:))&
                                        /(rho_temp(:this%N-i+1,i-1) - rho_temp(2:this%N-i+2,i-1)) &
                                        + rho_temp(2:this%N-i+2,i-2)
        end do

        this%rho_z = rho_temp(1,:)
    end subroutine thiele_init_cp_sp

    function thiele_eval_scalar_cp_sp(this,x) result(res)
        class(thiele_interp(sp,cp)), intent(in) :: this
        complex(sp), intent(in) :: x
        complex(sp) :: res

        integer :: i
        res = 0

        do i = this%N,3,-1
            res = (x - this%x_r(i-1))/(this%rho_z(i) - this%rho_z(i-2) + res)
        end do

        res = this%rho_z(1) + (x - this%x_z(1))/(this%rho_z(2) + res)

    end function thiele_eval_scalar_cp_sp

    function thiele_eval_vec_cp_sp(this,x) result(res)
        class(thiele_interp(sp,cp)), intent(in) :: this
        complex(sp), intent(in) :: x(:)
        complex(sp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function thiele_eval_vec_cp_sp

    subroutine thiele_init_re_dp(this,x,y)
        class(thiele_interp(dp,re)), intent(out) :: this
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:)

        real(dp), allocatable :: rho_temp(:,:)
        integer :: i

        this%x_r = x
        this%y_r = y
        this%N = size(x)

        allocate(rho_temp(this%N,this%N))
        rho_temp = 0
        rho_temp(:,1) = this%y_r
        rho_temp(:this%N-1,2) = (this%x_r(:this%N-1) - this%x_r(2:))/(rho_temp(:this%N-1,1) - rho_temp(2:,1))

        do i = 3,this%N
            rho_temp(:this%N-i+1,i) = (this%x_r(:this%N-i+1) - this%x_r(i:))&
                                        /(rho_temp(:this%N-i+1,i-1) - rho_temp(2:this%N-i+2,i-1)) &
                                        + rho_temp(2:this%N-i+2,i-2)
        end do

        this%rho_r = rho_temp(1,:)
    end subroutine thiele_init_re_dp

    function thiele_eval_scalar_re_dp(this,x) result(res)
        class(thiele_interp(dp,re)), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: res

        integer :: i
        res = 0

        do i = this%N,3,-1
            res = (x - this%x_r(i-1))/(this%rho_r(i) - this%rho_r(i-2) + res)
        end do

        res = this%rho_r(1) + (x - this%x_r(1))/(this%rho_r(2) + res)

    end function thiele_eval_scalar_re_dp

    function thiele_eval_vec_re_dp(this,x) result(res)
        class(thiele_interp(dp,re)), intent(in) :: this
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function thiele_eval_vec_re_dp

subroutine thiele_init_cp_dp(this,x,y)
        class(thiele_interp(dp,cp)), intent(out) :: this
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(in) :: y(:)

        complex(dp), allocatable :: rho_temp(:,:)
        integer :: i

        this%x_r = x
        this%y_r = y
        this%N = size(x)

        allocate(rho_temp(this%N,this%N))
        rho_temp = 0
        rho_temp(:,1) = this%y_z
        rho_temp(:this%N-1,2) = (this%x_z(:this%N-1) - this%x_z(2:))/(rho_temp(:this%N-1,1) - rho_temp(2:,1))

        do i = 3,this%N
            rho_temp(:this%N-i+1,i) = (this%x_z(:this%N-i+1) - this%x_z(i:))&
                                        /(rho_temp(:this%N-i+1,i-1) - rho_temp(2:this%N-i+2,i-1)) &
                                        + rho_temp(2:this%N-i+2,i-2)
        end do

        this%rho_z = rho_temp(1,:)
    end subroutine thiele_init_cp_dp

    function thiele_eval_scalar_cp_dp(this,x) result(res)
        class(thiele_interp(dp,cp)), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: res

        integer :: i
        res = 0

        do i = this%N,3,-1
            res = (x - this%x_r(i-1))/(this%rho_z(i) - this%rho_z(i-2) + res)
        end do

        res = this%rho_z(1) + (x - this%x_z(1))/(this%rho_z(2) + res)

    end function thiele_eval_scalar_cp_dp

    function thiele_eval_vec_cp_dp(this,x) result(res)
        class(thiele_interp(dp,cp)), intent(in) :: this
        complex(dp), intent(in) :: x(:)
        complex(dp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function thiele_eval_vec_cp_dp

    subroutine thiele_init_re_qp(this,x,y)
        class(thiele_interp(qp,re)), intent(out) :: this
        real(qp), intent(in) :: x(:)
        real(qp), intent(in) :: y(:)

        real(qp), allocatable :: rho_temp(:,:)
        integer :: i

        this%x_r = x
        this%y_r = y
        this%N = size(x)

        allocate(rho_temp(this%N,this%N))
        rho_temp = 0
        rho_temp(:,1) = this%y_r
        rho_temp(:this%N-1,2) = (this%x_r(:this%N-1) - this%x_r(2:))/(rho_temp(:this%N-1,1) - rho_temp(2:,1))

        do i = 3,this%N
            rho_temp(:this%N-i+1,i) = (this%x_r(:this%N-i+1) - this%x_r(i:))&
                                        /(rho_temp(:this%N-i+1,i-1) - rho_temp(2:this%N-i+2,i-1)) &
                                        + rho_temp(2:this%N-i+2,i-2)
        end do

        this%rho_r = rho_temp(1,:)
    end subroutine thiele_init_re_qp

    function thiele_eval_scalar_re_qp(this,x) result(res)
        class(thiele_interp(qp,re)), intent(in) :: this
        real(qp), intent(in) :: x
        real(qp) :: res

        integer :: i
        res = 0

        do i = this%N,3,-1
            res = (x - this%x_r(i-1))/(this%rho_r(i) - this%rho_r(i-2) + res)
        end do

        res = this%rho_r(1) + (x - this%x_r(1))/(this%rho_r(2) + res)

    end function thiele_eval_scalar_re_qp

    function thiele_eval_vec_re_qp(this,x) result(res)
        class(thiele_interp(qp,re)), intent(in) :: this
        real(qp), intent(in) :: x(:)
        real(qp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function thiele_eval_vec_re_qp

subroutine thiele_init_cp_qp(this,x,y)
        class(thiele_interp(qp,cp)), intent(out) :: this
        complex(qp), intent(in) :: x(:)
        complex(qp), intent(in) :: y(:)

        complex(qp), allocatable :: rho_temp(:,:)
        integer :: i

        this%x_r = x
        this%y_r = y
        this%N = size(x)

        allocate(rho_temp(this%N,this%N))
        rho_temp = 0
        rho_temp(:,1) = this%y_z
        rho_temp(:this%N-1,2) = (this%x_z(:this%N-1) - this%x_z(2:))/(rho_temp(:this%N-1,1) - rho_temp(2:,1))

        do i = 3,this%N
            rho_temp(:this%N-i+1,i) = (this%x_z(:this%N-i+1) - this%x_z(i:))&
                                        /(rho_temp(:this%N-i+1,i-1) - rho_temp(2:this%N-i+2,i-1)) &
                                        + rho_temp(2:this%N-i+2,i-2)
        end do

        this%rho_z = rho_temp(1,:)
    end subroutine thiele_init_cp_qp

    function thiele_eval_scalar_cp_qp(this,x) result(res)
        class(thiele_interp(qp,cp)), intent(in) :: this
        complex(qp), intent(in) :: x
        complex(qp) :: res

        integer :: i
        res = 0

        do i = this%N,3,-1
            res = (x - this%x_r(i-1))/(this%rho_z(i) - this%rho_z(i-2) + res)
        end do

        res = this%rho_z(1) + (x - this%x_z(1))/(this%rho_z(2) + res)

    end function thiele_eval_scalar_cp_qp

    function thiele_eval_vec_cp_qp(this,x) result(res)
        class(thiele_interp(qp,cp)), intent(in) :: this
        complex(qp), intent(in) :: x(:)
        complex(qp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function thiele_eval_vec_cp_qp

end module rational_interpolation