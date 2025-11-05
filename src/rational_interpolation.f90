module rational_interpolation
    use kinds
    use stdlib_math, only: is_close
    use ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
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
        generic :: init => init_re_sp, init_cp_sp, init_re_dp, init_cp_dp
        generic :: eval => eval_vec_re_sp, eval_scalar_re_sp,eval_vec_cp_sp, eval_scalar_cp_sp, eval_vec_re_dp,&
            & eval_scalar_re_dp,eval_vec_cp_dp, eval_scalar_cp_dp

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
    end type thiele_interp

    type, extends(interp_base), public :: MTT_interp
        real(wp), allocatable :: a_r(:)
        complex(wp), allocatable :: a_z(:)
        real(wp), allocatable :: x_p_r(:)
        complex(wp), allocatable :: x_p_z(:)
        real(wp), allocatable :: y_p_r(:)
        complex(wp), allocatable :: y_p_z(:)
        logical :: is_constant
        integer :: t
    contains
        procedure :: init_re_sp => MTT_init_re_sp
        procedure :: init_cp_sp => MTT_init_cp_sp
        procedure :: eval_scalar_re_sp => MTT_eval_scalar_re_sp
        procedure :: eval_vec_re_sp => MTT_eval_vec_re_sp
        procedure :: eval_scalar_cp_sp => MTT_eval_scalar_cp_sp
        procedure :: eval_vec_cp_sp => MTT_eval_vec_cp_sp
        procedure :: init_re_dp => MTT_init_re_dp
        procedure :: init_cp_dp => MTT_init_cp_dp
        procedure :: eval_scalar_re_dp => MTT_eval_scalar_re_dp
        procedure :: eval_vec_re_dp => MTT_eval_vec_re_dp
        procedure :: eval_scalar_cp_dp => MTT_eval_scalar_cp_dp
        procedure :: eval_vec_cp_dp => MTT_eval_vec_cp_dp
    end type MTT_interp

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

        this%x_z = x
        this%y_z = y
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
            res = (x - this%x_z(i-1))/(this%rho_z(i) - this%rho_z(i-2) + res)
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

        this%x_z = x
        this%y_z = y
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
            res = (x - this%x_z(i-1))/(this%rho_z(i) - this%rho_z(i-2) + res)
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

    subroutine MTT_init_re_sp(this,x,y)
        class(MTT_interp(sp,re)), intent(out) :: this
        real(sp), intent(in) :: x(:)
        real(sp), intent(in) :: y(:)

        real(sp), allocatable :: g_j(:),q_i(:),q_i_m(:),q_i_p(:)
        integer :: j,i
        logical :: found_x_j
        real(sp), parameter :: zero = 0.0_sp

        this%x_r = x
        this%y_r = y
        this%N = size(x)

        do i = 1,this%N
            do j = i+1,this%N
                if (is_close(x(i),x(j))) then
                    call interp_distinct_error('MTT')
                end if
            end do
        end do

        allocate(this%x_p_r(this%N),this%y_p_r(this%N),g_j(2:this%N),this%a_r(2:this%N))
        this%x_p_r = x
        this%y_p_r = y
        g_j = y(2:)-y(1)
        this%a_r = 0

        ! Stage (a)
        this%t = 1
        do j = 2,this%N
            if (.not.ieee_is_finite(g_j(j)).or.is_close(g_j(j),zero)) then
                found_x_j = .false.
                do i = j+1,this%N
                    if (ieee_is_finite(g_j(i)).and..not.is_close(g_j(i),zero)) then
                        this%x_p_r([i,j]) = this%x_p_r([j,i])
                        g_j([i,j]) = g_j([j,i])
                        this%y_p_r([i,j]) = this%y_p_r([j,i])
                        found_x_j = .true.
                        exit
                    end if
                end do
                if (.not.found_x_j) then
                    if (all(ieee_is_finite(g_j(j:)))) then
                        ! Stage (b)
                        exit
                    else
                        ! Stage (c)
                        call interp_exist_error('MTT, Stage (c)')
                    end if
                end if
            end if
            this%a_r(j) = g_j(j)/(this%x_p_r(j)-this%x_p_r(j-1))
            if (j<this%N) then
                g_j(j+1:) = this%a_r(j)*(this%x_p_r(j+1:)-this%x_p_r(j-1))/g_j(j+1:) - 1.0_sp
            end if
            this%t = j
        end do

        ! Termination stage
        if (this%t == 1) then
            this%is_constant = .true.
        else
            allocate(q_i(this%N),q_i_m(this%N),q_i_p(this%N))
            q_i_m = 1.0_sp
            q_i = q_i_m + this%a_r(3)*(this%x_r-this%x_p_r(2))
            if (any(is_close(q_i,zero))) call interp_exist_error('MTT, q_2')
            do i = 3,this%t-1
                q_i_p = q_i + this%a_r(i+1)*(this%x_r-this%x_p_r(i))*q_i_m
                if (any(is_close(q_i_p,zero))) call interp_exist_error('MTT, q_i')
                q_i_m = q_i
                q_i = q_i_p
            end do
        end if

    end subroutine MTT_init_re_sp

    function MTT_eval_scalar_re_sp(this,x) result(res)
        class(MTT_interp(sp,re)), intent(in) :: this
        real(sp), intent(in) :: x
        real(sp) :: res

        integer :: i
        real(sp) :: frac

        res = this%y_r(1)

        if (this%t > 1) then

            frac = 1.0_sp

            do i = this%t,3,-1
                frac = 1.0_sp + this%a_r(i)*(x - this%x_p_r(i-1))/frac
            end do

            frac = this%a_r(2)*(x-this%x_p_r(1))/frac
            res = res + frac
        end if

    end function MTT_eval_scalar_re_sp

    function MTT_eval_vec_re_sp(this,x) result(res)
        class(MTT_interp(sp,re)), intent(in) :: this
        real(sp), intent(in) :: x(:)
        real(sp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function MTT_eval_vec_re_sp

    subroutine MTT_init_cp_sp(this,x,y)
        class(MTT_interp(sp,cp)), intent(out) :: this
        complex(sp), intent(in) :: x(:)
        complex(sp), intent(in) :: y(:)

        complex(sp), allocatable :: g_j(:),q_i(:),q_i_m(:),q_i_p(:)
        integer :: j,i
        logical :: found_x_j
        complex(sp), parameter :: zero = (0.0_sp,0.0_sp)

        this%x_z = x
        this%y_z = y
        this%N = size(x)

        do i = 1,this%N
            do j = i+1,this%N
                if (is_close(x(i),x(j))) then
                    call interp_distinct_error('MTT')
                end if
            end do
        end do

        allocate(this%x_p_z(this%N),this%y_p_z(this%N),g_j(2:this%N),this%a_z(2:this%N))
        this%x_p_z = x
        this%y_p_z = y
        g_j = y(2:)-y(1)
        this%a_z = 0

        ! Stage (a)
        this%t = 1
        do j = 2,this%N
            if (.not.is_finite_sp(g_j(j)).or.is_close(g_j(j),zero)) then
                found_x_j = .false.
                do i = j+1,this%N
                    if (is_finite_sp(g_j(i)).and..not.is_close(g_j(i),zero)) then
                        this%x_p_z([i,j]) = this%x_p_z([j,i])
                        g_j([i,j]) = g_j([j,i])
                        this%y_p_z([i,j]) = this%y_p_z([j,i])
                        found_x_j = .true.
                        exit
                    end if
                end do
                if (.not.found_x_j) then
                    if (all(is_finite_sp(g_j(j:)))) then
                        ! Stage (b)
                        exit
                    else
                        ! Stage (c)
                        call interp_exist_error('MTT, Stage (c)')
                    end if
                end if
            end if
            this%a_z(j) = g_j(j)/(this%x_p_z(j)-this%x_p_z(j-1))
            if (j<this%N) then
                g_j(j+1:) = this%a_z(j)*(this%x_p_z(j+1:)-this%x_p_z(j-1))/g_j(j+1:) - 1.0_sp
            end if
            this%t = j
        end do

        ! Termination stage
        if (this%t == 1) then
            this%is_constant = .true.
        else
            allocate(q_i(this%N),q_i_m(this%N),q_i_p(this%N))
            q_i_m = 1.0_sp
            q_i = q_i_m + this%a_z(3)*(this%x_z-this%x_p_z(2))
            if (any(is_close(q_i,zero))) call interp_exist_error('MTT, q_2')
            do i = 3,this%t-1
                q_i_p = q_i + this%a_z(i+1)*(this%x_z-this%x_p_z(i))*q_i_m
                if (any(is_close(q_i_p,zero))) call interp_exist_error('MTT, q_i')
                q_i_m = q_i
                q_i = q_i_p
            end do
        end if
    end subroutine MTT_init_cp_sp

    function MTT_eval_scalar_cp_sp(this,x) result(res)
        class(MTT_interp(sp,cp)), intent(in) :: this
        complex(sp), intent(in) :: x
        complex(sp) :: res

        integer :: i
        complex(sp) :: frac

        res = this%y_z(1)

        if (this%t > 1) then

            frac = 1.0_sp

            do i = this%t,3,-1
                frac = 1.0_sp + this%a_z(i)*(x - this%x_p_z(i-1))/frac
            end do

            frac = this%a_z(2)*(x-this%x_p_z(1))/frac
            res = res + frac
        end if

    end function MTT_eval_scalar_cp_sp

    function MTT_eval_vec_cp_sp(this,x) result(res)
        class(MTT_interp(sp,cp)), intent(in) :: this
        complex(sp), intent(in) :: x(:)
        complex(sp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function MTT_eval_vec_cp_sp

    elemental function is_finite_sp(z) result(res)
        complex(sp), intent(in) :: z
        logical :: res

        res = ieee_is_finite(real(z,kind=sp)).and.ieee_is_finite(aimag(z))
    end function is_finite_sp

    subroutine MTT_init_re_dp(this,x,y)
        class(MTT_interp(dp,re)), intent(out) :: this
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:)

        real(dp), allocatable :: g_j(:),q_i(:),q_i_m(:),q_i_p(:)
        integer :: j,i
        logical :: found_x_j
        real(dp), parameter :: zero = 0.0_dp

        this%x_r = x
        this%y_r = y
        this%N = size(x)

        do i = 1,this%N
            do j = i+1,this%N
                if (is_close(x(i),x(j))) then
                    call interp_distinct_error('MTT')
                end if
            end do
        end do

        allocate(this%x_p_r(this%N),this%y_p_r(this%N),g_j(2:this%N),this%a_r(2:this%N))
        this%x_p_r = x
        this%y_p_r = y
        g_j = y(2:)-y(1)
        this%a_r = 0

        ! Stage (a)
        this%t = 1
        do j = 2,this%N
            if (.not.ieee_is_finite(g_j(j)).or.is_close(g_j(j),zero)) then
                found_x_j = .false.
                do i = j+1,this%N
                    if (ieee_is_finite(g_j(i)).and..not.is_close(g_j(i),zero)) then
                        this%x_p_r([i,j]) = this%x_p_r([j,i])
                        g_j([i,j]) = g_j([j,i])
                        this%y_p_r([i,j]) = this%y_p_r([j,i])
                        found_x_j = .true.
                        exit
                    end if
                end do
                if (.not.found_x_j) then
                    if (all(ieee_is_finite(g_j(j:)))) then
                        ! Stage (b)
                        exit
                    else
                        ! Stage (c)
                        call interp_exist_error('MTT, Stage (c)')
                    end if
                end if
            end if
            this%a_r(j) = g_j(j)/(this%x_p_r(j)-this%x_p_r(j-1))
            if (j<this%N) then
                g_j(j+1:) = this%a_r(j)*(this%x_p_r(j+1:)-this%x_p_r(j-1))/g_j(j+1:) - 1.0_dp
            end if
            this%t = j
        end do

        ! Termination stage
        if (this%t == 1) then
            this%is_constant = .true.
        else
            allocate(q_i(this%N),q_i_m(this%N),q_i_p(this%N))
            q_i_m = 1.0_dp
            q_i = q_i_m + this%a_r(3)*(this%x_r-this%x_p_r(2))
            if (any(is_close(q_i,zero))) call interp_exist_error('MTT, q_2')
            do i = 3,this%t-1
                q_i_p = q_i + this%a_r(i+1)*(this%x_r-this%x_p_r(i))*q_i_m
                if (any(is_close(q_i_p,zero))) call interp_exist_error('MTT, q_i')
                q_i_m = q_i
                q_i = q_i_p
            end do
        end if

    end subroutine MTT_init_re_dp

    function MTT_eval_scalar_re_dp(this,x) result(res)
        class(MTT_interp(dp,re)), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: res

        integer :: i
        real(dp) :: frac

        res = this%y_r(1)

        if (this%t > 1) then

            frac = 1.0_dp

            do i = this%t,3,-1
                frac = 1.0_dp + this%a_r(i)*(x - this%x_p_r(i-1))/frac
            end do

            frac = this%a_r(2)*(x-this%x_p_r(1))/frac
            res = res + frac
        end if

    end function MTT_eval_scalar_re_dp

    function MTT_eval_vec_re_dp(this,x) result(res)
        class(MTT_interp(dp,re)), intent(in) :: this
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function MTT_eval_vec_re_dp

    subroutine MTT_init_cp_dp(this,x,y)
        class(MTT_interp(dp,cp)), intent(out) :: this
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(in) :: y(:)

        complex(dp), allocatable :: g_j(:),q_i(:),q_i_m(:),q_i_p(:)
        integer :: j,i
        logical :: found_x_j
        complex(dp), parameter :: zero = (0.0_dp,0.0_dp)

        this%x_z = x
        this%y_z = y
        this%N = size(x)

        do i = 1,this%N
            do j = i+1,this%N
                if (is_close(x(i),x(j))) then
                    call interp_distinct_error('MTT')
                end if
            end do
        end do

        allocate(this%x_p_z(this%N),this%y_p_z(this%N),g_j(2:this%N),this%a_z(2:this%N))
        this%x_p_z = x
        this%y_p_z = y
        g_j = y(2:)-y(1)
        this%a_z = 0

        ! Stage (a)
        this%t = 1
        do j = 2,this%N
            if (.not.is_finite_dp(g_j(j)).or.is_close(g_j(j),zero)) then
                found_x_j = .false.
                do i = j+1,this%N
                    if (is_finite_dp(g_j(i)).and..not.is_close(g_j(i),zero)) then
                        this%x_p_z([i,j]) = this%x_p_z([j,i])
                        g_j([i,j]) = g_j([j,i])
                        this%y_p_z([i,j]) = this%y_p_z([j,i])
                        found_x_j = .true.
                        exit
                    end if
                end do
                if (.not.found_x_j) then
                    if (all(is_finite_dp(g_j(j:)))) then
                        ! Stage (b)
                        exit
                    else
                        ! Stage (c)
                        call interp_exist_error('MTT, Stage (c)')
                    end if
                end if
            end if
            this%a_z(j) = g_j(j)/(this%x_p_z(j)-this%x_p_z(j-1))
            if (j<this%N) then
                g_j(j+1:) = this%a_z(j)*(this%x_p_z(j+1:)-this%x_p_z(j-1))/g_j(j+1:) - 1.0_dp
            end if
            this%t = j
        end do

        ! Termination stage
        if (this%t == 1) then
            this%is_constant = .true.
        else
            allocate(q_i(this%N),q_i_m(this%N),q_i_p(this%N))
            q_i_m = 1.0_dp
            q_i = q_i_m + this%a_z(3)*(this%x_z-this%x_p_z(2))
            if (any(is_close(q_i,zero))) call interp_exist_error('MTT, q_2')
            do i = 3,this%t-1
                q_i_p = q_i + this%a_z(i+1)*(this%x_z-this%x_p_z(i))*q_i_m
                if (any(is_close(q_i_p,zero))) call interp_exist_error('MTT, q_i')
                q_i_m = q_i
                q_i = q_i_p
            end do
        end if
    end subroutine MTT_init_cp_dp

    function MTT_eval_scalar_cp_dp(this,x) result(res)
        class(MTT_interp(dp,cp)), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: res

        integer :: i
        complex(dp) :: frac

        res = this%y_z(1)

        if (this%t > 1) then

            frac = 1.0_dp

            do i = this%t,3,-1
                frac = 1.0_dp + this%a_z(i)*(x - this%x_p_z(i-1))/frac
            end do

            frac = this%a_z(2)*(x-this%x_p_z(1))/frac
            res = res + frac
        end if

    end function MTT_eval_scalar_cp_dp

    function MTT_eval_vec_cp_dp(this,x) result(res)
        class(MTT_interp(dp,cp)), intent(in) :: this
        complex(dp), intent(in) :: x(:)
        complex(dp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        res = 0

        do i = 1,size(x)
            res(i) = this%eval(x(i))
        end do

    end function MTT_eval_vec_cp_dp

    elemental function is_finite_dp(z) result(res)
        complex(dp), intent(in) :: z
        logical :: res

        res = ieee_is_finite(real(z,kind=dp)).and.ieee_is_finite(aimag(z))
    end function is_finite_dp


    subroutine interp_exist_error(error_message)
        character(len=*), intent(in) :: error_message

        write(stderr,'(a,a)') 'Error! Interpolant does not exist: ', error_message
        error stop
    end subroutine interp_exist_error

    subroutine interp_distinct_error(error_message)
        character(len=*), intent(in) :: error_message

        write(stderr,'(a,a)') 'Error! Interpolation points are not distinct in: ', error_message
        error stop
    end subroutine interp_distinct_error
end module rational_interpolation