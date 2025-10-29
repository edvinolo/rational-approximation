module rational_interpolation
    use kinds
    implicit none

    type, abstract, public :: real_interp_base(wp)
        integer, kind :: wp = dp
        real(wp), allocatable :: x(:)
        real(wp), allocatable :: y(:)
    contains
        procedure(real_init), deferred :: init
        generic :: eval => eval_vec, eval_scalar
        procedure(real_eval_scalar), private, deferred :: eval_scalar
        procedure(real_eval_vec), private, deferred :: eval_vec
    end type real_interp_base

    type, abstract, public :: cp_interp_base(wp)
        integer, kind :: wp = dp
        complex(wp), allocatable :: x(:)
        complex(wp), allocatable :: y(:)
    contains
        procedure(cp_init), deferred :: init
        generic :: eval => eval_vec, eval_scalar
        procedure(cp_eval_scalar), private, deferred :: eval_scalar
        procedure(cp_eval_vec), private, deferred :: eval_vec
    end type cp_interp_base

    type, extends(real_interp_base), public :: real_thiele_interp
    contains
        procedure :: init => real_thiele_init
        procedure :: eval_scalar => real_thiele_eval_scalar
        procedure :: eval_vec => real_thiele_eval_vec
    end type real_thiele_interp

    type, extends(cp_interp_base), public :: cp_thiele_interp
    contains
        procedure :: init => cp_thiele_init
        procedure :: eval_scalar => cp_thiele_eval_scalar
        procedure :: eval_vec => cp_thiele_eval_vec
    end type cp_thiele_interp

    abstract interface
        subroutine real_init(this,x,y)
            import :: real_interp_base
            class(real_interp_base), intent(out) :: this
            real(this%wp), intent(in) :: x(:)
            real(this%wp), intent(in) :: y(:)
        end subroutine real_init

        subroutine cp_init(this,x,y)
            import :: cp_interp_base
            class(cp_interp_base), intent(out) :: this
            complex(this%wp), intent(in) :: x(:)
            complex(this%wp), intent(in) :: y(:)
        end subroutine cp_init

        function real_eval_vec(this,x) result(res)
            import :: real_interp_base
            class(real_interp_base), intent(in) :: this
            real(this%wp), intent(in) :: x(:)
            real(this%wp), allocatable :: res(:)
        end function real_eval_vec

        function cp_eval_vec(this,x) result(res)
            import :: dp
            import :: cp_interp_base
            class(cp_interp_base), intent(in) :: this
            complex(this%wp), intent(in) :: x(:)
            complex(this%wp), allocatable :: res(:)
        end function cp_eval_vec

        function real_eval_scalar(this,x) result(res)
            import :: real_interp_base
            class(real_interp_base), intent(in) :: this
            real(this%wp), intent(in) :: x
            real(this%wp) :: res
        end function real_eval_scalar

        function cp_eval_scalar(this,x) result(res)
            import :: dp
            import :: cp_interp_base
            class(cp_interp_base), intent(in) :: this
            complex(this%wp), intent(in) :: x
            complex(this%wp) :: res
        end function cp_eval_scalar
    end interface

contains
    subroutine real_thiele_init(this,x,y)
        class(real_thiele_interp), intent(out) :: this
        real(this%wp), intent(in) :: x(:)
        real(this%wp), intent(in) :: y(:)

        this%x = x
        this%y = y

    end subroutine real_thiele_init

    function real_thiele_eval_scalar(this,x) result(res)
        class(real_thiele_interp), intent(in) :: this
        real(this%wp), intent(in) :: x
        real(this%wp) :: res

        res = 0

    end function real_thiele_eval_scalar

    function real_thiele_eval_vec(this,x) result(res)
        class(real_thiele_interp), intent(in) :: this
        real(this%wp), intent(in) :: x(:)
        real(this%wp), allocatable :: res(:)

        allocate(res(size(x)))
        res = 0

    end function real_thiele_eval_vec

    subroutine cp_thiele_init(this,x,y)
        class(cp_thiele_interp), intent(out) :: this
        complex(this%wp), intent(in) :: x(:)
        complex(this%wp), intent(in) :: y(:)

        this%x = x
        this%y = y

    end subroutine cp_thiele_init

    function cp_thiele_eval_scalar(this,x) result(res)
        class(cp_thiele_interp), intent(in) :: this
        complex(this%wp), intent(in) :: x
        complex(this%wp) :: res

        res = 0

    end function cp_thiele_eval_scalar

    function cp_thiele_eval_vec(this,x) result(res)
        class(cp_thiele_interp), intent(in) :: this
        complex(this%wp), intent(in) :: x(:)
        complex(this%wp), allocatable :: res(:)

        allocate(res(size(x)))
        res = 0

    end function cp_thiele_eval_vec
end module rational_interpolation
