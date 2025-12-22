module robust_pade
    use kinds
    use polynomial_eval, only: rat_eval
    use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
    use stdlib_optval, only: optval
    use stdlib_linalg, only: svd
    use stdlib_linalg_blas, only: gemv
    implicit none

        real(sp), parameter :: DEFAULT_TOL_sp = 1.0e2_sp*epsilon(1.0_sp)
        real(dp), parameter :: DEFAULT_TOL_dp = 1.0e2_dp*epsilon(1.0_dp)

    type, public :: pade_re(wp)
        integer, kind :: wp
        integer :: m
        integer :: n
        integer :: mu
        integer :: nu
        real(wp) :: tol
        real(wp), allocatable :: c(:)
        real(wp), allocatable :: p(:)
        real(wp), allocatable :: q(:)

    contains
        generic :: init => init_re_sp, init_re_dp
        generic :: eval => eval_vec_re_sp, eval_scalar_re_sp, eval_vec_re_dp, eval_scalar_re_dp

        procedure :: init_re_sp => pade_init_re_sp
        procedure :: eval_scalar_re_sp => pade_eval_scalar_re_sp
        procedure :: eval_vec_re_sp => pade_eval_vec_re_sp
        procedure :: init_re_dp => pade_init_re_dp
        procedure :: eval_scalar_re_dp => pade_eval_scalar_re_dp
        procedure :: eval_vec_re_dp => pade_eval_vec_re_dp
    end type pade_re

    type, public :: pade_cp(wp)
        integer, kind :: wp
        integer :: m
        integer :: n
        integer :: mu
        integer :: nu
        real(wp) :: tol
        complex(wp), allocatable :: c(:)
        complex(wp), allocatable :: p(:)
        complex(wp), allocatable :: q(:)

    contains
        generic :: init => init_cp_sp, init_cp_dp
        generic :: eval => eval_vec_cp_sp, eval_scalar_cp_sp, eval_vec_cp_dp, eval_scalar_cp_dp

        procedure :: init_cp_sp => pade_init_cp_sp
        procedure :: eval_scalar_cp_sp => pade_eval_scalar_cp_sp
        procedure :: eval_vec_cp_sp => pade_eval_vec_cp_sp
        procedure :: init_cp_dp => pade_init_cp_dp
        procedure :: eval_scalar_cp_dp => pade_eval_scalar_cp_dp
        procedure :: eval_vec_cp_dp => pade_eval_vec_cp_dp
    end type pade_cp


contains

    ! Try to compute a Padé approximant of type [m/n] for a function with Taylor coeffs. c.
    ! The SVD-based robust Padé algorithm of Gonnet, Güttel, and Trefethen (SIAM Review 55, 1, 101-117, 2013),
    ! is used. This means that the final approximant could be of lower order [mu,nu], mu <= m, nu <= n,
    ! depending on noise and tolerance levels.
    subroutine pade_init_re_sp(this,m,n,c,tol)
        class(pade_re(sp)), intent(out) :: this
        integer, intent(in) :: m
        integer, intent(in) :: n
        real(sp), intent(in) :: c(:)
        real(sp), intent(in), optional :: tol

        real(sp) :: tau
        real(sp), allocatable :: G_upper(:,:)
        real(sp), allocatable :: G_tilde(:,:)
        real(sp), allocatable :: G(:,:)
        real(sp), allocatable :: s(:)
        real(sp), allocatable :: U(:,:)
        real(sp), allocatable :: V_t(:,:)
        integer :: rho
        integer :: i,j
        real(sp) :: one,zero

        if (m + n + 1 /= size(c)) then
            write(stderr,*) 'Error in pade_init_re_sp, len(c) not equal to m + n + 1.'
            error stop
        end if

        this%tol = optval(tol, DEFAULT_TOL_sp)
        this%c = c
        this%m = m
        this%n = n
        this%mu = m
        this%nu = n

        one = 1.0_sp
        zero = 0.0_sp
        allocate(this%p(m+1),this%q(n+1),source = zero)
        allocate(G_tilde(n,n+1),source = zero)
        allocate(G(m+n+1,n+1),source = zero)
        allocate(s(n+1),U(n,n),V_t(n+1,n+1))


        tau = norm2(abs(c))*this%tol

        if (all(abs(c(:m+1)) <= tau)) then ! Paper says that it shuold also be c(1) = ... = c(m+1), but not sure if that matters if they are all less than tau.
            this%q(1) = 1.0_sp
            this%mu = 0
            this%nu = 0
            return
        end if

        if(n > 0) then
            ! Fill G
            do i = 1,n+1
                G(i:n+m+1,i) = c(1:m+n+2-i)
            end do

            ! G_upper = G(:m+1,:) ! Maybe this matrix is unnecessary, just use a slice of G when needed
            G_tilde = G(m+2:,:)
        end if

        do while(this%nu > 0)
            call svd(G_tilde(:this%nu,:this%nu+1),s(:this%nu+1),u=U(:this%nu,:this%nu),vt=V_t(:this%nu+1,:this%nu+1))
            rho = 0
            do i = 1,this%nu
                if (s(i) > tau) rho = rho + 1
            end do

            if (this%nu > rho) then
                this%nu = rho
                this%mu = this%mu - (this%nu-rho)
            else
                exit
            end if

        end do

        if (this%nu == 0) then
            this%p = c
            this%q = 1
        else
            this%q(:this%nu+1) = V_t(this%nu+1,:this%nu+1)
            G_upper = G(:this%mu+1,:this%nu+1)
            call gemv('N',this%mu+1,this%nu+1,one,G_upper,this%mu+1,this%q(:this%nu+1),&
                    1,zero,this%p(:this%mu+1),1)
        end if

        ! count b < tol
        rho = 0 ! use rho for lambda
        do i = 1,this%nu+1
            if (abs(this%q(i)) <= this%tol) then
                rho = rho + 1
            else
                exit
            end if
        end do

        if (rho > 0) then
            this%p(:this%mu+1-rho) = this%p(rho+1:this%mu+1)
            this%q(:this%nu+1-rho) = this%q(rho+1:this%nu+1)

            this%p(rho+1:) = zero
            this%q(rho+1:) = zero

            this%mu = this%mu-rho
            this%nu = this%nu-rho
        end if

        ! count b < tol
        rho = 0 ! use rho for lambda
        do i = this%nu+1,1,-1
            if (abs(this%q(i)) <= this%tol) then
                rho = rho + 1
            else
                j = i+1
                exit
            end if
        end do

        if (rho > 0) then
            this%q(j:) = zero
            this%nu = this%nu-rho
        end if

        ! count a < tau
        rho = 0 ! use rho for lambda
        do i = this%mu+1,1,-1
            if (abs(this%p(i)) <= tau) then
                rho = rho + 1
            else
                j = i + 1
                exit
            end if
        end do

        if (rho > 0) then
            this%p(j:) = zero
            this%mu = this%mu-rho
        end if

        this%p = this%p(:this%mu+1)/this%q(1)
        this%q = this%q(:this%nu+1)/this%q(1)

    end subroutine pade_init_re_sp

    subroutine pade_init_re_dp(this,m,n,c,tol)
        class(pade_re(dp)), intent(out) :: this
        integer, intent(in) :: m
        integer, intent(in) :: n
        real(dp), intent(in) :: c(:)
        real(dp), intent(in), optional :: tol

        real(dp) :: tau
        real(dp), allocatable :: G_upper(:,:)
        real(dp), allocatable :: G_tilde(:,:)
        real(dp), allocatable :: G(:,:)
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: U(:,:)
        real(dp), allocatable :: V_t(:,:)
        integer :: rho
        integer :: i,j
        real(dp) :: one,zero

        if (m + n + 1 /= size(c)) then
            write(stderr,*) 'Error in pade_init_re_dp, len(c) not equal to m + n + 1.'
            error stop
        end if

        this%tol = optval(tol, DEFAULT_TOL_dp)
        this%c = c
        this%m = m
        this%n = n
        this%mu = m
        this%nu = n

        one = 1.0_dp
        zero = 0.0_dp
        allocate(this%p(m+1),this%q(n+1),source = zero)
        allocate(G_tilde(n,n+1),source = zero)
        allocate(G(m+n+1,n+1),source = zero)
        allocate(s(n+1),U(n,n),V_t(n+1,n+1))


        tau = norm2(abs(c))*this%tol

        if (all(abs(c(:m+1)) <= tau)) then ! Paper says that it shuold also be c(1) = ... = c(m+1), but not sure if that matters if they are all less than tau.
            this%q(1) = 1.0_dp
            this%mu = 0
            this%nu = 0
            return
        end if

        if(n > 0) then
            ! Fill G
            do i = 1,n+1
                G(i:n+m+1,i) = c(1:m+n+2-i)
            end do

            ! G_upper = G(:m+1,:) ! Maybe this matrix is unnecessary, just use a slice of G when needed
            G_tilde = G(m+2:,:)
        end if

        do while(this%nu > 0)
            call svd(G_tilde(:this%nu,:this%nu+1),s(:this%nu+1),u=U(:this%nu,:this%nu),vt=V_t(:this%nu+1,:this%nu+1))
            rho = 0
            do i = 1,this%nu
                if (s(i) > tau) rho = rho + 1
            end do

            if (this%nu > rho) then
                this%nu = rho
                this%mu = this%mu - (this%nu-rho)
            else
                exit
            end if

        end do

        if (this%nu == 0) then
            this%p = c
            this%q = 1
        else
            this%q(:this%nu+1) = V_t(this%nu+1,:this%nu+1)
            G_upper = G(:this%mu+1,:this%nu+1)
            call gemv('N',this%mu+1,this%nu+1,one,G_upper,this%mu+1,this%q(:this%nu+1),&
                    1,zero,this%p(:this%mu+1),1)
        end if

        ! count b < tol
        rho = 0 ! use rho for lambda
        do i = 1,this%nu+1
            if (abs(this%q(i)) <= this%tol) then
                rho = rho + 1
            else
                exit
            end if
        end do

        if (rho > 0) then
            this%p(:this%mu+1-rho) = this%p(rho+1:this%mu+1)
            this%q(:this%nu+1-rho) = this%q(rho+1:this%nu+1)

            this%p(rho+1:) = zero
            this%q(rho+1:) = zero

            this%mu = this%mu-rho
            this%nu = this%nu-rho
        end if

        ! count b < tol
        rho = 0 ! use rho for lambda
        do i = this%nu+1,1,-1
            if (abs(this%q(i)) <= this%tol) then
                rho = rho + 1
            else
                j = i+1
                exit
            end if
        end do

        if (rho > 0) then
            this%q(j:) = zero
            this%nu = this%nu-rho
        end if

        ! count a < tau
        rho = 0 ! use rho for lambda
        do i = this%mu+1,1,-1
            if (abs(this%p(i)) <= tau) then
                rho = rho + 1
            else
                j = i + 1
                exit
            end if
        end do

        if (rho > 0) then
            this%p(j:) = zero
            this%mu = this%mu-rho
        end if

        this%p = this%p(:this%mu+1)/this%q(1)
        this%q = this%q(:this%nu+1)/this%q(1)

    end subroutine pade_init_re_dp

    subroutine pade_init_cp_sp(this,m,n,c,tol)
        class(pade_cp(sp)), intent(out) :: this
        integer, intent(in) :: m
        integer, intent(in) :: n
        complex(sp), intent(in) :: c(:)
        real(sp), intent(in), optional :: tol

        real(sp) :: tau
        complex(sp), allocatable :: G_upper(:,:)
        complex(sp), allocatable :: G_tilde(:,:)
        complex(sp), allocatable :: G(:,:)
        real(sp), allocatable :: s(:)
        complex(sp), allocatable :: U(:,:)
        complex(sp), allocatable :: V_t(:,:)
        integer :: rho
        integer :: i,j
        complex(sp) :: one,zero

        if (m + n + 1 /= size(c)) then
            write(stderr,*) 'Error in pade_init_cp_sp, len(c) not equal to m + n + 1.'
            error stop
        end if

        this%tol = optval(tol, DEFAULT_TOL_sp)
        this%c = c
        this%m = m
        this%n = n
        this%mu = m
        this%nu = n

        one = 1.0_sp
        zero = 0.0_sp
        allocate(this%p(m+1),this%q(n+1),source = zero)
        allocate(G_tilde(n,n+1),source = zero)
        allocate(G(m+n+1,n+1),source = zero)
        allocate(s(n+1),U(n,n),V_t(n+1,n+1))


        tau = norm2(abs(c))*this%tol

        if (all(abs(c(:m+1)) <= tau)) then ! Paper says that it shuold also be c(1) = ... = c(m+1), but not sure if that matters if they are all less than tau.
            this%q(1) = 1.0_sp
            this%mu = 0
            this%nu = 0
            return
        end if

        if(n > 0) then
            ! Fill G
            do i = 1,n+1
                G(i:n+m+1,i) = c(1:m+n+2-i)
            end do

            ! G_upper = G(:m+1,:) ! Maybe this matrix is unnecessary, just use a slice of G when needed
            G_tilde = G(m+2:,:)
        end if

        do while(this%nu > 0)
            call svd(G_tilde(:this%nu,:this%nu+1),s(:this%nu+1),u=U(:this%nu,:this%nu),vt=V_t(:this%nu+1,:this%nu+1))
            rho = 0
            do i = 1,this%nu
                if (s(i) > tau) rho = rho + 1
            end do

            if (this%nu > rho) then
                this%nu = rho
                this%mu = this%mu - (this%nu-rho)
            else
                exit
            end if

        end do

        if (this%nu == 0) then
            this%p = c
            this%q = 1
        else
            this%q(:this%nu+1) = V_t(this%nu+1,:this%nu+1)
            G_upper = G(:this%mu+1,:this%nu+1)
            call gemv('N',this%mu+1,this%nu+1,one,G_upper,this%mu+1,this%q(:this%nu+1),&
                    1,zero,this%p(:this%mu+1),1)
        end if

        ! count b < tol
        rho = 0 ! use rho for lambda
        do i = 1,this%nu+1
            if (abs(this%q(i)) <= this%tol) then
                rho = rho + 1
            else
                exit
            end if
        end do

        if (rho > 0) then
            this%p(:this%mu+1-rho) = this%p(rho+1:this%mu+1)
            this%q(:this%nu+1-rho) = this%q(rho+1:this%nu+1)

            this%p(rho+1:) = zero
            this%q(rho+1:) = zero

            this%mu = this%mu-rho
            this%nu = this%nu-rho
        end if

        ! count b < tol
        rho = 0 ! use rho for lambda
        do i = this%nu+1,1,-1
            if (abs(this%q(i)) <= this%tol) then
                rho = rho + 1
            else
                j = i+1
                exit
            end if
        end do

        if (rho > 0) then
            this%q(j:) = zero
            this%nu = this%nu-rho
        end if

        ! count a < tau
        rho = 0 ! use rho for lambda
        do i = this%mu+1,1,-1
            if (abs(this%p(i)) <= tau) then
                rho = rho + 1
            else
                j = i + 1
                exit
            end if
        end do

        if (rho > 0) then
            this%p(j:) = zero
            this%mu = this%mu-rho
        end if

        this%p = this%p(:this%mu+1)/this%q(1)
        this%q = this%q(:this%nu+1)/this%q(1)

    end subroutine pade_init_cp_sp

    subroutine pade_init_cp_dp(this,m,n,c,tol)
        class(pade_cp(dp)), intent(out) :: this
        integer, intent(in) :: m
        integer, intent(in) :: n
        complex(dp), intent(in) :: c(:)
        real(dp), intent(in), optional :: tol

        real(dp) :: tau
        complex(dp), allocatable :: G_upper(:,:)
        complex(dp), allocatable :: G_tilde(:,:)
        complex(dp), allocatable :: G(:,:)
        real(dp), allocatable :: s(:)
        complex(dp), allocatable :: U(:,:)
        complex(dp), allocatable :: V_t(:,:)
        integer :: rho
        integer :: i,j
        complex(dp) :: one,zero

        if (m + n + 1 /= size(c)) then
            write(stderr,*) 'Error in pade_init_cp_dp, len(c) not equal to m + n + 1.'
            error stop
        end if

        this%tol = optval(tol, DEFAULT_TOL_dp)
        this%c = c
        this%m = m
        this%n = n
        this%mu = m
        this%nu = n

        one = 1.0_dp
        zero = 0.0_dp
        allocate(this%p(m+1),this%q(n+1),source = zero)
        allocate(G_tilde(n,n+1),source = zero)
        allocate(G(m+n+1,n+1),source = zero)
        allocate(s(n+1),U(n,n),V_t(n+1,n+1))


        tau = norm2(abs(c))*this%tol

        if (all(abs(c(:m+1)) <= tau)) then ! Paper says that it shuold also be c(1) = ... = c(m+1), but not sure if that matters if they are all less than tau.
            this%q(1) = 1.0_dp
            this%mu = 0
            this%nu = 0
            return
        end if

        if(n > 0) then
            ! Fill G
            do i = 1,n+1
                G(i:n+m+1,i) = c(1:m+n+2-i)
            end do

            ! G_upper = G(:m+1,:) ! Maybe this matrix is unnecessary, just use a slice of G when needed
            G_tilde = G(m+2:,:)
        end if

        do while(this%nu > 0)
            call svd(G_tilde(:this%nu,:this%nu+1),s(:this%nu+1),u=U(:this%nu,:this%nu),vt=V_t(:this%nu+1,:this%nu+1))
            rho = 0
            do i = 1,this%nu
                if (s(i) > tau) rho = rho + 1
            end do

            if (this%nu > rho) then
                this%nu = rho
                this%mu = this%mu - (this%nu-rho)
            else
                exit
            end if

        end do

        if (this%nu == 0) then
            this%p = c
            this%q = 1
        else
            this%q(:this%nu+1) = V_t(this%nu+1,:this%nu+1)
            G_upper = G(:this%mu+1,:this%nu+1)
            call gemv('N',this%mu+1,this%nu+1,one,G_upper,this%mu+1,this%q(:this%nu+1),&
                    1,zero,this%p(:this%mu+1),1)
        end if

        ! count b < tol
        rho = 0 ! use rho for lambda
        do i = 1,this%nu+1
            if (abs(this%q(i)) <= this%tol) then
                rho = rho + 1
            else
                exit
            end if
        end do

        if (rho > 0) then
            this%p(:this%mu+1-rho) = this%p(rho+1:this%mu+1)
            this%q(:this%nu+1-rho) = this%q(rho+1:this%nu+1)

            this%p(rho+1:) = zero
            this%q(rho+1:) = zero

            this%mu = this%mu-rho
            this%nu = this%nu-rho
        end if

        ! count b < tol
        rho = 0 ! use rho for lambda
        do i = this%nu+1,1,-1
            if (abs(this%q(i)) <= this%tol) then
                rho = rho + 1
            else
                j = i+1
                exit
            end if
        end do

        if (rho > 0) then
            this%q(j:) = zero
            this%nu = this%nu-rho
        end if

        ! count a < tau
        rho = 0 ! use rho for lambda
        do i = this%mu+1,1,-1
            if (abs(this%p(i)) <= tau) then
                rho = rho + 1
            else
                j = i + 1
                exit
            end if
        end do

        if (rho > 0) then
            this%p(j:) = zero
            this%mu = this%mu-rho
        end if

        this%p = this%p(:this%mu+1)/this%q(1)
        this%q = this%q(:this%nu+1)/this%q(1)

    end subroutine pade_init_cp_dp


    pure function pade_eval_scalar_re_sp(this,x) result(res)
        class(pade_re(sp)), intent(in) :: this
        real(sp), intent(in) :: x
        real(sp) :: res

        res = rat_eval(this%p,this%q,x)
    end function pade_eval_scalar_re_sp

    pure function pade_eval_scalar_re_dp(this,x) result(res)
        class(pade_re(dp)), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: res

        res = rat_eval(this%p,this%q,x)
    end function pade_eval_scalar_re_dp

    pure function pade_eval_scalar_cp_sp(this,x) result(res)
        class(pade_cp(sp)), intent(in) :: this
        complex(sp), intent(in) :: x
        complex(sp) :: res

        res = rat_eval(this%p,this%q,x)
    end function pade_eval_scalar_cp_sp

    pure function pade_eval_scalar_cp_dp(this,x) result(res)
        class(pade_cp(dp)), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: res

        res = rat_eval(this%p,this%q,x)
    end function pade_eval_scalar_cp_dp


    pure function pade_eval_vec_re_sp(this,x) result(res)
        class(pade_re(sp)), intent(in) :: this
        real(sp), intent(in) :: x(:)
        real(sp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        do i = 1, size(x)
            res(i) = rat_eval(this%p,this%q,x(i))
        end do
    end function pade_eval_vec_re_sp

    pure function pade_eval_vec_re_dp(this,x) result(res)
        class(pade_re(dp)), intent(in) :: this
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        do i = 1, size(x)
            res(i) = rat_eval(this%p,this%q,x(i))
        end do
    end function pade_eval_vec_re_dp

    pure function pade_eval_vec_cp_sp(this,x) result(res)
        class(pade_cp(sp)), intent(in) :: this
        complex(sp), intent(in) :: x(:)
        complex(sp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        do i = 1, size(x)
            res(i) = rat_eval(this%p,this%q,x(i))
        end do
    end function pade_eval_vec_cp_sp

    pure function pade_eval_vec_cp_dp(this,x) result(res)
        class(pade_cp(dp)), intent(in) :: this
        complex(dp), intent(in) :: x(:)
        complex(dp), allocatable :: res(:)

        integer :: i

        allocate(res(size(x)))
        do i = 1, size(x)
            res(i) = rat_eval(this%p,this%q,x(i))
        end do
    end function pade_eval_vec_cp_dp

end module robust_pade