program thiele_sin
    use kinds
    use rational_interpolation, only: thiele_interp,re,cp
    implicit none

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
    type(thiele_interp(wp,re)) :: real_thiele

    integer :: i,unit

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = sin(x_int(i))
    end do

    call real_thiele%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = real_thiele%eval(x(i))
    end do


    open(file='example/thiele_sin/int_data.dat',newunit=unit,action='write')
    do i = 1,N_int
        write(unit,'(es25.17e3,a,es25.17e3)') x_int(i), ' ', y_int(i)
    end do
    close(unit)

    open(file='example/thiele_sin/data.dat',newunit=unit,action='write')
    do i = 1,N
        write(unit,'(es25.17e3,a,es25.17e3)') x(i), ' ', y(i)
    end do
    close(unit)

end program thiele_sin