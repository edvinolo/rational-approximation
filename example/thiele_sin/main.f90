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

    complex(wp), parameter :: dx_int_z = (0.4_wp,0.1_wp)
    complex(wp), parameter :: dx_z = (0.13_wp,0.1_wp)
    complex(wp), parameter :: x_int_start_z = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_start_z = (-0.2_wp,-0.2_wp)

    complex(wp), allocatable :: x_int_z(:)
    complex(wp), allocatable :: x_z(:)
    complex(wp), allocatable :: y_int_z(:)
    complex(wp), allocatable :: y_z(:)
    type(thiele_interp(wp,cp)) :: cmplx_thiele

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


    allocate(x_int_z(N_int),y_int_z(N_int),x_z(N),y_z(N))

    do i = 1,N_int
        x_int_z(i) = x_int_start_z + (i-1)*dx_int_z
        y_int_z(i) = sin(x_int_z(i))
    end do

    call cmplx_thiele%init(x_int_z,y_int_z)

    do i = 1,N
        x_z(i) = x_start_z + (i-1)*dx_z
    end do

    do i = 1,N
        y_z(i) = cmplx_thiele%eval(x_z(i))
    end do

    y_int_z = cmplx_thiele%eval(x_int_z)

    ! print *, maxval(abs(sin(x_int_z)-y_int_z)/abs(y_int_z))

end program thiele_sin