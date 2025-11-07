program MTT_cos
    use kinds
    use rational_interpolation, only: MTT_interp_re,MTT_interp_cp
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
    type(MTT_interp_re(wp)) :: real_MTT

    complex(wp), parameter :: dx_int_z = (0.4_wp,0.1_wp)
    complex(wp), parameter :: dx_z = (0.13_wp,0.1_wp)
    complex(wp), parameter :: x_int_start_z = (0.0_wp,0.0_wp)
    complex(wp), parameter :: x_start_z = (-0.2_wp,-0.2_wp)

    complex(wp), allocatable :: x_int_z(:)
    complex(wp), allocatable :: x_z(:)
    complex(wp), allocatable :: y_int_z(:)
    complex(wp), allocatable :: y_z(:)
    type(MTT_interp_cp(wp)) :: cmplx_MTT

    integer :: i,unit,t
    real(wp) :: abs_re,abs_cp

    allocate(x_int(N_int),y_int(N_int),x(N),y(N))

    do i = 1,N_int
        x_int(i) = x_int_start + (i-1)*dx_int
        y_int(i) = cos(x_int(i))
    end do

    call real_MTT%init(x_int,y_int)

    do i = 1,N
        x(i) = x_start + (i-1)*dx
    end do

    do i = 1,N
        y(i) = real_MTT%eval(x(i))
    end do

    open(file='example/MTT_cos/int_data.dat',newunit=unit,action='write')
    do i = 1,N_int
        write(unit,'(es25.17e3,a,es25.17e3)') x_int(i), ' ', y_int(i)
    end do
    close(unit)

    open(file='example/MTT_cos/data.dat',newunit=unit,action='write')
    do i = 1,N
        write(unit,'(es25.17e3,a,es25.17e3)') x(i), ' ', y(i)
    end do
    close(unit)

    y_int = real_MTT%eval(x_int)
    abs_re = maxval(abs(cos(x_int)-y_int)/abs(y_int))

    ! Check that stage (b) of MTT works
    x_int = [-1.0_wp,0.0_wp,1.0_wp,2.0_wp,3.0_wp]
    y_int = (x_int**2-1.0_wp)/(x_int+2.0_wp)
    call real_MTT%init(x_int,y_int)
    t = real_MTT%t ! Should be 4

    allocate(x_int_z(N_int),y_int_z(N_int),x_z(N),y_z(N))

    do i = 1,N_int
        x_int_z(i) = x_int_start_z + (i-1)*dx_int_z
        y_int_z(i) = cos(x_int_z(i))
    end do

    call cmplx_MTT%init(x_int_z,y_int_z)

    do i = 1,N
        x_z(i) = x_start_z + (i-1)*dx_z
    end do

    do i = 1,N
        y_z(i) = cmplx_MTT%eval(x_z(i))
    end do

    y_int_z = cmplx_MTT%eval(x_int_z)
    abs_cp = maxval(abs(cos(x_int_z)-y_int_z)/abs(y_int_z))

    open(file='example/MTT_cos/checks.dat',newunit=unit,action='write')
    write(unit,'(a,es25.17e3,a)')  'Re rel. abs. interp. error ', abs_re
    write(unit,'(a,es25.17e3,a)')  'Cp rel. abs. interp. error ', abs_cp
    write(unit,'(a,i0,a)')  'Check t for (z^2-1)/(z+2): ', t, '. Should be 4.'
    close(unit)

end program MTT_cos