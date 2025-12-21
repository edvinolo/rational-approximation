program pade_exp
    use kinds
    use polynomial_eval, only: poly_eval
    use robust_pade, only: pade_re,pade_cp
    implicit none

    integer, parameter :: N_coeff = 6
    integer, parameter :: m = 3
    integer, parameter :: n = 2
    integer, parameter :: wp = dp

    real(wp), parameter :: c_re(N_coeff) = 1.0_wp/[1.0_wp,1.0_wp,2.0_wp,6.0_wp,24.0_wp,120.0_wp]
    type(pade_re(wp)) :: real_pade

    complex(wp), parameter :: c_cp(N_coeff) = cmplx(c_re,y=0,kind=wp)
    type(pade_cp(wp)) :: cmplx_pade

    call real_pade%init(m,n,c_re)
    print *, real_pade%mu,real_pade%nu
    print *, real_pade%p
    print *, real_pade%q
    print *, 'e_pade =  ', real_pade%eval(1.0_wp), ' e = ', exp(1.0_wp)
    print *, 'e_pade =  ', real_pade%eval([1.0_wp,log(2.0_wp)]), ' e = ', exp([1.0_wp,log(2.0_wp)])
    print *, 'e_taylor = ', poly_eval(c_re,1.0_wp), poly_eval(c_re,log(2.0_wp))

    call cmplx_pade%init(m,n,c_cp)
    print *, cmplx_pade%mu,cmplx_pade%nu
    print *, cmplx_pade%p
    print *, cmplx_pade%q
    print *, 'e_pade =  ', cmplx_pade%eval((1.0_wp,0.0_wp)), ' e = ', exp((1.0_wp,0.0_wp))

end program pade_exp