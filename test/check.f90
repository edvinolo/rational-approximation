program check
use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
use testdrive, only: run_testsuite, new_testsuite, testsuite_type
use test_interpolation, only: collect_interpolation_suite
implicit none
integer :: stat,is
type(testsuite_type), allocatable :: testsuites(:)
character(len=*), parameter :: fmt = '("#", *(1x, a))'

stat = 0
testsuites = [&
    new_testsuite("Interpolation test suite",collect_interpolation_suite)&
    ]

write(stdout,'(a)') 'Running tests!'

do is = 1, size(testsuites)
    write(stderr,fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, stdout, stat)
end do

if (stat > 0) then
    write(stderr,'(i0,1x,a)') stat, "test(s) failed!"
    error stop
else if (stat == 0) then
    write(stdout,'(a)') 'All tests passed succesfully!'
end if

end program check
