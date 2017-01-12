module test_module
    implicit none

    interface test
        module procedure test_d
        module procedure test_c
    end interface test

    contains
        subroutine test_d(args)
        implicit none
        real *8 :: args
        print *, args
    end subroutine test_d

    subroutine test_c(args)
        implicit none
        complex :: args
        print *, args
    end subroutine test_c
end module test_module

