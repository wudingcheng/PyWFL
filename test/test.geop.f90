program main
    use wfl
    ! implicit none
    parameter(n=10)
    real x(n)
    real first, factor
    first=1.0
    factor=2.0
    x=geop(first,factor,n)
    print *, x
end program main