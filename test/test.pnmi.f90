program main
    use wfl, only: pnmi
    ! implicit none
    integer n
    parameter(n=3)
    real *8 :: ts, tn, pinmsn(0:n, 0:n)
    ts = 30.0
    tn = 40.0
    call pnmi(n, ts, tn, pinmsn)
    print *, pinmsn

end program main