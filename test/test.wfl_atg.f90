program main
    use wfl
    real*8 :: s,t,p0
    ! check value: atg=3.255976d-4 c/dbar for s=40 (pss-78), t=40 deg celsius, p=10000 decibars
    s=40.d0
    t=40.d0
    p0=1.d8
    print *, wfl_atg(s,t,p0), p0

end program main