program main
    use wfl
    ! implicit none
    integer n, m
    real w(n, n), r_lon, r_lat
    n=60
    m=60
    r_lon=300.*1000.
    r_lat=300.*1000.
    call wfl_fan_filter(n, 30, r_lon, r_lat, w)
    print *, w
end program main