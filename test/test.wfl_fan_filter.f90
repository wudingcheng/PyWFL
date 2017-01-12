program main
    use wfl
    ! implicit none
    integer, parameter :: n=5, m=5
    real w(0:n, 0:n), r_lon, r_lat
    integer i
    ! n=60
    ! m=60
    r_lon=300.*1000.
    r_lat=300.*1000.
    call wfl_fan_filter(n, n, r_lon, r_lat, w)
    do i = 0, n
        print *, w(:, i)
    enddo
end program main