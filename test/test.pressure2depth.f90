program main
    use wfl
    implicit none
    real *8 :: depth, pressure, lat

    depth = 7321.45
    pressure = 7500.004
    lat = 30
    pressure = wfl_depth2pressure(depth, lat)
    print *, "depth -> pressure", pressure
    print *, "pressure -> depth", wfl_pressure2depth(pressure, lat)

end program main