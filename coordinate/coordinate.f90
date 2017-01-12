subroutine py_wfl_2lam(alon, alat, cm, ps, pn, ol, x, y, ea, ef)
    use wfl, only: wfl_2lam
    implicit none
    real*8, intent(in) :: alon, alat, cm, ps, pn, ol
    real*8, intent(out) :: x, y
    real*8, optional, intent(in) :: ea, ef
    call wfl_2lam(alon, alat, cm, ps, pn, ol, x, y, ea, ef)
end subroutine py_wfl_2lam

subroutine py_wfl_blh2xyz(b, l, h, x, y, z, ea, ef)
    use wfl, only: wfl_blh2xyz
    real *8, intent(in) :: b, l, h
    real *8, intent(out) :: x, y, z
    real *8, optional, intent(in) :: ea, ef
    call wfl_blh2xyz(b, l, h, x, y, z, ea, ef)
end subroutine py_wfl_blh2xyz

subroutine py_wfl_map_xyz2blh(x, y, z, b, l, h, ea, ef)
    use wfl, only: wfl_map_xyz2blh
    real *8, intent(in) :: x, y, z
    real *8, intent(out) :: b, l, h
    real *8, optional, intent(in) :: ea, ef
    call wfl_map_xyz2blh(x, y, z, b, l, h, ea, ef)
end subroutine py_wfl_map_xyz2blh

subroutine py_wfl_invlam(x, y, cm, ps, pn, ol, alon, alat, ea, ef)
    use wfl, only: wfl_invlam
    implicit none
    real*8, intent(in) :: x, y, cm, ps, pn, ol
    real*8, intent(out) :: alon, alat
    real*8, optional, intent(in) :: ea, ef
    call wfl_invlam(x, y, cm, ps, pn, ol, alon, alat, ea, ef)
end subroutine py_wfl_invlam

subroutine py_wfl_xyz2enu(nn, x, y, z, e, n, u, ea, ef, s_lat, s_lon)
    use wfl, only: wfl_xyz2enu
    integer, intent(in) :: nn
    real *8, dimension(nn), intent(in) :: x, y, z
    real *8, dimension(nn), intent(out) :: e, n, u
    real *8, optional, intent(in) :: ea, ef, s_lat, s_lon
    call wfl_xyz2enu(nn, x, y, z, e, n, u, ea, ef, s_lat, s_lon)
end subroutine py_wfl_xyz2enu

subroutine  py_wfl_map_gc2gd( x, y, z, phi, elong, height, ea, ef)
    use wfl, only: wfl_map_gc2gd
    real *8, intent(in) ::  x, y, z
    real *8, intent(out) ::  elong, phi, height
    real *8, optional, intent(in) :: ea, ef
    call wfl_map_gc2gd( x, y, z, phi, elong, height, ea, ef)
end subroutine py_wfl_map_gc2gd

subroutine  py_wfl_map_datum(index_reference_name, a, f)
    use wfl, only: wfl_map_datum
    integer, intent(in) ::  index_reference_name
    real*8, intent(out) :: a, f
    call wfl_map_datum(index_reference_name, a, f)
end subroutine py_wfl_map_datum

subroutine py_wfl_map_geo2ups(lat,lon,x,y, ea, ef)
    use wfl, only: wfl_map_geo2ups
    real *8,intent(in) ::  lat,lon
    real *8,intent(out) :: x, y
    real *8,optional,intent(in) :: ea,ef
    call wfl_map_geo2ups(lat,lon,x,y, ea, ef)
end subroutine py_wfl_map_geo2ups