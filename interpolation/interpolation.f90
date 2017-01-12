subroutine py_bicubic(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y)
    use wfl, only: bicubic
    integer, intent(in):: m1, n1, m2, n2
    real, intent(in):: x1a(m1), x2a(n1), x1(m2), x2(n2), ya(m1, n1)
    real, intent(out):: y(m2, n2)
    call bicubic(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y)
end subroutine py_bicubic

subroutine py_bicubic_circle(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y)
    use wfl, only: bicubic
    integer, intent(in):: m1, n1, m2, n2
    real, intent(in):: x1a(m1), x2a(n1), x1(m2), x2(n2), ya(m1, n1)
    real, intent(out):: y(m2, n2)
    call bicubic_circle(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y)
end subroutine py_bicubic_circle

subroutine py_bicubic_spline(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y)
    use wfl, only: bicubic_spline
    integer, intent(in):: m1, n1, m2, n2
    real, intent(in):: x1a(m1), x2a(n1), ya(m1, n1), x1(m2), x2(n2)
    real, intent(out):: y(m2, n2)
    call bicubic_spline(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y)
end subroutine py_bicubic_spline

subroutine py_bicubic_spline_circle(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y)
    use wfl, only: bicubic_spline_circle
    integer, intent(in):: m1, n1, m2, n2
    real, intent(in):: x1a(m1), x2a(n1), ya(m1, n1), x1(m2), x2(n2)
    real, intent(out):: y(m2, n2)
    call bicubic_spline_circle(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y)
end subroutine py_bicubic_spline_circle

subroutine py_bilinear(x1a, x2a, ya, x1, x2, y)
    use wfl, only: bilinear
    real*8, intent(in):: x1a(2), x2a(2), ya(4), x1, x2
    real*8, intent(out):: y
    call bilinear(x1a, x2a, ya, x1, x2, y)
end subroutine py_bilinear

! subroutine py_bilinear_circle(n1, m1, x1a, x2a, ya, n2, m2, x1, x2, y, blankvalue)
!     use wfl, only: bilinear_circle
!     integer, intent(in):: n1, m1, n2, m2
!     real, intent(in):: x1a(n1), x2a(m1), ya(n1, m1), x1(n2), x2(m2)
!     real, intent(out):: y(n2, m2)
!     real, optional, intent(in):: blankvalue
! end subroutine py_bilinear_circle

subroutine py_bilinear_square(n1, m1, x1a, x2a, ya, n2, m2, x1, x2, y)
    use wfl, only: bilinear_square
    integer, intent(in):: n1, m1, n2, m2
    real, intent(in):: x1a(n1), x2a(m1), ya(n1, m1), x1(n2), x2(m2)
    real, intent(out):: y(n2, m2)
    call bilinear_square(n1, m1, x1a, x2a, ya, n2, m2, x1, x2, y)
end subroutine py_bilinear_square

subroutine py_linear_int(n, t, x, m, tint, xout, mask)
    use wfl, only: linear_int
    implicit none
    integer, intent(in):: n, m
    real, intent(in):: x(n), t(n), tint(m)
    logical, optional, intent(in):: mask(n)
    real, intent(out):: xout(m)
    call linear_int(n, t, x, m, tint, xout, mask)
end subroutine py_linear_int

subroutine py_linear_int_miss(n, x, xout, mask, method)
    use wfl, only: linear_int_miss
    integer, intent(in):: n
    real, intent(in):: x(n)
    logical, intent(in):: mask(n)
    real, intent(out):: xout(n)
    character(len=*), optional, intent(in):: method
    call linear_int_miss(n, x, xout, mask, method)
end subroutine py_linear_int_miss

subroutine py_points9_ave_circle(nlon, nlat, x, xout, weight, blankvalue, flag)
    use wfl, only: points9_ave_circle
    implicit none
    integer, intent(in):: nlon, nlat
    real *8, dimension(nlon, nlat), intent(in):: x
    real *8, dimension(nlon, nlat), intent(out):: xout
    real *8, optional, intent(in):: blankvalue, weight(nlon, nlat)
    logical, optional, intent(in) :: flag
    if(present(flag) .and. flag) then
        call points9_ave_circle(nlon, nlat, x, xout, weight, blankvalue)
    else
        call points9_ave_circle(nlon, nlat, x, xout, weight)
    endif
end subroutine py_points9_ave_circle

subroutine py_polint(xa, ya, n, x, y, dy)
    use wfl, only: polint
    integer, intent(in):: n
    real, intent(in)::x, xa(n), ya(n)
    real, intent(out):: y, dy
    call polint(xa, ya, n, x, y, dy)
end subroutine py_polint

subroutine py_polint2(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y, morder, norder)
    use wfl, only: polint2
    integer, intent(in):: m1, n1, m2, n2, morder, norder
    real, intent(in):: x1a(m1), x2a(n1), ya(m1, n1), x1(m2), x2(n2)
    real, intent(out):: y(m2, n2)
    call polint2(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y, morder, norder)
end subroutine py_polint2

subroutine py_polint2_circle(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y, morder, norder)
    use wfl, only: polint2_circle
    integer, intent(in):: m1, n1, m2, n2, morder, norder
    real, intent(in):: x1a(m1), x2a(n1), ya(m1, n1), x1(m2), x2(n2)
    real, intent(out):: y(m2, n2)
    call polint2_circle(m1, n1, x1a, x2a, ya, m2, n2, x1, x2, y, morder, norder)
end subroutine py_polint2_circle

subroutine py_spline_nr(x, y, n, yp1, ypn, y2)
    use wfl, only: spline_nr
    integer, intent(in):: n
    real, intent(in):: yp1, ypn, x(n), y(n)
    real, intent(out):: y2(n)
    call spline_nr(x, y, n, yp1, ypn, y2)
end subroutine py_spline_nr

subroutine py_spline3(n1, t1, x1, n2, t2, x2, d1, d2, c)
    use wfl, only: spline3
    integer, intent(in)::n1, n2
    real*8, intent(in):: t1(n1), x1(n1), t2(n2)
    real*8, dimension(n2), intent(out)::x2, d1, d2, c
    call spline3(n1, t1, x1, n2, t2, x2, d1, d2, c)
end subroutine py_spline3

subroutine py_spline3_miss(n1, t1, x1, n2, t2, x2, mask, method)
    use wfl, only: spline3_miss
    integer, intent(in) :: n1, n2
    real *8, intent(in) :: t1(n1), x1(n1), t2(n2)
    real *8, intent(out):: x2(n2)
    logical, intent(in):: mask(n1)
    character(len=*), optional, intent(in):: method
    call spline3_miss(n1, t1, x1, n2, t2, x2, mask, method)
end subroutine py_spline3_miss

subroutine py_splint(xa, ya, y2a, n, x, y)
    use wfl, only: splint
    integer, intent(in):: n
    real, intent(in) :: x, xa(n), y2a(n), ya(n)
    real, intent(out):: y
    call splint(xa, ya, y2a, n, x, y)
end subroutine py_splint