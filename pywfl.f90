subroutine py_factorial(n, result)
    use wfl, only: factorial
    real *8, intent(in) :: n
    real *8, intent(out) :: result
    result = factorial(n)
end subroutine py_factorial


subroutine py_bilinear(x1a,x2a,ya,x1,x2,y)
    use wfl, only: bilinear
    real*8,intent(in):: x1a(2),x2a(2),ya(4),x1,x2
    real*8,intent(out):: y
    call bilinear(x1a,x2a,ya,x1,x2,y)
end subroutine py_bilinear

subroutine py_bilinear_circle(n1,m1,x1a,x2a,ya,n2,m2,x1,x2,y,blankvalue)
    use wfl, only: bilinear_circle
    integer,intent(in):: n1,m1,n2,m2
    real,intent(in):: x1a(n1),x2a(m1),ya(n1,m1),x1(n2),x2(m2)
    real,intent(out):: y(n2,m2)
    real,optional,intent(in):: blankvalue
end subroutine py_bilinear_circle

subroutine py_bilinear_square(n1,m1,x1a,x2a,ya,n2,m2,x1,x2,y)
    use wfl, only: bilinear_square
    integer,intent(in):: n1,m1,n2,m2
    real,intent(in):: x1a(n1),x2a(m1),ya(n1,m1),x1(n2),x2(m2)
    real,intent(out):: y(n2,m2)
    call bilinear_square(n1,m1,x1a,x2a,ya,n2,m2,x1,x2,y)
end subroutine py_bilinear_square


subroutine py_bicubic(m1,n1,x1a,x2a,ya,m2,n2,x1,x2,y)
    use wfl, only: bicubic
    integer,intent(in):: m1,n1,m2,n2
    real,intent(in):: x1a(m1),x2a(n1),x1(m2),x2(n2),ya(m1,n1)
    real,intent(out):: y(m2,n2)
    call bicubic(m1,n1,x1a,x2a,ya,m2,n2,x1,x2,y)
end subroutine py_bicubic

subroutine py_bicubic_circle(m1,n1,x1a,x2a,ya,m2,n2,x1,x2,y)
    use wfl, only: bicubic
    integer,intent(in):: m1,n1,m2,n2
    real,intent(in):: x1a(m1),x2a(n1),x1(m2),x2(n2),ya(m1,n1)
    real,intent(out):: y(m2,n2)
    call bicubic_circle(m1,n1,x1a,x2a,ya,m2,n2,x1,x2,y)
end subroutine py_bicubic_circle

subroutine py_bicubic_spline(m1,n1,x1a,x2a,ya,m2,n2,x1,x2,y)
    use wfl, only: bicubic_spline
    integer,intent(in):: m1,n1,m2,n2
    real,intent(in):: x1a(m1),x2a(n1),ya(m1,n1),x1(m2),x2(n2)
    real,intent(out):: y(m2,n2)
    call bicubic_spline(m1,n1,x1a,x2a,ya,m2,n2,x1,x2,y)
end subroutine py_bicubic_spline

subroutine py_bicubic_spline_circle(m1,n1,x1a,x2a,ya,m2,n2,x1,x2,y)
    use wfl, only: bicubic_spline_circle
    integer,intent(in):: m1,n1,m2,n2
    real,intent(in):: x1a(m1),x2a(n1),ya(m1,n1),x1(m2),x2(n2)
    real,intent(out):: y(m2,n2)
    call bicubic_spline_circle(m1,n1,x1a,x2a,ya,m2,n2,x1,x2,y)
end subroutine py_bicubic_spline_circle

