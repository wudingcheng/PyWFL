subroutine py_wfl_svd_c(a,m,n,u,w,v)
    use wfl, only: wfl_svd
    integer,intent(in) :: m,n
    complex *8,intent(in) :: a(m,n)
    complex *8,intent(out) :: u(m,n),w(n),v(n,n)
    call wfl_svd(a,m,n,u,w,v)
end subroutine py_wfl_svd_c

subroutine py_wfl_svd_d(a,m,n,u,w,v)
    use wfl, only: wfl_svd
    integer,intent(in) :: m,n
    real *8,intent(in) :: a(m,n)
    real *8,intent(out) :: u(m,n),w(n),v(n,n)
    call wfl_svd(a,m,n,u,w,v)
end subroutine py_wfl_svd_d