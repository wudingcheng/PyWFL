subroutine py_aic(m_order_max,m,rsdus,rs,m_order)
    use wfl, only:aic
    integer,intent(in):: m_order_max,m
    real*8,intent(in):: rsdus(m_order_max)
    real*8,intent(out)::rs(m_order_max)
    integer,intent(out):: m_order
    call aic(m_order_max,m,rsdus,rs,m_order)
end subroutine py_aic

! function py_arth(first,increment,n)
!     use wfl, only:arth
!     real *8, intent(in) :: first,increment
!     integer, intent(in) :: n
!     real *8 :: result_arth
!     result_arth = arth(first,increment,n)
! end function py_arth

subroutine  py_bic(m_order_max,m,rsdus,rs,m_order)
    use wfl, only: bic
    integer,intent(in):: m_order_max,m
    real*8,intent(in):: rsdus(m_order_max)
    real*8,intent(out)::rs(m_order_max)
    integer,intent(out):: m_order
    call bic(m_order_max,m,rsdus,rs,m_order)
end subroutine py_bic

function py_chinv(p,dof) result(chinv_result)
    use wfl, only: chinv
    real*8,intent(in):: p,dof
    real*8 :: chinv_result
    chinv_result = chinv(p,dof)
end function py_chinv

subroutine py_cor2(n,a,b,co,colev,dof)
    use wfl, only: cor2
    integer,intent(in)::n
    real,dimension(n),intent(in):: a,b
    real,intent(out)::co,colev
    integer,optional,intent(in):: dof
    call cor2(n, a, b, co, colev, dof)
end subroutine

subroutine py_corfre(n,dt,a,b,out)
    use wfl, only:corfre
    integer,intent(in)::n
    real,dimension(n),intent(in):: a,b
    real,intent(out):: out(5,n-1)
    call corfre(n, dt, a, b, out)
end subroutine py_corfre