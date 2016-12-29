subroutine py_arco(n,m,x,v,r,ch,best_order)
    use wfl
    integer,intent(in):: n,m
    real *8,intent(in):: x(n)
    real *8,intent(out):: v(0:m),r(0:m)
    character,optional, intent(in):: ch*3
    integer,optional,intent(out):: best_order
    if(present(ch) .and. ch .eq. 'fpe') then
        call arco(n,m,x,v,r,ch)
    else
        call arco(n,m,x,v,r,best_order=best_order)
    endif
end subroutine py_arco

subroutine py_emd(n,t,x,m,imf,std)
    use wfl, only: emd
    integer,intent(in):: n,m
    real,intent(in):: x(n),t(n)
    real,intent(out) :: imf(n,m)
    real,optional,intent(in):: std
    if(present(std) .and. std .eq. 0) then
        call emd(n, t, x, m, imf)
    else
        call emd(n, t, x, m, imf, std)
    endif
end subroutine py_emd
