subroutine py_fft_filter(n, x, dt, per1, per2, method, xout, window)
    use wfl
    integer , intent(in)          :: n
    real *8, intent(in)          :: x(n)
    real *8, intent(in)          :: dt, per1, per2
    character(len=*), intent(in) :: method
    real *8, intent(out)         :: xout(n)
    real *8, intent(in)          :: window
    call fft_filter(n, x, dt, per1, per2, method, xout, window)
end subroutine py_fft_filter

subroutine py_vondrak(n,e,t,x,u,em)
    use wfl
    integer,intent(in)::n
    real*8,intent(in)::e,t(n),x(n)
    real*8,intent(out)::u(n),em
    call vondrak(n,e,t,x,u,em)
end subroutine py_vondrak

function py_gauss_jekeli(r,distance_from_center, ch)  result(gauss_result)
    use wfl
    implicit none
    real *8, intent(in)          :: r, distance_from_center
    character(len=*), intent(in) :: ch
    real *8 :: gauss_result
    if(ch .eq. 'distance') then
        write(*, *) ch
        gauss_result=gauss_jekeli(r, distance_from_center)
    else
        write(*, *) 'no', ch
        gauss_result=gauss_jekeli(r, distance_from_center, ch)
    endif
end function py_gauss_jekeli

subroutine py_gauss_jekeli_fre(n,r,w)
    use wfl
    implicit none
    integer, intent(in) :: n
    real *8, intent(in) :: r
    real *8, intent(out):: w(0:n)
end subroutine py_gauss_jekeli_fre

subroutine py_gauss_jekeli_fre_nm(n,m_max,r_lon, r_lat,w)
    use wfl
    integer, intent(in) :: n, m_max
    real *8, intent(in):: r_lon, r_lat
    real *8, intent(out):: w(0:n,0:n)
    call gauss_jekeli_fre_nm(n,m_max,r_lon, r_lat,w)
end subroutine py_gauss_jekeli_fre_nm

subroutine py_weight_aver(n,x,npoints,xout,weight_no)
    use wfl
    integer, intent(in) :: n,npoints
    real,    intent(in) ::  x(n)
    real,   intent(out) :: xout(n)
    logical, intent(in) :: weight_no
    if(weight_no) then
        call weight_aver(n, x, npoints, xout, 1.0)
    else
        call weight_aver(n, x, npoints, xout)
    endif
end subroutine py_weight_aver

subroutine py_wfl_butterworth(l, fl, fh, iband0, n, dt, x0, xout, method, edge_effect, fre, amp, ph)
    use wfl
    integer,intent(in) :: l,n,iband0
    real *8,intent(in) :: fl,fh,dt
    real *8,intent(in),dimension(n) :: x0
    real *8,intent(out),dimension(n) :: xout
    character(len=*), intent(in):: method, edge_effect
    real *8,intent(out), dimension(n/2-1) :: fre, amp, ph
    if( method .eq. 'no_phase_shift') then
        if(edge_effect .eq. 'no') then
            call wfl_butterworth(l,fl,fh,iband0,n,dt,x0,xout,method=method,edge_effect=edge_effect,fre=fre,amp=amp,ph=ph)
        else
            call wfl_butterworth(l, fl, fh, iband0, n, dt, x0, xout, method=method, fre=fre, amp=amp, ph=ph)
        endif
    else
        call wfl_butterworth(l, fl, fh, iband0, n, dt, x0, xout, fre=fre, amp=amp, ph=ph)
    endif
end subroutine py_wfl_butterworth

function py_wfl_butterworth_order(fl,fh,dt,iband) result(l)
    use wfl
    real *8, intent(in) :: fl,fh,dt
    integer, intent(in) :: iband
    integer :: l
    l = wfl_butterworth_order(fl,fh,dt,iband)
end function py_wfl_butterworth_order

subroutine  py_wfl_running_mean(n, x, m, ans, weight)
    use wfl
    integer, intent(in) :: n, m
    real *8, intent(in) :: x(n)
    real *8, intent(in) :: weight(m)
    real *8, dimension(n), intent(out) :: ans
    call wfl_running_mean(n, x, m, ans, weight)
end subroutine py_wfl_running_mean