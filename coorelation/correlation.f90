subroutine py_cor2(n, a, b, co, colev, dof)
    use wfl, only: cor2
    integer, intent(in)::n
    real, dimension(n), intent(in):: a, b
    real, intent(out)::co, colev
    integer, optional, intent(in):: dof
    if(present(dof) .and. dof .ne. 0) then
        call cor2(n, a, b, co, colev, dof)
    else
        call cor2(n, a, b, co, colev)
    endif
end subroutine py_cor2

subroutine py_corfre(n, dt, a, b, out)
    use wfl, only: corfre
    integer, intent(in)::n
    real, dimension(n), intent(in):: a, b
    real, intent(out):: out(5, n-1)
    call corfre(n, dt, a, b, out)
end subroutine py_corfre

subroutine py_corlag(n, dt, a, b, dof, out)
    use wfl, only: corlag
    integer, intent(in)::n
    real, intent(in)::dt, dof
    real, dimension(n), intent(in):: a, b
    real, intent(out):: out(4, -n/2+2:n/2-2)
    call corlag(n, dt, a, b, dof, out)
end subroutine py_corlag

subroutine py_wfl_arcov(n, x, x_arcov, ch)
    use wfl, only: wfl_arcov
    integer, intent(in) :: n
    real*8, intent(in) :: x(n)
    real*8, intent(out):: x_arcov(0:n-1)
    character(len=1), optional, intent(in) :: ch
    if(present(ch) .and. ch .eq. 'n') then
        call wfl_arcov(n, x, x_arcov, ch)
    else
        call wfl_arcov(n, x, x_arcov)
    endif
end subroutine py_wfl_arcov

function py_wfl_edof(n, x, y) result(n_edof)
    use wfl, only: wfl_edof
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: x, y
    integer n_edof
    n_edof = wfl_edof(n, x, y)
end function py_wfl_edof

subroutine  py_wfl_fre_response(n, dt, x, y, nw, fre_res, p, amp_ci, ph_ci)
    use wfl, only: wfl_fre_response
    integer, intent(in) :: n
    real *8, dimension(n), intent(in) :: x, y
    real *8, intent(in) :: nw, dt
    real *8, dimension(n/2, 4), intent(out) :: fre_res
    real *8, optional, dimension(n/2, 2), intent(out) :: amp_ci, ph_ci
    real *8, optional, intent(in) :: p
    call wfl_fre_response(n, dt, x, y, nw, fre_res, p, amp_ci, ph_ci)
end subroutine py_wfl_fre_response

subroutine  py_wfl_mtm_coh(n, dt, x, y, nw, coh_spec, cross_spec, p, f_test)
    use wfl, only: wfl_mtm_coh
    integer, intent(in) :: n
    real *8, dimension(n), intent(in) :: x, y
    real *8, intent(in) :: nw, dt
    real *8, dimension(n/2, 3), intent(out) :: coh_spec
    real *8, optional, dimension(n/2, 3), intent(out) :: cross_spec, f_test(n/2)
    real *8, optional, intent(in) :: p
    call wfl_mtm_coh(n, dt, x, y, nw, coh_spec, cross_spec, p, f_test)
end subroutine py_wfl_mtm_coh

subroutine  py_wfl_multicor(n, m, x, x_out, out_type, fft, remove_mean)
    use wfl, only: wfl_multicor, wfl_multicor_fft
    integer, intent(in) :: n, m
    real *8, dimension(n, m), intent(in) :: x
    real *8, dimension(0:n-1, m, m), intent(out) :: x_out
    character(len=*), intent(in) :: out_type
    character(len=*), optional, intent(in) :: remove_mean
    logical, optional, intent(in) :: fft
    print *, x
    if(present(fft) .and. fft) then
        print *, "FFT"
        call wfl_multicor_fft(n, m, x, x_out, out_type, remove_mean)
    else
        print *, "No FFT"
        call wfl_multicor(n, m, x, x_out, out_type, remove_mean)
    endif
    do i=0,n-1
        write(*,'(1x,i3,3x,f20.6,2f15.6)') i, x_out(i,1,2)
    enddo
end subroutine py_wfl_multicor


