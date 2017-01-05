subroutine py_arco(n, m, x, v, r, ch, best_order)
    use wfl
    integer, intent(in):: n, m
    real *8, intent(in):: x(n)
    real *8, intent(out):: v(0:m), r(0:m)
    character, optional, intent(in):: ch*3
    integer, optional, intent(out):: best_order
    if(present(ch) .and. ch .eq. 'fpe') then
        call arco(n, m, x, v, r, ch)
    else
        call arco(n, m, x, v, r, best_order=best_order)
    endif
end subroutine py_arco

subroutine py_emd(n, t, x, m, imf, std)
    use wfl, only: emd
    integer, intent(in):: n, m
    real, intent(in):: x(n), t(n)
    real, intent(out) :: imf(n, m)
    real, optional, intent(in):: std
    if(present(std) .and. std .eq. 0) then
        call emd(n, t, x, m, imf)
    else
        call emd(n, t, x, m, imf, std)
    endif
end subroutine py_emd

subroutine py_emd_tol(n, t, x, m, imf)
    use wfl, only: emd_tol
    integer, intent(in):: n, m
    real, intent(in):: x(n), t(n)
    real, intent(out) :: imf(n, m)
    call emd_tol(n, t, x, m, imf)
end subroutine py_emd_tol

subroutine py_fft_spec(n, dt, x, p, spec, lag1_co)
    use wfl, only: fft_spec
    integer, intent(in) :: n
    real *8, intent(in) :: dt, x(n), p
    real *8, intent(out) :: spec(n/2, 4)
    real *8, optional :: lag1_co
    if(present(lag1_co) .and. lag1_co .eq. 0) then
        call fft_spec(n, dt, x, p, spec)
    else
        call fft_spec(n, dt, x, p, spec, lag1_co)
    endif
end subroutine py_fft_spec

subroutine py_fft_fre(n, dt, x, pe, am, ph, ch)
    use wfl, only: fft_fre
    integer, intent(in):: n
    real, intent(in):: dt, x(n)
    real, intent(out):: pe(n-1), am(n-1), ph(n-1)
    character(len=*), optional, intent(in):: ch
    call fft_fre(n, dt, x, pe, am, ph, ch)
end subroutine py_fft_fre

subroutine py_fft_period(n, dt, x, pe, am, ch)
    use wfl, only: fft_period
    integer, intent(in):: n
    real, intent(in):: dt, x(n)
    real, intent(out):: pe(n-1), am(n-1)
    character(len=1), optional, intent(in):: ch
    if(present(ch) .and. ch .eq. 'y') then
        call fft_period(n, dt, x, pe, am, ch)
    else
        call fft_period(n, dt, x, pe, am)
    endif
end subroutine py_fft_period

subroutine py_fns(n, dt, lag1_cor, per, pk)
    use wfl, only: fns
    integer, intent(in):: n
    real *8, intent(in):: dt, lag1_cor
    real *8, dimension(1:n/2), intent(out):: per, pk
    call fns(n, dt, lag1_cor, per, pk)
end subroutine py_fns

subroutine py_fre_wavenum(nt, dt, n, c, s, pe, fw, f, w)
    use wfl, only: fre_wavenum
    integer, intent(in):: n, nt
    real*8, intent(in):: c(0:n, 0:n, nt), s(0:n, 0:n, nt), dt
    real*8, intent(out):: fw(nt/2, 0:n), f(nt/2), w(0:n), pe(nt/2)
    call fre_wavenum(nt, dt, n, c, s, pe, fw, f, w)
end subroutine py_fre_wavenum

subroutine py_wfl_hilbert(n, x, xhil)
    use wfl, only: wfl_hilbert
    integer, intent(in)::n
    real*8, intent(in)::x(n)
    real*8, intent(out)::xhil(n)
    call wfl_hilbert(n, x, xhil)
end subroutine py_wfl_hilbert

subroutine py_hilbert_amp(n, dt, x, amp, ph, fre)
    use wfl, only: hilbert_amp
    integer, intent(in)::n
    real, intent(in)::dt, x(n)
    real, dimension(n), intent(out)::amp, ph, fre
    call hilbert_amp(n, dt, x, amp, ph, fre)
end subroutine py_hilbert_amp

subroutine py_pg(n, dt, x, mw, per_low, per_up, per_out, sp)
    use wfl, only:pg
    integer, intent(in) :: n, mw
    real *8, intent(in) :: x(n), dt, per_low, per_up
    real *8, intent(out) :: per_out(mw), sp(mw)
    call pg(n, dt, x, mw, per_low, per_up, per_out, sp)
end subroutine py_pg

subroutine py_wave_linear(n, dt, s0, step, jtot, param, sigma, y, scale, amp, coi)
    use wfl, only:wave_linear
    integer, intent(in)::n, jtot
    real, intent(in)::dt, s0, step, param, sigma, y(n)
    real, intent(out)::scale(jtot), amp(n, jtot), coi(n)
    call wave_linear(n, dt, s0, step, jtot, param, sigma, y, scale, amp, coi)
end subroutine py_wave_linear

subroutine py_wavelet_max_amp(ntime, nper, period, amp, per, per_bands, ampmax)
    use wfl, only: wavelet_max_amp
    integer, intent(in):: ntime, nper
    real, intent(in):: period(nper), amp(ntime, nper), per, per_bands
    real, intent(out):: ampmax(ntime)
    call wavelet_max_amp(ntime, nper, period, amp, per, per_bands, ampmax)
end subroutine py_wavelet_max_amp

subroutine py_wfl_arspec(n, dt, np, co, spec)
    use wfl, only: wfl_arspec
    integer, intent(in) :: n, np
    real*8, intent(in) :: dt, co(np)
    real*8, intent(out) :: spec(n/2, 2)
    call wfl_arspec(n, dt, np, co, spec)
end subroutine py_wfl_arspec

subroutine py_wfl_ceof( data0, nx, nt, nmodes, icovcor, &
    eigenvalues, eigenvectors, princomp, variance, cumvariance, &
    sqrootweight, error_eof, method)
    use wfl, only: wfl_ceof
    integer, intent(in) :: nx, nt, nmodes, icovcor
    complex*8, intent(in) :: data0(nx, nt)
    real, intent(out) :: eigenvalues(nmodes), variance(nmodes), cumvariance(nmodes)
    complex*8, intent(out) :: eigenvectors(nx, nmodes), princomp(nt, nmodes)
    real, optional, intent(in):: sqrootweight(nx)
    real, optional, intent(out):: error_eof(nmodes)
    character(len=*), optional, intent(in) :: method  !method for solving the eigenvalue system
    call wfl_ceof( data0, nx, nt, nmodes, icovcor, &
        eigenvalues, eigenvectors, princomp, variance, cumvariance, &
        sqrootweight, error_eof, method)
end subroutine py_wfl_ceof
