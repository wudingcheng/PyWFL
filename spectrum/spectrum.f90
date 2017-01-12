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

subroutine py_wfl_conv(n, x, m, respns, ans, method, isconv)
    use wfl, only: wfl_conv
    integer :: n, m
!f2py integer, intent(hide) :: n, m
    complex *8, intent(in) :: x(n), respns(m)
    complex *8, dimension(n+m), intent(out) :: ans
    character(len=*), optional, intent(in) :: method
    logical, optional, intent(in) :: isconv
    call wfl_conv(n, x, m, respns, ans, method, isconv)
end subroutine py_wfl_conv

subroutine py_wfl_window(window_length, window_type, window_parameter, wfl_window_result, bandwidth)
    use wfl, only: wfl_window
    implicit none
    integer, intent(in) :: window_length
    character(len=*),intent(in) :: window_type, window_parameter
    real *8, dimension(window_length),intent(out) :: wfl_window_result
    real *8,optional,intent(out) :: bandwidth
    call wfl_window(window_length, window_type, window_parameter, wfl_window_result, bandwidth)
end subroutine py_wfl_window

!-----------------------------------------------------------------------------


subroutine py_hilbert_emd_fre(iunit,n,m,t,fre_low,fre_up,emd_imf,nfreline,filename,tlow,tup)
    use wfl, only: hilbert_emd_fre
    integer,intent(in):: iunit,n,m,nfreline
    real,intent(in):: t(n),emd_imf(n,m),fre_low,fre_up
    character(len=*),intent(in)::filename
    real,optional,intent(in):: tlow,tup
    call hilbert_emd_fre(iunit,n,m,t,fre_low,fre_up,emd_imf,nfreline,filename,tlow,tup)
end subroutine py_hilbert_emd_fre

subroutine py_hilbert_emd_per(iunit,n,m,t,per_low,per_up,emd_imf,nfreline,filename,tlow,tup)
    use wfl, only: hilbert_emd_per
    integer,intent(in):: iunit,n,m,nfreline
    real,intent(in):: t(n),emd_imf(n,m),per_low,per_up
    character(len=*),intent(in)::filename
    real,optional,intent(in):: tlow,tup
    call hilbert_emd_per(iunit,n,m,t,per_low,per_up,emd_imf,nfreline,filename,tlow,tup)
end subroutine py_hilbert_emd_per

subroutine py_hilbert_spec_fre(iunit,n,t,fre_low,fre_up,x,nfreline,filename)
    use wfl, only: hilbert_spec_fre
    integer,intent(in):: iunit,n,nfreline
    real,intent(in):: t(n),x(n),fre_low,fre_up
    character(len=*),intent(in)::filename
    call hilbert_spec_fre(iunit,n,t,fre_low,fre_up,x,nfreline,filename)
end subroutine py_hilbert_spec_fre

subroutine py_hilbert_spec_per(iunit,n,t,per_low,per_up,x,nperline,filename)
    use wfl, only: hilbert_spec_per
    integer,intent(in):: iunit,n,nperline
    real,intent(in):: t(n),x(n),per_low,per_up
    character(len=*),intent(in)::filename
    call hilbert_spec_per(iunit,n,t,per_low,per_up,x,nperline,filename)
end subroutine py_hilbert_spec_per

function py_wfl_ceof_ith_eigen(nt,nmodes,eigenvalues) result(ith)
    use wfl, only: wfl_ceof_ith_eigen
    implicit none
    integer,intent(in) :: nt,nmodes
    real*8, intent(in) :: eigenvalues(nmodes)
    integer :: ith
    ith = wfl_ceof_ith_eigen(nt,nmodes,eigenvalues)
end function py_wfl_ceof_ith_eigen

subroutine py_wfl_fft_am_window(n,dt,x,m,window,overlap,p,am_window,am,lag1_co)
    use wfl, only: wfl_fft_am_window
    implicit none
    integer,intent(in):: n,m
    real *8,intent(in):: dt,x(n),p
    character(len=*),intent(in) :: window
    logical,intent(in) :: overlap
    real *8,intent(out) :: am_window(m/2,2)
    real *8,intent(out):: am(n/2,3)
    real *8,optional,intent(in) :: lag1_co
    if(present(lag1_co) .and. lag1_co .lt. 0) then
        call wfl_fft_am_window(n,dt,x,m,window,overlap,p,am_window,am)
    else
        call wfl_fft_am_window(n,dt,x,m,window,overlap,p,am_window,am,lag1_co)
    endif
end subroutine py_wfl_fft_am_window

subroutine py_wfl_fft_omiga(n,dt,omiga)
    use wfl, only: wfl_fft_omiga
    implicit none
    integer,intent(in):: n
    real *8,intent(in):: dt
    real *8,intent(out) :: omiga(n)
    call wfl_fft_omiga(n,dt,omiga)
end subroutine py_wfl_fft_omiga

subroutine py_wfl_fft_spec_window(n,dt,x,m,window,overlap,p,spec_window,spec,lag1_co)
    use wfl, only: wfl_fft_spec_window
    implicit none
    integer,intent(in):: n,m
    real *8,intent(in):: dt,x(n),p
    character(len=*),intent(in) :: window
    logical,intent(in) :: overlap
    real *8,intent(out) :: spec_window(m/2,2)
    real *8,intent(out):: spec(n/2,3)
    real *8,optional,intent(in) :: lag1_co
    if(present(lag1_co) .and. lag1_co .lt. 0) then
        call wfl_fft_spec_window(n,dt,x,m,window,overlap,p,spec_window,spec)
    else
        call wfl_fft_spec_window(n,dt,x,m,window,overlap,p,spec_window,spec,lag1_co)
    endif
end subroutine py_wfl_fft_spec_window

subroutine py_wfl_mtm(n, dt, x, nw, x_spec, spec_type, f_test, x_edof, p, red_noise_spec, lag1_co)
    use wfl, only: wfl_mtm
    integer,intent(in) :: n
    real *8,dimension(n),intent(in) :: x
    real *8,intent(in) :: nw,dt
    real *8,dimension(n/2,2),intent(out) :: x_spec
    character(len=*),optional,intent(in) :: spec_type
    real *8,optional,dimension(n/2),intent(out) :: x_edof,f_test
    real *8,optional,dimension(n/2,3),intent(out) :: red_noise_spec
    real *8,optional,intent(in) :: p,lag1_co
    call wfl_mtm(n, dt, x, nw, x_spec, spec_type, f_test, x_edof, p, red_noise_spec, lag1_co)
end subroutine py_wfl_mtm

subroutine py_wfl_mtm_overlap( n, dt, x, m, overlap, nw, x_spec, spec_type, p, red_noise_spec, lag1_co )
    use wfl, only: wfl_mtm_overlap
    integer,intent(in) :: n,m
    real *8,intent(in):: dt,nw,x(n),overlap
    real *8,intent(out) :: x_spec(m/2,2)
    character(len=*),optional :: spec_type
    real *8,optional,intent(in) :: p,lag1_co
    real *8,optional,intent(out) :: red_noise_spec(m/2,3)
    call wfl_mtm_overlap( n, dt, x, m, overlap, nw, x_spec, spec_type, p, red_noise_spec, lag1_co )
end subroutine py_wfl_mtm_overlap

subroutine py_wfl_multicor_spec_mtm(n, m, dt, x, nw, per, cross_spec)
    use wfl, only: wfl_multicor_spec_mtm
    integer,intent(in) :: n,m
    real *8,dimension(n,m),intent(in) :: x
    real *8,intent(in) :: nw,dt
    real *8,dimension(n/2),intent(out) :: per
    complex(8),dimension(n/2,m,m),intent(out) :: cross_spec
    call wfl_multicor_spec_mtm(n, m, dt, x, nw, per, cross_spec)
end subroutine py_wfl_multicor_spec_mtm

subroutine py_wfl_multicor_spectra(n,m,dt,x,max_lag,window_type,pe,x_out,x_out2,x_out3,out_type, p, cll, clu)
    use wfl, only: wfl_multicor_spectra
    integer, intent(in) :: n,m,max_lag
    real *8,intent(in) :: dt
    real *8,dimension(n,m),intent(in) :: x
    real *8,dimension(n/2),intent(out) :: pe
    real *8,dimension(n/2,m,m),intent(out) :: x_out,x_out2,x_out3
    character(len=*),intent(in) :: out_type,window_type
    real *8,optional,intent(in) :: p
    real *8,optional,dimension(n/2,m,m),intent(out) :: cll,clu
    call wfl_multicor_spectra(n,m,dt,x,max_lag,window_type,pe,x_out,x_out2,x_out3,out_type, p, cll, clu)
end subroutine py_wfl_multicor_spectra

subroutine py_wfl_svd_coupled(n,mp,mq,s,p,kmode,lk,uk,vk,ak,bk,scfk,cumk,sk,pk)
    use wfl, only: wfl_svd_coupled
    integer,intent(in):: n,mp,mq,kmode
    real *8,intent(in):: s(n,mp),p(n,mq)
    real *8,intent(out):: lk(kmode),uk(mp,kmode),vk(mq,kmode),ak(n,kmode),bk(n,kmode),scfk(kmode),cumk(kmode)
    real *8,optional,intent(out):: sk(n,mp),pk(n,mq)
    call wfl_svd_coupled(n,mp,mq,s,p,kmode,lk,uk,vk,ak,bk,scfk,cumk,sk,pk)
end subroutine py_wfl_svd_coupled

subroutine py_wfl_wavelet(n,dt,s0,dj,jtot,y,wave,period,scale,param,coi)
    use wfl, only:wfl_wavelet
    integer,intent(in) :: n,jtot
    real ,intent(in) :: dt,s0,dj,y(n)
    real ,intent(out) :: period(jtot)
    complex(4),intent(out) :: wave(n,jtot)
    real ,optional,intent(in) :: param
    real ,optional,intent(out) :: coi(n),scale(jtot)
    call wfl_wavelet(n,dt,s0,dj,jtot,y,wave,period,scale,param,coi)
end subroutine py_wfl_wavelet

subroutine py_wfl_wavelet_coherency(nt,time,nj,scale,wave1,wave2,wave_coher,wave_phase,smooth,m,w,tolerance,global_coher,global_phase,power1,power2)
    use wfl, only: wfl_wavelet_coherency
    integer,intent(in):: nt,nj
    real ,intent(in) :: time(nt),scale(nj)
    complex(4),intent(in):: wave1(nt,nj),wave2(nt,nj)
    real ,intent(out):: wave_coher(nt,nj),wave_phase(nt,nj)
    logical,optional,intent(in) :: smooth
    integer,optional,intent(in) :: m,w
    real ,optional,intent(in):: tolerance
    real ,optional,intent(out):: global_coher(nj),global_phase(nj),power1(nt,nj),power2(nt,nj)
    call wfl_wavelet_coherency(nt,time,nj,scale,wave1,wave2,wave_coher,wave_phase,smooth,m,w,tolerance,global_coher,global_phase,power1,power2)
end subroutine py_wfl_wavelet_coherency

subroutine py_wfl_wavelet_coherency_lag(nt,time,nj,scale,wave1,wave2,wave_coher_lag,wave_phase_lag,time_idx_lag,smooth,m,w,tolerance)
    use wfl, only: wfl_wavelet_coherency_lag
    integer,intent(in):: nt,nj
    real ,intent(in) :: time(nt),scale(nj)
    complex(4),intent(in):: wave1(nt,nj),wave2(nt,nj)
    real ,intent(out):: wave_coher_lag(-nt/2:nt/2,nj),wave_phase_lag(-nt/2:nt/2,nj)
    real ,optional,intent(out):: time_idx_lag(-nt/2:nt/2)
    logical,optional,intent(in) :: smooth
    integer,optional,intent(in) :: m,w
    real ,optional,intent(in) :: tolerance
    call wfl_wavelet_coherency_lag(nt,time,nj,scale,wave1,wave2,wave_coher_lag,wave_phase_lag,time_idx_lag,smooth,m,w,tolerance)
end subroutine py_wfl_wavelet_coherency_lag

subroutine py_wfl_wavelet_signif(n,dt,y,s0,dj,jtot, &
    scale,period,dof,fft_theory,signif,param,siglvl,lag1,isigtest)
    use wfl, only:wfl_wavelet_signif
    integer,intent(in) :: n,jtot
    real ,intent(in) :: dt,y(n),s0,dj,scale(jtot),period(jtot)
    real ,intent(inout) :: dof(jtot)
    real ,intent(out) :: fft_theory(jtot),signif(jtot)
    real ,optional :: param,siglvl,lag1
    integer,optional :: isigtest
    call wfl_wavelet_signif(n,dt,y,s0,dj,jtot, &
        scale,period,dof,fft_theory,signif,param,siglvl,lag1,isigtest)
end subroutine py_wfl_wavelet_signif