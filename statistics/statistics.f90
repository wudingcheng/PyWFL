subroutine py_aic(m_order_max,m,rsdus,rs,m_order)
    use wfl, only:aic
    integer,intent(in):: m_order_max,m
    real*8,intent(in):: rsdus(m_order_max)
    real*8,intent(out)::rs(m_order_max)
    integer,intent(out):: m_order
    call aic(m_order_max,m,rsdus,rs,m_order)
end subroutine py_aic

! function py_arth(first,increment,n) result(arth_r)
!     use wfl, only:arth
!     real *8, intent(in) :: first,increment
!     integer, intent(in) :: n
!     real *8, dimension(n) :: arth_r
!     arth_r = arth(first,increment,n)
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

recursive function py_cumsum(arr,seed) result(ans)
    use wfl, only: cumsum
    real*8, dimension(:), intent(in) :: arr
    real*8, optional, intent(in) :: seed
    real*8, dimension(size(arr)) :: ans
    ans=cumsum(arr,seed)
end function py_cumsum

subroutine py_diff(n, xin, xout)
    use wfl, only: diff
    integer,intent(in):: n
    real,intent(in):: xin(n)
    real,intent(out):: xout(n-1)
    call diff(n, xin, xout)
end subroutine py_diff

function py_geop(first,factor,n) result(res)
    use wfl, only:geop
    real *8, intent(in) :: first,factor
    integer, intent(in) :: n
    real *8, dimension(n) :: res
    res = geop(first, factor, n)
    end function py_geop

subroutine py_leap_rm(n, x, zleap_criterion, xout)
    use wfl, only:leap_rm
    integer,intent(in):: n
    real,intent(in):: x(n),zleap_criterion
    real,intent(out):: xout(n)
    call leap_rm(n, x, zleap_criterion, xout)
end subroutine py_leap_rm

! function py_mean(x, dim, mask) result(res)
!     use wfl, only: mean
!     real *8, intent(in) :: x
!     integer, optional, intent(in) :: dim
!     logical, dimension(size(x)), optional, intent(in) :: mask
!     real, dimension(size(x,dim=dim)) :: res
!     res = mean(x, dim, mask)
! end function py_mean

function py_monotonic(x) result(iflag)
    use wfl, only: monotonic
    real,dimension(:),intent(in):: x
    integer iflag
    iflag = monotonic(x)
end function py_monotonic

subroutine py_montecarlo(n,ivars,co90,co90s,co95,co95s,co99,co99s)
    use wfl, only: montecarlo, montecarlomulti
    integer,intent(in)::n
    real,intent(out)::co90,co90s,co95,co95s,co99,co99s
    integer, optional, intent(in) :: ivars
    if(present(ivars) .and. ivars .gt. 2) then
        call montecarlomulti(n,ivars,co90,co90s,co95,co95s,co99,co99s)
    else
        call montecarlo(n,co90,co90s,co95,co95s,co99,co99s)
    endif
end subroutine py_montecarlo

subroutine py_month_mean(n,x,yyyy,mm,m,xout)
    use wfl, only: month_mean
    integer,intent(in):: n,yyyy,mm,m
    real,intent(in):: x(n)
    real,intent(out):: xout(m)
    call month_mean(n,x,yyyy,mm,m,xout)
end subroutine py_month_mean

subroutine py_month_mean_idx(n,t,nmon,iy_begin,im_begin,iy_end,im_end,idx_begin,idx_end)
    use wfl, only: month_mean_idx
    integer,intent(in):: n,nmon,iy_begin,im_begin,iy_end,im_end
    real,intent(in) :: t(n)
    integer,intent(out):: idx_begin(nmon),idx_end(nmon)
    call month_mean_idx(n,t,nmon,iy_begin,im_begin,iy_end,im_end,idx_begin,idx_end)
end subroutine py_month_mean_idx

function py_varpercent(a, b) result (var_percent)
    use wfl, only: varpercent
    real *8, dimension(:):: a,b
    real *8 var_percent
    var_percent = varpercent(a, b)
end function py_varpercent

function py_wfl_cdff(p,dof1,dof2) result(cdff_result)
    use wfl, only: wfl_cdff
    real*8,intent(in):: p,dof1,dof2
    real*8 cdff_result
    cdff_result = wfl_cdff(p,dof1,dof2)
end function py_wfl_cdff

function py_wfl_cdfnor(p) result(cdfnor_result)
    use wfl, only: wfl_cdfnor
    real*8,intent(in):: p
    real*8 cdfnor_result
    cdfnor_result = wfl_cdfnor(p)
end function py_wfl_cdfnor

function py_wfl_cdft(p,dof) result(cdft_result)
    use wfl, only: wfl_cdft
    real*8,intent(in):: p,dof
    real*8 cdft_result
    cdft_result = wfl_cdft(p,dof)
end function py_wfl_cdft

subroutine py_wfl_taylor_diagram(n,r,f,out,isnormal)
    use wfl, only: wfl_taylor_diagram
    integer,intent(in) :: n
    real *8,dimension(n),intent(in):: r, f
    real *8,dimension(6),intent(out) :: out
    character(len=1),optional,intent(in) :: isnormal
    call wfl_taylor_diagram(n,r,f,out,isnormal)
end subroutine py_wfl_taylor_diagram