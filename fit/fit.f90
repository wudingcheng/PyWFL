subroutine py_fitline(x,y,ndata,a,b,siga,sigb)
    use wfl, only: fitline
    integer, intent(in):: ndata
    real *8, intent(in):: x(ndata),y(ndata)
    real *8, intent(out):: a,b,siga,sigb
    call fitline(x,y,ndata,a,b,siga,sigb)
end subroutine py_fitline

subroutine py_fitmulti(x,y,ndat,ma,a,covar,chisq,s,r,v)
    use wfl, only: fitmulti
    integer*4, intent(in) :: ndat,ma
    real*8,intent(in) :: x(ndat,ma-1),y(ndat)
    real*8,intent(out) :: a(ma),covar(ma,ma),chisq,s,r,v(ma)
    call fitmulti(x,y,ndat,ma,a,covar,chisq,s,r,v)
end subroutine py_fitmulti

subroutine py_fitpoly(x,y,ndat,ma,a,covar,chisq)
    use wfl, only:fitpoly
    integer,intent(in):: ma,ndat
    real*8,intent(in):: x(ndat),y(ndat)
    real*8,intent(out):: chisq,a(ma),covar(ma,ma)
    call fitpoly(x,y,ndat,ma,a,covar,chisq)
end subroutine py_fitpoly

function py_fitpoly_best_degree(x, y, n, max_order, ch )  result(best_degree)
    use wfl, only:fitpoly_best_degree
    integer,intent(in):: n,max_order
    real*8,intent(in):: x(n),y(n)
    character(len=*),optional,intent(in):: ch
    integer best_degree
    best_degree = fitpoly_best_degree(x, y, n, max_order, ch )
end function py_fitpoly_best_degree

subroutine py_detrend(n,x,y,yout)
    use wfl, only:detrend
    integer,intent(in):: n
    real*8,intent(in):: x(n),y(n)
    real*8,intent(out):: yout(n)
    call detrend(n,x,y,yout)
end subroutine py_detrend

subroutine py_wfl_harmonic_fix(n, t, x, m_per, pe, amph, zper, itrend, timebegin, poly_co)
    use wfl, only:wfl_harmonic
    integer,intent(in) :: n, m_per, itrend
    real *8,intent(in) :: t(n),x(n),pe(m_per)
    real *8,intent(out) :: amph(4,m_per),zper(n,m_per+3)
    integer,optional,intent(in) :: timebegin
    real *8,optional,dimension(:,:),intent(out) :: poly_co
    if(present(timebegin) .and. timebegin .ne. x(0)) then
        call wfl_harmonic(n, t, x, m_per, pe, amph, zper, itrend, poly_co=poly_co)
    else
        call wfl_harmonic(n, t, x, m_per, pe, amph, zper, itrend, timebegin=timebegin, poly_co=poly_co)
    endif
end subroutine py_wfl_harmonic_fix

subroutine py_wfl_harmonic_polys(n, t, x, m_per, pe, amph, zper, timebegin, method, best_fit_order, poly_co)
    use wfl, only:wfl_harmonic
    integer,intent(in) :: n,m_per
    real *8,intent(in) :: t(n),x(n),pe(m_per)
    real *8,intent(out) :: amph(4,m_per),zper(n,m_per+3)
    integer,optional,intent(in) :: timebegin
    integer,intent(out) :: best_fit_order
    character(len=3),optional,intent(in) :: method
    real *8,optional,dimension(:,:),intent(out) :: poly_co
    if(present(timebegin) .and. timebegin .ne. x(0)) then
        call wfl_harmonic(n=n, t=t, x=x, m_per=m_per, pe=pe, amph=amph, &
            & zper=zper, method=method, best_fit_order=best_fit_order, &
            & poly_co=poly_co)
    else
        call wfl_harmonic(n=n, t=t, x=x, m_per=m_per, pe=pe, amph=amph, &
            & zper=zper, timebegin=timebegin, method=method, &
            & best_fit_order=best_fit_order, poly_co=poly_co)
    endif
end subroutine py_wfl_harmonic_polys

subroutine py_polynomial(x,y,n,m,yout)
    use wfl, only: polynomial
    integer,intent(in):: n,m
    real*8,dimension(n),intent(in):: x,y
    real*8,dimension(n),intent(out):: yout
    call polynomial(x,y,n,m,yout)
end subroutine py_polynomial

subroutine py_quadratic(n,x,y,yout)
    use wfl, only:quadratic
    integer,intent(in):: n
    real*8,intent(in):: x(n),y(n)
    real*8,intent(out):: yout(n)
    call quadratic(n,x,y,yout)
end subroutine py_quadratic

subroutine py_wfl_depoly(n, t, x, xout, m_order, max_order, method, nbest_fit_order)
    use wfl, only:wfl_depoly
    integer,intent(in) :: n
    real *8,dimension(n),intent(in) :: t,x
    real *8,dimension(n),intent(out) :: xout
    integer,optional,intent(in) :: m_order,max_order
    character(len=*),optional,intent(in) :: method
    integer,optional,intent(out) :: nbest_fit_order
    call wfl_depoly(n, t, x, xout, m_order, max_order, method, nbest_fit_order)
end subroutine py_wfl_depoly

subroutine py_wfl_harmonic_more(n,t,x,m_per,pe,idef,pdef,am,ph,dam,dph,co_self,zper,itrend,timebegin,method,best_fit_order, poly_co)
    use wfl, only:wfl_harmonic_more
    integer,intent(in):: n,m_per,idef
    real*8,intent(in):: t(n),x(n),pe(m_per),pdef(n,idef)
    real*8,dimension(m_per),intent(out):: am,ph,dam,dph
    real*8,intent(out):: zper(n,m_per+4),co_self(idef,2)
    integer,optional,intent(in):: itrend,timebegin   ! timebegin=20000101
    character(len=*),optional,intent(in):: method
    integer,optional,intent(out) :: best_fit_order
    real*8,optional,dimension(:,:),intent(out) :: poly_co
    call wfl_harmonic_more(n,t,x,m_per,pe,idef,pdef,am,ph,dam,dph,co_self,zper,itrend,timebegin,method,best_fit_order, poly_co)
end subroutine py_wfl_harmonic_more


! subroutine py_wfl_lsfit(ndat,x,y,ma,a,covar,funcs,chisq,sig)
!     use wfl, only:wfl_lsfit
!     integer,intent(in) :: ndat,ma
!     real *8,intent(in),dimension(ndat) :: x,y
!     real *8,intent(out) :: a(ma),covar(ma,ma)
!     real *8,optional,intent(out) :: chisq
!     real *8,optional,intent(in) :: sig(ndat)
!     external funcs
!     call wfl_lsfit(ndat,x,y,ma,a,covar,funcs,chisq,sig)
! end subroutine py_wfl_lsfit

subroutine py_wfl_miss_data_recover(n,t,x,m_per,pe,n_miss,nt_miss,xout,itrend)
    use wfl, only:wfl_miss_data_recover
    integer,intent(in) ::  n,m_per,n_miss
    real *8,dimension(n),intent(in) :: t,x,pe(m_per),nt_miss(n_miss,2)
    real *8,dimension(n),intent(out) :: xout
    integer,optional,intent(in) :: itrend
    call wfl_miss_data_recover(n,t,x,m_per,pe,n_miss,nt_miss,xout,itrend)
end subroutine py_wfl_miss_data_recover

subroutine py_wfl_rm_jump( n, t, x, m_per, pe, n_jump, nt_jump, xout, itrend, amph)
    use wfl, only:wfl_rm_jump
    integer,intent(in) :: n,n_jump,m_per
    real *8,dimension(n),intent(in) :: t,x
    real *8,dimension(m_per),intent(in) :: pe
    integer,dimension(n_jump),intent(in) :: nt_jump
    real *8,dimension(n,6),intent(out) :: xout
    integer,optional,intent(in) :: itrend
    real *8,optional,dimension(4,m_per),intent(out):: amph
    call wfl_rm_jump( n, t, x, m_per, pe, n_jump, nt_jump, xout, itrend, amph)
end subroutine py_wfl_rm_jump



