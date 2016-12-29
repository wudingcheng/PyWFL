module harmonics
    implicit none
    
    contains
        subroutine py_pnm(n,theta,p)
            use wfl, only: pnm
            integer,intent(in)::n
            real*8,intent(in)::theta
            real*8,intent(out):: p(0:n,0:n)
            call pnm(n, theta, p)
        end subroutine py_pnm

        subroutine py_pnmi(n,ts,tn,pinmsn)
            use wfl, only: pnmi
            integer,intent(in)::n
            real*8, intent(in)::ts,tn
            real*8, intent(out):: pinmsn(0:n,0:n)
            call pnmi(n, ts,tn,pinmsn)
        end subroutine py_pnmi

        subroutine py_wfl_pn(n,x,pnx)
            use wfl, only: wfl_pn
            integer, intent(in) :: n
            real *8, intent(in) :: x
            real *8, intent(out):: pnx(0:n)
            call wfl_pn(n, x, pnx)
        end subroutine py_wfl_pn

        subroutine py_wfl_pn_diff(n,x,pn_diff)
            use wfl, only: wfl_pn_diff
            integer, intent(in) :: n
            real *8, intent(in) :: x
            real *8, intent(out):: pn_diff(0:n)
            call wfl_pn_diff(n, x, pn_diff)
        end subroutine py_wfl_pn_diff

        ! subroutine py_wfl_pnm_derivative(n,theta,dp)
        !     use wfl
        !     integer,intent(in)::n
        !     real*8,intent(in)::theta
        !     real*8,intent(out):: dp(0:n,0:n)
        !     call wfl_pnm_derivative(n, theta, dp)
        ! end subroutine py_wfl_pnm_derivative

        ! function py_cs2kai(cs,cha)  result(kai)
        !     use wfl
        !     implicit none
        !     real*8,intent(in) :: cs
        !     character(len=1),intent(in):: cha
        !     real*8 kai
        !     kai = cs2kai(cs, cha)
        ! end function py_cs2kai

        subroutine py_cs2k(n,c,s,k,cha)
            use wfl, only: cs2k
            integer, intent(in):: n
            real *8, intent(in),dimension(0:n,0:n) :: c,s
            complex*16, intent(out) :: k(0:n,-n:n)
            character(len=7),optional, intent(in)  :: cha
            if(present(cha) .and. (cha .eq. 'quantum')) then
                call cs2k(n, c, s, k, cha)
            else
                call cs2k(n, c, s, k)
            endif
        end subroutine py_cs2k

        subroutine py_wfl_func_exp_fft(lon,lat,func,cord_lon,cord_lat,n,c,s,cha,diff_percent,iterative)
            use wfl, only: wfl_func_exp_fft 
            integer, intent(in) :: lon,lat,n
            real*8, intent(in) :: func(lon,lat-1),cord_lon(lon),cord_lat(lat)
            real*8, intent(out) :: c(0:n,0:n),s(0:n,0:n)
            character*6, optional, intent(in) :: cha
            real*8, optional,intent(out) :: diff_percent
            integer,optional,intent(in) :: iterative
            if(present(cha) .and. (cha .eq. 'smooth')) then
                call wfl_func_exp_fft(lon,lat,func,cord_lon,cord_lat,n,c,s,cha,diff_percent,iterative)
            else
                call wfl_func_exp_fft(lon,lat,func,cord_lon,cord_lat,n,c,s,diff_percent=diff_percent,iterative=iterative)
            endif
        end subroutine py_wfl_func_exp_fft

        subroutine py_wfl_dpnm(n, theta, dp)
            use wfl, only: wfl_dpnm
            integer, intent(in) :: n
            real *8, intent(in) :: theta
            real *8, intent(out):: dp(0:n, 0:n)
            call wfl_dpnm(n, theta, dp)
        end subroutine py_wfl_dpnm
    
        subroutine py_wfl_dpnm2(n, theta, dp, ddp)
            use wfl, only: wfl_dpnm
            integer, intent(in) :: n
            real *8, intent(in) :: theta
            real *8, intent(out):: dp(0:n, 0:n)
            real *8, intent(out):: ddp(0:n, 0:n)
            call wfl_dpnm(n, theta, dp, ddp)
        end subroutine py_wfl_dpnm2

end module harmonics
