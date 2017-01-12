subroutine py_factorial(n, result)
    use wfl, only: factorial
    real *8, intent(in) :: n
    real *8, intent(out) :: result
    result = factorial(n)
end subroutine py_factorial

! subroutine py_cfftpack(z,n)
!     use wfl, only: cfftpack
!     integer,intent(in):: n
!     complex,intent(inout):: z(n)
!     call cfftpack(z,n)
! end subroutine py_cfftpack

! subroutine py_cfftpackb(z,n)
!     use wfl, only: cfftpackb
!     integer,intent(in):: n
!     complex,intent(inout):: z(n)
!     call cfftpackb(z,n)
! end subroutine py_cfftpackb

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

