subroutine func(x, m, xout)
    real,intent(in) :: x
    integer,intent(in) ::  m
    real, intent(out) :: xout(m)
    !functions defined by user
    xout(1)=1.0
    do i=2,m
        xout(i)=xout(i-1)*x
    enddo
end subroutine func

subroutine pass_func(x, m, func, xout)
    integer, intent(in) :: m
    real, intent(in) :: x(m)
    real, intent(out) :: xout(m)
    ! integer i
    external func
    ! in the loop, f2py cannot handle this, and return default 0
    ! do i = 0, m
    !     call func(x(i), xout(i))
    ! enddo

    ! call this directly, the result is true
    call func(x, m, xout)
end subroutine pass_func


subroutine lsfit(ndat, x, y, ma, a, covar, funcs, chisq, sig)
    use wfl, only: gaussj8, covsrt8
    implicit none
    integer,intent(in) :: ndat,ma
    real*8,intent(in),dimension(ndat) :: x,y
    real*8,intent(out) :: a(ma),covar(ma,ma)
    real*8,optional,intent(out) :: chisq
    real*8,optional,intent(in) :: sig(ndat)
    external funcs

    ! uses covsrt,gaussj  !see "fitpoly.f90" for these subroutines

    integer :: i,j,k,l,m,mfit
    real*8 :: sum,wt,ym,afunc(ma),beta(ma),sig2i,sig0(ndat)
    real*8 :: afunc0(ma,ndat)
    !integer :: ia(ma)
    !real*8 :: sig(ndat)
    !npc=ma

    !mfit=0
    !
    !!sig=1.0d0 !initialize the weight variance to be equal
    !ia=1.0d0  !initialize all input perameters to be fitted

    mfit=ma

    if(mfit.eq.0) then
    print *,  'error: wfl_lsfit: no parameters to be fitted'
    endif

    covar=0.d0
    beta=0.d0

    sig0=1.d0
    if(present(sig)) sig0=sig

    call funcs(ndat,x,ma,afunc0)
    do i=1,ndat !loop over data to accumulate coefficients of the normal
    !call funcs(x(i),ma,afunc) !equations.
    afunc(:)=afunc0(:,i)
    ym=y(i)
    sig2i=1.d0/sig0(i)/sig0(i)
    j=0
        do l=1,ma
        j=j+1
        wt=afunc(l)*sig2i
        k=0
            do m=1,l
                k=k+1
                covar(j,k)=covar(j,k)+wt*afunc(m)
            enddo
        beta(j)=beta(j)+ym*wt
        enddo
    enddo


    do j=2,mfit !fill in above the diagonal from symmetry.
        do k=1,j-1
        covar(k,j)=covar(j,k)
        enddo
    enddo


    call gaussj8(covar,mfit,ma,beta,1,1) !matrix solution. defined in fitpoly.f90


    a(1:ma)=beta(1:ma)  !partition solution to appropriate coefficients a.

    if(present(chisq)) then
        chisq=0. !evaluate ¦ö2  of the fit.
        do i=1,ndat
        !call funcs(x(i),ma,afunc)
        afunc(:)=afunc0(:,i)
        sum=0.
            do j=1,ma
            sum=sum+a(j)*afunc(j)
            enddo
        chisq=chisq+((y(i)-sum)/sig0(i))*((y(i)-sum)/sig0(i))
        enddo
    endif

    call covsrt8(covar,ma,mfit) !sort covariance matrix to true order of fitting,defined in fitpoly.f90

end subroutine lsfit

subroutine fpoly6(x,n,p)

    integer,intent(in) :: n
    real*8,intent(in) :: x
    real*8,intent(out) :: p(n)
    !Fitting routine for a polynomial of degree n-1, with n coe.cients.
    integer j

    p(1)=1.
    do j=2,n
        p(j)=p(j-1)*x
    enddo

end subroutine fpoly6


! program main
!     ! implicit none
!     use wfl
!     integer,parameter:: n=1000,norder=3
!     real*8 :: x(n), y(n)
!     real*8 :: a(norder), co(norder,norder)
!     real*8 chisq

!     do i=1,n
!         x(i)=real(i)/1000.d0
!         y(i)=6.7d0+0.00045*x(i)+0.00000394*x(i)*x(i)
!     enddo
!     call wfl_lsfit(n,x,y,norder,a,co,fpoly6,chisq=chisq)

! end program main