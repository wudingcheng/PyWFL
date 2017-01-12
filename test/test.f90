subroutine pass_func(x, m, func, xout)
    integer, intent(in) :: m
    real, intent(in) :: x(m)
    real, intent(out) :: xout(m)
    real temp(m)
    external func

    print *, 'x=', x
    call func(x, m, temp)
    print *, temp
    xout(:) = temp(:) + 1.0
    print *, xout
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

    integer :: i,j,k,l,m,mfit
    real*8 :: sum,wt,ym,afunc(ma),beta(ma),sig2i,sig0(ndat)
    real*8 :: afunc0(ma,ndat)

    mfit=ma

    if(mfit.eq.0) then
        print *,  'error: wfl_lsfit: no parameters to be fitted'
        return
    endif

    covar=0.d0
    beta=0.d0

    sig0=1.d0
    if(present(sig) .and. sig(0) .gt. 0) sig0=sig

    call funcs(ndat,x,ma,afunc0)
    do i=1, ndat !loop over data to accumulate coefficients of the normal
        afunc(:)=afunc0(:,i)
        ym=y(i)
        sig2i=1.d0/sig0(i)**2
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

subroutine lstsq(H, ndat, y, ma, a, covar, chisq, sig)
    use wfl, only: gaussj8, covsrt8
    implicit none
    integer,intent(in) :: ndat,ma
    real*8,intent(in),dimension(ndat) :: y
    real*8, intent(in), dimension(ma, ndat) :: H
    real*8,intent(out) :: a(ma),covar(ma,ma)
    real*8,optional,intent(out) :: chisq
    real*8,optional,intent(in) :: sig(ndat)


    integer :: i,j,k,l,m,mfit
    real*8 :: sum,wt,ym,afunc(ma),beta(ma),sig2i,sig0(ndat)
    real*8 :: afunc0(ma,ndat)

    mfit=ma

    if(mfit.eq.0) then
        print *,  'error: wfl_lstsq: no parameters to be fitted'
    endif

    covar=0.d0
    beta=0.d0

    sig0=1.d0
    if(present(sig)) sig0=sig

    ! call funcs(ndat,x,ma,afunc0)
    ! print *, 'afunc0', afunc0
    do i=1,ndat !loop over data to accumulate coefficients of the normal
    !call funcs(x(i),ma,afunc) !equations.
        afunc(:)=H(:,i)
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
    print *, covar
    call gaussj8(covar,mfit,ma,beta,1,1) !matrix solution. defined in fitpoly.f90


    a(1:ma)=beta(1:ma)  !partition solution to appropriate coefficients a.

    if(present(chisq)) then
        chisq=0. !evaluate ¦ö2  of the fit.
        do i=1,ndat
        !call funcs(x(i),ma,afunc)
        afunc(:)=H(:,i)
        sum=0.
            do j=1,ma
            sum=sum+a(j)*afunc(j)
            enddo
        chisq=chisq+((y(i)-sum)/sig0(i))*((y(i)-sum)/sig0(i))
        enddo
    endif

    call covsrt8(covar,ma,mfit) !sort covariance matrix to true order of fitting,defined in fitpoly.f90

end subroutine lstsq



