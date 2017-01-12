program main
    use wfl
    integer :: n
    integer,parameter:: m=10
    real*8,allocatable,dimension(:,:) :: p
    integer,dimension(m) :: i_n
    real*8,dimension(m) :: res
    real*8,dimension(0:180*12+1) :: res1,res0(7300)
    real*8 s1,s2, theta

    i_n(1:m)=[360,720,1080,1440,1800,1900,2000,2100,3600,7200]

    write(*,*) 'n(degree) max_res time_used(s)'
    do j=1,m
    call cpu_time(s1)
    ! print *, 'n=',i_n(j)
    n=i_n(j)
    allocate(p(0:n,0:n) )

    res0=0.d0; res1=0.d0
    do i=0,180*12+1
        theta=dble(i-1)/12.d0
    p=0.0d0
    call pnm( n, theta, p )

    kk=0
    do in=0,n
    kk=kk+1
    res0(kk)=abs( sum(p(in,0:in)*p(in,0:in))-2.d0*dble(in)-1.0d0 )
    enddo

    res1(i)=maxval( res0(1:kk) )
    enddo

    res(j)=maxval(res1(0:180*12+1) )
    deallocate(p)
    call cpu_time(s2)
    write(*,'(i10,d15.7,f15.0)') n,res(j),s2-s1

    enddo

end program main