program main
    use wfl
    integer n
    parameter(n=64)
    integer i,j
    real x(n,2),x_out(0:n-1,2,2)
    real*8 out(-n:n-1,2,2),y(n,2)

    x=0.
    x(1:n,1)=[47,64,23,71,38,65,55,41,59,48,71,35,56,40,58,44,80,55,&
             37,74,51,58,50,60,44,57,50,45,25,59,50,71,56,74,50,58,&
             45,54,36,54,48,55,45,57,50,62,44,64,43,52,38,60,55,41,&
             53,49,34,35,54,45,68,38,50,60]
    x(1:n,2)=x(1:n,1)
    print *, x
    call cpu_time(s1)
    call wfl_multicor (n, 2, x, x_out,'cor',remove_mean='y' )
    call cpu_time(s2)
    y=dble(x)
    call wfl_multicor_fft(n, 2, y, out, 'cor',remove_mean='y' )
    call cpu_time(s3)

    print *, 'time: ', s2-s1,s3-s2
    print *, 'same out'
    do i=0,n-1
        write(*,'(1x,i3,3x,f20.6,2f15.6)') i,out(i,1,2),x_out(i,1,2)
    enddo

    aa=0.
end program main