program main
    use wfl
    implicit none
    real :: a(5), ans(5)
    integer :: i
    do i = 1, 5
        a(i) = i
    enddo
    write(*, *) a
    call wfl_running_mean(5, a, 3, ans)
    write(*, *) ans

    
end program main