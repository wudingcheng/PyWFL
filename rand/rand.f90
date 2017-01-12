subroutine  py_wfl_random_mt(n, x, iseed, ran_a, ran_b, ran_type )
    use wfl, only: wfl_random_mt
    integer,intent(in) :: n
    real*8,dimension(n),intent(out) :: x
    integer,optional,intent(in) :: iseed,ran_type
    real*8,optional,intent(in) :: ran_a,ran_b
    call wfl_random_mt(n, x, iseed, ran_a, ran_b, ran_type )
end subroutine  py_wfl_random_mt

subroutine  py_wfl_random_mt_norm(n, x, iseed, ran_mean,ran_sig, ran_type )
    use wfl, only: wfl_random_mt_norm
    integer,intent(in) :: n
    real*8,dimension(n),intent(out) :: x
    integer,optional,intent(in) :: iseed,ran_type
    real*8,optional,intent(in) :: ran_mean,ran_sig
    call wfl_random_mt_norm(n, x, iseed, ran_mean,ran_sig, ran_type )
end subroutine  py_wfl_random_mt_norm

subroutine py_wfl_random_norm(n, x, iseed, ran_mean, ran_sig, ran_type)
    use wfl, only: wfl_random_norm
    integer, intent(in) :: n
    real, intent(out) :: x(n)
    integer, optional, intent(in) :: iseed
    real, optional, intent(in) :: ran_mean, ran_sig
    integer, optional, intent(in) :: ran_type
    call wfl_random_norm(n, x, iseed, ran_mean, ran_sig, ran_type)
end subroutine py_wfl_random_norm

! subroutine py_wfl_random_uni(n, x, iseed, ran_a, ran_b)
!     use wfl, only: wfl_random_uni
!     integer,intent(in) :: n
!     real*8,intent(out) :: x(n)
!     real*8,optional,intent(in) :: ran_a,ran_b
!     call wfl_random_uni(n, x, iseed, ran_a, ran_b)
! end subroutine py_wfl_random_uni

subroutine  py_wfl_noise_color(n, beta, x, iseed)
    use wfl, only: wfl_noise_color
    integer,intent(in) :: n
    real*8,intent(in) :: beta !0<=beta<=3
    real*8,dimension(n),intent(out) :: x
    integer,optional,intent(in) :: iseed
    if(present(iseed) .and. iseed .ne. 0) then
        call wfl_noise_color(n, beta, x, iseed)
    else
        call wfl_noise_color(n, beta, x)
    endif
end subroutine  py_wfl_noise_color

subroutine  py_wfl_noise_multi_color(n, m, beta, pieces, x, iseed)
    use wfl, only: wfl_noise_multi_color
    integer,intent(in) :: n,m
    real *8,dimension(m),intent(in) :: beta !0<=beta<=3
    integer,dimension(m-1),intent(in) :: pieces
    real *8,dimension(n),intent(out) :: x
    integer,optional,intent(in) :: iseed
    if(present(iseed) .and. iseed .ne. 0) then
        call wfl_noise_multi_color(n, m, beta, pieces, x, iseed)
    else
        call wfl_noise_multi_color(n, m, beta, pieces, x)
    endif
end subroutine  py_wfl_noise_multi_color