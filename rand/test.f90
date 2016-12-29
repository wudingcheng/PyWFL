program main
    use wfl
    implicit none
    integer,parameter:: n=2048,m=3
    real,dimension(n):: x
    real:: beta(m)
    integer:: piece(m-1), iseed

    beta=[0.5,1.,2.]
    piece=[100,300]
    iseed=12345
    call wfl_noise_multi_color(n, m, beta, piece, x, iseed)
    print *, x
end program main