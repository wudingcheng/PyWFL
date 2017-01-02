program main
    use wfl
    implicit none
    real x(2), y(2), val(4), x1, y1, v
    integer i
    do i=0,3
        if(i<=1) then
            x(i) = i * 1.0
            y(i) = i * 1.0
            val(i) = i * 1.0
        else
            val(i) = i * 1.0
        endif
    enddo 
    x1 = 1.2
    y1 = 1.2
    call bilinear(x, y, val, x1, y1, v)
    print *, v
    
end program main