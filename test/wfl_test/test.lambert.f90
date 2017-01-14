program main
    use wfl
    real*8 alon,alat,cm,ps,pn,ol,x,y
    alat=36.d0
    alon=-84.d0
    cm=-81.d0
    ps=30.d0
    pn=40.d0
    ol=35.d0
    call wfl_2lam8(alon,alat,cm,ps,pn,ol,x,y)
    print *, x,y

    call wfl_invlam8(x,y,cm,ps,pn,ol,alon,alat)
    print *, alon,alat
end program main