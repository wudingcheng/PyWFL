program main
    use wfl
    parameter(n=1024,m_per=2,idef=2,idef2=0)
    real*8 t(n),x(n),zper(n,m_per+3),pdef(n,idef),co(idef,2),pdef2(n,idef2),co2(idef2,2),poly_co(2,2)
    real*8,dimension(m_per):: am,pe,dam,dph,ph
    real*8 tmp(n)
    integer itrend,timebegin,n_order
    pe(1)=365.2422; pe(2)=182.6211

    print *, 'test 1 ... '
    print *, 'full model !!!'

    !first form the model use true data
    do i=1,n
        t(i)=51544.0d0+dble(i-1)*5.d0
        dt=dble(i-1)*5.0d0
        x(i)=435.5d0+0.5d0*dt+1.0*dcos(2.d0*pi_dp*dt/pe(1)-60.d0*pi_dp/180.d0)+1.5d0*dcos(2.d0*pi_dp*dt/pe(2)-135.d0*pi_dp/180.d0)+3.7d0*sind(dble(i))+51.5d0*(sind(dble(i))+cosd(dble(2*i)) )
        pdef(i,1)=sind(dble(i))
        pdef(i,2)=sind(dble(i))+cosd(dble(2*i))
    enddo

    !call wfl_harmonic_more(n,t,x,m_per,pe,idef,pdef,am,ph,dam,dph,co,zper,timebegin=20000101,method='aic')
    !to get the coefficients using least squares fit method
    call wfl_harmonic_more(n,t,x,m_per,pe,idef,pdef,am,ph,dam,dph,co,zper, &
        & itrend=2,timebegin=20000101,best_fit_order=n_order,poly_co=poly_co)
    !after the above call, you will find the polynomial coefficients on the screen
    print *, 'harmonic amplitude and phase:'
    print *, am,ph !for harmonic amplitude and phase
    print *, 'coefficients of user defined parameters: '
    print *, co(1,1),co(1,2),co(2,1),co(2,2) !for user defined parameters's coefficients and std of coefficients


    !test for no polynomial fit and no harmonics
    print *, ''
    print *, 'test 2 ...'
    print *, 'only user defined parameters !!!'
    print *, ''
    print *, 'coefficients of user defined parameters: '
    print *, co

    do i=1,n
        t(i)=51544.0d0+dble(i-1)*5.d0
        dt=dble(i-1)*5.0d0
        x(i)=3.7d0*sind(dble(i))+51.5d0*(sind(dble(i))+cosd(dble(2*i)) )
        pdef(i,1)=sind(dble(i))
        pdef(i,2)=sind(dble(i))+cosd(dble(2*i))
    enddo

    call wfl_harmonic_more(n,t,x,0,pe,idef,pdef,am,ph,dam,dph,co,zper, &
        & itrend=0,timebegin=20000101,best_fit_order=n_order, poly_co=poly_co)


    !test for no user defined parameters, also general ls fits
    print *, ''
    print *, 'test 3 ...'
    print *, 'polynomial + harmonic model !!!'
    do i=1,n
        t(i)=51544.0d0+dble(i-1)*5.d0
        dt=dble(i-1)*5.0d0
        x(i)=435.5d0+0.5d0*dt+1.0*dcos(2.*pi*dt/pe(1)-60.0*pi/180.)+1.5*dcos(2.*pi*dt/pe(2)-135.0*pi/180.) !+3.7d0*sind(dble(i))+51.5d0*(sind(dble(i))+cosd(dble(2*i)) )
        pdef(i,1)=sind(dble(i))
        pdef(i,2)=sind(dble(i))+cosd(dble(2*i))
    enddo

    call wfl_harmonic_more(n,t,x,m_per,pe,idef2,pdef2,am,ph,dam,dph,co2,zper, &
        & itrend=2,timebegin=20000101,best_fit_order=n_order, poly_co=poly_co)
    print *, 'harmonic amplitude and phase:'
    print *, am,ph

end program main
