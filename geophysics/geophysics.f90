function py_am2eam(am, cha) result(eam)
    use wfl
    real *8 , intent(in):: am
    character(len=*)    :: cha
    real *8 :: eam
    eam = am2eam(am, cha)
end function py_am2eam

subroutine geo_distance(phi_1, lat_1, phi_2, lat_2, ch, dis)
    use wfl
    implicit none
    real, intent(in) :: phi_1, lat_1, phi_2, lat_2
    character(8), intent(in) :: ch
    real, intent(out) :: dis
    if(ch .eq. 'angle') then
        dis = angle_dis(phi_1, lat_1, phi_2, lat_2)
    else
        dis = angle_dis(phi_1, lat_1, phi_2, lat_2, ch)
    endif
end subroutine geo_distance

function py_area(dlon, lat1, lat2) result(area_result)
    use wfl
    real *8, intent(in):: dlon, lat1, lat2
    real *8 :: area_result
    area_result = area(dlon, lat1, lat2)
end function py_area

function py_cs2kai(cs, cha)  result(kai)
    use wfl
    implicit none
    real *8, intent(in) :: cs
    character(len=1), intent(in):: cha
    real *8 kai
    kai = cs2kai(cs, cha)
end function py_cs2kai

function py_kai2cs(kai, cha)  result(cs)
    use wfl
    implicit none
    real *8, intent(in) :: kai
    character(len=1), intent(in):: cha
    real *8 cs
    cs = kai2cs(kai, cha)
end function py_kai2cs

subroutine py_degree_var(n, c, s, degvar)
    use wfl
    integer, intent(in):: n
    real *8, dimension(0:n, 0:n), intent(in):: c, s
    real *8, dimension(0:n), intent(out):: degvar
    call degree_var(n, c, s, degvar)
end subroutine py_degree_var

subroutine j2j3(colat, dlat, dlon, hi, rho, dj2, dj3)
    use wfl
    real *8, intent(in):: colat, dlat, dlon, hi, rho
    real *8, intent(out):: dj2, dj3
    call dj2j3(colat, dlat, dlon, hi, rho, dj2, dj3)
end subroutine j2j3

subroutine py_eaam_xyz_j2j3(lon, lat, lay1, ltrop, slp, pl, u, v, z, hm, &
          fp, fpn, f3p, f3pn, fw, fwl, f3w, f3wl, pmn, pm, xypn, xyp, zpn, zp, &
          j2n, j2, j3n, j3, method)
    use wfl
    integer, intent(in):: lon, lat, lay1, ltrop
    real *4, dimension(lon, lat+1, lay1), intent(inout)::z, u, v, slp(lon, lat+1)
    real *4, intent(in):: hm(lon, lat), pl(lay1)
    complex *8, dimension(lon, lat), intent(out):: fpn, fp, fw, fwl, xyp, xypn
    real *8 , dimension(lon, lat), intent(out):: f3pn, f3p, f3w, f3wl, zp, zpn, pm, pmn, &
              j2, j3, j2n, j3n
    character(len=*), intent(in) :: method
    call eaam_xyz_j2j3(lon, lat, lay1, ltrop, slp, pl, u, v, z, hm, &
          fp, fpn, f3p, f3pn, fw, fwl, f3w, f3wl, pmn, pm, xypn, xyp, zpn, zp, &
          j2n, j2, j3n, j3, method)
end subroutine py_eaam_xyz_j2j3

subroutine py_load_vert_def(nlon, nlat, sigma, u, depart_lon, zlat)
    !purpose: get the radial elastic deformation of the earth's surface estimated from a given surface load globally.
    use wfl
    implicit none
    integer, intent(in):: nlon, nlat
    real *8, intent(in):: sigma(nlon, nlat) !unit: kg/m/m
    real *8, intent(out):: u(nlon, nlat)
    real *8, optional, intent(in):: depart_lon, zlat(nlat)  !must be presented at the same time
    call load_vert_def(nlon, nlat, sigma, u, depart_lon, zlat)
end subroutine py_load_vert_def

subroutine py_load_vert_exp(p, nlon, nlat, u, pbar)
    use wfl
    implicit none
    integer, intent(in):: nlon, nlat
    real *8, dimension(nlon, nlat), intent(in):: p
    real *8, dimension(nlon, nlat), intent(out):: u, pbar
    call load_vert_exp(p, nlon, nlat, u, pbar)
end subroutine py_load_vert_exp

subroutine py_love(n, ch, h, l, k)
    use wfl
    integer, intent(in) :: n
    character(len=*), intent(in) :: ch
    real *8, dimension(0:n), intent(out):: h, l, k
    select case(ch)
        case('wang')
            call love(n, h, l, k)
        case('farrell')
            call wfl_love_farrell(n, h, l, k)
        case default
            write(*, *) "Unkown"
            stop
    end select
end subroutine py_love

subroutine py_wfl_cs2cs(n, c, s, cout, sout, cha)
    use wfl
    integer, intent(in) ::  n
    real *8, dimension(0:n, 0:n), intent(in)  :: c, s
    real *8, dimension(0:n, 0:n), intent(out) :: cout, sout
    character(len=3), optional, intent(in) :: cha
    call wfl_cs2cs(n, c, s, cout, sout, cha)
end subroutine py_wfl_cs2cs

subroutine py_wfl_cs2vertical(n, c, s, u, depart_lon, alat, water_density, frame)
    use wfl
    integer, intent(in) :: n !, nlon, nlat
    real *8, dimension(0:n, 0:n), intent(in) :: c, s
    real *8, dimension(2*n, n), intent(out) :: u
    real *8, optional, intent(in) :: depart_lon, alat(n), water_density
    character(len=*), optional, intent(in) :: frame
    if(present(depart_lon) .and. depart_lon > 0) then
        call wfl_cs2vertical(n, c, s, u, depart_lon, alat, water_density, frame)
    else
        call wfl_cs2vertical(n, c, s, u, water_density=water_density, frame=frame)
    endif
end subroutine py_wfl_cs2vertical

subroutine py_wfl_var_reduction(n, x, x_correction, var_out, var_reduction, var_inversion)
   use wfl
   integer, intent(in) :: n
   real *8, dimension(n), intent(in) :: x, x_correction
   real *8, dimension(3), intent(out) :: var_out, var_reduction(2)
   real *8, optional, dimension(3), intent(out) :: var_inversion
   call wfl_var_reduction(n, x, x_correction, var_out, var_reduction, var_inversion)
end subroutine py_wfl_var_reduction

function py_wfl_depth2pressure(depth, lat) result (pressure)
    use wfl, only: wfl_depth2pressure
    real *8, intent(in) :: depth, lat
    real *8 :: pressure
    pressure = wfl_depth2pressure(depth, lat)
end function py_wfl_depth2pressure

function py_wfl_dm2lod(dm) result (dm2lod)
    use wfl, only: wfl_dm2lod
    real *8, intent(in) :: dm
    real *8 :: dm2lod
    dm2lod = wfl_dm2lod(dm)
end function py_wfl_dm2lod

subroutine py_wfl_fan_filter(nmax, mmax, r_lon, r_lat, w)
    use wfl, only: wfl_fan_filter
    integer, intent(in) :: nmax, mmax
    real *8, intent(in) :: r_lon, r_lat
    real *8, dimension(0:nmax, 0:nmax), intent(out) :: w
    call wfl_fan_filter(nmax, mmax, r_lon, r_lat, w)
end subroutine py_wfl_fan_filter

subroutine py_wfl_green(n, theta, n_max_love, green_value, cha, cha_norm)
    use wfl, only: wfl_green
    integer, intent(in) :: n, n_max_love
    real *8, intent(in) :: theta(n)
    real *8, intent(out) :: green_value(n)
    character(len=*), intent(in) :: cha
    character(len=*), optional, intent(in) :: cha_norm
    call wfl_green(n, theta, n_max_love, green_value, cha, cha_norm)
end subroutine py_wfl_green

subroutine py_wfl_green_func_theta(n, theta)
    use wfl, only: wfl_green_func_theta
    integer, intent(in)::n
    real *8, dimension(n), intent(out):: theta
    call wfl_green_func_theta(n, theta)
end subroutine py_wfl_green_func_theta

subroutine py_wfl_grace_destripes(n, c, s, i_segment, n_max, m_max, c_out, s_out, method, fit_order)
    use wfl, only: wfl_grace_destripes
    integer, intent(in) :: n, i_segment
    real *8, dimension(0:n, 0:n), intent(in) :: c, s
    integer, dimension(i_segment), intent(in) :: n_max, m_max
    real *8, dimension(0:n, 0:n), intent(out) :: c_out, s_out
    character(len=3), optional, intent(in) :: method
    integer, optional, dimension(i_segment), intent(in) :: fit_order
    call wfl_grace_destripes(n, c, s, i_segment, n_max, m_max, c_out, s_out, method, fit_order)
end subroutine py_wfl_grace_destripes

function py_wfl_in_situ_temp(s, theta, p0, pr0) result (t0)
    use wfl, only: wfl_in_situ_temp
    real *8, intent(in) :: s, theta, p0
    real *8, optional, intent(in) :: pr0
    real *8 :: t0
    t0 = wfl_in_situ_temp(s, theta, p0, pr0)
end function py_wfl_in_situ_temp

function py_wfl_potential_temp(s, t0, p0, pr0) result (theta)
    use wfl, only: wfl_potential_temp
    real *8, intent(in) :: s, t0, p0
    real *8, optional, intent(in) :: pr0
    real *8 :: theta
    theta = wfl_potential_temp(s, t0, p0, pr0)
end function py_wfl_potential_temp

subroutine py_wfl_j2j4j6(r, gm, omiga, inv_f, j, c)
    use wfl, only: wfl_j2j4j6
    real *8, intent(in) :: r, gm, omiga, inv_f
    real *8, intent(out) :: j(3), c(3)
    call wfl_j2j4j6(r, gm, omiga, inv_f, j, c)
end subroutine py_wfl_j2j4j6

subroutine py_wfl_lod2torque(n, dt, lod, torque_lod)
    use wfl, only: wfl_lod2torque
    integer, intent(in) ::  n
    real *8, intent(in) :: dt
    real *8, dimension(n), intent(in) :: lod
    real *8, dimension(n), intent(out) :: torque_lod
    call wfl_lod2torque(n, dt, lod, torque_lod)
end subroutine py_wfl_lod2torque

subroutine py_wfl_mountain_torque(lon, lat, sp, hm, te, ten, tm_dp, tmn_dp, tm_dh, tmn_dh, &
    tm3_dp, tm3n_dp, tm3_dh, tm3n_dh)
    use wfl, only: wfl_mountain_torque
    integer, intent(in) :: lon, lat
    real *8, dimension(lon, lat), intent(in) :: sp, hm !centered data
    complex(8), dimension(lon, lat), intent(out):: te, ten, tm_dp, tmn_dp, tm_dh, tmn_dh
    real *8, dimension(lon, lat), intent(out) :: tm3_dp, tm3n_dp, tm3_dh, tm3n_dh
    call wfl_mountain_torque(lon, lat, sp, hm, te, ten, tm_dp, tmn_dp, tm_dh, tmn_dh, &
        tm3_dp, tm3n_dp, tm3_dh, tm3n_dh)
end subroutine py_wfl_mountain_torque

subroutine py_wfl_prem(n, prem_out, prem_parameter)
    use wfl, only: wfl_prem
    integer, intent(in) :: n
    real *8, intent(out) :: prem_out(n, 2)
    character(len=*), optional, intent(in) :: prem_parameter
    call wfl_prem(n, prem_out, prem_parameter)
end subroutine py_wfl_prem

function py_wfl_pressure2depth(pressure, lat) result (depth)
    use wfl, only: wfl_pressure2depth
    real *8, intent(in) :: pressure, lat
    real *8 :: depth
    depth = wfl_pressure2depth(pressure, lat)
end function py_wfl_pressure2depth

subroutine py_wfl_rm_sp_from_sg(n, fre, fre_idx, spec_sp, spec_sg, spec_out_sg, m_order)
    use wfl, only: wfl_rm_sp_from_sg
    integer, intent(in) :: n, fre_idx(2)
    real *8, dimension(n), intent(in) :: fre, spec_sp, spec_sg
    real *8, dimension(n), intent(out) :: spec_out_sg
    integer, optional, intent(in) :: m_order
    if(present(m_order) .and. m_order .lt. 0) then
        call wfl_rm_sp_from_sg(n, fre, fre_idx, spec_sp, spec_sg, spec_out_sg, m_order)
    else
        call wfl_rm_sp_from_sg(n, fre, fre_idx, spec_sp, spec_sg, spec_out_sg)
    endif
end subroutine py_wfl_rm_sp_from_sg

subroutine py_wfl_search_resonance(n, fre, spec, n_window, q, resonance_parameter)
    use wfl, only: wfl_search_resonance
    integer, intent(in) :: n, n_window
    real *8, intent(in) :: spec(n), q, fre(n)
    real *8, intent(out) :: resonance_parameter(n)
    call wfl_search_resonance(n, fre, spec, n_window, q, resonance_parameter)
end subroutine py_wfl_search_resonance

subroutine py_wfl_sg_multistation( ns, alon, alat, n, sg, sg_out )
    use wfl, only: wfl_sg_multistation
    integer, intent(in) :: n, ns
    real *8, intent(in) :: alon(n), alat(n), sg(n, ns)
    complex(8), intent(out):: sg_out(n, 3)
    call wfl_sg_multistation( ns, alon, alat, n, sg, sg_out )
end subroutine py_wfl_sg_multistation

subroutine py_wfl_soiltemp(amp, period, phase, average_t, depth, year_begin, nday, soil_t, td, unit)
    use wfl, only: wfl_soiltemp
    real *8, intent(in) :: amp, period, phase, average_t, depth
    integer, intent(in) :: year_begin, nday
    real *8, intent(out) :: soil_t(nday, 3)
    real *8, optional, intent(in) :: td
    character(len=*), optional, intent(in) :: unit
    call wfl_soiltemp(amp, period, phase, average_t, depth, year_begin, nday, soil_t, td, unit)
end subroutine py_wfl_soiltemp

subroutine py_wfl_solar_heat_tide(n, fre, spec, m, index_fre, n_window, q, spec_out)
    use wfl, only: wfl_solar_heat_tide
    integer, intent(in) :: n, m, n_window
    real *8, intent(in) :: spec(n), q, fre(n)
    integer, intent(in) :: index_fre(m)
    real *8, intent(out) :: spec_out(n)
    call wfl_solar_heat_tide(n, fre, spec, m, index_fre, n_window, q, spec_out)
end subroutine py_wfl_solar_heat_tide

subroutine py_wfl_stress_torque(lon, lat, ets, nts, tf, tf3)
    use wfl, only: wfl_stress_torque
    integer, intent(in) :: lon, lat
    real *8, dimension(lon, lat), intent(in) :: ets, nts
    complex(8), dimension(lon, lat), intent(out):: tf
    real *8, dimension(lon, lat), intent(out) :: tf3
    call wfl_stress_torque(lon, lat, ets, nts, tf, tf3)
end subroutine py_wfl_stress_torque

subroutine py_wfl_thermalv(amp, period, phase, year_begin, nday, dh, pr, ltec, td, unit)
    use wfl, only: wfl_thermalv
    real *8, intent(in) :: amp, period, phase
    integer, intent(in) :: year_begin, nday
    real *8, intent(out) :: dh(nday, 3)
    real *8, optional, intent(in) :: pr, ltec, td
    character(len=*), optional, intent(in) :: unit
    call wfl_thermalv(amp, period, phase, year_begin, nday, dh, pr, ltec, td, unit)
end subroutine py_wfl_thermalv

subroutine py_wfl_torque(lon, lat, sp, hm, ets, nts, egws, ngws, te, ten, tm_dp, tmn_dp, tm_dh, tmn_dh, tf, tg, tm3_dp, tm3n_dp, tm3_dh, tm3n_dh, tf3, tg3)
    use wfl, only: wfl_torque
    integer, intent(in) :: lon, lat
    real *8, dimension(lon, lat+1), intent(in) :: sp, hm(lon, lat), ets, nts, egws, ngws
    complex(8), dimension(lon, lat), intent(out):: te, ten, tm_dp, tmn_dp, tm_dh, tmn_dh, tf, tg
    real *8, dimension(lon, lat), intent(out) :: tm3_dp, tm3n_dp, tm3_dh, tm3n_dh, tf3, tg3
    call wfl_torque(lon, lat, sp, hm, ets, nts, egws, ngws, te, ten, tm_dp, tmn_dp, tm_dh, tmn_dh, tf, tg, tm3_dp, tm3n_dp, tm3_dh, tm3n_dh, tf3, tg3)
end subroutine py_wfl_torque

subroutine py_wfl_torque2lod(n, dt, torque, dlod_torque)
    use wfl, only: wfl_torque2lod
    integer, intent(in) ::  n
    real *8 :: dt
    real *8, dimension(n), intent(in) :: torque
    real *8, dimension(n), intent(out) :: dlod_torque
    call wfl_torque2lod(n, dt, torque, dlod_torque)
end subroutine py_wfl_torque2lod

subroutine py_wfl_torque2lod_adams(n, dt, torque, lod0, lod, t, lod_poly, flag)
    use wfl, only: wfl_torque2lod_adams
    integer, intent(in) ::  n
    real *8, intent(in) :: dt
    real *8, dimension(n), intent(in) :: torque, lod0(3)
    real *8, dimension(n), intent(out) :: lod
    real *8, optional, dimension(n), intent(in) :: t
    real *8, optional, dimension(n), intent(out) :: lod_poly
    logical, optional, intent(in) :: flag
    if(present(flag) .and. flag) then
        call wfl_torque2lod_adams(n, dt, torque, lod0, lod, t, lod_poly)
    else
        call wfl_torque2lod_adams(n, dt, torque, lod0, lod)
    endif
end subroutine py_wfl_torque2lod_adams

subroutine py_wfl_vd_station_index(nlat, site_lon, site_lat, depart_lon, alat, idx_lon, idx_lat)
    use wfl, only: wfl_vd_station_index
    integer, intent(in) :: nlat
    real *8, intent(in) :: site_lon, site_lat
    real *8, intent(out) :: depart_lon
    real *8, dimension(nlat), intent(out) :: alat
    integer, intent(out) :: idx_lon, idx_lat
    call wfl_vd_station_index(nlat, site_lon, site_lat, depart_lon, alat, idx_lon, idx_lat)
end subroutine py_wfl_vd_station_index

subroutine py_wfl_wind2tao(windu, windv, stressu, stressv, air_density)
    use wfl, only: wfl_wind2tao
    real *8, intent(in) :: windu, windv
    real *8, intent(out) :: stressu, stressv
    real *8, optional, intent(in) :: air_density
    call wfl_wind2tao(windu, windv, stressu, stressv, air_density)
end subroutine py_wfl_wind2tao

function py_windv2tao(v) result (tao)
    use wfl, only: windv2tao
    real *8, intent(in) :: v
    real *8 tao
    tao = windv2tao(v)
end function py_windv2tao

subroutine py_x2prograd(a1, fai1, a2, fai2, a_pro, fai_pro, a_retro, fai_retro)
    use wfl, only: x2prograd
    real, intent(in):: a1, fai1, a2, fai2
    real, intent(out):: a_pro, a_retro, fai_pro, fai_retro
    call x2prograd(a1, fai1, a2, fai2, a_pro, fai_pro, a_retro, fai_retro)
end subroutine py_x2prograd

subroutine py_x2prograd_0(a1, fai1, a2, fai2, f1, f2, g1, g2)
    use wfl, only: x2prograd_0
    real, intent(in):: a1, fai1, a2, fai2
    real, intent(out):: f1, f2, g1, g2
    call x2prograd_0(a1, fai1, a2, fai2, f1, f2, g1, g2)
end subroutine py_x2prograd_0

function py_angle_dis(fai1,lambda1,fai2,lambda2,ch) result(angle_dis_result)
    use wfl, only: angle_dis
    real *8,intent(in) :: fai1,lambda1,fai2,lambda2
    character(len=*),optional,intent(in):: ch
    real *8 :: angle_dis_result
    if(present(ch) .and. ch .eq. 'distance') then
        angle_dis_result = angle_dis(fai1,lambda1,fai2,lambda2,ch)
    else
        angle_dis_result = angle_dis(fai1,lambda1,fai2,lambda2)
    endif
end function py_angle_dis


subroutine py_eaam_xyz_slp(lon,lat,slp,hm,fp,fpn,f3p,f3pn,pmn,pm,xypn,xyp,zpn,zp,j2n,j2,j3n,j3,is_sp)
    use wfl, only: eaam_xyz_slp
    implicit none
    integer,intent(in):: lon,lat
    real,intent(in):: slp(lon,lat),hm(lon,lat)
    complex*8,dimension(lon,lat),intent(out):: fpn,fp,xyp,xypn
    real*8,dimension(lon,lat),intent(out):: f3pn,f3p,zp,zpn,pm,pmn,j2,j3,j2n,j3n
    logical,optional,intent(in) :: is_sp
    call eaam_xyz_slp(lon,lat,slp,hm,fp,fpn,f3p,f3pn,pmn,pm,xypn,xyp,zpn,zp,j2n,j2,j3n,j3,is_sp)
end subroutine py_eaam_xyz_slp

subroutine py_eaam_xyz_wind(lon,lat,lay,ltrop,pl,u,v,fw,fwl,f3w,f3wl,flag, slp,hm,z)
    use wfl, only: eaam_xyz_wind
    integer,intent(in):: lon,lat,lay,ltrop
    real,intent(in):: pl(lay),u(lon,lat,lay),v(lon,lat,lay)
    complex*8,intent(out):: fw(lon,lat),fwl(lon,lat)
    real*8,intent(out):: f3w(lon,lat),f3wl(lon,lat)
    real,optional,intent(in):: slp(lon,lat),hm(lon,lat),z(lon,lat,lay)
    logical, optional, intent(in) :: flag
    if(present(flag) .and. flag) then
        call eaam_xyz_wind(lon,lat,lay,ltrop,pl,u,v,fw,fwl,f3w,f3wl,slp,hm,z)
    else
        call eaam_xyz_wind(lon,lat,lay,ltrop,pl,u,v,fw,fwl,f3w,f3wl)
    endif
end subroutine py_eaam_xyz_wind

subroutine py_eaam_xyz_wind_en(lon,lat,lay,ltrop,pl,u,v,fw,fwl,f3w,f3wl,flag,slp,hm,z)
    use wfl, only: eaam_xyz_wind_en
    integer,intent(in):: lon,lat,lay,ltrop
    real,intent(in):: pl(lay),u(lon,lat,lay),v(lon,lat,lay)
    complex*8,intent(out):: fw(lon,lat,3),fwl(lon,lat,3)
    real*8,intent(out):: f3w(lon,lat),f3wl(lon,lat)
    real,optional,intent(in):: slp(lon,lat),hm(lon,lat),z(lon,lat,lay)
    logical, optional, intent(in) :: flag
    if(present(flag) .and. flag) then
        call eaam_xyz_wind_en(lon,lat,lay,ltrop,pl,u,v,fw,fwl,f3w,f3wl,slp,hm,z)
    else
        call eaam_xyz_wind_en(lon,lat,lay,ltrop,pl,u,v,fw,fwl,f3w,f3wl)
    endif
end subroutine py_eaam_xyz_wind_en

! subroutine py_eaam_xyz_wind_en_wahr(lon,lat,lay,ltrop,pl,u0,v0,fw,fwl,f3w,f3wl,flag,slp,hm,z)
!     use wfl, only: eaam_xyz_wind_en_wahr
!     integer,intent(in):: lon,lat,lay,ltrop
!     real,intent(in):: pl(lay),u0(lon,lat,lay),v0(lon,lat,lay)
!     complex*8,dimension(lon,lat,3),intent(out):: fw,fwl
!     real*8,dimension(lon,lat),intent(out):: f3w,f3wl
!     real,optional,intent(in):: slp(lon,lat),hm(lon,lat),z(lon,lat,lay)
!     logical, optional, intent(in) :: flag
!      if(present(flag) .and. flag) then
!         call eaam_xyz_wind_en_wahr(lon,lat,lay,ltrop,pl,u0,v0,fw,fwl,f3w,f3wl,flag,slp,hm,z)
!     else
!         call eaam_xyz_wind_en_wahr(lon,lat,lay,ltrop,pl,u0,v0,fw,fwl,f3w,f3wl)
!     endif
! end subroutine py_eaam_xyz_wind_en_wahr

function py_eam_barnes2eubanks(eam,cha) result(eam2)
    use wfl, only: eam_barnes2eubanks
    real*8,intent(in):: eam
    character(len=*):: cha
    real*8 eam2
    eam2 = eam_barnes2eubanks(eam, cha)
end function py_eam_barnes2eubanks

function py_eam2masms(xin,ch) result(xyz2masms_result)
    use wfl, only: eam2masms
    real*8,intent(in):: xin
    character(len=1),intent(in):: ch
    real*8 xyz2masms_result
    xyz2masms_result = eam2masms(xin,ch)
end function py_eam2masms

subroutine py_ewam_xyz_j2j3(lon,lat,pre,evp,runoff,hm,fp,fpn,f3p,f3pn,pmn,pm,xypn,xyp,zpn,zp,j2n,j2,j3n,j3,grid_status,method)
    use wfl, only: ewam_xyz_j2j3
    integer,intent(in):: lon,lat
    real*4,dimension(lon,lat+1),intent(inout)::pre,evp,runoff
    real*4,intent(in)::hm(lon,lat)
    complex*8,dimension(lon,lat),intent(out)::fpn,fp,xyp,xypn
    real*8,dimension(lon,lat),intent(out)::f3pn,f3p,zp,zpn,pm,pmn,j2,j3,j2n,j3n
    character(len=*) :: grid_status
    character(len=*),optional:: method
    call ewam_xyz_j2j3(lon,lat,pre,evp,runoff,hm,fp,fpn,f3p,f3pn,pmn,pm,xypn,xyp,zpn,zp,j2n,j2,j3n,j3,grid_status,method)
end subroutine py_ewam_xyz_j2j3

subroutine py_gauss_grid(gauss_lon,gauss_lat,gauss_lat_edges)
    use wfl, only: gauss_grid
    real,intent(out):: gauss_lon(192),gauss_lat(94),gauss_lat_edges(95)
    call gauss_grid(gauss_lon,gauss_lat,gauss_lat_edges)
end subroutine py_gauss_grid

subroutine py_geocenter(zlat_s,zlat_n,zlon_w,zlon_e,rho,r1,r2,x,y,z)
    use wfl, only: geocenter
    real*8,intent(in):: zlat_s,zlat_n,zlon_w,zlon_e,rho,r1,r2
    real*8,intent(out):: x,y,z
    call geocenter(zlat_s,zlat_n,zlon_w,zlon_e,rho,r1,r2,x,y,z)
end subroutine py_geocenter

function py_masms2eam(xin,ch) result(masms2xyz_result)
    use wfl, only: masms2eam
    real*8,intent(in):: xin
    character(len=1),intent(in):: ch
    real*8 masms2xyz_result
    masms2xyz_result = masms2eam(xin,ch)
end function py_masms2eam

subroutine py_oam_int(rho,r1,r2,zlat_s,zlat_n,zlon_w,zlon_e,u,v,xoamp,yoamp,zoamp,xoamc,yoamc,zoamc)
    use wfl, only: oam_int
    real*8,intent(in):: rho,r1,r2,zlat_s,zlat_n,zlon_w,zlon_e,u,v
    real*8,intent(out)::xoamp,yoamp,zoamp,xoamc,yoamc,zoamc
    call oam_int(rho,r1,r2,zlat_s,zlat_n,zlon_w,zlon_e,u,v,xoamp,yoamp,zoamp,xoamc,yoamc,zoamc)
end subroutine py_oam_int

 subroutine py_oam_simple(rho,r,zlat,zlon,dlat,dlon,dr,u,v,xoamp,yoamp,zoamp,xoamc,yoamc,zoamc)
    use wfl, only: oam_simple
    real*8,intent(in):: rho,r,zlat,zlon,dlat,dlon,dr,u,v
    real*8,intent(out)::xoamp,yoamp,zoamp,xoamc,yoamc,zoamc
    call oam_simple(rho,r,zlat,zlon,dlat,dlon,dr,u,v,xoamp,yoamp,zoamp,xoamc,yoamc,zoamc)
end subroutine py_oam_simple

subroutine py_ocnden(s,c68,depth,den,lat)
    use wfl, only: ocnden
    real*8,intent(in)::s,c68,depth
    real*8,optional,intent(in) :: lat
    real*8,intent(out)::den
    call ocnden(s,c68,depth,den,lat)
end subroutine py_ocnden

subroutine py_pm2xy(n,dt,pm,chi_xy,cha)
    use wfl, only: pm2xy
    integer,intent(in)::n
    real,intent(in)::dt
    complex,intent(in):: pm(n)
    complex,intent(out)::chi_xy(n)
    character*3 cha
    call pm2xy(n,dt,pm,chi_xy,cha)
end subroutine py_pm2xy

subroutine py_prograd2x_0(f1,f2,g1,g2,a1,fai1,a2,fai2)
    use wfl, only: prograd2x_0
    real,intent(in):: f1,f2,g1,g2
    real,intent(out):: a1,fai1,a2,fai2
    call prograd2x_0(f1,f2,g1,g2,a1,fai1,a2,fai2)
end subroutine py_prograd2x_0

subroutine py_soda_xyz_st(x360,y128,yedges,z20,zedges)
    use wfl, only: soda_xyz_st
    real*8,intent(out):: z20(20),y128(128),x360(360),yedges(129),zedges(21)
    call soda_xyz_st(x360,y128,yedges,z20,zedges)
end subroutine py_soda_xyz_st

subroutine py_soda_xyz_uvh(uvhx,uvhy,uvhyedges)
    use wfl, only: soda_xyz_uvh
    real*8,intent(out):: uvhx(360),uvhy(128),uvhyedges(129)
    call soda_xyz_uvh(uvhx,uvhy,uvhyedges)
end subroutine py_soda_xyz_uvh

subroutine py_sphere_len(b1,l1,b2,l2,slen)
    use wfl, only: sphere_len
    real *8,intent(in):: b1,l1,b2,l2
    real *8,intent(out):: slen
    call sphere_len(b1,l1,b2,l2,slen)
end subroutine py_sphere_len

subroutine py_sphere_len_inv(b,l,slen,delta_b,delta_l)
    use wfl, only: sphere_len_inv
    real *8,intent(in):: b,l,slen
    real *8,intent(out):: delta_b,delta_l
    call sphere_len_inv(b,l,slen,delta_b,delta_l)
end subroutine py_sphere_len_inv

function py_volume(dlon,lat1,lat2,z1,z2) result(volume_result)
    use wfl, only: volume
    real*8,intent(in):: dlon,lat1,lat2,z1,z2
    real*8 volume_result
    volume_result = volume(dlon,lat1,lat2,z1,z2)
end function py_volume

function py_wfl_atg(s,t,p0) result(atg)
    use wfl, only: wfl_atg
    real*8,intent(in) :: s,t,p0
    real*8 :: atg
    atg = wfl_atg(s,t,p0)
end function py_wfl_atg






