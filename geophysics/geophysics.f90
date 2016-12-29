module geophysics
    implicit none
    
    contains
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

        function py_area(dlon,lat1,lat2) result(area_result)
            use wfl
            real*8, intent(in):: dlon,lat1,lat2
            real*8 :: area_result
            area_result = area(dlon, lat1, lat2)
        end function py_area

        function py_cs2kai(cs,cha)  result(kai)
            use wfl
            implicit none
            real*8,intent(in) :: cs
            character(len=1),intent(in):: cha
            real*8 kai
            kai = cs2kai(cs, cha)
        end function py_cs2kai

        function py_kai2cs(kai,cha)  result(cs)
            use wfl
            implicit none
            real*8,intent(in) :: kai
            character(len=1),intent(in):: cha
            real*8 cs
            cs = kai2cs(kai, cha)
        end function py_kai2cs

        subroutine py_degree_var(n,c,s,degvar)
            use wfl
            integer, intent(in):: n
            real *8, dimension(0:n,0:n),intent(in):: c,s
            real *8, dimension(0:n),intent(out):: degvar
            call degree_var(n,c,s,degvar)
        end subroutine py_degree_var

        subroutine j2j3(colat,dlat,dlon,hi,rho,dj2,dj3)
            use wfl
            real*8,intent(in):: colat,dlat,dlon,hi,rho
            real*8,intent(out):: dj2,dj3
            call dj2j3(colat,dlat,dlon,hi,rho,dj2,dj3)
        end subroutine j2j3

        subroutine py_eaam_xyz_j2j3(lon,lat,lay1,ltrop,slp,pl,u,v,z,hm,&
                  fp,fpn,f3p,f3pn,fw,fwl,f3w,f3wl,pmn,pm,xypn,xyp,zpn,zp,&
                  j2n,j2,j3n,j3,method)
            use wfl
            integer, intent(in):: lon,lat,lay1,ltrop
            real *4, dimension(lon,lat+1,lay1), intent(inout)::z,u,v,slp(lon,lat+1)
            real *4, intent(in):: hm(lon,lat), pl(lay1)
            complex *8, dimension(lon,lat),intent(out):: fpn,fp,fw,fwl,xyp,xypn
            real *8 , dimension(lon,lat),intent(out):: f3pn,f3p,f3w,f3wl,zp,zpn,pm,pmn,&
                      j2,j3,j2n,j3n
            character(len=*), intent(in) :: method
            call eaam_xyz_j2j3(lon,lat,lay1,ltrop,slp,pl,u,v,z,hm,&
                  fp,fpn,f3p,f3pn,fw,fwl,f3w,f3wl,pmn,pm,xypn,xyp,zpn,zp,&
                  j2n,j2,j3n,j3,method)
        end subroutine py_eaam_xyz_j2j3

        subroutine py_load_vert_def(nlon, nlat, sigma, u, depart_lon, zlat)
            !purpose: get the radial elastic deformation of the earth's surface estimated from a given surface load globally.
            use wfl
            implicit none
            integer,intent(in):: nlon,nlat
            real *8,intent(in):: sigma(nlon,nlat) !unit: kg/m/m
            real *8,intent(out):: u(nlon,nlat)
            real *8,optional,intent(in):: depart_lon,zlat(nlat)  !must be presented at the same time
            call load_vert_def(nlon, nlat, sigma, u, depart_lon, zlat)
        end subroutine py_load_vert_def

        subroutine py_load_vert_exp(p, nlon, nlat, u, pbar)
            use wfl
            implicit none
            integer,intent(in):: nlon,nlat
            real *8,dimension(nlon,nlat),intent(in):: p
            real *8,dimension(nlon,nlat),intent(out):: u, pbar
            call load_vert_exp(p, nlon, nlat, u, pbar)
        end subroutine py_load_vert_exp

        subroutine py_love(n, ch, h, l, k)
            use wfl
            integer, intent(in) :: n
            character(len=*), intent(in) :: ch
            real*8, dimension(0:n),intent(out):: h,l,k
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

        subroutine py_wfl_cs2cs(n,c,s,cout,sout,cha)
            use wfl
            integer, intent(in) ::  n
            real *8, dimension(0:n,0:n),intent(in)  :: c,s
            real *8, dimension(0:n,0:n),intent(out) :: cout,sout
            character(len=3),optional,intent(in) :: cha
            call wfl_cs2cs(n,c,s,cout,sout,cha)
        end subroutine py_wfl_cs2cs

        subroutine  py_wfl_cs2vertical(n,c,s,u,depart_lon,alat,water_density,frame)
            use wfl
            integer, intent(in) :: n !,nlon,nlat
            real *8, dimension(0:n,0:n),intent(in) :: c, s
            real *8, dimension(2*n,n),intent(out) :: u
            real *8, optional,intent(in) :: depart_lon,alat(n),water_density
            character(len=*),optional,intent(in) :: frame
            call wfl_cs2vertical(n,c,s,u,depart_lon,alat,water_density,frame)
        end subroutine py_wfl_cs2vertical

        subroutine py_wfl_var_reduction(n, x, x_correction, var_out, var_reduction, var_inversion)
           use wfl
           integer, intent(in) :: n
           real *8, dimension(n),intent(in) :: x, x_correction
           real *8, dimension(3),intent(out) :: var_out, var_reduction(2)
           real *8, optional,dimension(3),intent(out) :: var_inversion
           call wfl_var_reduction(n, x, x_correction, var_out, var_reduction, var_inversion)
        end subroutine py_wfl_var_reduction

        
end module geophysics