module harmonic
    real *8, parameter :: pi = 3.1415926535897932384626d0

    integer polys, n, num_par
    real *8 timebegin, std
    real *8, allocatable, dimension(:) :: t, x, periods, residuals, out, out_fit
    real *8, allocatable, dimension(:, :) :: pdef, design_matrix, cov, covar

    contains
        subroutine gen_design_matrix(poly)
            integer, intent(in) :: poly
            real *8, allocatable, dimension(:) :: t_temp
            integer num, i, j, m, p, n
            integer period_length
            n = size(t)
            periods_length = size(periods)
            num_par = 1 + poly
            periods_flag = .False.

            if(periods_length > 0) then
                num_par = num_par + size(periods) * 2
                periods_flag = .True.
            endif

            if(size(pdef, 1) > 0) then
                num_par = num_par + size(pdef, 2)
            endif

            allocate(t_temp(n))
            t_temp(:) = t(:) - timebegin
            allocate(design_matrix(n, num_par))
            do i = 1, n
                num = 1
                design_matrix(i, num) = 1
                do j = 1, poly
                    design_matrix(i, num+j) = t_temp(i) ** j
                enddo
                num = num + poly

                do m = 1, periods_length
                    design_matrix(i, num+m*2-1) = sin(2.0 * pi * t_temp(i) / periods(m))
                    design_matrix(i, num+m*2) = cos(2.0 * pi * t_temp(i) / periods(m))
                enddo
                num = num + 2*periods_length
            enddo

            deallocate(t_temp)
            if(size(pdef, 1) > 0) then
                design_matrix(:, num+1:num_par) = pdef(:, :)
            endif

        end subroutine gen_design_matrix


        subroutine fit
            implicit none
            integer i, j, flag
            real *8 chisq
            real *8, allocatable, dimension(:, :) :: A, AP, Q
            real *8, allocatable, dimension(:) :: Y, res, std_temp

            ! check input
            if(.not.(allocated(t) .and. allocated(x))) then
                print *, "t or x are not initialized!"
                stop
            else if(size(t) .ne. size(x)) then
                print *, "t and x are not the same length!"
                stop
            endif

            if(.not. allocated(design_matrix)) then
                print *, "Design matrix is not generated!"
                stop
            endif
            n = size(t)

            ! if(polys >= 0) then
            !     call gen_design_matrix(polys)
            ! endif

            if(.not. allocated(cov)) then
                allocate(cov(n, n))
                cov = 0
                forall(i=1:n) cov(i, i) = 1
            endif

            ! do condition adjustment
            ! y = Ax --> x = (A^{T}PA)^{-1}(A^{T}Py), here AP -> A^{T}P, A -> A^{T}PA, Y -> A^{T}Py
            allocate(A(num_par, num_par), Y(num_par), AP(num_par, n))
            AP = matmul(transpose(design_matrix), cov)
            A = matmul(AP, design_matrix)
            Y = matmul(AP, x)
            allocate(Q(num_par, num_par))
            call inverse(A, Q, num_par)
            out = matmul(Q, Y)
            out_fit = matmul(design_matrix, out)
            ! allocate(residuals(n))
            residuals(:) = out_fit(:) - x(:)

            allocate(std_temp(n))
            do i = 1, n
                std_temp(i) = dot_product(residuals, cov(:, i))
            enddo
            std = dot_product(std_temp, residuals) / (n - num_par)
            covar = Q(:, :) * std
            deallocate(A, Y, AP, Q, design_matrix, std_temp)

        end subroutine fit


end module harmonic

subroutine inverse(a,c,n)
    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! Alex G. December 2009
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - dimension
    ! output ...
    ! c(n,n) - inverse matrix of A
    ! comments ...
    ! the original matrix a(n,n) will be destroyed
    ! during the calculation
    !===========================================================
    implicit none
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do
end subroutine inverse

!c************************************************************************
!c***                      subroutine qr(a,m,n,q,l)                   *
!c************************************************************************
!c  this subroutine is to calculate the q-r matrix using the householder *
!c  method                                                                *
!c  *  a(m,n)=the matrix to be decomposed to q-r matrix(m>=n)   when input*   *
!c         when output is the up triangle matrix r                    *
!c     q(m*m)=retun the q matrix                                          *
!c     l=0 or not o,if =0 then fail,else right.                       *
!c  this subroutine copied from qinghua program                          *
!c-----------------------------------------------------------------------*
subroutine qr(a,m,n,q,l)
    implicit none
    integer,intent(in):: m,n
    real*8,intent(inout):: a(m,n)
    real*8,intent(out):: q(m,m)
    integer,intent(out):: l

    integer i,j,k,nn
    real*8 alpha,t,u

    if (m.lt.n) then
      l=0
      write(*,40)
      return
    end if
40  format(1x,'  fail')
    do 10 i=1,m
    do 10 j=1,m
      q(i,j)=0.0
      if (i.eq.j) q(i,j)=1.0
10  continue
    nn=n
    if (m.eq.n) nn=m-1
    do 200 k=1,nn
      u=0.0
      do 20 i=k,m
        if (abs(a(i,k)).gt.u) u=abs(a(i,k))
20    continue
      alpha=0.0
      do 30 i=k,m
        t=a(i,k)/u
        alpha=alpha+t*t
30    continue
      if (a(k,k).gt.0.0) u=-u
      alpha=u*sqrt(alpha)
      if (abs(alpha)+1.0.eq.1.0) then
        l=0
        write(*,40)
        return
      end if
      u=sqrt(2.0*alpha*(alpha-a(k,k)))
      if (u+1.0.ne.1.0) then
        a(k,k)=(a(k,k)-alpha)/u
        do 50 i=k+1,m
50      a(i,k)=a(i,k)/u
        do 80 j=1,m
          t=0.0
          do 60 l=k,m
60        t=t+a(l,k)*q(l,j)
          do 70 i=k,m
70        q(i,j)=q(i,j)-2.0*t*a(i,k)
80      continue
        do 110 j=k+1,n
          t=0.0
          do 90 l=k,m
90        t=t+a(l,k)*a(l,j)
          do 100 i=k,m
100       a(i,j)=a(i,j)-2.0*t*a(i,k)
110     continue
        a(k,k)=alpha
        do 120 i=k+1,m
120     a(i,k)=0.0
      end if
200 continue
    l=1
    do 210 i=1,m-1
    do 210 j=i+1,m
      t=q(i,j)
      q(i,j)=q(j,i)
      q(j,i)=t
210 continue
    return
end