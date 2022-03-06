
MODULE Tridiagonal

                implicit none

        public:: SLAE3diag, matrix3diag_mult

        private

contains

        !A three-diagonal matrix solver for systems of linear algebraic equations
        !Solves simultaneously m systems with 3-diagonal matrices: A^{(m)}x^{(m)} = b^{(m)}
        !The calling syntax is:
        !call SLAE3diag(N,     !intent(in) size of the system (number of variables/equations)
        !               M,     !intent(in) number of systems to be solved simultaneuosly
        !               Ldiag, !intent(in), dimension(n,m) set of lower diagonals of m matrices: a_{i,i-1}^{(m)}; note that Ldiag(1,:) are ignored.
        !                diag, !intent(in), dimension(n,m) set of main diagonals of m matrices: a_{i,i}^{(m)}.
        !               Udiag, !intent(in), dimension(n,m) set of upper diagonals of m matrices: a_{i,i+1}^{(m)}; note that Udiag(n,:) are ignored.
        !                 RHS, !intent(in), dimension(n,m) set of right-hand sides of m systems: b_{i}^{(m)}.
        !                   x, !intent(out), dimension(n,m) set of solutions returned in this parameter: x_{i}^{(m)}.
        subroutine SLAE3diag(n, m, Ldiag, diag, Udiag, rhs, x)
                implicit none
                integer, intent(in):: n, m                                        !order of the systems (n x n matrix) and number of simultaneous systems to solve
                double precision, dimension(n,m), intent(in):: Ldiag, diag, Udiag !lower diagonal, main diagonal, upper diagonal
                                                                                  !note that ?diag(n,?) is an element of the line n; so Ldiag(1,:) and Udiag(n,:) are ignored
                double precision, dimension(n,m), intent(in):: rhs                !right-hand side of the system(s)
                double precision, dimension(n,m), intent(out):: x                 !the solutions to be returned
                double precision, dimension(n,m):: a,b                            !sweeping coefficients
                double precision :: den                                           !denominator
                double precision :: den_r                                         !1/denominator
!               double precision, dimension(m):: den                              !denominator
                integer:: i,j

                do j=1,m               
!                  write(*,*) j
                   a(1,j) = -Udiag(1,j)/diag(1,j)                                    !express x_1 via x_2 (x1=a1*x2+b1) using the first row AKA equation 1
                   b(1,j) = rhs(1,j)/diag(1,j)                                       !express x_1 via x_2 (x1=a1*x2+b1) using the first row AKA equation 1                                   

                   do i=2,n-1                                                        !express x_i via x_{i+1} in the similar way
                      den = diag(i,j) + a(i-1,j)*Ldiag(i,j)
                      if (den .GE. 0.0d0) then 
                         den_r = 1.0D0/(den + 0.000001D0)
                      else
                         den_r = 1.0d0/(den - 0.000001D0)
                      endif
                      a(i,:) = -Udiag(i,:)*den_r
                      b(i,:) = ( rhs(i,j) - b(i-1,j)*Ldiag(i,j) ) * den_r
                   enddo
                   den =  Ldiag(n,j)*a(n-1,j) + diag(n,j) 
                   if (den .GE. 0.0d0) then 
                         den_r = 1.0D0/(den + 0.000001D0)
                   else
                         den_r = 1.0d0/(den - 0.000001D0)
                   endif
                   x(n,j) = ( rhs(n,j) - Ldiag(n,j)*b(n-1,j) ) * den_r  !exclude x_{n-1} and solve the last equation to get x_n
                   do i=n-1,1,-1                                        !use the connections to get all x_i one by one
                      x(i,j) = a(i,j)*x(i+1,j) + b(i,j)
                   enddo

                enddo


        end subroutine SLAE3diag

        !matrix multiplicator for 3-diagonal matrices; puts the results into the y matrix. If RHS is provided, returns the error in the err (that must be also provided then)
        subroutine matrix3diag_mult(n, m, Ldiag, diag, Udiag, x, y, rhs, err)
                implicit none
                integer, intent(in):: n, m !order of the system (n x n matrix) and number of simultaneous systems to solve
                double precision, dimension(n,m), intent(in):: Ldiag, diag, Udiag
                double precision, dimension(n,m), intent(in), optional:: rhs
                double precision, dimension(m), intent(out), optional:: err
                double precision, dimension(n,m), intent(in):: x
                double precision, dimension(n,m), intent(out):: y
                double precision, dimension(n,m):: a,b
                double precision, dimension(m):: den
                integer:: i, sys
                y = 0.
                do sys = 1,m
                        do i=1,n
                                    y(i,sys) = y(i,sys) + diag(i,sys)*x(i,sys)
                            if(i>1) y(i,sys) = y(i,sys) + Ldiag(i,sys)*x(i-1,sys)
                            if(i<n) y(i,sys) = y(i,sys) + Udiag(i,sys)*x(i+1,sys)
                        enddo
                enddo
                if(present(rhs)) err = sqrt(sum(((y-rhs)/(rhs+0.0000001))**2,dim=1))
        end subroutine matrix3diag_mult
END MODULE

