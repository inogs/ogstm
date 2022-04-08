subroutine OASIM_CALLER()

use myalloc
use OPT_mem
use oasim, only: calc_unit, oasim_lib

IMPLICIT NONE

! local

INTEGER ji, jj, jl

integer year,julianday
INTEGER points(1)
logical ERR

double precision :: sec_b, sec_e
double precision, dimension(1,33) :: local_taua, local_asymp, local_ssalb, Edout, Esout
double precision, dimension(1) :: l_sp, l_msl, l_w10, l_tco3, l_t2m, l_d2m, l_tcc, l_tclw, l_cdrem


points(1)=1
year=2019
julianday=30
sec_b = 43200.0 ! 12:00
sec_e = 46800.0   ! 13:00

do ji = 1,1
do jl = 1,33
local_taua(ji, jl )  = taua(jl,1,ji)
local_asymp(ji, jl ) = asymp(jl,1,ji)
local_ssalb(ji, jl ) = ssalb(jl,1,ji)
enddo
enddo


do ji=1,1
l_sp(ji) = sp(1,ji)
l_msl(ji) = msl(1,ji)
l_w10(ji) = w10(1,ji)
l_tco3(ji)= tco3(1,ji)
l_t2m(ji) = t2m(1,ji)
l_d2m(ji) = d2m(1,ji)
l_tcc(ji) = tcc(1,ji)
l_tclw(ji) = tclw(1,ji)
l_cdrem(ji) = cdrem(1,ji)

enddo




write(*,*) 'sp = ', l_sp
write(*,*) 'msl = ', l_msl
write(*,*) 'w10 = ', l_w10
write(*,*) 'tco3 = ', l_tco3
write(*,*) 't2m = ', l_t2m
write(*,*) 'd2m = ',  l_d2m
write(*,*) 'tcc = ', l_tcc
write(*,*) 'tclw = ', l_cdrem
write(*,*) 'cdrem = ', l_cdrem





call calc%monrad(points,year,julianday, sec_b, sec_e, &
                   l_sp, l_msl, l_w10, l_tco3, l_t2m, l_d2m, &
                   l_tcc, l_tclw, l_cdrem, &
                   local_taua, local_asymp, local_ssalb, &
                   Edout, Esout, ERR)


write(*,*) 'ERR = ', ERR

do jl=1,33
write(*,*) Edout(1,jl), Esout(1,jl)
enddo




end subroutine OASIM_CALLER
