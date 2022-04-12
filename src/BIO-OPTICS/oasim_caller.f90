subroutine OASIM_CALLER(datestring)

use myalloc
use OPT_mem
use TIME_MANAGER
use mpi
use oasim, only: calc_unit, oasim_lib

IMPLICIT NONE

character(LEN=17), INTENT(IN) ::  datestring

! local

INTEGER ji, jj, jl
integer year, month, day, julianday, counter, seconds
INTEGER points(jpj*jpi)
logical ERR
integer INTERVAL_TO_CALL_OASIM ! should be defined in namelist


double precision :: sec, sec_b, sec_e, s_interval
double precision, dimension(jpi*jpj,33) :: local_taua, local_asymp, local_ssalb, Edout, Esout
double precision, dimension(jpj*jpi) :: l_sp, l_msl, l_w10, l_tco3, l_t2m, l_d2m, l_tcc, l_tclw, l_cdrem




INTERVAL_TO_CALL_OASIM=3600
call read_date_string(datestring, year, month, day, sec)
seconds = int(sec)

if (mod(seconds, INTERVAL_TO_CALL_OASIM).gt.0 ) RETURN

trcoptparttime = MPI_WTIME() ! cronometer-start


CALL forcings_KEXT(datestring)
CALL forcings_atm_clim(datestring)
CALL forcings_atm_aero(datestring)

counter=1
do ji=1,jpi
do jj=1,jpj
points(counter)=counter
counter=counter+1
enddo
enddo


call tau2julianday(datestringToTAU(datestring), deltaT, julianday)

s_interval=real(INTERVAL_TO_CALL_OASIM,8)

sec_b = sec - s_interval/2.
sec_e = sec + s_interval/2.



call flatten_2d(sp  ,l_sp)
call flatten_2d(msl ,l_msl)
call flatten_2d(w10 ,l_w10)
call flatten_2d(tco3,l_tco3)
call flatten_2d(t2m, l_t2m)
call flatten_2d(d2m, l_d2m)
call flatten_2d(tcc, l_tcc)
call flatten_2d(tclw, l_tclw)
call flatten_2d(cdrem, l_cdrem)

call flatten_33(taua,local_taua)
call flatten_33(asymp,local_asymp)
call flatten_33(ssalb,local_ssalb)




call calc%monrad(points,year,julianday, sec_b, sec_e, &
                   l_sp, l_msl, l_w10, l_tco3, l_t2m, l_d2m, &
                   l_tcc, l_tclw, l_cdrem, &
                   local_taua, local_asymp, local_ssalb, &
                   Edout, Esout, ERR)


call unflatten_33(Edout, Ed_0m)
call unflatten_33(Esout, Es_0m)


call trc3streams(datestring) ! 3-stream radiative model


trcoptparttime = MPI_WTIME() - trcoptparttime ! cronometer-stop
trcopttottime = trcopttottime + trcoptparttime


end subroutine OASIM_CALLER
