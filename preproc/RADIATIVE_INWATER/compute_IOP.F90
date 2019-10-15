program compute

USE myalloc, ONLY: jpi,jpj,jpk,glamt,gphit
USE OPT_mem
USE TIME_MANAGER

IMPLICIT NONE

integer            :: uni1,uni2, uni3, uni4
CHARACTER(len=128) :: INPUT_OASIM_FILE, INPUT_PFT_FILE, INPUT_CDOM_FILE, INPUT_NAP_FILE
CHARACTER(len=128) :: OUTPUT_FILE
CHARACTER(len=32)  :: time_string
character(LEN=17)  ::  datestring
integer :: i,k,chl,time_2hr,nc, wl
INTEGER :: year, month, day, ihr
INTEGER :: day_of_year
INTEGER :: it_actual
INTEGER :: MODE ! 0-exact, 1-approx
INTEGER :: bottom
CHARACTER(LEN=20) :: V_POSITION
double precision :: dT
double precision :: solz(1,1), rmud(1,1)
double precision :: Ed_OASIM(nlt),Es_OASIM(nlt)
double precision :: sec
double precision, allocatable :: depthz(:)
double precision, allocatable :: Edz(:,:),Esz(:,:),Euz(:,:),PARz(:,:)
double precision, allocatable :: CHLz(:,:),CDOMz(:,:),NAPz(:,:)
double precision :: Eu_0m(nlt)

OUTPUT_FILE='pippo.nc'

ifst = .TRUE. ! to initialize local vector

dT =1800.0D0

jpi  = 1
jpj  = 1

uni1 = 20
uni2 = 21

call getarg(1, INPUT_OASIM_FILE)
open(uni1, file= INPUT_OASIM_FILE, status="old",action="read")

do i=1,nlt
       read(uni1,*) Ed_OASIM(i), Es_OASIM(i)
       write(*,*)  Ed_OASIM(i), Es_OASIM(i)
enddo

close (uni1)


call getarg(2, INPUT_PFT_FILE)

write (*,*) " INPUT_PFT_FILE= ", INPUT_PFT_FILE

open(uni2, file= INPUT_PFT_FILE, status="old",action="read")

read(uni2,*) datestring

write(*,*) 'datestring ', datestring

allocate(glamt(jpj,jpi),gphit(jpj,jpi))

read(uni2,*) gphit(jpj,jpi)

write(*,*) 'latitude ', gphit(jpj,jpi)

read(uni2,*) jpk

write(*,*) 'Number of vertical levels ', jpk

allocate(e3t(jpk,1,1))
allocate(depthz(jpk))
allocate(CHLz(jpk,nchl),CDOMz(jpk, nlt),NAPz(jpk, nlt))

do i =1,jpk
       read(uni2,*) depthz(i),CHLz(i,1), CHLz(i,2), CHLz(i,3), CHLz(i,4)
       write(*,*)  depthz(i),CHLz(i,1), CHLz(i,2), CHLz(i,3), CHLz(i,4)
enddo

if ( depthz(1) > 0.0D0) then

     e3t(1,1,1) = depthz(1) ! First data at zero depth (OASIM Ed0-)

     do i=2,jpk
        e3t(i,1,1) = depthz(i) - depthz(i-1) 
     enddo    
else

     STOP ' BGC Argo data at 0 depth! '

endif

close (uni2)

call getarg(3, INPUT_CDOM_FILE)
write (*,*) " INPUT_CDOM_FILE= ", INPUT_CDOM_FILE
open(uni3, file= INPUT_CDOM_FILE, status="old",action="read")
do i =1,jpk
       read(uni3,*) depthz(i), (CDOMz(i,wl),wl=1,nlt)
       write(*,*)   depthz(i), (CDOMz(i,wl),wl=1,nlt)
enddo
close (uni3)


call getarg(4, INPUT_NAP_FILE)
write (*,*) " INPUT_NAP_FILE= ", INPUT_NAP_FILE
open(uni4, file= INPUT_NAP_FILE, status="old",action="read")
do i =1,jpk
       read(uni4,*) depthz(i), (NAPz(i,wl),wl=1,nlt)
       write(*,*)   depthz(i), (NAPz(i,wl),wl=1,nlt)
enddo
close (uni4)


call getarg(5, OUTPUT_FILE)

write (*,*) " OUTPUT_FILE= ", OUTPUT_FILE

call OPEN_AVE_edeseu(TRIM(OUTPUT_FILE), jpk,  'Edz', 'Esz', 'Euz', nc)

call  myalloc_OPT()

do i=1,nlt
   write(*,*) lam(i), aw(i), bw(i)
enddo


allocate(Edz(jpk,nlt),Esz(jpk,nlt),Euz(jpk,nlt),PARz(jpk,nchl+1))


call read_date_string(datestring, year, month, day, sec)

call tau2julianday(datestringToTAU(datestring), dT, day_of_year)


ihr  =  int(sec/3600.d0) !from 0 to 23


write(*,*) "year", year
write(*,*) "day_of_year", day_of_year
write(*,*) "ihr", ihr

call sfcsolz(year, day_of_year, ihr, solz)

!call sfcsolz10_14(year, day_of_year, ihr, solz)

write(*,*) 'solz', solz

call getrmud(solz,rmud)

write(*,*) 'rmud', rmud

MODE = 0

V_POSITION = "AVERAGE"

Ed_0m(:,1,1) = Ed_OASIM(:)
Es_0m(:,1,1) = Es_OASIM(:)

bottom = jpk

do i=1,33 ! PAR RANGE
     write(*,*) lam(i), Ed_0m(i,1,1),Es_0m(i,1,1)
enddo

call edeseu_IOP(MODE,V_POSITION,bottom,e3t(:,1,1),Ed_0m(:,1,1),Es_0m(:,1,1),CHLz,CDOMz,NAPz,rmud,Edz,Esz,Euz,Eu_0m(:),PARz)

do i=1,33 ! PAR RANGE
     write(*,*) lam(i), Edz(1,i),Esz(1,i),Euz(1,i), Eu_0m(i)
enddo

call WRITE_AVE_edeseu(TRIM(OUTPUT_FILE), nc, jpk, 'Edz', 'Esz', 'Euz', Ed_0m, Es_0m, Eu_0m, Edz, Esz, Euz)


call CLOSE_AVE_edeseu(TRIM(OUTPUT_FILE),nc)

end program
