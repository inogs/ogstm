program compute

USE myalloc, ONLY: jpi,jpj,jpk,glamt,gphit
USE OPT_mem
USE TIME_MANAGER

IMPLICIT NONE

logical            :: first
CHARACTER(len=128) :: INPUT_OASIM_FILE, INPUT_ACDOM_FILE, INPUT_APHY_FILE, INPUT_ANAP_FILE, INPUT_PFT_FILE
CHARACTER(len=128) :: OUTPUT_FILE
CHARACTER(len=32)  :: time_string
character(LEN=17)  ::  datestring
integer :: i,k,chl,time_2hr,nc
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
double precision, allocatable :: Edz(:,:),Esz(:,:),Euz(:,:),PARz(:,:)
double precision, allocatable :: CHLz(:,:),CDOMz(:),NAPz(:)
double precision :: Eu_0m(nlt)



datestring='20000101-12:00:00'

dT =1800.0D0

jpi  = 1
jpj  = 1

call getarg(1, INPUT_OASIM_FILE)

write (*,*) "INPUT_OASIM_FILE= ", INPUT_OASIM_FILE

call getarg(2, INPUT_ACDOM_FILE)

write (*,*) " INPUT_ACDOM_FILE= ", INPUT_ACDOM_FILE

call getarg(3, INPUT_APHY_FILE)

write (*,*) " INPUT_APHY_FILE= ", INPUT_APHY_FILE

call getarg(4, INPUT_ANAP_FILE)

write (*,*) " INPUT_ANAP_FILE= ", INPUT_ANAP_FILE

call getarg(5, INPUT_PFT_FILE)

write (*,*) " INPUT_PFT_FILE= ", INPUT_PFT_FILE

call getarg(6, OUTPUT_FILE)

write (*,*) " OUTPUT_FILE= ", OUTPUT_FILE

!call OPEN_AVE_edeu(OUTPUT_FILE, 'Edz', 'Esz', 'Euz', nc)



jpk = 7

call  myalloc_OPT()

do i=1,nlt
   write(*,*) lam(i), aw(i), bw(i)
enddo

! H thikcness of layers !!!

allocate(e3t(jpk,1,1))
allocate(Edz(jpk,nlt),Esz(jpk,nlt),Euz(jpk,nlt),PARz(jpk,nchl+1))
allocate(CHLz(jpk,nchl),CDOMz(jpk),NAPz(jpk))
allocate(glamt(jpj,jpi),gphit(jpj,jpi))

e3t(1,1,1) = 4.0
e3t(2,1,1) = 5.0
e3t(3,1,1) = 10.0
e3t(4,1,1) = 20.0
e3t(5,1,1) = 10.0
e3t(6,1,1) = 20.0
e3t(7,1,1) = 10.0

gphit(1,1) = 35.0D0

!call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"diatoms", P(:,1))
!call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"chlorophytes", P(:,2))
!call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"cyanobacteria", P(:,3))
!call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"coccolitophores", P(:,4))
!call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"dinoflagellates", P(:,5))

CHLz(:,:) = 1.
CDOMz(:) = 1.
NAPz(:) = 1.


!call readnc_boussole_acdom(TRIM(INPUT_ACDOM_FILE),'acdom', 1, acdom)

!call readnc_boussole_acdom(TRIM(INPUT_APHY_FILE), 'aphy', 1, aphy)

!call readnc_boussole_acdom(TRIM(INPUT_ANAP_FILE), 'anap', 1, anap)

!call readnc_boussole(trim(INPUT_OASIM_FILE),'Ed_0m', time_2hr, Ed)

!call readnc_boussole(trim(INPUT_OASIM_FILE),'Es_0m', time_2hr, Es)


call read_date_string(datestring, year, month, day, sec)

call tau2julianday(datestringToTAU(datestring), dT, day_of_year)

call sfcsolz(year, day_of_year, ihr, solz)

call getrmud(solz,rmud)

MODE = 0

V_POSITION = "AVERAGE"

bottom = 1

Ed_0m(:,1,1) = 1.0d0
Es_0m(:,1,1) = 1.0D0

do i=1,33 ! PAR RANGE
     write(*,*) lam(i), Ed_0m(i,1,1),Es_0m(i,1,1)
enddo

call edeseu(MODE,V_POSITION,bottom,e3t(:,1,1),Ed_0m(:,1,1),Es_0m(:,1,1),CHLz,CDOMz,NAPz,rmud,Edz,Esz,Euz,Eu_0m(:),PARz)

do i=1,33 ! PAR RANGE
     write(*,*) lam(i), Edz(1,i),Esz(1,i),Euz(1,i), Eu_0m(i)
enddo


!call CLOSE_AVE_edeu(OUTPUT_FILE,nc)

end program
