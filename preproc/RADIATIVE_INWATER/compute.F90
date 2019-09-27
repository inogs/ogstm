program compute

IMPLICIT NONE

integer, parameter :: nlt=33
integer, parameter :: nchl=5
integer, parameter :: km=7
logical            :: first
CHARACTER(len=128) :: INPUT_OASIM_FILE, INPUT_ACDOM_FILE, INPUT_APHY_FILE, INPUT_ANAP_FILE, INPUT_PFT_FILE
CHARACTER(len=128) :: OUTPUT_FILE
CHARACTER(len=32)  :: time_string
integer :: i,k,chl,time_2hr,nc
integer :: lam(nlt)
real(8) :: aw(nlt),bw(nlt)
real(8) :: ac(nchl,nlt),bc(nchl,nlt)
real(8) :: rmud
real(8) :: Ed(nlt),Es(nlt)
real(8) :: H(km)
real(8) :: P(km,nchl)
real(8) :: acdom(km,nlt), aphy(km,nlt), anap(km,nlt)
real(8) :: CHLT(km)
real(8) :: zenith_angle

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

call OPEN_AVE_edeu(OUTPUT_FILE, 'Edz', 'Esz', 'Euz', nc)

call lidata(lam,aw,bw,ac,bc)

do i=1,nlt
   write(*,*) lam(i), aw(i), bw(i)
enddo

! H thikcness of layers !!!

H(1) = 4.0
H(2) = 5.0
H(3) = 10.0
H(4) = 20.0
H(5) = 10.0
H(6) = 20.0
H(7) = 10.0

P(:,:) = 0.0
call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"diatoms", P(:,1))
call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"chlorophytes", P(:,2))
call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"cyanobacteria", P(:,3))
call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"coccolitophores", P(:,4))
call readnc_boussole_PFT(TRIM(INPUT_PFT_FILE),"dinoflagellates", P(:,5))

CHLT(:) = 0.

do chl=1,nchl
        CHLT(:) = P(:,chl) + CHLT(:)
enddo

write(*,*) 'CHLT', CHLT

call readnc_boussole_acdom(TRIM(INPUT_ACDOM_FILE),'acdom', 1, acdom)

call readnc_boussole_acdom(TRIM(INPUT_APHY_FILE), 'aphy', 1, aphy)

call readnc_boussole_acdom(TRIM(INPUT_ANAP_FILE), 'anap', 1, anap)

first = .TRUE.
do time_2hr=1,12

    call readnc_boussole(trim(INPUT_OASIM_FILE),'Ed_0m', time_2hr, Ed)

    call readnc_boussole(trim(INPUT_OASIM_FILE),'Es_0m', time_2hr, Es)

    zenith_angle = -90.+ 180./24.* (2.0*real(time_2hr,8) - 1.0)

    write(*,*) "============time 2hr ", time_2hr, " ================"
    write(*,*) "zenith_angle ", zenith_angle

    do i=1,33 ! PAR RANGE
        write(*,*) lam(i), Ed(i),Es(i)
    enddo

    call getrmud(zenith_angle,rmud)

    call edeu(first,aw,bw,ac,bc,Ed,Es,H,P,aphy,anap,acdom,rmud,OUTPUT_FILE,nc,time_2hr)

enddo

call CLOSE_AVE_edeu(OUTPUT_FILE,nc)

end program
