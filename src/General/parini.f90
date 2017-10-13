      SUBROUTINE parini()
!---------------------------------------------------------------------
!
!                       ROUTINE parini
!                     ******************
!
!  Purpose :
!  --------
!     Reads Domain Decomposition from domdec.txt file
!           Domain size          from meshmask.nc


       USE myalloc
       USE modul_param
       IMPLICIT NONE

! local declarations
! ==================
      INTEGER ji,jj,nn
      INTEGER nproci, nprocj
      integer, allocatable, dimension(:,:) :: domdec

      call COUNTLINE ('domdec.txt', nn)
      if (nn.ne.mpi_glcomm_size) then
         if (lwp) write(*,*) 'domdec.txt is not compliant with MPI ranks '
         STOP
      endif

        narea = myrank+1
        allocate(domdec(mpi_glcomm_size,13))

        open(3333,file='domdec.txt', form='formatted')
        read(3333,*) ((domdec(ji,jj), jj=1,13),ji=1,mpi_glcomm_size)
        close(3333)

        if (domdec(narea,1).ne.myrank) write(*,*) 'ERROR'
!        if (lwp) THEN
!        do ji=1,mpi_glcomm_size
!        write(*,'(13I3)') domdec(ji,:)
!        enddo
!        endif
        nproci = maxval(domdec(:,2))+1
        nprocj = maxval(domdec(:,3))+1

        ji     = domdec(narea, 2)+1
        jj     = domdec(narea, 3)+1
        jpi    = domdec(narea, 4)
        jpj    = domdec(narea, 5)
        nimpp  = domdec(narea, 6)
        njmpp  = domdec(narea, 7)
        nbondi = domdec(narea, 8)
        nbondj = domdec(narea, 9)
        nowe   = domdec(narea,10)
        noea   = domdec(narea,11)
        nono   = domdec(narea,12)
        noso   = domdec(narea,13)

      jpreci = 1
      jprecj = 1

! -----  from inimpp -----------------------------
        nlci = jpi
        nlcj = jpj
        nldi= 1  +jpreci
        nlei=nlci-jpreci
        IF(ji.eq.1) nldi=1 ! western boundary without ghost cell
        IF(ji.eq.nproci) nlei=nlci
        nldj= 1  +jprecj
        nlej=nlcj-jprecj
        IF(jj.eq.1) nldj=1 ! south boundary without ghost cell
        IF(jj.eq.nprocj) nlej=nlcj
        nperio=0
! ------------------------------------------------
      IF(lwp) THEN
      WRITE(numout,*) ' '
      WRITE(numout,*) 'Dom_Size'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' number of lines for overlap  jpreci   = ',jpreci
      WRITE(numout,*) ' number of lines for overlap  jprecj   = ',jprecj
      WRITE(numout,*) ' '

      ENDIF


      call getDimension('meshmask.nc','x',jpiglo)
      call getDimension('meshmask.nc','y',jpjglo)
      call getDimension('meshmask.nc','z',jpk)
      jpkb   = jpk



      IF(lwp) THEN
      WRITE(numout,*) 'Dimension_Med_Grid'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpiglo  : first  dimension of global domain --> i ',jpiglo
      WRITE(numout,*) ' jpjglo  : second dimension of global domain --> j ',jpjglo
      WRITE(numout,*) ' jpk     : number of levels           > or = jpk   ',jpk
      WRITE(numout,*) ' jpkb    : first vertical layers where biology is active > or = jpkb   ',jpkb
      WRITE(numout,*) ' '
      ENDIF



      jpim1=jpi-1
      jpjm1=jpj-1
      jpkm1=jpk-1
      jpij=jpi*jpj
      jpkbm1=jpkb-1




      CLOSE(numnam)
      CONTAINS
! **************************************************************
      SUBROUTINE COUNTLINE(FILENAME,LINES)
          implicit none
          character FILENAME*(*)
          integer lines
          integer TheUnit

          TheUnit = 326

          lines=0
          OPEN(UNIT=TheUnit,file=FILENAME,status='old')
          DO WHILE (.true.)
       read(TheUnit, *, END=21)
       lines = lines+1
          ENDDO

21        CLOSE(TheUnit)

      END SUBROUTINE COUNTLINE
! **************************************************************
      END SUBROUTINE PARINI
