      SUBROUTINE parini(ME)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE parini
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Read and print options for the run (namelist)
CC
CC   Input : 
CC   ------
CC      the namelist file ( unit numnam ) :
CC
CC   Output :
CC   -------

       USE myalloc
       USE myalloc_mpp
       IMPLICIT NONE

CC local declarations
CC ==================
      INTEGER ji,jj,nn,ME
      INTEGER ilu, iused(1,100)
      INTEGER ibloc, ilseq

      CHARACTER(LEN=32) clname
      CHARACTER(LEN=21) clold, clfor, clseq, clnew, cldir, clunf

      NAMELIST/MPI_config/ mpi_pack_size
      NAMELIST/Domain_Characteristic/  jperio
      NAMELIST/Number_Fluxes/ jpflx, jptaux, jptauy, jpwind, jpemp,jpkef, jpice, jpqsr
CC----------------------------------------------------------------------

C 0. Initializations
C ------------------

C ... Open specifier
      clold = 'OLD'
      clnew = 'NEW'
      clfor = 'FORMATTED'
      clunf = 'UNFORMATTED'
      clseq = 'SEQUENTIAL'
      cldir = 'DIRECT'

C ... Sequential value
      ilseq = 1

C ... computation of the record length for direct access file
C     this length depend of 512 for the t3d machine
      ibloc = 512
      ilu = 0
      DO ji = 1, 100
        iused(1,ji) = 0
      END DO

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' routine parini'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' initialization of Domain'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist.init'
          WRITE(numout,*) ' ********'
          WRITE(numout,*) ' '
      ENDIF


      numnam = 1004
      clname = 'namelist.init'
      CALL ctlopn( numnam, clname, clold, clfor, clseq,ilseq, ilu, iused, numout, lwp, 1)


      call COUNTLINE ('Dom_Dec_jpi.ascii', jpni)
      call COUNTWORDS('Dom_Dec_jpi.ascii', jpnj)
      jpnij = jpni*jpnj
      jpreci = 1
      jprecj = 1



      IF(lwp) THEN
      WRITE(numout,*) ' '
      WRITE(numout,*) 'Dom_Size'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' number of processors following i : jpni   = ', jpni
      WRITE(numout,*) ' number of processors following j : jpnj   = ', jpnj
      WRITE(numout,*) ' '
      WRITE(numout,*) ' local domains : < or = jpni x jpnj number of processors   = ', jpnij
      WRITE(numout,*) ' number of lines for overlap  jpreci   = ',jpreci
      WRITE(numout,*) ' number of lines for overlap  jprecj   = ',jprecj
      WRITE(numout,*) ' '

      ENDIF


      call getDimension('meshmask.nc','x',jpiglo)
      call getDimension('meshmask.nc','y',jpjglo)
      call getDimension('meshmask.nc','z',jpk)
      jpidta = jpiglo
      jpjdta = jpjglo
      jpkdta = jpk
      jpkb   = jpk
      jpilec = jpiglo
      jpjlec = jpjglo
      jpklec = jpk



      IF(lwp) THEN
      WRITE(numout,*) 'Dimension_Med_Grid'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpiglo  : first  dimension of global domain --> i ',jpiglo
      WRITE(numout,*) ' jpjglo  : second dimension of global domain --> j ',jpjglo
      WRITE(numout,*) ' jpk     : number of levels           > or = jpk   ',jpk
      WRITE(numout,*) ' jpidta  : first  dimension of global domain --> i ',jpidta
      WRITE(numout,*) ' jpjdta  : second dimension of global domain --> j ',jpjdta
      WRITE(numout,*) ' jpkdta  : number of levels           > or = jpk   ',jpkdta
      WRITE(numout,*) ' jpkb    : first vertical layers where biology is active > or = jpkb   ',jpkb
      WRITE(numout,*) ' '
      ENDIF

      jpkmod=(1+(jpk-1)/jpnij)



      allocate(ilcit(jpni, jpnj)) ; ilcit = huge(ilcit(1,1))
      allocate(ilcjt(jpni, jpnj)) ; ilcjt = huge(ilcjt(1,1))

      open(3333,file='Dom_Dec_jpi.ascii', form='formatted')
      open(3334,file='Dom_Dec_jpj.ascii', form='formatted')
      
      read(3333,*) ((ilcit(ji,jj), jj=1,jpnj),ji=1,jpni)
      read(3334,*) ((ilcjt(ji,jj), jj=1,jpnj),ji=1,jpni)
      
      close(3333)
      close(3334)

      do nn =1, jpni*jpnj
        if(ME+1 .EQ. nn) then
          ji = 1 + mod(nn -1, jpni)
          jj = 1 + (nn -1)/jpni
          jpi =  ilcit(ji,jj) 
          jpj =  ilcjt(ji,jj)
        endif
      enddo

CC Other dimension parameters
CC --------------------------
CC      jpim1   :  jpi - 1
CC      jpjm1   :  jpj - 1
CC      jpkm1   :  jpk - 1
CC      jpij    :  jpi x jpj

      jpim1=jpi-1
      jpjm1=jpj-1
      jpkm1=jpk-1
      jpij=jpi*jpj
      jpkbm1=jpkb-1

C ... Namelist :MPI_config
C     Deafault value
      mpi_pack_size = 1

      REWIND( numnam )
      READ  ( numnam, MPI_config )

      IF(lwp) THEN
      WRITE(numout,*) 'MPI_config'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' mpi_pack_size  :size of mpi packet -->  ',mpi_pack_size
      WRITE(numout,*) ' '
      ENDIF

C
C ... namelist.init: Domain_Characteristic
C
C
      REWIND( numnam )
      READ  ( numnam, Domain_Characteristic )


      IF(lwp) THEN
      WRITE(numout,*) 'Domain_Characteristic'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' lateral cond. type for the global domain (4, 3, 2, 1 or 0)',jperio
      ENDIF

C ..  Namelist :  Number of fluxes


      REWIND( numnam )
      READ  ( numnam, Number_Fluxes )

      IF(lwp) THEN
      WRITE(numout,*) 'Number of fluxes'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpflx   : number of fluxes',jpflx
      WRITE(numout,*) ' jptaux  :  zonal wind stress',jptaux
      WRITE(numout,*) ' jptauy  :  meridional wind stress',jptauy
      WRITE(numout,*) ' jpemp   :  E - P in mm/day',jpemp
      WRITE(numout,*) ' jpqsr   :  solar radiation',jpqsr
      WRITE(numout,*) ' jpkef   :  Extinction Coefficient',jpkef
      WRITE(numout,*) ' jpilec  : first horizontal dimension > or = jpi ',jpilec
      WRITE(numout,*) ' jpjlec  : second                     > or = jpj ',jpjlec
      WRITE(numout,*) ' jpklec  : number of levels           > or = jpk ',jpklec
      WRITE(numout,*) ' '
      ENDIF


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
      SUBROUTINE COUNTWORDS(filename,n)
        IMPLICIT NONE
        CHARACTER*(*) filename
        INTEGER N
        ! local
        INTEGER I
        CHARACTER(LEN=1024) str, str_blank


        open(unit=21,file=filename, form='formatted')
        read(21,'(A)') str
        close(21)

        str_blank=' '//trim(str)
        N=0
        do i = 1,len(trim(str))
         if ((str_blank(i:i).eq.' ').and.(str_blank(i+1:i+1).ne.' ') )  N=N+1
        enddo

      END SUBROUTINE COUNTWORDS


! ***************************************************************
      END SUBROUTINE PARINI
