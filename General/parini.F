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
      INTEGER ji,jj,nn,ME!,jk
!      INTEGER iost, ivide1, ivide2
      INTEGER ilu, iused(1,100)
      INTEGER ilglo, ildta, ibloc, ilseq

      CHARACTER(LEN=32) clname
      CHARACTER(LEN=21) clold, clfor, clseq, clnew, cldir, clunf

      NAMELIST/Dom_Size/ jpni, jpnj, jpnij, jpreci, jprecj

      NAMELIST/Dimension_Med_Grid/ jpiglo, jpjglo, jpk, jpidta,jpjdta, jpkdta, jpkb
      NAMELIST/Domain_Characteristic/  jperio, jpisl, jpnisl
      NAMELIST/Number_Fluxes/ jpflx, jptaux, jptauy, jpwind, jpemp,
     &                        jpkef, jpice, jpqsr, jpilec, jpjlec,jpklec
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
C
C ... Sequential value
      ilseq = 1
C
C ... computation of the record length for direct access file
C     this length depend of 512 for the t3d machine
      ibloc = 512
      ildta = ibloc*((jpidta*jpjdta)/ibloc+1)*jpbyt
      ilglo = ibloc*((jpiglo*jpjglo)/ibloc+1)*jpbyt
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
      REWIND( numnam )
      clname = 'namelist.init'
      CALL ctlopn( numnam, clname, clold, clfor, clseq,ilseq, ilu, iused, numout, lwp, 1)

C 1. Read namlist files
C ---------------------

C ... Namelist Dom_Size :


      REWIND( numnam )
      READ  ( numnam, Dom_Size )

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
C
C ... Namelist :Dimension_Med_Grid


      REWIND( numnam )
      READ  ( numnam, Dimension_Med_Grid )

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
C
      jpkmod=(1+(jpk-1)/jpnij)

C ... Namelist.init :Matrix_Size
C
C
C

      allocate(ilcit(jpni, jpnj))
      allocate(ilcjt(jpni, jpnj))

      open(3333,file='Dom_Dec_jpi.ascii', form='formatted')
      open(3334,file='Dom_Dec_jpj.ascii', form='formatted')
      
      read(3333,*) ((ilcit(ji,jj), ji=1,jpni),jj=1,jpnj) 
      read(3334,*) ((ilcjt(ji,jj), ji=1,jpni),jj=1,jpnj) 
      
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
C
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

C
C ... namelist.init: Domain_Characteristic
C
C
      REWIND( numnam )
      READ  ( numnam, Domain_Characteristic )
C
C
      IF(lwp) THEN
      WRITE(numout,*) 'Domain_Characteristic'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' lateral cond. type for the global domain (4, 3, 2, 1 or 0)',jperio
      WRITE(numout,*) ' jpisl  : number of islands', jpisl
      WRITE(numout,*) ' jpnisl  : maximum number of points per island', jpnisl
      ENDIF
C
C ..  Namelist :  Number of fluxes
C
C
      REWIND( numnam )
      READ  ( numnam, Number_Fluxes )
C
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


      END SUBROUTINE PARINI
