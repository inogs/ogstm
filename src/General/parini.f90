      SUBROUTINE parini(ME)
!---------------------------------------------------------------------
!
!                       ROUTINE parini
!                     ******************
!
!  Purpose :
!  --------
!     Reads Domain Decomposition from Dom_Dec_jp?_ascii files
!           Domain size          from meshmask.nc


       USE myalloc
       USE myalloc_mpp
       IMPLICIT NONE

! local declarations
! ==================
      INTEGER ji,jj,nn,ME

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
      jpkb   = jpk


! epascolo warning
      IF(lwp) THEN
      WRITE(numout,*) 'Dimension_Med_Grid'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpjglo  : first dimension of global domain --> j ',jpjglo
      WRITE(numout,*) ' jpiglo  : second  dimension of global domain --> i ',jpiglo
      WRITE(numout,*) ' jpk     : number of levels           > or = jpk   ',jpk
      WRITE(numout,*) ' jpkb    : first vertical layers where biology is active > or = jpkb   ',jpkb
      WRITE(numout,*) ' '
      ENDIF


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
