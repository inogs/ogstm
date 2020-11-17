      MODULE DA_MEM
      USE myalloc

      IMPLICIT NONE
      public

      CHARACTER(LEN=35) CHLSUP_FOR_DA
      INTEGER DA_Nprocs
      INTEGER TREd_procs_per_node
      INTEGER max_procs_per_one_node
      INTEGER AssimilationLevels
      character(LEN=200) satfile_suffix
      character(LEN=40 ) satvarname

      CHARACTER(LEN=3) ::varlistDA(17)
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: CHL_SUP
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:):: CHLtot
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: CHLsat,VAR2D
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: ERRsat, MISFIT,cMISFIT
      REAL(4), allocatable, dimension (:,:,:) :: CORR, CORR_c, FACTOR
      REAL(4), allocatable, dimension(:,:,:) :: tottrnDA

      !REAL(4), allocatable, dimension (:,:,:) :: FACTORN, FactorNN,FactorNP,FactorNC,FactorNS
      REAL(4), allocatable, dimension (:,:,:) :: CORRN,CORRP,CORRC,CORRS
      REAL(4), allocatable,dimension(:,:,:) :: Nfraction, Pfraction,Cfraction, CHLfraction, Sfraction ! variabili 3D
      REAL(4), allocatable,dimension(:,:,:) :: N_CHL, P_CHL, N_C, P_C ! rapporti  3D
      REAL(4), allocatable,dimension(:,:,:) :: CHL_C, S_C             ! rapporti  3D


      CONTAINS

      SUBROUTINE DA_INIT
      varlistDA( 1)='P1l'
      varlistDA( 2)='P2l'
      varlistDA( 3)='P3l'
      varlistDA( 4)='P4l'

      varlistDA( 5)='P1c'
      varlistDA( 6)='P2c'
      varlistDA( 7)='P3c'
      varlistDA( 8)='P4c'

      varlistDA( 9)='P1n'
      varlistDA(10)='P2n'
      varlistDA(11)='P3n'
      varlistDA(12)='P4n'

      varlistDA(13)='P1p'
      varlistDA(14)='P2p'
      varlistDA(15)='P3p'
      varlistDA(16)='P4p'

      varlistDA(17)='P1s'


      if (myrank==0) then
        ALLOCATE (CHL_SUP(jpiglo,jpjglo))
        ALLOCATE ( CHLsat(jpiglo,jpjglo),  VAR2D(jpiglo,jpjglo) )
        ALLOCATE ( ERRsat(jpiglo,jpjglo), MISFIT(  jpiglo,jpjglo), cMISFIT(  jpiglo,jpjglo)  )
        ALLOCATE (CHLtot(jpiglo, jpjglo, jpk))

        ALLOCATE ( CORR(jpiglo, jpjglo, jpk), CORR_c(jpiglo, jpjglo, jpk), FACTOR(jpiglo, jpjglo, jpk))
        !ALLOCATE ( FACTORN (jpiglo, jpjglo, jpk))
        ALLOCATE ( CORRN(jpiglo, jpjglo, jpk))
        ALLOCATE ( CORRP(jpiglo, jpjglo, jpk))
        ALLOCATE ( CORRC(jpiglo, jpjglo, jpk))
        ALLOCATE ( CORRS(jpiglo, jpjglo, jpk))

        ALLOCATE (  N_CHL(jpiglo, jpjglo, jpk), P_CHL(jpiglo, jpjglo, jpk))
        ALLOCATE (  N_C(jpiglo, jpjglo, jpk), P_C(jpiglo, jpjglo, jpk))
        ALLOCATE (  CHL_C(jpiglo, jpjglo, jpk), S_C(jpiglo, jpjglo, jpk))
        ALLOCATE (   Nfraction(jpiglo, jpjglo, jpk), Pfraction(jpiglo, jpjglo, jpk), Cfraction(jpiglo, jpjglo, jpk))
        ALLOCATE ( CHLfraction(jpiglo, jpjglo, jpk), Sfraction(jpiglo, jpjglo, jpk))

        !allocate(tottrnDA(jpiglo, jpjglo, jpk))      ;  tottrnDA = huge(tottrnDA(1,1,1))

      endif
      END SUBROUTINE DA_INIT




        SUBROUTINE readNetCDF_2dvar(fileNetCDF,varname,im,jm,MATRIX)
        use netcdf
        implicit none
        character fileNetCDF*(*) ,varname*(*)
        integer ncid, stat, VARid, mycount
        integer im,jm
        real MATRIX(im,jm)

        mycount = 0

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid); call handle_err1(stat,mycount,fileNetCDF )
        stat = nf90_inq_varid (ncid, varname, VARid)
        call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat,mycount,fileNetCDF )
        stat = nf90_get_var (ncid,VARid,MATRIX)
        call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat,mycount,fileNetCDF )

        stat = nf90_close(ncid)                         ; call handle_err1(stat,mycount,fileNetCDF )


        end SUBROUTINE readNetCDF_2dvar

      SUBROUTINE readNetCDF_3dvar(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      real(4) M(jpiglo,jpjglo,jpk)

      counter = 0;
      start    = (/1,       1,       1,  1/)
      thecount = (/jpiglo,  jpjglo, jpk, 1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      END SUBROUTINE readNetCDF_3dvar


      SUBROUTINE read_corr(fileNetCDF,varname, M,jpk_corr)
      USE myalloc
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer jpk_corr
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      real(4) M(jpiglo,jpjglo,jpk)
      real(4) Mr(jpiglo,jpjglo,jpk_corr)

      counter = 0;
      start    = (/1,       1,       1,       1/)
      thecount = (/jpiglo,  jpjglo, jpk_corr, 1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,Mr,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)


      M = 0.0;
      M(:,:,1:jpk_corr) = Mr;
      END SUBROUTINE read_corr



      ! *************************************************
      ! ** true if the datestring is in the restart list
      ! *************************************************
      LOGICAL FUNCTION IsaDAvar(string)
      IMPLICIT NONE
      CHARACTER(LEN=3), INTENT(IN) :: string
          ! LOCAL
          INTEGER I

          IsaDAvar = .false.


          DO I=1,17
          if (varlistDA(I).eq.string) THEN
            IsaDAvar = .true.
            CYCLE
          endif
          ENDDO

      END FUNCTION IsaDAvar


      LOGICAL FUNCTION IsaCHLvar(string)
      IMPLICIT NONE
      CHARACTER(LEN=3), INTENT(IN) :: string
          ! LOCAL
          INTEGER I

          IsaCHLvar = .false.


          DO I=1,4
          if (varlistDA(I).eq.string) THEN
            IsaCHLvar = .true.
            CYCLE
          endif
          ENDDO

      END FUNCTION IsaCHLvar



      END MODULE
