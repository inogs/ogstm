      MODULE DIA_mem
          USE myalloc


          implicit none

          INTEGER FsizeGlo, Fsize, FsizeMax
          INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: INDFluxGlo, INDflxDUMP, INDflxDUMPglo, INDflxBuff, INDflxDUMPZERO
          INTEGER(4), allocatable, DIMENSION(:,:)   :: flx_ridxt

          REAL(8),    ALLOCATABLE, DIMENSION(:,:,:) :: diaflx
          REAL(8),    ALLOCATABLE, DIMENSION(:,:)   :: MflxDumpGlo, diaflxBuff
          REAL(4),    ALLOCATABLE, DIMENSION(:,:)   :: MflxDumpGlo_Float
          logical existFileFluxes

          CONTAINS

! **********************************************************
          SUBROUTINE alloc_dia
          USE NETCDF
            IMPLICIT NONE

            INQUIRE(FILE='Fluxes.nc', EXIST=existFileFluxes)
            if (existFileFluxes) then
               call getDimension('Fluxes.nc','n',FsizeGlo)
            else
               FsizeGlo = 0
            endif

            if (FsizeGlo.gt.0) then

               allocate(INDFluxGlo(FsizeGlo)) 
       
               INDFluxGlo = huge(INDFluxGlo(1))

            endif
            mem_all = FsizeGlo*4


          END SUBROUTINE alloc_dia


! **********************************************************

       SUBROUTINE alloc_DIA_local_flx()

       allocate(flx_ridxt     (Fsize,4           )) 
       flx_ridxt  = huge(flx_ridxt(1,1))
       allocate(INDflxDUMP    (Fsize             ))  
       INDflxDUMP = huge(INDflxDUMP(1))
       allocate(diaflx        (Fsize, jptra, 7   )) 
       diaflx     = huge(diaflx(1,1,1))
       diaflx = 0
      END SUBROUTINE alloc_DIA_local_flx



      SUBROUTINE alloc_DIA_GLOBAL_flx()
      allocate(INDflxDUMPZERO(   Fsize             )) 
      INDflxDUMPZERO    = huge(INDflxDUMPZERO(1))
      allocate(MflxDumpGlo      (7, FsizeGlo       ))  
      MflxDumpGlo       = huge(MflxDumpGlo(1,1))
      allocate(MflxDumpGlo_FLOAT(7, FsizeGlo       ))  
      MflxDumpGlo_FLOAT = huge(MflxDumpGlo_FLOAT(1,1))
      allocate(INDflxDUMPglo    (FsizeGlo          )) 
      INDflxDUMPglo     = huge(INDflxDUMPglo(1))

      END SUBROUTINE alloc_DIA_GLOBAL_flx


      SUBROUTINE alloc_DIA_MPI_flx()
       allocate(INDflxBuff    (FsizeMax          )) 
       INDflxBuff = huge(INDflxBuff(1))
       allocate(diaflxBuff    (FsizeMax, 7))        
       diaflxBuff = huge(diaflxBuff(1,1))
      END SUBROUTINE alloc_DIA_MPI_flx

      END MODULE DIA_mem
