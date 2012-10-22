      MODULE DIA_mem
          USE myalloc


          implicit none

          INTEGER FsizeGlo, Fsize
          INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: INDFluxGlo, INDflxDUMP
          INTEGER(4), allocatable, DIMENSION(:,:)   :: flx_ridxt

          REAL(8),    ALLOCATABLE, DIMENSION(:,:,:) :: diaflx
          REAL(8),    ALLOCATABLE, DIMENSION(:,:)   :: MflxDump

          CONTAINS

! **********************************************************
          SUBROUTINE alloc_dia
          USE NETCDF
            IMPLICIT NONE
            logical existFile
            INQUIRE(FILE='Fluxes.nc', EXIST=existFile)

            if (existFile) then
               call getDimension('Fluxes.nc','n',FsizeGlo)
            else
               FsizeGlo = 0
            endif

            if (FsizeGlo.gt.0) then
               allocate(INDFluxGlo(FsizeGlo))

            endif
            mem_all = FsizeGlo*4


          END SUBROUTINE alloc_dia


! **********************************************************

       SUBROUTINE alloc_DIA_local_flx()

       allocate(flx_ridxt (Fsize,4        ))
       allocate(INDflxDUMP(Fsize          ))
       allocate(MflxDump  (7, Fsize        ))
       allocate(diaflx    (Fsize, jptra, 7))
       diaflx = 0
      END SUBROUTINE alloc_DIA_local_flx


      END MODULE DIA_mem
