      MODULE DIA_mem
          USE myalloc


          implicit none

          INTEGER FsizeGlo, Fsize, FsizeMax
          INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: INDFluxGlo, INDflxDUMP, INDflxDUMPglo, INDflxBuff, INDflxDUMPZERO
          INTEGER(4), allocatable, DIMENSION(:,:)   :: flx_ridxt

          double precision,    ALLOCATABLE, DIMENSION(:,:,:) :: diaflx
          double precision,    ALLOCATABLE, DIMENSION(:,:)   :: MflxDumpGlo, diaflxBuff
          real,    ALLOCATABLE, DIMENSION(:,:)   :: MflxDumpGlo_Float
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
       allocate(diaflx        (7, Fsize, jptra   ))
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



      subroutine clean_memory_dia()
          
          if (allocated(INDFluxGlo)) then
          deallocate(INDFluxGlo)
          endif

          ! conditional deallocation is chosen only in case of
          ! conditional allocation

          if (allocated(flx_ridxt)) then
              deallocate(flx_ridxt)
          endif

          if (allocated(INDflxDUMP)) then
              deallocate(INDflxDUMP)
          endif

          if (allocated(diaflx)) then
              deallocate(diaflx)
          endif

          if (allocated(INDflxDUMPZERO)) then
              deallocate(INDflxDUMPZERO)
          endif

          if (allocated(MflxDumpGlo)) then
              deallocate(MflxDumpGlo)
          endif

          if (allocated(MflxDumpGlo_FLOAT)) then
              deallocate(MflxDumpGlo_FLOAT)
          endif

          if (allocated(INDflxDUMPglo)) then
              deallocate(INDflxDUMPglo)
          endif

          deallocate(INDflxBuff)
          deallocate(diaflxBuff)

      end subroutine clean_memory_dia



      END MODULE DIA_mem
