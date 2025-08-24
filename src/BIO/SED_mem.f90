       MODULE SED_mem

     USE modul_param 
       USE myalloc
       USE DIA_mem

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif


       IMPLICIT NONE

       public


      INTEGER :: dimen_jvsed
      INTEGER :: nsed=26
      INTEGER, allocatable :: sed_idx(:)
      INTEGER, allocatable :: jarr_sed(:,:),jarr_sed_flx(:,:)
      double precision, allocatable :: ztra(:,:)
      double precision, allocatable :: zwork(:,:,:)


!!!CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_SED()

#ifdef __OPENMP1
      INTEGER :: ntids, omp_get_max_threads
      EXTERNAL :: omp_get_max_threads
#else
      INTEGER :: ntids = 1
#endif
      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
#endif
       dimen_jvsed=0

       allocate(sed_idx(nsed))  
       sed_idx = huge(sed_idx(1))

       sed_idx(1)  = ppR6c
       sed_idx(2)  = ppR6n
       sed_idx(3)  = ppR6p
       sed_idx(4)  = ppR6s

       sed_idx(5)  = ppR8c
       sed_idx(6)  = ppR8n
       sed_idx(7)  = ppR8p
       sed_idx(8)  = ppR8s

       sed_idx(9)  = ppP1c
       sed_idx(10)  = ppP1n
       sed_idx(11)  = ppP1p
       sed_idx(12)  = ppP1s
       sed_idx(13)  = ppP1l

       sed_idx(14) = ppP2c
       sed_idx(15) = ppP2n
       sed_idx(16) = ppP2p
       sed_idx(17) = ppP2l

       sed_idx(18) = ppP3c
       sed_idx(19) = ppP3n
       sed_idx(20) = ppP3p
       sed_idx(21) = ppP3l

       sed_idx(22) = ppP4c
       sed_idx(23) = ppP4n
       sed_idx(24) = ppP4p
       sed_idx(25) = ppP4l
       sed_idx(26) = ppO5c
       allocate(jarr_sed(2, jpi*jpj))        
       jarr_sed     = huge(jarr_sed(1,1))
       allocate(jarr_sed_flx(jpk,jpi*jpj)) 
       jarr_sed_flx = huge(jarr_sed_flx(1,1)) 
       allocate( ztra(nsed,ntids))         
       ztra         = huge(ztra(1,1))
       allocate(zwork(jpk,nsed, ntids))   
       zwork        = huge(zwork(1,1,1))


#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_SED
      
      
      
      subroutine clean_memory_sed

          deallocate(sed_idx)
          deallocate(jarr_sed)
          deallocate(jarr_sed_flx)
          deallocate(ztra)
          deallocate(zwork)

      end subroutine clean_memory_sed



      END MODULE 
