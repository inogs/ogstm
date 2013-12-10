       MODULE myalloc_mpp

       USE modul_param_mpp

       IMPLICIT NONE

       public

CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/modulo_mpp.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC---------------------------------------------------------------------
CCC                      COMMON MPP
CCC                     *************
CCC  PURPOSE :
CCC  ---------
CCC     Include mpp file
CCC 
CC   MODIFICATIONS :
CC   -------------
CC      original : 94 (j. Escobar, M. Imbard)
CC      modification for MPI : 98 (l. Colombet)
CC
CCC---------------------------------------------------------------------
CCC  OPA8, IDRIS
CCC---------------------------------------------------------------------
CC
#ifdef key_mpp_mpi
CC
CC MPI variable definition
CC -----------------------
CC 
#include <mpif.h>
CC 
CC    size : number of process
CC    rank : process number  [ 0 - size-1 ]
CC
      INTEGER  size,rank
      INTEGER  :: mpi_pack_size =18
#  else
CC
CC PVM and SHMEM VERSION
CC
#include  <fpvm3.h>
CC
CC PVM variable definition
CC -----------------------
CC      npvm_ipas     : pvm initialization flag 
CC      npvm_mytid    : pvm tid
CC      npvm_me       : node number [ 0 - nproc-1 ]
CC      npvm_nproc    : real(8) number of nodes
CC      npvm_tids     : tids array [ 0 - nproc-1 ]
CC
      INTEGER npvm_ipas,npvm_mytid,npvm_me
      INTEGER npvm_nproc,npvm_inum
      INTEGER npvm_tids(0:nprocmax-1)
CC
      COMMON/compvm/npvm_ipas,npvm_mytid,npvm_me,
     $              npvm_nproc,npvm_inum,
     $              npvm_tids
CC
CC T3D variable definition
CC -----------------------
CC      nt3d_ipas     : pvm initialization flag
CC      nt3d_mytid    : pvm tid
CC      nt3d_me       : node number [ 0 - nproc-1 ]
CC      nt3d_nproc    : real(8) number of nodes
CC      nt3d_tids     : tids array [ 0 - nproc-1 ]
CC
      INTEGER nt3d_ipas,nt3d_mytid,nt3d_me
      INTEGER nt3d_nproc
      INTEGER nt3d_tids(0:nprocmax-1)
CC
      COMMON/comt3d/nt3d_ipas,nt3d_mytid,nt3d_me,
     $              nt3d_nproc,
     $              nt3d_tids
CC
#    ifdef key_mpp_shmem
CC
CC SHMEM VERSION
CC
#include <mpp/shmem.fh>
CC
CC real(8) sum reduction
CC
      INTEGER nrs1sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER nrs2sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      REAL(8) wrs1wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      REAL(8) wrs2wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      REAL(8) wrstab_shmem(jpmppsum)
      COMMON /mpprsum_local / wrstab_shmem, wrs1wrk_shmem,
     $    nrs1sync_shmem, wrs2wrk_shmem, nrs2sync_shmem
CC
CC minimum and maximum reduction
CC
      INTEGER ni1sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER ni2sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      REAL(8) wi1wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,jpmppsum/2+1))
      REAL(8) wi2wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,jpmppsum/2+1))
      REAL(8) wintab_shmem(jpmppsum), wi1tab_shmem(jpmppsum),
     $    wi2tab_shmem(jpmppsum)
      COMMON /mpprmimx_local /wintab_shmem, wi1tab_shmem, wi2tab_shmem  
     $    ,wi1wrk_shmem, ni1sync_shmem, wi2wrk_shmem, ni2sync_shmem
CC
CC value not equal zero for barotropic stream function around islands
CC
      INTEGER ni11sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER ni12sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER ni21sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER ni22sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      REAL(8) wi11wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      REAL(8) wi12wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      REAL(8) wi21wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      REAL(8) wi22wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      REAL(8) wiltab_shmem(jpmppsum), wi11tab_shmem(jpmppsum),
     $    wi12tab_shmem(jpmppsum), wi21tab_shmem(jpmppsum),
     $    wi22tab_shmem(jpmppsum)
CC
      INTEGER ni11wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      INTEGER ni12wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      INTEGER ni21wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      INTEGER ni22wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      INTEGER  niitab_shmem(jpmppsum), ni11tab_shmem(jpmppsum),
     $    ni12tab_shmem(jpmppsum)
CC
      COMMON /mpprisl_local /wiltab_shmem, wi11tab_shmem, wi12tab_shmem,  
     $    wi21tab_shmem, wi22tab_shmem,wi11wrk_shmem, ni11sync_shmem,
     $    wi12wrk_shmem, ni12sync_shmem,wi21wrk_shmem, ni21sync_shmem,
     $    wi22wrk_shmem, ni22sync_shmem,
     $    niitab_shmem , ni11tab_shmem , ni12tab_shmem,  
     $    ni11wrk_shmem, ni12wrk_shmem , ni21wrk_shmem
CC
      INTEGER nis1sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER nis2sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER nis1wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      INTEGER nis2wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      INTEGER nistab_shmem(jpmppsum)
      COMMON /mppisum_local / nistab_shmem, nis1wrk_shmem,
     $    nis1sync_shmem, nis2wrk_shmem, nis2sync_shmem
CC
CC integer sum reduction
CC
      INTEGER nil1sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER nil2sync_shmem(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER nil1wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      INTEGER nil2wrk_shmem(max(SHMEM_REDUCE_MIN_WRKDATA_SIZE,
     $    jpmppsum/2+1))
      INTEGER niltab_shmem(jpmppsum)
      COMMON /mppimin_local / niltab_shmem, nil1wrk_shmem,
     $    nil2sync_shmem, nil2wrk_shmem, nil1sync_shmem
#    endif
#endif

        END MODULE
