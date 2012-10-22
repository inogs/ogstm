       MODULE modul_param_mpp

       USE modul_param
       
       IMPLICIT NONE

       public

CCC---------------------------------------------------------------------
CCC                      PARAMETER MPP
CCC                     ****************
CCC  PURPOSE :
CCC  ---------
CCC     Include mpp parameter
CCC 
CC   MODIFICATIONS :
CC   -------------
CC      original : 94 (j. Escobar, M. Imbard)
CC      modification for SHMEM : (A.M. Treguier)
CC      modification for MPI : 98 (l. Colombet)
CC
CCC---------------------------------------------------------------------
CCC  OPA8, IDRIS
CCC---------------------------------------------------------------------
CC
CC The processor number is a required power of two :
CC                       1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,...
CC MPP dimension
CC -------------
CC      ndimmax       : maximun dimension
CC      ndim_mpp      : dimension for this simulation
CC
CC-CC      INTEGER, parameter :: nprocmax=2**10,
CC-CC     &                      ndim_mpp=jpnij
      INTEGER :: nprocmax,
     &           ndim_mpp
C      PARAMETER(nprocmax=2**10)
C      PARAMETER(ndim_mpp=jpnij)
CC
#ifdef key_mpp_mpi
CC
CC MPI variable definition
CC -----------------------
CC  
CC not  PARAMETER
CC
#  else
CC
CC PVM and SHMEM VERSION
CC
CC PVM file definition
CC -------------------
CC      executable    : executable name 
CC      opaall        : group name 
CC      simfile       : file name 
CC

      CHARACTER(LEN=47) executable
      CHARACTER opaall*(*)
      CHARACTER(LEN=80) simfile
      PARAMETER(executable='opa')
      PARAMETER(simfile='pvm3_ndim')
#ifdef key_mpp_shmem
      PARAMETER(opaall="")
#else
      PARAMETER(opaall="opaall")
#endif
CC
CC PVM control
CC -----------
CC      mynode_print  : flag for print, mynode   routine
CC      mpprecv_print : flag for print, mpprecv  routine
CC      mppsend_print : flag for print, mppsend  routine
CC      mppsync_print : flag for print, mppsync  routine
CC      mpprsum_print : flag for print, mpprsum  routine
CC      mppisum_print : flag for print, mppisum  routine
CC      mpprisl_print : flag for print, mpprisl  routine
CC      mppimin_print : flag for print, mppimin  routine
CC      mpparent_print  : flag for print, mpparent routine
CC
      INTEGER mynode_print
      INTEGER mpprecv_print,mppsend_print,mppsync_print
      INTEGER mpprsum_print,mpprisl_print,mppimin_print 
      INTEGER mppisum_print,mpprmax_print,mpprmin_print
      INTEGER mpparent_print
      PARAMETER ( mynode_print  =0)
      PARAMETER ( mpprecv_print =0)
      PARAMETER ( mppsend_print =0)
      PARAMETER ( mppsync_print =0)
      PARAMETER ( mpprsum_print =0)
      PARAMETER ( mpprisl_print =0)
      PARAMETER ( mppimin_print =0)
      PARAMETER ( mppisum_print =0)
      PARAMETER ( mpprmax_print =0)
      PARAMETER ( mpprmin_print =0)
      PARAMETER ( mpparent_print=0)
CC
CC Variable definition
CC -------------------
CC      jpvmreal(8) : codage propre a pvm3 des real(8)4 et real(8)8 
CC                ( real(8)4 = 4, real(8)8 = 6)                      
CC      jpvmint  : codage propre a pvm3 des integer4 et integer8 
CC                ( integer4 = 3, integer8 = 21 : ext Cray)   
CC
#    ifdef key_mpp_pvm
      INTEGER jpvmreal(8),jpvmint
!      PARAMETER (jpvmreal(8)=6,jpvmint=3)
#      else
      INTEGER jpvmreal(8),jpvmint
!      PARAMETER (jpvmreal(8)=6,jpvmint=21)
#    endif
#    ifdef key_mpp_shmem
C
C SHMEM VERSION
C
C   Maximum  dimension of array to sum on the processors
C ----------------------
C
C
      INTEGER jpmppsum,jpmsec,jpmpplat
      PARAMETER(jpmsec=50000,jpmpplat=30)
      PARAMETER (jpmppsum=max(jpisl*jpisl,jpmpplat*jpk,jpmsec))
C
#    endif
#endif

        END MODULE
