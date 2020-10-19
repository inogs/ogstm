MODULE NODE_NAME
        !---------------------------------------------------------------------
        !
        !                       MODULE MATRIX_VARS
        !                     ******************
        !
        !
        !  gcoidessa@inogs.it  developing 
        !  Purpose :
        !  ---------
        !  OBTAIN name of node of each rank
        !  -------
        !  subroutines: - NODE_NAME_FILL(): found node name
        !


        USE calendar
        USE myalloc
        USE IO_mem
        USE FN_mem
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module


        IMPLICIT NONE
        

        PUBLIC
        
        INTEGER :: lengt
        
        CHARACTER*(MPI_MAX_PROCESSOR_NAME) local_array

        CONTAINS

!------------------------------------------------------------------------------
        SUBROUTINE NODE_NAME_FILL()
        
        INTEGER :: IERROR
        
        CALL MPI_GET_PROCESSOR_NAME(local_array, lengt, IERROR)
        
        END SUBROUTINE NODE_NAME_FILL        
!---------------------------------------------------------------------------------
END MODULE
                


















                
