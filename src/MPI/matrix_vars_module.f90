MODULE MATRIX_VARS
        !---------------------------------------------------------------------
        !
        !                       MODULE MATRIX_VARS
        !                     ******************
        !
        !
        !  gcoidessa@inogs.it  developing 
        !  Purpose :
        !  ---------
        !  OBTAIN MATRIX OF VARIABLES TO DUMP
        !  ---------
        !  Subroutines: 
        !       -ALLOCATE_MATRIX_VARS(): allocate matrix to use
        !       -POPULATE_MATRIX_VARS(): fill matrix to use with vars       
        !       -CLEAN_MATRIX_VARS()   


        USE calendar
        USE myalloc
        USE IO_mem
        USE FN_mem
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module

        USE dtype_procs_string_module        
        USE nodes_module

        IMPLICIT NONE
        
        CHARACTER(len=12) :: novars

        PUBLIC
        
        type(processor_string), allocatable :: matrix_state_1(:,:)
        type(processor_string), allocatable :: matrix_state_2(:,:)
        type(processor_string), allocatable :: matrix_diag_2d_1(:,:)
        type(processor_string), allocatable :: matrix_diag_2d_2(:,:)
        type(processor_string), allocatable :: matrix_diag_1(:,:)
        type(processor_string), allocatable :: matrix_diag_2(:,:)
        INTEGER :: matrix_state_1_row, matrix_state_2_row
        INTEGER :: matrix_diag_2d_1_row, matrix_diag_2d_2_row
        INTEGER :: matrix_diag_1_row, matrix_diag_2_row
        INTEGER :: matrix_col != nodes

        CONTAINS

!------------------------------------------------------------
        SUBROUTINE ALLOCATE_MATRIX_VARS()
        
        ALLOCATE (matrix_state_1(matrix_state_1_row, matrix_col))
        ALLOCATE (matrix_state_2(matrix_state_2_row, matrix_col))
        ALLOCATE (matrix_diag_2d_1(matrix_diag_2d_1_row, matrix_col))
        ALLOCATE (matrix_diag_2d_2(matrix_diag_2d_2_row, matrix_col))
        ALLOCATE (matrix_diag_1(matrix_diag_1_row, matrix_col))
        ALLOCATE (matrix_diag_2(matrix_diag_2_row, matrix_col))
        END SUBROUTINE ALLOCATE_MATRIX_VARS
!------------------------------------------------------------

!------------------------------------------------------------
        SUBROUTINE POPULATE_MATRIX_VARS()
        INTEGER :: i,j
        INTEGER :: counter,counter_high,trans_high,counter_dia2d,counter_dia2d_high,trans_dia2d_high,counter_dia,counter_dia_high,trans_dia_high
        counter = 0
        counter_high = 0
        counter_dia2d = 0
        counter_dia2d_high = 0
        counter_dia = 0
        counter_dia_high = 0

        ! define how many rows will be in the procs matrix
        
        matrix_col = nodes
        novars="novars_input"
        

        !IF (FREQ_GROUP.eq.2) THEN

        IF (MOD(jptra,nodes)==0)THEN
                matrix_state_2_row = (jptra/nodes)
        ELSE
                matrix_state_2_row = (jptra/nodes) + 1
        END IF
        
        !ELSE
        
        IF (MOD(jptra_high,nodes)==0)THEN
                matrix_state_1_row = (jptra_high/nodes)
        ELSE
                matrix_state_1_row = (jptra_high/nodes) + 1
        END IF

        !END IF
        IF (MOD(jptra_dia_2d,nodes)==0)THEN
                matrix_diag_2d_2_row = (jptra_dia_2d/nodes)
        ELSE
                matrix_diag_2d_2_row = (jptra_dia_2d/nodes) + 1
        END IF

        !ELSE

        IF (MOD(jptra_dia2d_high,nodes)==0)THEN
                matrix_diag_2d_1_row = (jptra_dia2d_high/nodes)
        ELSE
                matrix_diag_2d_1_row = (jptra_dia2d_high/nodes) + 1
        END IF

        IF (MOD(jptra_dia,nodes)==0)THEN
                matrix_diag_2_row = (jptra_dia/nodes)
        ELSE
                matrix_diag_2_row = (jptra_dia/nodes) + 1
        END IF

        !ELSE

        IF (MOD(jptra_dia_high,nodes)==0)THEN
                matrix_diag_1_row = (jptra_dia_high/nodes)
        ELSE
                matrix_diag_1_row = (jptra_dia_high/nodes) + 1
        END IF

        !allocate matrix

        CALL ALLOCATE_MATRIX_VARS()

        !POPULATE matrix_procs with string values;
        !will be used next, when will count the number of variable for each assignments

        !IF (FREQ_GROUP.eq.2) THEN
        DO i=1,matrix_state_2_row
                DO j=1,matrix_col
                        matrix_state_2(i,j)%var_name = novars
                END DO
        END DO
        !ELSE
        DO i=1,matrix_state_1_row
                DO j=1,matrix_col
                        matrix_state_1(i,j)%var_name = novars
                END DO
        END DO
        !END IF
        DO i=1,matrix_diag_2d_2_row
                DO j=1,matrix_col
                        matrix_diag_2d_2(i,j)%var_name = novars
                END DO
        END DO
        !ELSE
        DO i=1,matrix_diag_2d_1_row
                DO j=1,matrix_col
                        matrix_diag_2d_1(i,j)%var_name = novars
                END DO
        END DO
        
        DO i=1,matrix_diag_2_row
                DO j=1,matrix_col
                        matrix_diag_2(i,j)%var_name = novars
                END DO
        END DO
        !ELSE
        DO i=1,matrix_diag_1_row
                DO j=1,matrix_col
                        matrix_diag_1(i,j)%var_name = novars
                END DO
        END DO
        
        !insert variable string inside matrix of procs
        !control on how many variables are inserted by counter==jptra for matrix that are not fully covered

        !IF (FREQ_GROUP.eq.2) THEN
        DO i=1,matrix_state_2_row
                DO j=1,matrix_col
                        IF (counter==jptra)THEN
                                EXIT
                        ELSE
                                matrix_state_2(i,j)%var_name = ctrcnm(counter+1)
                                counter=counter + 1
                        END IF
                END DO
        END DO
        !ELSE
        DO i=1,matrix_state_1_row
                DO j=1,matrix_col
                        IF (counter_high==jptra_high)THEN
                                EXIT
                        ELSE    
                                trans_high = highfreq_table(counter_high+1)                             
                                matrix_state_1(i,j)%var_name =ctrcnm(trans_high)
                                counter_high=counter_high + 1
                        END IF
                END DO
        END DO
        !END IF
        DO i=1,matrix_diag_2d_2_row
                DO j=1,matrix_col
                        IF (counter_dia2d==jptra_dia_2d)THEN
                                EXIT
                        ELSE
                                matrix_diag_2d_2(i,j)%var_name = dianm_2d(counter_dia2d+1)
                                counter_dia2d=counter_dia2d + 1
                        END IF
                END DO
        END DO
        !ELSE
        DO i=1,matrix_diag_2d_1_row
                DO j=1,matrix_col
                        IF (counter_dia2d_high==jptra_dia2d_high)THEN
                                EXIT
                        ELSE
                                trans_dia2d_high = highfreq_table_dia2d(counter_dia2d_high+1)                     
                                matrix_diag_2d_1(i,j)%var_name =dianm_2d(trans_dia2d_high)
                                counter_dia2d_high=counter_dia2d_high + 1
                        END IF
                END DO
        END DO

        DO i=1,matrix_diag_2_row
                DO j=1,matrix_col
                        IF (counter_dia==jptra_dia)THEN
                                EXIT
                        ELSE
                                matrix_diag_2(i,j)%var_name = dianm(counter_dia+1)
                                counter_dia=counter_dia + 1
                        END IF
                END DO
        END DO
        !ELSE
        DO i=1,matrix_diag_1_row
                DO j=1,matrix_col
                        IF (counter_dia_high==jptra_dia_high)THEN
                                EXIT
                        ELSE
                                trans_dia_high = highfreq_table_dia(counter_dia_high+1)
                                matrix_diag_1(i,j)%var_name =dianm(trans_dia_high)
                                counter_dia_high=counter_dia_high + 1
                        END IF
                END DO
        END DO

        END SUBROUTINE POPULATE_MATRIX_VARS
!------------------------------------------------------------
        SUBROUTINE CLEAN_MATRIX_VARS()

        DEALLOCATE (matrix_state_1)
        DEALLOCATE (matrix_state_2)
        DEALLOCATE (matrix_diag_1)
        DEALLOCATE (matrix_diag_1)
        DEALLOCATE (matrix_diag_2d_1)
        DEALLOCATE (matrix_diag_2d_1)
        
        END SUBROUTINE CLEAN_MATRIX_VARS


!--------------------------------------------------




END MODULE

