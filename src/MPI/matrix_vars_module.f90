MODULE MATRIX_VARS
        !---------------------------------------------------------------------
        !
        !                       MODULE MATRIX_VARS
        !                     ******************
        !
        !
        !  gcoidess developing 
        !  Purpose :
        !  ---------
        !  OBTAIN MATRIX OF VARIABLES TO DUMP
        !  ---------
        !  Subroutines: 
        !       -ALLOCATE_MATRIX_VARS()
        !       -POPULATE_MATRIX_VARS()       
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
        
        !INTEGER IERROR
        !INTEGER :: nodes = 3 !declaration only to see if everything is correct
        !INTEGER :: jptra = 10!"           "
        !INTEGER :: matrix_row
        !INTEGER :: matrix_col != nodes
        !INTEGER :: i,j
        !INTEGER :: couniter
        !CHARACTER :: ctrcnm(10)
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

        INTEGER, allocatable, dimension(:) :: DA_table
        type(processor_string), allocatable :: matrix_DA(:,:)
        type(processor_string), allocatable :: PX_matrix(:)
        INTEGER :: matrix_DA_row,num_DA_vars,PX_DA

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
#ifdef ExecDA
        CALL DA_VARS()
#endif
        END SUBROUTINE POPULATE_MATRIX_VARS
!------------------------------------------------------------
        SUBROUTINE DA_VARS()

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

        !how many DA?
        INTEGER :: tmp_var_num,i,j,DA_vars_parallel
        !INTEGER, allocatable, dimension(:) :: DA_table
        !type(processor_string), allocatable :: matrix_DA(:,:)
        !type(processor_string), allocatable :: PX_matrix(:)
        INTEGER :: counter_DA
        INTEGER::trans_DA


        num_DA_vars = size(varlistDA)       
        write(*,*) 'DA vars are', num_DA_vars
        allocate(DA_table(num_DA_vars))
        
        !create DA_table
        !table of index transformtion from varlist_da and ctrcnm
        tmp_var_num = 0
        do i = 1, num_DA_vars
                do j=1, jptra      
                        write(*,*)'valistda', varlistDA(i)
                        IF (varlistDA(i).eq.trim(ctrcnm(j))) then
                                write(*,*)'var list is',varlistDA(i)
                                write(*,*)'ctrcnm is', ctrcnm(j)
                                write(*,*)'tmp_var_num is',tmp_var_num
                                tmp_var_num=tmp_var_num + 1
                                write(*,*)'tmp_var_num is',tmp_var_num
                                DA_table(tmp_var_num) = j
                        end if
                end do
        end do

        do i=1, num_DA_vars
                write(*,*) 'Da table num', i, 'is', DA_table(i)
        end do
        
        !define hard coded how many variables of varlist_DA will be printed
        !separately
        !PX=p1l+p2l+p3l+p4l -> for chl 
        PX_DA=4
        allocate(PX_matrix(PX_DA))
       
        DA_vars_parallel=num_DA_vars-PX_DA
 
        !define the DA_matrix for printing in parallel
        IF (MOD(DA_vars_parallel,nodes)==0)THEN
                matrix_DA_row = (DA_vars_parallel/nodes)
        ELSE
                matrix_DA_row = (DA_vars_parallel/nodes) + 1
        END IF
        
        allocate(matrix_DA(matrix_DA_row,matrix_col))

        DO i=1,matrix_DA_row
                DO j=1,matrix_col
                        matrix_DA(i,j)%var_name = novars
                END DO
        END DO
        
        counter_DA=0
        DO i=1,matrix_DA_row
                DO j=1,matrix_col
                        IF (counter_DA==num_DA_vars)THEN
                                EXIT
                        ELSE IF(counter_DA < PX_DA) then
                                write(*,*) 'counter da is', counter_DA
                                trans_DA = DA_table(counter_DA+1)
                                write(*,*) 'trans da is', trans_DA
                                PX_matrix(j)%var_name = ctrcnm(trans_DA)
                                counter_DA=counter_DA + 1
                        ELSE
                                trans_DA = DA_table(counter_DA+1)
                                matrix_DA(i,j)%var_name = ctrcnm(trans_DA)
                                counter_DA=counter_DA + 1                       
                        END IF
                END DO
        END DO

        END SUBROUTINE DA_VARS
!--------------------------------------------------
        SUBROUTINE CLEAN_MATRIX_VARS()

        DEALLOCATE (matrix_state_1)
        DEALLOCATE (matrix_state_2)
        DEALLOCATE (matrix_diag_1)
        DEALLOCATE (matrix_diag_1)
        DEALLOCATE (matrix_diag_2d_1)
        DEALLOCATE (matrix_diag_2d_1)
        DEALLOCATE (matrix_DA)
        END SUBROUTINE CLEAN_MATRIX_VARS


!--------------------------------------------------




END MODULE

