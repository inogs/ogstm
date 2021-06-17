MODULE DA_VARS_module
        !---------------------------------------------------------------------
        !
        !                       MODULE MATRIX_VARS
        !                     ******************
        !
        !
        !  gcoidessa@inogs.it developing 
        !
        !  Purpose :
        !  ---------
        !  OBTAIN DA MATRIX OF VARIABLES TO DUMP
        !  ---------
        !  Subroutines: 
        !       -DA_VARS(): define matrix PX and DA_matrix that are used for
        !       dumping variables
        !       -CLEAN_DA_VARS() : deallocate main matrix  

        USE calendar
        USE myalloc
        USE IO_mem
        USE FN_mem
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module

        USE dtype_procs_string_module
        USE nodes_module
        USE DA_MEM



        IMPLICIT NONE

        PUBLIC

        INTEGER, allocatable, dimension(:) :: DA_table
        type(processor_string), allocatable :: matrix_DA(:,:)
        type(processor_string), allocatable :: PX_matrix(:)
        INTEGER :: matrix_DA_row,num_DA_vars,PX_DA

        
        CONTAINS

!-------------------------------------------------------------------------------------------------------
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
        USE MATRIX_VARS
        USE DA_MEM

        IMPLICIT NONE

        INTEGER :: tmp_var_num,i,j,DA_vars_parallel
        INTEGER :: counter_DA
        INTEGER :: trans_DA


        num_DA_vars = size(DA_varlist)
        allocate(DA_table(num_DA_vars))

        !create DA_table
        !table of index transformtion from DA_varlist and ctrcnm
        tmp_var_num = 0
        do i = 1, num_DA_vars
                do j=1, jptra
                        IF (DA_varlist(i).eq.trim(ctrcnm(j))) then
                                tmp_var_num=tmp_var_num + 1
                                DA_table(tmp_var_num) = j
                        end if
                end do
        end do


        !define hard coded how many variables of DA_varlist will be printed separately
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

        DO i=1, PX_DA
                trans_DA = DA_table(counter_DA + 1)
                PX_matrix(i)%var_name = ctrcnm(trans_DA)
                counter_DA=counter_DA + 1
        END DO

        DO i=1,matrix_DA_row
                DO j=1, matrix_col
                        IF (counter_DA==num_DA_vars)THEN
                                EXIT
                        ELSE
                                trans_DA = DA_table(counter_DA+1)
                                matrix_DA(i,j)%var_name = ctrcnm(trans_DA)
                                counter_DA=counter_DA + 1
                        END IF
                END DO
        END DO


        
        END SUBROUTINE DA_VARS
!----------------------------------------------------------------------------------------------------------------------

        SUBROUTINE CLEAN_DA_VARS()

        DEALLOCATE (matrix_DA)
        DEALLOCATE (PX_matrix)
        DEALLOCATE (DA_table)

        END SUBROUTINE CLEAN_DA_VARS

!--------------------------------------------------------------------------------------------------------------------------


END MODULE 
