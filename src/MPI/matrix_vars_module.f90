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
        type(processor_string), allocatable :: matrix_phys_2d_1(:,:)
        type(processor_string), allocatable :: matrix_phys_2d_2(:,:)
        type(processor_string), allocatable :: matrix_phys_1(:,:)
        type(processor_string), allocatable :: matrix_phys_2(:,:)
        INTEGER :: matrix_state_1_row, matrix_state_2_row
        INTEGER :: matrix_diag_2d_1_row, matrix_diag_2d_2_row
        INTEGER :: matrix_diag_1_row, matrix_diag_2_row
        INTEGER :: matrix_phys_2d_1_row, matrix_phys_2d_2_row
        INTEGER :: matrix_phys_1_row, matrix_phys_2_row
        INTEGER :: matrix_col != nodes

        INTEGER :: jptra_dia_2d_wri, jptra_dia_2d_high_wri
        INTEGER :: jptra_dia_wri, jptra_dia_high_wri
        INTEGER :: low_freq_jptra_dia, high_freq_jptra_dia ! counters

        INTEGER :: jptra_phys_2d_high_wri, jptra_phys_2d_wri, jptra_phys_high_wri, jptra_phys_wri
        
        INTEGER, allocatable :: highfreq_table_dia_wri(:)
        INTEGER, allocatable :: highfreq_table_dia_2d_wri(:)
        INTEGER, allocatable :: lowfreq_table_dia_wri(:)
        INTEGER, allocatable :: lowfreq_table_dia_2d_wri(:)
        INTEGER, allocatable :: highfreq_table_phys_wri(:)
        INTEGER, allocatable :: highfreq_table_phys_2d_wri(:)
        INTEGER, allocatable :: lowfreq_table_phys_wri(:)
        INTEGER, allocatable :: lowfreq_table_phys_2d_wri(:)

        CONTAINS

!------------------------------------------------------------
        SUBROUTINE ALLOCATE_MATRIX_VARS()
        
        ALLOCATE (matrix_state_1(matrix_state_1_row, matrix_col))
        ALLOCATE (matrix_state_2(matrix_state_2_row, matrix_col))
        ALLOCATE (matrix_diag_2d_1(matrix_diag_2d_1_row, matrix_col))
        ALLOCATE (matrix_diag_2d_2(matrix_diag_2d_2_row, matrix_col))
        ALLOCATE (matrix_diag_1(matrix_diag_1_row, matrix_col))
        ALLOCATE (matrix_diag_2(matrix_diag_2_row, matrix_col))
        ALLOCATE (matrix_phys_2d_1(matrix_phys_2d_1_row, matrix_col))
        ALLOCATE (matrix_phys_2d_2(matrix_phys_2d_2_row, matrix_col))
        ALLOCATE (matrix_phys_1(matrix_phys_1_row, matrix_col))
        ALLOCATE (matrix_phys_2(matrix_phys_2_row, matrix_col))
        END SUBROUTINE ALLOCATE_MATRIX_VARS
!------------------------------------------------------------

!------------------------------------------------------------
        SUBROUTINE POPULATE_MATRIX_VARS()
        INTEGER :: i,j
        INTEGER :: counter,counter_high,trans_high,counter_dia2d,trans_dia2d_low,counter_dia2d_high
        INTEGER :: trans_dia2d_high,trans_dia_low,counter_dia,counter_dia_high,trans_dia_high
        INTEGER :: counter_phys2d,trans_phys2d_low,counter_phys2d_high,trans_phys2d_high,trans_phys_low
        INTEGER :: counter_phys,counter_phys_high,trans_phys_high

        counter = 0
        counter_high = 0
        counter_dia2d = 0
        counter_dia2d_high = 0
        counter_dia = 0
        counter_dia_high = 0
        counter_phys2d = 0
        counter_phys2d_high = 0
        counter_phys = 0
        counter_phys_high = 0

        ! define how many rows will be in the procs matrix
        
        matrix_col = nodes
        novars="novars_input"
        
        CALL DIA_MATRIX_VARS()


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
        IF (MOD(jptra_dia_2d_wri,nodes)==0)THEN
                matrix_diag_2d_2_row = (jptra_dia_2d_wri/nodes)
        ELSE
                matrix_diag_2d_2_row = (jptra_dia_2d_wri/nodes) + 1
        END IF

        !ELSE

        IF (MOD(jptra_dia_2d_high_wri,nodes)==0)THEN
                matrix_diag_2d_1_row = (jptra_dia_2d_high_wri/nodes)
        ELSE
                matrix_diag_2d_1_row = (jptra_dia_2d_high_wri/nodes) + 1
        END IF

        IF (MOD(jptra_dia_wri,nodes)==0)THEN
                matrix_diag_2_row = (jptra_dia_wri/nodes)
        ELSE
                matrix_diag_2_row = (jptra_dia_wri/nodes) + 1
        END IF

        !ELSE

        IF (MOD(jptra_dia_high_wri,nodes)==0)THEN
                matrix_diag_1_row = (jptra_dia_high_wri/nodes)
        ELSE
                matrix_diag_1_row = (jptra_dia_high_wri/nodes) + 1
        END IF

        IF (MOD(jptra_phys_2d_wri,nodes)==0)THEN
                matrix_phys_2d_2_row = (jptra_phys_2d_wri/nodes)
        ELSE
                matrix_phys_2d_2_row = (jptra_phys_2d_wri/nodes) + 1
        END IF

        !ELSE

        IF (MOD(jptra_phys_2d_high_wri,nodes)==0)THEN
                matrix_phys_2d_1_row = (jptra_phys_2d_high_wri/nodes)
        ELSE
                matrix_phys_2d_1_row = (jptra_phys_2d_high_wri/nodes) + 1
        END IF

        IF (MOD(jptra_phys_wri,nodes)==0)THEN
                matrix_phys_2_row = (jptra_phys_wri/nodes)
        ELSE
                matrix_phys_2_row = (jptra_phys_wri/nodes) + 1
        END IF

        !ELSE

        IF (MOD(jptra_phys_high_wri,nodes)==0)THEN
                matrix_phys_1_row = (jptra_phys_high_wri/nodes)
        ELSE
                matrix_phys_1_row = (jptra_phys_high_wri/nodes) + 1
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

        DO i=1,matrix_phys_2d_2_row
                DO j=1,matrix_col
                        matrix_phys_2d_2(i,j)%var_name = novars
                END DO
        END DO
        !ELSE
        DO i=1,matrix_phys_2d_1_row
                DO j=1,matrix_col
                        matrix_phys_2d_1(i,j)%var_name = novars
                END DO
        END DO

        DO i=1,matrix_phys_2_row
                DO j=1,matrix_col
                        matrix_phys_2(i,j)%var_name = novars
                END DO
        END DO
        !ELSE
        DO i=1,matrix_phys_1_row
                DO j=1,matrix_col
                        matrix_phys_1(i,j)%var_name = novars
                END DO
        END DO

        !OK
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
        !diagnostic 2d low freq
        DO i=1,matrix_diag_2d_2_row
                DO j=1,matrix_col
                        IF (counter_dia2d ==jptra_dia_2d_wri)THEN
                                EXIT
                        ELSE
                                trans_dia2d_low = lowfreq_table_dia_2d_wri(counter_dia2d + 1)
                                matrix_diag_2d_2(i,j)%var_name = dianm_2d(trans_dia2d_low)
                                counter_dia2d=counter_dia2d + 1
                        END IF
                END DO
        END DO

        !OK        !ELSE

        DO i=1,matrix_diag_2d_1_row
                DO j=1,matrix_col
                        IF (counter_dia2d_high==jptra_dia_2d_high_wri)THEN
                                EXIT
                        ELSE
                                trans_dia2d_high = highfreq_table_dia_2d_wri(counter_dia2d_high+1)                     
                                matrix_diag_2d_1(i,j)%var_name =dianm_2d(trans_dia2d_high)
                                counter_dia2d_high=counter_dia2d_high + 1
                        END IF
                END DO
        END DO

!OK
        DO i=1,matrix_diag_2_row
                DO j=1,matrix_col
                        IF (counter_dia==jptra_dia_wri)THEN
                                EXIT
                        ELSE
                                trans_dia_low = lowfreq_table_dia_wri(counter_dia +1)
                                matrix_diag_2(i,j)%var_name = dianm(trans_dia_low)
                                counter_dia=counter_dia + 1
                        END IF
                END DO
        END DO

        DO i=1,matrix_diag_1_row
                DO j=1,matrix_col
                        IF (counter_dia_high==jptra_dia_high_wri)THEN
                                EXIT
                        ELSE
                                trans_dia_high = highfreq_table_dia_wri(counter_dia_high+1)
                                matrix_diag_1(i,j)%var_name =dianm(trans_dia_high)
                                counter_dia_high=counter_dia_high + 1
                        END IF
                END DO
        END DO
!--------------------------------------- PHYS
        if(freq_ave_phys .eq. 2) then  
        DO i=1,matrix_phys_2d_2_row
                DO j=1,matrix_col
                        IF (counter_phys2d ==jptra_phys_2d_wri)THEN
                                EXIT
                        ELSE
                                trans_phys2d_low = lowfreq_table_phys_2d_wri(counter_phys2d + 1)
                                matrix_phys_2d_2(i,j)%var_name = physnm_2d(trans_phys2d_low)
                                counter_phys2d=counter_phys2d + 1
                        END IF
                END DO
        END DO

        DO i=1,matrix_phys_2_row
                DO j=1,matrix_col
                        IF (counter_phys==jptra_phys_wri)THEN
                                EXIT
                        ELSE
                                trans_phys_low = lowfreq_table_phys_wri(counter_phys +1)
                                matrix_phys_2(i,j)%var_name = physnm(trans_phys_low)
                                counter_phys=counter_phys + 1
                        END IF
                END DO
        END DO

        else

        !OK        !ELSE

        DO i=1,matrix_phys_2d_1_row
                DO j=1,matrix_col
                        IF (counter_phys2d_high==jptra_phys_2d_high_wri)THEN
                                !write(*,*) 'exit'
                                EXIT
                        ELSE
                                !write(*,*) '368',i ,j,highfreq_table_phys_2d_wri(counter_phys2d_high+1)
                                trans_phys2d_high = highfreq_table_phys_2d_wri(counter_phys2d_high+1)
                                !write(*,*) '368', trans_phys2d_high,highfreq_table_phys_2d_wri(counter_phys2d_high+1)
                                matrix_phys_2d_1(i,j)%var_name=physnm_2d(trans_phys2d_high)
                                counter_phys2d_high=counter_phys2d_high + 1
                        END IF
                END DO
        END DO

!OK
        DO i=1,matrix_phys_1_row
                DO j=1,matrix_col
                        IF (counter_phys_high==jptra_phys_high_wri)THEN
                                EXIT
                        ELSE
                                trans_phys_high = highfreq_table_phys_wri(counter_phys_high+1)
                                !write(*,*) '384', highfreq_table_phys_wri(counter_phys_high+1)
                                matrix_phys_1(i,j)%var_name=physnm(trans_phys_high)
                                counter_phys_high=counter_phys_high + 1
                        END IF
                END DO
        END DO
        
        end if 

        END SUBROUTINE POPULATE_MATRIX_VARS
!----------------------------------------------------------
        SUBROUTINE DIA_MATRIX_VARS()

        INTEGER :: i
        
        !high freq dia 3d
        jptra_dia_high_wri = 0
        DO i =1, jptra_dia
                IF (diahf(i).eq.1 .and. diaWR(i) == 1) then
                        jptra_dia_high_wri = jptra_dia_high_wri + 1
                END IF
        ENDDO
        ALLOCATE (highfreq_table_dia_wri(jptra_dia_high_wri))
        !high freq phys 3d
        jptra_phys_high_wri = 0
        DO i=1,jptra_phys
                if (freq_ave_phys.eq.1 .and. physWR(i) == 1) then
                        jptra_phys_high_wri = jptra_phys_high_wri + 1
                end if
        ENDDO

        ALLOCATE (highfreq_table_phys_wri(jptra_phys_high_wri))

        jptra_dia_high_wri = 0
        DO i =1, jptra_dia
                IF (diahf(i).eq.1 .and. diaWR(i) == 1) then
                        jptra_dia_high_wri = jptra_dia_high_wri + 1
                        highfreq_table_dia_wri(jptra_dia_high_wri) = i
                END IF
        ENDDO

        jptra_phys_high_wri = 0
        DO i =1, jptra_phys
                IF (freq_ave_phys.eq.1 .and. physWR(i) == 1) then
                        jptra_phys_high_wri = jptra_phys_high_wri + 1
                        highfreq_table_phys_wri(jptra_phys_high_wri) = i
                END IF
        ENDDO

        !high freq dia 2d
        jptra_dia_2d_high_wri = 0
        DO i =1, jptra_dia_2d
                IF (diahf_2d(i).eq.1 .and. diaWR_2d(i) == 1) then
                        jptra_dia_2d_high_wri = jptra_dia_2d_high_wri + 1
                END IF
        ENDDO
        
        ALLOCATE (highfreq_table_dia_2d_wri(jptra_dia_2d_high_wri))

        jptra_phys_2d_high_wri = 0
        !write(*,*) '443', jptra_phys_2d, jptra_phys_2d_high_wri
        DO i =1, jptra_phys_2d
                IF (freq_ave_phys.eq.1 .and. physWR_2d(i) == 1) then
                        jptra_phys_2d_high_wri = jptra_phys_2d_high_wri + 1
                        !write(*,*) '447', jptra_phys_2d, jptra_phys_2d_high_wri
                END IF
        ENDDO


        ALLOCATE (highfreq_table_phys_2d_wri(jptra_phys_2d_high_wri))

        jptra_dia_2d_high_wri = 0
        DO i =1, jptra_dia_2d
                IF (diahf_2d(i).eq.1 .and. diaWR_2d(i) == 1) then
                        jptra_dia_2d_high_wri = jptra_dia_2d_high_wri + 1
                        highfreq_table_dia_2d_wri(jptra_dia_2d_high_wri) = i
                END IF
        ENDDO

        jptra_phys_2d_high_wri = 0
        DO i =1, jptra_phys_2d
                IF (freq_ave_phys.eq.1 .and. physWR_2d(i) == 1) then
                        !write(*,*) '465', jptra_phys_2d, jptra_phys_2d_high_wri
                        jptra_phys_2d_high_wri = jptra_phys_2d_high_wri + 1
                        !write(*,*) '467', jptra_phys_2d, jptra_phys_2d_high_wri
                        highfreq_table_phys_2d_wri(jptra_phys_2d_high_wri) = i
                        !write(*,*) '469',i, jptra_phys_2d, jptra_phys_2d_high_wri,highfreq_table_phys_2d_wri(jptra_phys_2d_high_wri)
                END IF
        ENDDO

        !low freq dia 3d
        
        jptra_dia_wri = 0
        DO i =1, jptra_dia
                IF (diaWR(i) == 1) then
                        jptra_dia_wri = jptra_dia_wri + 1
                END IF
        ENDDO

        ALLOCATE (lowfreq_table_dia_wri(jptra_dia_wri))

        jptra_phys_wri = 0
        DO i =1, jptra_phys
                IF (freq_ave_phys.eq.2 .and. physWR(i) == 1) then
                        jptra_phys_wri = jptra_phys_wri + 1
                END IF
        ENDDO


        ALLOCATE (lowfreq_table_phys_wri(jptra_phys_wri))

        jptra_dia_wri = 0
        DO i =1, jptra_dia
                IF (diaWR(i) == 1) then
                        jptra_dia_wri = jptra_dia_wri + 1
                        lowfreq_table_dia_wri(jptra_dia_wri) = i
                END IF
        ENDDO

        jptra_phys_wri = 0
        DO i =1, jptra_phys
                IF (freq_ave_phys.eq.2 .and. physWR(i) == 1) then
                        jptra_phys_wri = jptra_phys_wri + 1
                        lowfreq_table_phys_wri(jptra_phys_wri) = i
                END IF
        ENDDO


        !lowfreq dia 2d

        jptra_dia_2d_wri = 0
        DO i =1, jptra_dia_2d
                IF (diaWR_2d(i) == 1) then
                        jptra_dia_2d_wri = jptra_dia_2d_wri + 1
                END IF
        ENDDO
        
        ALLOCATE (lowfreq_table_dia_2d_wri(jptra_dia_2d_wri))

        jptra_phys_2d_wri = 0
        DO i =1, jptra_phys_2d
                IF (freq_ave_phys.eq.2 .and. physWR_2d(i) == 1) then
                        jptra_phys_2d_wri = jptra_phys_2d_wri + 1
                END IF
        ENDDO


        ALLOCATE (lowfreq_table_phys_2d_wri(jptra_phys_2d_wri))

        jptra_dia_2d_wri = 0
        DO i =1, jptra_dia_2d
                IF (diaWR_2d(i) == 1) then
                        jptra_dia_2d_wri = jptra_dia_2d_wri + 1
                        lowfreq_table_dia_2d_wri(jptra_dia_2d_wri) = i
                END IF
        ENDDO

        jptra_phys_2d_wri = 0
        DO i =1, jptra_phys_2d
                IF (freq_ave_phys.eq.2 .and. physWR_2d(i) == 1) then
                        jptra_phys_2d_wri = jptra_phys_2d_wri + 1
                        lowfreq_table_phys_2d_wri(jptra_phys_2d_wri) = i
                END IF
        ENDDO



        END SUBROUTINE DIA_MATRIX_VARS
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

