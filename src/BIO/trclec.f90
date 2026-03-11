 SUBROUTINE trclec
!---------------------------------------------------------------------
!
!                       ROUTINE trclec
!                     ******************
!
!  PURPOSE :
!  ---------
!     READ and PRINT options for the passive tracer run (namelist)
!     ADDED READ options for physic tracers through namelist.phys
!   INPUT :
!   -----
!      the namelist FILE ( UNIT numnat ) :
!            &nattrc           : general
!            &nattrc_diag      : general
!      namelist phys (UNIT numphys):
!            &PHYS_FREQ
!            &PHYS>_3D and 2D


       USE myalloc
       USE BIO_mem 
#ifdef key_trc_fabm
       USE fabm
       USE ogstm_yaml_reader
#endif

       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================

      INTEGER :: i,j, ji
      INTEGER :: nmatch
#ifdef key_trc_fabm
      type(ConfigYAML) :: cfg
#endif

!----------------------------------------------------------------------
! statement functions
! ===================

!passive tracers
#ifdef key_trc_bfm

      namelist /NATTRC/           ctrcnm, ctrcun, ctrmax, ctr_hf
      namelist /NATTRC_DIAG/      dianm, diaun, diahf, diaWR
      namelist /NATTRC_DIAG_2d/   dianm_2d, diaun_2d, diahf_2d ,diaWR_2d

#elif  key_trc_fabm

      call read_config_yaml("ogstm.yaml", cfg)
      call print_config(cfg)

! Interior state variables features and dump frequency
      ctrmax(:)   = 100000.0
      ctr_hf(:)   = 0
      do i = 1, size(model_fabm%interior_state_variables)
      ctrcnm(i) = model_fabm%interior_state_variables(i)%name
      ctrcun(i) = model_fabm%interior_state_variables(i)%units
      end do
      do j = 1, size(cfg%interior_state)
          nmatch = 0
          do i = 1, size(model_fabm%interior_state_variables)
              if  ( trim(model_fabm%interior_state_variables(i)%name) .EQ. trim(cfg%interior_state(j)%name)) then
                  nmatch = nmatch + 1
                  write(*,'(a)') " - "//trim(cfg%interior_state(j)%name)
                  write(*,'(a,1x,es12.5)') "    ctrmax:", cfg%interior_state(j)%ctrmax
                  ctrmax(i)=cfg%interior_state(j)%ctrmax
                  write(*,'(a,1x,i0)')     "    ctrhf: ", cfg%interior_state(j)%ctrhf
                  ctr_hf(i)=cfg%interior_state(j)%ctrhf
                  write(*,'(a,1x,i0)')     "    relax: TO BE COMPLETED", cfg%interior_state(j)%relax
              endif
          end do
          if (nmatch == 0) then
              write(*,*) "ERROR: variable not found ..." , &
              trim(cfg%interior_state(j)%name)       
              stop 1
          else if (nmatch > 1) then
              write(*,*) "ERROR: duplicate FABM variable name: ", &
              trim(cfg%interior_state(j)%name)
              stop 2
          end if
      end do
! Interior diagnostic variables features and dump frequency
      diahf(:)    = 0
      diaWR(:)    = 0
      do i = 1, size(model_fabm%interior_diagnostic_variables)
          dianm(i) = model_fabm%interior_diagnostic_variables(i)%name
          diaun(i) = model_fabm%interior_diagnostic_variables(i)%units
      end do
      do j = 1, size(cfg%interior_diagnostic)
          nmatch = 0
          do i = 1, size(model_fabm%interior_diagnostic_variables)
              if  ( trim(model_fabm%interior_diagnostic_variables(i)%name) .EQ. trim(cfg%interior_diagnostic(j)%name)) then
                  nmatch = nmatch + 1
                  diahf(i)=cfg%interior_diagnostic(j)%diahf
                  diaWR(i)=cfg%interior_diagnostic(j)%diaWR
                  IF (diaWR(i)>0) THEN
                    write(*,'(a)') " - "//trim(cfg%interior_diagnostic(j)%name)
                    write(*,'(a,1x,i0)') "    diahf:", cfg%interior_diagnostic(j)%diahf
                    write(*,'(a,1x,i0)') "    diaWR:", cfg%interior_diagnostic(j)%diaWR
                  END IF
                      
              endif
          end do
          if (nmatch == 0) then
              write(*,*) "ERROR: variable not found ..." , &
              trim(cfg%interior_diagnostic(j)%name)       
              stop 1
          else if (nmatch > 1) then
              write(*,*) "ERROR: duplicate FABM variable name: ", &
              trim(cfg%interior_diagnostic(j)%name)
              stop 2
          end if
      end do
! Horizontal diagnostic variables features and dump frequency
      diahf_2d(:) = 0
      diaWR_2d(:) = 0
      do i = 1, size(model_fabm%horizontal_diagnostic_variables)
          dianm_2d(i) = model_fabm%horizontal_diagnostic_variables(i)%name
          diaun_2d(i) = model_fabm%horizontal_diagnostic_variables(i)%units
      end do
      do j = 1, size(cfg%horizontal_diagnostic)
          nmatch = 0
          do i = 1, size(model_fabm%horizontal_diagnostic_variables)
              if  ( trim(model_fabm%horizontal_diagnostic_variables(i)%name) .EQ. trim(cfg%horizontal_diagnostic(j)%name)) then
                  nmatch = nmatch + 1
                  diahf_2d(i)=cfg%horizontal_diagnostic(j)%diahf_2d
                  diaWR_2d(i)=cfg%horizontal_diagnostic(j)%diaWR_2d
                  if (diaWR_2d(i)>0) then
                      write(*,'(a)') " - "//trim(cfg%horizontal_diagnostic(j)%name)
                      write(*,'(a,1x,i0)') "    diahf_2d:", cfg%horizontal_diagnostic(j)%diahf_2d
                      write(*,'(a,1x,i0)') "    diaWR_2d:", cfg%horizontal_diagnostic(j)%diaWR_2d
                  end if
              endif
          end do
          if (nmatch == 0) then
              write(*,*) "ERROR: variable not found ..." , &
              trim(cfg%horizontal_diagnostic(j)%name)       
              stop 1
          else if (nmatch > 1) then
              write(*,*) "ERROR: duplicate FABM variable name: ", &
              trim(cfg%horizontal_diagnostic(j)%name)
              stop 2
          end if
      end do

!     namelist /NATTRC/           ctrmax, ctr_hf
!     namelist /NATTRC_DIAG/      diahf, diaWR
!     namelist /NATTRC_DIAG_2d/   diahf_2d ,diaWR_2d

#else

      namelist /NATTRC/           ctrcnm, ctrcun, ctrmax, ctr_hf
      namelist /NATTRC_DIAG/      dianm, diaun, diahf, diaWR
      namelist /NATTRC_DIAG_2d/   dianm_2d, diaun_2d, diahf_2d ,diaWR_2d

#endif

!physics tracers

namelist /PHYS_num/   jptra_phys, jptra_phys_2d
      namelist /PHYS_freq/  freq_ave_phys
      namelist /PHYS_3D/    physnm, physun, physWR
      namelist /PHYS_2D/    physnm_2d, physun_2d, physWR_2d

!-----------------------
      OPEN(unit=numphys, file='namelist.phys', status= 'OLD')!'FORMATTED','SEQUENTIAL')

      REWIND(numphys)
      READ(numphys,phys_num)

      CLOSE(numphys)

      allocate(physnm(jptra_phys))
      allocate(physun(jptra_phys))
      allocate(physnm_2d(jptra_phys))
      allocate(physun_2d(jptra_phys))
      allocate(physWR(jptra_phys))
      allocate(physWR_2d(jptra_phys))

!------------------------

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' ROUTINE trclec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist for passive tracers'
          WRITE(numout,*) ' ****************************'
          WRITE(numout,*) ' '
      ENDIF

!----------------------- READING PASSIVE TRACERS NAMELIST
#if !defined(key_trc_fabm)
      OPEN(unit=numnat, file='namelist.passivetrc', status= 'OLD') !'FORMATTED', 'SEQUENTIAL')

!      *****  namelist nattrc STATE VARIABLES:

      REWIND(numnat)
      READ(numnat,nattrc)

!      *****  namelist nattrc_diag DIAGNOSTIC VARIABLES:

      REWIND(numnat)
      READ(numnat,nattrc_diag)

!      *****  namelist nattrc_diag_2d DIAGNOSTIC VARIABLES surface or
!      bottom:

      REWIND(numnat)
      READ(numnat,nattrc_diag_2d)

!      *****

      CLOSE(numnat)
#endif
!---------------------- READING PHYSICS TRACERS NAMELIST

      OPEN(unit=numphys, file='namelist.phys', status= 'OLD') !'FORMATTED', 'SEQUENTIAL')

!      ***** namelist PHYS_freq:

      REWIND(numphys)
      READ(numphys,phys_freq)

!      ***** namelist PHYS_3D

      REWIND(numphys)
      READ(numphys,phys_3d)

!      ***** namelist PHYS_2D

      REWIND(numphys)
      READ(numphys,phys_2d)

!      *****

      CLOSE(numphys)

!----------------------
!--------------------- PRINTING
      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) 'nattrc'
          WRITE(numout,*) ' '
          DO ji=1,jptra
            WRITE(numout,250) 'tracer nb: ',ji,' name = ',ctrcnm(ji),&
                 ' in unit = ',ctrcun(ji), ' max ', &
                 ctrmax(ji), ' highFreq ', ctr_hf(ji)
250   FORMAT (A,I3, A8,A15,A,A,A,ES20.7,A,I1)
            WRITE(numout,*) ' '
          END DO
          WRITE(numout,*) ' '
      ENDIF

!--------------------- STATE VARIABLE 3D

      jptra_high = 0
      do ji =1,jptra
          if (ctr_hf(ji).eq.1) jptra_high = jptra_high + 1
      enddo
      allocate(highfreq_table(jptra_HIGH))
       highfreq_table = huge(highfreq_table(1))

      jptra_high = 0
      do ji =1,jptra
          if (ctr_hf(ji).eq.1) then
              jptra_high = jptra_high + 1
              highfreq_table(jptra_high) = ji
              if (lwp) WRITE(numout,*) ctrcnm(ji),&
                 ' belongs also to high freq group'
           else
               if (lwp) WRITE(numout,*) ctrcnm(ji),&
                 ' belongs only to low freq group'
           endif
      enddo

!------------------- DIAGNOSTIC VARIABLES 3D

      jptra_dia_high= 0
      do ji =1, jptra_dia
          IF (diahf(ji).eq.1) jptra_dia_high = jptra_dia_high + 1
      ENDDO

      if (lwp) write(*,*) 'High freq diagnostics number :', jptra_dia_HIGH
      allocate(highfreq_table_dia(jptra_dia_HIGH))

      jptra_dia_high = 0

      do ji =1, jptra_dia
          IF (diahf(ji).eq.1) then
            jptra_dia_high = jptra_dia_high + 1
            highfreq_table_dia(jptra_dia_high) = ji
            if (lwp) WRITE(numout,*) dianm(ji),&
               ' belongs also to high freq group'
          ELSE
            if (lwp) WRITE(numout,*) dianm(ji),&
               ' belongs only to low freq group'
          ENDIF
      enddo

!---------------- DIAGNOSTIC VARIABLE 2D
        jptra_dia2d_high= 0
      do ji =1, jptra_dia_2d
          IF (diahf_2d(ji).eq.1) jptra_dia2d_high = jptra_dia2d_high + 1
      ENDDO

      if (lwp) write(*,*) 'High freq diagnostics number 2d:', jptra_dia2d_HIGH
      allocate(highfreq_table_dia2d(jptra_dia2d_HIGH))

      jptra_dia2d_high = 0

      do ji =1, jptra_dia_2d
          IF (diahf_2d(ji).eq.1) then
            jptra_dia2d_high = jptra_dia2d_high + 1
            highfreq_table_dia2d(jptra_dia2d_high) = ji
            if (lwp) WRITE(numout,*) dianm_2d(ji),&
               ' belongs also to high freq group'
          ELSE
            if (lwp) WRITE(numout,*) dianm_2d(ji),&
               ' belongs only to low freq group'
          ENDIF
      enddo

!------------------------------------
      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** number of passive tracer jptra = ',jptra
          WRITE(numout,*) ' '
      ENDIF

      END SUBROUTINE trclec

