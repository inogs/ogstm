      SUBROUTINE SNUTELL(datestr, ISLOG, ApplyConditions)
      use myalloc
      use DA_MEM
      use filenames

      IMPLICIT NONE
      character*17 datestr
      logical ISLOG, ApplyConditions


      ! LOCAL
      character(LEN=3) Var
      character ctype
      integer i,j,k, Itype
      integer, parameter :: jpk_200=26



      Var= 'chl'; call read_corr(CORR_FILE,Var, CORR, jpk_200 ) ! CORR can be logarithmic

      if (ISLOG) then
       ! CHLnew = exp(LOG(CHLtot)+CORR)
       ! CORR_c = CHLnew - CHLtot
      else
        CORR_c = CORR
      endif


      Factor = 1.0 + CORR_c/CHLtot
      write(*,*) 'Corr Factor Calculated'
      ! CorrFactor.nc ha questi
!        float C1(depth, latitude, longitude) ; --> ci scrivo Factor dopo la correzione "Factor limitation to 10.0"
!        float P1(depth, latitude, longitude) ;
!        float P2(depth, latitude, longitude) ;
!        float P3(depth, latitude, longitude) ;
!        float P4(depth, latitude, longitude) ;
!        float chl(depth, latitude, longitude) ; --> ci scrivo Factor appena calcolato
        if(ApplyConditions) then
          do k=1,jpk_200
          do j=1,jpjglo
          do i=1,jpiglo
            if (Factor(i,j,k).gt.10.0) Factor(i,j,k)=  10.0
          enddo
          enddo
          enddo
        write(*,*) 'Factor limitation to 10.0 ... done '
       endif

        DO Itype=1,4

           write(ctype, '(i1)') Itype
           write(*,*) 'chiamata NEWCORRFACTOR', Itype
           CALL NEWCORRFACTOR(datestr,ctype, ApplyConditions)

        ENDDO




      END SUBROUTINE SNUTELL


!****************************************************************************
!****************************************************************************
!****************************************************************************

      SUBROUTINE NEWCORRFACTOR(datestr, CTYPE, ApplyConditions)
        USE myalloc
        USE DA_MEM
        USE TIME_MANAGER
        implicit none

        character CTYPE
        character*17 datestr
        ! local
        integer, parameter :: jpk_200=26
        integer i,j,k, SOMMA

        LOGICAL,    DIMENSION(jpiglo,jpjglo,jpk) :: COND2, COND3, COND4, Factor_GT1, Factor_LT0, COND
        LOGICAL,    DIMENSION(jpiglo,jpjglo,jpk) :: COND5, COND6, COND7, COND8
        INTEGER(1), DIMENSION(jpiglo,jpjglo,jpk) :: ICOND2, ICOND3, ICOND4, IFactor_LT0
        INTEGER(1), DIMENSION(jpiglo,jpjglo,jpk) :: ICOND5, ICOND6, ICOND7,ICOND8
        character*3 PVar, NVar, CVar, IVAR
        character*3 SVar
        CHARACTER(LEN=39) AVEFILE  ! DA/ave.20091231-12:00:00.P1n.nc
        logical B
        real(4) MAX_N_CHL, MAX_P_CHL, MAX_P_C, MAX_N_C
        real(4) OPT_N_C, OPT_P_C, OPT_S_C, LIM_THETA
        real(4) ONE
        CHARACTER(LEN=37) filename
        real(8) julian
        logical ApplyConditions
!        real fillValueAVE
!        integer(1) fillVAlueCOND
!        fillValueAVE  = 1.0e+20
!        fillValueCOND = 2
         julian=datestring2sec(datestr)
        
        ONE = 1.0
        MAX_N_CHL = 150. ! Derived from max chl:c=0.02 (BFMconsortium)
        MAX_P_CHL =  10.
        MAX_P_C   =  7.86e-4*2 ! values from BFMconsortium parametrs document (P.Lazzari)
        OPT_P_C   =  7.86e-4
        MAX_N_C   =  1.26e-2*2 ! values from BFMconsortium parametrs document (P.Lazzari)
        OPT_N_C   =  1.26e-2
        OPT_S_C   =  0.01 ! values from BFMconsortium parametrs document (P.Lazzari)



        LIM_THETA =  0.01

        PVar='P'//CTYPE//'p'
        NVar='P'//CTYPE//'n'
        IVar='P'//CTYPE//'l'
        CVar='P'//CTYPE//'c'
        if (CTYPE.eq.'1')then
        SVar='P'//CTYPE//'s'
        endif


         write(*,*) Pvar, ' ' , Nvar, ' ', Cvar,  ' ', Ivar
        AVEFILE = 'DA__FREQ_1/ave.'//datestr//'.P1l.nc'
        AVEFILE(34:36) = PVar ; call readNetCDF_3dvar(AVEFILE,PVar, Pfraction    )
        AVEFILE(34:36) = NVar ; call readNetCDF_3dvar(AVEFILE,NVar, Nfraction    )
        AVEFILE(34:36) = IVar ; call readNetCDF_3dvar(AVEFILE,IVar, CHLfraction  )
        AVEFILE(34:36) = CVar ; call readNetCDF_3dvar(AVEFILE,CVar, Cfraction    )
        if (CTYPE.eq.'1')then
           AVEFILE(34:36) = SVar ; call readNetCDF_3dvar(AVEFILE,SVar, Sfraction    )
        endif



        CORR  = Factor
        CORRC = Factor
        CORRN = Factor
        CORRP = Factor
        CORRS = Factor

        N_CHL = Nfraction / CHLfraction
        P_CHL = Pfraction / CHLfraction
        CHL_C = CHLfraction / Cfraction
        N_C   = Nfraction / Cfraction
        P_C   = Pfraction / Cfraction
        if (CTYPE.eq.'1')then
        S_C   = Sfraction / Cfraction
        endif

        ! Booleans 3D
        Factor_GT1 = (Factor.gt.1.)
        Factor_LT0 = (Factor.le.0.)
        COND2      = ((N_CHL.gt.MAX_N_CHL) .OR. (P_CHL.gt.MAX_P_CHL )) .AND. Factor_GT1
        COND3      = (N_C.gt.(4*MAX_N_C)) .AND. Factor_GT1
        COND4      = (P_C.gt.(4*MAX_P_C)) .AND. Factor_GT1

        ! from boolean to integers
        ICOND2      = 0
        ICOND3      = 0
        ICOND4      = 0
        IFactor_LT0 = 0

        if(ApplyConditions) then
          do k=1,jpk_200
          do j=1,jpjglo
          do i=1,jpiglo
            if (COND2     (i,j,k))    ICOND2  (i,j,k)=1
            if (COND3     (i,j,k))    ICOND3  (i,j,k)=1
            if (COND4     (i,j,k))    ICOND4  (i,j,k)=1
            if (Factor_LT0(i,j,k)) IFactor_LT0(i,j,k)=1
          enddo
          enddo
          enddo
        else
          do k=1,jpk_200
          do j=1,jpjglo
          do i=1,jpiglo
            if (Factor_LT0(i,j,k)) IFactor_LT0(i,j,k)=1
          enddo
          enddo
          enddo
          COND2 = .false.
          COND3 = .false.
          COND4 = .false.
        endif

!        if (CTYPE.eq.'1') then
!          call setFillValue_3D_byte(jpiglo,jpjglo,jpk,CHL,fillValueAVE, IFactor_LT0, fillValueCond   )
!          call  MODIFY_NC_3D_byte('conditionsForCorr2.nc','FactNeg', jpiglo,jpjglo,jpk, IFactor_LT0)
!        endif
!
!
!
!        call setFillValue_3D_byte(jpiglo,jpjglo,jpk,CHL,fillValueAVE, ICOND2, fillValueCond  )
!        call setFillValue_3D_byte(jpiglo,jpjglo,jpk,CHL,fillValueAVE, ICOND3, fillValueCond  )
!        call setFillValue_3D_byte(jpiglo,jpjglo,jpk,CHL,fillValueAVE, ICOND4, fillValueCond  )
!
!        call  MODIFY_NC_3D_byte('conditionsForCorr2.nc','cond_2'//CTYPE, jpiglo,jpjglo,jpk, ICOND2)
!        call  MODIFY_NC_3D_byte('conditionsForCorr2.nc','cond_3'//CTYPE, jpiglo,jpjglo,jpk, ICOND3)
!        call  MODIFY_NC_3D_byte('conditionsForCorr2.nc','cond_4'//CTYPE, jpiglo,jpjglo,jpk, ICOND4)

        COND = COND2 .OR. COND3 .OR. COND4
        ! Forcing CORR = 1

        SOMMA=0
        do k=1,jpk_200
        do j=1,jpjglo
        do i=1,jpiglo
           if (COND(i,j,k)) then
             CORR(i,j,k) = 1.0;
             SOMMA=SOMMA+1
           endif
           if (FACTOR_LT0(i,j,k)) then
              CORR(i,j,k) = .01
           endif
        enddo
        enddo
        enddo

        CORRC = CORR
        CORRN = CORR
        CORRP = CORR
        CORRS = CORR
        
        !Check on positive CORR and on internal quotas
        if(ApplyConditions) then
          Factor_GT1 = (CORR.gt.1.)
          COND5      = (CHL_C.lt.LIM_THETA) .AND. Factor_GT1
          COND6      = (N_C.gt.OPT_N_C) .AND. Factor_GT1
          COND7      = (P_C.gt.OPT_P_C) .AND. Factor_GT1
          if (CTYPE.eq.'1')then
          COND8      = (S_C.gt.OPT_S_C) .AND. Factor_GT1
          endif

          ICOND5      = 0
          ICOND6      = 0
          ICOND7      = 0
          if (CTYPE.eq.'1')then
          ICOND8      = 0
          endif

          do k=1,jpk_200
          do j=1,jpjglo
          do i=1,jpiglo
            if (COND5     (i,j,k))    ICOND5  (i,j,k)=1
            if (COND6     (i,j,k))    ICOND6  (i,j,k)=1
            if (COND7     (i,j,k))    ICOND7  (i,j,k)=1
            if (CTYPE.eq.'1')then
                if (COND8     (i,j,k))    ICOND8  (i,j,k)=1
            endif
          enddo
          enddo
          enddo
        
          do k=1,jpk_200
          do j=1,jpjglo
          do i=1,jpiglo
            if (COND5(i,j,k)) then
              CORRC(i,j,k) = max(ONE , (CORR(i,j,k)*CHL_C(i,j,k)/LIM_THETA) )
            endif
            if (COND6(i,j,k)) then
              CORRN(i,j,k) = max(ONE,CORRC(i,j,k)/N_C(i,j,k)*OPT_N_C)
            endif
            if (COND7(i,j,k)) then
              CORRP(i,j,k) = max(ONE,CORRC(i,j,k)/P_C(i,j,k)*OPT_P_C)
            endif
            if (CTYPE.eq.'1')then
            if (COND8(i,j,k)) then
              CORRS(i,j,k) = max(ONE,CORRC(i,j,k)/S_C(i,j,k)*OPT_S_C)
            endif
            endif
          enddo
          enddo
          enddo
        endif

        write(*,*) 'In', SOMMA , 'points Correction factor forced to avoid numerical side effects '





        filename = 'RESTARTS/RST.'//datestr//'.'//IVAR//'.nc'
        write(*,*) 'SNUTELL writes ', trim(filename), ' and all the group'

       ! qui dovrei scrivere il RST.after

        filename(32:34)=IVAR ; tottrn = real(CHLfraction* CORR ,8) ; CALL write_restartDA(filename,julian)
        filename(32:34)=NVAR ; tottrn = real(  Nfraction* CORRN,8) ; CALL write_restartDA(filename,julian)
        filename(32:34)=PVAR ; tottrn = real(  Pfraction* CORRP,8) ; CALL write_restartDA(filename,julian)
        filename(32:34)=CVAR ; tottrn = real(  Cfraction* CORRC,8) ; CALL write_restartDA(filename,julian)
        if (CTYPE.eq.'1')then
        filename(32:34)=SVAR ; tottrn = real(  Sfraction* CORRS,8) ; CALL write_restartDA(filename,julian)
        endif

        END SUBROUTINE NEWCORRFACTOR

