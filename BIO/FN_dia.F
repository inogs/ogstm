       MODULE FN_dia

       USE myalloc
       USE FN_mem
       USE myalloc_mpp

       IMPLICIT NONE

       public

!----------------------------------------------------------------------

!----------------------------------------------------------------------
      CONTAINS

   
      subroutine count_FN(kt)
       REAL(8) cou_tot 
       REAL(8) b,FN_local(jptra)
       INTEGER a,FN_winners(jptra)
       INTEGER ji,jj,jk,jn,jnn,kt

           DO jn=1,jptra
            DO jk=1,jpk
             DO jj=2,jpjm1
              DO ji=2,jpim1
                IF(tra(ji,jj,jk,jn).LT.0.) THEN
!C                write(*,*) 'tra(',ji,',',jj,',',jk,',',jn,')=',
!C     &                      tra(ji,jj,jk,jn)
                 cou_FN(ji,jj,jk) = cou_FN(ji,jj,jk) + 1.
                 FN_ranking(jn) = FN_ranking(jn) + 1.
                END IF
               END DO
              END DO
            END DO
           END DO

           IF(mod( kt , nwritetrc ).EQ.0) THEN

           cou_tot = 0.
           DO jn=1,jptra
             FN_winners(jn)=jn
             FN_local(jn)=FN_ranking(jn)
           END DO

           DO jk=1,jpk
             DO jj=1,jpj
              DO ji=1,jpi
                 cou_tot = cou_tot +cou_FN(ji,jj,jk) 
              END DO
             END DO
           END DO
          
           DO jnn=1,jptra-1
            DO jn=1,jptra-1
              IF(FN_local(jn).LT.FN_local(jn+1)) THEN

               a = FN_winners(jn)
               FN_winners(jn) = FN_winners(jn+1)
               FN_winners(jn+1) = a 

               b = FN_local(jn)
               FN_local(jn) = FN_local(jn+1)
               FN_local(jn+1) = b

!C                write(*,*) 'FN_winners(',jn,')=',FN_winners(jn)
              ENDIF 
            END DO
           END DO
           write(*,*) 'Proc N',rank,'FN_count',cou_tot     
           write(*,*) 'First FN:',ctrcnm(FN_winners(1)),'FN_count',FN_ranking(FN_winners(1))
           write(*,*) 'Second FN:',ctrcnm(FN_winners(2)),'FN_count',FN_ranking(FN_winners(2))
           write(*,*) 'Third FN:',ctrcnm(FN_winners(3)),'FN_count',FN_ranking(FN_winners(3))
           write(*,*) 'Fourth FN:',ctrcnm(FN_winners(4)),'FN_count',FN_ranking(FN_winners(4))
           write(*,*) 'Fifth FN:',ctrcnm(FN_winners(5)),'FN_count',FN_ranking(FN_winners(5))
!C           DO jn=1,jptra
!C                 write(*,*) 'FN_ranking(',jn,')=',FN_ranking(jn)
!C                 write(*,*) 'FN_winners(',jn,')=',ctrcnm(FN_winners(jn))
!C           END DO
           END IF
      
      END subroutine count_FN

      END MODULE 
