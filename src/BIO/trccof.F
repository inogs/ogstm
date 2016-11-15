      SUBROUTINE trccof
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trccof
!!!                     ******************
!!!
!!!  Purpose :
!!!  --------
!!!     Compute the arrays used for newtonian damping on passive tracer fields
!!!     ( if "key_trc_dmp" key activated )
!!!
CC   Method :
CC   -------
CC      Arrays defining the damping are computed for each grid
CC      point for each passive tracer (restotr)
CC      Damping depends on distance to coast, depth and latitude



       USE myalloc
       ! epascolo USE myalloc_mpp
       USE BC_mem

       IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================

#if defined key_trc_dmp
       INTEGER jn, jk, jj, ji, tra_idx


C===================================================
C===================================================

      DO jn=1,jptra
            DO ji=1,jpi
          DO jj=1,jpj
        DO jk=1,jpk
              restotr(jk,jj,ji,jn)=0.
            END DO
          END DO
        END DO
      END DO

      !print *,allocated(restotr),allocated(restocorr),allocated(resto)

      DO jn=1,jn_gib
        tra_idx=tra_matrix_gib(jn)
            DO ji=1,jpi
          DO jj=1,jpj
        DO jk=1,jpk
              !print *,jk,jj,ji,tra_idx,jn,restocorr(jn),resto(jk,jj,ji,jn)
              restotr(jk,jj,ji,tra_idx)=restocorr(jn)*resto(jk,jj,ji,jn)
            END DO
          END DO
        END DO
      END DO
#  else
C     No evaluation of restoring coefficients
#endif

      END SUBROUTINE trccof
