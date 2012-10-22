      SUBROUTINE trcsbc

! Computes surface boundary conditions on passive tracers

      USE myalloc
      USE stafun

      IMPLICIT NONE

      INTEGER  :: ji,jj,jn
      REAL(8)  :: ztra,zse3t

!      zsrau = 1. / rau0

! Conc/dilution process

      DO jn = 1, jptra
         DO jj = 1, jpj
            DO ji = 1, jpi

                  zse3t = 1. / fse3t(ji,jj,1)

                  ztra = 1./ rhopn(ji,jj,1) * zse3t * tmask(ji,jj,1) * emp(ji,jj) * trn(ji,jj,1,jn) ! original emps(ji,jj)

                  if (jn .GT. 49) then
                     tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + ztra
                     tra_pp(ji,jj,jn-49,jptra_dia) = ztra ! save evaporation rate
                  endif
            END DO
         END DO
      END DO

! Contribution of surface elevation virtual flux

      END SUBROUTINE trcsbc
