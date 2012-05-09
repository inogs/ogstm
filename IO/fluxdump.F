      SUBROUTINE fluxdump(datenow)
!     ******************
      USE calendar
      USE myalloc
      USE myalloc_mpp
      USE IO_mem
      USE FN_mem
      USE TIME_MANAGER
      USE DIA_mem
#ifdef key_mpp
      USE myalloc_mpp
#endif

      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datenow

      INTEGER jf,js,jn
      INTEGER jig,jjg,jkg
      INTEGER flux_unit

      CHARACTER(LEN=34) flux_file
      CHARACTER(LEN=2) rankstring



      if (Fsize .GT. 0) THEN

      write(rankstring,1414) rank

      flux_file= 'flux.'//datenow//'.COR'//rankstring//'.ascii'

      flux_unit=5028

          OPEN(UNIT=flux_unit,FILE=flux_file,FORM='FORMATTED')


          do js=1,7
             do jn=1,jptra
                do jf=1,Fsize
                   jig = flx_ridxt(jf,2)
                   jjg = flx_ridxt(jf,3)
                   jkg = flx_ridxt(jf,4)
                   write(flux_unit,*) jig, jjg, jkg, diaflx(jf,jn,js)
                enddo
             enddo
          enddo

          CLOSE(flux_unit)

      ENDIF
1414  FORMAT(I2.2)

      end SUBROUTINE fluxdump
