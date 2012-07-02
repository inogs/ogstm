      SUBROUTINE fluxdump(datemean, datefrom, dateend)
!     ******************

      USE DIA_mem
      USE netcdf
      USE myalloc_mpp

      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, datefrom, dateend

      INTEGER jf,js,jn
!      INTEGER jkk!,jjg,jkg
!      INTEGER flux_unit

      CHARACTER(LEN=31) flux_file
      CHARACTER(LEN=2) rankstring

      integer s, nc,INDid, junk!, counter
      integer nid, tid
      INTEGER IDS(jptra)


      if (Fsize .GT. 0) THEN

      write(rankstring,1414) rank

!      flux_file= 'flux.'//datenow//'.COR'//rankstring//'.as'
!
!      flux_unit=5028
!
!          OPEN(UNIT=flux_unit,FILE=flux_file,FORM='FORMATTED')
!
!          do js=1,7
!             do jn=1,jptra
!                do jf=1,Fsize
!                   jkk = INDFluxGlo(flx_ridxt(jf,1))
!                   write(flux_unit,*) jkk, diaflx(jf,jn,js)
!                enddo
!             enddo
!          enddo
!
!          CLOSE(flux_unit)
1414  FORMAT(I2.2)


      flux_file= 'flux.'//datemean//'.COR'//rankstring//'.nc'

      s = nf90_create(flux_file, NF90_CLOBBER, nc)

      s = nf90_put_att(nc, nf90_global, 'Time_Start'     , datefrom)
      s = nf90_put_att(nc, nf90_global, 'Time___End'     ,  dateend)


      s= nf90_def_dim(nc,'n'           , Fsize,  nid)
      s= nf90_def_dim(nc,'Type'        ,     7,  tid)

      ! ******************   var definition
      s = nf90_def_var(nc,'index' ,   nf90_int, (/nid/),   IndID)

      DO jn=1,jptra
         s = nf90_def_var(nc,ctrcnm(jn) ,   nf90_double, (/tid,nid/),   IDS(jn))
      ENDDO

      s =nf90_enddef(nc)

      !  ************************ now, put var
      s = nf90_put_var(nc, IndID,      INDflxDUMP)

      DO jn=1,jptra

          do js=1,7
               do jf=1,Fsize
                 MflxDump(js,jf)=       diaflx(jf,jn,js)
               enddo
          enddo

          s = nf90_put_var(nc,IDS(jn) ,  MflxDump)

      ENDDO



      s= nf90_close(nc)

! reset integral of flux

      diaflx = 0

      ENDIF



      end SUBROUTINE fluxdump
