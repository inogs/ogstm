      SUBROUTINE domrea
!---------------------------------------------------------------------

!                       ROUTINE DOMREA
!                     ******************

!  PURPOSE :
!  ---------
!       Reads files:
!                   meshmask.nc
!                   bounmask.nc
!                   BC/ATM_yyyymmdd-HH:MM:SS.nc
!                   BC/GIB_yyyymmdd-HH:MM:SS.nc
!                   BC/TIN_yyyymmdd-HH:MM:SS.nc


! parameters and commons
! ======================
      USE calendar
      USE myalloc
      USE myalloc_mpp
      USE BC_mem
      USE DIA_mem
      USE TIME_MANAGER

      IMPLICIT NONE

! local declarations
! ==================

      INTEGER ii, jj, kk, jn, iinew, jjnew, iiend, jjend
      INTEGER ierr

      CHARACTER(LEN=11) maskfile
      CHARACTER(LEN=50) filename
      CHARACTER(LEN=3), DIMENSION(7) :: var_nc
      CHARACTER(LEN=5) nomevar01
      LOGICAL B
      
! -------------------
! To read only one meshmask.nc

      maskfile        = 'meshmask.nc'
      filename='BC/TI_yyyy0215-00:00:00.nc' ! 26 chars


      iiend = MIN(jpi+nimpp-1, jpiglo)
      jjend = MIN(jpj+njmpp-1, jpjglo)
      do kk=1, jpk
         do jj=njmpp, jjend
            do ii=nimpp, iiend
               iinew = ii - nimpp + 1
               jjnew = jj - njmpp + 1
               idxt2glo(iinew, jjnew, kk,1)=ii !
               idxt2glo(iinew, jjnew, kk,2)=jj ! matrix to go from local to global
               idxt2glo(iinew, jjnew, kk,3)=kk !
            enddo
         enddo
      enddo

! 1. Horzontal grid-point position
! --------------------------------
      call readnc_global_double_2d(maskfile,'glamt', totglamt)
!     call readnc_slice_double_2d(maskfile,'glamt', glamt)
      call readnc_slice_double_2d(maskfile,'glamu', glamu)
      call readnc_slice_double_2d(maskfile,'glamv', glamv)
      call readnc_slice_double_2d(maskfile,'glamf', glamf)

      call readnc_global_double_2d(maskfile,'gphit', totgphit)
      call readnc_slice_double_2d(maskfile,'gphit', gphit)
      call readnc_slice_double_2d(maskfile,'gphiu', gphiu)
      call readnc_slice_double_2d(maskfile,'gphiv', gphiv)
      call readnc_slice_double_2d(maskfile,'gphif', gphif)


! 2. Horizontal scale factors
! ---------------------------
      call readnc_slice_double_2d(maskfile,'e1t', e1t)
      call readnc_slice_double_2d(maskfile,'e1u', e1u)
      call readnc_slice_double_2d(maskfile,'e1v', e1v)
      call readnc_slice_double_2d(maskfile,'e1f', e1f)

      call readnc_slice_double_2d(maskfile,'e2t', e2t)
      call readnc_slice_double_2d(maskfile,'e2u', e2u)
      call readnc_slice_double_2d(maskfile,'e2v', e2v)
      call readnc_slice_double_2d(maskfile,'e2f', e2f)



! 3. masks
! --------

      CALL readnc_slice_double (maskfile,'umask', umask )
      CALL readnc_slice_double (maskfile,'vmask', vmask )
      CALL readnc_slice_double (maskfile,'fmask', fmask )
      CALL readnc_slice_double (maskfile,'tmask', tmask )
!      CALL readnc_global_double(maskfile,'tmask', tmaskglo)


!      Initialization of mbathy
      mbathy(:,:) = 0
      NWATERPOINTS=0
      do jj=1, jpj
       do ii=1, jpi          
        do kk=1, jpk
         if (tmask(ii, jj, kk).NE.0.) then
            mbathy(ii,jj) = mbathy(ii,jj) +1
            NWATERPOINTS = NWATERPOINTS +1
         endif
        enddo
       enddo
      enddo

      CALL readnc_slice_double_2d(maskfile,'ff', ff )


! 4. depth and vertical scale factors
! -----------------------------------

       CALL readmask_double_1d(maskfile,'gdept',gdept)
       CALL readmask_double_1d(maskfile,'gdepw',gdepw)
       CALL readmask_double_1d(maskfile,'e3t',    e3t)
       CALL readmask_double_1d(maskfile,'e3w',    e3w)


!       Restoration Mask ****************


      var_nc(1) = 'O2o'
      var_nc(2) = 'N1p'
      var_nc(3) = 'N3n'
      var_nc(4) = 'N5s'
      var_nc(5) = 'O3c'
      var_nc(6) = 'O3h'
      var_nc(7) = 'N6r'

      IF (NWATERPOINTS.GT.0) THEN
      do jn=1,jn_gib

         nomevar01='re'//var_nc(jn)
         call readnc_slice_double('bounmask.nc',nomevar01,resto(:,:,:,jn));

      enddo
      ELSE
        resto=0.0
      ENDIF
      call readnc_slice_int   ('bounmask.nc','index',idxt);



! ************************************ BFM points re-indexing *******

      NBFMPOINTS = BFM_count()
      call myalloc_BFM()
      B=BFM_Indexing();

! *********************************   Gibraltar area
      filename  ='BC/GIB_'//TC_GIB%TimeStrings(1)//'.nc'


       if (lwp) write(*,*) 'domrea->filename: ', filename, '    '

      CALL readnc_int_1d(filename, 'gib_idxt_N1p', Gsizeglo, gib_idxtglo)

      if (lwp) write(*,*) 'domrea->readnc_int_1d  finita'
      if (lwp) write(*,*) 'domrea->Gsizeglo', Gsizeglo

      Gsize = COUNT_InSubDomain(Gsizeglo, gib_idxtglo)

      write(*,*) 'domrea->Gsize   : ', Gsize, 'rank=', rank


      if (Gsize.NE.0) then
          if (lwp) write(*,*) 'domrea-> lancio alloc_DTATRC_local_gib'
          call alloc_DTATRC_local_gib

          B=GIBRe_indexing()

      endif

! ********************************  Rivers ******
      filename       ='BC/TIN_'//TC_TIN%TimeStrings(1)//'.nc';
      CALL readnc_int_1d(filename, 'riv_idxt', Rsizeglo, riv_idxtglo)
      Rsize = COUNT_InSubDomain(Rsizeglo,riv_idxtglo)

      if (Rsize.NE. 0) then
          call alloc_DTATRC_local_riv

          B=RIVRe_Indexing()

      endif


       if(lwp) write(*,*) 'RIV finiti'

! ******************************************* Atmospherical inputs
      filename       = 'BC/ATM_'//TC_ATM%TimeStrings(1)//'.nc'
      CALL readnc_int_1d(filename, 'atm_idxt',Asizeglo,atm_idxtglo)
      Asize = COUNT_InSubDomain(Asizeglo,atm_idxtglo)

      if (Asize.NE. 0) then
         call alloc_DTATRC_local_atm
         write(*,*) 'domrea->ATMRE_Indexing ATM iniziata, rank=', rank
         B=ATMRe_Indexing()
         write(*,*) 'domrea->ATMRE_Indexing ATM finita, rank=', rank
      endif



! **************************************** FLUXES
      if (FSizeGlo.GT.0) then ! file could not exist
         CALL readnc_int_1d('Fluxes.nc','index',FSizeGlo,INDFluxGlo)
         Fsize = COUNT_InSubDomain(FSizeGlo, INDFluxGlo)
      else
         Fsize=0
      endif

      call MPI_ALLREDUCE(Fsize, FsizeMax, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)

      write(*,*) 'rank=', rank, ' Fsize = ' , Fsize,' FsizeMax = ' , FsizeMax

      if (Fsize.NE.0) then
         call alloc_DIA_local_flx
         if (rank ==0) call alloc_DIA_GLOBAL_flx()
         B=FLXRe_Indexing()
         write(*,*) 'domrea->FLXRE_Indexing finita, rank=', rank

      endif

      call alloc_DIA_MPI_flx()

      write(*,*) 'DOMREA finita, rank = ', rank


      CONTAINS

! *****************************************************************
!     FUNCTION COUNT_InSubDomain
!     RETURNS the number of points of a specific boundary condition
!     in the subdomain of the current processor
! *****************************************************************
      INTEGER FUNCTION COUNT_InSubDomain(sizeGLO,idxtGLOBAL)
          USE modul_param , ONLY: jpi, jpj, jpk
          USE myalloc     , ONLY: idxt
          USE myalloc_mpp , ONLY: rank

          IMPLICIT NONE
          INTEGER, INTENT(IN) :: sizeGLO
          INTEGER, INTENT(IN) :: idxtGLOBAL(sizeGLO)

          ! local
          INTEGER ii,jj,kk,jv
          INTEGER counter,junk

           counter = 0
           do kk =1, jpk
            do jj =1, jpj
             do ii =1, jpi
                junk = idxt(ii,jj,kk)
                do jv =1, sizeGLO
                  if (junk.EQ.idxtGLOBAL(jv))  counter = counter + 1
                enddo
             enddo
            enddo
           enddo

          COUNT_InSubDomain = counter

      END FUNCTION COUNT_InSubDomain

! *************************************************************************
      LOGICAL FUNCTION RIVRE_Indexing()
          IMPLICIT NONE
          ! local
          INTEGER ii,jj,kk,jv
          INTEGER counter,junk


          counter=0
           do kk =1, jpk
            do jj =1, jpj
             do ii =1, jpi
                junk = idxt(ii,jj,kk)
                do jv =1, RsizeGLO
                   if ( junk.EQ.riv_idxtglo(jv) )  then
                      counter = counter + 1
                      riv_ridxt(1,counter) = jv
                      riv_ridxt(2,counter) = ii
                      riv_ridxt(3,counter) = jj
                      riv_ridxt(4,counter) = kk
                   endif
              enddo
             enddo
            enddo
           enddo


      RIVRE_Indexing = .true.
      END FUNCTION RIVRE_Indexing
! *************************************************************************
      LOGICAL FUNCTION GIBRE_Indexing()
          IMPLICIT NONE

          ! local
          INTEGER ii,jj,kk,jv
          INTEGER counter,junk


          counter=0
           do kk =1, jpk
            do jj =1, jpj
             do ii =1, jpi
                junk = idxt(ii,jj,kk)
                do jv =1, Gsizeglo
                   if ( junk.EQ.gib_idxtglo(jv) )  then
                      counter = counter + 1
                      gib_ridxt(1,counter) = jv
                      gib_ridxt(2,counter) = ii
                      gib_ridxt(3,counter) = jj
                      gib_ridxt(4,counter) = kk
                   endif
              enddo
             enddo
            enddo
           enddo


      GIBRE_Indexing= .true.
      END FUNCTION GIBRE_Indexing



! *************************************************************************

          LOGICAL FUNCTION  ATMRE_Indexing()

          IMPLICIT NONE

          ! local
          INTEGER ii,jj,kk,jv
          INTEGER counter,junk


          counter=0
           do kk =1, jpk
            do jj =1, jpj
             do ii =1, jpi
                junk = idxt(ii,jj,kk)
                do jv =1, AsizeGLO
                   if ( junk.EQ.atm_idxtglo(jv) )  then
                      counter = counter + 1
                      atm_ridxt(1,counter) = jv
                      atm_ridxt(2,counter) = ii
                      atm_ridxt(3,counter) = jj
                      atm_ridxt(4,counter) = kk
                   endif
              enddo
             enddo
            enddo
           enddo
        ATMRE_Indexing =.true.
      END FUNCTION ATMRE_Indexing

! *************************************************************************
      LOGICAL FUNCTION  FLXRE_Indexing()

          IMPLICIT NONE

          ! local
          INTEGER ii,jj,kk,jv
          INTEGER counter,junk


          counter=0
           do kk =1, jpk
            do jj =1, jpj
             do ii =1, jpi
                junk = idxt(ii,jj,kk)
                do jv =1, FsizeGLO
                   if ( junk.EQ.INDFluxGlo(jv) )  then
                      counter = counter + 1
                      flx_ridxt(counter,1) = jv
                      flx_ridxt(counter,2) = ii
                      flx_ridxt(counter,3) = jj
                      flx_ridxt(counter,4) = kk
                   endif
              enddo
             enddo
            enddo
           enddo


      do ii=1,Fsize
         INDflxDUMP(ii) = INDFluxGlo(flx_ridxt(ii,1))
      enddo

        FLXRE_Indexing =.true.
      END FUNCTION FLXRE_Indexing



          LOGICAL FUNCTION  BFM_Indexing()

          IMPLICIT NONE

          ! local
          INTEGER ii,jj,kk
          INTEGER counter


          counter=0
          if (atlantic_bfm) then

               do kk =1, jpkb-1
                do jj =1, jpj-1
                 do ii =2, jpi-1

                   if (tmask(ii,jj,kk).EQ.1.0 ) then
                      counter = counter + 1
                      BFMpoints(1,counter) = ii
                      BFMpoints(2,counter) = jj
                      BFMpoints(3,counter) = kk

                   endif

                 enddo
                enddo
               enddo
           else
           ! NO ACTIVATION IN ATLANTIC BUFFER
               do kk =1, jpkb-1
                do jj =1, jpj-1
                 do ii =2, jpi-1

                   if ( (tmask(ii,jj,kk).EQ.1.0 ) .and. (resto(ii,jj,kk,1).eq.0.0 )) then
                      counter = counter + 1
                      BFMpoints(1,counter) = ii
                      BFMpoints(2,counter) = jj
                      BFMpoints(3,counter) = kk

                   endif

                 enddo
                enddo
               enddo
           endif


        BFM_Indexing =.true.
      END FUNCTION BFM_Indexing



! ***************************************************************
          INTEGER FUNCTION  BFM_count()
          USE myalloc, only:  NBFMPOINTS_SUP
          IMPLICIT NONE

          ! local
          INTEGER ii,jj,kk
          INTEGER :: counter = 0


       if (atlantic_bfm) then
           do kk =1, jpkb-1
            if (kk.eq.2)  NBFMPOINTS_SUP = counter
            do jj =1, jpj-1
             do ii =2, jpi-1

               if (tmask(ii,jj,kk).EQ.1.0 ) counter = counter + 1

             enddo
            enddo
           enddo
        else
           ! NO ACTIVATION IN ATLANTIC BUFFER
           do kk =1, jpkb-1
            if (kk.eq.2)  NBFMPOINTS_SUP = counter
            do jj =1, jpj-1
             do ii =2, jpi-1

               if ( (tmask(ii,jj,kk).EQ.1.0 ) .and. (resto(ii,jj,kk,1).eq.0.0) )  counter = counter + 1

             enddo
            enddo
           enddo

        endif
        BFM_count =counter
      END FUNCTION BFM_count





      END SUBROUTINE domrea
