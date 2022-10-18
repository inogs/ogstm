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
      USE BC_mem
      USE DIA_mem
      USE TIME_MANAGER
      USE mpi
      USE MPI_GATHER_INFO
      USE nodes_module
      USE MATRIX_VARS

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      IMPLICIT NONE

! local declarations
! ==================

      INTEGER kk,jj,ii, jn, iinew, jjnew, iiend, jjend
      INTEGER ierr

      CHARACTER(LEN=11) maskfile

      character(len=10) bfmmask_file

      CHARACTER(LEN=50) filename
      CHARACTER(LEN=3), DIMENSION(7) :: var_nc
      CHARACTER(LEN=5) nomevar01
      LOGICAL B

! -------------------
! To read only one meshmask.nc

      maskfile        = 'meshmask.nc'

      bfmmask_file = 'bfmmask.nc'

      filename='BC/TI_yyyy0215-00:00:00.nc' ! 26 chars


      iiend = MIN(jpi+nimpp-1, jpiglo)
      jjend = MIN(jpj+njmpp-1, jpjglo)

            do ii=nimpp, iiend
         do jj=njmpp, jjend
      do kk=1, jpk
               iinew = ii - nimpp + 1
               jjnew = jj - njmpp + 1
               idxt2glo(kk,jjnew,iinew,3)=ii !
               idxt2glo(kk,jjnew,iinew,2)=jj ! matrix to go from local to global
               idxt2glo(kk,jjnew,iinew,1)=kk !
            enddo
         enddo
      enddo

! 1. Horzontal grid-point position
! --------------------------------
      DO wr_procs=1, nodes
                if (myrank==writing_procs(wr_procs))then
                        call readnc_global_double_2d(maskfile,'glamt',totglamt)
                        call readnc_global_double_2d(maskfile,'gphit',totgphit)
                end if
      end do

      call readnc_slice_double_2d(maskfile,'glamt', glamt)
      call readnc_slice_double_2d(maskfile,'glamu', glamu)
      call readnc_slice_double_2d(maskfile,'glamv', glamv)
      !call readnc_slice_double_2d(maskfile,'glamf', glamf)

      call readnc_slice_double_2d(maskfile,'gphit', gphit)
      call readnc_slice_double_2d(maskfile,'gphiu', gphiu)
      call readnc_slice_double_2d(maskfile,'gphiv', gphiv)
      !call readnc_slice_double_2d(maskfile,'gphif', gphif)


! 2. Horizontal scale factors
! ---------------------------
      call readnc_slice_double_2d(maskfile,'e1t', e1t)
      call readnc_slice_double_2d(maskfile,'e1u', e1u)
      call readnc_slice_double_2d(maskfile,'e1v', e1v)
      !call readnc_slice_double_2d(maskfile,'e1f', e1f)

      call readnc_slice_double_2d(maskfile,'e2t', e2t)
      call readnc_slice_double_2d(maskfile,'e2u', e2u)
      call readnc_slice_double_2d(maskfile,'e2v', e2v)
      !call readnc_slice_double_2d(maskfile,'e2f', e2f)



! 3. masks
! --------

      CALL readnc_slice_int1 (maskfile,'umask', umask )
      CALL readnc_slice_int1 (maskfile,'vmask', vmask )
!      CALL readnc_slice_double (maskfile,'fmask', fmask )
      CALL readnc_slice_int1 (maskfile,'tmask', tmask )
!      CALL readnc_global_double(maskfile,'tmask', tmaskglo)

      call readnc_slice_logical(bfmmask_file, 'bfmmask', bfmmask)

!      Initialization of mbathy
      mbathy(:,:) = 0
      NWATERPOINTS=0
     do ii=1, jpi
       do jj=1, jpj
        do kk=1, jpk
         if (tmask(kk,jj,ii).eq.1) then
            mbathy(jj,ii) = mbathy(jj,ii) +1
            NWATERPOINTS = NWATERPOINTS +1
         endif
        enddo
       enddo
      enddo

!      CALL readnc_slice_double_2d(maskfile,'ff', ff )


! 4. depth and vertical scale factors
! -----------------------------------

#ifdef gdept1d
       CALL readmask_double_1d(maskfile,'gdept', gdept)
#else
      CALL readnc_slice_double(maskfile,'gdept', gdept)
#endif
       CALL readmask_double_1d(maskfile,'gdepw', gdepw)
       jpk_eu = 0
       do kk=1,jpk
#ifdef gdept1d
          if (gdept(kk).lt.Euphotic_lev)  jpk_eu=kk
#else
          if (gdept(kk,1,1).lt.Euphotic_lev)  jpk_eu=kk
#endif

       enddo
       if (lwp) write(*,*) 'Euphotic level at k = ', jpk_eu

      CALL readnc_slice_double (maskfile,'e3t_0', e3t_0 )
      CALL readnc_slice_double (maskfile,'e3u_0', e3u_0 )
      CALL readnc_slice_double (maskfile,'e3v_0', e3v_0 )
      CALL readnc_slice_double (maskfile,'e3w_0', e3w_0 )

      IF (.not.IS_FREE_SURFACE) then
         e3t = e3t_0
         e3t_back = e3t
         e3u = e3u_0
         e3v = e3v_0
         e3w = e3w_0
      ENDIF


      h_column = 0.0
      DO ii= 1,jpi
      DO jj= 1,jpj
      DO kk=1,mbathy(jj,ii)
           h_column(jj,ii) = h_column(jj,ii) + e3t_0(kk,jj,ii)
      ENDDO
      ENDDO
      ENDDO



!       Restoration Mask ****************

      ! resto is kept just to provide compliance with bfmv2, but should be removed with bfmv5
      var_nc(1) = 'O2o'
      var_nc(2) = 'N1p'
      var_nc(3) = 'N3n'
      var_nc(4) = 'N5s'
      var_nc(5) = 'O3c'
      var_nc(6) = 'O3h'

      IF (NWATERPOINTS.GT.0) THEN
      do jn=1,jn_gib

         nomevar01='re'//var_nc(jn)
         call readnc_slice_float('bounmask.nc',nomevar01,resto(:,:,:,jn),0)

      enddo
      ELSE
        resto=0.0
      ENDIF

      call readnc_slice_int   ('bounmask.nc','index',idxt)
      


! ************************************ BFM points re-indexing *******

      NBFMPOINTS = BFM_count()
      call myalloc_BFM()
      B=BFM_Indexing()

!************************************ North,South,East,West Boundaries ****
#ifdef key_mpp
      B = SENDRECV_count()
      call myalloc_sendrecv()
      B = SENDRECV_Indexing()
#endif


! *********************************   Gibraltar area
      !filename  ='BC/GIB_'//TC_GIB%TimeStrings(1)//'.nc'


      !if (lwp) write(*,*) 'domrea->filename: ', filename, '    '

      !CALL readnc_int_1d(filename, 'gib_idxt_N1p', Gsizeglo, gib_idxtglo)

      !if (lwp) write(*,*) 'domrea->readnc_int_1d  finita'
      !if (lwp) write(*,*) 'domrea->Gsizeglo', Gsizeglo

      !Gsize = COUNT_InSubDomain_GIB(Gsizeglo, gib_idxtglo)

      !write(*,*) 'domrea->Gsize   : ', Gsize, 'myrank=', myrank


      !if (Gsize.NE.0) then
          !if (lwp) write(*,*) 'domrea-> lancio alloc_DTATRC_local_gib'
          !call alloc_DTATRC_local_gib

          !B=GIBRe_indexing()

      !endif

! ********************************  Rivers ******
      !filename       ='BC/TIN_'//TC_TIN%TimeStrings(1)//'.nc'
      !print *,"---",Rsizeglo
      
      !CALL readnc_int_1d(filename, 'riv_idxt', Rsizeglo, riv_idxtglo)

      !Rsize = COUNT_InSubDomain(Rsizeglo,riv_idxtglo)

      !if (Rsize.NE. 0) then
          !call alloc_DTATRC_local_riv

          !B=RIVRe_Indexing()

      !endif
      !print *,Rsize,Rsizeglo


      !if(lwp) write(*,*) 'RIV finiti'

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      allocate(boundaries)
      boundaries = bc_set("boundaries.nml")

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

! ******************************************* Atmospherical inputs
      filename       = 'BC/ATM_'//TC_ATM%TimeStrings(1)//'.nc'
      ! CALL readnc_int_1d(filename, 'atm_idxt',Asizeglo,atm_idxtglo)
      ! Asize = COUNT_InSubDomain(Asizeglo,atm_idxtglo)

      ! if (Asize.NE. 0) then
      call alloc_DTATRC_local_atm
      !    write(*,*) 'domrea->ATMRE_Indexing ATM iniziata, myrank=', myrank
      !    B=ATMRe_Indexing()
      !    write(*,*) 'domrea->ATMRE_Indexing ATM finita, myrank=', myrank
      ! endif



! **************************************** FLUXES
      if (FSizeGlo.GT.0) then ! file could not exist
         CALL readnc_int_1d('Fluxes.nc','index',FSizeGlo,INDFluxGlo)
         Fsize = COUNT_InSubDomain(FSizeGlo, INDFluxGlo)
      else
         Fsize=0
      endif

      call MPI_ALLREDUCE(Fsize, FsizeMax, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)

      ! LEVEL1 write(*,*) 'myrank=', myrank, ' Fsize = ' , Fsize,' FsizeMax = ' , FsizeMax

      if (Fsize.NE.0) then
         call alloc_DIA_local_flx
      end if

      if (myrank ==0) call alloc_DIA_GLOBAL_flx()

      if (Fsize.NE.0) then
         B=FLXRe_Indexing()
         write(*,*) 'domrea->FLXRE_Indexing finita, myrank=', myrank

      endif

      call alloc_DIA_MPI_flx()

      write(*,*) 'DOMREA finita, myrank = ', myrank

   !   DEALLOCATE(idxt)
    !  DEALLOCATE(idxt2glo)
     ! DEALLOCATE(resto)


      CONTAINS

! *****************************************************************
!     FUNCTION COUNT_InSubDomain
!     RETURNS the number of points of a specific boundary condition
!     in the subdomain of the current processor
! *****************************************************************
      INTEGER FUNCTION COUNT_InSubDomain(sizeGLO,idxtGLOBAL)
          USE modul_param , ONLY: jpk,jpj,jpi
          USE myalloc     , ONLY: idxt

          IMPLICIT NONE
          INTEGER, INTENT(IN) :: sizeGLO
          INTEGER, INTENT(IN) :: idxtGLOBAL(sizeGLO)

          ! local
          INTEGER kk,jj,ii,jv
          INTEGER counter,junk

           counter = 0
             do ii =1, jpi
            do jj =1, jpj
           do kk =1, jpk
                if (tmask(kk,jj,ii).eq.1) then
                  junk = idxt(kk,jj,ii)
                  do jv =1, sizeGLO
                    if (junk.EQ.idxtGLOBAL(jv)) then
                       counter = counter + 1
                    endif
                  enddo
                endif
             enddo
            enddo
           enddo

          COUNT_InSubDomain = counter

      END FUNCTION COUNT_InSubDomain

! *************************************************************************

! *****************************************************************
!     FUNCTION COUNT_InSubDomain_GIB
!     RETURNS the number of points of a specific boundary condition
!     in the subdomain of the current processor

!     It is identical to COUNT_InSubDomain, apart from the
!     EXIT command, put in order to save time about GIB count at 1/24


!     This EXIT command implies that we can apply an only BC to each cell.
!     Without it we can have e.g. more that one river condition on each cell,
!     useful in case of degratated model

! *****************************************************************
      INTEGER FUNCTION COUNT_InSubDomain_GIB(sizeGLO,idxtGLOBAL)
          USE modul_param , ONLY: jpk,jpj,jpi
          USE myalloc     , ONLY: idxt

          IMPLICIT NONE
          INTEGER, INTENT(IN) :: sizeGLO
          INTEGER, INTENT(IN) :: idxtGLOBAL(sizeGLO)

          ! local
          INTEGER kk,jj,ii,jv
          INTEGER counter,junk

           counter = 0
             do ii =1, jpi
            do jj =1, jpj
           do kk =1, jpk
                if (tmask(kk,jj,ii).eq.1) then
                  junk = idxt(kk,jj,ii)
                  do jv =1, sizeGLO
                    if (junk.EQ.idxtGLOBAL(jv)) then
                       counter = counter + 1
                       EXIT
                    endif
                  enddo
                endif
             enddo
            enddo
           enddo

          COUNT_InSubDomain_GIB = counter

      END FUNCTION COUNT_InSubDomain_GIB

! *************************************************************************

      LOGICAL FUNCTION RIVRE_Indexing()
          IMPLICIT NONE
          ! local
          INTEGER kk,jj,ii,jv
          INTEGER counter,junk


          counter=0
          kk=1
             do ii =1, jpi
            do jj =1, jpj
                junk = idxt(kk,jj,ii)
                do jv =1, RsizeGLO
                   if ( junk.EQ.riv_idxtglo(jv) )  then
                      counter = counter + 1
                      riv_ridxt(1,counter) = jv
                      riv_ridxt(2,counter) = kk
                      riv_ridxt(3,counter) = jj
                      riv_ridxt(4,counter) = ii
                   endif
              enddo
            enddo
           enddo


      RIVRE_Indexing = .true.
      END FUNCTION RIVRE_Indexing
! *************************************************************************
      LOGICAL FUNCTION GIBRE_Indexing()
          IMPLICIT NONE

          ! local
          INTEGER kk,jj,ii,jv
          INTEGER counter,junk


          counter=0
             do ii =1, jpi
            do jj =1, jpj
           do kk =1, jpk
                junk = idxt(kk,jj,ii)
                do jv =1, Gsizeglo
                   if ( junk.EQ.gib_idxtglo(jv) )  then
                      counter = counter + 1
                      gib_ridxt(1,counter) = jv
                      gib_ridxt(2,counter) = kk
                      gib_ridxt(3,counter) = jj
                      gib_ridxt(4,counter) = ii
                   endif
              enddo
             enddo
            enddo
           enddo


      GIBRE_Indexing= .true.
      END FUNCTION GIBRE_Indexing


! *************************************************************************
      LOGICAL FUNCTION  FLXRE_Indexing()

          IMPLICIT NONE

          ! local
          INTEGER kk,jj,ii,jv
          INTEGER counter,junk


          counter=0
             do ii =1, jpi
            do jj =1, jpj
           do kk =1, jpk
                junk = idxt(kk,jj,ii)
                do jv =1, FsizeGLO
                   if ( junk.EQ.INDFluxGlo(jv) )  then
                      counter = counter + 1
                      flx_ridxt(counter,1) = jv
                      flx_ridxt(counter,2) = kk
                      flx_ridxt(counter,3) = jj
                      flx_ridxt(counter,4) = ii
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
          INTEGER kk,jj,ii
          INTEGER counter


          counter=0
          if (atlantic_bfm) then

                 do ii =2, jpi-1
                do jj =1, jpj-1
               do kk =1, jpkb-1

                   if (tmask(kk,jj,ii).EQ.1.0 ) then
                      counter = counter + 1
                      BFMpoints(1,counter) = kk
                      BFMpoints(2,counter) = jj
                      BFMpoints(3,counter) = ii

                   endif

                 enddo
                enddo!idxt
               enddo
           else
           ! NO ACTIVATION IN ATLANTIC BUFFER
               do kk =1, jpkb-1
                do jj =1, jpj-1
                 do ii =2, jpi-1

                   if ( (tmask(kk,jj,ii).EQ.1.0 ) .and. (resto(kk,jj,ii,1).eq.0.0 )) then
                      counter = counter + 1
                      BFMpoints(1,counter) = kk
                      BFMpoints(2,counter) = jj
                      BFMpoints(3,counter) = ii

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
          INTEGER kk,jj,ii
          INTEGER :: counter = 0


       if (atlantic_bfm) then
             do ii =2, jpi-1
            do jj =1, jpj-1
           do kk =1, jpkb-1
            if (kk.eq.2)  NBFMPOINTS_SUP = counter

               if (tmask(kk,jj,ii).EQ.1.0 ) counter = counter + 1

             enddo
            enddo
           enddo
        else
           ! NO ACTIVATION IN ATLANTIC BUFFER
             do ii =2, jpi-1
            do jj =1, jpj-1
           do kk =1, jpkb-1
            if (kk.eq.2)  NBFMPOINTS_SUP = counter

               if ( (tmask(kk,jj,ii).EQ.1.0 ) .and. (resto(kk,jj,ii,1).eq.0.0) )  counter = counter + 1

             enddo
            enddo
           enddo

        endif
        BFM_count =counter
      END FUNCTION BFM_count


! ***************************************************************
#ifdef key_mpp
          LOGICAL FUNCTION  SENDRECV_count()
          USE myalloc
          IMPLICIT NONE

          ! local
          INTEGER ii,jj,kk
          INTEGER :: counter

      counter=0
        do jj =1, jpj
           do kk =1, jpk-1
               if (tmask(kk,jj,2).EQ.1 ) counter = counter + 1
           enddo
        enddo
      WEST_count_send =counter

      counter=0
         do jj =1, jpj
           do kk =1, jpk-1
               if (tmask(kk,jj,1).EQ.1 ) counter = counter + 1
           enddo
         enddo
      WEST_count_recv =counter




      counter=0
         do jj =1, jpj
           do kk =1, jpk-1
               if (tmask(kk,jj,jpi-1).EQ.1 ) counter = counter + 1
           enddo
         enddo
      EAST_count_send = counter
      counter=0
         do jj =1, jpj
           do kk =1, jpk-1
               if (tmask(kk,jj,jpi).EQ.1 ) counter = counter + 1
           enddo
         enddo
      EAST_count_recv = counter


      counter=0
         do ii =1, jpi
           do kk =1, jpk-1
               if (tmask(kk,2,ii).EQ.1 ) counter = counter + 1
           enddo
         enddo
      SOUTH_count_send = counter
      counter=0
         do ii =1, jpi
           do kk =1, jpk-1
               if (tmask(kk,1,ii).EQ.1 ) counter = counter + 1
           enddo
         enddo
      SOUTH_count_recv = counter


      counter=0
         do ii =1, jpi
           do kk =1, jpk-1
               if (tmask(kk,jpj-1,ii).EQ.1 ) counter = counter + 1
           enddo
         enddo
      NORTH_count_send = counter
      counter=0
         do ii =1, jpi
           do kk =1, jpk-1
               if (tmask(kk,jpj,ii).EQ.1 ) counter = counter + 1
           enddo
         enddo
      NORTH_count_recv = counter



      SENDRECV_count = .true.
      END FUNCTION SENDRECV_count



     LOGICAL FUNCTION  SENDRECV_Indexing()

          IMPLICIT NONE

          ! local
          INTEGER ii,jj,kk
          INTEGER counter


      counter=0
        do jj =1, jpj
           do kk =1, jpk-1
               if (tmask(kk,jj,2).EQ.1 ) then
                  counter = counter + 1
                  WESTpoints_send(1,counter) = jj
                  WESTpoints_send(2,counter) = kk
              endif
           enddo
        enddo
      WEST_count_send =counter

      counter=0
         do jj =1, jpj
           do kk =1, jpk-1
               if (tmask(kk,jj,1).EQ.1 ) then
                   counter = counter + 1
                   WESTpoints_recv(1,counter) = jj
                   WESTpoints_recv(2,counter) = kk
               endif
           enddo
         enddo
      WEST_count_recv =counter


      counter=0
         do jj =1, jpj
           do kk =1, jpk-1
               if (tmask(kk,jj,jpi-1).EQ.1 ) then
                   counter = counter + 1
                   EASTpoints_send(1,counter) = jj
                   EASTpoints_send(2,counter) = kk
               endif
           enddo
         enddo
      EAST_count_send = counter
      counter=0
         do jj =1, jpj
           do kk =1, jpk-1
               if (tmask(kk,jj,jpi).EQ.1 ) then
                   counter = counter + 1
                   EASTpoints_recv(1,counter) = jj
                   EASTpoints_recv(2,counter) = kk
               endif
           enddo
         enddo
      EAST_count_recv = counter


      counter=0
         do ii =1, jpi
           do kk =1, jpk-1
               if (tmask(kk,2,ii).EQ.1 ) then
                   counter = counter + 1
                   SOUTHpoints_send(1,counter) = ii
                   SOUTHpoints_send(2,counter) = kk
               endif
           enddo
         enddo
      SOUTH_count_send = counter
      counter=0
         do ii =1, jpi
           do kk =1, jpk-1
               if (tmask(kk,1,ii).EQ.1 ) then
                   counter = counter + 1
                   SOUTHpoints_recv(1,counter) = ii
                   SOUTHpoints_recv(2,counter) = kk
               endif
           enddo
         enddo
      SOUTH_count_recv = counter


      counter=0
         do ii =1, jpi
            do kk =1, jpk-1
               if (tmask(kk,jpj-1,ii).EQ.1 ) then
                   counter = counter + 1
                   NORTHpoints_send(1,counter) = ii
                   NORTHpoints_send(2,counter) = kk
               endif
           enddo
         enddo
      NORTH_count_send = counter
      counter=0
         do ii =1, jpi
            do kk =1, jpk-1
               if (tmask(kk,jpj,ii).EQ.1 ) then
                   counter = counter + 1
                   NORTHpoints_recv(1,counter) = ii
                   NORTHpoints_recv(2,counter) = kk
               endif
           enddo
         enddo
      NORTH_count_recv = counter


      !$acc update device( EASTpoints_send( 1:2, 1:EAST_count_send ) )
      !$acc update device( WESTpoints_send( 1:2, 1:WEST_count_send ) )
      !$acc update device( NORTHpoints_send( 1:2, 1:NORTH_count_send ) )
      !$acc update device( SOUTHpoints_send( 1:2, 1:SOUTH_count_send ) )
      
      !$acc update device( EASTpoints_recv( 1:2, 1:EAST_count_recv ) )
      !$acc update device( WESTpoints_recv( 1:2, 1:WEST_count_recv ) )
      !$acc update device( NORTHpoints_recv( 1:2, 1:NORTH_count_recv ) )
      !$acc update device( SOUTHpoints_recv( 1:2, 1:SOUTH_count_recv ) )
      


        SENDRECV_Indexing =.true.
      END FUNCTION SENDRECV_Indexing

#endif
      END SUBROUTINE domrea
