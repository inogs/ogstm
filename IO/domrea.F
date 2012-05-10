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

      call readnc_slice_double_2d(maskfile,'glamt', glamt)
      call readnc_slice_double_2d(maskfile,'glamu', glamu)
      call readnc_slice_double_2d(maskfile,'glamv', glamv)
      call readnc_slice_double_2d(maskfile,'glamf', glamf)


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


! F79 Initialization of mbathy
      mbathy(:,:) = 0
      do jj=1, jpj
       do ii=1, jpi          
        do kk=1, jpk
         if (tmask(ii, jj, kk).NE.0.) then
            mbathy(ii,jj) = mbathy(ii,jj) +1
         endif
        enddo
       enddo
      enddo

      CALL readnc_slice_double_2d(maskfile,'ff', ff )


! 4. depth and vertical scale factors
! -----------------------------------

#ifdef key_s_coord
      CALL restget(nummsh,'hbatt',jpiglo,jpjglo,1,0,clog,hbattglo)
      CALL restget(nummsh,'hbatu',jpiglo,jpjglo,1,0,clog,hbatuglo)
      CALL restget(nummsh,'hbatv',jpiglo,jpjglo,1,0,clog,hbatvglo)
      CALL restget(nummsh,'hbatf',jpiglo,jpjglo,1,0,clog,hbatfglo)

      do jj=jjsta, jjend
         do ii=iista, iiend
            iinew = ii - nimpp + 1
            jjnew = jj - njmpp + 1
            hbatt(iinew, jjnew) = hbattglo(ii, jj)
            hbatu(iinew, jjnew) = hbatuglo(ii, jj)
            hbatv(iinew, jjnew) = hbatvglo(ii, jj)
            hbatf(iinew, jjnew) = hbatfglo(ii, jj)
         enddo
      enddo


      CALL restget(nummsh,'gsigt',1,1,jpk,0,clog,gsigt)
      CALL restget(nummsh,'gsigw',1,1,jpk,0,clog,gsigw)
      CALL restget(nummsh,'gsi3w',1,1,jpk,0,clog,gsi3w)
      CALL restget(nummsh,'esigt',1,1,jpk,0,clog,esigt)
      CALL restget(nummsh,'esigw',1,1,jpk,0,clog,esigw)
#  else


       CALL readmask_double_1d(maskfile,'gdept',gdept)
       CALL readmask_double_1d(maskfile,'gdepw',gdepw)
       CALL readmask_double_1d(maskfile,'e3t',    e3t)
       CALL readmask_double_1d(maskfile,'e3w',    e3w)



#endif
#ifdef key_trc_dmp


!       Restoration Mask ****************


      var_nc(1) = 'O2o'
      var_nc(2) = 'N1p'
      var_nc(3) = 'N3n'
      var_nc(4) = 'N5s'
      var_nc(5) = 'O3c'
      var_nc(6) = 'O3h'
      var_nc(7) = 'N6r'

      do jn=1,jn_gib

         nomevar01='re'//var_nc(jn)
         call readnc_slice_double('bounmask.nc',nomevar01,resto(:,:,:,jn));

      enddo
      call readnc_slice_int   ('bounmask.nc','index',idxt);





! *********************************   Gibraltar area
      filename  ='BC/GIB_'//TC_GIB%TimeStrings(1)//'.nc'


       if (lwp) write(*,*) 'domrea->filename: ', filename, '    '

      CALL readnc_int_1d(filename, 'gib_idxt_N1p', Gsizeglo, gib_idxtglo)

      if (lwp) write(*,*) 'domrea->ioogsnc_idx_1d  finita'
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
      !CALL ioogsnc_idx_1d2(filename,nomevar04, Rsizeglo,riv_idxtglo)
      Rsize = COUNT_InSubDomain(Rsizeglo,riv_idxtglo)

      if (Rsize.NE. 0) then
          call alloc_DTATRC_local_riv

          B=RIVRe_Indexing()

      endif


       if(lwp) write(*,*) 'RIV finiti'

! ******************************************* Atmospherical inputs
      filename       = 'BC/ATM_'//TC_ATM%TimeStrings(1)//'.nc'
      CALL readnc_int_1d(filename, 'atm_idxt',Asizeglo,atm_idxtglo)
!      CALL ioogsnc_idx_1d2(filename, nomevar04,Asizeglo,atm_idxtglo)
      Asize = COUNT_InSubDomain(Asizeglo,atm_idxtglo)

      if (Asize.NE. 0) then
         call alloc_DTATRC_local_atm
         write(*,*) 'domrea->ATMRE_Indexing ATM iniziata, rank=', rank
         B=ATMRe_Indexing()
         write(*,*) 'domrea->ATMRE_Indexing ATM finita, rank=', rank
      endif

#endif

! **************************************** FLUXES

      CALL readnc_int_1d('Fluxes.nc','index',FSizeGlo,INDFluxGlo)



      Fsize = COUNT_InSubDomain(FSizeGlo, INDFluxGlo)
      write(*,*) 'rank=', rank, ' Fsize = ' , Fsize

      if (Fsize.NE.0) then
         call alloc_DIA_local_flx
         write(*,*) 'domrea->FLXRE_Indexing iniziata, rank=', rank
         B=FLXRe_Indexing()
         write(*,*) 'domrea->FLXRE_Indexing finita, rank=', rank

      endif

       write(*,*) 'DOMREA finita, rank = ', rank

!       call flush()

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
                      riv_ridxt(counter,1) = jv
                      riv_ridxt(counter,2) = ii
                      riv_ridxt(counter,3) = jj
                      riv_ridxt(counter,4) = kk
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
                      gib_ridxt(counter,1) = jv
                      gib_ridxt(counter,2) = ii
                      gib_ridxt(counter,3) = jj
                      gib_ridxt(counter,4) = kk
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
                      atm_ridxt(counter,1) = jv
                      atm_ridxt(counter,2) = ii
                      atm_ridxt(counter,3) = jj
                      atm_ridxt(counter,4) = kk
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



      END SUBROUTINE domrea
