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
      USE DR_mem
      USE BC_mem
      USE TIME_MANAGER
      IMPLICIT NONE

! local declarations
! ==================

      INTEGER ii, jj, kk, jn, iinew, jjnew, iiend, jjend


      CHARACTER*11 maskfile
      CHARACTER*50 filename
      CHARACTER(LEN=3), DIMENSION(7) :: var_nc
      CHARACTER*5 nomevar01
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

!      nomevar01(1:5)='glamt'; CALL ioogsnc_mesh_2d_d(nomefile,nomevar01,1,jpjglo,jpiglo,glamtglo)
!      nomevar01(1:5)='glamu'; CALL ioogsnc_mesh_2d_d(nomefile,nomevar01,1,jpjglo,jpiglo,glamuglo)
!      nomevar01(1:5)='glamv'; CALL ioogsnc_mesh_2d_d(nomefile,nomevar01,1,jpjglo,jpiglo,glamvglo)
!      nomevar01(1:5)='glamf'; CALL ioogsnc_mesh_2d_d(nomefile,nomevar01,1,jpjglo,jpiglo,glamfglo)
!
!      do jj=jjsta, jjend
!         do ii=iista, iiend
!            iinew = ii - nimpp + 1
!            jjnew = jj - njmpp + 1
!            glamt(iinew, jjnew) = glamtglo(ii, jj)
!            glamu(iinew, jjnew) = glamuglo(ii, jj)
!            glamv(iinew, jjnew) = glamvglo(ii, jj)
!            glamf(iinew, jjnew) = glamfglo(ii, jj)
!         enddo
!      enddo


      call readnc_slice_double_2d(maskfile,'gphit', gphit)
      call readnc_slice_double_2d(maskfile,'gphiu', gphiu)
      call readnc_slice_double_2d(maskfile,'gphiv', gphiv)
      call readnc_slice_double_2d(maskfile,'gphif', gphif)

!      nomevar01(1:5) = 'gphit'; CALL ioogsnc_mesh_2d_d(nomefile,nomevar01,1,jpjglo,jpiglo,gphitglo)
!      nomevar01(1:5) = 'gphiu'; CALL ioogsnc_mesh_2d_d(nomefile,nomevar01,1,jpjglo,jpiglo,gphiuglo)
!      nomevar01(1:5) = 'gphiv'; CALL ioogsnc_mesh_2d_d(nomefile,nomevar01,1,jpjglo,jpiglo,gphivglo)
!      nomevar01(1:5) = 'gphif'; CALL ioogsnc_mesh_2d_d(nomefile,nomevar01,1,jpjglo,jpiglo,gphifglo)
!
!      do jj=jjsta, jjend
!         do ii=iista, iiend
!            iinew = ii - nimpp + 1
!            jjnew = jj - njmpp + 1
!            gphit(iinew, jjnew) = gphitglo(ii, jj)
!            gphiu(iinew, jjnew) = gphiuglo(ii, jj)
!            gphiv(iinew, jjnew) = gphivglo(ii, jj)
!            gphif(iinew, jjnew) = gphifglo(ii, jj)
!         enddo
!      enddo

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


!      nomevar02(1:3)='e1t'; CALL ioogsnc_mesh_2d_d(nomefile, nomevar02,1,jpjglo,jpiglo,e1tglo)
!      nomevar02(1:3)='e1u'; CALL ioogsnc_mesh_2d_d(nomefile, nomevar02,1,jpjglo,jpiglo,e1uglo)
!      nomevar02(1:3)='e1v'; CALL ioogsnc_mesh_2d_d(nomefile, nomevar02,1,jpjglo,jpiglo,e1vglo)
!      nomevar02(1:3)='e1f'; CALL ioogsnc_mesh_2d_d(nomefile, nomevar02,1,jpjglo,jpiglo,e1fglo)
!
!      do jj=jjsta, jjend
!         do ii=iista, iiend
!            iinew = ii - nimpp + 1
!            jjnew = jj - njmpp + 1
!            e1t(iinew, jjnew) = e1tglo(ii, jj)
!            e1u(iinew, jjnew) = e1uglo(ii, jj)
!            e1v(iinew, jjnew) = e1vglo(ii, jj)
!            e1f(iinew, jjnew) = e1fglo(ii, jj)
!         enddo
!      enddo

!      nomevar02(1:3)='e2t'; CALL ioogsnc_mesh_2d_d(nomefile, nomevar02,1,jpjglo,jpiglo,e2tglo)
!      nomevar02(1:3)='e2u'; CALL ioogsnc_mesh_2d_d(nomefile, nomevar02,1,jpjglo,jpiglo,e2uglo)
!      nomevar02(1:3)='e2v'; CALL ioogsnc_mesh_2d_d(nomefile, nomevar02,1,jpjglo,jpiglo,e2vglo)
!      nomevar02(1:3)='e2f'; CALL ioogsnc_mesh_2d_d(nomefile, nomevar02,1,jpjglo,jpiglo,e2fglo)
!
!      do jj=jjsta, jjend
!         do ii=iista, iiend
!            iinew = ii - nimpp + 1
!            jjnew = jj - njmpp + 1
!            e2t(iinew, jjnew) = e2tglo(ii, jj)
!            e2u(iinew, jjnew) = e2uglo(ii, jj)
!            e2v(iinew, jjnew) = e2vglo(ii, jj)
!            e2f(iinew, jjnew) = e2fglo(ii, jj)
!         enddo
!      enddo

! 3. masks
! --------
      write(*,*) 'reading 3D'

      CALL readnc_slice_double (maskfile,'umask', umask )
      CALL readnc_slice_double (maskfile,'vmask', vmask )
      CALL readnc_slice_double (maskfile,'fmask', fmask )
      CALL readnc_slice_double (maskfile,'tmask', tmask )
      CALL readnc_global_double(maskfile,'tmask', tmaskglo)
!      nomevar01(1:5) = 'tmask'; CALL ioogsnc_mesh_3d_d(nomefile, nomevar01,1,jpk,jpjglo,jpiglo,tmaskglo)
!      nomevar01(1:5) = 'umask'; CALL ioogsnc_mesh_3d_d(nomefile, nomevar01,1,jpk,jpjglo,jpiglo,umaskglo)
!      nomevar01(1:5) = 'vmask'; CALL ioogsnc_mesh_3d_d(nomefile, nomevar01,1,jpk,jpjglo,jpiglo,vmaskglo)
!      nomevar01(1:5) = 'fmask'; CALL ioogsnc_mesh_3d_d(nomefile, nomevar01,1,jpk,jpjglo,jpiglo,fmaskglo)

!      do kk=1, jpk
!         do jj=jjsta, jjend
!            do ii=iista, iiend
!               iinew = ii - nimpp + 1
!               jjnew = jj - njmpp + 1
!               tmask(iinew, jjnew, kk) = tmaskglo(ii, jj, kk)
!               umask(iinew, jjnew, kk) = umaskglo(ii, jj, kk)
!               vmask(iinew, jjnew, kk) = vmaskglo(ii, jj, kk)
!               fmask(iinew, jjnew, kk) = fmaskglo(ii, jj, kk)
!            enddo
!         enddo
!      enddo

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
!      nomevar03 = 'ff' ; nomevar03(3:3) = CHAR(0)
!      CALL ioogsnc_mesh_2d_d(nomefile, nomevar03,1,jpjglo,jpiglo,ffglo)
!
!      do jj=jjsta, jjend
!         do ii=iista, iiend
!            iinew = ii - nimpp + 1
!            jjnew = jj - njmpp + 1
!            ff(iinew, jjnew) = ffglo(ii, jj)
!         enddo
!      enddo


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


!      nomevar01(1:5)='gdept'; CALL ioogsnc_mesh_1d_d(nomefile, nomevar01,1,jpk,gdept)
!      nomevar01(1:5)='gdepw'; CALL ioogsnc_mesh_1d_d(nomefile, nomevar01,1,jpk,gdepw)
!      nomevar02(1:3) = 'e3t'; CALL ioogsnc_mesh_1d_d(nomefile, nomevar02,1,jpk,e3t)
!      nomevar02(1:3) = 'e3w'; CALL ioogsnc_mesh_1d_d(nomefile, nomevar02,1,jpk,e3w)
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

!         CALL ioogsnc_mesh_3d_d(nomefile, nomevar01,1,jpk,jpjglo,jpiglo,restoglo)
!
!         do kk=1, jpk
!            do jj=jjsta, jjend
!               do ii=iista, iiend
!                  iinew = ii - nimpp + 1
!                  jjnew = jj - njmpp + 1
!                  resto(iinew, jjnew, kk,jn) = restoglo(ii, jj, kk)
!               enddo
!            enddo
!         enddo
      enddo
      call readnc_slice_int   ('bounmask.nc','index',idxt);
!      nomevar01(1:5)='index'; CALL ioogsnc_idx_3d(nomefile, nomevar01,1,jpk,jpjglo,jpiglo,idxtglo)
!
!      do kk=1, jpk
!         do jj=jjsta, jjend
!            do ii=iista, iiend
!               iinew = ii - nimpp + 1
!               jjnew = jj - njmpp + 1
!               idxt(iinew, jjnew, kk) = idxtglo(ii, jj, kk)
!            enddo
!         enddo
!      enddo




! *********************************   Gibraltar area
      filename  ='BC/GIB_'//TC_GIB%TimeStrings(1)//'.nc'


       if (lwp) write(*,*) 'domrea->filename: ', filename, '    '


      !CALL ioogsnc_idx_1d2(filename, nomevar05,Gsizeglo, gib_idxtglo)
      CALL readnc_int_1d(filename, 'gib_idxt_N1p', Gsizeglo, gib_idxtglo)

      if (lwp) write(*,*) 'domrea->ioogsnc_idx_1d  finita'
      if (lwp) write(*,*) 'domrea->Gsizeglo', Gsizeglo

      Gsize = COUNT_InSubDomain(Gsizeglo, gib_idxtglo)

      write(*,*) 'domrea->Gsize   : ', Gsize, 'rank=', rank


      if (Gsize.NE.0) then
          if (lwp) write(*,*) 'domrea-> lancio alloc_DTATRC_local_gib'
          call alloc_DTATRC_local_gib

          write(*,*) 'domrea->GIBRE_Indexing GIB iniziata, rank=', rank
          B=GIBRe_indexing()
          write(*,*) 'domrea->GIBRE_Indexing GIB finita, rank=', rank
      endif

! ********************************  Rivers ******
      filename       ='BC/TIN_'//TC_TIN%TimeStrings(1)//'.nc';
      CALL readnc_int_1d(filename, 'riv_idxt', Rsizeglo, riv_idxtglo)
      !CALL ioogsnc_idx_1d2(filename,nomevar04, Rsizeglo,riv_idxtglo)
      Rsize = COUNT_InSubDomain(Rsizeglo,riv_idxtglo)

      if (Rsize.NE. 0) then
          call alloc_DTATRC_local_riv
          write(*,*) 'domrea->RIVRE_Indexing RIV iniziata, rank=', rank
          B=RIVRe_Indexing()
          write(*,*) 'domrea->RIVRE_Indexing RIV finita, rank=', rank
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

       if(lwp) write(*,*) 'DOMREA finita, wait for flushing...'

       call flush()

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


      LOGICAL FUNCTION RIVRE_Indexing()
!          USE modul_param , ONLY: jpi, jpj, jpk
!          USE myalloc     , ONLY: idxt
!          USE myalloc_mpp , ONLY: rank
!          USE BC_mem      , ONLY: Rsizeglo,riv_idxtglo,riv_ridxt
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

      LOGICAL FUNCTION GIBRE_Indexing()
!          USE modul_param , ONLY: jpi, jpj, jpk
!          USE myalloc     , ONLY: idxt
!          USE myalloc_mpp , ONLY: rank
!          USE BC_mem      , ONLY: Gsizeglo,gib_idxtglo,gib_ridxt
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



          LOGICAL FUNCTION  ATMRE_Indexing()
!          USE modul_param , ONLY: jpi, jpj, jpk
!          USE myalloc     , ONLY: idxt
!          USE myalloc_mpp , ONLY: rank
!          USE BC_mem      , ONLY: Asizeglo,atm_idxtglo,atm_ridxt
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



      END SUBROUTINE domrea

! ************************************************************
!      SUBROUTINE RE_Indexing
!      global indexes ---> local (processor subdomain) indexes
! ************************************************************
      SUBROUTINE RE_Indexing(sizeGLO, idxtGLOBAL, N, ridxt)
          USE modul_param , ONLY: jpi, jpj, jpk
          USE myalloc     , ONLY: idxt
          USE myalloc_mpp , ONLY: rank
          IMPLICIT NONE

          INTEGER, INTENT(IN) :: sizeGLO
          INTEGER, INTENT(IN) :: idxtGLOBAL(sizeGLO)
          INTEGER, INTENT(IN) :: N
          INTEGER, INTENT(OUT):: ridxt(N,4)
          ! local
          INTEGER ii,jj,kk,jv
          INTEGER counter,junk


          counter=0
           do kk =1, jpk
            do jj =1, jpj
             do ii =1, jpi
                junk = idxt(ii,jj,kk)
                do jv =1, sizeGLO
                   if ( junk.EQ.idxtGLOBAL(jv) )  then
                      counter = counter + 1
                      ridxt(counter,1) = jv
                      ridxt(counter,2) = ii
                      ridxt(counter,3) = jj
                      ridxt(counter,4) = kk
                   endif
              enddo
             enddo
            enddo
           enddo


      END SUBROUTINE RE_Indexing

