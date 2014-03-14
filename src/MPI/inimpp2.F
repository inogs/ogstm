CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/inimpp2.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE inimpp2
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE INIMPP2
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Lay out the global domain over processors.
CCC
CCC     FOR USING THIS VERSION, A PREPROCESSING TRAITMENT IS RECOMMENDED
CCC     FOR DEFINING BETTER CUTTING OUT.
CCC
CCC     This routine is used with a the bathymetry file.
CCC
CCC     In this version, the land processors are avoided and the adress
CCC     processor (nproc, narea,noea, ...) are calculated again.
CCC     The jpnij parameter can be lesser than jpni x jpnj
CCC     and this jpnij parameter must be calculated before with an
CCC     algoritmic preprocessing program.
CCC
CC   METHOD :
CC   -------
CC      Global domain is distributed in smaller local domains.
CC      Periodic condition is a function of the local domain position
CC      (global boundary or neighbouring domain) and of the global
CC      periodic
CC      Type :
CC                     jperio global periodic condition
CC                     nperio local  periodic condition
CC
CC
CC   INPUT :
CC   -----
CC      argument		: no
CC
CC   OUTPUT :
CC   ------
CC      common
CC            /COMDOM/          : domain parameters
CC                    nimpp     : longitudinal index 
CC                    njmpp     : latitudinal  index
CC                    nperio    : lateral condition type 
CC                    narea     : number for local area
CC                    nlci      : first dimension
CC                    nlcj      : second dimension
CC                    nproc     : number for local processor
CC                    noea      : number for local neighboring processor
CC                    nowe      : number for local neighboring processor
CC                    noso      : number for local neighboring processor
CC                    nono      : number for local neighboring processor
CC
CC   WORKSPACE :
CC   ---------
CC      local
CC              ji, jj, jn, jproc
CC              ii, ij, ifreq, il1, il2, iresti, irestj, zidom, zjdom
CC              nlcit(), nlcjt(), nimppt(), njmppt()

CCC      USE calendarmpp
      USE calendar
      USE myalloc
        IMPLICIT NONE
CC----------------------------------------------------------------------

      INTEGER nbne, nbnw, nbse, nbsw, npne, npnw, npse, npsw
      INTEGER, PARAMETER :: nsolv = 3

#ifdef key_mpp
      CHARACTER(LEN=25) clexp
C
      INTEGER ji,jj,jn,jproc,jarea
      INTEGER ii,ij,ifreq,il1,il2,iim,ijm,icont,ili,ilj
      INTEGER isurf,ijm1,imil
      INTEGER iino,ijno,iiso,ijso,iiea,ijea,iiwe,ijwe
      INTEGER iinw,ijnw,iine,ijne,iisw,ijsw,iise,ijse
      INTEGER iresti,irestj,iproc
      INTEGER imask(jpiglo,jpjglo),idata(jpidta,jpjdta)
      INTEGER iimppt(jpni,jpnj),ijmppt(jpni,jpnj)
      INTEGER ilci  (jpni,jpnj),ilcj  (jpni,jpnj)
      INTEGER ipproc(jpni,jpnj)
      INTEGER ibondj(jpni,jpnj),ibondi(jpni,jpnj)
      INTEGER ipolj(jpni,jpnj)
      INTEGER ilei(jpni,jpnj),ilej(jpni,jpnj)
      INTEGER ildi(jpni,jpnj),ildj(jpni,jpnj)
      INTEGER ioea(jpni,jpnj),iowe(jpni,jpnj)
      INTEGER ioso(jpni,jpnj),iono(jpni,jpnj)
      INTEGER ione(jpni,jpnj),ionw(jpni,jpnj)
      INTEGER iose(jpni,jpnj),iosw(jpni,jpnj)
      INTEGER ibne(jpni,jpnj),ibnw(jpni,jpnj)
      INTEGER ibse(jpni,jpnj),ibsw(jpni,jpnj)
      INTEGER iin(jpnij),ijn(jpnij)
C
      INTEGER iglo(4),iloc(4),iabsf(4),iabsl(4),ihals(4),ihale(4)
C
      REAL(8)    zidom ,zjdom

C 0. initialisation
C -----------------
C
#ifdef key_mpp_shmem
C Initialisation of shmem array
C
      CALL mppshmem
#endif
C
C read bathymetry file
C
      REWIND numbat
      READ(numbat,9001) clexp,iim,ijm
      READ(numbat,'(/)')
      ifreq=40
      il1=1
      DO jn=1,jpidta/ifreq+1
        READ(numbat,'(/)')
        il2=min(jpidta,il1+ifreq-1)
        READ(numbat,9002) (ii,ji=il1,il2,5)
        READ(numbat,'(/)')
        DO jj=jpjdta,1,-1
          READ(numbat,9003) ij,(idata(ji,jj),ji=il1,il2)
        END DO
        il1=il1+ifreq
      END DO
 9001 FORMAT(1x,a15,2i8)
 9002 FORMAT(3x,13(i3,12x))
 9003 FORMAT(i3,41i3)
      DO jj=1,jpjglo
        DO ji=1,jpiglo
          imask(ji,jj)=min(1,max(0,idata(ji+nizoom-1,jj+njzoom-1)))
        END DO
      END DO
C
C
C  1. Dimension arrays for subdomains
C -----------------------------------
C
C  Computation of local domain sizes ilci() ilcj()
C  These dimensions depend on global sizes jpni,jpnj and jpiglo,jpjglo
C  The subdomains are squares leeser than or equal to the global
C  dimensions divided by the number of processors minus the overlap
C  array.
C
      nreci=2*jpreci
      nrecj=2*jprecj
      iresti = mod ( jpiglo - nreci , jpni )
      irestj = mod ( jpjglo - nrecj , jpnj )
C
      IF (iresti.eq.0) iresti = jpni
      DO jj=1,jpnj
        DO ji=1,iresti
          ilci(ji,jj) = jpi
        END DO
        DO ji=iresti+1,jpni
          ilci(ji,jj) = jpi -1
        END DO
      END DO

      IF (irestj.eq.0) irestj = jpnj
      DO ji=1,jpni
        DO jj=1,irestj
          ilcj(ji,jj) = jpj
        END DO
        DO jj=irestj+1,jpnj
          ilcj(ji,jj) = jpj -1
        END DO
      END DO
C
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' inimpp: defines mpp subdomains'
      IF(lwp)WRITE(numout,*) ' ******  ----------------------'
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)"iresti=",iresti," irestj=",irestj
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)"jpni=",jpni," jpnj=",jpnj
C
      zidom=nreci
      DO ji=1,jpni
        zidom = zidom + ilci(ji,1) - nreci
      END DO
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)" sum ilci(i,1)=",zidom," jpiglo=",jpiglo
C
      zjdom=nrecj
      DO jj=1,jpnj
        zjdom = zjdom + ilcj(1,jj) - nrecj
      END DO
      IF(lwp)WRITE(numout,*)" sum ilci(1,j)=",zjdom," jpjglo=",jpjglo
      IF(lwp)WRITE(numout,*)
C
C
C  2. Index arrays for subdomains
C -------------------------------
C
      DO jj=1,jpnj
        DO ji=1,jpni
          iimppt(ji,jj)=1
          ijmppt(ji,jj)=1
          ipproc(ji,jj)=-1
        END DO
      END DO
C
      IF(jpni.gt.1)THEN
          DO jj=1,jpnj
            DO ji=2,jpni
              iimppt(ji,jj)=iimppt(ji-1,jj)+ilci(ji-1,jj)-nreci
            END DO
          END DO
      ENDIF
C
      IF(jpnj.gt.1)THEN
          DO jj=2,jpnj
            DO ji=1,jpni
              ijmppt(ji,jj)=ijmppt(ji,jj-1)+ilcj(ji,jj-1)-nrecj
            END DO
          END DO
      ENDIF
C
C
C 3. Subdomain description in the Regular Case
C --------------------------------------------
C
      nperio=0
      icont=-1
      DO jarea=1,jpni*jpnj
        ii=1+mod(jarea-1,jpni)
        ij=1+   (jarea-1)/jpni
        ili=ilci(ii,ij)
        ilj=ilcj(ii,ij)
C
        ibondj(ii,ij)=-1
        IF(jarea.gt.jpni) ibondj(ii,ij)=0
        IF(jarea.gt.(jpnj-1)*jpni) ibondj(ii,ij)=1
        IF(jpnj.eq.1) ibondj(ii,ij)=2
C
        ibondi(ii,ij)=0
        IF(mod(jarea,jpni).eq.1) ibondi(ii,ij)=-1
        IF(mod(jarea,jpni).eq.0) ibondi(ii,ij)= 1
        IF(jpni.eq.1) ibondi(ii,ij)= 2
C
C 2.4 Subdomain neighbors
C
        iproc=jarea-1
        ioso(ii,ij)=iproc-jpni
        iowe(ii,ij)=iproc-1
        ioea(ii,ij)=iproc+1
        iono(ii,ij)=iproc+jpni
        ildi(ii,ij)= 1 +jpreci
        ilei(ii,ij)=ili-jpreci
        ionw(ii,ij)=iono(ii,ij)-1
        ione(ii,ij)=iono(ii,ij)+1
        iosw(ii,ij)=ioso(ii,ij)-1
        iose(ii,ij)=ioso(ii,ij)+1
        ibsw(ii,ij)=1
        ibnw(ii,ij)=1
        IF(mod(iproc,jpni) .EQ. 0) THEN
            ibsw(ii,ij)=0
            ibnw(ii,ij)=0
        ENDIF
        ibse(ii,ij)=1
        ibne(ii,ij)=1
        IF(mod(iproc,jpni) .EQ. jpni-1) THEN
            ibse(ii,ij)=0
            ibne(ii,ij)=0
        ENDIF
        IF(iproc .LT. jpni) THEN
            ibsw(ii,ij)=0
            ibse(ii,ij)=0
        ENDIF
        IF(iproc .GE. (jpnj-1)*jpni) THEN
            ibnw(ii,ij)=0
            ibne(ii,ij)=0
        ENDIF
        IF(ibondi(ii,ij).eq.-1.or.ibondi(ii,ij).eq.2) ildi(ii,ij)=1
        IF(ibondi(ii,ij).eq. 1.or.ibondi(ii,ij).eq.2) ilei(ii,ij)=ili
        ildj(ii,ij)= 1  +jprecj
        ilej(ii,ij)=ilj-jprecj
        IF(ibondj(ii,ij).eq.-1.or.ibondj(ii,ij).eq.2) ildj(ii,ij)=1
        IF(ibondj(ii,ij).eq. 1.or.ibondj(ii,ij).eq.2) ilej(ii,ij)=ilj
C
C warning ii*ij (zone) .ne. nproc (processors)!
C
        IF(jperio.eq.1.OR.jperio.eq.4) THEN
            IF(jpni.eq.1)THEN
                ibondi(ii,ij)=2
            ELSE
                ibondi(ii,ij)=0
            ENDIF
            IF(mod(jarea,jpni).eq.0) THEN
                ioea(ii,ij)=iproc-(jpni-1)
                ione(ii,ij)=ione(ii,ij)-jpni
                iose(ii,ij)=iose(ii,ij)-jpni
            ENDIF
            IF(mod(jarea,jpni).eq.1) THEN
                iowe(ii,ij)=iproc+ jpni-1
                ionw(ii,ij)=ionw(ii,ij)+jpni
                iosw(ii,ij)=iosw(ii,ij)+jpni 
            ENDIF
            ibsw(ii,ij)=1
            ibnw(ii,ij)=1
            ibse(ii,ij)=1
            ibne(ii,ij)=1
            IF(iproc .LT. jpni) THEN
                ibsw(ii,ij)=0
                ibse(ii,ij)=0
            ENDIF
            IF(iproc .GE. (jpnj-1)*jpni) THEN
                ibnw(ii,ij)=0
                ibne(ii,ij)=0
            ENDIF
        ENDIF
        ipolj(ii,jj)=0
        IF(jperio.eq.3 .OR. jperio.eq.4) THEN
            ijm1=jpni*(jpnj-1)
            imil=ijm1+(jpni+1)/2
            IF(jarea.gt.ijm1) ipolj(ii,ij)=3
            IF(mod(jpni,2).eq.1.and.jarea.eq.imil) ipolj(ii,ij)=4
            IF(ipolj(ii,ij).eq.3) THEN
                iono(ii,ij)=jpni*jpnj-jarea+ijm1
            ENDIF
        ENDIF
C
        isurf=0
        DO jj=1,ilj
          DO  ji=1,ili
            IF(imask(ji+iimppt(ii,ij)-1,
     $         jj+ijmppt(ii,ij)-1).eq.1) isurf=isurf+1
          END DO
        END DO
        IF(isurf.ne.0) THEN
            icont=icont+1
            ipproc(ii,ij)=icont
            iin(icont+1)=ii
            ijn(icont+1)=ij
        ENDIF
      END DO
C
C... Control
      IF(icont+1.ne.jpnij) THEN
          IF(lwp) THEN 
              WRITE(numout,*) ' Eliminate land processors algorithm'
              WRITE(numout,*)
              WRITE(numout,*) ' jpni =',jpni,' jpnj =',jpnj
              WRITE(numout,*) ' jpnij =',jpnij, '< jpni x jpnj' 
              WRITE(numout,*)
              WRITE(numout,*) ' E R R O R '
              WRITE(numout,*) ' ***********, inimpp finds jpnij=',icont
     $            +1
              WRITE(numout,*) ' we stop'
          ENDIF 
          STOP 'inimpp'
      ENDIF 
C
C
C 4. Subdomain print
C ------------------
C
      IF(lwp) THEN
          ifreq=4
          il1=1
          DO jn=1,(jpni-1)/ifreq+1
            il2=min(jpni,il1+ifreq-1)
            WRITE(numout,*)
            WRITE(numout,9401) (ji,ji=il1,il2)
            WRITE(numout,9400) ('***',ji=il1,il2-1)
            DO jj=1,jpnj
C             WRITE(numout,9400)
              WRITE(numout,9403) ('   ',ji=il1,il2-1)
              WRITE(numout,9402) jj,
     $            (ilci(ji,jj),ilcj(ji,jj),ji=il1,il2)
              WRITE(numout,9404) (ipproc(ji,jj),ji=il1,il2)
              WRITE(numout,9403) ('   ',ji=il1,il2-1)
              WRITE(numout,9400) ('***',ji=il1,il2-1)
C             WRITE(numout,9400)
            END DO
            il1=il1+ifreq
          END DO
 9400     FORMAT('     ***',20('*************',a3))
 9403     FORMAT('     *     ',20('         *   ',a3))
 9401     FORMAT('        ',20('   ',i3,'          '))
 9402     FORMAT(' ',i3,' *  ',20(i3,'  x',i3,'   *   '))
 9404     FORMAT('     *  ',20('      ',i3,'   *   '))
      ENDIF
C
C
C 5. neighbour treatment
C ----------------------
C
      DO jarea=1,jpni*jpnj
        iproc=jarea-1
        ii=1+mod(jarea-1,jpni)
        ij=1+   (jarea-1)/jpni
        IF(    ipproc(ii,ij).EQ.-1
     $      .AND.iono(ii,ij).GE.0
     $      .AND.iono(ii,ij).LE.jpni*jpnj-1) THEN
            iino=1+mod(iono(ii,ij),jpni)
            ijno=1+   (iono(ii,ij))/jpni
            IF(ibondj(iino,ijno).EQ.1) ibondj(iino,ijno)=2
            IF(ibondj(iino,ijno).EQ.0) ibondj(iino,ijno)=-1
        ENDIF
        IF(    ipproc(ii,ij).EQ.-1
     $      .AND.ioso(ii,ij).GE.0
     $      .AND.ioso(ii,ij).LE.jpni*jpnj-1) THEN
            iiso=1+mod(ioso(ii,ij),jpni)
            ijso=1+   (ioso(ii,ij))/jpni
            IF(ibondj(iiso,ijso).EQ.-1) ibondj(iiso,ijso)=2
            IF(ibondj(iiso,ijso).EQ.0 ) ibondj(iiso,ijso)=1
        ENDIF
        IF(    ipproc(ii,ij).EQ.-1
     $      .AND.ioea(ii,ij).GE.0
     $      .AND.ioea(ii,ij).LE.jpni*jpnj-1) THEN
            iiea=1+mod(ioea(ii,ij),jpni)
            ijea=1+   (ioea(ii,ij))/jpni
            IF(ibondi(iiea,ijea).EQ.1) ibondi(iiea,ijea)=2
            IF(ibondi(iiea,ijea).EQ.0) ibondi(iiea,ijea)=-1
        ENDIF
        IF(    ipproc(ii,ij).EQ.-1
     $      .AND.iowe(ii,ij).GE.0
     $      .AND.iowe(ii,ij).LE.jpni*jpnj-1) THEN
            iiwe=1+mod(iowe(ii,ij),jpni)
            ijwe=1+   (iowe(ii,ij))/jpni
            IF(ibondi(iiwe,ijwe).EQ.-1) ibondi(iiwe,ijwe)=2
            IF(ibondi(iiwe,ijwe).EQ.0 ) ibondi(iiwe,ijwe)=1
        ENDIF
        IF(    ipproc(ii,ij).EQ.-1
     $      .AND.ibne(ii,ij).EQ.1) THEN
            iine=1+mod(ione(ii,ij),jpni)
            ijne=1+   (ione(ii,ij))/jpni
            IF(ibsw(iine,ijne).EQ.1) ibsw(iine,ijne)=0
        ENDIF
        IF(    ipproc(ii,ij).EQ.-1
     $      .AND.ibsw(ii,ij).EQ.1) THEN
            iisw=1+mod(iosw(ii,ij),jpni)
            ijsw=1+   (iosw(ii,ij))/jpni
            IF(ibne(iisw,ijsw).EQ.1) ibne(iisw,ijsw)=0
        ENDIF
        IF(    ipproc(ii,ij).EQ.-1
     $      .AND.ibnw(ii,ij).EQ.1) THEN
            iinw=1+mod(ionw(ii,ij),jpni)
            ijnw=1+   (ionw(ii,ij))/jpni
            IF(ibse(iinw,ijnw).EQ.1) ibse(iinw,ijnw)=0
        ENDIF
        IF(    ipproc(ii,ij).EQ.-1
     $      .AND.ibse(ii,ij).EQ.1) THEN
            iise=1+mod(iose(ii,ij),jpni)
            ijse=1+   (iose(ii,ij))/jpni
            IF(ibnw(iise,ijse).EQ.1) ibnw(iise,ijse)=0
        ENDIF
      ENDDO
C
C
C 6. Change processor name
C ------------------------
C
      nproc=narea-1
      ii=iin(narea)
      ij=ijn(narea)
      IF (ioso(ii,ij).ge.0 .and. 
     $    ioso(ii,ij).le.(jpni*jpnj-1)) THEN 
          iiso=1+mod(ioso(ii,ij),jpni)
          ijso=1+   (ioso(ii,ij))/jpni
          noso=ipproc(iiso,ijso)
      ENDIF 
      IF (iowe(ii,ij).ge.0 .and. 
     $    iowe(ii,ij).le.(jpni*jpnj-1)) THEN 
          iiwe=1+mod(iowe(ii,ij),jpni)
          ijwe=1+   (iowe(ii,ij))/jpni
          nowe=ipproc(iiwe,ijwe)
      ENDIF 
      IF (ioea(ii,ij).ge.0 .and. 
     $    ioea(ii,ij).le.(jpni*jpnj-1)) THEN 
          iiea=1+mod(ioea(ii,ij),jpni)
          ijea=1+   (ioea(ii,ij))/jpni
          noea=ipproc(iiea,ijea)
      ENDIF
      IF (iono(ii,ij).ge.0 .and. 
     $    iono(ii,ij).le.(jpni*jpnj-1)) THEN 
          iino=1+mod(iono(ii,ij),jpni)
          ijno=1+   (iono(ii,ij))/jpni
          nono=ipproc(iino,ijno)
      ENDIF
      IF (iose(ii,ij).ge.0 .and. 
     $    iose(ii,ij).le.(jpni*jpnj-1)) THEN 
          iise=1+mod(iose(ii,ij),jpni)
          ijse=1+   (iose(ii,ij))/jpni
          npse=ipproc(iise,ijse)
      ENDIF
      IF (iosw(ii,ij).ge.0 .and. 
     $    iosw(ii,ij).le.(jpni*jpnj-1)) THEN 
          iisw=1+mod(iosw(ii,ij),jpni)
          ijsw=1+   (iosw(ii,ij))/jpni
          npsw=ipproc(iisw,ijsw)
      ENDIF
      IF (ione(ii,ij).ge.0 .and. 
     $    ione(ii,ij).le.(jpni*jpnj-1)) THEN 
          iine=1+mod(ione(ii,ij),jpni)
          ijne=1+   (ione(ii,ij))/jpni
          npne=ipproc(iine,ijne)
      ENDIF
      IF (ionw(ii,ij).ge.0 .and. 
     $    ionw(ii,ij).le.(jpni*jpnj-1)) THEN 
          iinw=1+mod(ionw(ii,ij),jpni)
          ijnw=1+   (ionw(ii,ij))/jpni
          npnw=ipproc(iinw,ijnw)
      ENDIF
      nbnw=ibnw(ii,ij)
      nbne=ibne(ii,ij)
      nbsw=ibsw(ii,ij)
      nbse=ibse(ii,ij)
      nlcj=ilcj(ii,ij)  
      nlci=ilci(ii,ij)  
      nldi=ildi(ii,ij)
      nlei=ilei(ii,ij)
      nldj=ildj(ii,ij)
      nlej=ilej(ii,ij)
      nbondi=ibondi(ii,ij)
      nbondj=ibondj(ii,ij)
      nimpp=iimppt(ii,ij)  
      njmpp=ijmppt(ii,ij)  

      DO jproc=1,jpnij
        ii=iin(jproc)
        ij=ijn(jproc)
        nimppt(jproc)=iimppt(ii,ij)  
        njmppt(jproc)=ijmppt(ii,ij)  
        nlcjt(jproc)=ilcj(ii,ij)
        nlcit(jproc)=ilci(ii,ij)
        nldit(jproc)=ildi(ii,ij)
        nleit(jproc)=ilei(ii,ij)
        nldjt(jproc)=ildj(ii,ij)
        nlejt(jproc)=ilej(ii,ij)
      END DO
C
      DO ji=1,jpi
        mindi(ji)=ji+nizoom-1+nimpp-1
      END DO
C
      DO jj=1,jpj
        mindj(jj)=jj+njzoom-1+njmpp-1
      END DO
C
C FETI method
C
      IF ((nperio .EQ. 1 .OR. nperio.eq.4).AND. nsolv .EQ. 3) THEN 
C
C...general case : Earth == infinite tube
C
          nbnw = 1
          npnw = narea
          nbne = 1
          npne = narea
          nbsw = 1
          npsw = (narea-2)
          nbse = 1
          npse = (narea-2)
C
C...real(8) boundary condition
C
          IF (nbondj .EQ. -1 .OR. nbondj .EQ. 2) THEN 
              nbsw = 0
              nbse = 0
          ENDIF 
C
          IF (nbondj .EQ. -1 .OR. nbondj .EQ. 2) THEN 
              nbsw = 0
              nbse = 0
          ENDIF 
C
          IF (nbondj .EQ.  1 .OR. nbondj .EQ. 2) THEN 
              nbnw = 0
              nbne = 0
          ENDIF 
      ENDIF 
C
C  Periodicity : no corner if nbondi = 2 and nperio != 1
C
      IF(lwp)WRITE(numout,*) ' nproc=  ',nproc
      IF(lwp)WRITE(numout,*) ' nowe=   ',nowe
      IF(lwp)WRITE(numout,*) ' noea=   ',noea
      IF(lwp)WRITE(numout,*) ' nono=   ',nono
      IF(lwp)WRITE(numout,*) ' noso=   ',noso
      IF(lwp)WRITE(numout,*) ' nbondi= ',nbondi
      IF(lwp)WRITE(numout,*) ' nbondj= ',nbondj
      IF(lwp)WRITE(numout,*) ' npolj=  ',npolj
      IF(lwp)WRITE(numout,*) ' nperio= ',nperio
      IF(lwp)WRITE(numout,*) ' nlci=   ',nlci
      IF(lwp)WRITE(numout,*) ' nlcj=   ',nlcj
      IF(lwp)WRITE(numout,*) ' nimpp=  ',nimpp
      IF(lwp)WRITE(numout,*) ' njmpp=  ',njmpp
      IF(lwp)WRITE(numout,*) ' nbse=   ',nbse,' npse= ',npse
      IF(lwp)WRITE(numout,*) ' nbsw=   ',nbsw,' npsw= ',npsw
      IF(lwp)WRITE(numout,*) ' nbne=   ',nbne,' npne= ',npne
      IF(lwp)WRITE(numout,*) ' nbnw=   ',nbnw,' npnw= ',npnw
C
      IF((nperio .EQ. 1 .OR. nperio.eq.4).and.jpni.ne.1) THEN
          STOP 'inimpp'
      ENDIF
C
C     7. calendar
C
      iglo(1)=jpiglo
      iglo(2)=jpjglo
      iglo(3)=jpk
      iglo(4)=1
      iloc(1)=nlci
      iloc(2)=nlcj
      iloc(3)=jpk
      iloc(4)=1
      iabsf(1)=nimppt(narea)
      iabsf(2)=njmppt(narea)
      iabsf(3)=1
      iabsf(4)=1
      iabsl(:)=iabsf(:)+iloc(:)-1
      ihals(1)=jpreci
      ihals(2)=jprecj
      ihals(3)=0
      ihals(4)=0
      ihale(1)=jpreci
      ihale(2)=jprecj
      ihale(3)=0
      ihale(4)=0
      if (nbondi .eq. -1 .or. nbondi .eq. 2) ihals(1)=0
      if (nbondi .eq. 1 .or. nbondi .eq. 2) ihale(1)=0
      if (nbondj .eq. -1 .or. nbondj .eq. 2) ihals(2)=0
      if (nbondj .eq. 1 .or. nbondj .eq. 2) ihale(2)=0
      IF(lwp) THEN
          WRITE(numout,*) 'iloc',iloc(1),iloc(2),iloc(3)
          WRITE(numout,*) 'iabsf',iabsf(1),iabsf(2),iabsf(3)
          WRITE(numout,*) 'ihals',ihals(1),ihals(2),ihals(3)
          WRITE(numout,*) 'ihale',ihale(1),ihale(2),ihale(3)
      ENDIF
CCC      CALL calendar_inimpp(jpnij,nproc,iglo,iloc,iabsf,iabsl,ihals,ihale)
C
#  else
C
C      No mpp computation 
C
      nimpp=1
      njmpp=1
      nlci=jpi
      nlcj=jpj
      nldi=1
      nlei=jpi
      nldj=1
      nlej=jpj
      nperio=jperio
      nbondi=2
      nbondj=2
      IF(lwp)WRITE(numout,*) 
      IF(lwp)WRITE(numout,*) ' *** inimpp'
      IF(lwp)WRITE(numout,*) ' nperio= ',nperio
C
#endif
C
C
C
      RETURN
      END
