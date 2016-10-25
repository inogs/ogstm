      SUBROUTINE inimpp
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE INIMPP
!!!                     ******************
!!!
!!!  PURPOSE :
!!!  ---------
!!!     Lay out the global domain over processors.
!!!
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
CC                    nbondi    : mark for "east-west local boundary"
CC                    nbondj    : mark for "north-south local boundary"
CC                    nproc     : number for local processor
CC                    noea      : number for local neighboring processor
CC                    nowe      : number for local neighboring processor
CC                    noso      : number for local neighboring processor
CC                    nono      : number for local neighboring processor
CC
CC   WORKSPACE :
CC   ---------
CC      local
CC              jj,ji, jn, jproc
CC              ij,ii, ifreq, il1, il2, iresti, irestj, zidom, zjdom
CC              nlcit(), nlcjt(), nimppt(), njmppt()
CC
CC   MODIFICATIONS:
CC   --------------
CC       original  : 94-11 (M. Guyon)
CC       additions : 95-04 (J. Escobar, M. Imbard)
CC       additions : 98-02 (M. Guyon) FETI method
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
CC----------------------------------------------------------------------

      USE calendar
      USE myalloc
        IMPLICIT NONE


      INTEGER nbne, nbnw, nbse, nbsw, npne, npnw, npse, npsw
      INTEGER, PARAMETER :: nsolv = 3

#ifdef key_mpp_mpi
      INTEGER jj,ji,jn,jproc
      INTEGER ij,ii
      INTEGER iresti,irestj,ijm1,imil
      INTEGER iimppt(jpnj,jpni),ijmppt(jpnj,jpni)
      INTEGER ibonit(jpnij),ibonjt(jpnij)
      INTEGER iglo(4),iloc(4),iabsf(4),iabsl(4),ihals(4),ihale(4)
      REAL(8)    zidom ,zjdom

C
C  1. Dimension arrays for subdomains
C -----------------------------------
C
C  Computation of local domain sizes ilcit() ilcjt()
C  These dimensions depend on global sizes jpnj,jpni and jpjglo,jpiglo
C  The subdomains are squares leeser than or equal to the global
C  dimensions divided by the number of processors minus the overlap
C  array.
C
      nreci=2*jpreci
      nrecj=2*jprecj

           iresti = 0
           irestj = 0
C

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
        zidom = zidom + ilcit(1,ji) - nreci
      END DO
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)" sum ilcit(i,1)=",zidom," jpiglo=",jpiglo
C
      zjdom=nrecj
      DO jj=1,jpnj
        zjdom = zjdom + ilcjt(jj,1) - nrecj
      END DO
      IF(lwp)WRITE(numout,*)" sum ilcit(1,j)=",zjdom," jpjglo=",jpjglo
      IF(lwp)WRITE(numout,*)
C
C
C  2. Index arrays for subdomains
C -------------------------------
C
      DO jj=1,jpnj
        DO ji=1,jpni
          iimppt(jj,ji)=1
          ijmppt(jj,ji)=1
        END DO
      END DO
C
      IF(jpni.gt.1)THEN
          DO jj=1,jpnj
            DO ji=2,jpni
              iimppt(jj,ji)=iimppt(jj,ji-1)+ilcit(jj,ji-1)-nreci
            END DO
          END DO
      ENDIF
C
      IF(jpnj.gt.1)THEN
          DO jj=2,jpnj
            DO ji=1,jpni
              ijmppt(jj,ji)=ijmppt(jj,ji-1)+ilcjt(jj,ji-1)-nrecj
            END DO
          END DO
      ENDIF
C
C 3. Subdomain description
C ------------------------
C
      DO jproc=1,jpnij
        ii=1+mod(jproc-1,jpni)
        ij=1+   (jproc-1)/jpni
        nimppt(jproc)=iimppt(ij,ii)
        njmppt(jproc)=ijmppt(ij,ii)
        nlcit(jproc)=ilcit(ij,ii)     
        nlci=nlcit(jproc)     
        nlcjt(jproc)=ilcjt(ij,ii)     
        nlcj=nlcjt(jproc)
        nbondj=-1
        IF(jproc.gt.jpni) nbondj=0
        IF(jproc.gt.(jpnj-1)*jpni) nbondj=1
        IF(jpnj.eq.1) nbondj=2
        ibonjt(jproc) = nbondj
C
        nbondi=0
        IF(mod(jproc,jpni).eq.1) nbondi=-1
        IF(mod(jproc,jpni).eq.0) nbondi= 1
        IF(jpni.eq.1) nbondi= 2
        ibonit(jproc) = nbondi
C
        nldi= 1  +jpreci
        nlei=nlci-jpreci
        IF(nbondi.eq.-1.or.nbondi.eq.2) nldi=1
        IF(nbondi.eq. 1.or.nbondi.eq.2) nlei=nlci
        nldj= 1  +jprecj
        nlej=nlcj-jprecj
        IF(nbondj.eq.-1.or.nbondj.eq.2) nldj=1
        IF(nbondj.eq. 1.or.nbondj.eq.2) nlej=nlcj
        nldit(jproc)=nldi
        nleit(jproc)=nlei
        nldjt(jproc)=nldj
        nlejt(jproc)=nlej
      END DO
C
C
C 4. From global to local
C -----------------------
C
      nperio=0
      IF(jperio.eq.2.and.nbondj.eq.-1) nperio=2
C
C
C 5. Subdomain neighbours
C ----------------------
C
      nproc=narea-1
      noso=nproc-jpni
      nowe=nproc-1
      noea=nproc+1
      nono=nproc+jpni
C great neighbours
      npnw=nono-1
      npne=nono+1
      npsw=noso-1
      npse=noso+1
      nbsw=1
      nbnw=1
      IF(mod(nproc,jpni) .EQ. 0) THEN
          nbsw=0
          nbnw=0
      ENDIF
      nbse=1
      nbne=1
      IF(mod(nproc,jpni) .EQ. jpni-1) THEN
          nbse=0
          nbne=0
      ENDIF
      IF(nproc .LT. jpni) THEN
          nbsw=0
          nbse=0
      ENDIF
      IF(nproc .GE. (jpnj-1)*jpni) THEN
          nbnw=0
          nbne=0
      ENDIF
      nlcj=nlcjt(narea)  
      nlci=nlcit(narea)  
      nldi=nldit(narea)
      nlei=nleit(narea)
      nldj=nldjt(narea)
      nlej=nlejt(narea)
      nbondi=ibonit(narea)
      nbondj=ibonjt(narea)
      nimpp=nimppt(narea)  
      njmpp=njmppt(narea)  

C
      DO ji=1,jpi
        mindi(ji)=ji+nizoom-1+nimpp-1
      END DO
C
      DO jj=1,jpj
        mindj(jj)=jj+njzoom-1+njmpp-1
      END DO
C
C w a r n i n g  narea (zone) .ne. nproc (processors)!
C
      IF(jperio.eq.1.or.jperio.eq.4) THEN
          IF(jpni.eq.1)THEN
              nbondi=2
              nperio=1
          ELSE
              nbondi=0
          ENDIF
          IF(mod(narea,jpni).eq.0) THEN
              noea=nproc-(jpni-1)
              npne=npne-jpni
              npse=npse-jpni
          ENDIF
          IF(mod(narea,jpni).eq.1) THEN
              nowe=nproc+(jpni-1)
              npnw=npnw+jpni
              npsw=npsw+jpni
          ENDIF
          nbsw=1
          nbnw=1
          nbse=1
          nbne=1
          IF(nproc .LT. jpni) THEN
              nbsw=0
              nbse=0
          ENDIF
          IF(nproc .GE. (jpnj-1)*jpni) THEN
              nbnw=0
              nbne=0
          ENDIF
      ENDIF
      npolj=0
      IF(jperio.eq.3 .OR. jperio.eq.4) THEN
          ijm1=jpni*(jpnj-1)
          imil=ijm1+(jpni+1)/2
          IF(narea.gt.ijm1) npolj=3
          IF(mod(jpni,2).eq.1.and.narea.eq.imil) npolj=4
          IF(npolj.eq.3) nono=jpni*jpnj-narea+ijm1
      ENDIF
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
          IF (nbondj .EQ. 1  .OR. nbondj .EQ. 2) THEN 
              nbnw = 0
              nbne = 0
          ENDIF 
      ENDIF 
C
C Periodicity : no corner if nbondi = 2 and nperio != 1


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
      IF((nperio.eq.1.or.nperio.eq.4).and.jpni.ne.1)THEN
          STOP 'inimpp'
      ENDIF


C     6. calendar
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
