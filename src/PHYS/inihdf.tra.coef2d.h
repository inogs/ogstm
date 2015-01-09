CCC
CCC             inihdf.tra.coef2d.h
CCC           ***********************
CCC
CCC   defined key :  'key_trahdfcoef2d'
CCC   ============
CCC
CC      2D eddy diffusivity coefficients ( longitude, latitude )
CC      --------------------------------
CC
CC       biharmonic operator    : ahtt = defined at T-level
CC                                ahtu,ahtv,ahtw never used
CC
CC       harmonic operator (ahtt never used)
CC       iso-model level   : ahtu, ahtv defined at u-, v-points
CC         isopycnal         : ahtu, ahtv, ahtw defined at u-, v-, w-pts
CC         or geopotential   
CC
CC     eddy induced velocity
CC       always harmonic   : aeiu, aeiv, aeiw defined at u-, v-, w-pts


      IF(lwp)WRITE(numout,*)
#ifdef key_trahdfeiv
      IF(lwp)WRITE(numout,*) ' inihdf: 2D eddy diffusivity and eddy',
     $                       ' induced velocity coefficients'
#  else
      IF(lwp)WRITE(numout,*) ' inihdf: 2D eddy diffusivity coefficient'
#endif
      IF(lwp)WRITE(numout,*) ' ======  --'
      IF(lwp)WRITE(numout,*)
C
#if defined key_trahdfbilap
C
C biharmonic operator :
C ====================
C
C ... set ahtt at T-point (here no space variation)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ahtt(ji,jj) = aht0
        END DO
      END DO
C
C ... Lateral boundary conditions on ( ahtt )
#ifdef key_mpp
C   ... Mpp: export boundary values to neighbouring processors
      CALL mpplnk2( ahtt, 1, 1 )
# else
C   ... mono- or macro-tasking: T-point, >0, 2D, no-slab
      CALL lbc( ahtt, 1, 1, 1, 1, 1, 1 )
#endif
C
C ... Control print
      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'inihdf: 2D ahtt array'
          CALL prihre(ahtt,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
      ENDIF
C
#  else
C
C harmonic operator :
C ==================
C
C ... set ahtu = ahtv at u- and v-points, and ahtw at w-point (idem T-point)
C     (here no space variation)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ahtu(ji,jj) = aht0
          ahtv(ji,jj) = aht0
          ahtw(ji,jj) = aht0
        END DO
      END DO
C
C ... Lateral boundary conditions on ( ahtu, ahtv, ahtw )
#  ifdef key_mpp
C   ... Mpp: export boundary values to neighbouring processors
      CALL mpplnk2( ahtu, 2, 1 )
      CALL mpplnk2( ahtv, 3, 1 )
      CALL mpplnk2( ahtw, 1, 1 )
#   else
C   ... mono- or macro-tasking: U- and V-point, >0, 2D, no-slab
      CALL lbc( ahtu, 2, 1, 1, 1, 1, 1 )
      CALL lbc( ahtv, 3, 1, 1, 1, 1, 1 )
      CALL lbc( ahtw, 1, 1, 1, 1, 1, 1 )
#  endif
C
C ... Control print
      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'inihdf: ahtu array'
          CALL prihre(ahtu,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
          WRITE(numout,*)
          WRITE(numout,*) 'inihdf: ahtv array'
          CALL prihre(ahtv,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
          WRITE(numout,*)
          WRITE(numout,*) 'inihdf: ahtw array'
          CALL prihre(ahtw,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
      ENDIF
C
#    ifdef key_trahdfeiv
C ... set aeiu = aeiv at u- and v-points, and aeiw at w-point (idem T-point)
C     (here no space variation)
      DO jj = 1, jpj
        DO ji = 1, jpi
          aeiu(ji,jj) = aeiv0
          aeiv(ji,jj) = aeiv0
          aeiw(ji,jj) = aeiv0
        END DO
      END DO
C
C ... Lateral boundary conditions on ( aeiu, aeiv, aeiw )
#  ifdef key_mpp
C   ... Mpp: export boundary values to neighbouring processors
      CALL mpplnk2( aeiu, 2, 1 )
      CALL mpplnk2( aeiv, 3, 1 )
      CALL mpplnk2( aeiw, 1, 1)
#   else
C   ... mono- or macro-tasking: U- and V-point, >0, 2D, no-slab
      CALL lbc( aeiu, 2, 1, 1, 1, 1, 1 )
      CALL lbc( aeiv, 4, 1, 1, 1, 1, 1 )
      CALL lbc( aeiw, 1, 1, 1, 1, 1, 1 )
#  endif
C
C ... Control print
      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'inihdf: aeiu array'
          CALL prihre(aeiu,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
          WRITE(numout,*)
          WRITE(numout,*) 'inihdf: aeiv array'
          CALL prihre(aeiv,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
          WRITE(numout,*)
          WRITE(numout,*) 'inihdf: aeiw array'
          CALL prihre(aeiw,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
      ENDIF
C
#    endif
C
#endif
