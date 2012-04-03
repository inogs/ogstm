C $Id: trchdf.isopycnal.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC        trchdf.isopycnal.h
CCC      **********************
CCC
CCC   defined key: 'key_trahdfiso' or ('key_s_coord' & 'key_trc_hdfgeop')
CCC   ============
CC
CC   Method :
CC   -------
CC         The horizontal component of the lateral diffusive trends is
CC	provided by a 2nd order operator rotated along neural or geopo-
CC	tential surfaces to which an eddy induced advection can be added
CC	It is computed using before fields (forward in time) and isopyc-
CC	nal or geopotential slopes computed in routine hdfslp.
CC
CC	horizontal fluxes associated with the rotated lateral mixing:
CC         zftu = (ahtt+ahtrb0) e2u*e3u/e1u di[ trb ]
CC	         - ahtt       e2u*uslp    dk[ mi(mk(trb)) ]
CC         zftv = (ahtt+ahtrb0) e1v*e3v/e2v dj[ trb ]
CC		 - ahtt       e2u*vslp    dk[ mj(mk(trb)) ]
CC	add horizontal Eddy Induced advective fluxes ('key_trahdfeiv'):
CC         zftu = zftu - dk-1[ ahtt e2u mi(wslpi) ] mi( trb ) 
CC         zftv = zftv - dk-1[ ahtt e1v mj(wslpj) ] mj( trb ) 
CC	take the horizontal divergence of the fluxes:
CC         difftr = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
CC      Add this trend to the general trend (tra):
CC         tra = tra + difftr
CC
CC      'key_trc_diatrd' defined: the trend is saved for diagnostics.
CC
CC      macro-tasked on tracer slab (jn-loop).
CC
CC   Input :
CC   ------
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC            /comcoh/          : scale factors
CC            /comask/          : masks
CC            /cottrc/          : previous fields (before) and next
CC                                fields (after) for passive tracer
CC            /comiso/          : isopycnal slopes
CC            /comtsk/          : multitasking
CC
CC   Output :
CC   -------
CC      common
CC            /cottrc/ tra      : passive tracer trend increased by the
CC                                before horizontal component of the
CC				  lateral diffusive trend 
CC            /cottrd/ trtrd    : now horizontal tracer diffusive trend
CC                                ('key_trc_diatrd' defined)
CC
CC   Modifications :
CC   --------------
CC       original :  94-08 (G. Madec, M. Imbard)
CC       addition :  97-05 (G. Madec) split into trchdf and trczdf
CC       addition :  98-03 (L. Bopp, MA Foujols) passive tracer generalisation
CC      modification : 00-10 (MA Foujols E Kestenare) USE passive tracer
CC                            coefficient
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

      USE myalloc

        IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ktask, kt

#if defined key_passivetrc

      INTEGER ji, jj, jk, jn
C
      REAL(8) zabe1, zabe2, zcof1, zcof2, zmsku, zmskv, zbtr, zta, zsa
      REAL(8) zcg1, zcg2, zuwk, zvwk, zuwk1, zvwk1
      REAL(8) ztagu, ztagv
      REAL(8) zftu(jpi,jpj), zftv(jpi,jpj), zdkt(jpi,jpj), zdk1t(jpi,jpj)
      REAL(8) zftug(jpi,jpj), zftvg(jpi,jpj)
CC----------------------------------------------------------------------
CC statement functions
CC ===================

!                          #include "stafun.h"
!                          #include "stafun.passivetrc.h"

CCC---------------------------------------------------------------------
CCC  OPA8.1, LODYC (1998)
CCC---------------------------------------------------------------------
C
      ztagu = 0.e0
      ztagv = 0.e0
C
C Passive tracer slab
C ==================
C
      DO 1000 jn = ktask, jptra, ncpu
C
C Horizontal slab
C ===============
C
        DO jk = 1, jpkm1
C
C
C 1. Vertical tracer gradient at level jk and jk+1
C ------------------------------------------------
C surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
C
          DO jj = 1, jpj
            DO ji = 1, jpi
              zdk1t(ji,jj)=(trb(ji,jj,jk,jn)-trb(ji,jj,jk+1,jn))
     $            *tmask(ji,jj,jk+1)
            END DO
          END DO
          IF (jk.EQ.1) THEN
              DO jj = 1, jpj
                DO ji = 1, jpi
                  zdkt(ji,jj) = zdk1t(ji,jj)
                END DO
              END DO
          ELSE
              DO jj = 1, jpj
                DO ji = 1, jpi
                  zdkt(ji,jj)=(trb(ji,jj,jk-1,jn)-trb(ji,jj,jk,jn))
     $                *tmask(ji,jj,jk)
                END DO
              END DO
          ENDIF
C
C
C 2. Horizontal fluxes
C --------------------
C
          DO jj = 1 , jpjm1
            DO ji = 1 , jpim1
              zabe1 = ( fsahtru(ji,jj,jk) + ahtrb0 )
     $            * e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
              zabe2 = ( fsahtrv(ji,jj,jk) + ahtrb0 )
     $            * e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)
C
              zmsku=1./max( tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)
     $            +tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1. )
              zmskv=1./max( tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)
     $            +tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1. )
C
              zcof1= -fsahtru(ji,jj,jk)*e2u(ji,jj)*uslp(ji,jj,jk)*zmsku
              zcof2= -fsahtrv(ji,jj,jk)*e1v(ji,jj)*vslp(ji,jj,jk)*zmskv
C
              zftu(ji,jj)= umask(ji,jj,jk) *
     $            (  zabe1 *( trb(ji+1,jj,jk,jn) - trb(ji,jj,jk,jn) )
     $            +zcof1 *( zdkt (ji+1,jj) + zdk1t(ji,jj)
     $            +zdk1t(ji+1,jj) + zdkt (ji,jj) )  )
              zftv(ji,jj)= vmask(ji,jj,jk) *
     $            (  zabe2 *( trb(ji,jj+1,jk,jn) - trb(ji,jj,jk,jn) )
     $            +zcof2 *( zdkt (ji,jj+1) + zdk1t(ji,jj)
     $            +zdk1t(ji,jj+1) + zdkt (ji,jj) )  )
            END DO
          END DO
C
#ifdef key_trahdfeiv
C
C II.3. Eddy induced horizontal advective fluxes
C ----------------------------------------------
C
          DO jj=1,jpjm1
            DO ji=1,jpim1
#if defined key_off_tra 
              zcg1= -ugm(ji,jj,jk)*e2u(ji,jj)*fse3u(ji,jj,jk)*0.5
              zcg2= -vgm(ji,jj,jk)*e1v(ji,jj)*fse3v(ji,jj,jk)*0.5
#else
              zuwk = ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj,jk) )
     $            * fsaeitru(ji,jj,jk  ) * umask(ji,jj,jk)
              zuwk1= ( wslpi(ji,jj,jk+1) + wslpi(ji+1,jj,jk+1) )
     $            * fsaeitru(ji,jj,jk+1) * umask(ji,jj,jk+1)
              zvwk = ( wslpj(ji,jj,jk  ) + wslpj(ji,jj+1,jk) )
     $            * fsaeitrv(ji,jj,jk  ) * vmask(ji,jj,jk)
              zvwk1= ( wslpj(ji,jj,jk+1) + wslpj(ji,jj+1,jk+1) )
     $            * fsaeitrv(ji,jj,jk+1) * vmask(ji,jj,jk+1)
C
              zcg1= -0.25 * e2u(ji,jj) * umask(ji,jj,jk) * ( zuwk-zuwk1
     $            )
              zcg2= -0.25 * e1v(ji,jj) * vmask(ji,jj,jk) * ( zvwk-zvwk1
     $            )
#endif
C
              zftug(ji,jj) = zcg1 * ( trb(ji+1,jj,jk,jn) + trb(ji,jj,jk
     $            ,jn) ) 
              zftvg(ji,jj) = zcg2 * ( trb(ji,jj+1,jk,jn) + trb(ji,jj,jk
     $            ,jn) ) 
C
              zftu(ji,jj) = zftu(ji,jj) + zftug(ji,jj)
              zftv(ji,jj) = zftv(ji,jj) + zftvg(ji,jj)

#    if defined key_diaeiv
              ugm(ji,jj,jk) = ugm(ji,jj,jk) -
     $            2.* zcg1/ (e2u(ji,jj)*fse3u(ji,jj,jk))
              vgm(ji,jj,jk) = vgm(ji,jj,jk) -
     $            2.* zcg2/ (e1v(ji,jj)*fse3v(ji,jj,jk))
#    endif
            END DO
          END DO
C
#  else
C
C II.3. NO  Eddy induced advective fluxes
C ------==-------------------------------
C
#endif
C
C II.4 Second derivative (divergence) and add to the general trend
C ----------------------------------------------------------------
C
          DO jj = 2 , jpjm1
            DO ji = 2 , jpim1
              zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
              zta =( zftu (ji,jj) - zftu (ji-1,jj  ) +
     $            zftv (ji,jj) - zftv (ji  ,jj-1)  ) * zbtr
              tra (ji,jj,jk,jn) = tra (ji,jj,jk,jn) + zta
#if defined key_trc_diatrd
#    if defined key_trahdfeiv
              ztagu = ( zftug(ji,jj) - zftug(ji-1,jj  ) ) * zbtr
              ztagv = ( zftvg(ji,jj) - zftvg(ji  ,jj-1) ) * zbtr
              trtrd(ji,jj,jk,jn,7) = ztagu
              trtrd(ji,jj,jk,jn,8) = ztagv
#    endif
              trtrd(ji,jj,jk,jn,4)=(zftu(ji,jj)-zftu(ji-1,jj))*zbtr
     $            - ztagu
              trtrd(ji,jj,jk,jn,5)=(zftv(ji,jj)-zftv(ji,jj-1))*zbtr
     $            - ztagv
#endif
            END DO
          END DO
        END DO
C
C
C End of task
C ===========
C
 1000 CONTINUE

#else
C
C no passive tracers
C
#endif
