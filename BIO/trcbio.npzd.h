c $Id: trcbio.npzd.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC      trcbio.npzd.h
CCC      *****************
CCC
CC   defined key : key_trc_npzd
CC   ===========
CC
CC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      COMMON
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC                                depths
CC            /cottrp/          : present and next fields for passive
CC                              : tracers
CC            /comtsk/          : multitasking
CC            /cotbio/          : biological parameters
CC
CC   OUTPUT :
CC   ------
CC      COMMON
CC            /cottrp/ tra      : general tracer trend increased by the
CC                                now horizontal tracer advection trend
CC            /cottbd/ trbio    : now horizontal tracer advection trend
CC                                (IF 'key_trc_diabio' is activated)
CC
CC   WORKSPACE :
CC   ---------
CC      local
CC               zdet,zzoo,zphy,znut              : now concentrations
CC               zlt,zlnut,zlpe                   : limitation terms for phyto
CC                                                 
CC               zflxnp,zflxpn,zflxpz,zflxdz      : fluxes between bio
CC                                                  boxes
CC               zflxpd,zflxzd,zflxdn
CC               zphya,zzooa,znuta,zdeta          : after bio trends
CC               zphimp, zmp, zphimz, zmz         : mortality terms
CC               zppz, zpdz, zpppz, zppdz, zfood  : preferences terms
CC               zfilpz, zfilpd                   : filtration terms
CC
CC   EXTERNAL :                   no
CC   --------
CC
CC   REFERENCES :                 no
CC   ----------
CC
CC   MODIFICATIONS:
CC   --------------
CC       original : 95-02 (M. Levy)
CC                  99-07 (M. Levy) version .h
CC                  99-09 (M. Levy) version with no deep mixing limitation
CC       adaptations : 00-12 (E. Kestenare) 
CC                     assign a parameter to name individual tracers
CC                     01-02 (E. Kestenare)
CC                     introduce jpno3 instead of jpnut 
CC                     01-02 (E. Kestenare) add sediments 
CC----------------------------------------------------------------------


#include "parameter.h"
#include "common.h"

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER kt,ktask

      INTEGER ji,jj,jk,jn
      REAL(8) ztot(jpi)
      REAL(8) zdet,zzoo,zphy,znut,zflxnp,zflxpn,zppz,zpdz,zpppz,
     $    zppdz,zfood,zfilpz,zfildz,zflxpz,zflxdz,zflxzl,zflxzn,zflxpd,
     $    zflxzd,zflxdn,zphya,zzooa,znuta,zdeta,ztra
      REAL(8) zle,zlt,zlnut
CCC 28/4/2004 Paolo aux variables
      REAL(8) db,zb,pb
CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

ccC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C   | --------------|
C   | NPZD MODEL| 
C   | --------------|

C
C vertical slab
C =============
C
      DO 1000 jj = ktask+1,jpjm1,ncpu
C
C 1. biological level
C ===================
C
C       set fbod to 0 for sediments
        DO ji = 2,jpim1
          fbod(ji,jj)=0.
        END DO
 
        DO jk=1,jpkbm1
          DO ji = 2,jpim1
C
          IF(Tmask(ji,jj,jk).eq.1) THEN
C
C 1.1 trophic variables( det, zoo, phy, nut)
C ------------------------------------------
C
C negative trophic variables DO not contribute to the fluxes
C
CCC Paolo 23/4/2004 minimum concentration (NAMELIST)
CCC            zdet = max(0.,trn(ji,jj,jk,jpdet))
CCC            zzoo = max(0.,trn(ji,jj,jk,jpzoo))
CCC            zphy = max(0.,trn(ji,jj,jk,jpphy))
CCC            znut = max(0.,trn(ji,jj,jk,jpno3))

 
            zdet = max(admin,trn(ji,jj,jk,jpdet))
            zzoo = max(azmin,trn(ji,jj,jk,jpzoo))
            zphy = max(apmin,trn(ji,jj,jk,jpphy))
            znut = max(0.,trn(ji,jj,jk,jpno3))
            
            db=admin-min(admin,trn(ji,jj,jk,jpdet));
            zb=azmin-min(azmin,trn(ji,jj,jk,jpzoo));
            pb=apmin-min(apmin,trn(ji,jj,jk,jpphy));
            znut = znut- db -zb -pb;


C
C
C 1.2 Limitations
C -----------------------
C
C 16/6/2004 Paolo i consider also the temperature
           zlt = 1.
            zlt = min(exp( log(2.)*( tn(ji,jj,jk)- 20.)/ 20.),1.)
            zle = 1. - exp( -xpar(ji,jj,jk)/aki/zlt) 
            zle = 1.
            zlnut = znut / ( znut + aknut )
 
C
C
C 1.3 sinks and sources
C ---------------------
C
C
C 4.1 phytoplankton production and exsudation
C
            zflxnp = tmumax * zle * zlt * zlnut  * zphy
            zflxpn = rgamma * zflxnp
C
C 4.2 zoolplankton production
C
C preferences
C
            zppz = rppz
            zpdz = 1. - rppz
            zpppz = ( zppz * zphy ) /
     $          ( ( zppz * zphy + zpdz * zdet ) + 1.e-13 )
            zppdz = ( zpdz * zdet ) /
     $          ( ( zppz * zphy + zpdz * zdet ) + 1.e-13 )
            zfood = zpppz * zphy + zppdz * zdet
C
C filtration
C
            zfilpz = taus * zpppz / (aks + zfood)
            zfildz = taus * zppdz / (aks + zfood)
C
C grazing
C

            zflxpz = zfilpz * zphy * zzoo
            zflxdz = zfildz * zdet * zzoo
C
C 3. fecal pellets production
C
            zflxzl = rpnaz * zflxpz + rdnaz * zflxdz
C
C 4. zooplankton liquide excretion IF zzoo greater THEN eggzoo
C
            zflxzn = tauzn * zzoo * (1. + sign(1., zzoo - eggzoo))/2.
C
C 5. mortality
C
C phytoplankton mortality
C
            zflxpd = tmminp * zphy
C
C
C zooplankton mortality if zzoo greater then eggzoo
C
C
            zflxzd = tmminz * zzoo * (1. + sign(1., zzoo - eggzoo))/2.
C
C
C 6. detritus breakdomn
C
            zflxdn = taudn * zdet
C
C
C 1.4 determination of trends
C ---------------------------
C
C total trend for each biological tracer
C
            zphya = zflxnp - zflxpn - zflxpz - zflxpd
            zzooa = zflxpz + zflxdz - zflxzd - zflxzl - zflxzn
            znuta = zflxpn + zflxzn + zflxdn - zflxnp
            zdeta = zflxpd + zflxzd + zflxzl - zflxdn - zflxdz
C
#if defined key_trc_diabio
            trbio(ji,jj,jk,1) = zflxnp
            trbio(ji,jj,jk,2) = zflxpz
            trbio(ji,jj,jk,3) = zflxpd
            trbio(ji,jj,jk,4) = zflxdz
            trbio(ji,jj,jk,5) = zflxzd + zflxzl
            trbio(ji,jj,jk,6) = zflxzn
            trbio(ji,jj,jk,7) = zflxdn
            trbio(ji,jj,jk,9) = zflxpn + zflxzn + zflxdn
#endif
C
C tracer flux at totox-point added to the general trend
C
CCC Paolo 26/4/2004 added time step rdt and consider  trn ?
CCC Paolo time conversion days->seconds
CCC Here tra is trend matrix
CCC            tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + zdeta
CCC            tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + zzooa
CCC            tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zphya
CCC            tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + znuta

            
             tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + qsr(ji,jj)/24./3600.
C            tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + zdeta/24./3600.
C     $                            *tmask(ji,jj,jk)    
            tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + zzooa/24./3600.
     $                            *tmask(ji,jj,jk)
            tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zphya/24./3600.
     $                            *tmask(ji,jj,jk)
            tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + znuta/24./3600.
     $                            *tmask(ji,jj,jk)
C
           END IF
          END DO
        END DO

C
C 2. under biological level
C =========================
C
        DO jk=jpkb,jpk
C
C 2.1 compute the remineralisation of all quantities towards nutrients
C --------------------------------------------------------------------
C
          DO ji=2,jpim1
            ztot(ji) = 0.
          END DO 
          DO jn=1,jptra
            IF (ctrcnm(jn).NE.'NO3') THEN 
                DO ji=2,jpim1
CCC trn is a concentration tra is a trend ????
                  ztra = remdmp(jk,jn) * trn(ji,jj,jk,jn)

CCC Paolo 3/5/2004 tra is a trend negative values admitted

                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) - ztra
                  ztot(ji) = ztot(ji) + ztra


                END DO 
            ENDIF
          END DO 
          DO jn=1,jptra
            IF (ctrcnm(jn).EQ.'NUT') THEN 
                DO ji=2,jpim1
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztot(ji)
                END DO
#if defined key_trc_diabio
                trbio(ji,jj,jk,1)=ztot(ji)
#endif 
            ENDIF
          END DO
        END DO 
C
#if defined key_trc_diabio
        DO jk=jpkb,jpkm1
          DO ji=2,jpim1
            trbio(ji,jj,jk,9)=tra(ji,jj,jk,jptra)
          END DO
        END DO
#endif
C
C
C END of slab
C ===========
C
 1000 CONTINUE
