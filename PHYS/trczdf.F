C $Id: trczdf.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE trczdf ( ktask, kt )
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trczdf
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Compute the trend due to the vertical tracer diffusion inclu-
CCC     ding the vertical component of lateral mixing (only for second
CCC     order operator, for fourth order it is already computed and
CCC     add to the general trend in trchdf.F) and the surface forcing
CCC     and add it to the general trend of the tracer equations.
CC
CC   Method :
CC   -------
CC         The vertical component of the lateral diffusive trends is
CC      provided by a 2nd order operator rotated along neural or geopo-
CC      tential surfaces to which an eddy induced advection can be added
CC      It is computed using before fields (forward in time) and isopyc-
CC      nal or geopotential slopes computed in routine hdfslp.
CC
CC	First part: vertical trends associated with the lateral mixing
CC	==========  (excluding the vertical flux proportional to dk[tr] )
CC	            ('key_trahdfiso' or 'key_trahdfgeop')
CC	vertical fluxes associated with the rotated lateral mixing:
CC	   zftw =-ahtt {  e2t*wslpi di[ mi(mk(trb)) ]
CC	               + e1t*wslpj dj[ mj(mk(trb)) ]  }
CC	save avt coef. resulting from vertical physics alone in zavt:
CC	   zavt = avt
CC	update and save in zavt the vertical eddy viscosity coefficient:
CC	   avt = avt + wslpi^2+wslj^2
CC	add vertical Eddy Induced advective fluxes ('key_trahdfeiv'):
CC	   zftw = zftw + { di[ahtt e2u mi(wslpi)]
CC			  +dj[ahtt e1v mj(wslpj)] } mk(trb)
CC	take the horizontal divergence of the fluxes:
CC         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
CC      Add this trend to the general trend (tra):
CC         tra = tra + difft
CC
CC	Second part: vertical trend associated with the vertical physics
CC	===========  (including the vertical flux proportional to dk[tr]
CC		      associated with the lateral mixing, through the
CC		      update of avt)
CC	The vertical diffusion of passive tracers (tr) is given by:
CC             difft = dz( avt dz(tr) ) = 1/e3t dk+1( avt/e3w dk(tr) )
CC         'key_zdfexplicit' defined: forward  time scheme (tr=trb), using
CC	   a time splitting technique (NOT YET IMPLEMENTED).
CC          default option          : backward time scheme, tr=tra
CC                                  : use ndttrc or 2 for time stepping
CC      Bottom boundary conditions:
CC         no flux on passive tracers: applied through the masked field avt
CC      Add this trend to the general trend tra :
CC         tra = tra + dz( avt dz(tr) )
CC
CC	Third part: recover avt resulting from the vertical physics
CC	==========  alone, as it can be used to compute the turbocline
CC		    depth in diamld ('key_trahdfiso' or 'key_trahdfgeop')
CC	   avt = zavt
CC
CC      'key_trc_diatrd' defined: trend saved for futher diagnostics.
CC
CC      macro-tasked on tracer slab (jn-loop)
CC
CC   Input :
CC   ------
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC            /comcoo/          : scale factors
CC            /comask/          : masks, bathymetry
CC            /cottrc/          : previous fields (before) and next
CC                                fields (after)  for passive tracer
CC
CC   Output :
CC   -------
CC      common
CC            /cottrc/ tra      : passive tracer (tr) trend increased
CC                                by the after vertical diffusion trend
CC            /cottrd/ trtrd    : after vertical tracer diffusion trend
CC                                (IF 'key_trc_diatrd' is defined)

CC parameters and commons
CC ======================


       USE myalloc
       USE myalloc_mpp
       USE ZDF_mem

        IMPLICIT NONE

CC---------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ktask, kt

#if defined key_passivetrc

      INTEGER ji, jj, jk, jn, jv
! omp variables
      INTEGER :: mytid, ntids

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

C
      REAL(8) ztavg, zdt
C
      INTEGER ikst, ikenm2, ikstp1
C
C only IMPLICIT scheme
C
      REAL(8) z2dtt, z2dtf, ztra, zcoef
#if ! defined key_trahdfbilap
#    if defined key_trahdfiso || defined key_trahdfgeop
      REAL(8) zcoef0, zcoef3, zcoef4, zcoeg3
      REAL(8) zbtr, zmku, zmkv, zuwki, zvwki, zuwk, zvwk, ztav
      REAL(8) ztfw(jpi,jpk), zdit(jpi,jpk), zdjt(jpi,jpk), zdj1t(jpi,jpk)
      REAL(8) zavt(jpi,jpj,jpk), ztfwg(jpi,jpk)
#    endif
#endif
CC---------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"


      trczdfparttime = MPI_WTIME()
#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif


      IF (dimen_jvzdf .EQ. 0) THEN
         DO jj = 2,jpjm1
           DO  ji = 2,jpim1
              IF(tmask(ji,jj,1) .NE. 0) THEN
                 dimen_jvzdf = dimen_jvzdf + 1
                 jarr_zdf(1,dimen_jvzdf) = ji
                 jarr_zdf(2,dimen_jvzdf) = jj
              ENDIF
            END DO
         END DO
      ENDIF

C
#    if ! defined key_trahdfbilap
#        if defined key_trahdfiso || defined key_trahdfgeop
C
C

C
C 0. SAVE avt in zavt 
C ----------------------
C
          DO jk = 2, jpkm1
            DO jj = 2, jpjm1
              DO ji = 2, jpim1
                zavt(ji,jj,jk) = avt(ji,jj,jk)
              END DO 
            END DO 
          END DO 
C

C
C
#        endif
#    endif
C
C passive tracer slab
C ===================
C

      TRACER_LOOP: DO  jn = 1, jptra, ntids
C
C 0. Local constant initialization
C --------------------------------
C
!$omp    parallel default(none) private(jv,ji,jj,jk,mytid,ztavg,zdt,z2dtt,ikst,ikstp1,ikenm2)
!$omp&                           shared(jn,dimen_jvzdf,jarr_zdf,ndttrc,rdttra,jpk,jpkm1,avt,zwi,
!$omp&                                   zws,zwd,zwy,trb,tra,zwt,zwz,zwx,tmask,e3t,e3w)

#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
      if( mytid + jn <= jptra ) then

        ztavg = 0.e0
C
C
C vertical slab
C =============
C
        DO jv = 1, dimen_jvzdf

           ji = jarr_zdf(1,jv)
           jj = jarr_zdf(2,jv)
C
C
C I. Vertical trends associated with lateral mixing
C -------------------------------------------------
C    (excluding the vertical flux proportional to dk[t] )
C
#if ! defined key_trahdfbilap
#    if defined key_trahdfiso || defined key_trahdfgeop
#    include "trczdf.isopycnal.h"
#      else
C    No vertical component associated with isopycnal diffussion
C    and with geopotential diffusion in s-coordinates
#    endif
#endif
C
C
C II. Vertical trend associated with the vertical physics
C -------------------------------------------------------
C     (including the vertical flux proportional to dk[t] associated
C      with the lateral mixing, through the avt update)
C
CCC
CCC        trczdf.implicit.h
CCC      *********************
CCC
CCC   define key :  no (default option)
CCC   ==========
CC

C
C
C II. dk[ avt dk[ (tr) ] ] diffusive trends
C ==========================================
C
C
C 0. initialization
C -----------------
C
#if defined key_trc_arakawa || defined key_trc_tvd || defined key_trc_upstream
      zdt= 2.*float(ndttrc)
#else
      zdt= float(ndttrc)
#endif
C
C II.0 Matrix construction
C ------------------------
C
C Diagonal, inferior, superior
C (including the bottom boundary condition via avt masked)
        DO jk = 1, jpkm1
          z2dtt = zdt * rdttra(jk)
C   ... Euler time stepping when starting from rest


            zwi(jk,mytid+1) = - z2dtt * avt(ji,jj,jk  )
     $                        /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
            zws(jk,mytid+1) = - z2dtt * avt(ji,jj,jk+1)
     $                        /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
            zwd(jk,mytid+1) = 1. - zwi(jk,mytid+1) - zws(jk,mytid+1)

        END DO
C
C Surface boundary conditions

          zwi(1,mytid+1) = 0.e0
          zwd(1,mytid+1) = 1. - zws(1,mytid+1)

C
C
C II.1. Vertical diffusion on tr
C ------------------------------
C
C Second member construction
        DO jk = 1, jpkm1
          z2dtt = zdt * rdttra(jk)
C   ... Euler time stepping when starting from rest


CCC Paolo tra is tracer or tracer trend?
CCC 31 10 2005 F79 We introduce a module for time advancing with diagnostic
CCC for false negatives
            zwy(jk,mytid+1) = trb(ji,jj,jk,jn+mytid) + z2dtt * tra(ji,jj,jk,jn+mytid)
CCC 31 10 2005 F79 we pass zwy to the F neg diagnostic
CCC            tra(ji,jj,jk,jn) = zwy(jk) * tmask(ji,jj,jk)

        END DO

CCC 31 10 2005 F79 we pass zwy to the F neg diagnostic
C       IF(jn .EQ. jptra) THEN
C        IF(jj .EQ. jpjm1) THEN
C         CALL count_FN(kt)
C        ENDIF
C       ENDIF
C
C Matrix inversion from the first level
        ikst = 1
C
C It is possible to change the maximum bottom value
C for matrix inversion (see Marina Levy s thesis)
C That will be done with some more documentation
C
CCC
CCC        ZDF.MATRIXSOLVER
CCC      ********************
CCC
CC Matrix inversion
CC----------------------------------------------------------------------
C   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
C
C        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
C        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
C        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
C        (        ...               )( ...  ) ( ...  )
C        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
C
C   m is decomposed in the product of an upper and lower triangular
C   matrix
C   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
C   The second member is in 2d array zwy
C   The solution is in 2d array zwx
C   The 2d arry zwt and zwz are work space arrays
C
C   N.B. the starting vertical index (ikst) is equal to 1 except for
C   the resolution of tke matrix where surface tke value is prescribed
C   so that ikstrt=2.

C
        ikstp1=ikst+1
        ikenm2=jpk-2

        zwt(ikst,mytid+1)=zwd(ikst,mytid+1)

        DO jk=ikstp1,jpkm1

            zwt(jk,mytid+1)=zwd(jk,mytid+1)-zwi(jk,mytid+1)*zws(jk-1,mytid+1)
     $                             /zwt(jk-1,mytid+1)
        END DO

          zwz(ikst,mytid+1)=zwy(ikst,mytid+1)

        DO jk=ikstp1,jpkm1

            zwz(jk,mytid+1)=zwy(jk,mytid+1)-zwi(jk,mytid+1)/zwt(jk-1,mytid+1)
     $                             *zwz(jk-1,mytid+1)
        END DO

        zwx(jpkm1,mytid+1)=zwz(jpkm1,mytid+1)/zwt(jpkm1,mytid+1)

        DO jk=ikenm2,ikst,-1

            zwx(jk,mytid+1)=( zwz(jk,mytid+1)-zws(jk,mytid+1)*zwx(jk+1,mytid+1) )
     $                 /zwt(jk,mytid+1)

        END DO

C
#if defined key_trc_diatrd
C Compute and save the vertical diffusive tracer trends
         DO jk = 1, jpkm1
           z2dtt = zdt * rdttra(jk)
C   ... Euler time stepping when starting from rest


            ztra = ( zwx(jk,mytid+1) - trb(ji,jj,jk,jn+mytid) ) / z2dtt
            trtrd(ji,jj,jk,jn+mytid,6) = ztra - tra(ji,jj,jk,jn+mytid)
#    if defined key_trc_hdfiso
     $                       + trtrd(ji,jj,jk,jn+mytid,6)
#    endif

         END DO
#endif
C
C Save the masked passive tracer after in tra
C (c a u t i o n:  tracer not its trend, Leap-frog scheme done
C                  it will not be done in trcnxt)
         DO jk = 1, jpkm1

            tra(ji,jj,jk,jn+mytid) = zwx(jk,mytid+1) * tmask(ji,jj,jk)

         END DO
C
C END of jv loop
C
        END DO 
#  if defined key_trc_diatrd
#         ifdef key_mpp
        CALL mpplnk( trtrd(1,1,1,jn,6), 1, 1 )
#         else      
        CALL lbc( trtrd(1,1,1,jn,6), 1, 1, 1, 1, jpk, 1 )
#         endif
#  endif
C
C End of task
C ===========
C
C
      end if

!$omp    end parallel
       END DO TRACER_LOOP


#    if ! defined key_trahdfbilap
#        if defined key_trahdfiso || defined key_trahdfgeop
C
C
C
C III. recover the avt resulting from vertical physics only
C ---------------------------------------------------------
C
          DO jk = 2, jpkm1
            DO jj = 2, jpjm1
              DO ji = 2, jpim1
                avt(ji,jj,jk) = zavt(ji,jj,jk)
              END DO 
            END DO 
          END DO 
C


#        endif
#    endif
C
C
CCC 10 11 2004  F79 cronometer-stop

       trczdfparttime = MPI_WTIME() - trczdfparttime
       trczdftottime = trczdftottime + trczdfparttime

CC-CC       write(*,*) "F79T:trczdfparttime", trczdfparttime
CC-CC       write(*,*) "F79T:trczdftottime", trczdftottime
CCC
C
#else
C
C no passive tracers
C
#endif
C
      RETURN
      END
