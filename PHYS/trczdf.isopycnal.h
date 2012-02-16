C $Id: trczdf.isopycnal.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC        trczdf.isopycnal.h
CCC      **********************
CCC
CCC   defined key : 'key_trahdfiso' or ('key_s_coord' & 'key_trahdfgeop')
CCC   ============
CCC
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (1997)
CCC---------------------------------------------------------------------
C
C
C I. vertical trends associated with the lateral mixing
C =====================================================
C  (excluding the vertical flux proportional to dk[t]
C
C
C I.1 horizontal tracer gradient
C ------------------------------
C
        DO jk = 1, jpkm1
          DO ji = 1, jpim1
C   ... i-gradient of passive tracer at jj
            zdit (ji,jk)=(trb(ji+1,jj,jk,jn)-trb(ji,jj,jk,jn))
     $          *umask(ji,jj,jk) 
C   ... j-gradient of passive tracer at jj
            zdjt (ji,jk)=(trb(ji,jj+1,jk,jn)-trb(ji,jj,jk,jn))
     $          *vmask(ji,jj,jk) 
C   ... j-gradient of passive tracer at jj+1
            zdj1t(ji,jk)=(trb(ji,jj,jk,jn)-trb(ji,jj-1,jk,jn))
     $          *vmask(ji,jj-1,jk)  
          END DO
        END DO
C
C
C I.2 Vertical fluxes
C -------------------
C
C Surface and bottom vertical fluxes set to zero
        DO ji = 1, jpi
          ztfw(ji, 1 ) = 0.e0
          ztfw(ji,jpk) = 0.e0
        END DO
C
C interior (2=<jk=<jpk-1)
        DO jk = 2, jpkm1
          DO ji = 2, jpim1
            zcoef0 = - fsahtw(ji,jj,jk) * tmask(ji,jj,jk)
            zmku = 1./max( umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)
     $                    +umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1. )
            zmkv = 1./max( vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)
     $                    +vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1. )
            zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
            zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)
C
            ztfw(ji,jk) = zcoef3 * ( zdit (ji  ,jk-1) + zdit (ji-1,jk)
     $                              +zdit (ji-1,jk-1) + zdit (ji  ,jk) )
     $                  + zcoef4 * ( zdjt (ji  ,jk-1) + zdj1t(ji  ,jk)
     $                              +zdj1t(ji  ,jk-1) + zdjt (ji  ,jk) )
          END DO
        END DO
C
C
C I.3  update avt
C ---------------
C
        DO jk = 2, jpkm1
          DO ji = 2, jpim1
C   ... add isopycnal vertical coeff. to avt
            avt(ji,jj,jk) = zavt(ji,jj,jk)
     $        + fsahtw(ji,jj,jk) * ( wslpi(ji,jj,jk)*wslpi(ji,jj,jk)
     $                              +wslpj(ji,jj,jk)*wslpj(ji,jj,jk) )
          END DO
        END DO
C
#    ifdef key_trahdfeiv
C
C I.4 Eddy induced vertical advective fluxes
C ------------------------------------------
C
        DO jk = 2, jpkm1
          DO ji = 2, jpim1
#if defined key_off_tra 
            zcoeg3 = -wgm(ji,jj,jk)* e1t(ji,jj)*e2t(ji,jj)*0.5
#else
            zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )
     $                              * e2u(ji-1,jj)*umask(ji-1,jj,jk)
            zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )
     $                              * e2u(ji  ,jj)*umask(ji  ,jj,jk)
            zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )
     $                              * e1v(ji,jj-1)*vmask(ji,jj-1,jk)
            zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )
     $                              * e1v(ji  ,jj)*vmask(ji  ,jj,jk)
C
            zcoeg3 = + 0.25 * tmask(ji,jj,jk) * fsaeiw(ji,jj,jk)
     $                                 * ( zuwk - zuwki + zvwk - zvwki )
#endif
C
            ztfwg(ji,jk) = + zcoeg3 *
     $          ( trb(ji,jj,jk,jn) + trb(ji,jj,jk-1,jn) )
C
            ztfw(ji,jk) = ztfw(ji,jk) + ztfwg(ji,jk)
#  if defined key_diaeiv
            wgm(ji,jj,jk) = wgm(ji,jj,jk) -
     $              2. *  zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
#  endif
          END DO
        END DO
C
#    else
C
C I.4 NO Eddy induced advective fluxes
C ----==------------------------------
C
#endif
C
C
C I.5 Divergence of vertical fluxes added to the general tracer trend
C -------------------------------------------------------------------
C
        DO jk = 1, jpkm1
          DO ji = 2, jpim1
            zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
            ztav = (  ztfw(ji,jk) - ztfw(ji,jk+1)  ) * zbtr
            tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztav
#if defined key_trc_diatrd
#        if defined key_trahdfeiv
            ztavg = ( ztfwg(ji,jk) - ztfwg(ji,jk+1) ) * zbtr
            trtrd(ji,jj,jk,jn,9) = ztavg
#        endif
            trtrd(ji,jj,jk,jn,6) = ztav - ztavg
#endif
          END DO
        END DO
