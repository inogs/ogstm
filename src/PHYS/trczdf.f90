      SUBROUTINE trczdf
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trczdf
!!!                     ******************
!!!
!!!  Purpose :
!!!  --------
!!!     Compute the trend due to the vertical tracer diffusion inclu-
!!!     ding the vertical component of lateral mixing (only for second
!!!     order operator, for fourth order it is already computed and
!!!     add to the general trend in trchdf.F) and the surface forcing
!!!     and add it to the general trend of the tracer equations.
!!
!!   Method :
!!   -------
!!         The vertical component of the lateral diffusive trends is
!!      provided by a 2nd order operator rotated along neural or geopo-
!!      tential surfaces to which an eddy induced advection can be added
!!      It is computed using before fields (forward in time) and isopyc-
!!      nal or geopotential slopes computed in routine hdfslp.
!!
!!    First part: vertical trends associated with the lateral mixing
!!    ==========  (excluding the vertical flux proportional to dk[tr] )
!!                ('key_trahdfiso' or 'key_trahdfgeop')
!!    vertical fluxes associated with the rotated lateral mixing:
!!       zftw =-ahtt {  e2t*wslpi di[ mi(mk(trb)) ]
!!                   + e1t*wslpj dj[ mj(mk(trb)) ]  }
!!    save avt coef. resulting from vertical physics alone in zavt:
!!       zavt = avt
!!    update and save in zavt the vertical eddy viscosity coefficient:
!!       avt = avt + wslpi^2+wslj^2
!!    add vertical Eddy Induced advective fluxes ('key_trahdfeiv'):
!!       zftw = zftw + { di[ahtt e2u mi(wslpi)]
!!              +dj[ahtt e1v mj(wslpj)] } mk(trb)
!!    take the horizontal divergence of the fluxes:
!!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
!!      Add this trend to the general trend (tra):
!!         tra = tra + difft
!!
!!    Second part: vertical trend associated with the vertical physics
!!    ===========  (including the vertical flux proportional to dk[tr]
!!              associated with the lateral mixing, through the
!!              update of avt)
!!    The vertical diffusion of passive tracers (tr) is given by:
!!             difft = dz( avt dz(tr) ) = 1/e3t dk+1( avt/e3w dk(tr) )
!!         'key_zdfexplicit' defined: forward  time scheme (tr=trb), using
!!       a time splitting technique (NOT YET IMPLEMENTED).
!!          default option          : backward time scheme, tr=tra
!!                                  : use ndttrc or 2 for time stepping
!!      Bottom boundary conditions:
!!         no flux on passive tracers: applied through the masked field avt
!!      Add this trend to the general trend tra :
!!         tra = tra + dz( avt dz(tr) )
!!



       USE myalloc
       USE myalloc_mpp
       USE ZDF_mem
       USE DIA_mem
        IMPLICIT NONE

!!---------------------------------------------------------------------
!! local declarations
!! ==================


      LOGICAL l1,l2,l3
      INTEGER ji, jj, jk, jn, jv, jf
! omp variables
      INTEGER :: mytid, ntids

#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      REAL(8) ztavg, zdt
      INTEGER ikst, ikenm2, ikstp1
!! only IMPLICIT scheme

      REAL(8) z2dtt
      REAL(8) delta_tra(jpk),int_tra(jpk)
      REAL(8) Aij

!!---------------------------------------------------------------------
!! statement functions
!! ===================


      trczdfparttime = MPI_WTIME() ! Cronometer start
#ifdef __OPENMP1
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
         ! Cross matrix between transect and local optimized indexing
         jarr_zdf_flx=0

         DO jf=1,Fsize
            DO jv=1, dimen_jvzdf
               DO jk=1,jpk
                  l1 = flx_ridxt(jf,2) .EQ. jarr_zdf(1,jv)
                  l2 = flx_ridxt(jf,3) .EQ. jarr_zdf(2,jv)
                  l3 = flx_ridxt(jf,4) .EQ. jk
                  IF ( l1 .AND. l2 .AND. l3) THEN
                     jarr_zdf_flx(jv,jk)= jf
                  END IF
               END DO
            END DO
         END DO

      ENDIF


!! passive tracer slab
!! ===================

!!!$omp    parallel do default(none) private(jv,ji,jj,jk,mytid,ztavg,zdt,z2dtt,ikst,ikstp1,ikenm2,jf,delta_tra,int_tra,Aij,jn)
!!!$omp&                           shared(dimen_jvzdf,jarr_zdf,ndttrc,rdt,jpk,jpkm1,avt,zwi,
!!!$omp&                                   zws,zwd,zwy,trb,tra,zwt,zwz,zwx,tmask,e1t,e2t,e3t,e3w,jarr_zdf_flx,Fsize,
!!!$omp&                                   diaflx)
      TRACER_LOOP: DO  jn = 1, jptra

!!!$         mytid = omp_get_thread_num()  ! take the thread ID

!      if( mytid + jn <= jptra ) then
        ztavg = 0.e0
!! vertical slab

        DO jv = 1, dimen_jvzdf

           ji  = jarr_zdf(1,jv)
           jj  = jarr_zdf(2,jv)
           Aij = e1t(ji,jj) * e2t(ji,jj)

!! I. Vertical trends associated with lateral mixing
!! -------------------------------------------------
!!    (excluding the vertical flux proportional to dk[t] )
!! II. Vertical trend associated with the vertical physics
!! -------------------------------------------------------
!!     (including the vertical flux proportional to dk[t] associated
!!      with the lateral mixing, through the avt update)
!! II. dk[ avt dk[ (tr) ] ] diffusive trends
!! ==========================================

!! 0. initialization
      zdt= ndttrc


!! II.0 Matrix construction
!! Diagonal, inferior, superior
!! (including the bottom boundary condition via avt masked)
!!   ... Euler time stepping when starting from rest
        DO jk = 1, jpkm1
          z2dtt = zdt * rdt
          zwi(jk,mytid+1) = - z2dtt * avt(ji,jj,jk  )/( e3t(ji,jj,jk) * e3w(ji,jj,jk  ) )
          zws(jk,mytid+1) = - z2dtt * avt(ji,jj,jk+1)/( e3t(ji,jj,jk) * e3w(ji,jj,jk+1) )
          zwd(jk,mytid+1) = 1. - zwi(jk,mytid+1) - zws(jk,mytid+1)
        END DO

!! Surface boundary conditions
          zwi(1,mytid+1) = 0.e0
          zwd(1,mytid+1) = 1. - zws(1,mytid+1)

!! II.1. Vertical diffusion on tr
!! ------------------------------
!! Second member construction
!!   ... Euler time stepping when starting from rest
        DO jk = 1, jpkm1
          z2dtt = zdt * rdt
          zwy(jk,mytid+1) = trb(ji,jj,jk,jn) + z2dtt * tra(ji,jj,jk,jn)
        END DO

!! Matrix inversion from the first level
        ikst = 1
!! It is possible to change the maximum bottom value
!! for matrix inversion (see Marina Levy s thesis)
!! That will be done with some more documentation



!!!        ZDF.MATRIXSOLVER
!!!      ********************
!!!
!! Matrix inversion
!!----------------------------------------------------------------------
!!   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
!!
!!        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
!!        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
!!        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
!!        (        ...               )( ...  ) ( ...  )
!!        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
!!
!!   m is decomposed in the product of an upper and lower triangular
!!   matrix
!!   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
!!   The second member is in 2d array zwy
!!   The solution is in 2d array zwx
!!   The 2d arry zwt and zwz are work space arrays
!!
!!   N.B. the starting vertical index (ikst) is equal to 1 except for
!!   the resolution of tke matrix where surface tke value is prescribed
!!   so that ikstrt=2.

        ikstp1=ikst+1
        ikenm2=jpk-2

        zwt(ikst,mytid+1)=zwd(ikst,mytid+1)

        DO jk=ikstp1,jpkm1
            zwt(jk,mytid+1)=zwd(jk,mytid+1)-zwi(jk,mytid+1)*zws(jk-1,mytid+1)/zwt(jk-1,mytid+1)
        END DO

          zwz(ikst,mytid+1)=zwy(ikst,mytid+1)

        DO jk=ikstp1,jpkm1
            zwz(jk,mytid+1)=zwy(jk,mytid+1)-zwi(jk,mytid+1)/zwt(jk-1,mytid+1)*zwz(jk-1,mytid+1)
        END DO

        zwx(jpkm1,mytid+1)=zwz(jpkm1,mytid+1)/zwt(jpkm1,mytid+1)

        DO jk=ikenm2,ikst,-1
            zwx(jk,mytid+1)=( zwz(jk,mytid+1)-zws(jk,mytid+1)*zwx(jk+1,mytid+1) )/zwt(jk,mytid+1)
        END DO

! calculate flux due to vertical diffusion (on top face of the grid cell jk)
! the flux is computed as an integral (int_tra) starting from the surface

      DO jk=1,jpkm1

         z2dtt = zdt * rdt
         delta_tra(jk) = ( zwx(jk,mytid+1) - zwy(jk,mytid+1) ) / z2dtt * Aij * e3t(ji,jj,jk)! or trn(ji,jj,jk,jn+mytid)

         IF (jk .EQ. 1) THEN
             int_tra(1)  = 0
         ELSE
             int_tra(jk) = delta_tra(jk) + int_tra(jk-1)
         ENDIF

      jf = jarr_zdf_flx(jv,jk)

! jf = 0 are the points not included in the transect they are excluded

         IF ((Fsize .GT. 0) .AND. (jf .GT. 0)) THEN

                  diaflx(jf,jn,7) = int_tra(jk)

         ENDIF

      ENDDO



!! Save the masked passive tracer after in tra
!! (c a u t i o n:  tracer not its trend, Leap-frog scheme done
!!                  it will not be done in trcnxt)
         DO jk = 1, jpkm1
            tra(ji,jj,jk,jn) = zwx(jk,mytid+1) * tmask(ji,jj,jk)
         END DO

        END DO ! jv

!      end if

       END DO TRACER_LOOP
!!!$omp    end parallel do



       trczdfparttime = MPI_WTIME() - trczdfparttime   !cronometer-stop
       trczdftottime = trczdftottime + trczdfparttime

      END SUBROUTINE trczdf
