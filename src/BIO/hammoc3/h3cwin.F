C $Id: h3cwin.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
      SUBROUTINE h3cwin
#if defined key_passivetrc
c -------------------------------------------------------------------
c
c sous programme de calcul du terme xKw du coefficient d echange de gaz
c si c est le premier passage : lecture du fichier 
c
c si pas de temps kt
c 
c  interpolation des champs mensuels
C
C  xkw=Xconv*a*(u^2+v) en m/s
c  references: guidelines ocmip
c
c -------------------------------------------------------------------
c

       USE myalloc

        IMPLICIT NONE

#    if defined key_trc_hamocc3
      INTEGER i, j, l, mo,nvxk
      
      nvxk   = 96
      WRITE(numout,*) 'coefficients echange de gaz: initialisation'
      WRITE(numout,*)
      WRITE(numout,1111) nvxk
 1111 FORMAT ('unite logique pour xkw :nvxk :',i7)
      WRITE(numout,'('' ---------------------------------- '')')
c
C Lecture  du fichier pression atmospherique
C--------------------------------
      DO mo=1,12
#        if defined key_vpp
        CALL read2s(nvxk,wsmo(1,1,mo),jpi,jpj)
#        else
        READ (nvxk) wsmo(:,:,mo)
#        endif
      END DO
#    else
C
C no wind reading
C
#    endif
#endif
      RETURN
      END

