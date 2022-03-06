      subroutine getrmud(sunz,rmud)
      USE myalloc, ONLY:jpj,jpi
      USE OPT_mem, ONLY: refrac_idx, rad
      IMPLICIT NONE 
!  Computes average cosine for direct irradiance just below
!  sea surface
      double precision,  intent(in) :: sunz(jpj,jpi)
      double precision, intent(out) :: rmud(jpj,jpi)
      double precision              ::  rsza(jpj,jpi), sinszaw(jpj,jpi)
      double precision              ::  szaw(jpj,jpi), rmudl(jpj,jpi)
 
!  Compute average cosine for direct irradiance in the water 
!  column given solar zenith angle (in degrees) at surface.
      rsza = sunz/rad
      sinszaw = sin(rsza)/refrac_idx
      szaw = asin(sinszaw)
      rmudl = 1.0D0/cos(szaw)   !avg cosine direct (1 over)
      rmud = min(rmudl,1.5D0)
      rmud = max(rmud,0.0D0)
 
      return
      end
