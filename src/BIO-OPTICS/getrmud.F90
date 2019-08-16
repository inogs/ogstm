      subroutine getrmud(sunz,rmud)
      IMPLICIT NONE 
!  Computes average cosine for direct irradiance just below
!  sea surface
      real(8), intent(in)  :: sunz
      real(8), intent(out) :: rmud
      real(8) rn,pi 
      real(8) rad, rsza, sinszaw
      real(8) szaw, rmudl
      data rn /1.341/  !refractive index of seawater
      data pi /3.14159265359/
 
!  Compute average cosine for direct irradiance in the water 
!  column given solar zenith angle (in degrees) at surface.
      rad = 180.0D0/pi 
      rsza = sunz/rad
      sinszaw = sin(rsza)/rn
      szaw = asin(sinszaw)
      rmudl = 1.0/cos(szaw)   !avg cosine direct (1 over)
      rmud = min(rmudl,1.5)
      rmud = max(rmud,0.0)
 
      return
      end
