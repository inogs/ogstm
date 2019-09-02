       subroutine ocalbedo(rad,sunz,ws,rod,ros)
c
c  Computes ocean surface albedo from solar zenith angle (sunz) 
c  and wind speed (ws, m/s).
c  Albedo is provided as direct (albd) and diffuse (albs).     
c
      parameter(nlt=33)
      real rod(nlt),ros(nlt)
c
c  Derive surface reflectance as a function of sunz and ws
c  for OASIM Bands
      call sfcrfl(rad,sunz,ws,rod,ros)
c
      return
      end
