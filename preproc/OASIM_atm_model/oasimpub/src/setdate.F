      subroutine setdate(idmo,imo,iyr,iday,ndmo,lpyr)
c
c  Sets dates for the run.
c
      real daymo
      common /bdaymo/ daymo(12)
c
c  Leap year check
      lpyr = 0
      lpchk = mod(iyr,4)
      if (lpchk .eq. 0)lpyr = 1
c
c  Number of days in a month
      daymo(1)  = 31.0
      daymo(2)  = 28.0
      daymo(3)  = 31.0
      daymo(4)  = 30.0
      daymo(5)  = 31.0
      daymo(6)  = 30.0
      daymo(7)  = 31.0
      daymo(8)  = 31.0
      daymo(9)  = 30.0
      daymo(10) = 31.0
      daymo(11) = 30.0
      daymo(12) = 31.0
      if (lpyr .eq. 1)daymo(2) = 29.0
c
c  New run dates
      if (idmo .gt. daymo(imo))then
       idmo = 1
       imo = imo+1
       if (imo .gt. 12)then
        imo = 1
        iyr = iyr+1
       endif
      endif
      call doy(idmo,imo,iyr,iday,lpyr)
      if (lpyr .eq. 1)daymo(2) = 29.0
      istday = iday
c
c  End of month day
      ndmo = nint(daymo(imo))
c
      return
      end
c
c ******************************************************************
      subroutine doy(idmo,imo,iyr,jday,lpyr)
c
c  Computes day of year from gregorian date 
c  Includes a factor for leap year.
c
      integer iyr,imo,idmo,jday
      integer feb,apr,aug,sep,oct
      integer feblp,aprlp,auglp,seplp,octlp
      data jan,feb,mar,apr,may,jun/31,59,90,120,151,181/
      data jul,aug,sep,oct,nov/212,243,273,304,334/
      data feblp,marlp,aprlp,maylp,junlp/60,91,121,152,182/
      data jullp,auglp,seplp,octlp,novlp/213,244,274,305,335/
c
      lpchk = mod(iyr,4)
      if (lpchk .eq. 0)then
       lpyr = 1
       go to 20
      else
       lpyr = 0
      endif
      if (imo .eq. 1)jday = idmo
      if (imo .eq. 2)jday = jan+idmo
      if (imo .eq. 3)jday = feb+idmo
      if (imo .eq. 4)jday = mar+idmo
      if (imo .eq. 5)jday = apr+idmo
      if (imo .eq. 6)jday = may+idmo
      if (imo .eq. 7)jday = jun+idmo
      if (imo .eq. 8)jday = jul+idmo
      if (imo .eq. 9)jday = aug+idmo
      if (imo .eq. 10)jday = sep+idmo
      if (imo .eq. 11)jday = oct+idmo
      if (imo .eq. 12)jday = nov+idmo
      go to 30
c  Leap year
20    if (imo .eq. 1)jday = idmo
      if (imo .eq. 2)jday = jan+idmo
      if (imo .eq. 3)jday = feblp+idmo
      if (imo .eq. 4)jday = marlp+idmo
      if (imo .eq. 5)jday = aprlp+idmo
      if (imo .eq. 6)jday = maylp+idmo
      if (imo .eq. 7)jday = junlp+idmo
      if (imo .eq. 8)jday = jullp+idmo
      if (imo .eq. 9)jday = auglp+idmo
      if (imo .eq. 10)jday = seplp+idmo
      if (imo .eq. 11)jday = octlp+idmo
      if (imo .eq. 12)jday = novlp+idmo
30    continue
c
      return
      end
