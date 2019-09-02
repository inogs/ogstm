      subroutine rstfle (idmo, imo, iyr, ndmo)
c
c  Makes a file containing month, day, and year representing
c  the last simulated date to be used for restarting model
c  runs.
c
      character*1 ctmp1
      character*2 cday, cmo, cyr
      character*6 cdate
      character*4 cdate2
      common /bdaymo/ daymo(12)
c
c  Make year last two digits
      if (iyr .lt. 2000)then
       iy = iyr-1900
      else
       iy = iyr-2000
      endif
c
c  Assign file name extensions and write to data file
      if (idmo .lt. 10)then
       write(ctmp1,'(i1)')idmo
       cday = '0'//ctmp1
      else
       write(cday,'(i2)')idmo
      endif
      if (imo .lt. 10)then
       write(ctmp1,'(i1)')imo
       cmo = '0'//ctmp1
      else
       write(cmo,'(i2)')imo
      endif
      write(cyr,'(i2)')iy
      cdate = cmo//cday//cyr
      open(4,file='File',status='unknown',form='formatted')
      write(4,'(a6)') cdate
      close(4)
      cdate2 = cmo//cyr
      open(4,file='OUT',status='unknown',form='formatted')
      write(4,'(a4)') cdate2
      close(4)
c
c  Check on end of month
      if (idmo .ge. ndmo)then
       open(4, file='MONTH', status='unknown', form='formatted')
       write(4,'(a2)')cmo
       close(4)
       open(4, file='YEAR', status='unknown', form='formatted')
       write(4,'(a2)')cyr
       close(4)
      endif
c
      return
      end
