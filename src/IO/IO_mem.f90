       MODULE IO_mem 

       USE modul_param 
       USE myalloc

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
! Common/comcoh/  : IO matrix
! ---------------------------------------------------------------------


      INTEGER :: jpi_rec, jpj_rec
      INTEGER :: jpi_max, jpj_max
      double precision :: elapsed_time_1=0.0, elapsed_time_2 = 0.0
      LOGICAL :: existFilebkp = .false.
      double precision, allocatable :: buffglamt(:) 
      double precision, allocatable :: buffgphit(:)
      double precision, allocatable :: bufftrn(:)
      double precision, allocatable :: bufftrb(:)
      double precision, allocatable :: buffsn(:)
      double precision, allocatable :: bufftn(:)
      double precision, allocatable :: buffvatm(:)
      double precision, allocatable :: buffemp(:)
      double precision, allocatable :: buffqsr(:)
      double precision, allocatable :: buffun(:)
      double precision, allocatable :: buffvn(:)
      double precision, allocatable :: buffwn(:)
      double precision, allocatable :: buffavt(:)
      double precision, allocatable :: buffe3t(:)
      double precision, allocatable :: bufftma(:)
      double precision, allocatable :: bufftrIO(:)
      double precision, allocatable :: buffDIA(:)
      double precision, allocatable :: buffDIA2d(:)
      real, allocatable :: d2f2d(:,:)
!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_IO()
      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

       allocate(buffglamt (jpi_max* jpj_max))    
        buffglamt = huge(buffglamt(1))
       allocate(buffgphit (jpi_max* jpj_max))      
        buffgphit = huge(buffgphit(1))
       allocate(bufftrn   (jpi_max* jpj_max* jpk)) 
        bufftrn   = huge(bufftrn(1))
       allocate(bufftrb   (jpi_max* jpj_max* jpk)) 
        bufftrb   = huge(bufftrb(1))
       allocate(buffsn    (jpi_max *jpj_max* jpk)) 
        buffsn    = huge(buffsn(1))
       allocate(bufftn    (jpi_max* jpj_max* jpk)) 
        bufftn    = huge(bufftn(1))
       allocate(buffvatm  (jpi_max* jpj_max))      
        buffvatm  = huge(buffvatm(1))
       allocate(buffemp   (jpi_max* jpj_max))      
        buffemp   = huge(buffemp(1))
       allocate(buffqsr   (jpi_max* jpj_max))      
        buffqsr   = huge(buffqsr(1))
       allocate(buffun    (jpi_max* jpj_max* jpk)) 
        buffun    = huge(buffun(1))
       allocate(buffvn    (jpi_max* jpj_max* jpk)) 
        buffvn    = huge(buffvn(1))
       allocate(buffwn    (jpi_max* jpj_max* jpk)) 
        buffwn    = huge(buffwn(1))
       allocate(buffavt   (jpi_max* jpj_max* jpk)) 
        buffavt   = huge(buffavt(1))
       allocate(buffe3t   (jpi_max* jpj_max* jpk)) 
        buffe3t   = huge(buffe3t(1))
       allocate(bufftma   (jpi_max* jpj_max* jpk)) 
        bufftma   = huge(bufftma(1))
       allocate(bufftrIO  (jpi_max* jpj_max* jpk)) 
        bufftrIO  = huge(bufftrIO(1))
       allocate(buffDIA   (jpi_max* jpj_max* jpk)) 
        buffDIA   = huge(buffDIA(1))
       allocate(buffDIA2d (jpi_max* jpj_max     )) 
        buffDIA2d = huge(buffDIA2d(1))

       if (lwp) then
       allocate(d2f2d     (jpjglo,jpiglo))         
        d2f2d     = huge(d2f2d(1,1))
       endif

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_IO



      subroutine clean_memory_io()

          deallocate(buffglamt)
          deallocate(buffgphit)
          deallocate(bufftrn)
          deallocate(bufftrb)
          deallocate(buffsn)
          deallocate(bufftn)
          deallocate(buffvatm)
          deallocate(buffemp)
          deallocate(buffqsr)
          deallocate(buffun)
          deallocate(buffvn)
          deallocate(buffwn)
          deallocate(buffavt)
          deallocate(buffe3t)
          deallocate(bufftma)
          deallocate(bufftrIO)
          deallocate(buffDIA)
          deallocate(buffDIA2d)
          
          if (lwp) then
              deallocate(d2f2d)
          endif

      end subroutine clean_memory_io



      END MODULE 
