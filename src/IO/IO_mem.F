       MODULE IO_mem 

       USE modul_param 
       USE myalloc

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif Mem_Monitor

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
! Common/comcoh/  : IO matrix
! ---------------------------------------------------------------------


      INTEGER :: jpi_rec, jpj_rec
      INTEGER :: jpi_max, jpj_max
      INTEGER :: ave_counter_1=0, ave_counter_2=0
      LOGICAL :: existFilebkp = .false.
      REAL(8), allocatable :: buffglamt(:) 
      REAL(8), allocatable :: buffgphit(:)
      REAL(8), allocatable :: bufftrn(:)
      REAL(8), allocatable :: bufftrb(:)
      REAL(8), allocatable :: buffsn(:)
      REAL(8), allocatable :: bufftn(:)
      REAL(8), allocatable :: buffvatm(:)
      REAL(8), allocatable :: buffemp(:)
      REAL(8), allocatable :: buffqsr(:)
      REAL(8), allocatable :: buffun(:)
      REAL(8), allocatable :: buffvn(:)
      REAL(8), allocatable :: buffwn(:)
      REAL(8), allocatable :: buffavt(:)
      REAL(8), allocatable :: bufftma(:)
      REAL(8), allocatable :: bufftrIO(:)
      REAL(8), allocatable :: buffpp(:)
      REAL(4), allocatable :: d2f3d(:,:,:)
      REAL(4), allocatable :: d2f2d(:,:)
!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_IO()
      INTEGER  :: err
      REAL(8)  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif Mem_Monitor

       allocate(buffglamt (jpi_max* jpj_max))      ; buffglamt = huge(buffglamt(1))
       allocate(buffgphit (jpi_max* jpj_max))      ; buffgphit = huge(buffgphit(1))
       allocate(bufftrn   (jpi_max* jpj_max* jpk)) ; bufftrn   = huge(bufftrn(1))
       allocate(bufftrb   (jpi_max* jpj_max* jpk)) ; bufftrb   = huge(bufftrb(1))
       allocate(buffsn    (jpi_max *jpj_max* jpk)) ; buffsn    = huge(buffsn(1))
       allocate(bufftn    (jpi_max* jpj_max* jpk)) ; bufftn    = huge(bufftn(1))
       allocate(buffvatm  (jpi_max* jpj_max))      ; buffvatm  = huge(buffvatm(1))
       allocate(buffemp   (jpi_max* jpj_max))      ; buffemp   = huge(buffemp(1))
       allocate(buffqsr   (jpi_max* jpj_max))      ; buffqsr   = huge(buffqsr(1))
       allocate(buffun    (jpi_max* jpj_max* jpk)) ; buffun    = huge(buffun(1))
       allocate(buffvn    (jpi_max* jpj_max* jpk)) ; buffvn    = huge(buffvn(1))
       allocate(buffwn    (jpi_max* jpj_max* jpk)) ; buffwn    = huge(buffwn(1))
       allocate(buffavt   (jpi_max* jpj_max* jpk)) ; buffavt   = huge(buffavt(1))
       allocate(bufftma   (jpi_max* jpj_max* jpk)) ; bufftma   = huge(bufftma(1))
       allocate(bufftrIO  (jpi_max* jpj_max* jpk)) ; bufftrIO  = huge(bufftrIO(1))
       allocate(buffpp    (jpi_max* jpj_max* jpk)) ; buffpp    = huge(buffpp(1))

       if (lwp) then
       allocate(d2f3d     (jpiglo,jpjglo,jpk))     ; d2f3d     = huge(d2f3d(1,1,1))
       allocate(d2f2d     (jpiglo,jpjglo))         ; d2f2d     = huge(d2f2d(1,1))
       endif

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif Mem_Monitor

      END subroutine myalloc_IO



      END MODULE 
