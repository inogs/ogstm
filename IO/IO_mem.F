       MODULE IO_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
! Common/comcoh/  : IO matrix
! ---------------------------------------------------------------------


      INTEGER :: jpi_rec, jpj_rec
      INTEGER :: jpi_max, jpj_max
      INTEGER :: ave_counter=0
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
      REAL(8), allocatable :: buffbblx(:)
      REAL(8), allocatable :: buffvn(:)
      REAL(8), allocatable :: buffbbly(:)
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

       allocate(buffglamt (jpi_max* jpj_max))
       allocate(buffgphit (jpi_max* jpj_max))
       allocate(bufftrn   (jpi_max* jpj_max* jpk))
       allocate(bufftrb   (jpi_max* jpj_max* jpk))
       allocate(buffsn    (jpi_max *jpj_max* jpk))
       allocate(bufftn    (jpi_max* jpj_max* jpk))
       allocate(buffvatm  (jpi_max* jpj_max))
       allocate(buffemp   (jpi_max* jpj_max))
       allocate(buffqsr   (jpi_max* jpj_max))
       allocate(buffun    (jpi_max* jpj_max* jpk))
       allocate(buffbblx  (jpi_max* jpj_max))
       allocate(buffvn    (jpi_max* jpj_max* jpk))
       allocate(buffbbly  (jpi_max* jpj_max))
       allocate(buffwn    (jpi_max* jpj_max* jpk))
       allocate(buffavt   (jpi_max* jpj_max* jpk))
       allocate(bufftma   (jpi_max* jpj_max* jpk))
       allocate(bufftrIO  (jpi_max* jpj_max* jpk))
       allocate(buffpp    (jpi_max* jpj_max* jpk))

       if (lwp) then
       allocate(d2f3d     (jpiglo,jpjglo,jpk))
       allocate(d2f2d     (jpiglo,jpjglo))
       endif

      mem_all = rea_len *(7*jpi_max*jpj_max + (10+jptra_dia)*jpi_max* jpj_max*jpk) +
     &    4* (jpiglo*jpjglo*jpk + jpiglo*jpjglo )


      END subroutine myalloc_IO



      END MODULE 
