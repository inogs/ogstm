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
      REAL(8), allocatable :: bufftrn_FN(:)
      REAL(8), allocatable :: buffpp(:)
      REAL(4), allocatable :: d2f3d(:,:,:)
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
       allocate(bufftrn_FN(jpi_max* jpj_max* jpk))
       allocate(buffpp    (jpi_max* jpj_max* jpk*jptra_dia))
       allocate(d2f3d     (jpiglo,jpjglo,jpk))

      mem_all = rea_len *(7*jpi_max* jpj_max + 11*jpi_max* jpj_max* jpk +
     &           jpi_max* jpj_max* jpk*3) + 4 * (jpiglo + jpjglo + jpk  + jpiglo* jpjglo* jpk )
      END subroutine myalloc_IO

      subroutine mydealloc_IO()

      INTEGER :: I,error
       error = 0
       deallocate(buffglamt,STAT=I)
       error = error + I 
       deallocate(buffgphit,STAT=I)
       error = error + I
       deallocate(bufftrn,STAT=I)
       error = error + I
       deallocate(bufftrb,STAT=I)
       error = error + I
       deallocate(buffsn,STAT=I)
       error = error + I
       deallocate(bufftn,STAT=I)
       error = error + I
       deallocate(buffvatm,STAT=I)
       error = error + I
       deallocate(buffemp,STAT=I)
       error = error + I
       deallocate(buffqsr,STAT=I)
       error = error + I
       deallocate(buffun,STAT=I)
       error = error + I
       deallocate(buffbblx,STAT=I)
       error = error + I
       deallocate(buffvn,STAT=I)
       error = error + I
       deallocate(buffbbly,STAT=I)
       error = error + I
       deallocate(buffwn,STAT=I)
       error = error + I
       deallocate(buffavt,STAT=I)
       error = error + I
       deallocate(bufftma,STAT=I)
       error = error + I
       deallocate(bufftrIO,STAT=I)
       error = error + I
       deallocate(bufftrn_FN,STAT=I)
       error = error + I
       deallocate(buffpp,STAT=I)
       error = error + I

       IF(error .NE. 0) THEN
           write(*,*) 'IO_mem.F:Error in deallocate_IO buffers'
           STOP
       ENDIF 
      END subroutine mydealloc_IO

      END MODULE 
