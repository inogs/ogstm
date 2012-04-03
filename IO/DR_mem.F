       MODULE DR_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
!       Domain Read matrix
! ---------------------------------------------------------------------

!      REAL(8), allocatable ::glamtglo(:,:)    , glamuglo(:,:)
!      REAL(8), allocatable ::glamvglo(:,:)    , glamfglo(:,:)
!      REAL(8), allocatable ::gphitglo(:,:)    , gphiuglo(:,:)
!      REAL(8), allocatable ::gphivglo(:,:)    , gphifglo(:,:)
!      REAL(8), allocatable ::  e1tglo(:,:)    ,   e1uglo(:,:)
!      REAL(8), allocatable ::  e1vglo(:,:)    ,   e1fglo(:,:)
!      REAL(8), allocatable ::  e2tglo(:,:)    ,   e2uglo(:,:)
!      REAL(8), allocatable ::  e2vglo(:,:)    ,   e2fglo(:,:)
!      REAL(8), allocatable ::   ffglo(:,:)
      REAL(8), allocatable ::tmaskglo(:,:,:)
!      REAL(8), allocatable ::umaskglo(:,:,:), vmaskglo(:,:,:), fmaskglo(:,:,:)
      REAL(8), allocatable ::hbattglo(:,:)    , hbatuglo(:,:)
      REAL(8), allocatable ::hbatvglo(:,:)    , hbatfglo(:,:)


!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_DR()
!      allocate(glamtglo(jpiglo,jpjglo))
!      allocate(glamuglo(jpiglo,jpjglo))
!      allocate(glamvglo(jpiglo,jpjglo))
!      allocate(glamfglo(jpiglo,jpjglo))
!      allocate(gphitglo(jpiglo,jpjglo))
!      allocate(gphiuglo(jpiglo,jpjglo))
!      allocate(gphivglo(jpiglo,jpjglo))
!      allocate(gphifglo(jpiglo,jpjglo))
!      allocate(e1tglo(jpiglo,jpjglo))
!      allocate(e1uglo(jpiglo,jpjglo))
!      allocate(e1vglo(jpiglo,jpjglo))
!      allocate(e1fglo(jpiglo,jpjglo))
!      allocate(e2tglo(jpiglo,jpjglo))
!      allocate(e2uglo(jpiglo,jpjglo))
!      allocate(e2vglo(jpiglo,jpjglo))
!      allocate(e2fglo(jpiglo,jpjglo))
!      allocate(ffglo(jpiglo,jpjglo))
      allocate(tmaskglo(jpiglo,jpjglo,jpk))
!      allocate(fmaskglo(jpiglo,jpjglo,jpk))
!      allocate(umaskglo(jpiglo,jpjglo,jpk))
!      allocate(vmaskglo(jpiglo,jpjglo,jpk))
      allocate(hbattglo(jpiglo,jpjglo))    
      allocate(hbatuglo(jpiglo,jpjglo))
      allocate(hbatvglo(jpiglo,jpjglo))
      allocate(hbatfglo(jpiglo,jpjglo))
      mem_all = rea_len*(4*jpiglo*jpjglo+ 1*jpiglo*jpjglo*jpk)
      
      END subroutine myalloc_DR

      subroutine mydealloc_DR()

      INTEGER :: I,error
      error = 0
!      deallocate(glamtglo,STAT=I)
!      error = error + I
!      deallocate(glamuglo,STAT=I)
!      error = error + I
!      deallocate(glamvglo,STAT=I)
!      error = error + I
!      deallocate(glamfglo,STAT=I)
!      error = error + I
!      deallocate(gphitglo,STAT=I)
!      error = error + I
!      deallocate(gphiuglo,STAT=I)
!      error = error + I
!      deallocate(gphivglo,STAT=I)
!      error = error + I
!      deallocate(gphifglo,STAT=I)
!      error = error + I
!      deallocate(e1tglo,STAT=I)
!      error = error + I
!      deallocate(e1uglo,STAT=I)
!      error = error + I
!      deallocate(e1vglo,STAT=I)
!      error = error + I
!      deallocate(e1fglo,STAT=I)
!      error = error + I
!      deallocate(e2tglo,STAT=I)
!      error = error + I
!      deallocate(e2uglo,STAT=I)
!      error = error + I
!      deallocate(e2vglo,STAT=I)
!      error = error + I
!      deallocate(e2fglo,STAT=I)
!      error = error + I
!      deallocate(ffglo,STAT=I)
!      error = error + I
      deallocate(tmaskglo,STAT=I)
      error = error + I
!      deallocate(fmaskglo,STAT=I)
!      error = error + I
!      deallocate(umaskglo,STAT=I)
!      error = error + I
!      deallocate(vmaskglo,STAT=I)
!      error = error + I
      deallocate(hbattglo,STAT=I)
      error = error + I
      deallocate(hbatuglo,STAT=I)
      error = error + I
      deallocate(hbatvglo,STAT=I)
      error = error + I
      deallocate(hbatfglo,STAT=I)
      error = error + I

       IF(error .NE. 0) THEN
           write(*,*) 'DR_mem.F:Error in deallocate_DR buffers'
           STOP
       ENDIF 
      END subroutine mydealloc_DR

      END MODULE 
