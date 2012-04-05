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
      !REAL(8), allocatable ::tmaskglo(:,:,:)
!      REAL(8), allocatable ::umaskglo(:,:,:), vmaskglo(:,:,:), fmaskglo(:,:,:)
#ifdef key_s_coord
      REAL(8), allocatable ::hbattglo(:,:)    , hbatuglo(:,:)
      REAL(8), allocatable ::hbatvglo(:,:)    , hbatfglo(:,:)
#endif

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
!      allocate(tmaskglo(jpiglo,jpjglo,jpk))
!      allocate(fmaskglo(jpiglo,jpjglo,jpk))
!      allocate(umaskglo(jpiglo,jpjglo,jpk))
!      allocate(vmaskglo(jpiglo,jpjglo,jpk))
#ifdef key_s_coord
      allocate(hbattglo(jpiglo,jpjglo))    
      allocate(hbatuglo(jpiglo,jpjglo))
      allocate(hbatvglo(jpiglo,jpjglo))
      allocate(hbatfglo(jpiglo,jpjglo))
#endif
      mem_all = rea_len*(1*jpiglo*jpjglo*jpk)
#ifdef key_s_coord
     & + rea_len*4*jpiglo*jpjglo
#endif
      
      END subroutine myalloc_DR



      END MODULE 
