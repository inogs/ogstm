       MODULE DR_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
!       Domain Read matrix
! ---------------------------------------------------------------------

#ifdef key_s_coord
      REAL(8), allocatable ::hbattglo(:,:)    , hbatuglo(:,:)
      REAL(8), allocatable ::hbatvglo(:,:)    , hbatfglo(:,:)
#endif

!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_DR()

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
