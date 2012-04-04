       MODULE DYN_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
!     Dynamic data read from file matrix
! ---------------------------------------------------------------------

#ifdef key_trahdfeiv
      REAL(8), allocatable, dimension(:,:,:) :: zaeiu,zaeiv, zaeiw
#endif
#ifdef key_oper
      INTEGER dynrea_flag, kenr_oper 
#endif

!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_DYN()

#ifdef key_trahdfeiv
       allocate(zaeiu   (jpiglo,jpjglo,jpk))
       allocate(zaeiv   (jpiglo,jpjglo,jpk))
       allocate(zaeiw   (jpiglo,jpjglo,jpk))
#endif
#ifdef key_oper
      dynrea_flag = 1
#endif
      mem_all = 0

#ifdef key_trahdfeiv
     &          +8*(3*jpiglo*jpjglo*jpk)
#endif


      END subroutine myalloc_DYN

      END MODULE



!      INTEGER, allocatable :: imel(:,:)
!      REAL(4), allocatable, dimension(:,:,:) :: zu_f, zv_f, zw_f,zt_f,zs_f
!      REAL(4), allocatable, dimension(:,:,:) :: zflx, zavt_f,

!      REAL(4), allocatable, dimension(:,:)   :: zbblx_f,zbbly_f,zemp_f,zqsr_f
!      REAL(4), allocatable, dimension(:,:)   ::  zlon,zlat zwind_f,zice_f,
!      REAL(8), allocatable, dimension(:)     :: zlev
!       allocate(imel    (jpiglo,jpjglo    ))
!       allocate(zlon    (jpiglo,jpjglo    ))
!       allocate(zlat    (jpiglo,jpjglo    ))
!       allocate(zlev    (              jpk))
!       allocate(zu_f    (jpiglo,jpjglo,jpk))
!       allocate(zv_f    (jpiglo,jpjglo,jpk))
!       allocate(zw_f    (jpiglo,jpjglo,jpk))
!       allocate(zt_f    (jpiglo,jpjglo,jpk))
!       allocate(zs_f    (jpiglo,jpjglo,jpk))
!       allocate(zbblx_f (jpiglo,jpjglo    ))
!       allocate(zbbly_f (jpiglo,jpjglo    ))
!       allocate(zemp_f  (jpiglo,jpjglo    ))
!       allocate(zqsr_f  (jpiglo,jpjglo    ))
!       allocate(zavt_f  (jpiglo,jpjglo,jpk))
!       allocate(zwind_ f(jpiglo,jpjglo    ))
!       allocate(zice_f  (jpiglo,jpjglo    ))
!       allocate(zflx    (jpiglo,jpjglo,jpflx))

!       mem_all = 4*(0*jpiglo*jpjglo*jpk
!     &           + 7* jpiglo*jpjglo )+
!     &          8*(jpiglo*jpjglo*jpflx+
!     &            2*jpiglo*jpjglo + jpk )
