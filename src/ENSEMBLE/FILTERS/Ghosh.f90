 
! developed by Simone Spada (sspada@ogs.it) at OGS

module Ghosh
    use mpi
    
    use myalloc, &
        only: lwp, myrank, mysize, mycomm
    USE ogstm_mpi_module, &
        ONLY: glcomm
    use Ens_Mem, &
        only: EnsRankZero, myrankZero, EnsIOUnit, &
            EnsDebug, EnsRank, EnsSize, &
            UseLocalObsDumping, &
            ForgettingFactor, &
            Ens_shared_alloc, Ens_shared_alloc_base, &
            EnsWeights
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAstate_kjit, win_DAstate, n_DAstate, gl_DAstate, gl_DAstate_kjitn, &
            DAMask
    use Ens_Utilities, &
        only: Ens_ReduceMeanAndBase
    use Ens_Obs, &
        only: nObs, n_SqrtR1, &
            Misfit, H, SqrtR1
    use LocalOperations, &
        only: LocalSpace, &
            Ens_Init_Local, Ens_Finalize_Local
    use Algebra, &
        only: Init_Algebra, Finalize_Algebra, SymChangeBase, OrtMatrix, det
    
    
    implicit none
    
    logical, parameter :: LocalAnalysis=.true.
    
    double precision, pointer, contiguous, dimension(:,:,:) :: Hstate, SqrtR1HLi, HLTR1HLi_loc
    integer :: n_Hstate, n_SqrtR1HLi, n_HLTR1HLi_loc
    integer :: win_Hstate, win_SqrtR1HLi, win_HLTR1HLi_loc
    double precision, pointer, contiguous, dimension(:,:) :: gl_Hstate
    double precision, pointer, contiguous, dimension(:,:,:,:) :: gl_SqrtR1HLi, gl_HLTR1HLi_loc
    
    double precision, dimension(:,:), allocatable :: EnsCov1, ForecastCov1, TTW1T
    
    double precision, dimension(:,:), allocatable :: OrtMatrixSampling
    
    double precision, pointer, contiguous, dimension(:) :: ChangeBase
    integer :: n_ChangeBase
    integer :: win_ChangeBase
    double precision, pointer, contiguous, dimension(:,:) :: gl_ChangeBase
    
    double precision, pointer, contiguous, dimension(:,:,:,:) :: NewBase
    integer :: n_NewBase
    integer :: win_NewBase
    double precision, pointer, contiguous, dimension(:,:,:,:,:) :: gl_NewBase
    
    type(LocalSpace) :: LocalPatch
    
    integer :: GhoshOrder
    
    
    
contains

subroutine GHOSH_Namelist(filename)
        
    character(len=*) :: filename
    
    NAMELIST/GHOSH_setup/ GhoshOrder
    
    if (lwp) then
        
    end if
    
    OPEN(unit=EnsIOUnit, file=filename, status='OLD')
    
        GhoshOrder=5

        REWIND(EnsIOUnit)
        READ(EnsIOUnit, GHOSH_setup)

        IF(lwp) THEN
            write(*,*) 'Namelist GHOSH parameters:'
            WRITE(*,*) ' '
            WRITE(*,*) ' GhoshOrder: ', GhoshOrder
        END IF

    CLOSE(EnsIOUnit)
    
end subroutine

subroutine Ghosh_Init

    double precision, dimension(:), POINTER, contiguous :: member_pointer
    double precision, dimension(:,:), POINTER, contiguous :: global_pointer
    
    integer ierr
    
    call GHOSH_Namelist('namelist.init')
    
    if (EnsDebug>1) then
        call mpi_barrier(glcomm,ierr)
        if (lwp) write(*,*) "Ghosh_Init start" 
    end if
    
    n_HLTR1HLi_loc=nj_DAstate*ni_DAstate*(EnsSize-1)
    call Ens_shared_alloc(n_HLTR1HLi_loc, member_pointer, global_pointer, win_HLTR1HLi_loc)
    HLTR1HLi_loc(1:EnsSize-1, 1:nj_DAstate, 1:ni_DAstate)=>member_pointer
    gl_HLTR1HLi_loc(1:EnsSize-1, 1:nj_DAstate, 1:ni_DAstate, 0:EnsSize-1)=>global_pointer
    HLTR1HLi_loc  = huge(HLTR1HLi_loc(1,1,1))
    
    allocate(EnsCov1(EnsSize-1, EnsSize-1))
    EnsCov1=Huge(EnsCov1(1,1))
    
    allocate(ForecastCov1(EnsSize-1, EnsSize-1))
    ForecastCov1=Huge(ForecastCov1(1,1))
    
    allocate(TTW1T(EnsSize-1,EnsSize-1))
    TTW1T=Huge(TTW1T(1,1))
    call TTW1T_builder
    
    allocate(OrtMatrixSampling(EnsSize, EnsSize))
    OrtMatrixSampling=Huge(OrtMatrixSampling(1,1))
    
    if (LocalAnalysis) then
        call Ens_Init_Local
        call LocalPatch%Init(EnsSize-1)
    end if
    
    if (EnsDebug>1) then
        call mpi_barrier(glcomm,ierr)
        if (lwp) write(*,*) "Ghosh_Init: before Init_Algebra" 
    end if
    
    call Init_Algebra
    
    if (EnsDebug>1) then
        call mpi_barrier(glcomm,ierr)
        if (lwp) write(*,*) "Ghosh_Init: Init_Algebra done" 
    end if
    
    n_ChangeBase=EnsSize-1
    call Ens_shared_alloc(n_ChangeBase, ChangeBase, global_pointer, win_ChangeBase)
    gl_ChangeBase(1:EnsSize-1, 0:EnsSize-1)=>global_pointer
    ChangeBase  = huge(ChangeBase(1))
    
    if (EnsDebug>1) then
        call mpi_barrier(glcomm,ierr)
        if (lwp) write(*,*) "Ghosh_Init: before Ens_shared_alloc_base" 
    end if
    
    n_NewBase=n_DAstate
    if (.false.) then
        call Ens_shared_alloc_base(n_NewBase, member_pointer, global_pointer, win_NewBase)
        if (EnsRank/=EnsRankZero) then
            NewBase(1:nk_DAstate, 1:nj_DAstate, 1:ni_DAstate, 1:ntra_DAstate)=>member_pointer
            NewBase  = huge(NewBase(1,1,1,1))
        end if
        gl_NewBase(1:nk_DAstate, 1:nj_DAstate, 1:ni_DAstate, 1:ntra_DAstate, 1:EnsSize-1)=>global_pointer
    else
        call Ens_shared_alloc(n_NewBase, member_pointer, global_pointer, win_NewBase)
        NewBase(1:nk_DAstate, 1:nj_DAstate, 1:ni_DAstate, 1:ntra_DAstate)=>member_pointer
        NewBase  = huge(NewBase(1,1,1,1))
        gl_NewBase(1:nk_DAstate, 1:nj_DAstate, 1:ni_DAstate, 1:ntra_DAstate, 0:EnsSize-1)=>global_pointer
    end if
    
    Select Case (GhoshOrder)
        
        Case (2,3)
            
        Case (5)
            !EnsWeights=gnappo!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
        Case default
            if (lwp) write(*,*) 'GhoshOrder = ', GhoshOrder, ' not supported yet. Aborting.'
            call mpi_barrier(glcomm,ierr)
            call MPI_abort(glcomm, 1, ierr)
        
    end select
    
    
    if (EnsDebug>1) then
        call mpi_barrier(glcomm,ierr)
        if (lwp) write(*,*) "Ghosh_Init done" 
    end if
    
end subroutine

subroutine Ghosh_Finalize
    integer ierror
    
    CALL MPI_Win_free(win_HLTR1HLi_loc, ierror)
    
    deallocate(EnsCov1)
    deallocate(ForecastCov1)
    deallocate(OrtMatrixSampling)
    
    if (LocalAnalysis) then
        call Ens_Finalize_Local
        call LocalPatch%Free
    end if
    
    call Finalize_Algebra
    
    CALL MPI_Win_free(win_ChangeBase, ierror)
    CALL MPI_Win_free(win_NewBase, ierror)
    
    deallocate(TTW1T)
    
end subroutine

subroutine Ghosh_Forecast
    !call Ens_ReduceMeanAndBase(........) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
end subroutine

subroutine Ghosh_Analysis
    
    double precision, dimension(:), POINTER, contiguous :: member_pointer
    double precision, dimension(:,:), POINTER, contiguous :: global_pointer
    
    integer :: indexi, indexj, indexk, indext, indexd, c, rankc
    integer ierror
    
    ForecastCov1=TTW1T*ForgettingFactor
    
    n_Hstate=nj_DAstate*ni_DAstate*nObs
    call Ens_shared_alloc(n_Hstate, member_pointer, global_pointer, win_Hstate)
    Hstate(1:nObs, 1:nj_DAstate, 1:ni_DAstate)=>member_pointer
    gl_Hstate(1:n_Hstate, 0:EnsSize-1)=>global_pointer
    Hstate  = huge(Hstate(1,1,1))
    
    n_SqrtR1HLi=nj_DAstate*ni_DAstate*n_SqrtR1
    call Ens_shared_alloc(n_SqrtR1HLi, member_pointer, global_pointer, win_SqrtR1HLi)
    SqrtR1HLi(1:n_SqrtR1, 1:nj_DAstate, 1:ni_DAstate)=>member_pointer
    gl_SqrtR1HLi(1:n_SqrtR1, 1:nj_DAstate, 1:ni_DAstate, 0:EnsSize-1)=>global_pointer
    SqrtR1HLi  = huge(SqrtR1HLi(1,1,1))
    
    Hstate=H(DAstate_kjit)
    
    call Ens_ReduceMeanAndBase(win_Hstate, n_Hstate, gl_Hstate)
    
    call Ens_ReduceMeanAndBase(win_DAstate, n_DAstate, gl_DAstate)
    
    if (EnsRank==EnsRankZero) then
        Hstate=H(DAstate_kjit)
        Hstate=Misfit(Hstate)
    end if
    
    SqrtR1HLi=SqrtR1(Hstate)
    
    if (LocalAnalysis) then
    
        CALL MPI_Win_fence(0, win_SqrtR1HLi, ierror)
            HLTR1HLi_loc=0.0d0
            c=0
            do indexd=0, EnsSize-1
                if (indexd/=EnsRankZero) then
                    c=c+1
                    do indexi=1,ni_DAstate
                        do indexj=1,nj_DAstate
                            HLTR1HLi_loc(c,indexj,indexi)=HLTR1HLi_loc(c,indexj,indexi)+DOT_PRODUCT(SqrtR1HLi(:, indexj, indexi),gl_SqrtR1HLi(:, indexj, indexi, indexd))
                        end do
                    end do
                end if
            end do
        CALL MPI_Win_fence(0, win_SqrtR1HLi, ierror)
        
        call LocalPatch%ComputeAll(HLTR1HLi_loc, UseLocalObsDumping)
        
        CALL MPI_Win_fence(0, win_HLTR1HLi_loc, ierror)
        
            if (EnsDebug>0.and.EnsRank==EnsRankZero) then
                if (myrank==myrankZero) then
                    write(*,*) '-------------------------------------------------------'
                    write(*,*) '-------------------------------------------------------'
                    write(*,*) ''
                    write(*,*) 'DEBUG: parameters estimation (332, 284).'
                    write(*,*) ''
                    write(*,*) '-------------------------------------------------------'
                    write(*,*) '-------------------------------------------------------'
                end if
                do indexj=0, mysize-1
                    call mpi_barrier(mycomm, ierror)
                    if (myrank==indexj) then
                        write(*,*) ''
                        write(*,*) 'myrank: ', myrank
                        write(*,*) ''
                        write(*,*) '-------------------------------------------------------'
!                             write(*,*) ''
!                             write(*,*) 'gl_DAstate_kjitn'
!                             write(*,*) ''
!                             do indexi=0,EnsSize-1
!                                 write(*,*) gl_DAstate_kjitn(1, 14, 28, (/14,19,23,27/), indexi)
!                             end do
!                             write(*,*) ''
!                             write(*,*) '-------------------------------------------------------'
!                             write(*,*) ''
!                             write(*,*) 'gl_Hstate'
!                             write(*,*) ''
!                             do indexi=0,EnsSize-1
!                                 write(*,*) gl_Hstate(14+27*nj_DAstate,indexi)
!                             end do
!                             write(*,*) ''
!                             write(*,*) '-------------------------------------------------------'
!                             write(*,*) ''
!                             write(*,*) 'ForecastCov1'
!                             write(*,*) ''
!                             do indexi=1,EnsSize-1
!                                 write(*,*) ForecastCov1(:,indexi)
!                             end do
!                             write(*,*) ''
!                             write(*,*) '-------------------------------------------------------'
                        write(*,*) ''
                        write(*,*) 'gl_HLTR1HLi_loc diagonal'
                        write(*,*) ''
                        EnsCov1(:,1:EnsRankZero)=gl_HLTR1HLi_loc(:, 14, 28, 0:EnsRankZero-1)
                        EnsCov1(:,EnsRankZero+1:EnsSize-1)=gl_HLTR1HLi_loc(:, 14, 28, EnsRankZero+1:EnsSize-1)
                        do indexi=1,EnsSize-3,3
                            write(*,*) EnsCov1(indexi, indexi), EnsCov1(indexi+1, indexi+1), EnsCov1(indexi+2, indexi+2)
                        end do
                        if (mod(EnsSize-1,3)==1) then
                            write(*,*) EnsCov1(EnsSize-1, EnsSize-1)
                        elseif (mod(EnsSize-1,3)==2) then
                            write(*,*) EnsCov1(EnsSize-2, EnsSize-2), EnsCov1(EnsSize-1, EnsSize-1)
                        endif
                        write(*,*) ''
                        write(*,*) '-------------------------------------------------------'
                        write(*,*) ''
                        write(*,*) 'determinant gl_HLTR1HLi_loc: ', det(EnsSize-1,EnsCov1)
                        write(*,*) ''
                        write(*,*) '-------------------------------------------------------'
                        write(*,*) ''
                        EnsCov1=ForecastCov1
                        write(*,*) 'determinant ForecastCov1: ', det(EnsSize-1,EnsCov1)
                        write(*,*) ''
                        write(*,*) ''
                        write(*,*) '-------------------------------------------------------'
                        write(*,*) ''
                        EnsCov1(:,1:EnsRankZero)=ForecastCov1(:,1:EnsRankZero)+gl_HLTR1HLi_loc(:, 14, 28, 0:EnsRankZero-1)
                        EnsCov1(:,EnsRankZero+1:EnsSize-1)=ForecastCov1(:,EnsRankZero+1:EnsSize-1)+gl_HLTR1HLi_loc(:, 14, 28, EnsRankZero+1:EnsSize-1)
                        write(*,*) 'determinant EnsCov1: ', det(EnsSize-1,EnsCov1)
                        write(*,*) ''
                        write(*,*) ''
                        write(*,*) '-------------------------------------------------------'
                        write(*,*) ''
                        write(*,*) 'Misfit'
                        write(*,*) ''
                        write(*,*) gl_HLTR1HLi_loc(:, 14, 28, EnsRankZero)
                        write(*,*) ''
                        write(*,*) '-------------------------------------------------------'
                        write(*,*) '-------------------------------------------------------'
                    end if
                end do
            end if
            
            if (EnsDebug>1) then
                call mpi_barrier(glcomm,ierror)
                if (lwp) write(*,*) "Ghosh analysis: before computing analysis cov" 
            end if
            
            rankc=-1
            do indexi=1,ni_DAstate
                do indexj=1,nj_DAstate
                    if (DAMask(1, indexj,indexi)==1) then
                        rankc=rankc+1
                        
                        if (Mod(rankc, EnsSize)==EnsRank) then
                        
                            EnsCov1(:,1:EnsRankZero)=ForecastCov1(:,1:EnsRankZero)+gl_HLTR1HLi_loc(:, indexj, indexi, 0:EnsRankZero-1)
                            EnsCov1(:,EnsRankZero+1:EnsSize-1)=ForecastCov1(:,EnsRankZero+1:EnsSize-1)+gl_HLTR1HLi_loc(:, indexj, indexi, EnsRankZero+1:EnsSize-1)
                            
!                                 if (EnsDebug>1.and.indexi==28.and.indexj==14) then
!                                     write(*,*) ''
!                                     write(*,*) 'EnsCov1'
!                                     write(*,*) ''
!                                     do indexd=1,EnsSize-1
!                                         write(*,*) EnsCov1(:, indexd)
!                                     end do
!                                     write(*,*) ''
!                                     write(*,*) '-------------------------------------------------------'
!                                     write(*,*) '-------------------------------------------------------'
!                                 end if                                
                            
                            call SymChangeBase(EnsCov1,ierror)
                            if (ierror/=0) then
                                write(*,*) 'myrank: ', myrank, ', EnsRank: ', EnsRank, ', i_DA: ', indexi, ', j_DA: ', indexj
                                call MPI_abort(glcomm, 1, ierror)
                            end if
                            
!                                 if (EnsDebug>1.and.indexi==28.and.indexj==14) then
!                                     write(*,*) ''
!                                     write(*,*) 'EnsCov1 SymChangeBase'
!                                     write(*,*) ''
!                                     do indexd=1,EnsSize-1
!                                         write(*,*) EnsCov1(:, indexd)
!                                     end do
!                                     write(*,*) ''
!                                     write(*,*) '-------------------------------------------------------'
!                                     write(*,*) '-------------------------------------------------------'
!                                 end if        
                            
                            gl_HLTR1HLi_loc(:, indexj, indexi, 0:EnsRankZero-1)=EnsCov1(:,1:EnsRankZero)
                            gl_HLTR1HLi_loc(:, indexj, indexi, EnsRankZero)=matmul(matmul(gl_HLTR1HLi_loc(:, indexj, indexi, EnsRankZero),EnsCov1),EnsCov1)
                            gl_HLTR1HLi_loc(:, indexj, indexi, EnsRankZero+1:EnsSize-1)=EnsCov1(:,EnsRankZero+1:EnsSize-1)
                            
                        end if                        
                    end if
                end do
            end do
        CALL MPI_Win_fence(0, win_HLTR1HLi_loc, ierror)   
        
        if (EnsDebug>0.and.EnsRank==EnsRankZero) then
            do indexj=0, mysize-1
                call mpi_barrier(mycomm, ierror)
                if (myrank==indexj) then
                    write(*,*) ''
                    write(*,*) 'myrank: ', myrank
                    write(*,*) ''
!                         write(*,*) '-------------------------------------------------------'
!                         write(*,*) ''
!                         write(*,*) 'gl_HLTR1HLi_loc changebase'
!                         write(*,*) ''
!                         do indexi=0,EnsRankZero-1
!                             write(*,*) gl_HLTR1HLi_loc(:, 14, 28, indexi)
!                         end do
!                         do indexi=EnsRankZero+1,EnsSize-1
!                             write(*,*) gl_HLTR1HLi_loc(:, 14, 28, indexi)
!                         end do
!                         write(*,*) ''
                    write(*,*) '-------------------------------------------------------'
                    write(*,*) ''
                    write(*,*) 'analysis coeff oldbase'
                    write(*,*) ''
                    write(*,*) gl_HLTR1HLi_loc(:, 14, 28, EnsRankZero)
                    write(*,*) ''
                    write(*,*) '-------------------------------------------------------'
                    write(*,*) ''
                    write(*,*) 'sum: ', sum(gl_HLTR1HLi_loc(:, 14, 28, EnsRankZero))
                    write(*,*) ''
                    write(*,*) '-------------------------------------------------------'
                    write(*,*) ''
                    write(*,*) 'abssum: ', sum(abs(gl_HLTR1HLi_loc(:, 14, 28, EnsRankZero)))
                    write(*,*) ''
                    write(*,*) '-------------------------------------------------------'
                    write(*,*) '-------------------------------------------------------'
                end if
            end do
        end if
        
        if (EnsDebug>1) then
            call mpi_barrier(glcomm,ierror)
            if (lwp) write(*,*) "Ghosh analysis: before changing base" 
        end if
        
        CALL MPI_Win_fence(0, win_DAstate, ierror)
        CALL MPI_Win_fence(0, win_ChangeBase, ierror)
            if (EnsRank==EnsRankZero) then
                c=0
                do indexd=0,EnsSize-1
                    if (indexd/=EnsRankZero) then
                        c=c+1
                        do indext=1, ntra_DAstate
                            do indexi=1,ni_DAstate
                                do indexj=1, nj_DAstate
                                    do indexk=1, nk_DAstate
                                        if (DAMask(indexk, indexj,indexi)==1) then
                                            DAstate_kjit(indexk, indexj, indexi, indext)=DAstate_kjit(indexk, indexj, indexi, indext) + &
                                                gl_DAstate_kjitn(indexk, indexj, indexi, indext, indexd) * HLTR1HLi_loc(c, indexj, indexi)
                                        else
                                            exit
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end if
                end do                
            
                EnsCov1=0.0d0
                do indexd=1,EnsSize-1
                    EnsCov1(indexd,indexd)=1.0d0
                end do
                
                CALL MPI_Win_fence(0, win_DAstate, ierror)
                
                if (EnsDebug>1) then
                    call mpi_barrier(glcomm,ierror)
                    if (lwp) write(*,*) "before 2nd order sampling" 
                end if
                
                if (myrank==myrankZero) then
                    call SecondOrdExactSampling(gl_ChangeBase) 
                    if (EnsDebug>1) write(*,*) 'computed 2nd order sampling transformation'
                end if
                call MPI_Bcast(gl_ChangeBase, EnsSize*(EnsSize-1), mpi_real8, 0, mycomm, ierror)
                
                if (EnsDebug>1) then
                    if (myrank==0) write(*,*) "broadcasted transformation" 
                end if
                
            else
                NewBase=0.0d0
                c=0
                do indexd=0,EnsSize-1
                    if (indexd/=EnsRankZero) then
                        c=c+1
                        do indext=1, ntra_DAstate
                            do indexi=1,ni_DAstate
                                do indexj=1, nj_DAstate
                                    do indexk=1, nk_DAstate
                                        if (DAMask(indexk, indexj,indexi)==1) then
                                            NewBase(indexk, indexj, indexi, indext)=NewBase(indexk, indexj, indexi, indext) + &
                                                gl_DAstate_kjitn(indexk, indexj, indexi, indext, indexd) * HLTR1HLi_loc(c, indexj, indexi)
                                        else
                                            exit
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end if
                end do
                
                CALL MPI_Win_fence(0, win_DAstate, ierror)
                
                if (EnsDebug>1) then
                    call mpi_barrier(glcomm,ierror)
                    if (lwp) write(*,*) "before 2nd order sampling" 
                end if
                
                do indext=1, ntra_DAstate
                    do indexi=1,ni_DAstate
                        do indexj=1, nj_DAstate
                            do indexk=1, nk_DAstate
                                if (DAMask(indexk, indexj,indexi)==0) exit
                                DAstate_kjit(indexk, indexj,indexi, indext) = gl_DAstate_kjitn(indexk, indexj,indexi, indext, EnsRankZero) 
                            end do
                        end do
                    end do
                end do
                
                if (EnsDebug>1) then
                    call mpi_barrier(mycomm,ierror)
                    if (myrank==0) write(*,*) "EnsRank=", EnsRank, ', DAstate_kjit updated' 
                end if
                
            end if
        CALL MPI_Win_fence(0, win_ChangeBase, ierror)
        CALL MPI_Win_fence(0, win_DAstate, ierror)
        
        if (EnsDebug>1) then
            call mpi_barrier(glcomm,ierror)
            if (lwp) write(*,*) "Ghosh analysis: before computing ensemble" 
        end if
        
        CALL MPI_Win_fence(0, win_NewBase, ierror)
        
            if (EnsDebug>1.and.EnsRank==EnsRankZero) then
                do indexj=0, mysize-1
                    call mpi_barrier(mycomm, ierror)
                    if (myrank==indexj) then
                        write(*,*) ''
                        write(*,*) 'myrank: ', myrank
                        write(*,*) ''
                        write(*,*) '-------------------------------------------------------'
                        write(*,*) ''
                        write(*,*) 'gl_NewBase'
                        write(*,*) ''
                        do indexi=0,EnsRankZero-1
                            write(*,*) gl_NewBase(1, 14, 28, (/14,19,23,27/), indexi)
                        end do
                        do indexi=EnsRankZero+1,EnsSize-1
                            write(*,*) gl_NewBase(1, 14, 28, (/14,19,23,27/), indexi)
                        end do
                        write(*,*) ''
                        write(*,*) '-------------------------------------------------------'
                        write(*,*) ''
                        write(*,*) 'DAstate_kjit'
                        write(*,*) ''
                        write(*,*) DAstate_kjit(1, 14, 28, (/14,19,23,27/))
                        write(*,*) ''
                        write(*,*) '-------------------------------------------------------'
                        write(*,*) '-------------------------------------------------------'
                    end if
                end do
            end if
        
            !do indexd=1,EnsSize-1
            c=0
            do indexd=0,EnsSize-1
                if (indexd/=EnsRankZero) then
                    c=c+1
                do indext=1, ntra_DAstate
                    do indexi=1,ni_DAstate
                        do indexj=1, nj_DAstate
                            do indexk=1, nk_DAstate
                                if (DAMask(indexk, indexj,indexi)==1) then
                                    DAstate_kjit(indexk, indexj, indexi, indext)=DAstate_kjit(indexk, indexj, indexi, indext) + &
                                        !gl_NewBase(indexk, indexj, indexi, indext, indexd) * ChangeBase(indexd)
                                        gl_NewBase(indexk, indexj, indexi, indext, indexd) * ChangeBase(c)
                                else
                                    exit
                                end if
                            end do
                        end do
                    end do
                end do
                end if
            end do            
        CALL MPI_Win_fence(0, win_NewBase, ierror)
        
    else
    
        if (lwp) write(*,*) 'Non-local analysis is not implemented at this moment. Aborting.'
        call MPI_abort(glcomm, 1, ierror)
        
    end if
    
    CALL MPI_Win_free(win_Hstate, ierror)
    CALL MPI_Win_free(win_SqrtR1HLi, ierror)
    
end subroutine

subroutine SecondOrdExactSampling(ChangeBaseMatrix)
    
    double precision, dimension(EnsSize-1,EnsSize), intent(out) :: ChangeBaseMatrix
    
    double precision :: AllWeightsSqrt, AllWeightsSqrt1
    integer indexi
    integer ierror
    
    AllWeightsSqrt1=sqrt(dble(EnsSize))
    AllWeightsSqrt=1.0d0/AllWeightsSqrt1
    
    call SymChangeBase(EnsCov1,ierror)
    if (ierror/=0) then
        write(*,*) 'SecondOrdExactSampling SymChangeBase failed. Aborting.'
        call MPI_abort(glcomm, 1, ierror)
    end if
    
    OrtMatrixSampling(:,1)=AllWeightsSqrt
    call OrtMatrix(OrtMatrixSampling, EnsSize, EnsSize, 1) 
    do indexi=2, EnsSize
        OrtMatrixSampling(:,indexi)=OrtMatrixSampling(:,indexi)*AllWeightsSqrt1
    end do
    
    ChangeBaseMatrix=transpose(OrtMatrixSampling(:,2:EnsSize))        
    ChangeBaseMatrix=MatMul(EnsCov1,ChangeBaseMatrix)

end subroutine

subroutine TTW1T_builder

    double precision, dimension(:,:), allocatable :: matrixT
    integer :: indexi
    
    if (.false.) then ! true if you want to build T matrix first, false if you want to make it faster.
        
        allocate(matrixT(0:EnsSize-1,EnsSize-1))
        
        do indexi=1, EnsRankZero
            matrixT(:,indexi)=-EnsWeights
            matrixT(indexi-1,indexi)=matrixT(indexi-1,indexi)+1.0d0
        end do
        do indexi=EnsRankZero+1, EnsSize-1
            matrixT(:,indexi)=-EnsWeights
            matrixT(indexi,indexi)=matrixT(indexi,indexi)+1.0d0
        end do
        do indexi=1, EnsSize-1
            TTW1T(:,indexi)=matmul(matrixT(:,indexi)/EnsWeights,matrixT)
        end do
        
        deallocate(matrixT)
        
    else
    
        do indexi=1, EnsRankZero
            TTW1T(:,indexi)=-1.0d0
            TTW1T(indexi,indexi)=1.0d0/EnsWeights(indexi-1)-1.0d0
        end do
        do indexi=EnsRankZero+1, EnsSize-1
            TTW1T(:,indexi)=-1.0d0
            TTW1T(indexi,indexi)=1.0d0/EnsWeights(indexi)-1.0d0
        end do
    
    end if
    
end subroutine

end module




