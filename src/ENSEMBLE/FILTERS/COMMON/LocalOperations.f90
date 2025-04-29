! developed by Simone Spada (sspada@ogs.it) at OGS

module LocalOperations
    use mpi
    
    use myalloc, &
        only: lwp, myrank, mycomm, &
            domdec, nono, noso, noea, nowe
    use Ens_Mem, &
        only: EnsSize, &
            Namelist_LocalRange => LocalRange
    use Ens_Custom, &   
        only: nj_DAstate, ni_DAstate, DAMask, &
            LocalWeightFunction
        
    implicit none
    
    integer xRank, yRank, nonoea, nonowe, nosoea, nosowe, LocalMaxRange
    
    type :: LocalSpace
        integer EnsDim
        integer LocalRange
        double precision, dimension(:,:,:), allocatable :: LocalPatch, VerticalPatch, HorizontalPatch, DiagonalPatch
    contains
        procedure :: Init => Local_alloc
        procedure :: Free => Local_dealloc
        procedure :: SendRecive => LocalSendRecive
        procedure :: Sum => SummingLocalPatch
        procedure :: WeightedSum => SummingLocalPatchWeighted
        procedure :: ComputeAll => LocalComputeAll
    end type
    
contains

subroutine Ens_Init_Local
    
    xRank=domdec(myrank+1, 2)
    yRank=domdec(myrank+1, 3)
    
    if (noea==-1) then
        if (nono==-1) then 
            nonoea=-1
        else
            nonoea=domdec(nono+1, 11)
        end if
        
        if (noso==-1) then 
            nosoea=-1
        else
            nosoea=domdec(noso+1, 11)
        end if
    else
        nonoea=domdec(noea+1, 12)
        nosoea=domdec(noea+1, 13)
    end if
    
    if (nowe==-1) then
        if (nono==-1) then 
            nonowe=-1
        else
            nonowe=domdec(nono+1, 10)
        end if
        
        if (noso==-1) then 
            nosowe=-1
        else
            nosowe=domdec(noso+1, 10)
        end if
    else
        nonowe=domdec(nowe+1, 12)
        nosowe=domdec(nowe+1, 13)
    end if
    
    LocalMaxRange=minval(domdec(:,4:5))-2
    if (Namelist_LocalRange>LocalMaxRange) then
        if (lwp) write(*,*) 'WARNING! The local range (', Namelist_LocalRange, ') is bigger than the maximum local range (', LocalMaxRange, '). LocalRange reduced to LocalMaxRange.'
        Namelist_LocalRange=LocalMaxRange
    end if
    
end subroutine

subroutine Ens_Finalize_Local
    
end subroutine

subroutine Local_alloc(this, EnsDim, Optional_LocalRange)

    class(LocalSpace), intent(inout) :: this
    integer, intent(in) :: EnsDim
    integer, optional, intent(in) :: Optional_LocalRange
    
    integer :: LocalRange
    
    this%EnsDim = EnsDim
    
    if (present(Optional_LocalRange)) then
        LocalRange = Optional_LocalRange
        if (LocalRange>LocalMaxRange) then
            if (lwp) write(*,*) 'WARNING! The local range (', LocalRange, ') is bigger than the maximum local range (', LocalMaxRange, '). LocalRange reduced to LocalMaxRange.'
            LocalRange=LocalMaxRange
        end if
    else
        LocalRange = Namelist_LocalRange
    end if
    
    this%LocalRange=LocalRange
    
    if (LocalRange>0) then
    
        allocate(this%LocalPatch(this%EnsDim, 1-LocalRange:nj_DAstate+LocalRange, 1-LocalRange:ni_DAstate+LocalRange))
        this%LocalPatch=huge(this%LocalPatch(1,1,1))
        
        allocate(this%VerticalPatch(this%EnsDim, LocalRange, ni_DAstate))
        this%VerticalPatch=Huge(this%VerticalPatch(1,1,1))
        
        allocate(this%HorizontalPatch(this%EnsDim, nj_DAstate, LocalRange))
        this%HorizontalPatch=Huge(this%HorizontalPatch(1,1,1))
        
        allocate(this%DiagonalPatch(this%EnsDim, LocalRange, LocalRange))
        this%DiagonalPatch=Huge(this%DiagonalPatch(1,1,1))
    
    end if
    
end subroutine

subroutine Local_dealloc(this)
    class(LocalSpace), intent(inout) :: this
    
    if (this%LocalRange>0) then
    
        deallocate(this%LocalPatch)
        deallocate(this%VerticalPatch)
        deallocate(this%HorizontalPatch)
        deallocate(this%DiagonalPatch)
    
    end if
    
end subroutine

subroutine LocalComputeAll(this, patch, WeightedObs)

    class(LocalSpace), intent(inout) :: this
    double precision, dimension(this%EnsDim,nj_DAstate, ni_DAstate), intent(inout) :: patch
    logical, intent(in) :: WeightedObs
    
    if (this%LocalRange>0) then
        
            call this%SendRecive(patch)
    
            if (WeightedObs) then
                call this%WeightedSum(patch)
            else
                call this%Sum(patch)
            end if
        
        end if
    
end subroutine

subroutine LocalSendRecive(this, patch, opt_FillValue)
    
    class(LocalSpace), intent(inout) :: this
    double precision, dimension(this%EnsDim,nj_DAstate, ni_DAstate), intent(in) :: patch
    double precision, optional, intent(in) :: opt_FillValue
    
    integer :: ierr, LocalRange
    double precision :: FillValue
    
!               3     4
!               |     ^
!               |     |
!               v     |
!           ________________
!          |                |
!     1<-- |                | 1 <--
!     2--> |                | 2 -->
!          |________________|
!               3     4
!               |     ^
!               |     |
!               v     |    
    
    if (present(opt_FillValue)) then
        FillValue=opt_FillValue
    else
        FillValue=0.0d0
    end if
    
    LocalRange=this%LocalRange
    
    this%LocalPatch=FillValue
    this%LocalPatch(:,1:nj_DAstate,1:ni_DAstate)=patch

    if (mod(xRank, 2)==0) then
        if (nowe/=-1) then
            this%HorizontalPatch=patch(:,:,1:LocalRange)            
            call MPI_Send(this%HorizontalPatch, nj_DAstate*LocalRange*this%EnsDim, MPI_real8, nowe, 10, mycomm, ierr)
            call MPI_Recv(this%HorizontalPatch, nj_DAstate*LocalRange*this%EnsDim, MPI_real8, nowe, 20, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1:nj_DAstate, 1-LocalRange:0)=this%HorizontalPatch
        end if
        
        if (noea/=-1) then
            this%HorizontalPatch=patch(:,:,ni_DAstate+1-LocalRange:ni_DAstate)
            call MPI_Send(this%HorizontalPatch, nj_DAstate*LocalRange*this%EnsDim, MPI_real8, noea, 20, mycomm, ierr)
            call MPI_Recv(this%HorizontalPatch, nj_DAstate*LocalRange*this%EnsDim, MPI_real8, noea, 10, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1:nj_DAstate, ni_DAstate+1:ni_DAstate+LocalRange)=this%HorizontalPatch
        end if
    else
        if (noea/=-1) then     
            call MPI_Recv(this%HorizontalPatch, nj_DAstate*LocalRange*this%EnsDim, MPI_real8, noea, 10, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1:nj_DAstate, ni_DAstate+1:ni_DAstate+LocalRange)=this%HorizontalPatch
            this%HorizontalPatch=patch(:,:,ni_DAstate+1-LocalRange:ni_DAstate)
            call MPI_Send(this%HorizontalPatch, nj_DAstate*LocalRange*this%EnsDim, MPI_real8, noea, 20, mycomm, ierr)
        end if
    
        if (nowe/=-1) then
            call MPI_Recv(this%HorizontalPatch, nj_DAstate*LocalRange*this%EnsDim, MPI_real8, nowe, 20, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1:nj_DAstate, 1-LocalRange:0)=this%HorizontalPatch
            this%HorizontalPatch=patch(:,:,1:LocalRange)
            call MPI_Send(this%HorizontalPatch, nj_DAstate*LocalRange*this%EnsDim, MPI_real8, nowe, 10, mycomm, ierr)   
        end if
    end if
    
    if (mod(yRank, 2)==0) then
        if (noso/=-1) then
            this%VerticalPatch=patch(:,1:LocalRange,:)
            call MPI_Send(this%VerticalPatch, LocalRange*ni_DAstate*this%EnsDim, MPI_real8, noso, 30, mycomm, ierr)
            call MPI_Recv(this%VerticalPatch, LocalRange*ni_DAstate*this%EnsDim, MPI_real8, noso, 40, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1-LocalRange:0,1:ni_DAstate)=this%VerticalPatch
        end if
        
        if (nono/=-1) then
            this%VerticalPatch=patch(:,nj_DAstate+1-LocalRange:nj_DAstate,:)
            call MPI_Send(this%VerticalPatch, LocalRange*ni_DAstate*this%EnsDim, MPI_real8, nono, 40, mycomm, ierr)
            call MPI_Recv(this%VerticalPatch, LocalRange*ni_DAstate*this%EnsDim, MPI_real8, nono, 30, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,nj_DAstate+1:nj_DAstate+LocalRange,1:ni_DAstate)=this%VerticalPatch
        end if
    else
        if (nono/=-1) then
            call MPI_Recv(this%VerticalPatch, LocalRange*ni_DAstate*this%EnsDim, MPI_real8, nono, 30, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,nj_DAstate+1:nj_DAstate+LocalRange,1:ni_DAstate)=this%VerticalPatch
            this%VerticalPatch=patch(:,nj_DAstate+1-LocalRange:nj_DAstate,:)
            call MPI_Send(this%VerticalPatch, LocalRange*ni_DAstate*this%EnsDim, MPI_real8, nono, 40, mycomm, ierr)
        end if
    
        if (noso/=-1) then
            call MPI_Recv(this%VerticalPatch, LocalRange*ni_DAstate*this%EnsDim, MPI_real8, noso, 40, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1-LocalRange:0,1:ni_DAstate)=this%VerticalPatch
            this%VerticalPatch=patch(:,1:LocalRange,:)
            call MPI_Send(this%VerticalPatch, LocalRange*ni_DAstate*this%EnsDim, MPI_real8, noso, 30, mycomm, ierr)
        end if
    end if
    
    if (mod((xRank+yRank)/2, 2)==0) then
        if (nosowe/=-1) then
            this%DiagonalPatch=patch(:,1:LocalRange,1:LocalRange)
            call MPI_Send(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nosowe, 50, mycomm, ierr)
            call MPI_Recv(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nosowe, 60, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1-LocalRange:0, 1-LocalRange:0)=this%DiagonalPatch
        end if
        
        if (nonoea/=-1) then
            this%DiagonalPatch=patch(:,nj_DAstate+1-LocalRange:nj_DAstate,ni_DAstate+1-LocalRange:ni_DAstate)
            call MPI_Send(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nonoea, 60, mycomm, ierr)
            call MPI_Recv(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nonoea, 50, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,nj_DAstate+1:nj_DAstate+LocalRange, ni_DAstate+1:ni_DAstate+LocalRange)=this%DiagonalPatch
        end if
    else
        if (nonoea/=-1) then
            call MPI_Recv(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nonoea, 50, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,nj_DAstate+1:nj_DAstate+LocalRange, ni_DAstate+1:ni_DAstate+LocalRange)=this%DiagonalPatch
            this%DiagonalPatch=patch(:,nj_DAstate+1-LocalRange:nj_DAstate,ni_DAstate+1-LocalRange:ni_DAstate)
            call MPI_Send(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nonoea, 60, mycomm, ierr)
        end if
    
        if (nosowe/=-1) then
            call MPI_Recv(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nosowe, 60, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1-LocalRange:0, 1-LocalRange:0)=this%DiagonalPatch
            this%DiagonalPatch=patch(:,1:LocalRange,1:LocalRange)
            call MPI_Send(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nosowe, 50, mycomm, ierr)
        end if
    end if
    
    if (mod(floor((xRank-yRank)/2.0), 2)==0) then
        if (nosoea/=-1) then
            this%DiagonalPatch=patch(:,1:LocalRange,ni_DAstate+1-LocalRange:ni_DAstate)
            call MPI_Send(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nosoea, 70, mycomm, ierr)
            call MPI_Recv(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nosoea, 80, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1-LocalRange:0, ni_DAstate+1:ni_DAstate+LocalRange)=this%DiagonalPatch
        end if

        if (nonowe/=-1) then
            this%DiagonalPatch=patch(:,nj_DAstate+1-LocalRange:nj_DAstate,1:LocalRange)
            call MPI_Send(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nonowe, 80, mycomm, ierr)
            call MPI_Recv(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nonowe, 70, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,nj_DAstate+1:nj_DAstate+LocalRange, 1-LocalRange:0)=this%DiagonalPatch
        end if
    else
        if (nonowe/=-1) then
            call MPI_Recv(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nonowe, 70, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,nj_DAstate+1:nj_DAstate+LocalRange, 1-LocalRange:0)=this%DiagonalPatch
            this%DiagonalPatch=patch(:,nj_DAstate+1-LocalRange:nj_DAstate,1:LocalRange)
            call MPI_Send(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nonowe, 80, mycomm, ierr)
        end if
    
        if (nosoea/=-1) then
            call MPI_Recv(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nosoea, 80, mycomm, MPI_STATUS_IGNORE, ierr)
            this%LocalPatch(:,1-LocalRange:0, ni_DAstate+1:ni_DAstate+LocalRange)=this%DiagonalPatch
            this%DiagonalPatch=patch(:,1:LocalRange,ni_DAstate+1-LocalRange:ni_DAstate)
            call MPI_Send(this%DiagonalPatch, LocalRange*LocalRange*this%EnsDim, MPI_real8, nosoea, 70, mycomm, ierr)
        end if
    end if

end subroutine

subroutine SummingLocalPatch(this, patch)
    
    class(LocalSpace), intent(inout) :: this
    double precision, dimension(this%EnsDim, nj_DAstate, ni_DAstate), intent(out) :: patch
    
    integer :: indexi, indexj, indexk, temp, LocalRange
    
    LocalRange=this%LocalRange
    
    patch=0.0d0
    do indexi=-LocalRange,LocalRange
        temp=floor(sqrt((0.5d0+LocalRange)**2-indexi**2))
        do indexj=-temp+1, temp+1
            patch(:,1,1)=patch(:,1,1)+this%LocalPatch(:,indexj, 1+indexi)
        end do
    end do
    do indexj=2, nj_DAstate
        patch(:,indexj,1)=patch(:,indexj-1,1)
        do indexi=-LocalRange,LocalRange
            temp=floor(sqrt((0.5d0+LocalRange)**2-indexi**2))
            patch(:,indexj,1)=patch(:,indexj,1) + &
                this%LocalPatch(:,indexj+temp, 1+indexi) - this%LocalPatch(:,indexj-1-temp, 1+indexi)
        end do
    end do
    do indexi=2, ni_DAstate
        do indexj=1,nj_DAstate
            patch(:,indexj,indexi)=patch(:,indexj,indexi-1)
            do indexk=-LocalRange,LocalRange
                temp=floor(sqrt((0.5d0+LocalRange)**2-indexk**2))
                patch(:,indexj,indexi)=patch(:,indexj,indexi) + &
                    this%LocalPatch(:,indexj+indexk, indexi+temp) - this%LocalPatch(:,indexj+indexk, indexi-1-temp)
            end do
        end do
    end do
    
end subroutine

subroutine SummingLocalPatchWeighted(this, patch)
    
    class(LocalSpace), intent(inout) :: this
    double precision, dimension(this%EnsDim, nj_DAstate, ni_DAstate), intent(out) :: patch
    
    integer :: indexi, indexj, indexk, indexl, temp, LocalRange
    
    LocalRange=this%LocalRange
    
    patch=0.0d0
    do indexi=1, ni_DAstate
        do indexj=1,nj_DAstate
            if (DAMask(1, indexj,indexi)==1) then
                do indexk=-LocalRange,LocalRange
                    temp=floor(sqrt((0.5d0+LocalRange)**2-indexk**2))
                    do indexl=-temp,temp
                        patch(:,indexj,indexi)=patch(:,indexj,indexi) + &
                            this%LocalPatch(:,indexj+indexl,indexi+indexk) * LocalWeightFunction(dble(indexk**2+ indexl**2))
                    end do
                end do
            end if
        end do
    end do
    
end subroutine
    
end module
