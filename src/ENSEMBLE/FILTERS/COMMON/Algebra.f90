! developed by Simone Spada (sspada@ogs.it) at OGS

module Algebra
    use Ens_Mem, &
        only: EnsSize
        
    implicit none
    
    integer EnsDim, lwork, liwork
    integer, dimension(:), allocatable :: isuppz, iwork
    double precision, dimension(:,:), allocatable :: eigenvectors
    double precision, dimension(:), allocatable :: eigenvalues, work
    
contains

    subroutine Init_Algebra
        
        EnsDim=EnsSize-1
        
        allocate(isuppz(2*EnsDim))
        isuppz=Huge(isuppz(1))
        
        if (EnsDim<20) then 
            lwork=26*EnsDim
            liwork=lwork
        else
            lwork=(EnsDim+6)*EnsDim
        end if
        liwork=lwork
        
        allocate(iwork(liwork))
        iwork=Huge(iwork(1))
        
        allocate(work(lwork))
        work=Huge(work(1))
        
        allocate(eigenvalues(EnsDim))
        eigenvalues=Huge(eigenvalues(1))
        
        allocate(eigenvectors(EnsDim, EnsDim))
        eigenvectors=Huge(eigenvectors(1,1))
        
    
    end subroutine
    
    subroutine Finalize_Algebra
        deallocate(isuppz)
        deallocate(iwork)
        deallocate(work)
        deallocate(eigenvalues)
        deallocate(eigenvectors)
        
    end subroutine

    subroutine SymChangeBase(matrix, ierr, debug_opt)
    
        double precision, dimension(EnsDim, EnsDim), intent(inout) :: matrix
        integer, intent(in), optional :: debug_opt
        integer, intent(out) :: ierr
        integer :: indexi, neigenvalues, debug
        double precision dlamch
        
        debug=0
        if (present(debug_opt)) debug=debug_opt
        ierr=0
        
        if (debug>0) write(*,*) 'before dsyevr'
        
        call dsyevr("V", "A", "U", EnsDim, matrix, EnsDim, 0.0d0, 0.0d0,0.0d0, 0.0d0, &
            dlamch('S'), neigenvalues, eigenvalues, eigenvectors, EnsDim, &
            isuppz, work, lwork, iwork, liwork, ierr)
            
        if (debug>0) write(*,*) 'after dsyevr'

        if (ierr/=0) then
            write(*,*) "something wrong with svd. ierr=", ierr
            return
        end if

        if (EnsDim/=neigenvalues) then
            ierr=1
            write(*,*) "something strange in the number of eigenvalues! neigenvalues=", neigenvalues
            return
        end if
        
        if (debug>0) write(*,*) eigenvalues
        
        if (eigenvalues(1)<10**(-10)) then
            ierr=1
            write(*,*) "Small eigenvalue:", eigenvalues
            return
        end if
        
        do indexi=1, EnsDim
            matrix(indexi,:)=eigenvectors(:,indexi)/sqrt(eigenvalues(indexi))
        end do
        matrix=MatMul(eigenvectors,matrix)

    end subroutine


    subroutine OrtMatrix(CurrentMatrix, nRows,nCols, nInitialCols)
        
        integer, intent(in) :: nRows, nInitialCols, nCols
        double precision, dimension (nRows,nCols), intent(inout) :: CurrentMatrix
        integer :: indexi, indexj
        
        call NormalNumber(nRows*(nCols-nInitialCols),CurrentMatrix(:,nInitialCols+1:nCols))
        
        do indexi=nInitialCols+1, nCols
            CurrentMatrix(:,indexi)=CurrentMatrix(:,indexi)/norm2(CurrentMatrix(:,indexi))
            do indexj=1, indexi-1
                CurrentMatrix(:,indexi)=CurrentMatrix(:,indexi) & 
                    -dot_product(CurrentMatrix(:,indexi),CurrentMatrix(:,indexj))*CurrentMatrix(:,indexj)
                CurrentMatrix(:,indexi)=CurrentMatrix(:,indexi)/norm2(CurrentMatrix(:,indexi))
            end do
        end do
        
    end subroutine

    subroutine NormalNumber(nsize, OutputArray)
        
        integer, intent(in) :: nsize
        double precision, dimension(nsize), intent(out) :: OutputArray
        double precision, dimension(2) :: workvector
        integer :: counter
        double precision :: circlerange
        
        counter=0
        do while (counter<nsize)
            call random_number(workvector)
            workvector=workvector*2.0d0
            workvector=workvector-1.0d0
            circlerange=workvector(1)*workvector(1)+workvector(2)*workvector(2)
            if ((circlerange >= 1.0d0).or.(circlerange==0.0d0)) cycle
            circlerange=sqrt(-2.0d0*log(circlerange)/circlerange)
            workvector=workvector*circlerange
            counter=counter+1
            OutputArray(counter)=workvector(1)
            if (counter<nsize) then
                counter=counter+1
                OutputArray(counter)=workvector(2)
            end if
        end do
    end subroutine

    double precision function det(n, matrix)
        
        integer, intent(in) :: n
        double precision, dimension(n,n), intent(inout) :: matrix
        
        integer, dimension(:), allocatable :: ipiv
        integer info
        integer indexi
        
        allocate(ipiv(n))
        
        call dgetrf(n,n,matrix,n,ipiv,info)
        
        deallocate(ipiv)
        
        det=1.0d0
        do indexi=1, n
            det=det*matrix(indexi,indexi)
        end do
        
    end function

end module
