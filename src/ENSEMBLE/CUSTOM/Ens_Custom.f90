! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Custom
    use modul_param, &
        only: jpi, jpj, jpk, jptra
    use myalloc, &
        only: nldj, nlej, nldi, nlei, bfmmask, tmask
    use Ens_Mem, &
        only: EnsSize, &
            LocalRange
        
    implicit none
    
    double precision, parameter :: small=exp(-20.0d0), logsmall=-20.0d0
    integer :: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate
    double precision, pointer, contiguous, dimension(:) :: DAstate 
    integer :: n_DAstate 
    integer :: win_DAstate
    double precision, pointer, contiguous, dimension(:,:) :: gl_DAstate
    double precision, pointer, contiguous, dimension(:,:,:,:) :: DAstate_kjit
    integer, dimension(:), allocatable :: DAVariablesIndex
    integer(1), dimension(:,:,:), allocatable :: DAMask
    
contains
    
    subroutine Ens_Init_DA
        use myalloc, &
            only: ctrcnm!, bfmmask, glamt
        use Ens_Mem, &
            only: Ens_shared_alloc
        
        !double precision, dimension(:), POINTER, contiguous :: member_pointer
        double precision, dimension(:,:), POINTER, contiguous :: global_pointer
        integer indexi, indexj
        
        nk_DAstate=jpk
        nj_DAstate=nlej-nldj+1
        ni_DAstate=nlei-nldi+1
        
        ntra_DAstate=0
        do indexi=1,jptra
            if (IsEnsDAVar(ctrcnm(indexi))) ntra_DAstate=ntra_DAstate+1
        end do
        allocate(DAVariablesIndex(ntra_DAstate))
        ntra_DAstate=0
        do indexi=1,jptra
            if (IsEnsDAVar(ctrcnm(indexi))) then
                ntra_DAstate=ntra_DAstate+1
                DAVariablesIndex(ntra_DAstate)=indexi
            end if
        end do
        
        n_DAstate=nk_DAstate * nj_DAstate * ni_DAstate * ntra_DAstate
        
        call Ens_shared_alloc(n_DAstate, DAstate, global_pointer, win_DAstate)
        !DAstate(1:jpk,1:jpj,1:jpi,1:jptra)=>member_pointer
        gl_DAstate(1:n_DAstate, 0:EnsSize-1)=>global_pointer
        DAstate_kjit(1:nk_DAstate, 1:nj_DAstate, 1:ni_DAstate, 1:ntra_DAstate)=>DAstate
        !DAstate  = huge(DAstate(1)) 
        DAstate  = 0.0d0
        
        allocate(DAMask(nk_DAstate, nj_DAstate, ni_DAstate))
        DAMask=bfmmask(:,nldj:nlej,nldi:nlei)
        !do indexi=1,ni_DAstate
        !    do indexj=1,nj_DAstate
        !        if (glamt(nldj-1+indexj,nldi-1+indexi)<-6.0d0) DAMask(:, indexj,indexi)=0
        !    end do
        !end do
        
    end
    
    subroutine Ens_Finalize_DA
        use mpi 
        
        integer ierror
        
        deallocate(DAVariablesIndex)
        CALL MPI_Win_free(win_DAstate, ierror)
        
        deallocate(DAMask)
        
    end subroutine
    
    subroutine Ens_state2DA
        use myalloc, &
            only: trn
            
        ! Transform and prepare the array for DA filter.
        ! e.g., DAstate_kjit=log(trn)
        
        integer indexi, indexj, indexk, indext
        
        !DAstate_kjit=log(trn(:,nldj:nlej,nldi:nlei,DAVariablesIndex))
        do indext=1, ntra_DAstate
            do indexi=1, ni_DAstate
                do indexj=1, nj_DAstate
                    do indexk=1, nk_DAstate
                        if (DAMask(indexk, indexj, indexi)==0) then
                            !DAstate_kjit(indexk:nk_DAstate, indexj, indexi, indext)=0.0d0
                            exit
                        end if
                        if (trn(indexk, nldj-1+indexj, nldi-1+indexi, DAVariablesIndex(indext))<small) then
                            DAstate_kjit(indexk, indexj, indexi, indext)=logsmall
                        else
                            DAstate_kjit(indexk, indexj, indexi, indext)=log(trn(indexk, nldj-1+indexj, nldi-1+indexi, DAVariablesIndex(indext)))
                        end if
                    end do
                end do
            end do
        end do        
            
    end subroutine
    
    subroutine Ens_DA2state
        use myalloc, &
            only: trn
        use ogstm_mpi_module, &
            only: mpplnk_my
        
        ! Transform back, from DA output to trn.
        ! e.g. trn=exp(DAstate_kjit)
        
        integer indexi, indexj, indexk, indext
        
        !trn(:,nldj:nlej,nldi:nlei,DAVariablesIndex)=exp(DAstate_kjit)
        do indext=1, ntra_DAstate
            do indexi=1, ni_DAstate
                do indexj=1, nj_DAstate
                    do indexk=1, nk_DAstate
                        if (DAMask(indexk, indexj, indexi)==0) exit
                        trn(indexk, nldj-1+indexj, nldi-1+indexi, DAVariablesIndex(indext))=exp(DAstate_kjit(indexk, indexj, indexi, indext))
                    end do
                end do
            end do
        end do        
        
        do indexi=1,ntra_DAstate
            CALL mpplnk_my(trn(:,:,:,DAVariablesIndex(indexi)))
        end do
        
    end subroutine
    
    subroutine Ens_RST2DA(tracer)
        use myalloc, &
            only: trn
        !Use Ens_IO, &
        !    only: trn_IO
        
        double precision, dimension(jpk, jpj, jpi, jptra), intent(out) :: tracer
            
        ! Transform for averaging, going from trn to tracer
        ! e.g., tracer=log(trn)
        
        integer indexi, indexj, indexk, indext
        
        !tracer=log(trn)
        do indext=1, jptra
            do indexi=1, jpi
                do indexj=1, jpj
                    do indexk=1, jpk
                        if (tmask(indexk, indexj, indexi)==0) exit
                        if (trn(indexk, indexj, indexi, indext)<small) then
                            tracer(indexk, indexj, indexi, indext)=logsmall
                        else
                            tracer(indexk, indexj, indexi, indext)=log(trn(indexk, indexj, indexi, indext))
                        end if
                    end do
                end do
            end do
        end do        
            
    end subroutine
    
    subroutine Ens_DA2RST(tracer)
        !Use Ens_IO, &
        !    only: trn_IO
        
        double precision, dimension(jpk, jpj, jpi, jptra), intent(inout) :: tracer
        
        ! Transform back, but in place
        ! e.g., tracer=exp(tracer)
        
        integer indexi, indexj, indexk, indext
        
        !tracer=exp(tracer)
        do indext=1, jptra
            do indexi=1, jpi
                do indexj=1, jpj
                    do indexk=1, jpk
                        if (tmask(indexk, indexj, indexi)==0) exit
                        tracer(indexk, indexj, indexi, indext)=exp(tracer(indexk, indexj, indexi, indext))
                    end do
                end do
            end do
        end do        
        
    end subroutine
    
    subroutine Ens_Ave2DA
        use modul_param, &
            only: jptra_dia, jptra_dia_2d
        use myalloc, &
            only: jptra_high, jptra_dia_high, jptra_dia2d_high, jptra_phys, jptra_phys_2d, &
                traIO, traIO_HIGH, tra_DIA_IO, tra_DIA_IO_HIGH, tra_DIA_2d_IO, tra_DIA_2d_IO_HIGH, &
                tra_PHYS_IO, tra_PHYS_IO_HIGH, tra_PHYS_2d_IO, tra_PHYS_2d_IO_HIGH
        use DIA_mem, &
            only: Fsize, diaflx     
            
        ! Transform for averaging, in place.
        ! e.g., traIO=log(traIO)
        
        integer indexi, indexj, indexk, indext, indexp
        double precision, dimension(:,:,:,:), pointer :: temparray
        
        !traIO=log(traIO)
        !traIO_HIGH=log(traIO_HIGH)
        do indext=1, jptra
            temparray => traIO
            do indexp=1,2
                do indexi=1, jpi
                    do indexj=1, jpj
                        do indexk=1, jpk
                            if (tmask(indexk, indexj, indexi)==0) exit
                            if (temparray(indexk, indexj, indexi, indext)<small) then
                                temparray(indexk, indexj, indexi, indext)=logsmall
                            else
                                temparray(indexk, indexj, indexi, indext)=log(temparray(indexk, indexj, indexi, indext))
                            end if
                        end do
                    end do
                end do
                if (indext>jptra_high) exit
                temparray => traIO_HIGH
            end do
        end do        
            
    end subroutine
    
    subroutine Ens_DA2Ave
        use modul_param, &
            only: jptra_dia, jptra_dia_2d
        use myalloc, &
            only: jptra_high, jptra_dia_high, jptra_dia2d_high, jptra_phys, jptra_phys_2d, &
                traIO, traIO_HIGH, tra_DIA_IO, tra_DIA_IO_HIGH, tra_DIA_2d_IO, tra_DIA_2d_IO_HIGH, &
                tra_PHYS_IO, tra_PHYS_IO_HIGH, tra_PHYS_2d_IO, tra_PHYS_2d_IO_HIGH
        use DIA_mem, &
            only: Fsize, diaflx      
            
        ! Transform back, in place
        ! e.g., traIO=exp(traIO)
        
        integer indexi, indexj, indexk, indext, indexp
        double precision, dimension(:,:,:,:), pointer :: temparray
        
        !traIO=exp(traIO)
        !traIO_HIGH=exp(traIO_HIGH)
        do indext=1, jptra
            temparray => traIO
            do indexp=1,2
                do indexi=1, jpi
                    do indexj=1, jpj
                        do indexk=1, jpk
                            if (tmask(indexk, indexj, indexi)==0) exit
                            temparray(indexk, indexj, indexi, indext)=exp(temparray(indexk, indexj, indexi, indext))
                        end do
                    end do
                end do
                if (indext>jptra_high) exit
                temparray => traIO_HIGH
            end do
        end do        
            
    end subroutine
    
    function IsEnsDAVar(name)
        implicit none
        
        character(LEN=*), intent(in) :: name
        logical :: IsEnsDAVar
        
        !if ((name(1:1).eq."P").or.(name(1:1).eq."N").or.(name(1:1).eq."O")) then
        if (.true.) then
            IsEnsDAVar=.true.
        else
            IsEnsDAVar=.false.
        end if
            
    end function
    
    function LocalWeightFunction(radius2)
        
        double precision, intent(in) :: radius2
        double precision :: LocalWeightFunction
        
        LocalWeightFunction=(radius2/((LocalRange+1)**2)-1.0d0)**2
        ! si puo' anche usare y= 1 + ( -3 + 2*x )*x^2 che e' di grado 3 anziche' 4

    end function 
    
end module
