! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Custom
    use modul_param, &
        only: jpi, jpj, jpk, jptra
    use Ens_Mem, &
        only: n_DAstate, gl_DAstate
        
implicit none
    
contains
    
    subroutine Ens_Init_Custom
        n_DAstate=jpk*jpj*jpi*jptra
    end
    
    subroutine Ens_Ave2DA
        use modul_param, &
            only: jptra_dia, jptra_dia_2d
        use myalloc, &
            only: jptra_high, jptra_dia_high, jptra_dia2d_high, jptra_phys, jptra_phys_2d, &
                traIO, traIO_HIGH, tra_DIA_IO, tra_DIA_IO_HIGH, tra_DIA_2d_IO, tra_DIA_2d_IO_HIGH, &
                tra_PHYS_IO, tra_PHYS_IO_HIGH, tra_PHYS_2d_IO, tra_PHYS_2d_IO_HIGH
        use DIA_mem, &
            only: Fsize, diaflx            
            
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
            
    end subroutine
    
    subroutine Ens_state2DA
        use myalloc, &
            only: trn
        Use Ens_Mem, &
            only: DAstate
            
        DAstate=reshape(trn, (/n_DAstate/))
            
    end subroutine
    
    subroutine Ens_DA2state_0
        
    end subroutine
    
    subroutine Ens_DA2state_all
        use myalloc, &
            only: trn
        Use Ens_Mem, &
            only: DAstate
            
        trn=reshape(DAstate, (/jpk, jpj, jpi, jptra/))
        
    end subroutine
    
end module
