      INTEGER, parameter :: jpk_ben = 100

      INTEGER, parameter :: jptra_b = 1

      INTEGER, parameter :: jptra_var_b    = 1

      INTEGER, parameter :: jptra_flux_b   = 1

      INTEGER, parameter :: jptra_dia_b    = jptra_var_b + jptra_flux_b
       
      INTEGER, parameter :: jptra_dia_b_2d = 1

 
C State variables

        integer,parameter :: ppQ11c=1 

C       diagnostic indexes
        integer,parameter :: ppDIA_b=1

c       flux indexes
        integer, parameter:: ppFLUX_b=1

c       variables 2d
        integer, parameter:: ppDIA2d_b=1
