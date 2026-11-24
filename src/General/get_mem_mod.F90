module check_mem
  use iso_c_binding
  implicit none
  ! get_mem is implemented in C (get_mem.c); this is only its interface.
  ! An empty Fortran function here would define a second get_mem symbol in the
  ! library, conflicting with the C one
  interface
     function get_mem(err) bind(c)
         integer :: err
         double precision :: get_mem
     end function get_mem
  end interface
end module check_mem
