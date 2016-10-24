module check_mem
  use iso_c_binding
  implicit none
contains
  function get_mem(err) bind(c)
      integer :: err
      double precision :: get_mem
  end function get_mem
end module check_mem

