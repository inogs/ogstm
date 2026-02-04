!===========================================================
! read_config_yaml_fortran_yaml.f90
!
! Single-file example for Fortran-YAML (Bolding & Bruggeman)
!===========================================================

module ogstm_yaml_reader
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use yaml
  use yaml_types
  implicit none
  private

  public :: ConfigYAML
  public :: read_config_yaml, print_config

  type :: InteriorStateItem
     character(len=:), allocatable :: name
     real(wp) :: ctrmax = 0.0_wp
     integer  :: ctrhf  = 0
     integer  :: relax  = 0
  end type

  type :: InteriorDiagnosticItem
     character(len=:), allocatable :: name
     integer :: diahf = 0
     integer :: diaWR = 0
  end type

  type :: HorizontalDiagnosticItem
     character(len=:), allocatable :: name
     integer :: diahf_2d = 0
     integer :: diaWR_2d = 0
  end type

  type :: ConfigYAML
     type(InteriorStateItem),        allocatable :: interior_state(:)
     type(InteriorDiagnosticItem),   allocatable :: interior_diagnostic(:)
     type(HorizontalDiagnosticItem), allocatable :: horizontal_diagnostic(:)
  end type

contains

  !---------------------------------------------------------
  ! Read YAML file -> ConfigYAML
  !---------------------------------------------------------
  subroutine read_config_yaml(filename, cfg)
    character(len=*), intent(in)  :: filename
    type(ConfigYAML), intent(out) :: cfg

    character(error_length) :: err
    class(type_node), pointer :: root
    class(type_dictionary), pointer :: root_dict
    type(type_error), pointer :: e

    err = ''
    nullify(e)

    root => parse(filename, unit=100, error=err)
    if (err /= '') then
       write(*,*) 'PARSE ERROR: ', trim(err)
       stop 1
    end if

    ! Root must be a dictionary for your YAML layout
    select type (root)
    type is (type_dictionary)
       root_dict => root
    class default
       write(*,*) 'ERROR: YAML root is not a dictionary.'
       stop 2
    end select

    call read_interior_state(root_dict, cfg%interior_state)
    call read_interior_diagnostic(root_dict, cfg%interior_diagnostic)
    call read_horizontal_diagnostic(root_dict, cfg%horizontal_diagnostic)

    ! (optional) free tree if your build uses finalize patterns elsewhere
    ! call root_dict%finalize()   ! only if you manage lifetime explicitly
  end subroutine read_config_yaml

  !---------------------------------------------------------
  ! Helpers: check if a key exists in a dictionary
  ! (Fortran-YAML stores pairs in a linked list: dict%first%next%...)
  !---------------------------------------------------------
  logical function has_key(dict, key)
    class(type_dictionary), intent(in) :: dict
    character(len=*), intent(in) :: key
    type(type_key_value_pair), pointer :: p
    has_key = .false.
    p => dict%first
    do while (associated(p))
       if (trim(p%key) == trim(key)) then
          has_key = .true.
          return
       end if
       p => p%next
    end do
  end function has_key

  integer function dict_size(dict)
    class(type_dictionary), intent(in) :: dict
    type(type_key_value_pair), pointer :: p
    dict_size = 0
    p => dict%first
    do while (associated(p))
       dict_size = dict_size + 1
       p => p%next
    end do
  end function dict_size

  !---------------------------------------------------------
  ! Read section: interior_state
  !---------------------------------------------------------
  subroutine read_interior_state(root, arr)
    class(type_dictionary), intent(in) :: root
    type(InteriorStateItem), allocatable, intent(out) :: arr(:)

    class(type_dictionary), pointer :: sec
    class(type_dictionary), pointer :: item
    type(type_error), pointer :: e
    type(type_key_value_pair), pointer :: p
    integer :: n, i

    nullify(e)

    if (.not. has_key(root, "interior_state")) then
       allocate(arr(0))
       return
    end if

    sec => root%get_dictionary("interior_state", required=.true., error=e)
    if (associated(e)) then
       write(*,*) "ERROR reading interior_state"
       stop 3
    end if

    n = dict_size(sec)
    allocate(arr(n))

    p => sec%first
    i = 0
    do while (associated(p))
       i = i + 1
       arr(i)%name = trim(p%key)

       select type (val => p%value)
       type is (type_dictionary)
          item => val
       class default
          write(*,*) "ERROR: interior_state.", trim(p%key), " is not a mapping"
          stop 4
       end select

       ! get_real/get_integer support defaults (optional)
       arr(i)%ctrmax = real(item%get_real("ctrmax", default=0.0_real_kind, error=e), wp)
       arr(i)%ctrhf  = item%get_integer("ctrhf", default=0, error=e)
       arr(i)%relax  = item%get_integer("relax", default=0, error=e)

       if (associated(e)) then
          write(*,*) "ERROR converting interior_state fields for ", trim(arr(i)%name)
          stop 5
       end if

       p => p%next
    end do
  end subroutine read_interior_state

  !---------------------------------------------------------
  ! Read section: interior_diagnostic
  !---------------------------------------------------------
  subroutine read_interior_diagnostic(root, arr)
    class(type_dictionary), intent(in) :: root
    type(InteriorDiagnosticItem), allocatable, intent(out) :: arr(:)

    class(type_dictionary), pointer :: sec
    class(type_dictionary), pointer :: item
    type(type_error), pointer :: e
    type(type_key_value_pair), pointer :: p
    integer :: n, i

    nullify(e)

    if (.not. has_key(root, "interior_diagnostic")) then
       allocate(arr(0))
       return
    end if

    sec => root%get_dictionary("interior_diagnostic", required=.true., error=e)
    if (associated(e)) then
       write(*,*) "ERROR reading interior_diagnostic"
       stop 6
    end if

    n = dict_size(sec)
    allocate(arr(n))

    p => sec%first
    i = 0
    do while (associated(p))
       i = i + 1
       arr(i)%name = trim(p%key)

       select type (val => p%value)
       type is (type_dictionary)
          item => val
       class default
          write(*,*) "ERROR: interior_diagnostic.", trim(p%key), " is not a mapping"
          stop 7
       end select

       arr(i)%diahf = item%get_integer("diahf", default=0, error=e)
       arr(i)%diaWR = item%get_integer("diaWR", default=0, error=e)

       if (associated(e)) then
          write(*,*) "ERROR converting interior_diagnostic fields for ", trim(arr(i)%name)
          stop 8
       end if

       p => p%next
    end do
  end subroutine read_interior_diagnostic

  !---------------------------------------------------------
  ! Read section: horizontal_diagnostic
  !---------------------------------------------------------
  subroutine read_horizontal_diagnostic(root, arr)
    class(type_dictionary), intent(in) :: root
    type(HorizontalDiagnosticItem), allocatable, intent(out) :: arr(:)

    class(type_dictionary), pointer :: sec
    class(type_dictionary), pointer :: item
    type(type_error), pointer :: e
    type(type_key_value_pair), pointer :: p
    integer :: n, i

    nullify(e)

    if (.not. has_key(root, "horizontal_diagnostic")) then
       allocate(arr(0))
       return
    end if

    sec => root%get_dictionary("horizontal_diagnostic", required=.true., error=e)
    if (associated(e)) then
       write(*,*) "ERROR reading horizontal_diagnostic"
       stop 9
    end if

    n = dict_size(sec)
    allocate(arr(n))

    p => sec%first
    i = 0
    do while (associated(p))
       i = i + 1
       arr(i)%name = trim(p%key)

       select type (val => p%value)
       type is (type_dictionary)
          item => val
       class default
          write(*,*) "ERROR: horizontal_diagnostic.", trim(p%key), " is not a mapping"
          stop 10
       end select

       arr(i)%diahf_2d = item%get_integer("diahf_2d", default=0, error=e)
       arr(i)%diaWR_2d = item%get_integer("diaWR_2d", default=0, error=e)

       if (associated(e)) then
          write(*,*) "ERROR converting horizontal_diagnostic fields for ", trim(arr(i)%name)
          stop 11
       end if

       p => p%next
    end do
  end subroutine read_horizontal_diagnostic

  !---------------------------------------------------------
  ! Debug print
  !---------------------------------------------------------
  subroutine print_config(cfg)
    type(ConfigYAML), intent(in) :: cfg
    integer :: i

    write(*,'(a)') "=== interior_state ==="
    do i = 1, size(cfg%interior_state)
       write(*,'(a)') " - "//trim(cfg%interior_state(i)%name)
       write(*,'(a,1x,es12.5)') "    ctrmax:", cfg%interior_state(i)%ctrmax
       write(*,'(a,1x,i0)')     "    ctrhf: ", cfg%interior_state(i)%ctrhf
       write(*,'(a,1x,i0)')     "    relax: ", cfg%interior_state(i)%relax
    end do

    write(*,'(a)') "=== interior_diagnostic ==="
    do i = 1, size(cfg%interior_diagnostic)
       write(*,'(a)') " - "//trim(cfg%interior_diagnostic(i)%name)
       write(*,'(a,1x,i0)') "    diahf:", cfg%interior_diagnostic(i)%diahf
       write(*,'(a,1x,i0)') "    diaWR:", cfg%interior_diagnostic(i)%diaWR
    end do

    write(*,'(a)') "=== horizontal_diagnostic ==="
    do i = 1, size(cfg%horizontal_diagnostic)
       write(*,'(a)') " - "//trim(cfg%horizontal_diagnostic(i)%name)
       write(*,'(a,1x,i0)') "    diahf_2d:", cfg%horizontal_diagnostic(i)%diahf_2d
       write(*,'(a,1x,i0)') "    diaWR_2d:", cfg%horizontal_diagnostic(i)%diaWR_2d
    end do
  end subroutine print_config

end module ogstm_yaml_reader


!program test_read_yaml
!  use config_yaml_reader
!  implicit none
!  type(ConfigYAML) :: cfg

!  call read_config_yaml("config.yaml", cfg)
!  call print_config(cfg)
!end program test_read_yaml

