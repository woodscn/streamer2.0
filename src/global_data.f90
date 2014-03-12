module global_data
  implicit none
  save
  type window_data
     real(8) :: xmin
     real(8) :: xmax
     real(8) :: ymin
     real(8) :: ymax
  end type window_data

  real(8), parameter :: gamma_const = 1.4d0
  real(8), parameter :: PI = 3.141592653589793d0
  type (window_data), protected :: window
  real(8), dimension(2), protected :: differentials = (/ 0., 0. /) ! (/ dxi, deta /)
  integer, dimension(2), protected :: dimensions = (/ 0, 0 /) !(/ nxi, neta /)
  real(8), protected :: h0
  logical, protected :: grid_preserving
  real(8), protected :: CFL
  real(8), protected :: time_step
  character(len=16), protected :: splitting_type
  character(len=128),protected :: input_file_name
  character(len=128),protected :: output_file_name
  integer, protected :: skip
  real(8), protected :: tmax
  integer, protected :: update_type
  logical, protected :: muscl_flag

  contains
    subroutine set_muscl_flag(in)
      implicit none
      logical, intent(in) :: in
      muscl_flag = in
    end subroutine set_muscl_flag

    subroutine set_grid_preserving(in)
      implicit none
      logical, intent(in) :: in
      grid_preserving = in
    end subroutine

    subroutine set_output_file_name(in)
      implicit none
      character(len=128), intent(in) :: in
      output_file_name = in
    end subroutine

    subroutine set_update_type(in)
      implicit none
      integer :: in
      update_type = in
    end subroutine

    subroutine set_h0(in)
      implicit none
      real(8) :: in
      h0 = in
    end subroutine

    subroutine set_CFL(in)
      implicit none
      real(8) :: in
      CFL = in
    end subroutine

    subroutine set_skip(in)
      implicit none
      integer :: in
      skip = in
    end subroutine

    subroutine set_tmax(in)
      implicit none
      real(8) :: in
      tmax = in
    end subroutine

    subroutine set_input_file_name(in)
      implicit none
      character(len=*) :: in
      input_file_name = in
    end subroutine

    subroutine set_splitting_type(in)
      implicit none
      character(len=*) :: in
      splitting_type = in
    end subroutine

    subroutine set_time_step(in)
      implicit none
      real(8), intent(in) :: in
      time_step = in
    end subroutine

    subroutine set_window(in)
      implicit none
      real(8), intent(in), dimension(:) :: in
      window%xmin = in(1) ; window%xmax = in(2)
      window%ymin = in(3) ; window%ymax = in(4)
    end subroutine

    subroutine set_differentials(in)
      implicit none
      real(8), intent(in), dimension(:) :: in
      differentials = in
    end subroutine

    subroutine set_dimensions(in)
      implicit none
      integer, intent(in), dimension(:) :: in
      dimensions = in
    end subroutine

end module global_data
