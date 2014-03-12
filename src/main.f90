program main
  use boundary_conditions_init, only: read_boundary
  use node_array, only: create_column, initialize_column, get_node, delete_column
  use create_grid!, only: create_single_column
  use timestep, only: time_step_func
  use grid_motion, only: move_grid
  use flow_update, only: update_in_place
  use output_data, only: write_files, write_files_matlab
  use global_data, only: set_time_step, dimensions, differentials, set_input_file_name, &
    tmax, skip, set_h0, h0
  use types, only: node_data
  use node_array, only: get_node, set_node
  use flow_update, only: flow_driver
  implicit none
  real(8) :: t = 0.0, dt = 0.0, xi
  integer :: nt = 0, ierror, i,j
  character(len=128) :: in
  logical :: change_h = .true., flag=.true.
  type(node_data) :: temp

  call get_command_argument(1,in,status=ierror)
!  if(ierror == 0)
  call set_input_file_name(in)
!  else
!    write(*,*) 'Error getting command line arguments'
!    stop
!  end if


  call read_boundary(verbose_flag=.false.)
  if( h0 /= 0 )then
    call create_single_column
  else
!    call create_cartesian_grid
     call create_duct_grid
  end if

  call write_files_matlab(0.d0, first_flag = .true.)
  do
     nt = nt + 1
     dt = min(time_step_func(), tmax - t)
     call set_time_step(dt)

     call flow_driver()

     call move_grid(dt)
     t = t + dt
     if(mod(nt,skip)==0)then
       write(*,*) "nt = ", nt, "t = ", t
       call write_files_matlab(t)
     end if
     if(dimensions(1) > 4 .and. change_h.and.flag)then
       call delete_column
       change_h = .false.
     elseif(dimensions(1) > 4 .and. flag)then
       call delete_column
       flag = .false.
     end if
!     if(t >  15.0 .and. change_h)then
!       change_h = .false.
!       do i = 1, dimensions(1)
!         do j = 1, dimensions(2)
!           temp = get_node([i,j])
!           temp%h = .999
!           call set_node([i,j],temp)
!         end do
!       end do
!       call set_h0(.999d0)
!     end if
     if(t >= tmax) exit
  end do
  call write_files_matlab(tmax)
  write(*,*) 'Reached tmax, exiting. t = ', t
end program main
