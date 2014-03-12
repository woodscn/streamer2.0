module grid_motion
  implicit none

contains
  subroutine move_grid(dt)
    use global_data, only: dimensions, window_data, window, differentials
    use node_array, only: get_node, set_node, create_column, delete_column
    use types, only: node_data
    implicit none
    real(8), intent(in) :: dt
    integer :: i, j
    type (node_data) :: p
    logical :: shortfall = .false., overrun = .false., first = .true.

    shortfall = .false. ; overrun = .false.
    do i = 1, dimensions(1)
       do j = 1, dimensions(2)
          p = get_node([i,j])
! Grid motion is now handled in the flow_update loop.
! This routine still handles the logic for creating
! and removing columns
!          p%x = p%x + p%h*p%u*dt
!          p%y = p%y + p%h*p%v*dt
          if( i == 1 .and. p%x >= window%xmin + differentials(1) ) &
               shortfall = .true.
!          if( p%x > window%xmax + differentials(1) ) &
          if( p%x > window%xmax ) &
               overrun = .true.
          call set_node([i,j],p)
       end do
    end do
    if(shortfall)then
!      write(*,*) 'Shortfall. Creating new column.'
      call create_column()
      if(first)then
        call delete_column()
        first = .false.
      end if
    end if
    if(overrun  )then
!      write(*,*) 'Overrun. Deleting end column.'
      call delete_column()
    end if

  end subroutine move_grid

end module grid_motion
