module create_grid
  implicit none
  contains
    subroutine create_single_column()
    use node_array, only: create_column
      implicit none
      call create_column(.false.,xi_in=0.d0)
    end subroutine

    subroutine create_cartesian_grid()
    use global_data, only: differentials
    use node_array, only: create_column
      implicit none
      real(8) :: xi
      xi=.6d0
      xi = xi - mod(xi,differentials(1))-1d-10
      do
        if(xi<0.0) exit
        call create_column(.false.,xi_in=xi)
        xi = xi - differentials(1)
      end do
      call create_column(.false.,xi_in=0.d0)
    end subroutine create_cartesian_grid

    subroutine create_duct_grid()
    use global_data, only: differentials, dimensions, PI
    use node_array, only: create_column, set_node
    use types, only: node_data
      implicit none
      real(8) :: xi, deta
      real(8), allocatable, dimension(:) :: eta
      integer :: n, i, j
      type (node_data) :: left, right
      xi=3.6d0
      xi = xi - mod(xi,differentials(1))-1d-10
      do
        if(xi<0.0) exit
        allocate(eta(dimensions(2)))
        if(xi<0.5)then
          eta(1) = 0.d0
        elseif(xi<1.0)then
          eta(1) = (xi-.5d0)*tan(PI*15.d0/180.d0)
        else
          eta(1) = (.5d0)*tan(PI*15.d0/180.d0)
        end if
        deta = (1.-eta(1))/dimensions(2)
        eta(1) = eta(1) + .5d0*deta
        do n = 2, dimensions(2)
          eta(n) = eta(n-1)+deta
        end do
        call create_column(.false.,xi_in=xi,eta_in=eta)
        deallocate(eta)
        xi = xi - differentials(1)
      end do
      call create_column(.false.,xi_in=0.d0)

      do i = 1, dimensions(1)
        do j = 1, dimensions(2)
          call set_node([i,j],geom_from_grid([i,j]))
        end do
      end do

    end subroutine create_duct_grid

    function geom_from_grid(indices)
    use types, only: node_data
    use node_array, only: get_node
    use global_data, only: differentials
      implicit none
      type (node_data) :: left, right, up, down, p
      type (node_data) :: geom_from_grid
      real(8) :: dxi, deta
      integer, dimension(:), intent(in) :: indices

      dxi = 2.d0*differentials(1)
      deta= 2.d0*differentials(2)
      p = get_node(indices)
      left = get_node(indices,'left' )
      right= get_node(indices,'right')
      up   = get_node(indices,'up'   )
      down = get_node(indices,'down' )
      if( left%boundary)then
        left = p
        dxi = differentials(1)
      end if
      if(right%boundary)then
        right= p
        dxi = differentials(1)
      end if
      if(up%boundary)then
        up   = p
        deta = differentials(2)
      end if
      if(down%boundary)then
        down = p
        deta = differentials(2)
      end if
      p%A = (right%x-left%x)/dxi
      p%B = (right%y-left%y)/dxi
      p%L = (   up%x-down%x)/deta
      p%M = (   up%y-down%y)/deta
      geom_from_grid = p
    end function geom_from_grid

end module create_grid
