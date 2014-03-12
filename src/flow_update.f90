module flow_update
use types, only: node_data
  implicit none
  save
    real(8) :: time, diff
    character(len=5) :: more, less, direction
    integer :: i, j
    type(node_data) :: old_node_value
!    type(node_data), dimension(400,200) :: old_node_array
    type(node_data), allocatable, dimension(:,:) :: old_node_array, interfaces, current_array
contains

  subroutine flow_driver(scheme)
  use global_data, only: dimensions, update_type, grid_preserving
  use node_array, only: set_node
  use h_update!, only: iterative_update
    implicit none
    character(len=*), intent(in), optional :: scheme
    allocate(old_node_array(dimensions(1),dimensions(2)), current_array(dimensions(1),dimensions(2)))
    select case(update_type)
    case(1)
      call stepping
    case(2)
      call splitting()
      if(grid_preserving .and. dimensions(1) > 1) call iteration(.false.,.false.)
      call move_cell_array
    case default
      write(*,*) 'Error in flow_driver -- invalid update_type'
      stop
    end select


    deallocate(old_node_array, current_array)

  end subroutine flow_driver

  subroutine stepping()
  use global_data, only: dimensions, update_type
    implicit none
    do i = 1, dimensions(1)
      do j = 1, dimensions(2)
        select case(update_type)
        case(1)
          call splitting()
          call move_cell()
        case(2)
          call update_in_place(verbose_flag = .false.)
        case default
          write(*,*) 'Error in stepping -- invalid update_type'
          stop
        end select
      end do
    end do
  end subroutine stepping

  subroutine splitting(verbose_flag)
  use global_data, only: time_step, differentials, splitting_type, update_type, dimensions
  use node_array, only: set_node
    implicit none
    logical, intent(in), optional :: verbose_flag
    logical :: verbose = .false.
    integer :: n, nmax = 0
    real(8) :: txi = 0.d0, teta = 0.d0

    if (present(verbose_flag)) verbose = verbose_flag
    select case(splitting_type)
    case('godunov')
      nmax = 2
      txi  = time_step
      teta = time_step
    case('strang')
      nmax = 3
      txi  = .5d0*time_step
      teta = time_step
    case default
      write(*,*) 'Error in splitting: Bad type flag'
      stop
    end select

    do n = 1, nmax
       if(mod(n,2)/=0)then
          direction = 'xi'
          diff = differentials(1)
          more = 'right'
          less = 'left'
          time = txi
          allocate(interfaces(dimensions(1)+1,dimensions(2)))
       else
          direction = 'eta'
          diff = differentials(2)
          more = 'up'
          less = 'down'
          time = teta
          allocate(interfaces(dimensions(1),dimensions(2)+1))
       end if
       select case(update_type)
       case(1)
         call update_in_place(verbose_flag = .false.)
       case(2)
         call stepping
         do i = 1, dimensions(1)
           do j = 1, dimensions(2)
             call set_node([i,j],current_array(i,j))
           end do
         end do
       case default
         write(*,*) 'Error in splitting -- Invalid update_type'
         stop
       end select

       deallocate(interfaces)

    end do
  end subroutine splitting

  subroutine update_in_place(verbose_flag)
  use types, only: node_data
  use node_array, only: get_node, set_node
  use global_data, only: update_type
    implicit none
    logical, intent(in), optional :: verbose_flag
    logical :: verbose = .false.
    type (node_data) :: interface_left, interface_right, current

    if(present(verbose_flag)) verbose = verbose_flag

    interface_left = interface_state('left')
    interface_right = interface_state('right')
    current = update_node(interface_left,interface_right)
    if(verbose) write(*,*) current
    select case(update_type)
    case(1)
      call set_node([i,j],current)
    case(2)
      current_array(i,j) = current
    case default
      write(*,*) 'Error in update_in_place -- Invalid update_type'
      stop
    end select

  end subroutine update_in_place

!  subroutine update_en_masse(verbose_flag)
!  use types, only: node_data
!  use node_array, only: set_node, set_node
!    implicit none
!    logical, intent(in), optional :: verbose_flag
!    logical :: verbose = .false.
!
!    do i = 1, dimensions(1)
!      do j = 1, dimensions(2)
!        interfaces(i,j) = interface_state('left')
!        current(i,j) = get_node([i,j])
!      end do
!    end do
!
!    select case(direction)
!    case('xi')
!      i = dimensions(1)
!      do j = 1, dimensions(2)
!        interfaces(i+1,j) = interface_state('right')
!      end do
!    case('eta')
!      j = dimensions(2)
!      do i = 1, dimensions(1)
!        interfaces(i,j+1) = interface_state('right')
!      end do
!    end select



  function interface_state(leftright)
  use node_array, only: get_node
  use muscl_update, only: muscl
  use riemann, only: riemann_solve, normal_vector
  use types, only: node_data
  use global_data, only: dimensions, muscl_flag
    implicit none
    character(len=*), intent(in) :: leftright
    type (node_data) :: left, right, interface_state
    real(8), dimension(2) :: normal

    select case(leftright)
    case('left')
       left = get_node([i,j],less)
       right= get_node([i,j])
       if( (less == 'left' .and. i > 2).or.(less == 'down' .and. j>2) )&
!       if(.not.(left%boundary) .and. muscl_flag)&
            call muscl(left,right,get_node([i,j],less,less),get_node([i,j],more))
    case('right')
       left = get_node([i,j])
       right= get_node([i,j],more)
       if( (more == 'right' .and. i < (dimensions(1)-1)).or.(more == 'up' .and. j<(dimensions(2)-1)))&
!       if(.not.right%boundary .and. muscl_flag)&
        call muscl(left,right,get_node([i,j],less),get_node([i,j],more,more))
!            call muscl(right,get_node([i,j],more),left,get_node([i,j],more,more))
    case default
      write(*,*) 'Error in interface_state -- Invalid leftright flag'
      stop
    end select

    normal = normal_vector(left,right,direction)
    interface_state = riemann_solve(left%rotate(normal,inverse_flag = .false.),right%rotate(normal,inverse_flag = .false.))
    interface_state = interface_state%rotate(normal,inverse_flag = .true.)
  end function interface_state




  function update_node(interface_left,interface_right)
  use types, only: node_data
  use node_array, only: get_node
    implicit none
    type (node_data), intent(in) :: interface_left, interface_right
    type (node_data) :: current, update_node
    type (node_data) :: flux_left, flux_right, flux_diff, change

    current   = get_node([i,j])
    old_node_value = current
    old_node_array(i,j) = old_node_value
    current   = current%primtocons()

    flux_left = geom_flux(current,interface_left)
    flux_right= geom_flux(current,interface_right)
    flux_diff = flux_right%subtract_node_data(flux_left)
    change    = flux_diff%scalar_node_data_multiply(time/diff)
    current   = current%subtract_node_data(change)

    flux_left = phys_flux(current,interface_left)
    flux_right= phys_flux(current,interface_right)
    flux_diff = flux_right%subtract_node_data(flux_left)
    change    = flux_diff%scalar_node_data_multiply(time/diff)
    current   = current%subtract_node_data(change)
    
    current   = current%constoprim()
    current%q = sqrt(current%u**2+current%v**2)
    current%theta = atan2(current%v,current%u)

!    call artificial_viscosity(current,2)

    update_node = current

  end function

  function total_flux(state,interface)
    use types, only: node_data
    implicit none
    type(node_data),intent(in) :: state, interface
    type(node_data) :: total_flux

    total_flux = geom_flux(state,interface)
    total_flux = total_flux%add_node_data(phys_flux(state,interface))

  end function total_flux

  function geom_flux(state,interface)
    use types, only: node_data
    implicit none
    type (node_data), intent(in) :: state, interface
    type (node_data) :: geom_flux

    select case(direction)
    case('xi')
       geom_flux%A = -state%h*interface%u
       geom_flux%B = -state%h*interface%v
    case('eta')
       geom_flux%L = -state%h*interface%u
       geom_flux%M = -state%h*interface%v
    case default
       write(*,*) 'Error in flux -- bad direction tag'
       stop
    end select

  end function geom_flux

  function phys_flux(state,interface)
    use types, only: node_data
    use global_data, only: gamma_const
    implicit none
    type (node_data), intent(in) :: state , interface
    type (node_data) :: phys_flux, out

    real(8), dimension(2) :: normal
    real(8) :: v_normal, e

    select case(direction)
    case('xi')
       normal = [state%M, -state%L]
!       out%A = -state%h*interface%u
!       out%B = -state%h*interface%v
    case('eta')
       normal = [-state%B, state%A]
!       out%L = -state%h*interface%u
!       out%M = -state%h*interface%v
    case default
       write(*,*) 'Error in flux -- bad direction tag'
       stop
    end select

    v_normal = dot_product([interface%u,interface%v],normal)
    e = .5d0*(interface%u**2+interface%v**2)+interface%p/((gamma_const-1.d0)*interface%rho)

    out%rho = interface%rho*(1.d0-state%h)*v_normal
    out%p   = interface%rho*(1.d0-state%h)*v_normal*e           + interface%p*v_normal
    out%u   = interface%rho*(1.d0-state%h)*v_normal*interface%u + interface%p*normal(1)
    out%v   = interface%rho*(1.d0-state%h)*v_normal*interface%v + interface%p*normal(2)
    phys_flux = out

  end function phys_flux

  subroutine metric_update(state, change)
  use types, only: node_data
    implicit none
    type (node_data), intent(inout) :: state
    type (node_data), intent(inout) :: change

    state%A = state%A - change%A
    state%B = state%B - change%B
    state%L = state%L - change%L
    state%M = state%M - change%M

    change%A = 0.d0 ; change%B = 0.d0 ; change%L = 0.d0 ; change%M = 0.d0

  end subroutine metric_update

  subroutine move_cell()
  use types, only: node_data
  use global_data, only: time_step
  use node_array, only: get_node, set_node
    implicit none
    type (node_data) :: new

    new = get_node([i,j])
    new%x = old_node_value%x + time_step*0.5d0*(old_node_value%h*old_node_value%u + new%h*new%u)
    new%y = old_node_value%y + time_step*0.5d0*(old_node_value%h*old_node_value%v + new%h*new%v)
    call set_node([i,j],new)
  end subroutine move_cell

! This uses incorrect values for old_node_array. It uses the values after the
! eta time step, rather than at the beginning of the whole time step.
  subroutine move_cell_array()
  use types, only: node_data
  use global_data, only: time_step, dimensions
  use node_array, only: get_node, set_node
    implicit none
    type (node_data) :: new, old

    do i = 1, dimensions(1)
      do j = 1, dimensions(2)
        new = get_node([i,j])
        old = old_node_array(i,j)
        new%x = old%x + time_step*0.5d0*(old%h*old%u + new%h*new%u)
        new%y = old%y + time_step*0.5d0*(old%h*old%v + new%h*new%v)
        call set_node([i,j],new)
      end do
    end do
  end subroutine move_cell_array

  subroutine artificial_viscosity(p,order)
    use node_array, only: get_node
    implicit none
    type (node_data) :: p,pm,pmm,pp,ppp,temp
    integer :: order

    if(order==4)then
       pm = get_node([i,j],less)
       if(.not. pm%boundary)then
          pmm= get_node([i,j],less,less)
       else
          pmm= pm
       end if
       pp = get_node([i,j],more)
       if(.not. pp%boundary)then
          ppp= get_node([i,j],more,more)
       else
          ppp=pp
       end if
       temp = ppp%subtract_node_data(pp%scalar_node_data_multiply(4.d0))
       temp = temp%add_node_data(p%scalar_node_data_multiply(6.d0))
       temp = temp%subtract_node_data(pm%scalar_node_data_multiply(4.d0))
       temp = temp%add_node_data(ppp)
       temp = temp%scalar_node_data_multiply(1.d0/diff**4)
    elseif(order==2)then
       pm = get_node([i,j],less)
       pp = get_node([i,j],more)
       temp = pm%subtract_node_data(p%scalar_node_data_multiply(2.d0))
       temp = temp%add_node_data(pp)
       temp = temp%scalar_node_data_multiply(1.d0/diff**2)
    end if
    
    temp = temp%scalar_node_data_multiply(.0001d0)
!write(*,*) temp%p, temp%rho, temp%u, temp%v
!read(*,*)

!    p = p%add_node_data(temp%scalar_node_data_multiply(-1d-1))
    p%rho = p%rho + temp%rho
    p%p = p%p + temp%p
    p%u = p%u + temp%u
    p%v = p%v + temp%v
    
    
  end subroutine artificial_viscosity
    

end module flow_update
