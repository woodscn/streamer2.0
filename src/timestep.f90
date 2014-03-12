module timestep
  implicit none
contains
  function spec_rad(p,direction)
    use types, only: node_data
    use global_data, only: gamma_const
    character(len=*), intent(in) :: direction
    type (node_data), intent(in) :: p
    real(8) spec_rad
    real(8), dimension(2) :: normal

    spec_rad = 0.d0
    select case(direction)
    case('xi')
       normal = [p%M,-p%L]
    case('eta')
       normal = [-p%B,p%A]
    case default
       write(*,*) 'Error in spec_rad -- Incompatible direction specified'
       stop
    end select
    spec_rad = ((1.0 - p%h)*abs( normal(1)*p%u + p%v*normal(2) ) + &
            sqrt(gamma_const*p%p/p%rho)*sqrt((normal(1)*p%u)**2+(normal(2)*p%v)**2))&
            /(p%A*p%M - p%B*p%L)
  end function spec_rad

  function time_step_func()
    use global_data, only: dimensions, CFL, differentials, window
    use types, only: node_data
    use node_array, only: get_node
    implicit none
    integer :: i, j, n
    integer, save :: counter = 0.d0
    type (node_data) :: p
    real(8) :: time_step_func
    real(8) :: rho_xi = 0.0 , rho_eta = 0.0 , dt = 1.0 , xnew, x, temp

    rho_xi = 0.d0
    rho_eta= 0.d0
    do i = 1 , dimensions(1)
       do j = 1 , dimensions(2)
          p = get_node([i,j])
          rho_xi = max(rho_xi ,spec_rad(p,'xi' ))
          rho_eta= max(rho_eta,spec_rad(p,'eta'))
       end do
    end do

    time_step_func = CFL*2.0/(rho_xi/differentials(1) + rho_eta/differentials(2))
!!$    time_step_func = 1d-5
    temp = 0.
    x = 0.
    xnew = 0.

    do j = 1 , dimensions(2)
      p = get_node([1,j])

      temp = p%x + p%h*p%u*time_step_func
      if(temp > xnew)then
        xnew = temp
        n = j
      end if
    end do
    if( xnew > window%xmin + differentials(1) )then
      p = get_node([1,n])
!      time_step_func = max((window%xmin + differentials(1) - p%x)/(p%h*p%u) , 1.d-14)
!write(*,*) window%xmin, differentials(1)
      time_step_func = (window%xmin + differentials(1) - p%x)/(p%h*p%u)
    end if
!    p = get_node([dimensions(1)-counter,1])
!    if( p%x + p%h*p%u*time_step_func > 0.5 )then
!      time_step_func = (0.5-p%x)/(p%h*p%u)
!      counter = counter + 1
!    end if


  end function time_step_func

end module timestep
