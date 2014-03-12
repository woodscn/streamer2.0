!__****h* main/boundary_conditions
!__ NAME
!__ 	Boundary_Conditions -- module containing routines for applying simulation
!__ 	boundary conditions
!__ AUTHOR
!__ 	C. Nathan Woods
!__ DESCRIPTION
!__ 	Boundary_conditions has access to the primary data structures from
!__ 	boundary_conditions_init on a read-only basis. It uses the boundary
!__ 	conditions described there to create boundary nodes (ghost cells) that
!__ 	enforce the actual simulation boundary conditions. It also contains
!__ 	a routine for initializing a new column of cells at a supersonic inflow
!__ 	boundary.
!__****

module boundary_conditions
  use boundary_conditions_init, only: top_boundary, bottom_boundary, left_boundary, right_boundary, leftcoords

contains
  function boundary_node(pcurrent, directionaltag, verbose_flag)
    use types, only: boundary_master, node_data, node
    use global_data, only: PI
    implicit none
    type (node_data), intent(in), target :: pcurrent
    character(len=*), intent(in) :: directionaltag
    logical, intent(in), optional :: verbose_flag
    type (boundary_master), dimension(:), allocatable :: masterp
    type (node_data) :: boundary_node
    real(8), pointer :: coordp => null()
    integer :: n, polyn, inda
    logical :: verbose = .false.
    real(8) :: poly_angle

    if(present(verbose_flag))verbose = verbose_flag

    select case(directionaltag)
    case('up')
       allocate( masterp(size(top_boundary)) )
       masterp = top_boundary
       coordp  => pcurrent%x
    case('down')
       allocate( masterp(size(bottom_boundary)) )
       masterp = bottom_boundary
       coordp  => pcurrent%x
    case('left')
       allocate( masterp(size(left_boundary)) )
       masterp = left_boundary
       coordp  => pcurrent%y
    case('right')
       allocate( masterp(size(right_boundary)) )
       masterp = right_boundary
       coordp  => pcurrent%y
    end select

    if(verbose)write(*,*) 'coordp = ',coordp
    do n = 1 , size(masterp,1)
       if( (n == 1 .and. coordp <= masterp(n)%coordmax) .or. &
           (masterp(n)%coordmin <= coordp .and. coordp <= masterp(n)%coordmax) .or. &
           (n == size(masterp,1) .and. masterp(n)%coordmin <= coordp) )then
!          if(verbose)write(*,*) 'masterp = ',masterp
          select case(masterp(n)%boundary_type)
          case('supersonicin')
             boundary_node%rho = masterp(n)%Q(1)
             boundary_node%p   = masterp(n)%Q(2)
             boundary_node%u   = masterp(n)%Q(3)
             boundary_node%v   = masterp(n)%Q(4)
             boundary_node%A   = masterp(n)%geom(1)
             boundary_node%B   = masterp(n)%geom(2)
             boundary_node%L   = masterp(n)%geom(3)
             boundary_node%M   = masterp(n)%geom(4)
          case('polywall')
             if( .not. allocated(masterp(n)%poly_coeffs) )then
                write(*,*) 'Error in boundary conditions -- polynomial wall called without allocated poly_coeffs array'
                stop
             end if
             polyn = size(masterp(n)%poly_coeffs,1)
             poly_angle = 0.d0
             do inda = 1, polyn
                poly_angle = poly_angle + (polyn - inda + 1)&
                  *masterp(n)%poly_coeffs(inda)*(coordp)**(polyn - inda)
             end do
             poly_angle = atan(poly_angle)
             call wall_reflection(poly_angle)
          case('solidwall')
             call wall_reflection(masterp(n)%wall_angle)
!             boundary_node%rho = pcurrent%rho
!             boundary_node%p   = pcurrent%p
!             boundary_node%u   = cos(2.*masterp(n)%wall_angle)*pcurrent%u + sin(2.*masterp(n)%wall_angle)*pcurrent%v
!             boundary_node%v   = sin(2.*masterp(n)%wall_angle)*pcurrent%u - cos(2.*masterp(n)%wall_angle)*pcurrent%v
!             if(abs(masterp(n)%wall_angle)<=.5*PI)then
!                boundary_node%A   = pcurrent%A
!                boundary_node%M   = (pcurrent%M*(pcurrent%A+boundary_node%A&
!                     *sin(masterp(n)%wall_angle)**2) - sin(masterp(n)%wall_angle)&
!                     *cos(masterp(n)%wall_angle)*(pcurrent%A*pcurrent%L+pcurrent%B&
!                     *pcurrent%M+boundary_node%A*pcurrent%L))&
!                     /(boundary_node%A*(cos(masterp(n)%wall_angle)**2 - &
!                     sin(masterp(n)%wall_angle)**2)-pcurrent%A&
!                     *sin(masterp(n)%wall_angle)**2+pcurrent%B*sin(masterp(n)%wall_angle)&
!                     *cos(masterp(n)%wall_angle))
!                boundary_node%L = tan(masterp(n)%wall_angle)*(pcurrent%M+boundary_node%M)&
!                     - pcurrent%L
!                boundary_node%B = tan(masterp(n)%wall_angle)*(pcurrent%A+boundary_node%A)&
!                     - pcurrent%B
!             elseif(abs(masterp(n)%wall_angle)>.5*PI)then
!                boundary_node%B   = pcurrent%B
!                boundary_node%L   = (pcurrent%L*(pcurrent%B+boundary_node%B&
!                     *cos(masterp(n)%wall_angle)**2) - sin(masterp(n)%wall_angle)&
!                     *cos(masterp(n)%wall_angle)*(pcurrent%A*pcurrent%L+pcurrent%B&
!                     *pcurrent%M+boundary_node%B*pcurrent%M))&
!                     /(boundary_node%B*(sin(masterp(n)%wall_angle)**2 - &
!                     cos(masterp(n)%wall_angle)**2)-pcurrent%B&
!                     *cos(masterp(n)%wall_angle)**2+pcurrent%A*sin(masterp(n)%wall_angle)&
!                     *cos(masterp(n)%wall_angle))
!                boundary_node%A = 1.d0/tan(masterp(n)%wall_angle)*(pcurrent%B+boundary_node%B)&
!                     - pcurrent%A
!                boundary_node%M = 1.d0/tan(masterp(n)%wall_angle)*(pcurrent%L+boundary_node%L)&
!                     - pcurrent%M
!             end if
!             write(*,*) masterp(n)%wall_angle
!             write(*,*) boundary_node
!             read(*,*)

          case('supersonicout')
!					boundary_node%rho = pcurrent%rho
!					boundary_node%p   = pcurrent%p
!					boundary_node%u   = pcurrent%u
!					boundary_node%v   = pcurrent%v
!					boundary_node%A   = pcurrent%A
!					boundary_node%B   = pcurrent%B
!					boundary_node%L   = pcurrent%L
!					boundary_node%M   = pcurrent%M
!					boundary_node%h   = pcurrent%h
             boundary_node = pcurrent
          case('pressure')
             boundary_node = pcurrent
             boundary_node%p = masterp(n)%free_pressure
          case('subsonicout')
             boundary_node = pcurrent
             if(sqrt(boundary_node%u**2+boundary_node%v**2)/sqrt(1.4*boundary_node%p/boundary_node%rho) < 1.d0)then
                write(*,*) "Subsonic" 
                boundary_node%p = masterp(n)%free_pressure
             end if
!             boundary_node%p = masterp(n)%free_pressure
          case('subsonicin')
             boundary_node = pcurrent
             boundary_node%rho = masterp(n)%Q(1)
             boundary_node%u = masterp(n)%Q(3)
             boundary_node%v = masterp(n)%Q(4)
             boundary_node%A   = masterp(n)%geom(1)
             boundary_node%B   = masterp(n)%geom(2)
             boundary_node%L   = masterp(n)%geom(3)
             boundary_node%M   = masterp(n)%geom(4)
          case('rootwall')
             if(directionaltag == 'left' .or. directionaltag == 'down')then
                call wall_reflection(masterp(n)%wall_angle+.5d0*0.1/sqrt(coordp+.1))
             else
                call wall_reflection(masterp(n)%wall_angle-.5d0*0.1/sqrt(coordp+.1))
             end if
          case('logwall')
             if(directionaltag == 'left'.or. directionaltag == 'down')then
                call wall_reflection(masterp(n)%wall_angle+(2.41778+.14*log(coordp+.1))/(18.2699+log(coordp+.01))**2)
             else
                call wall_reflection(masterp(n)%wall_angle-(2.41778+.14*log(coordp+.1))/(18.2699+log(coordp+.01))**2)
             end if
          end select
          exit
       end if
    end do
    deallocate( masterp )
    contains
      subroutine wall_reflection(angle)
      implicit none
         real(8) :: angle

         boundary_node%rho = pcurrent%rho
         boundary_node%p   = pcurrent%p
         boundary_node%u   = cos(2.*angle)*pcurrent%u + sin(2.*angle)*pcurrent%v
         boundary_node%v   = sin(2.*angle)*pcurrent%u - cos(2.*angle)*pcurrent%v
         if(abs(angle)<=.5*PI)then
            boundary_node%A   = pcurrent%A
            boundary_node%M   = (pcurrent%M*(pcurrent%A+boundary_node%A&
                 *sin(angle)**2) - sin(angle)&
                 *cos(angle)*(pcurrent%A*pcurrent%L+pcurrent%B&
                 *pcurrent%M+boundary_node%A*pcurrent%L))&
                 /(boundary_node%A*(cos(angle)**2 - &
                 sin(angle)**2)-pcurrent%A&
                 *sin(angle)**2+pcurrent%B*sin(angle)&
                 *cos(angle))
            boundary_node%L = tan(angle)*(pcurrent%M+boundary_node%M)&
                 - pcurrent%L
            boundary_node%B = tan(angle)*(pcurrent%A+boundary_node%A)&
                 - pcurrent%B
          elseif(abs(angle)>.5*PI)then
            boundary_node%B   = pcurrent%B
            boundary_node%L   = (pcurrent%L*(pcurrent%B+boundary_node%B&
                 *cos(angle)**2) - sin(angle)&
                 *cos(angle)*(pcurrent%A*pcurrent%L+pcurrent%B&
                 *pcurrent%M+boundary_node%B*pcurrent%M))&
                 /(boundary_node%B*(sin(angle)**2 - &
                 cos(angle)**2)-pcurrent%B&
                 *cos(angle)**2+pcurrent%A*sin(angle)&
                 *cos(angle))
             boundary_node%A = 1.d0/tan(angle)*(pcurrent%B+boundary_node%B)&
                  - pcurrent%A
             boundary_node%M = 1.d0/tan(angle)*(pcurrent%L+boundary_node%L)&
                  - pcurrent%M
          end if
        end subroutine wall_reflection
  end function boundary_node

end module boundary_conditions

!program boundary_conditions_tester
!	use boundary_conditions_init
!	use boundary_conditions
!	use types
!    implicit none
!    real(8) temp
!    type (node_data) :: test_node, boundnode
!    logical :: verbose = .true.
!    call read_boundary(verbose)
!
!
!
!end program boundary_conditions_tester
