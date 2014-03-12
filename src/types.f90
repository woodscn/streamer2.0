!__****h* main/types
!__ NAME
!__		Types -- module containing type definitions used in ucs2d
!__	AUTHOR
!__		C. Nathan Woods
!__	USAGE
!__		Load the desired type definitions using the use ___, only: ___
!__		construct.
!__	EXAMPLE
!__		use types, only: node_data
!__		The above will load the type definition for node_data only,
!__		leaving the types coord_init, etc. undefined.
!__****

module types
  type node_data
!__****t* node_array/node_data
!__ NAME
!__		Node_Data -- Primary data type. Stores all nodal information except for
!__ 	pointers.
!__ DESCRIPTION
!__		Node_data is a structure intended to hold all nodal information
!__		necessary to define the simulation. Node_data contains the
!__		following information:
!__		* Primitive flow variables
!__		  (density, pressure, x- and y-velocity components)
!__		* Metric components (A, B, L, M)
!__		* Node position, in Cartesian coordinates (x, y)
!__		* Grid-motion function (h, g)
!__             * Flow speed (q)
!__		* Flow angle, measured ccw from positive x-axis (theta)
!__ NOTES
!__		Rather than conservative variables, primitive flow variables are
!__		stored because they are used more frequently. The conservative
!__		variables are only used to actually update the solution in time.
!__
!__		The grid-motion function is stored as both g and h, in order
!__		to better accomodate orthogonality-preserving conditions,
!__		especially for stagnation flows. 
!__
!__ TODO
!__ 	
!__****
     real(8) :: rho = 0
     real(8) :: p = 0
     real(8) :: u = 0
     real(8) :: v = 0
     real(8) :: A = 0
     real(8) :: B = 0
     real(8) :: L = 0
     real(8) :: M = 0
     real(8) :: x = 0
     real(8) :: y = 0
     real(8) :: g = 0
     real(8) :: h = 0
     real(8) :: q = 0
     real(8) :: theta = 0
     logical :: boundary = .false.
   contains
     procedure, pass :: primtocons
     procedure, pass :: constoprim
     procedure, pass :: add_node_data
     procedure, pass :: subtract_node_data
     procedure, pass :: scalar_node_data_multiply
     procedure, nopass :: node_data_scalar_multiply
     procedure, pass :: node_data_node_data_equivalance
     procedure, pass :: node_data_integer_exponentiation
     procedure, pass :: rotate
     procedure, pass :: node_data_max
!     procedure, pass :: ndnd
     generic :: operator(+) => add_node_data
     generic :: operator(-) => subtract_node_data
!     generic :: operator(*) => scalar_node_data_multiply, node_data_scalar_multiply
     generic :: operator(==) => node_data_node_data_equivalance
     generic :: operator(**) => node_data_integer_exponentiation
     generic :: max => node_data_max
!     generic :: assignment(=) => ndnd
  end type  node_data
  
  interface operator(*)
     module procedure scalar_node_data_multiply
     module procedure node_data_scalar_multiply
  end interface operator(*)

  type, extends(node_data) :: node
!__****t*
!__ NAME
!__ 	Node -- Extension of node_data, including pointers to other nodes.
!__****
     type (node), pointer :: up => null()
     type (node), pointer :: down => null()
     type (node), pointer :: left => null()
     type (node), pointer :: right => null()
     type (node), pointer :: next => null()
  end type node

!__****t* types/coord_init
!__ NAME
!__		Coord_Init -- type definition for holding (x,y) coordinate data
!__****
  type coord_init
     real(8) :: x
     real(8) :: y
  end type coord_init

!__****t* types/boundary_data , types/boundary_master
!__ NAME
!__		Boundary_Master -- type definition used to define simulation boundaries
!__ DESCRIPTION
!__		The derived type boundary_data contains information about
!__      	boundary conditions at a node. Originally intended to be
!__     	contained in the node_data type, it is currently used only in
!__     	its extended form boundary_master as part of overall simulation
!__ 	        boundary conditions. The connection between nodes and their
!__      	boundary conditions is now described entirely using the
!__     	association status of node pointers and the (x,y) coordinate
!__     	position of the node.
!__ TODO
!__		Since the details of boundary conditions are no longer stored
!__		at simulation nodes, the definition for boundary_data should be
!__		rolled into boundary_master. Also, the 'corner' type is no
!__		longer necessary, nor are boundary_type1 and boundary_type2.
!__****
  type boundary_data
     character(len=16) :: boundary_type = ''
     !>  interior, supersonicout, solidwall, pressure, corner
     character(len=16) :: boundary_type1= ''! For corners only, since a corner will have 2 BC's
     character(len=16) :: boundary_type2= ''
     real(8) :: wall_angle = 0.
     real(8) :: free_pressure = 0.
     real(8) :: Q(4) = (/ 0.0, 0.0, 0.0, 0.0 /)
     real(8) :: geom(4) = [ 1.0, 0.0, 0.0, 1.0 ]
     real(8), allocatable :: poly_coeffs(:)
  end type boundary_data

  type, extends(boundary_data) :: boundary_master
  real(8) :: coordmin = 0.0
  real(8) :: coordmax = 0.0
end type boundary_master


contains

!  subroutine ndnd(out,in)
!    implicit none
!    type(node_data), intent(in) :: in
!    type(node_data), intent(out):: out
!
!    out%rho = in%rho
!    out%p   = in%p
!    out%u   = in%u
!    out%v   = in%v
!    out%A   = in%A
!    out%B   = in%B
!    out%L   = in%L
!    out%M   = in%M
!    out%x   = in%x
!    out%y   = in%y
!    out%g   = in%g
!    out%h   = in%h
!    out%q   = in%q
!    out%theta=in%theta
!    out%boundary=in%boundary
!  end subroutine ndnd

  function add_node_data(node1,node2)
    implicit none
    class(node_data), intent(in) :: node1, node2
    type(node_data) :: add_node_data
!    type(node_data) ::  out

    add_node_data%rho = node1%rho + node2%rho
    add_node_data%p   = node1%p   + node2%p
    add_node_data%u   = node1%u   + node2%u
    add_node_data%v   = node1%v   + node2%v
    add_node_data%A   = node1%A   + node2%A
    add_node_data%B   = node1%B   + node2%B
    add_node_data%L   = node1%L   + node2%L
    add_node_data%M   = node1%M   + node2%M
    add_node_data%h   = node1%h   + node2%h
    add_node_data%g   = node1%g   + node2%g
    add_node_data%x   = node1%x   + node2%x
    add_node_data%y   = node1%y   + node2%y
    add_node_data%q   = node1%q   + node2%q
    add_node_data%theta = node1%theta + node2%theta
    add_node_data%boundary = .false.

!    add_node_data = out

  end function add_node_data

  function subtract_node_data(node1,node2)
    implicit none
    class(node_data), intent(in) :: node1, node2
    type(node_data) :: subtract_node_data

    subtract_node_data%rho = node1%rho - node2%rho
    subtract_node_data%p   = node1%p   - node2%p
    subtract_node_data%u   = node1%u   - node2%u
    subtract_node_data%v   = node1%v   - node2%v
    subtract_node_data%A   = node1%A   - node2%A
    subtract_node_data%B   = node1%B   - node2%B
    subtract_node_data%L   = node1%L   - node2%L
    subtract_node_data%M   = node1%M   - node2%M
    subtract_node_data%h   = node1%h   - node2%h
    subtract_node_data%g   = node1%g   - node2%g
    subtract_node_data%x   = node1%x   - node2%x
    subtract_node_data%y   = node1%y   - node2%y
    subtract_node_data%q   = node1%q   - node2%q
    subtract_node_data%theta = node1%theta - node2%theta
    subtract_node_data%boundary = .false.

  end function subtract_node_data

  function scalar_node_data_multiply(node1,scalar)
    implicit none
!    type(node_data), intent(in) :: node1
    class(node_data), intent(in) :: node1
    type(node_data) :: scalar_node_data_multiply
    real(8), intent(in) :: scalar

    scalar_node_data_multiply%rho = node1%rho * scalar
    scalar_node_data_multiply%p   = node1%p   * scalar
    scalar_node_data_multiply%u   = node1%u   * scalar
    scalar_node_data_multiply%v   = node1%v   * scalar
    scalar_node_data_multiply%A   = node1%A   * scalar
    scalar_node_data_multiply%B   = node1%B   * scalar
    scalar_node_data_multiply%L   = node1%L   * scalar
    scalar_node_data_multiply%M   = node1%M   * scalar
    scalar_node_data_multiply%h   = node1%h   * scalar
    scalar_node_data_multiply%g   = node1%g   * scalar
    scalar_node_data_multiply%x   = node1%x   * scalar
    scalar_node_data_multiply%y   = node1%y   * scalar
    scalar_node_data_multiply%q   = node1%q   * scalar
    scalar_node_data_multiply%theta = node1%theta * scalar
    scalar_node_data_multiply%boundary = .false.

  end function scalar_node_data_multiply

  function node_data_scalar_multiply(scalar,node1)
    implicit none
!    type(node_data), intent(in) :: node1
    class(node_data), intent(in) :: node1
    type(node_data) :: node_data_scalar_multiply
    real(8), intent(in) :: scalar

    node_data_scalar_multiply%rho = node1%rho * scalar
    node_data_scalar_multiply%p   = node1%p   * scalar
    node_data_scalar_multiply%u   = node1%u   * scalar
    node_data_scalar_multiply%v   = node1%v   * scalar
    node_data_scalar_multiply%A   = node1%A   * scalar
    node_data_scalar_multiply%B   = node1%B   * scalar
    node_data_scalar_multiply%L   = node1%L   * scalar
    node_data_scalar_multiply%M   = node1%M   * scalar
    node_data_scalar_multiply%h   = node1%h   * scalar
    node_data_scalar_multiply%g   = node1%g   * scalar
    node_data_scalar_multiply%x   = node1%x   * scalar
    node_data_scalar_multiply%y   = node1%y   * scalar
    node_data_scalar_multiply%q   = node1%q   * scalar
    node_data_scalar_multiply%theta = node1%theta * scalar

  end function node_data_scalar_multiply

  function node_data_integer_exponentiation(node, scalar)
    implicit none
    class(node_data), intent(in) :: node
    integer, intent(in) :: scalar
    type(node_data) :: node_data_integer_exponentiation
    
    node_data_integer_exponentiation%p = node%p**scalar
    node_data_integer_exponentiation%rho = node%rho**scalar
    node_data_integer_exponentiation%u = node%u**scalar
    node_data_integer_exponentiation%v = node%v**scalar
    node_data_integer_exponentiation%A = node%A**scalar
    node_data_integer_exponentiation%B = node%B**scalar
    node_data_integer_exponentiation%L = node%L**scalar
    node_data_integer_exponentiation%M = node%M**scalar
    node_data_integer_exponentiation%h = node%h**scalar
    node_data_integer_exponentiation%g = node%g**scalar
    node_data_integer_exponentiation%q = node%q**scalar
    node_data_integer_exponentiation%theta = node%theta**scalar
    node_data_integer_exponentiation%x = node%x**scalar
    node_data_integer_exponentiation%y = node%y**scalar
    node_data_integer_exponentiation%boundary = .false.
  end function node_data_integer_exponentiation
    
  function node_data_node_data_equivalance(node1, node2)
    implicit none
    class(node_data), intent(in) :: node1, node2
    logical :: node_data_node_data_equivalance
    
    node_data_node_data_equivalance = (node1%rho == node2%rho) .and. &
         (node1%p == node2%p) .and. (node1%u == node2%u) .and. & 
         (node1%v == node2%v) .and. (node1%A == node2%A) .and. &
         (node1%B == node2%B) .and. (node1%L == node2%L) .and. &
         (node1%M == node2%M) .and. (node1%h == node2%h) .and. &
         (node1%g == node2%g) .and. (node1%x == node2%x) .and. &
         (node1%y == node2%y) .and. (node1%q == node2%q) .and. &
         (node1%theta == node2%theta) .and. (node1%boundary .eqv. node2%boundary)
  end function node_data_node_data_equivalance

  function node_data_max(node)
    implicit none
    class(node_data), intent(in) :: node
    real(8) :: node_data_max

    node_data_max = maxval([node%p, node%rho, node%u, node%v, &
         node%A, node%B, node%L, node%M, node%h, node%g, node%x, &
         node%y, node%q, node%theta])
  end function node_data_max

  function primtocons(p)
    use global_data, only: gamma_const
    implicit none
    type(node_data) :: primtocons
    class(node_data), intent(in) :: p
    real(8) :: delta, energy

    delta = p%A*p%M - p%B*p%L
    energy = 0.5d0*(p%u**2+p%v**2)+p%p/(p%rho*(gamma_const-1.d0))
    
    primtocons%rho = p%rho*delta
    primtocons%p   = p%rho*delta*energy
    primtocons%u   = p%rho*delta*p%u
    primtocons%v   = p%rho*delta*p%v
    primtocons%A   = p%A
    primtocons%B   = p%B
    primtocons%L   = p%L
    primtocons%M   = p%M
    primtocons%h   = p%h
    primtocons%g   = p%g
    primtocons%q   = p%q
    primtocons%theta=p%theta
    primtocons%x   = p%x
    primtocons%y   = p%y
    primtocons%boundary = p%boundary
  end function primtocons

  function constoprim(p)
    use global_data, only: gamma_const
    implicit none
    type(node_data) :: constoprim
    class(node_data), intent(in) :: p

    constoprim%rho = p%rho/(p%A*p%M-p%B*p%L)
    constoprim%u   = p%u/p%rho
    constoprim%v   = p%v/p%rho
    constoprim%p   =(p%p/p%rho - .5d0*(constoprim%u**2+constoprim%v**2))*(gamma_const - 1.d0)*constoprim%rho
    constoprim%A   = p%A
    constoprim%B   = p%B
    constoprim%L   = p%L
    constoprim%M   = p%M
    constoprim%h   = p%h
    constoprim%g   = p%g
    constoprim%q   = p%q
    constoprim%theta=p%theta
    constoprim%x   = p%x
    constoprim%y   = p%y
    constoprim%boundary = p%boundary

  end function constoprim

  function rotate(in,normal,inverse)
    implicit none
    type(node_data) :: rotate
    class(node_data), intent(in) :: in
    logical, intent(in) :: inverse
    real(8), dimension(2) :: normal

    rotate%p = in%p
    rotate%rho = in%rho
!    rotate%u = in%u
!    rotate%v = in%v
    rotate%A = in%A
    rotate%B = in%B
    rotate%L = in%L
    rotate%M = in%M
    rotate%h = in%h
    rotate%g = in%g
    rotate%q = in%q
    rotate%theta = in%theta
    rotate%x = in%x
    rotate%y = in%y
    rotate%boundary = in%boundary
    if(.not. inverse)then
       rotate%u = normal(1)*in%u + normal(2)*in%v
       rotate%v =-normal(2)*in%u + normal(1)*in%v
    else
       rotate%u = normal(1)*in%u - normal(2)*in%v
       rotate%v = normal(2)*in%u + normal(1)*in%v
    end if
  end function rotate

!!$  subroutine types_test()
!!$    implicit none
!!$    type (node_data), pointer :: head , tail
!!$    type (boundary_master) :: master_test
!!$    integer :: STAT ,istat
!!$
!!$    nullify( head, tail )
!!$    allocate( head , STAT=istat )
!!$    if( STAT /= 0 )then
!!$       write(*,*) 'Error allocating node_data in types_test'
!!$       stop
!!$    end if
!!$    tail => head
!!$    write(*,*) 'Is head assocated?' , associated(head)
!!$    head = node_data(1.0,1.0,1.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0)
!!$    write(*,*) 'Density = ' , head%rho
!!$    write(*,*) 'Pressure = ' , head%p
!!$    write(*,*) 'Velocity = ' , head%u , head%v
!!$    write(*,*) 'Metric   = ' , head%A, head%B, head%L, head%M
!!$    write(*,*) 'Position = ' , head%x, head%y
!!$    write(*,*) 'Grid motion: ', head%g, head%theta
!!$    write(*,*) 'Boundary type: ', head%boundary
!!$    allocate( head%down , STAT=istat )
!!$    if( STAT /= 0 )then
!!$       write(*,*) 'Error allocating node_data in types_test'
!!$       stop
!!$    end if
!!$    tail => head%down
!!$    tail = node_data(1.0,1.0,1.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,2.0,boundary_data('free'))
!!$    write(*,*) 'Density = ' , tail%rho
!!$    write(*,*) 'Pressure = ' , tail%p
!!$    write(*,*) 'Velocity = ' , tail%u , tail%v
!!$    write(*,*) 'Metric   = ' , tail%A, tail%B, tail%L, tail%M
!!$    write(*,*) 'Position = ' , tail%x, tail%y
!!$    write(*,*) 'Grid motion: ', tail%g, tail%theta
!!$    write(*,*) 'Boundary type: ', tail%boundary
!!$
!!$    write(*,*) 'Master_test = ' , master_test
!!$
!!$    write(*,*) head%primtocons()
!!$  end function types_test

  function node_tester()
!__****f*
!__ NAME
!__     Node_tester -- Exercise and verify node type and associated operators.
!__****
    implicit none
    type(node_data) :: p1, p2, p3, p3cons, temp, zero, one
    type(node) :: n1, n2
    real(8) :: angle1, angle2
    real(8), dimension(2) :: normal
    integer :: node_tester
    node_tester = 1
    zero = node_data(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,&
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .false.)
    one  = node_data(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,&
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, .false.)
    p1%rho = 1.0 ; p2%rho = 0.75        ; p3%rho = .25                    ; p3cons%rho = .155
    p1%p   = 1.0 ; p2%p   = 0.5         ; p3%p   = .5                     ; p3cons%p   = .81936875
    p1%u   = 1.0 ; p2%u   = 0.6         ; p3%u   = .75                    ; p3cons%u   = .11625
    p1%v   = 0.0 ; p2%v   = 0.8         ; p3%v   = .1                     ; p3cons%v   = .0155
    p1%A   = 1.0 ; p2%A   = 1./sqrt(2.) ; p3%A   = .8                     ; p3cons%A   = .8
    p1%B   = 0.0 ; p2%B   = 1./sqrt(2.) ; p3%B   = .3                     ; p3cons%B   = .3
    p1%L   = 0.0 ; p2%L   = 1./sqrt(2.) ; p3%L   =-.2                     ; p3cons%L   =-.2
    p1%M   = 1.0 ; p2%M   = 1./sqrt(2.) ; p3%M   = .7                     ; p3cons%M   = .7
    p1%x   = 0.0 ; p2%x   = 1.35        ; p3%x   = .5                     ; p3cons%x   = .5
    p1%y   = 0.0 ; p2%y   = 1.25        ; p3%y   = .6                     ; p3cons%y   = .6
    p1%h   = 1.0 ; p2%h   = 0.25        ; p3%h   = .25                    ; p3cons%h   = .25
    p1%g   = 1.0 ; p2%g   = 0.25        ; p3%g   = .25*sqrt(.75**2+.1**2) ; p3cons%g   = .25*sqrt(.75**2+.1**2)
    p1%q   = 1.0 ; p2%q   = 1.          ; p3%q   =     sqrt(.75**2+.1**2) ; p3cons%q   =     sqrt(.75**2+.1**2)
    p1%theta=0.0 ; p2%theta=0.          ; p3%theta=   atan2(.1,.75)       ; p3cons%theta=   atan2(.1,.75)
    p1%boundary = .false. ; p2%boundary = .false. ; p3%boundary = .false. ; p3cons%boundary = .false.
!!$    p1 = node_data(1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0,&
!!$         1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, .false.)
!!$    p2 = node_data(.75, .5, .6, .8, 1.0/sqrt(2.), 1.0/sqrt(2.),&
!!$         -1.0/sqrt(2.), 1.0/sqrt(2.), 1.35, 1.25, .25, .25, 1.0,&
!!$         0.0, .false.)
!!$    p3 = node_data(.25, .5, .75, .1, .8, .3, -.2, .7,&
!!$         .5, .6, .25, .25*sqrt(.75**2+.1**2), sqrt(.75**2+.1**2),&
!!$         atan2(.1,.75), .false.)
!!$    p3cons = node_data(.155, .775, .11625, .0155, .8, .3, -.2, .7&
!!$         , .5, .6, .25, .25*sqrt(.75**2+.1**2), &
!!$         sqrt(.75**2+.1**2), atan2(.1,.75), .false.)
    angle1 = .1d0
    angle2 = .33d0
    ! Primtocons/Constoprim inverse test
    temp = p2%constoprim()
    temp = p2 - temp%primtocons()
    if( abs(temp%max()) <1d-15 )then
       node_tester = node_tester * 1
    else 
       node_tester = node_tester * 0
    end if
    ! Primtocons verification test
    write(*,*) node_tester
    if( abs(temp%max()) < 1d-15 )then
       node_tester = node_tester * 1
    else
       node_tester = node_tester * 0
    end if
    ! Additive inverse test
    if( p2 - p2 == zero .and. p1 - p1 == zero )then
       node_tester = node_tester * 1
    else
       node_tester = node_tester * 0
    end if
    ! Additive verification test
    temp%rho = p3%rho + p3cons%rho
    temp%p   = p3%p   + p3cons%p
    temp%u   = p3%u   + p3cons%u
    temp%v   = p3%v   + p3cons%v
    temp%A   = p3%A   + p3cons%A
    temp%B   = p3%B   + p3cons%B
    temp%L   = p3%L   + p3cons%L
    temp%M   = p3%M   + p3cons%M
    temp%x   = p3%x   + p3cons%x
    temp%y   = p3%y   + p3cons%y
    temp%h   = p3%h   + p3cons%h
    temp%g   = p3%g   + p3cons%g
    temp%theta=p3%theta+p3cons%theta
    temp%q   = p3%q   + p3cons%q
    temp%boundary = p3%boundary
    temp = ( temp - ( p3 + p3cons ) )**2
    if( abs(temp%max()) < 1d-15 )then
       node_tester = node_tester * 1
    else
       node_tester = node_tester * 0
    end if
    ! Multiplicative inverse test
    temp = (p2 - (1.d0/1.992d0)*(p2*1.992d0))**2
    if( abs(temp%max()) < 1d-15 )then
       node_tester = node_tester * 1
    else
       node_tester = node_tester * 0
    end if
    ! Rotative inverse test
    normal = [cos(angle1), sin(angle1)]
    temp = p2%rotate(normal,.false.)
    temp = (p2-temp%rotate(normal,.true.))**2
    if( temp%max() < 1d-15 )then
       node_tester = node_tester * 1
    else
       node_tester = node_tester * 0
    end if
    
    write (*,*) "Node_tester = ",node_tester
  end function node_tester

  integer function types_tester()
    implicit none
    integer success
    success = 0
    success = success + node_tester()
    types_tester = success
  end function types_tester

end module types
