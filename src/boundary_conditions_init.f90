!__****h* init/boundary_conditions_init
!__ NAME
!__ 	Boundary_Conditions_Init -- Module containing routines for initializing
!__ 	boundary conditions
!__ AUTHOR
!__ 	C. Nathan Woods
!__****

module boundary_conditions_init
  use types, only: boundary_master, coord_init
  use global_data, only: gamma_const
  type (boundary_master), target, save, dimension(:), allocatable, protected ::    top_boundary
  type (boundary_master), target, save, dimension(:), allocatable, protected :: bottom_boundary
  type (boundary_master), target, save, dimension(:), allocatable, protected ::   left_boundary
  type (boundary_master), target, save, dimension(:), allocatable, protected ::  right_boundary
  type (coord_init), dimension(:), allocatable, protected :: leftcoords

contains
!__****f* boundary_conditions_init/read_boundary
!__ NAME
!__		Read_Boundary -- Reads boundary data from file and stores
!__ 	it for use in the simulation
!__ INPUTS
!__ 	boundary_input.in -- Data file describing boundary conditions
!__ OUTPUTS
!__ 	top_boundary -- Array of type boundary_master describing eta = neta+1 boundary
!__ 	bottom_boundary -- Array of type boundary_master describing eta = 0 boundary
!__ 	left_boundary -- Array of type boundary_master describing xi = 0 boundary
!__ 	right_boundary -- Array of type boundary_master describing xi = nxi+1 boundary
!__ 	leftcoords -- Array of type coord_init describing left-most (x,y) coordinates
!__ 	window -- Structure describing the extent of the computational "window"
!__ 	differentials -- Two-element array containing differential lengths (dxi, deta)
!__ 	dimensions -- Two-element array containing dimensions (nxi, neta)
!__ DESCRIPTION
!__ 	Read_boundary reads input data from boundary_input.in, located in
!__ 	the main program directory (in this case ../ucs2d2_0, one directory ABOVE the
!__ 	/src directory that holds source files), and writes that data into persistent
!__ 	data structures that can be used (read-only) throughout the simulation. Each
!__ 	of four boundary structures, labelled top_, bottom_, left_, and right_boundary,
!__ 	is an allocatable array of boundary_master elements. Each element contains:
!__ 	* A coordinate range of the form [xmin, xmax] to specify a boundary segment
!__ 	* A descriptor specifying the type of boundary condition to be enforced on
!__ 		the boundary segment, e.g. supersonic inflow, a solid wall, etc.
!__ 	* Any further data necessary to completely specify the particular boundary
!__ 		condition
!__
!__ 	The extent of the computational window is derived from the maximum extent of
!__ 	the prescribed boundary conditions.
!__
!__ 	The differential lengths (dxi, deta) and simulation dimensions (nxi, neta) are
!__ 	also read from boundary_input.in. Only one set of values should be specified;
!__ 	the other is computed based on either desired differential lengths or on a
!__ 	desired nominal grid resolution.
!__****
  subroutine read_boundary(verbose_flag)
    use global_data, only: set_window, set_differentials, set_dimensions, &
      differentials, dimensions, window, set_splitting_type, input_file_name, &
      set_h0, set_CFL, set_skip, set_tmax, set_update_type, set_output_file_name, &
      set_grid_preserving, set_muscl_flag
    implicit none
    integer :: ierror, nbounds, n, side, nxi, neta
    character(len=128) :: btype, error_message, output_file_name!, buffer
    real(8) :: coordmin, coordmax, Wbound(4), theta
    logical, optional :: verbose_flag
    logical :: verbose = .false.
    type (boundary_master), pointer, dimension(:) :: bmp
    real(8), dimension(4) :: mins, maxes
    real(8) :: dxi, deta, tmax, h0, CFL
    character(len=16) :: split_type
    integer :: skip, update_type
    logical :: grid_preserving, muscl_flag
    integer :: poly_order
    
    if( present(verbose_flag))verbose=verbose_flag
    open(UNIT=8, FILE=trim(input_file_name), STATUS='old', &
         ACTION='read', IOSTAT=ierror, IOMSG=error_message)
    if( ierror /= 0 )then;
       write(*,*) error_message;
       stop;
    end if;
    
    do side = 1 , 4

       read(UNIT=8,FMT='(I1)',IOSTAT=ierror,IOMSG=error_message) nbounds
       if( ierror /= 0 )then;
          write(*,*) error_message;
          stop;
       end if;
       if(verbose)write(*,*) 'nbounds = ',nbounds
       !			allocate( bmp(nbounds) )


       select case(side)
       case(1)
          allocate( top_boundary(nbounds), STAT=ierror )
          bmp => top_boundary
       case(2)
          allocate( bottom_boundary(nbounds), STAT=ierror )
          bmp => bottom_boundary
       case(3)
          allocate( left_boundary(nbounds), STAT=ierror )
          bmp => left_boundary
       case(4)
          allocate( right_boundary(nbounds), STAT=ierror )
          bmp => right_boundary
       end select
       if(ierror/=0)then
          write(*,*) 'Error in boundary_conditions_init -- allocating boundary_master arrays'
          stop
       end if

       do n = 1 , nbounds
          read(8,*) coordmin, coordmax, btype
          bmp(n)%coordmin = coordmin
          bmp(n)%coordmax = coordmax
          bmp(n)%boundary_type = trim(btype)
          select case (btype)
          case('supersonicin')
             read(8,*) bmp(n)%Q, bmp(n)%geom
             bmp(n)%Q(3) = bmp(n)%Q(3)*sqrt(gamma_const*bmp(n)%Q(2)/bmp(n)%Q(1))*cos(bmp(n)%Q(4))
             bmp(n)%Q(4) = bmp(n)%Q(3)*sqrt(gamma_const*bmp(n)%Q(2)/bmp(n)%Q(1))*sin(bmp(n)%Q(4))
          case('subsonicin')
             read(8,*) bmp(n)%Q, bmp(n)%geom
             bmp(n)%Q(3) = bmp(n)%Q(3)*sqrt(gamma_const*bmp(n)%Q(2)/bmp(n)%Q(1))*cos(bmp(n)%Q(4))
             bmp(n)%Q(4) = bmp(n)%Q(3)*sqrt(gamma_const*bmp(n)%Q(2)/bmp(n)%Q(1))*sin(bmp(n)%Q(4))
          case('solidwall')
             read(8,*) theta
             if(verbose)write(*,*) 'Wall angle theta = ',theta
             bmp(n)%wall_angle=theta
          case('rootwall')
             read(8,*) theta
             if(verbose)write(*,*) 'Wall angle theta = ',theta
             bmp(n)%wall_angle=theta
          case('logwall')
             read(8,*) theta
             if(verbose)write(*,*) 'Wall angle theta = ',theta
             bmp(n)%wall_angle=theta
          case('pressure')
             read(8,*) Wbound(2)
             if(verbose)write(*,*) 'Freestream pressure = ',Wbound(2)
             bmp(n)%free_pressure = Wbound(2)
          case('subsonicout')
             read(8,*) Wbound(2)
             if(verbose)write(*,*) 'Freestream pressure = ',Wbound(2)
             bmp(n)%free_pressure = Wbound(2)
          case('polywall')
             read(8,*) poly_order
             allocate(bmp(n)%poly_coeffs(poly_order))
             read(8,*) bmp(n)%poly_coeffs
          end select
!          if(verbose)write(*,*) 'bmp = ', bmp(n)
       end do
       if(nbounds>1)then
          mins (side) = minval(bmp%coordmin)
          maxes(side) = maxval(bmp%coordmax)
       else
          mins (side) = bmp(1)%coordmin
          maxes(side) = bmp(1)%coordmax
       end if
    end do
    if(mins(1)==mins(2) .and. maxes(1)==maxes(2) .and. mins(3)==mins(4) .and. maxes(3)==maxes(4) )then
!       window%xmin = mins(1)
!       window%xmax = maxes(1)
!       window%ymin = mins(3)
!       window%ymax = maxes(3)
       call set_window([mins(1),maxes(1),mins(3),maxes(3)])
       if(verbose)write(*,*) window
    else
       write(*,*) 'Boundary conditions are inconsistent'
       stop
    end if

!    write(*,*) top_boundary
!    write(*,*) bottom_boundary
!    write(*,*) left_boundary
!    write(*,*) right_boundary

    read(8,*) nxi,  dxi
    read(8,*) neta, deta
    if(dxi == 0.)dxi = (window%xmax-window%xmin)/real(nxi )
    if(deta== 0.)deta= (window%ymax-window%ymin)/real(neta)
    if(nxi == 0 )nxi = nint((window%xmax-window%xmin)/dxi )! Maybe round down?
    if(neta== 0 )neta= nint((window%ymax-window%ymin)/deta)

    if(verbose)write(*,*) 'dxi = ', dxi, 'deta = ', deta
    if(verbose)write(*,*) 'nxi = ', nxi, 'neta = ', neta

!    differentials(1) = dxi
!    differentials(2) = deta
    call set_differentials([dxi,deta])

!    dimensions(1) = 0
    !		dimensions(1) = nxi
!    dimensions(2) = neta

    call set_dimensions([0,neta])

    n = 0

    allocate( leftcoords(neta), STAT=ierror )
    if( ierror /= 0 )then;
       write(*,*) error_message;
       stop;
    end if;
    do
       n = n + 1
       if(window%ymin+(n-.5)*deta > window%ymax)exit
       leftcoords(n)%y = window%ymin + (n-.5)*deta
       leftcoords(n)%x = window%xmin
    end do
    if(size(leftcoords%y,1) /= neta)then;
       write(*,*) 'Error: neta and deta are incompatible';
       stop;
    end if;
    if(verbose) write(*,*) 'xleft = ', leftcoords%x
    if(verbose) write(*,*) 'yleft = ', leftcoords%y

    read(8,*) split_type
    call set_splitting_type(trim(split_type))
    read(8,*) h0, grid_preserving
    call set_h0(h0)
    call set_grid_preserving(grid_preserving)
    read(8,*) CFL
    call set_CFL(CFL)
    read(8,*) tmax
    call set_tmax(tmax)
    read(8,*) skip
    call set_skip(skip)
    read(8,*) update_type
    call set_update_type(update_type)
    read(8,*) output_file_name
    call set_output_file_name(output_file_name)
    read(8,*) muscl_flag
    call set_muscl_flag(muscl_flag)

  end subroutine read_boundary

  integer function boundary_test_case(verbose_flag)
    use types
    use global_data
    implicit none

    logical, optional, intent(in) :: verbose_flag
    logical :: verbose
    if(present(verbose_flag))then
       verbose = verbose_flag
    else
       verbose = .false.
    end if

    open(UNIT=88, FILE='test_boundaries.in', &
         ACTION='write', IOSTAT=ierror, IOMSG=error_message)
    write(88,*) "1"
    write(88,*) "0.0 0.6 supersonicout"
    write(88,*) "1"
    write(88,*) "0.0 0.6 supersonicout"
    write(88,*) "2"
    write(88,*) "0.0 0.5 supersonicin"
    write(88,*) "1.0 1.0 2.4 0.0 1.0 0.0 0.0 1.0"
    write(88,*) "0.5 1.0 supersonicin"
    write(88,*) "0.5 0.25 7.0 0.0 1.0 0.0 0.0 1.0"
    write(88,*) "1"
    write(88,*) "0.0 1.0 supersonicout"
    write(88,*) "22 0."
    write(88,*) "20 0."
    write(88,*) "strang"
    write(88,*) "0.000000 .false."
    write(88,*) "0.700000"
    write(88,*) "2.000000"
    write(88,*) "1000"
    write(88,*) "2"
    write(88,*) "test_boundaries.dat"
    write(88,*) "true"

  end function boundary_init_tester

end module boundary_conditions_init
