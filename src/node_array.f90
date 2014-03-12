!__****h* main/node_array
!__ Name
!__		Node array -- module for array-like operations for a
!__		linked list of node data elements.
!__ DESCRIPTION
!__		Node_array was originally conceived as a repository
!__		for procedures that would be used to handle array-
!__		like operation on a linked list of node_data elements.
!__		In its current state, Node_array contains only the
!__		subroutine for allocating a new column of nodes.
!__ TODO
!__		Eventually, this will also contain a routine to
!__		deallocate nodes.
!__****

module node_array
  use types, only: node, node_data
  type (node), pointer, private, save :: head => null()

  interface assignment (=)
     module procedure node_node_data_assign
  end interface
  interface operator (.copy.)
     module procedure ndnc
  end interface


contains

  subroutine create_column(verbose_flag,neta_in,xi_in,eta_in)

    implicit none
    logical, optional :: verbose_flag
    integer, optional :: neta_in
    real(8), optional :: xi_in
    real(8), optional, dimension(:) :: eta_in
    call allocate_column(verbose_flag,neta_in)
    if(present(xi_in).and.present(eta_in))then
      call initialize_column(verbose_flag,.true.,xi=xi_in,eta_in = eta_in)
    elseif(present(xi_in))then
      call initialize_column(verbose_flag,.true.,xi=xi_in)
    else
      call initialize_column(verbose_flag,.false.)
    end if
  end subroutine create_column

  subroutine allocate_column(verbose_flag,neta_in)
    use global_data, only: dimensions, set_dimensions

    implicit none
    logical, optional, intent(in) :: verbose_flag
    logical :: verbose = .false.
    integer, intent(in), optional :: neta_in
    type (node), pointer :: p1 => null(), p2 => null()
    integer :: istat , i, neta

    neta = dimensions(2)
    if(present(neta_in)) neta = neta_in
    if(present(verbose_flag)) verbose = verbose_flag

    if(verbose) write(*,*) 'Neta used for testing: ',neta

    if(neta < 1)then
       write(*,*) 'Error in create_column -- inappropriate dimension specification'
       stop
    end if

!		Case for creating new column to attach to already-extant array.
    if(associated(head))then
       if(verbose) write(*,*) 'Head pointer is already associated. Adding new column.'
       allocate(head%left,STAT=istat)
       if(istat > 0)then
          write(*,*) 'Error in create_column -- allocating new column base'
          stop
       end if
       head%left%right => head ! Aim reciprocal pointer back to head.
!			Aim appropriate pointers: p1 to bottom element of column to be
!			created, and p2 to bottom element of column to the right of p1.
       p1 => head%left ; p2 => head
!			Reassign head to point to new bottom-left element.
       head => p1
!			Loop to create remainder of new column and assign up-down and left-
!			right links.
       do
          allocate(p1%up,STAT=istat)
          if(istat > 0)then
             write(*,*) 'Error in create_column -- allocating new column'
             stop
          end if
          p1%up%down => p1
          p1 =>p1%up ; p2=> p2%up
          p1%right => p2 ; p2%left => p1
!    			End loop when new column has the same number of elements as the
!				previous column.
          if( .not. associated(p2%up) ) exit
       end do
    else
       if(verbose) write(*,*) 'Head pointer unassociated. Creating first column.'
       allocate(head,STAT=istat)
       if(istat > 0)then
          write(*,*) 'Error in create_column -- allocating head'
          stop
       end if
       p1=>head ; p2=> null()
       do i = 2 , neta
          allocate(p1%up,STAT=istat)
          if(istat > 0)then
             write(*,*) 'Error in create_column -- allocating initial column'
             stop
          end if
          p1%up%down => p1
          p1 =>p1%up
       end do
    end if
!    dimensions(1) = dimensions(1) + 1
    call set_dimensions([dimensions(1)+1,dimensions(2)])
    if(verbose) write(*,*) 'nxi = ', dimensions(1)

  end subroutine allocate_column

  subroutine delete_column(verbose_flag)
    use global_data, only: dimensions, set_dimensions
    implicit none
    logical, intent(in), optional :: verbose_flag
    logical :: verbose = .false.
    integer :: n, ierror = 0
    type (node), pointer :: p => null(), p2 => null(), p3 => null()

    if(present(verbose_flag)) verbose = verbose_flag

    if(verbose)write(*,*) 'Number of columns: ',dimensions(1)
    p => head
    if(dimensions(1) > 1)then
       do n = 2, dimensions(1)
          if(associated(p%right))then
             p => p%right
          else
             ierror = 1
          end if
       end do
    end if
    if(associated(p%right) .or. ierror/=0)then
       write(*,*) 'Error in delete_column -- inconsistent xi-dimensions'
       stop
    end if

    do n = 1, dimensions(2)
       p2 => p%up
       p3 => p%left
       if(verbose)write(*,*) 'Association status of current pointer: ', associated(p2)
       if(n/=dimensions(2) .and. .not. associated(p2))then
          ierror = 1
          write(*,*) n
       end if
       deallocate(p)
       if(associated(p3))p3%right => null()
       p => p2
    end do
    if(associated(p) .or. ierror/=0)then
       write(*,*) 'Error in delete_column -- inconsistent eta-dimensions'
       stop
    end if
!    dimensions(1) = dimensions(1) - 1
    call set_dimensions([dimensions(1)-1,dimensions(2)])
    if(dimensions(1)==0)then
       head => null()
       write(*,*) 'Deallocating head pointer'
    end if
  end subroutine delete_column

  function get_node(indices, direction, direction2)
    use boundary_conditions, only: boundary_node
    implicit none
    type (node_data) :: get_node, p2
    integer, dimension(2), intent(in) :: indices
    character(len=*), intent(in), optional :: direction, direction2
    integer :: i, j, n
    type (node), pointer :: p => null()

    i = indices(1) ; j = indices(2)
    p => head
    if(i>1)then
       do n = 2, i
          p => p%right
       end do
    end if
    if(j>1)then
       do n = 2, j
          p => p%up
       end do
    end if
    p2 = p
    if(present(direction))then
       select case(direction)
       case('left')
          p => p%left
       case('right')
          p => p%right
       case('up')
          p => p%up
       case('down')
          p => p%down
       case default
          write(*,*) 'Error in get_node -- illegal direction flag'
          stop
       end select
    end if
    if(present(direction2))then
       if(associated(p))then
          select case(direction2)
          case('left')
             p => p%left
          case('right')
             p => p%right
          case('up')
             p => p%up
          case('down')
             p => p%down
          case default
             write(*,*) 'Error in get_node -- illegal direction flag'
             stop
          end select
       else
          write(*,*) 'Error in get_node -- attempt to step away from unassociated node'

       end if
    end if

    if(associated(p))then
       get_node = p
    else
       get_node = boundary_node(p2, direction)
       get_node%boundary = .true.
!			write(*,*) 'boundary_node'
    end if
!		write(*,*) get_node
    
  end function get_node

  subroutine set_node(indices, node_value)
    implicit none
    type (node_data), intent(in) :: node_value
    integer, dimension(2), intent(in) :: indices
    integer :: i, j, n
    type (node), pointer :: p

    i = indices(1) ; j = indices(2)
    p => head
    if(i>1)then
       do n = 2, i
          p => p%right
       end do
    end if
    if(j>1)then
       do n = 2, j
          p => p%up
       end do
    end if

    if(associated(p))then
       p=p.copy.node_value
    else
       write(*,*) 'Error in set_node -- Cannot set out-of-bounds node'
       stop
    end if

  end subroutine set_node

	!__****f* boundary_conditions/initialize_column
!__ NAME
!__ 	Initialize_Column -- Apply Dirichlet boundary conditions to an uninitialized
!__ 	column of cells that has previously been allocated and linked.
!__ INPUTS
!__ 	head -- Pointer to the node_data element at the bottom-left (nxi, neta = 1,1)
!__ 		corner of the simulation region.
!__ 	h0 -- Nominal value of grid-motion-function h. Possible redundant if using g
!__ 	differentials -- Two-element array (dxi, deta) specifying dimensions of grid
!__ 		cells
!__ 	dimensions -- Two-element array (nxi, neta) specifying dimensions of grid
!__ 	top_, bottom_, left_, and right_boundary -- Boundary arrays (read-only)
!__ 	leftcoords -- Upstream grid specification (read-only)
!__ OUTPUTS
!__ 	Initialized column of cells starting from head, obeying supersonic upstream
!__ 	boundary condition
!__ DESCRIPTION
!__ 	Initialize_column applies Dirichlet (supersonic inflow) boundary conditions
!__ 	to an already allocated and linked column of cells based on the information
!__ 	obtained in boundary_conditions_init. Each cell is assigned coordinates
!__ 	from the array leftcoords, and those coordinates are then used to determine
!__ 	the flow state (density, pressure, velocity) in that cell. The metric is then
!__ 	either computed based on neighboring coordinate values or else set equal to
!__ 	the Cartesian metric. Finally, g is either extrapolated from within the
!__ 	simulation region or computed from h0.
!__****
  subroutine initialize_column(verbose_flag,xi_flag,xi,eta_in)
    use global_data, only: h0, differentials, dimensions
    use types, only: node_data, node
    use boundary_conditions_init, only: top_boundary, bottom_boundary, left_boundary, right_boundary, leftcoords
    use boundary_conditions, only: boundary_node
    implicit none
!		type (node), pointer, intent(in) :: head ! Pointer to bottom-left-most node
    logical, intent(in), optional :: verbose_flag ! Activates many intermediate write statements
    logical, intent(in) :: xi_flag
    real(8), intent(in), optional :: xi
    real(8), intent(in), optional, dimension(:) :: eta_in
    type (node_data) :: temp, p
    logical :: verbose = .false.
    integer :: neta, n
!    real(8) :: theta
    neta = dimensions(2)

    if( present(verbose_flag) )verbose = verbose_flag

    do n = 1 , neta
       if(verbose)write(*,*) 'n = ', n
       if(xi_flag)then
         p%x = xi
       else
         p%x = leftcoords(n)%x
       end if
       if(present(eta_in))then
         p%y = eta_in(n)
       else
         p%y = leftcoords(n)%y
       end if
       temp = boundary_node(p,'left')
       p%rho= temp%rho
       p%p  = temp%p
       p%u  = temp%u
       p%v  = temp%v
!       if(dimensions(1)>1 .and. .not. xi_flag)then
!          temp = get_node((/1,n/),'right')
!          p%g  =  temp%g
!          p%A  = (temp%x - p%x)/differentials(1)
!          p%L  = temp%L - differentials(1)/differentials(2)*(temp%A-p%A)
!          p%B  = (temp%y - p%y)/differentials(1)
!          p%M  = temp%M - differentials(1)/differentials(2)*(temp%B-p%B)
!       else
          p%g  =  log(h0*sqrt(p%u**2+p%v**2))
          p%A  =  1.0
          p%B  =  0.0
          p%L  = 0.0
          p%M  = 1.0
!       end if
       p%q  = sqrt(p%u**2+p%v**2)
       p%theta = atan2(p%v,p%u)
       p%h  = h0
       p%g  = log(p%h*p%q)
       if(verbose)write(*,*) 'y = ', p%y,'p = ', p%rho,p%u,p%p
       call set_node( (/1,n/) , p )
    end do
  end subroutine initialize_column

  subroutine node_node_data_assign(out,in)
    implicit none
    type (node_data), intent(out) :: out
    type (node), intent(in) :: in

    out%rho = in%rho
    out%p   = in%p
    out%u   = in%u
    out%v   = in%v
    out%A   = in%A
    out%B   = in%b
    out%L   = in%L
    out%M   = in%M
    out%x   = in%x
    out%y   = in%y
    out%h   = in%h
    out%g   = in%g
    out%q   = in%q
    out%theta = in%theta
  end subroutine node_node_data_assign

  function ndnc(n,nd)
    implicit none
    type (node), intent(in) :: n
    type (node_data), intent(in) :: nd
    type (node) :: ndnc

    ndnc%rho   = nd%rho
    ndnc%p     = nd%p
    ndnc%u     = nd%u
    ndnc%v     = nd%v
    ndnc%A     = nd%A
    ndnc%B     = nd%b
    ndnc%L     = nd%L
    ndnc%M     = nd%M
    ndnc%x     = nd%x
    ndnc%y     = nd%y
    ndnc%h     = nd%h
    ndnc%g     = nd%g
    ndnc%q     = nd%q
    ndnc%theta = nd%theta

    ndnc%up    => n%up
    ndnc%down  => n%down
    ndnc%left  => n%left
    ndnc%right => n%right
  end function ndnc
!! Broken -- overwrites initial pointers and nullifies them.
!	subroutine node_data_node_assign(out,in)
!		implicit none
!		type (node), intent(out) :: out
!		type (node_data), intent(in) :: in
!
!		out%rho = in%rho
!		out%p   = in%p
!		out%u   = in%u
!		out%v   = in%v
!		out%A   = in%A
!		out%B   = in%b
!		out%L   = in%L
!		out%M   = in%M
!		out%h   = in%h
!		out%g   = in%g
!		out%q   = in%q
!		out%theta = in%theta
!	end subroutine node_data_node_assign



  integer function create_column_tester(verbose_flag)
!__****f* create_column/create_column_tester
!__ NAME
!__ 	Create_column_tester
!__ DESCRIPTION
!__ 	Test routine for node_array/create_column. Allocates three ten-element
!__ 	columns, tests the association status of their pointers, and tests that
!__ 	each "inverse" pointer in fact links back to the original point. Testing
!__ 	has NOT been implemented for determining that each node is actually a
!__ 	distinct memory location--that is, that pointer address are distinct.
    use global_data, only: dimensions, set_dimensions
    implicit none
    logical, intent(in), optional :: verbose_flag
    integer :: neta = 10, n = 0
    type (node), pointer :: p
    logical, dimension(4) :: pattern = .false.
    logical, allocatable, dimension(:) :: success
    logical :: verbose = .false.
    integer :: out=0

    if(present(verbose_flag))verbose = verbose_flag
!    dimensions(2) = 10
    call set_dimensions([0,10])

! Test creation of initial array column
    call allocate_column(verbose_flag=.false., neta_in=neta)
    p=>head
!		write(*,*) 'Head association state: '
!		write(*,*) '    left -- ', temp(1)
!		write(*,*) '   right -- ', temp(2)
!		write(*,*) '      up -- ', temp(3)
!		write(*,*) '    down -- ', temp(4)
    allocate( success(neta) )
    pattern = (/ .false., .false., .true., .false. /)
    success(1) = compare_state(pattern=pattern, state=node_association_state(p))
    if(verbose)write(*,*) 'Head allocated successfully? Success = ', success(1)

    do n = 2, neta - 1
       p => p%up
!			write(*,*) 'Association state, node ',n
!			write(*,*) '    left -- ', associated(p%left)
!			write(*,*) '   right -- ', associated(p%right)
!			write(*,*) '      up -- ', associated(p%up)
!			write(*,*) '    down -- ', associated(p%down)
       pattern = (/ .false., .false., .true., .true. /)
       success(n) = compare_state(pattern=pattern, state=node_association_state(p))
       if(verbose)&
            write(*,*) 'Node ', n, 'allocated successfully? Success = ', success(n)
    end do
    p=> p%up
    pattern = (/ .false., .false., .false., .true. /)
    success(neta) = compare_state(pattern=pattern, state=node_association_state(p))
    if(verbose)then
       write(*,*) 'Node ', neta, 'allocated successfully? Success = ', success(neta)
       write(*,*) 'Column allocated successfully? Success = ', reduce_logical_array(success)
    end if
    if(reduce_logical_array(success))then
       out = out + 0
    else
       out = 1
    end if
! Test creation of additional column
    call allocate_column(verbose_flag=.false., neta_in=neta)
    p=>head
    pattern = (/ .false., .true., .true., .false. /)
    success(1) = compare_state(pattern=pattern, state=node_association_state(p))
    if(verbose)write(*,*) 'New head allocated successfully? Success = ', success(1)
    do n = 2, neta - 1
       p => p%up
!			write(*,*) 'Association state, node ',n
!			write(*,*) '    left -- ', associated(p%left)
!			write(*,*) '   right -- ', associated(p%right)
!			write(*,*) '      up -- ', associated(p%up)
!			write(*,*) '    down -- ', associated(p%down)
       pattern = (/ .false., .true., .true., .true. /)
       success(n) = compare_state(pattern=pattern, state=node_association_state(p))
       if(verbose)&
            write(*,*) 'Node ', n, 'allocated successfully? Success = ', success(n)
    end do
    p=> p%up
    pattern = (/ .false., .true., .false., .true. /)
    success(neta) = compare_state(pattern=pattern, state=node_association_state(p))
    if(verbose)then
       write(*,*) 'Node ', neta, 'allocated successfully? Success = ', success(neta)
       write(*,*) 'Column allocated successfully? Success = ', reduce_logical_array(success)
    end if
    if(reduce_logical_array(success))then
       out = out + 0
    else
       out = 1
    end if
! Test links in initial column
    p=>head%right
    pattern = (/ .true., .false., .true., .false. /)
    success(1) = compare_state(pattern=pattern, state=node_association_state(p))
    if(verbose)write(*,*) 'Old head linked successfully? Success = ', success(1)
    do n = 2, neta - 1
       p => p%up
!			write(*,*) 'Association state, node ',n
!			write(*,*) '    left -- ', associated(p%left)
!			write(*,*) '   right -- ', associated(p%right)
!			write(*,*) '      up -- ', associated(p%up)
!			write(*,*) '    down -- ', associated(p%down)
       pattern = (/ .true., .false., .true., .true. /)
       success(n) = compare_state(pattern=pattern, state=node_association_state(p))
       if(verbose)&
            write(*,*) 'Node ', n, 'linked successfully? Success = ', success(n)
    end do
    p=> p%up
    pattern = (/ .true., .false., .false., .true. /)
    success(neta) = compare_state(pattern=pattern, state=node_association_state(p))
    if(verbose)then
       write(*,*) 'Node ', neta, 'linked successfully? Success = ', success(neta)
       write(*,*) 'Column linked successfully? Success = ', reduce_logical_array(success)
    end if
    if(reduce_logical_array(success))then
       out = out + 0
    else
       out = 1
    end if
! Test links of internal column
    call allocate_column(verbose_flag=.false., neta_in=neta)
    p=>head%right
    pattern = (/ .true., .true., .true., .false. /)
    success(1) = compare_state(pattern=pattern, state=node_association_state(p))
    if(verbose)write(*,*) 'Old head linked successfully? Success = ', success(1)
    do n = 2, neta - 1
       p => p%up
!			write(*,*) 'Association state, node ',n
!			write(*,*) '    left -- ', associated(p%left)
!			write(*,*) '   right -- ', associated(p%right)
!			write(*,*) '      up -- ', associated(p%up)
!			write(*,*) '    down -- ', associated(p%down)
       pattern = (/ .true., .true., .true., .true. /)
       success(n) = compare_state(pattern=pattern, state=node_association_state(p))
       if(verbose)&
            write(*,*) 'Node ', n, 'linked successfully? Success = ', success(n)
    end do
    p=> p%up
    pattern = (/ .true., .true., .false., .true. /)
    success(neta) = compare_state(pattern=pattern, state=node_association_state(p))
    if(verbose)then
       write(*,*) 'Node ', neta, 'linked successfully? Success = ', success(neta)
       write(*,*) 'Column linked successfully? Success = ', reduce_logical_array(success)
    end if
    if(reduce_logical_array(success))then
       out = out + 0
    else
       out = 1
    end if

    deallocate(success)

    create_column_tester = out

  end function create_column_tester

  logical function reduce_logical_array(array)
!__****f* create_column_tester/reduce_logical_array
!__ NAME
!__ 	Reduce_logical_array
!__ DESCRIPTION
!__ 	Support function for create_column_tester. Takes a logical array and
!__ 	returns a logical .true. if all element of the array are true, and
!__ 	.false. otherwise.
    implicit none
    logical, intent(in), dimension(:) :: array
    integer :: n = 0
    reduce_logical_array = .true.
    do n = 2, size(array)
       reduce_logical_array = array(n) .and. array(n-1) .and. reduce_logical_array
    end do
  end function reduce_logical_array

  function node_association_state(p)
!__****f* create_column_tester/node_association_state
!__ NAME
!__ 	Node_association_state
!__ DESCRIPTION
!__ 	Support function for create_column_tester. Accepts node pointer and
!__ 	evaluates the association state of each of its link pointers. Also
!__ 	checks that the "inverse" pointer in fact returns the original node.
!__ 	Outputs a four-element logical array containing pass/fail information
!__ 	for each link.
    implicit none
    type (node), intent(in), pointer :: p
    logical, dimension(4) :: node_association_state, temp
    temp = (/ associated(p%left), associated(p%right), &
         associated(p%up), associated(p%down) /)
    !		Test that each associated link goes both ways.
    if(temp(1))temp(1)=temp(1).and.associated(p%left%right,p)
    if(temp(2))temp(2)=temp(2).and.associated(p%right%left,p)
    if(temp(3))temp(3)=temp(3).and.associated(p%up%down,p)
    if(temp(4))temp(4)=temp(4).and.associated(p%down%up,p)
    node_association_state = temp
  end function node_association_state

  logical function compare_state(pattern,state)
!__****f* create_column_tester/compare_state
!__ NAME
!__ 	Compare_state
!__ DESCRIPTION
!__ 	Support utility for Create_column_tester. Compares two logical arrays
!__ 	and returns a scalar T if they are identical.
    implicit none
    logical, intent(in), dimension(:) :: pattern, state
    integer :: n
    if(size(pattern) /= size(state))then
       write(*,*) 'Error in compare_state -- incompatible arrays'
       stop
    end if
    compare_state = .true. ! Initialize to .true. to enable .and. logic.
    do n = 1 , size(pattern)
       compare_state = compare_state .and. (pattern(n) .eqv. state(n))
    end do
  end function compare_state


end module node_array

!!$program node_array_tester
!!$  use node_array
!!$  use boundary_conditions_init, only: read_boundary
!!$  use global_data, only: dimensions
!!$  implicit none
!!$  integer :: n
!!$  type (node), pointer :: p => null()
!!$  
!!$  write(*,*) 'Create_column_tester',create_column_tester(verbose_flag=.true.)
!!$  write(*,*) dimensions
!!$  do n = 1, dimensions(1)
!!$     write(*,*) 'Calling delete_column'
!!$     call delete_column()
!!$  end do
!!$  call read_boundary(verbose_flag=.false.)
!!$  write(*,*) 'Finished read_boundary'
!!$  call create_column(verbose_flag=.true.,neta_in=10)
!!$  write(*,*) 'Finished create_column'
!!$  call initialize_column(verbose_flag=.true.)
!!$  write(*,*) 'Finished initialize_column'
!!$  read(*,*)
!!$  write(*,*) get_node((/1,1/))
!!$  write(*,*) 'Finished get_node'
!!$  read(*,*)
!!$end program node_array_tester
