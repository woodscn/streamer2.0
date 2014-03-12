module output_data
    implicit none
    contains
    subroutine write_files(time,first_flag)
      use global_data, only: dimensions, output_file_name
      use types, only: node_data
      use node_array, only: get_node
      implicit none
      save

      real(8), intent(in) :: time
      logical, intent(in), optional :: first_flag
      logical :: first = .false.
      real(8), allocatable, dimension(:,:) :: rho, p, u, v, x, y
      real(8), allocatable, dimension(:,:,:) :: geom
      integer :: i, j, neta, nxi, error
      type (node_data) :: node
      character(len=128) :: error_message
      real :: temp

      temp = time

      nxi = dimensions(1) ; neta = dimensions(2)
      first = .false.
      if(present(first_flag))first = first_flag
      allocate( rho(neta, nxi), p(neta, nxi), u(neta, nxi), v(neta, nxi)&
        , x(neta, nxi), y(neta, nxi), geom(neta,nxi,4), stat = error )
      if(error /= 0)then
        write(*,*) "Arrays not allocated!!"
        stop
      end if

      do i = 1, nxi
        do j = 1, neta
          node = get_node([i,j])
          rho(j,i) = node%rho
          p(j,i) = node%p
          u(j,i) = node%u
          v(j,i) = node%v
          x(j,i) = node%x
          y(j,i) = node%y
          geom(j,i,1) = node%A
          geom(j,i,2) = node%B
          geom(j,i,3) = node%L
          geom(j,i,4) = node%M
        end do
      end do

      if(first)then
        open( unit = 3141, file = output_file_name, iostat = error, iomsg = error_message, action = 'write' )
        write(3141,*) 'variables = "x", "y", "rho", "p", "u", "v", "A", "B", "L", "M"'
        if( error /= 0 )then
          write(*,*) 'Error opening file for output'
          write(*,*) error_message
          stop
        end if
      end if
      1001 format (10E24.15)
      write(3141,*) 'zone f = "block"  i = ', neta, ' j = ', nxi
      write(3141,1001) x
      write(3141,*) ' '
      write(3141,1001) y
      write(3141,*) ' '
      write(3141,1001) rho
      write(3141,*) ' '
      write(3141,1001) p
      write(3141,*) ' '
      write(3141,1001) u
      write(3141,*) ' '
      write(3141,1001) v
      write(3141,*) ' '
      write(3141,1001) geom(:,:,1)
      write(3141,*) ' '
      write(3141,1001) geom(:,:,2)
      write(3141,*) ' '
      write(3141,1001) geom(:,:,3)
      write(3141,*) ' '
      write(3141,1001) geom(:,:,4)
      write(3141,*) ' '
      flush(3141)


      deallocate( rho, p, u, v, x, y, geom, stat = error )
      if(error/=0)then
        write(*,*) "Error in deallocation"
        stop
      end if
    end subroutine write_files

    subroutine write_files_matlab(time,first_flag)!(x,y,E,h,t,dt)
! Inputs: x, y, E (array of conserved variables), t, dt
! Outputs: two_D.dat, an ASCII file that contains the primitive variables together
!          with their coordinate values and the associated time. Data is designed
!          to be read by the matlab file geom_data_reader.m. It should be noted
!          that the file is opened and headers are written in the main program.
    use global_data, only: dimensions, output_file_name
    use types, only: node_data
    use node_array, only: get_node

      implicit none
      real(8), intent(in) :: time
      logical, intent(in), optional :: first_flag
      logical :: first
      integer :: i,j
      character(len=128) :: error_message
      integer :: error
      type (node_data) :: current

      first = .false.
      if(present(first_flag))first = first_flag
      if(first)then
        open( unit = 3141, file = output_file_name, iostat = error, iomsg = error_message, action = 'write' )
        write(3141,*) ' nt= ' , 1000
!        write(3141,*) ' '
        if( error /= 0 )then
          write(*,*) 'Error opening file for output'
          write(*,*) error_message
          stop
        end if
      end if
! We output primitive variables
      write(3141,*) ' '
      write(3141,*) 't=',time,'nx=',dimensions(1),'ny=',dimensions(2),'dt=',0.
      write(3141,*) ' '
      write(3141,*) 'x= ','y= ','u= ','v= ','rho= ','P= '
      do i = 1 , dimensions(1)
         do j = 1 , dimensions(2)
            current = get_node([i,j])
            write(3141,*) ' ', current%x , current%y , current%u , current%v , current%rho , current%p &
                 , current%A , current%B , current%L , current%M , current%h
         end do
      end do
    end subroutine write_files_matlab


end module output_data
