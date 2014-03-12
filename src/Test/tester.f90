program tester
  use types
!!$  use node_array
  implicit none
  integer success
  success = types_tester()
  write(*,*) "Success = ",success
!!$  write(*,*) "Create_node_tester", create_node_tester(.true.)
end program tester
