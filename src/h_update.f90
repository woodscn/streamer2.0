module h_update
  use global_data, only: h0
  implicit none

contains
!!$  subroutine h_tester(verbose_flag)
!!$    implicit none
!!$    
!!$  end subroutine h_tester

  elemental function g_test(x,y,type)
    use global_data, only: PI
    implicit none
    integer, intent(in) :: type
    real(8), intent(in) :: x, y
    real(8) :: g_test
    select case(type)
    case(1)
       g_test = sin(x)*sin(y)+1.5
    case(2)
       g_test = x*(x+1)*exp(y)+1.
    case(3)
       g_test = log(x+1)*(y)+1.
    case(4)
       g_test = 1.5-x*y**3
    case(5)
       g_test = log(x+y+1.5)*cos(x)
    case(6)
       g_test = exp(y/(x+1.))+sin((x+1.5)-y)+1.
    case(7)
       g_test = tan(PI*(x-0)*1)*tan(PI*(y-0)*1)
    end select
  end function g_test

  elemental function dg_testdx(x,y,type)
    use global_data, only: PI
    implicit none
    integer, intent(in) :: type
    real(8), intent(in) :: x, y
    real(8) :: dg_testdx
    integer, parameter :: n=1, m=1
    real(8), parameter :: x0=0., y0=0.
    select case(type)
    case(1)
       dg_testdx = cos(x)*sin(y)
    case(2)
       dg_testdx = (2*x+1)*exp(y)
    case(3)
       dg_testdx = y/x
    case(4)
       dg_testdx = -y**3
    case(5)
       dg_testdx = cos(x)/(x+y+1.5)-log(x+y+1.5)*sin(x)
    case(6)
       dg_testdx = y/(x+1.)**2*(y/(x+1.))*exp(y/(x+1.)) + cos((x+1.5)-y)
    case(7)
       dg_testdx = n*PI/cos(n*PI*(x-x0))**2*tan(m*PI*(y-y0))
    end select
  end function dg_testdx

  elemental function dg_testdy(x,y,type)
    use global_data, only: PI
    implicit none
    integer, intent(in) :: type
    real(8), intent(in) :: x, y
    real(8) :: dg_testdy
    integer, parameter :: n=1, m=1
    real(8), parameter :: x0=0., y0=0.
    select case(type)
    case(1)
       dg_testdy =sin(x)*cos(y)
    case(2)
       dg_testdy = x*(x+1)*exp(y)
    case(3)
       dg_testdy = log(x+1)
    case(4)
       dg_testdy = -3.*x*y**2
    case(5)
       dg_testdy = 1./(x+y+1.5)*cos(x)
    case(6)
       dg_testdy = exp(y/(x+1.))/(x+1.)-cos((x+1.5)-y)
    case(7)
       dg_testdy = m*PI*tan(n*PI*(x-x0))/cos(m*PI*(y-y0))**2
    end select
  end function dg_testdy

    
  subroutine iteration(verbose_flag,verify_flag)
    use node_array, only: get_node, set_node
    use types, only: node_data
    use global_data, only: differentials, dimensions, h0
    implicit none
    integer :: i, j, nxi, neta
    logical, intent(in), optional :: verbose_flag, verify_flag
    logical :: verbose = .false., verify = .false.
    real(8), allocatable, dimension(:,:) :: h, g, theta, q, A, B, L, M, S2, T2, &
         co, si, alpha, beta, gamma, change, dtdxi, dtdeta, xi, eta
    !      volatile :: g

    real(8) :: tol = 1.d-14
    integer :: iter, inda, indb
    integer :: alphatest = 1, betatest = 1, gammatest = 1, gtest = 1
    type(node_data) :: p

    if(present(verbose_flag)) verbose = verbose_flag
    if(present(verify_flag))  verify  = verify_flag
    nxi = dimensions(1) ; neta = dimensions(2)
    !      if(dimensions(1) > 1)then

    allocate( h(nxi,neta), g(nxi,neta), theta(nxi,neta), q(nxi,neta), A(nxi,neta),&
         B(nxi,neta), L(nxi,neta), M(nxi,neta), S2(nxi,neta), T2(nxi,neta), co(nxi,neta),&
         si(nxi,neta), alpha(nxi,neta), beta(nxi,neta), gamma(nxi,neta), dtdxi(nxi,neta),&
         dtdeta(nxi,neta), change(nxi,neta) )

    do i = 1, nxi
       do j = 1, neta
          p = get_node([i,j])
          h(i,j) = p%h
          theta(i,j) = p%theta
          q(i,j) = p%q
          g(i,j) = p%g
          A(i,j) = p%A
          B(i,j) = p%B
          L(i,j) = p%L
          M(i,j) = p%M
       end do
    end do

    !      dtdxi(1,:)       = (theta(2,:) - theta(1,:))/differentials(1)
    !      dtdxi(nxi,:)     = (theta(nxi,:)-theta(nxi-1,:))/differentials(1)
    !      if(nxi>2)then
    !        dtdxi(2:nxi-1,:) = (theta(2:nxi-1,:) - theta(1:nxi-2,:))/(2.*differentials(1))
    !      end if
    !
    !      dtdeta(:,1)        = (theta(:,2) - theta(:,1))/differentials(2)
    !      dtdeta(:,neta)     = (theta(:,neta) - theta(:,neta-1))/differentials(2)
    !      if(neta>2)then
    !        dtdeta(:,2:neta-1) = (theta(:,2:neta-1) - theta(:,1:neta-2))/(2.*differentials(2))
    !      end if

    call array_gradient(theta,dtdxi,dtdeta)

    S2 = L**2 + M**2
    T2 = A**2 + B**2

    alpha = S2*(A*sin(theta) - B*cos(theta))
    beta  = T2*(M*cos(theta) - L*sin(theta))
    gamma = S2*(A*cos(theta) + B*sin(theta))*dtdxi &
         - T2*(L*cos(theta) + M*sin(theta))*dtdeta

    if(verify)then
       allocate(xi(dimensions(1),dimensions(2)),eta(dimensions(1),dimensions(2)))
       do inda = 1, dimensions(1)
          do indb = 1, dimensions(2)
             xi (inda,indb) = differentials(1)*(inda-1)
             eta(inda,indb) = differentials(2)*(indb-1)
          end do
       end do
       alphatest = 1 ; betatest = 1 ; gtest = 1
       alpha = g_test(xi,eta,alphatest)
       beta  = g_test(xi,eta, betatest)
       gamma = -(alpha*dg_testdx(xi,eta,gtest)+beta*dg_testdy(xi,eta,gtest))
    end if



    change = 0.
    iter = 0
    if(verify)then
       g=1.5
    else
       g=0.
    end if

    do
       iter = iter + 1

       g     (2:nxi,2:neta) =&
             ( alpha(2:nxi,2:neta)*g(1:nxi-1,2:neta)*differentials(1) &
            +   beta(2:nxi,2:neta)*g(2:nxi,1:neta-1)*differentials(2)  &
            -  gamma(2:nxi,2:neta)*differentials(1)*differentials(2) ) &
            / (alpha(2:nxi,2:neta)*differentials(2) +  beta(2:nxi,2:neta)*differentials(1) )

       change(2:nxi,2:neta) =&
             ( alpha(2:nxi,2:neta)*g(1:nxi-1,2:neta)*differentials(1) &
            +   beta(2:nxi,2:neta)*g(2:nxi,1:neta-1)*differentials(2)  &
            -  gamma(2:nxi,2:neta)*differentials(1)*differentials(2) ) &
            / (alpha(2:nxi,2:neta)*differentials(2) +  beta(2:nxi,2:neta)*differentials(1) )&
            - g(2:nxi,2:neta)
       
       if(maxval(abs(change(2:nxi,2:neta))) < tol)then
          if(verbose) write(*,*) 'iter = ', iter
          exit
       end if

       if(iter > 1000 )then
          write(*,*) 'did not converge'
          stop
       end if
    end do
    h = exp(g)/q
    h=h/maxval(h)*h0

    if(verify)then
       write(*,*) "Verification parameter = ", sqrt(1./(dimensions(1)*dimensions(2))*sum(&
            (g_test(xi,eta,gtest)-g)/g&
            **2))
       deallocate(xi,eta)
       read(*,*)
    end if

    do i = 1, nxi
       do j = 1, neta
          p = get_node([i,j])
          p%h = h(i,j)
          p%g = g(i,j)
          call set_node([i,j],p)
       end do
    end do
    
    deallocate( h, g, theta, q, A,&
         B, L, M, S2, T2, co,&
         si, alpha, beta, gamma, dtdxi,&
         dtdeta, change )

    !      end if
  end subroutine iteration

  subroutine variational_iteration(verbose_flag,verify_flag)
    use node_array, only: get_node, set_node
    use types, only: node_data
    use global_data, only: differentials, dimensions, h0
    implicit none
    integer :: i, j, nxi, neta
    integer :: alphatest= 1, betatest= 1, gammatest= 1, gtest = 2
    logical, intent(in), optional :: verbose_flag, verify_flag
    logical :: verbose = .false.
    real(8), allocatable, dimension(:,:) :: h, g, theta, q, A, B, L, M, S2, T2, &
         co, si, alpha, beta, gamma, change, dtdxi, dtdeta, dgdxi, dgdeta, temp, xi, eta
    !      volatile :: g

    real(8) :: tol = 1.d-3, tempmax
    integer :: iter
    type(node_data) :: p

    if(present(verbose_flag)) verbose = verbose_flag
    nxi = dimensions(1) ; neta = dimensions(2)
    !      if(dimensions(1) > 1)then

    allocate( h(nxi,neta), g(nxi,neta), theta(nxi,neta), q(nxi,neta), A(nxi,neta),&
         B(nxi,neta), L(nxi,neta), M(nxi,neta), S2(nxi,neta), T2(nxi,neta), co(nxi,neta),&
         si(nxi,neta), alpha(nxi,neta), beta(nxi,neta), gamma(nxi,neta), dtdxi(nxi,neta),&
         dtdeta(nxi,neta), change(nxi,neta), dgdxi(nxi,neta), dgdeta(nxi,neta), temp(1,neta) )

    do i = 1, nxi
       do j = 1, neta
          p = get_node([i,j])
          h(i,j) = p%h
          theta(i,j) = p%theta
          q(i,j) = p%q
          g(i,j) = p%g
          A(i,j) = p%A
          B(i,j) = p%B
          L(i,j) = p%L
          M(i,j) = p%M
       end do
    end do
    call array_gradient(theta,dtdxi,dtdeta)

    S2 = L**2 + M**2
    T2 = A**2 + B**2

    alpha = S2*(A*sin(theta) - B*cos(theta))!/differentials(1)
    beta  = T2*(M*cos(theta) - L*sin(theta))!/differentials(2)
    gamma = S2*(A*cos(theta) + B*sin(theta))*dtdxi &
         - T2*(L*cos(theta) + M*sin(theta))*dtdeta

    if(present(verify_flag))then
       if(verify_flag)then
          allocate(xi(dimensions(1),dimensions(2)),eta(dimensions(1),dimensions(2)))
          xi = reshape( (/ ((differentials(1)*(i-1)/real(dimensions(1))&
               ,i=1,dimensions(1)),j=1,dimensions(2))/),(/dimensions(1),dimensions(2)/))
          eta= reshape( (/ ((differentials(2)*(j-1)/real(dimensions(2))&
               ,j=1,dimensions(2)),i=1,dimensions(1))/),(/dimensions(1),dimensions(2)/))


          alpha = g_test(&
               reshape( [ ( (differentials(1)*(i-1)/real(dimensions(1)),i=1,dimensions(1))&
               , j=1,dimensions(2) ) ] ,[dimensions(1),dimensions(2)]),&
               reshape( [ ( (differentials(2)*(j-1)/real(dimensions(2)),j=1,dimensions(2))&
               , i=1,dimensions(1) ) ] ,[dimensions(1),dimensions(2)]),&
               alphatest)
          beta  = g_test(&
               reshape( [ ( (differentials(1)*(i-1)/real(dimensions(1)),i=1,dimensions(1))&
               , j=1,dimensions(2) ) ] ,[dimensions(1),dimensions(2)]),&
               reshape( [ ( (differentials(2)*(j-1)/real(dimensions(2)),j=1,dimensions(2))&
               , i=1,dimensions(1) ) ] ,[dimensions(1),dimensions(2)]),&
               betatest)
          gamma = g_test(&
               reshape( [ ( (differentials(1)*(i-1)/real(dimensions(1)),i=1,dimensions(1))&
               , j=1,dimensions(2) ) ] ,[dimensions(1),dimensions(2)]),&
               reshape( [ ( (differentials(2)*(j-1)/real(dimensions(2)),j=1,dimensions(2))&
               , i=1,dimensions(1) ) ] ,[dimensions(1),dimensions(2)]),&
               gammatest)
          
          
       end if
    end if

    change = 0.
    iter = 0
    temp = 0.
    do
       iter = iter + 1
       call array_gradient(g,dgdxi,dgdeta)
       tempmax = 0.
       do i = 2, nxi
          temp(1,:) = (dgdxi(i,:) + beta(i,:)/alpha(i,:)*dgdeta(i,:) + gamma(i,:)/alpha(i,:) &
               + temp(1,:))*differentials(1)
          g(i,:) = g(i,:) - temp(1,:)
          tempmax = max(tempmax,abs(maxval(temp)))
       end do

       if(tempmax < tol)then
          if(verbose) write(*,*) 'iter = ', iter
          exit
       end if

       if(iter > 1000 )then
          write(*,*) 'did not converge'
          stop
       end if
    end do
    h = exp(g)/q
    h=h/maxval(h)*h0

!!$    if(present(verify_flag))then
!!$       if(verify_flag)then
!!$          write(*,*) "Verification parameter = ", sqrt(1./(dimensions(1)*dimensions(2))*sum(sum(&
!!$               (alpha*dg_testdx(&
!!$               reshape( [ ( (differentials(1)*(i-1)/real(dimensions(1)),i=1,dimensions(1))&
!!$               , j=1,dimensions(2) ) ] ,[dimensions(1),dimensions(2)]),&
!!$               reshape( [ ( (differentials(2)*(j-1)/real(dimensions(2)),j=1,dimensions(2))&
!!$               , i=1,dimensions(1) ) ] ,[dimensions(1),dimensions(2)]),&
!!$               gtest)+beta*dg_testdy(&
!!$               reshape( [ ( (differentials(1)*(i-1)/real(dimensions(1)),i=1,dimensions(1))&
!!$               , j=1,dimensions(2) ) ] ,[dimensions(1),dimensions(2)]),&
!!$               reshape( [ ( (differentials(2)*(j-1)/real(dimensions(2)),j=1,dimensions(2))&
!!$               , i=1,dimensions(1) ) ] ,[dimensions(1),dimensions(2)]),&
!!$               gtest)+gamma - g)&
!!$               **2)))
!!$       end if
!!$    end if


    do i = 1, nxi
       do j = 1, neta
          p = get_node([i,j])
          p%h = h(i,j)
          p%g = g(i,j)
          call set_node([i,j],p)
       end do
    end do

    deallocate( h, g, theta, q, A,&
         B, L, M, S2, T2, co,&
         si, alpha, beta, gamma, dtdxi,&
         dtdeta, change, dgdxi, dgdeta )

    !      end if
  end subroutine variational_iteration

  subroutine array_gradient(in,dx,dy)
    use global_data, only: dimensions, differentials
    implicit none
    real(8), dimension(:,:), intent(in) :: in
    real(8), dimension(size(in,1),size(in,2)), intent(out) :: dx, dy
    integer :: nxi, neta

    nxi  = dimensions(1)
    neta = dimensions(2)

    dx(1,:)   = (in(2,:) - in(1,:))/differentials(1)
    dx(nxi,:) = (in(nxi,:)-in(nxi-1,:))/differentials(1)
    if(nxi>2)then
       dx(2:nxi-1,:) = (in(2:nxi-1,:) - in(1:nxi-2,:))/(2.*differentials(1))
    end if

    dy(:,1)    = (in(:,2) - in(:,1))/differentials(2)
    dy(:,neta) = (in(:,neta) - in(:,neta-1))/differentials(2)
    if(neta>2)then
       dy(:,2:neta-1) = (in(:,2:neta-1) - in(:,1:neta-2))/(2.*differentials(2))
    end if

  end subroutine array_gradient

end module h_update
