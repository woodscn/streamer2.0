module pressure_implicit
  use types, only: node_data
  use global_data, only: gamma_const
  implicit none
  type (node_data), save :: state_implicit
  real(8) :: x_implicit, h_implicit, plusminus

contains

  real(8) function Delta_func(M0,psi,h,plusminus)
    implicit none
    real(8), intent(in) :: M0, psi, h, plusminus
    Delta_func = (((gamma_const-2.d0*h+1.d0)-(1.d0-h)&
         *(2.d0-plusminus*M0*(gamma_const-1.d0)))/(psi**((gamma_const-1.d0)&
         /(2.d0*gamma_const))*(gamma_const-2.d0*h+1.d0)&
         -(1.d0-h)*(2.d0-plusminus*M0*(gamma_const-1.d0)))&
         )**(2.d0*h/(gamma_const-2.d0*h+1.d0))
  end function Delta_func

  real function pressure_func(p)
    use types, only: node_data
    implicit none
    real, intent(in) :: p
    real(8) :: a0, M, Delta, Delta0, S, psi

    a0 = sqrt(gamma_const*state_implicit%p/state_implicit%rho)
    M = state_implicit%u/a0
    Delta0 = state_implicit%A*state_implicit%M - state_implicit%B*state_implicit%L
    S = sqrt(state_implicit%L**2 + state_implicit%M**2)
    psi = p/state_implicit%p
    Delta = Delta0*Delta_func(M,psi,h_implicit,plusminus)
    pressure_func = (2.d0*(1.d0-h_implicit)/(gamma_const-2.d0*h_implicit+1.d0)&
         -plusminus*(gamma_const-1.d0)/(gamma_const-2.d0*h_implicit+1.d0)&
         *((1.d0-h_implicit)*M-Delta/S*x_implicit/a0)&
         )**(2.d0*gamma_const/(gamma_const-1.d0)) - psi
  end function pressure_func
end module pressure_implicit

module riemann
  implicit none
contains
  function normal_vector(left,right,direction)
    use types, only: node_data

    implicit none
    type(node_data), intent(in) :: left, right
    character(len=*), intent(in) :: direction

    real(8) :: deltaL, deltaR
    real(8), dimension(2) :: gradL, gradR, normal, normal_vector

    deltaL = ( left%A* left%M -  left%B* left%L)
    deltaR = (right%A*right%M - right%B*right%L)
    select case(direction)
    case('xi')
       gradL = [ left%M , - left%L]/deltaL
       gradR = [right%M , -right%L]/deltaR
    case('eta')
       gradL = [- left%B ,  left%A]/deltaL
       gradR = [-right%B , right%A]/deltaR
    case default
       write(*,*) 'Error in normal_vector -- Bad direction flag.'
       stop
    end select
    normal = gradL + gradR
    normal_vector = normal/sqrt(normal(1)**2 + normal(2)**2)
  end function normal_vector

  function riemann_solve( tempL, tempR, verbose_flag, t_out )
    use types, only: node_data
    use global_data, only: gamma_const
    implicit none
    type (node_data), intent(in) :: tempL, tempR
    type (node_data) :: riemann_solve, left, right
    logical, intent(in), optional :: verbose_flag
    real(8), intent(in), optional :: t_out
    real(8) :: tout = 1.d0
    logical :: verbose = .false.
    integer :: n
    integer, parameter :: nx = 1000

    real(8) :: DL, PL, UL, VL, AL
    real(8) :: DR, PR, UR, VR, AR
    real(8) :: Pstar, Ustar, DstarL, DstarR
    real(8) :: tol = 1.d-14, x
    real(8) :: PsiL, PsiR, temp=1.d0, fL, fR, dfL, dfR
    real(8), dimension(nx,5) :: data

!    if(present(verbose_flag))verbose=verbose_flag
    if(present(t_out)) tout = t_out
    left = tempL ; right = tempR
    DL =  left%rho ; PL =  left%p ; UL =  left%u ; VL =  left%v
    DR = right%rho ; PR = right%p ; UR = right%u ; VR = right%v

    AL = sqrt(gamma_const*PL/DL) ; AR = sqrt(gamma_const*PR/DR)
    Pstar = guessp(left,right)
    if(verbose)write(*,*) "Initial guess P = " , Pstar
    temp = 1.d0 ; n = 0
    do while(abs(temp) .gt. tol)
       if(verbose)write(*,*) n , "Pstar = " , Pstar; n = n + 1
       if(n .gt. 10)then ;  write(*,*) "Failed Convergence" ; stop ; end if

       PsiL = Pstar/PL
       call u_fun( left,Pstar,fL,dfL)

       PsiR = Pstar/PR
       call u_fun(right,Pstar,fR,dfR)

       temp = ( UR - UL + fR + fL )/( dfL + dfR )
       Pstar = max( Pstar - temp , tol )
       if(verbose)write(*,*) "fL , fR , Pstar = " , fL , fR , Pstar
    end do

    Ustar = .5*(UR+fR+UL-fL)
    DstarL = beta(Pstar/PL)*DL
    DstarR = beta(Pstar/PR)*DR
    if(verbose)write(*,*) Ustar , DstarL , DstarR
    if(verbose)then
       open(unit = 1, file = 'riemann_test.dat')
       do n = 1 , nx
          x = (1.d0/(real(nx-1,8))*real(n-1,8)-0.5d0)
          call sample(x/tout,left,right,Pstar,Ustar,DstarL,DstarR,riemann_solve,verbose)
          data(n,1)   = x
          data(n,2:5) = [riemann_solve%rho,riemann_solve%u,riemann_solve%v,riemann_solve%p]
          write(1,*) data(n,:)
       end do
    else
       x = 0.0d0
       call sample(x,left,right,Pstar,Ustar,DstarL,DstarR,riemann_solve,verbose)
    end if

  end function riemann_solve


  real(8) pure function guessp(left,right)
    use global_data, only: gamma_const
    use types, only: node_data
    implicit none
    type (node_data), intent(in) :: left, right
    real(8) :: aL, aR, gL, gR, tol

    aL = sqrt(gamma_const* left%p/ left%rho)
    aR = sqrt(gamma_const*right%p/right%rho)
    tol = 1d-8

    ! Linearised guess
    guessp = .5*(left%p+right%p)&
         -.125*(right%u-left%u)*(left%rho+right%rho)*(aL+aR)
    if(.not.( guessp .gt. min(left%p,right%p) .and. guessp .lt. max(left%p,right%p) &
         .and. max(left%p,right%p)/min(left%p,right%p) .le. 2.0))then
       if(guessp .lt. min(left%p,right%p))then
          ! Two-rarefaction solution
          guessp = (&
               (aL+aR-.5*(gamma_const-1.)*(right%u-left%u))&
               /(aL/left%p**((gamma_const-1.)/(2.*gamma_const))&
               +aR/right%p**((gamma_const-1.)/(2.*gamma_const)))&
               )**(2.*gamma_const/(gamma_const-1.))
       else
          ! Two-shock solution
          gL=sqrt( 2./((gamma_const+1.)* left%rho)/(guessp+ left%p*(gamma_const-1.)/(gamma_const+1.)))
          gR=sqrt( 2./((gamma_const+1.)*right%rho)/(guessp+right%p*(gamma_const-1.)/(gamma_const+1.)))
          guessp = max(tol,&
               (gL*left%p+gR*right%p-(right%u-left%u))/(gL+gR) )
       end if
    end if
    !          guessp = .5*( PL + PR )
  end function guessp

  pure subroutine u_fun(in,pstar,f,df)
    use global_data, only: gamma_const
    use types, only: node_data
    type (node_data), intent(in) :: in
    real(8), intent(in) :: pstar
    real(8) , intent(out) :: f , df
    real(8) :: A, B, psi, a0

    psi = pstar/in%p
    a0  = sqrt(gamma_const*in%p/in%rho)

    if( psi .gt. 1. )then
       A = 2.d0/((gamma_const+1.d0)*in%rho)
       B = (gamma_const-1.d0)/(gamma_const+1.d0)*in%p
       f = in%p*(psi-1.d0)*sqrt(A/(in%p*psi+B))
       df= sqrt(A/(B+in%p*psi))*(1.d0-in%p*(psi-1.d0)/(2.*(B+psi*in%p)))
    else
       f = 2.d0*a0/(gamma_const-1.d0)*(psi**((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
       df= 1.d0/(in%rho*a0)*psi**((gamma_const+1.d0)/(-2.d0*gamma_const))
    end if

  end subroutine u_fun


  pure function beta(psi)
    use global_data, only: gamma_const
    implicit none
    real(8), intent(in) :: psi
    real(8) :: beta
    if( psi .gt. 1 )then
       beta = ( (gamma_const+1.)*psi + gamma_const-1. )&
            /( gamma_const+1. + (gamma_const-1.)*psi )
    else
       beta = psi**(1.d0/gamma_const)
    end if
  end function beta

  subroutine sample(x,left,right,Pstar,Ustar,DstarL,DstarR,out,verbose)
    use types, only: node_data
    use global_data, only: gamma_const
    use pressure_implicit
    implicit none
    type (node_data), intent(in) :: left, right
    real(8) , intent(in)  :: Pstar, Ustar, DstarL, DstarR, x
    type (node_data) , intent(out) :: out
    logical, intent(in) :: verbose

    real(8) :: PsiL , SL , SR , DeltaL, DeltaR
    real(8) :: PsiR , aL , aR , betaL , betaR
    real(8) :: cL , cR , cLT , cLH , cRT , cRH , h
    logical :: test_flag = .false.
    real, external :: zbrent
    if(verbose) test_flag = .true.

    PsiL = Pstar/left%p
    PsiR = Pstar/right%p
    aL   = sqrt(gamma_const* left%p/ left%rho)
    aR   = sqrt(gamma_const*right%p/right%rho)
    betaL= DstarL/ left%rho
    betaR= DstarR/right%rho
    SL = sqrt( left%L**2 +  left%M**2)
    SR = sqrt(right%L**2 + right%M**2)
    DeltaL =  left%A* left%M -  left%B* left%L
    DeltaR = right%A*right%M - right%B*right%L
    h = .5d0*(left%h+right%h)

    if( Ustar .gt. x )then
       !       write(*,*) "The boundary lies to the left of the contact wave"

       out%v = left%v
       if( PsiL .gt. 1.d0 )then
          !          write(*,*) " Left shock"
          cL = (1.d0-h)*left%u-aL*sqrt((gamma_const+1.d0)/(2.d0*gamma_const)*(PsiL-1.d0)+1.d0)
          if(test_flag)cL = cL*SL/DeltaL
          !          write(*,*) "Left shock speed = " , cL
          if( cL .gt. x )then
             !             write(*,*) " The boundary lies to the left of the shock"
             out%p = left%p
             out%u = left%u
             out%rho = left%rho
          else
             !             write(*,*) " The boundary lies in the left central region"
             out%p = Pstar
             out%u = Ustar
             out%rho = DstarL
          end if
       else
          !          write(*,*) "Left rarefaction wave"
          cLT = (1.d0-h)*left%u - aL
          if(test_flag) cLT = cLT*SL/DeltaL
          !          write(*,*) "Left rarefaction tail speed = " , cLT
          cLH = (1.d0-h)*Ustar - sqrt(gamma_const*Pstar/DstarL)
          if(test_flag) cLH = cLH*SL/(DeltaL*Delta_func(left%u/aL,Pstar/left%p,h,-1.d0))
          !          write(*,*) "Left rarefaction head speed = " , cLH
          if( cLT .gt. x )then
             !             write(*,*) " The boundary lies to the left of the wave"
             out%p = left%p
             out%u = left%u
             out%rho = left%rho

          elseif( cLH .lt. x )then
             !             write(*,*) "The boundary lies in the left central region"
             out%p = Pstar
             out%u = Ustar
             out%rho = DstarL
          else
             !             write(*,*) "The boundary lies within the left expansion wave"
             if(test_flag)then
                state_implicit = left
                x_implicit = x
                h_implicit = h
                plusminus  =-1.d0
                out%p = zbrent(pressure_func,real(left%p),real(Pstar),1e-7)
             else
                out%p = left%p*(2.d0*(1.d0-h)/(gamma_const-2.d0*h+1.d0) &
                     + (gamma_const-1.d0)/(aL*(gamma_const-2.d0*h+1.d0))*((1.d0-h)*left%u-x))**(2.d0*gamma_const/(gamma_const-1.d0))
             end if
             out%u = left%u - 2.d0*aL/(gamma_const-1.d0)*((out%p/left%p)**((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
             out%rho = left%rho*(out%p/left%p)**(1.d0/gamma_const)
          end if
       end if
    else
       !       write(*,*) " The boundary lies to the right of the contact wave"
       out%v = right%v
       if( PsiR .gt. 1. )then
          !          write(*,*) " Right shock"
          cR = (1.d0-h)*right%u+aR*sqrt((gamma_const+1.d0)/(2.d0*gamma_const)*(PsiR-1.d0)+1.d0)
          if(test_flag) cR = cR*SR/DeltaR
          !          write(*,*) " Right shock speed = " , cR
          if( cR .lt. x )then
             !             write(*,*) " The boundary lies to the right of the shock"
             out%p = right%p
             out%u = right%u
             out%rho = right%rho
          else
             !             write(*,*) " The boundary lies in the right central region"
             out%p = Pstar
             out%u = Ustar
             out%rho = DstarR
          end if
       else
          !          write(*,*) " Right rarefaction wave"
          cRT = (1.d0-h)*right%u + aR
          if(test_flag) cRT = cRT*SR/DeltaR
          !          write(*,*) "Right rarefaction tail speed = " , cRT
          cRH = (1.d0-h)*Ustar + sqrt(gamma_const*Pstar/DstarR)
          if(test_flag)cRH = cRH*SR/(DeltaR*Delta_func(right%u/aR,Pstar/right%p,h,1.d0))
          !          write(*,*) "Right rarefaction head speed = " , cRH
          if( cRT .lt. x )then
             !             write(*,*) " The boundary lies to the right of the wave"
             out%p = right%p
             out%u = right%u
             out%rho = right%rho

          elseif( cRH .gt. x )then
             !             write(*,*) "The boundary lies in the right central region"
             out%p = Pstar
             out%u = Ustar
             out%rho = DstarR
          else
             !             write(*,*) " The boundary lies within the right expansion wave"
             if(test_flag)then
                state_implicit = right
                x_implicit = x
                h_implicit = h
                plusminus  = 1.d0
                out%p = zbrent(pressure_func,real(right%p),real(Pstar),1e-7)
             else
                out%p = right%p*(2.*(1.d0-h)/(gamma_const-2.d0*h+1.d0)-(gamma_const-1.d0)/(aR*(gamma_const-2.d0*h+1.d0))&
                     *(1.d0-h)*right%u)**((2.d0*gamma_const)/(gamma_const-1.d0))
             end if
             out%u = right%u + 2.d0*aR/(gamma_const-1.d0)*((out%p/right%p)**((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
             out%rho = right%rho*(out%p/right%p)**(1.d0/gamma_const)
          end if
       end if
    end if
  end subroutine sample
end module riemann

!program riemann_tester
!use riemann
!use types, only: node_data
!	implicit none
!	type (node_data) :: oneL, oneR, twoL, twoR, threeL, threeR ,&
!		fourL, fourR, fiveL, fiveR, test
!	real(8) :: h0
!	  oneL%rho = 1.0D0     ;   oneL%p =    1.0D0   ;   oneL%u = 0.0D0     ;
!	  oneR%rho = 0.125D0   ;   oneR%p =    0.1D0   ;   oneR%u = 0.0D0     ;
!	  twoL%rho = 1.0D0     ;   twoL%p =    0.4D0   ;   twoL%u =-2.0D0     ;
!	  twoR%rho = 1.0D0     ;   twoR%p =    0.4D0   ;   twoR%u = 2.0D0     ;
!	threeL%rho = 1.0D0     ; threeL%p = 1000.0D0   ; threeL%u = 0.0D0     ;
!	threeR%rho = 1.0D0     ; threeR%p =    0.01D0  ; threeR%u = 0.0D0     ;
!	 fourL%rho = 1.0D0     ;  fourL%p =    0.01D0  ;  fourL%u = 0.0D0     ;
!	 fourR%rho = 1.0D0     ;  fourR%p =  100.0D0   ;  fourR%u = 0.0D0     ;
!	 fiveL%rho = 5.99924D0 ;  fiveL%p = 460.894D0  ; fiveL%u = 19.5975D0  ;
!	 fiveR%rho = 5.99242D0 ;  fiveR%p =  46.0950D0 ; fiveR%u = -6.19633D0 ;
!
!      oneL%A = 1.0d0 ;   oneL%B = 0.0d0 ;   oneL%L = 0.0d0 ;   oneL%M = 1.0d0
!      oneR%A = 1.0d0 ;   oneR%B = 0.0d0 ;   oneR%L = 0.0d0 ;   oneR%M = 1.0d0
!      twoL%A = 1.0d0 ;   twoL%B = 0.0d0 ;   twoL%L = 0.0d0 ;   twoL%M = 1.0d0
!      twoR%A = 1.0d0 ;   twoR%B = 0.0d0 ;   twoR%L = 0.0d0 ;   twoR%M = 1.0d0
!    threeL%A = 1.0d0 ; threeL%B = 0.0d0 ; threeL%L = 0.0d0 ; threeL%M = 1.0d0
!    threeR%A = 1.0d0 ; threeR%B = 0.0d0 ; threeR%L = 0.0d0 ; threeR%M = 1.0d0
!     fourL%A = 1.0d0 ;  fourL%B = 0.0d0 ;  fourL%L = 0.0d0 ;  fourL%M = 1.0d0
!     fourR%A = 1.0d0 ;  fourR%B = 0.0d0 ;  fourR%L = 0.0d0 ;  fourR%M = 1.0d0
!     fiveL%A = 1.0d0 ;  fiveL%B = 0.0d0 ;  fiveL%L = 0.0d0 ;  fiveL%M = 1.0d0
!     fiveR%A = 1.0d0 ;  fiveR%B = 0.0d0 ;  fiveR%L = 0.0d0 ;  fiveR%M = 1.0d0
!
!      oneL%A = 1.0d0 ;   oneL%B = 0.0d0 ;   oneL%L = 0.0d0 ;   oneL%M = 1.0d0
!      oneR%A = 1.0d0 ;   oneR%B = 0.0d0 ;   oneR%L = 0.0d0 ;   oneR%M = 1.0d0
!
!	h0 = 0.5d0
!
!      oneL%h = h0
!      oneR%h = h0
!      twoL%h = h0
!      twoR%h = h0
!    threeL%h = h0
!    threeR%h = h0
!     fourL%h = h0
!     fourR%h = h0
!     fiveL%h = h0
!     fiveR%h = h0
!
!	test = riemann_solve(fiveL, fiveR, verbose_flag=.true., t_out=.015d0)
!
!end program riemann_tester
