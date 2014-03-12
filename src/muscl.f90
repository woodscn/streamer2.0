module muscl_update
  implicit none
contains
  subroutine muscl(f1,f2,f11,f22)!,fL,fR)
    use types, only: node_data
    implicit none
    type (node_data) :: f1, f2, f11, f22, fL, fR
    !	intent(in) :: f1, f2, f11, f22
    !	intent(out):: fL, fR

    fL = f1 ; fR = f2
    fL%rho = funcL(f1%rho,f2%rho,f11%rho) ; fR%rho = funcR(f1%rho,f2%rho,f22%rho)
    fL%p   = funcL(f1%p  ,f2%p  ,f11%p  ) ; fR%p   = funcR(f1%p  ,f2%p  ,f22%p  )
    fL%u   = funcL(f1%u  ,f2%u  ,f11%u  ) ; fR%u   = funcR(f1%u  ,f2%u  ,f22%u  )
    fL%v   = funcL(f1%v  ,f2%v  ,f11%v  ) ; fR%v   = funcR(f1%v  ,f2%v  ,f22%v  )

    fL%theta = 0.0d0 ; fL%q = 0.0d0 ; fL%g = 0.0d0 ; fL%x = 0.0d0 ; fL%y = 0.0d0
    fR%theta = 0.0d0 ; fR%q = 0.0d0 ; fR%g = 0.0d0 ; fR%x = 0.0d0 ; fR%y = 0.0d0

    f1=fL ; f2=fR

  end subroutine muscl

  real(8) function funcR(f1,f2,f22)
    implicit none
    real(8) :: f1, f2, f22, r
    r = (f2-f1)/(f22-f2)
    funcR = f2 - 0.5d0*(f22-f2)*phi(r)
  end function funcR

  real(8) function funcL(f1,f2,f11)
    implicit none
    real(8) :: f11, f1, f2, r
    r = (f2-f1)/(f1-f11)
    funcL = f1 + 0.5d0*(f1-f11)*phi(r)
  end function funcL


  real(8) function phi(r)
    implicit none
    real(8) :: r, b
    phi = max(0.d0,min(1.d0,r))
!    phi = (r+abs(r))/(1.d0+abs(r))
!    phi = maxval([0.d0,min(2.d0*r,1.d0),min(r,2.d0)])
!    phi = max(0.d0,minval([2.d0*r,(1.d0+2.d0*r)/3.d0,2.d0]))
!    phi = max(0.d0,minval([2.d0*r,.5d0*(1.d0+r),2.d0]))
!    phi = 1.5d0*(r**2+r)/(r**2+r+1.d0)
!    phi = (r**2+r)/(r**2+1.d0)
!    phi = 1.5d0*(r+abs(r))/(r+2.d0)
!    b = 1.5
!    phi = maxval([0.d0,min(b*r,1.d0),min(r,b)])
  end function phi
end module muscl_update
