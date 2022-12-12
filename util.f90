module util

  use glob_var
  
  implicit none

contains  

  function upcase(string)

!==============================================================
! Function to convert character string to uppercase
!==============================================================
    
    character(len=*), intent(in):: string
    character(len=len(string)) :: re_string
    character(len=len(string)) :: upcase
    integer :: j
    
    re_string=trim(adjustl(string))
    
    do j = 1,len(re_string)
       if(re_string(j:j) >= "a" .and. re_string(j:j) <= "z") then
          upcase(j:j) = achar(iachar(re_string(j:j)) - 32)
       else
          upcase(j:j) = re_string(j:j)
       end if
    end do
    
    
  end function upcase


  function atomcount(iat,iarr)
    

!===============================================================================
! Function to count the number of atomic species requested
!===============================================================================
    
    integer, intent(in) :: iat
    integer, dimension(:), intent(in) :: iarr
    integer :: atomcount

!---Count the number of atoms with specified identity

    atomcount = count(mask = iarr(:) .eq. iat) 

  end function atomcount
    


  subroutine atompick(iat,iarr,oarr)

!===============================================================================
! Routine to provide an index array for particular atomic species
!===============================================================================
    
    intent(in) :: iat,iarr
    intent(out) :: oarr

    integer :: iat
    integer, dimension(:) :: iarr,oarr
 
    integer :: i, seq

    seq = 0

!---Index relevant atoms to recall as necessary

    do i=1, size(iarr(:))
       if(iarr(i) .eq. iat) then
          seq = seq + 1
          oarr(seq) = i
       end if
    end do
    
  end subroutine atompick

  function cross(a,b)

!===============================================================================
!
! Provides the cross product of a and b
!
!===============================================================================

    real(dp), dimension(3) :: a, b, cross
    intent(in) :: a, b
    
    
    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
    
  end function cross
  
   
  
  function rotate(a,b,theta)

!===============================================================================
!
! Rotates a around b by the angles theta
!
!===============================================================================

    real(dp), dimension(3) :: a, b, rotate
    real(dp) :: theta
    intent(in) :: a, b
    
    
    real(dp) :: co, si
    
    si = sin(theta)
    co = cos(theta)
    
    
    !Initialise
    rotate = 0.0_dp
    
    
    rotate(1)=             (a(1)*a(1) + (1.0 - a(1)*a(1))*co) *b(1)    
    rotate(1)= rotate(1) + (a(1)*a(2)*(1.0 - co) - a(3)*si) *b(2)
    rotate(1)= rotate(1) + (a(1)*a(3)*(1.0 - co) + a(2)*si) *b(3)
    
    rotate(2)=             (a(1)*a(2)*(1.0 - co) + a(3)*si) *b(1)
    rotate(2)= rotate(2) + (a(2)*a(2) + (1.0 - a(2)*a(2))*co) *b(2)
    rotate(2)= rotate(2) + (a(2)*a(3)*(1.0 - co) - a(1)*si) *b(3)
    
    rotate(3)=             (a(1)*a(3)*(1.0 - co) - a(2)*si) *b(1)
    rotate(3)= rotate(3) + (a(2)*a(3)*(1.0 - co) + a(1)*si) *b(2)
    rotate(3)= rotate(3) + (a(3)*a(3) + (1.0 - a(3)*a(3))*co) *b(3)  
    
    
    
  end function rotate
  

  subroutine invert(a,b,d)

!===============================================================================
!     
!     Adapted from DL routine to invert a 3 * 3 matrix using cofactors
!     
!===============================================================================

    real(dp),dimension(9) :: a,b
    real(dp) :: d,r

    intent(in) :: a
    intent(out) :: b

    ! calculate adjoint matrix
    b(1)=a(5)*a(9)-a(6)*a(8)
    b(2)=a(3)*a(8)-a(2)*a(9)
    b(3)=a(2)*a(6)-a(3)*a(5)
    b(4)=a(6)*a(7)-a(4)*a(9)
    b(5)=a(1)*a(9)-a(3)*a(7)
    b(6)=a(3)*a(4)-a(1)*a(6)
    b(7)=a(4)*a(8)-a(5)*a(7)
    b(8)=a(2)*a(7)-a(1)*a(8)
    b(9)=a(1)*a(5)-a(2)*a(4)
     
    ! calculate determinant
    d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
    r=0.0_dp
    
    if(abs(d).gt.0.0_dp) r=1.0_dp/d

    !complete inverse matrix
    b(1)=r*b(1)
    b(2)=r*b(2)
    b(3)=r*b(3)
    b(4)=r*b(4)
    b(5)=r*b(5)
    b(6)=r*b(6)
    b(7)=r*b(7)
    b(8)=r*b(8)
    b(9)=r*b(9)
    
  end subroutine invert

end module util
