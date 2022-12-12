module cvatdefine
  
  use glob_var
  use error_close

  implicit none



contains
  subroutine cvatdef(natms,atlst,nacn,cvlst)
    
!===============================================================================
! Define the atoms identified for CV analysis
!===============================================================================

    integer :: natms,nacn
    type(atom_list),dimension(:) :: atlst 
    type(atom_list),dimension(:),allocatable :: cvlst 

    intent(in) :: nacn,natms
    intent(inout) :: atlst
    intent(out) :: cvlst

    integer :: i
    integer :: natmcl

!---Initialise 
    allocate(cvlst(nacn))
    natmcl=0


!---Complete the information in cvlst arrays

    do i=1, natms
       
       !Populate the cvlst arrays
       if(atlst(i)%iclu .gt. 0) then
          natmcl=natmcl+1
          
          cvlst(natmcl) = atlst(i)
          
       end if
    end do
       
  end subroutine cvatdef

end module cvatdefine
