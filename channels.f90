module channels

 implicit none

!------------------------------------------------------------------------------- 
! Subroutines to open input and output channels 
!------------------------------------------------------------------------------- 
 
 
contains


  subroutine opench(chanid,chname,stat,act,type)

    integer, intent(in) :: chanid
    character(len=*), intent(in) :: chname,stat,act,type
    integer :: ierr

    open(unit=chanid,file=chname,status=stat,action=act,form=type,iostat=ierr)
    if (ierr .ne. 0) then
       write(*,*) "Unable to open the input file:", chname
       write(*,*) "Killing the job"
       stop
    end if

  end subroutine opench




  subroutine closech(chanid)

    integer, intent(in) :: chanid

    close(chanid)

  end subroutine closech




end module channels
