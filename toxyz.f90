module toxyz

  use glob_var
  use channels
  use error_close
  implicit none

  contains
!===============================================================================
! Produces xyz appended files which can be read by VMD
!===============================================================================

    subroutine writexyz(natms,heading,atlab,xyz)

      intent(in) :: natms,heading,atlab,xyz


      integer :: natms
      character(len=*) :: heading
      character(len=8), dimension(:) :: atlab
      real(dp), dimension(:,:) :: xyz

      integer :: i, ierr
      logical :: conn

      inquire(unit=un_xyz,opened=conn)

      if(.not. conn) then
!---Open the output file         

         call opench(un_xyz,"CVanal.xyz","replace","write","formatted")

         write(un_log,'(/,a)') "Succesfully opened the XYZ output"
 

      end if
!---Write the frame in xyz format

      write(un_xyz,'(i6)',iostat=ierr) natms
      if(ierr .ne. 0) call error(un_log,12,1)

      write(un_xyz,'(a80)',iostat=ierr) adjustl(heading)
      if(ierr .ne. 0) call error(un_log,12,1)

      do i=1, natms
         write(un_xyz,'(a8,1x,3(f10.5,1x))',iostat=ierr) &
              atlab(i), xyz(:,i)
         if(ierr .ne. 0) call error(un_log,12,1)
      end do
      
      
      
    end subroutine writexyz

end module toxyz

















