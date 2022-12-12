module readinp

  use glob_var
  use channels
  use error_close
  use util


  implicit none
  
contains
  
  subroutine inputscan(nframes,nomit,nskip,numconn,cnatnam,&
       dcoord,xyzcl,trjfm)

!===============================================================================
! Subroutine to read in the input file, CVanal.inp
!===============================================================================

    integer :: nframes,nomit,nskip,numconn
    real(kind=dp) :: dcoord
    logical :: xyzcl
    type(histfrm) :: trjfm
    character(len=8), dimension(:), allocatable :: cnatnam
    
    intent(out) :: nframes,nomit,nskip,numconn
    intent(out) :: dcoord,cnatnam,xyzcl,trjfm
 


    integer :: i
    character(len=10) :: line
    character :: formc
    integer :: ierr
 
    
!---Open the input to read
    call opench(un_analinp,'CVanal.inp','old','read','formatted')
    write(un_log,'(/,a)') "Succesfully opened CVanal.inp"


!---Read the directives and parameters
    rdln: do
       read(un_analinp,*,iostat=ierr) line

       if(ierr .ne. 0) call error(un_log,3,1)
   
       line = adjustl(line)
       line = upcase(line)
   
       backspace(un_analinp,iostat=ierr)
       if(ierr .ne. 0) call error(un_log,3,1)    

       select case(line(1:5))
          
       case("CLOSE")
          !End if the EOF is found          
          exit rdln

       case ("NOMIT")
          !Number of records to omit
          read(un_analinp,*) line, nomit
          write(un_log,'(a,i7)') "Records to omit from the&
               & start of the trjectory:", nomit
   
       case ("NSKIP")
          !Number of records to skip
          read(un_analinp,*) line, nskip
          write(un_log,'(a,i7)') "Records to skip in between reads:", &
               nskip
             

       case("NREAD")
          !Number of records to read
          read(un_analinp,*) line, nframes
          write(un_log,'(a,i7)') "Total number of records to read:", &
               nframes


       case("CONNE")
          !Number of atom types for the connectivity lists
          read(un_analinp,*) line, numconn
          write(un_log,'(a,i7)') "Atomic species for CV analysis:",&
               numconn

          !Store the atomic labels
          allocate(cnatnam(numconn))

          do i=1, numconn
             read(un_analinp,*,iostat=ierr) cnatnam(i)
             cnatnam(i)=adjustl(cnatnam(i))
             cnatnam(i)=upcase(cnatnam(i))
             if(ierr .ne. 0) call error(un_log,3,1)
          end do

          write(un_log,'(7a8)') cnatnam(1:numconn)


       case("COORD")
          !Truncation distance for atom connections
          read(un_analinp,*,iostat=ierr) line, dcoord
          if(ierr .ne. 0) call error(un_log,3,1)
          write(un_log,'(a,f7.2)') "Distance cut-off for &
               &coordinating atoms:", dcoord
          

       case("P_XYZ")
          !Chose whether to print the xyz file
          read(un_analinp,*,iostat=ierr) line, xyzcl
          if(ierr .ne. 0) call error(un_log,3,1)


       case("FORMA")
          !Trajectory format
          read(un_analinp,*,iostat=ierr) line,trjfm%vrsn,formc
          if(ierr .ne. 0) call error(un_log,3,1)


          formc=upcase(formc)
          if(formc .eq. 'F') then
             trjfm%frm=.true.
          elseif(formc .eq. 'U') then
             trjfm%frm=.false.
          else
             if(ierr .ne. 0) call error(un_log,3,1)
          end if

  
       case default
          read(un_analinp,*,iostat=ierr) line
          if(ierr .ne. 0) call error(un_log,3,1)

       end select
    end do rdln

        
    call closech(un_analinp)


    
  end subroutine inputscan
  
end module readinp
