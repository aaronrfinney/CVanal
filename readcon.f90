module readcon

  use glob_var
  use error_close
  use channels
  use util

  implicit none
  
contains

  subroutine atcnt(natms)

!===============================================================================
! Subroutine to count the number of atoms in CONFIG
!===============================================================================

    integer, intent(out) :: natms
  
    character(len=80) :: header    
    character(len=8) :: name
    integer :: j, ierr, levcfg, imcon


!---Count the number of atoms in CONFIG       

    !Initialise
    natms=0
    
    !---Open CONFIG
    call opench(un_config,'CONFIG','old','read','formatted')
    write(un_log,'(/,a)') "Succesfully opened CONFIG"    
    
    read(un_config,'(a80)',iostat=ierr) header
    if(ierr .ne. 0) call error(un_log,5,1)
    read(un_config,'(2i10)',iostat=ierr) levcfg,imcon
    if(ierr .ne. 0) call error(un_log,5,1)
    
    
    if(imcon .gt. 0) then
       read(un_config,'(//)',iostat=ierr) 
    end if
    
    !---Count the number of atoms
    scancfg: do
       
       read(un_config,'(a8)',iostat=ierr) name
       if(ierr .lt. 0) then
          exit scancfg
       elseif(ierr .gt. 0) then
          call error(un_log,5,1)
       else
          natms = natms + 1
       end if
       
       read(un_config,*)
       
       do j=1,levcfg
          read(un_config,*)
       end do
       
    end do scancfg
    
    write(un_log,'(a,i7)') "Total number of atoms &
         &in CONFIG:", natms

    call closech(un_config)

  end subroutine atcnt




  
  subroutine configscan(natms,nuniq,auniq,ain,atlab,atind,clind,aicl,&
       nacn,cnatnam)
    
    
!===============================================================================
! Subroutine to read CONFIG and extract unique atom types.
!
! And identify atoms for the connectivity lists.
!===============================================================================


    integer :: natms,nuniq,nacn
    character(len=8), dimension(:) :: atlab,auniq
    integer, dimension(:) :: ain,atind,clind,aicl
    character(len=8), dimension(:) :: cnatnam
    
    intent(in) :: natms,cnatnam
    intent(out) :: nuniq,auniq,ain,atlab,atind,clind,aicl,nacn
    

    character(len=80) :: header    
    character(len=8) :: name
    integer :: i, j, k, ierr
    integer :: levcfg, imcon
    integer :: nconty

    logical :: uniquetype

    !Number of atom types which can be connected
    nconty=size(cnatnam)

    !Initialise
    auniq="        "
    nuniq=0
    atind=0
    ain=0
    clind=0
    aicl=0
    nacn=0
     
    
!---Open CONFIG
    call opench(un_config,'CONFIG','old','read','formatted')
     
    read(un_config,'(a80)',iostat=ierr) header 
    if(ierr .ne. 0) call error(un_log,5,1)
    read(un_config,'(2i10)',iostat=ierr) levcfg,imcon
    if(ierr .ne. 0) call error(un_log,5,1)
    
       
    if(imcon .gt. 0) then
       read(un_config,'(//)',iostat=ierr) 
    end if
    

!---Get atom names and indices

    do i=1, natms
       read(un_config,'(a8)',iostat=ierr) name
       if(ierr .ne. 0) call error(un_log,5,1) 

       atlab(i) = adjustl(name)
       ain(i)=i

!---Idenify unique atoms

       uniquetype = .true.

       do j=1,limit
          if(name .eq. auniq(j)) then 
             uniquetype = .false.
          end if

       end do

       if(uniquetype) then

          nuniq = nuniq + 1
          auniq(nuniq) = name

          atind(i)=nuniq

          if(nuniq .eq. limit) call error(un_log,6,1)

       else
          
          do k=1, nuniq
             if(name .eq. auniq(k)) then            
                atind(i)=k
             end if
          end do
       end if

!---Check if this is a CV atom and initialise log numbers
 
       uncon: do j=1, nconty

          if(upcase(name) .eq. cnatnam(j)) then

             clind(i)=1 !All CV atoms logged as 1
             aicl(i)=j  
             nacn=nacn+1 

             exit uncon

          end if
       end do uncon



!---Skip the velocities and forces
       
       do j=1, (levcfg+1)
          read(un_config,*)
       end do

    end do

!---Write the number of unique atom types

    write(un_log,'(a32,i7)') "Number of atom types&
         & identified:", nuniq
    !write(un_log,'(7a8)') auniq(1:nuniq)

    write(un_log,'(a32,i7)') "Number of CV atoms:", nacn


!---Close the CONFIG file

    call closech(un_config)

  end subroutine configscan

        
end module readcon
