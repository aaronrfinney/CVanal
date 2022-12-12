module readhist

  use glob_var
  use error_close
  use channels


  implicit none
  
contains

  subroutine readtraj_f(ird,natms,atlab,rxyz,nstp,tstp,&
       cell,imcon,nomit,nskip,skip)

!===============================================================================
! Subroutine to read a formatted DL 2 trajectory file
!===============================================================================

    intent(in) :: natms,atlab,nskip,nomit,skip
    intent(out) :: rxyz,nstp,tstp,cell,imcon
    intent(inout) :: ird

    integer, save :: frcount
    integer :: ird,natms,nstp,nomit,nskip,imcon
    real(dp) :: mass,chrg
    real(dp), dimension(:,:) :: rxyz
    real(dp), dimension(:) :: cell
    real(dp) :: tstp
    character(len=*), dimension(:) :: atlab
    logical :: skip
     
    integer :: ierr, i, j, ai
    integer :: keytrj,atoms
    character(len=80) :: trjheader
    character(len=8) :: string
    character(len=8), dimension(natms) :: label
    logical :: readrec


  
    readrec=.true.

    if(ird .eq. 0) then
!---Read the preamble
       frcount=0
       call opench(un_hist,'HISTORY','old','read','formatted')
       write(un_log,'(/,a)') "Succesfully reading HISTORY:"
       
       read(un_hist,'(a80)',iostat=ierr) trjheader
       if(ierr .ne. 0) call error(un_log,8,1)
       write(un_log,'(3x,a77)') adjustl(trjheader)
       
       read(un_hist,'(3i10)',iostat=ierr) keytrj,imcon,atoms
       if(ierr .ne. 0) call error(un_log,8,1)

       !Check to ensure # atoms is expected
       if(atoms .ne. natms) call error(un_log,9,1)
  
       ird=1

       !Omit records if necessary
       call fframeskip(ird,nomit)
 
    end if

    
!---Read the data

    read(un_hist,'(a8,4i10,f12.6)',iostat=ierr) &
         string, nstp, atoms, keytrj, imcon, tstp
    
    
    
    !Check for error/EOF
    if(ierr .gt. 0) then
       call error(un_log,8,1)
    elseif(ierr .lt. 0) then
       ird = -1
       readrec = .false.
    end if
    
    if(readrec) then
       frcount=frcount+1
       
       if(imcon .gt. 0) then
          read(un_hist,'(3g12.4)',iostat=ierr) cell(:)
          if(ierr .ne. 0) call error(un_log,8,1)
       end if
       

       do i=1, natms
          read(un_hist,'(a8,i10,2f12.6)',iostat=ierr) label(i),&
               ai, mass, chrg
          if(ierr .ne. 0) call error(un_log,8,1)
          
          !Check that atom type is expected
          if(frcount .eq. 1) then
             if(label(i) .ne. atlab(i)) then
                call error(un_log,7,1)
             end if
          end if
          
          read(un_hist,'(1p,3e12.4)',iostat=ierr) rxyz(:,i)
          if(ierr .ne. 0) call error(un_log,8,1)
                      
          do j=1, keytrj
             read(un_hist,*,iostat=ierr)
             if(ierr .ne. 0) call error(un_log,8,1)
          end do
          
       end do
       
       !---Skip frames
       if(skip) then
          call fframeskip(ird,nskip)
       end if
       
    end if
 
  end subroutine readtraj_f



  subroutine fframeskip(ird,nskip)

!-------------------------------------------------------------------------------
! Skips frames in the formatted DL 2 trajectory
!-------------------------------------------------------------------------------

    intent(inout) :: ird
    intent(in) :: nskip

    integer :: ird, nskip

    character(len=8) :: string
    real(dp) :: tstp
    integer :: keytrj, imcon, atoms,nstp
    integer :: i, j, ierr, lnskp
    logical :: readf

    readf=.true.
  

!---Skip frames
   
    do i=1, nskip
       read(un_hist,'(a8,4i10,f12.6)',iostat=ierr) &
            string, nstp, atoms, keytrj, imcon, tstp
       
       !Check for error/EOF
       if(ierr .gt. 0) then
          call error(un_log,8,1)
       elseif(ierr .lt. 0) then
          ird=-1
          readf=.false.
       end if
       
       if(readf) then
            
          if(imcon .gt. 0) then
             read(un_hist,'(//)') 
          end if
          lnskp=atoms*(keytrj+2)
          do j=1,lnskp
             read(un_hist,*) 
          end do
       end if
    end do
    
  end subroutine fframeskip




  subroutine readtraj_u(ird,natms,atlab,rxyz,nstp,tstp,&
       cell,imcon,nomit,nskip,skip)
   
!===============================================================================
! Subroutine to read the unformatted DL 2 trajectory file
!===============================================================================
      
    intent(in) :: natms,atlab,nskip,nomit,skip
    intent(out) :: rxyz,nstp,tstp,cell,imcon
    intent(inout) :: ird

    integer :: ird,natms,nstp,nomit,nskip,imcon
    real(dp), dimension(natms) :: mass, chrg
    real(dp), dimension(:,:) :: rxyz
    real(dp), dimension(:) :: cell
    real(dp) :: tstp
    character(len=*), dimension(:) :: atlab
    logical :: skip
    
    integer :: ierr, i, j
    integer :: keytrj, atoms
    real(dp) :: rnstp, ratoms, rkeytrj, rimcon
    character(len=80) :: trjheader
    character(len=8), dimension(natms) :: label
    logical :: readrec
    integer, save :: frcount
    
 
    readrec=.true.
    
    if(ird .eq. 0) then
!---Read the preamble
       frcount=0      
       call opench(un_hist,'HISTORY','old','read','unformatted')
       write(un_log,'(/,a)') "Succesfully reading HISTORY"
       
       read(un_hist,iostat=ierr) trjheader
       if(ierr .ne. 0) call error(un_log,8,1)
       write(un_log,'(3x,a80)') adjustl(trjheader)

       
       read(un_hist,iostat=ierr) ratoms
       if(ierr .ne. 0) call error(un_log,8,1)

       atoms=int(ratoms)
       
       !Check to ensure # atoms is expected
       if(atoms .ne. natms) call error(un_log,9,1)

       read(un_hist,iostat=ierr) (label(i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)

       !Check atom type is expected
       do i=1,natms
          if(label(i) .ne. atlab(i)) then
             call error(un_log,7,1)
          end if
       end do
       
       read(un_hist,iostat=ierr) (mass(i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)
       read(un_hist,iostat=ierr) (chrg(i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)
       
         
       ird=1

       !Omit records if necessary
       call uframeskip(ird,nomit)


    end if


!---Read the data

    read(un_hist,iostat=ierr) rnstp, ratoms, &
         rkeytrj, rimcon, tstp
    
    nstp=int(rnstp)
    keytrj=int(rkeytrj)
    imcon=int(rimcon)
    
    
    !Check for error/EOF
    if(ierr .gt. 0) then
       call error(un_log,8,1)
    elseif(ierr .lt. 0) then
       ird=-1
       readrec=.false.
    end if
    
    if(readrec) then
       frcount=frcount+1
       
       if(imcon .gt. 0) then
          read(un_hist,iostat=ierr) cell(:)
          if(ierr .ne. 0) call error(un_log,8,1)
       end if
       
       read(un_hist,iostat=ierr) (rxyz(1,i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)
       read(un_hist,iostat=ierr) (rxyz(2,i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)
       read(un_hist,iostat=ierr) (rxyz(3,i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)
       
       do j=1, 3*keytrj
          read(un_hist,iostat=ierr) 
          if(ierr .ne. 0) call error(un_log,8,1)
       end do
       
       
       !---Skip frames
       if(skip) then
          call uframeskip(ird,nskip)
       end if
       
    end if
 
  end subroutine readtraj_u


  subroutine uframeskip(ird,nskip)

!-------------------------------------------------------------------------------
! Skips frames in the unformatted DL 2 trajectory
!-------------------------------------------------------------------------------    
    intent(inout) :: ird
    intent(in) :: nskip

    integer :: ird, nskip
    
    real(dp) :: tstp, nstp, atoms, rkeytrj, imcon
    integer :: i, j, lnskp, ierr, keytrj 
    logical :: readf

    readf=.true.
    
!---Skip frames


    
    do i=1, nskip
       
       read(un_hist,iostat=ierr) &
            nstp, atoms, rkeytrj, imcon, tstp

       keytrj=int(rkeytrj)
       
       !Check for error/EOF
       if(ierr .gt. 0) then
          call error(un_log,8,1)
       elseif(ierr .lt. 0) then
          ird=-1
          readf=.false.
       end if
       
       if(readf) then
          
          if(imcon .gt. 0) then
             read(un_hist) 
          end if
          lnskp=(keytrj+1)*3
          do j=1,lnskp
             read(un_hist)
          end do
       end if
    end do
    
  end subroutine uframeskip



  
  subroutine readtraj_dl4_f(ird,natms,atlab,rxyz,nstp,tstp,&
       cell,imcon,nomit,nskip,skip)

!===============================================================================
! Subroutine to read a DL 4 trajectory file
!===============================================================================

    intent(in) :: natms,atlab,nskip,nomit,skip
    intent(out) :: rxyz,nstp,tstp,cell,imcon
    intent(inout) :: ird

    integer, save :: frcount
    integer :: ird,natms,nstp,nomit,nskip,imcon
    real(dp) :: mass,chrg
    real(dp), dimension(:,:) :: rxyz
    real(dp), dimension(:) :: cell
    real(dp) :: tstp
    character(len=*), dimension(:) :: atlab
    logical :: skip
     
    integer :: ierr, i, j, ai
    integer :: keytrj,atoms
    character(len=80) :: trjheader
    character(len=8) :: string
    character(len=8), dimension(natms) :: label
    logical :: readrec


  
    readrec=.true.

    if(ird .eq. 0) then
!---Read the preamble
       frcount=0
       call opench(un_hist,'HISTORY','old','read','formatted')
       write(un_log,'(/,a)') "Succesfully reading HISTORY:"
       
       read(un_hist,'(a72)',iostat=ierr) trjheader
       if(ierr .ne. 0) call error(un_log,8,1)
       write(un_log,'(3x,a77)') adjustl(trjheader)

              
       read(un_hist,'(3i10)',iostat=ierr) keytrj,imcon,atoms
       if(ierr .ne. 0) call error(un_log,8,1)

       
       !Check to ensure # atoms is expected
       if(atoms .ne. natms) call error(un_log,9,1)
  
       ird=1

       !Omit records if necessary
       
       call fframeskip_dl4(ird,nomit)

       
 
    end if

    
!---Read the data

    read(un_hist,'(a8,2i10,2i2,f20.6)',iostat=ierr) &
         string, nstp, atoms, keytrj, imcon, tstp

        
    
    !Check for error/EOF
    if(ierr .gt. 0) then
       call error(un_log,8,1)
    elseif(ierr .lt. 0) then
       ird = -1
       readrec = .false.
    end if

        
    if(readrec) then
       frcount=frcount+1

       if(imcon .gt. 0) then
          read(un_hist,'(3g20.10)',iostat=ierr) cell(:)
          if(ierr .ne. 0) call error(un_log,8,1)

          
       end if
       

       do i=1, natms
          read(un_hist,'(a8,i10,2f12.6)',iostat=ierr) label(i),&
               ai, mass, chrg
          if(ierr .ne. 0) call error(un_log,8,1)

          !print *, label(i),ai, mass, chrg
          
          !Check that atom type is expected
          if(frcount .eq. 1) then
             if(label(i) .ne. atlab(i)) then
                call error(un_log,7,1)
             end if
          end if
          
          read(un_hist,'(3f20.8)',iostat=ierr) rxyz(:,i)
          if(ierr .ne. 0) call error(un_log,8,1)

          !write(*,'(3f20.8)') rxyz(:,i)
                      
          do j=1, keytrj
             read(un_hist,*,iostat=ierr)
             if(ierr .ne. 0) call error(un_log,8,1)
          end do
          
       end do
       
       !---Skip frames
       if(skip) then
          call fframeskip_dl4(ird,nskip)
       end if
       
    end if
 
  end subroutine readtraj_dl4_f



  subroutine fframeskip_dl4(ird,nskip)

!-------------------------------------------------------------------------------
! Skips frames in the DL 4 trajectory
!-------------------------------------------------------------------------------

    intent(inout) :: ird
    intent(in) :: nskip

    integer :: ird, nskip

    character(len=8) :: string
    real(dp) :: tstp
    integer :: keytrj, imcon, atoms,nstp
    integer :: i, j, ierr, lnskp
    logical :: readf

    readf=.true.
    

!---Skip frames
   
    do i=1, nskip
       read(un_hist,'(a8,2i10,2i2,f20.6)',iostat=ierr) &
            string, nstp, atoms, keytrj, imcon, tstp

       !Check for error/EOF
       if(ierr .gt. 0) then
          call error(un_log,8,1)
       elseif(ierr .lt. 0) then
          ird=-1
          readf=.false.
       end if
       
       if(readf) then
            
          if(imcon .gt. 0) then
             read(un_hist,'(//)') 
          end if
          lnskp=atoms*(keytrj+2)
          do j=1,lnskp
             read(un_hist,*) 
          end do
       end if
    end do
    
  end subroutine fframeskip_dl4




  subroutine readtraj_dlp4_u(ird,natms,atlab,rxyz,nstp,tstp,&
       cell,imcon,nomit,nskip,skip)
   
!===============================================================================
! Subroutine to read an unformatted DL 4 trajectory file
!===============================================================================
      
    intent(in) :: natms,atlab,nskip,nomit,skip
    intent(inout) :: rxyz,nstp,tstp,cell,imcon
    intent(inout) :: ird

    integer :: ird,natms,nstp,nomit,nskip,imcon
    real(dp), dimension(natms) :: mass
    real(dp),dimension(natms) :: rsd
    real(dp), dimension(3,natms) :: rxyz
    real(dp), dimension(9) :: cell
    real(dp) :: tstp,time
    character(len=8), dimension(natms) :: atlab
    logical :: skip
    
    integer :: ierr,i,j
    integer :: keytrj,atoms,records
    character(len=72) :: trjheader
    character(len=8), dimension(:),allocatable :: label
    real(dp), dimension(:),allocatable :: chrg
    logical :: readrec
    integer, save :: frcount
    

    allocate(label(natms))
    allocate(chrg(natms))
    readrec=.true.
    
    if(ird .eq. 0) then
!---Read the preamble
       frcount=0      
       call opench(un_hist,'HISTORY','old','read','unformatted')
       write(un_log,'(/,a)') "Succesfully reading HISTORY"


       read(un_hist,iostat=ierr) trjheader
       if(ierr .ne. 0) call error(un_log,8,1)
       write(un_log,'(3x,a72)') adjustl(trjheader)
  
       read(un_hist,iostat=ierr) keytrj,imcon,atoms,nstp,records
       if(ierr .ne. 0) call error(un_log,8,1)
  
       !Check to ensure # atoms is expected
       if(atoms .ne. natms) call error(un_log,9,1)


       read(un_hist,iostat=ierr) (label(i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)


       !Check atom type is expected
       do i=1,natms
          if(label(i) .ne. atlab(i)) then
             call error(un_log,7,1)
          end if
       end do
  
       read(un_hist,iostat=ierr) (mass(i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)


       read(un_hist,iostat=ierr) (chrg(i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)

  
       read(un_hist,iostat=ierr) (rsd(i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)
        
         
       ird=1

       !Omit records if necessary
       call uframeskip_dlp4(ird,nomit)
  

    end if
 
!---Read the data

    read(un_hist,iostat=ierr) nstp,atoms,&
         keytrj,imcon,tstp,time

    
    
    !Check for error/EOF
    if(ierr .gt. 0) then
       call error(un_log,8,1)
    elseif(ierr .lt. 0) then
       ird=-1
       readrec=.false.
    end if
    
    if(readrec) then
       frcount=frcount+1
       
       if(imcon .gt. 0) then
          read(un_hist,iostat=ierr) cell(:)
          if(ierr .ne. 0) call error(un_log,8,1)
       end if
       
       read(un_hist,iostat=ierr) (rxyz(1,i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)
       read(un_hist,iostat=ierr) (rxyz(2,i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)
       read(un_hist,iostat=ierr) (rxyz(3,i),i = 1,natms)
       if(ierr .ne. 0) call error(un_log,8,1)

       
       do j=1, 3*keytrj
          read(un_hist,iostat=ierr) 
          if(ierr .ne. 0) call error(un_log,8,1)
       end do
       
       
       !---Skip frames
       if(skip) then
          call uframeskip_dlp4(ird,nskip)
       end if
       
    end if

    deallocate(label,chrg)

 
  end subroutine readtraj_dlp4_u


  subroutine uframeskip_dlp4(ird,nskip)

!-------------------------------------------------------------------------------
! Skips frames in the unformatted DL 4 trajectory
!-------------------------------------------------------------------------------    
    intent(inout) :: ird
    intent(in) :: nskip

    integer :: ird, nskip
    
    real(dp) :: tstp,time
    integer :: nstp,atoms,keytrj,imcon
    integer :: i, j, lnskp, ierr 
    logical :: readf

    readf=.true.
    
!---Skip frames


    
    do i=1, nskip
       
       read(un_hist,iostat=ierr) &
            nstp,atoms,keytrj,imcon,tstp,time

       !Check for error/EOF
       if(ierr .gt. 0) then
          call error(un_log,8,1)
       elseif(ierr .lt. 0) then
          ird=-1
          readf=.false.
       end if
       
       if(readf) then
          
          if(imcon .gt. 0) then
             read(un_hist) 
          end if
          lnskp=(keytrj+1)*3
          do j=1,lnskp
             read(un_hist)
          end do
       end if
    end do
    
  end subroutine uframeskip_dlp4


end module readhist
