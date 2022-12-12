module readfld
  
  use glob_var
  use error_close
  use channels
  use util
  

contains  
  subroutine fieldscan(mdef,atlst,mlst,numol,natms,numconn,cnatnam,nmcn)

!===============================================================================
! A routine to read through a DL_POLY FIELD file and extract moleular
! information. 
!
!===============================================================================

    type(moldefine) :: mdef
    type(atom_list),dimension(:) :: atlst
    type(mol_list),dimension(:),allocatable :: mlst
    integer :: numol,natms,numconn,nmcn
    character(len=8), dimension(:) :: cnatnam

    intent(in) :: natms,numconn,cnatnam
    intent(out) :: mdef,atlst,mlst,numol,nmcn

    character(len=20) :: line
    character(len=8) :: formspec
    character(len=8), dimension(:), allocatable :: tmp_nam
    
    real(dp),dimension(:),allocatable :: tmp_mass, tmp_chrg
    integer,dimension(:,:,:), allocatable :: tmp_lbnd
    integer,dimension(:),allocatable :: molcntyp

    integer :: iat,iatom,mxbnd,jprev,nmax,irpt
    integer :: mxbnd_old,natm_f
    integer :: ierr,i,j,k
    integer :: im,nfirst
    logical :: cont

    iatom = 0
    mxbnd = 0


!---Open the FIELD file

    call opench(un_fld,"FIELD","old","read","formatted")
    write(un_log,'(/,a)') "Succesfully opened FIELD"

!---Define molecules

    !Find the number of molecules
    do
       read(un_fld,*,iostat=ierr) line
       if(ierr .ne. 0) call error(un_log,5,1)

       line = adjustl(line)
       line = upcase(line)
       if(line(1:4) .eq. "MOLE") then
          exit
       end if
    end do
    
    !Allocate memory for molecular information and initialise

    backspace(un_fld,iostat=ierr)
    if(ierr .ne. ierr) call error(un_log,5,1) 
    
    read(un_fld,*,iostat=ierr) line, mdef%nmtp
    if(ierr .ne. ierr) call error(un_log,5,1)

    allocate(mdef%nmol(mdef%nmtp)) 
    allocate(mdef%nmat(mdef%nmtp))
    allocate(mdef%nbnd(mdef%nmtp)) 
    allocate(mdef%mname(mdef%nmtp))
    allocate(mdef%nadd(mdef%nmtp))
    allocate(mdef%iref(3,mdef%nmtp))
    
    mdef%nmol=0
    mdef%nmat=0
    mdef%nbnd=0
    mdef%mname=""
    mdef%iref=0
    mdef%nadd=0

    allocate(molcntyp(mdef%nmtp))
    molcntyp = 0


    !Obtain definitions for each molecule type

    do i=1,mdef%nmtp

       !Get the name, # of molecules and # of atoms for this molecular type
       read(un_fld,*,iostat=ierr) mdef%mname(i)
       if(ierr .ne. 0) call error(un_log,5,1)

       mdef%mname(i)=adjustl(mdef%mname(i))
       mdef%mname(i)=upcase(mdef%mname(i))

       read(un_fld,*,iostat=ierr) line, mdef%nmol(i)
       if(ierr .ne. 0) call error(un_log,5,1)

       read(un_fld,*,iostat=ierr) line, mdef%nmat(i)
       if(ierr .ne. 0) call error(un_log,5,1)

       if(mdef%nmat(i) .ge. 3) mdef%iref(:,i)=(/ 1 , 2 , 3 /)

       !Read atomic weights
       
       allocate(tmp_mass(mdef%nmat(i)))
       allocate(tmp_chrg(mdef%nmat(i)))
       allocate(tmp_nam(mdef%nmat(i)))

       iat=0

       rpts:do 
          
          if(iat .ge. mdef%nmat(i) ) then
             exit
          end if
          
          iat=iat+1
          read(un_fld,*,iostat=ierr) tmp_nam(iat), tmp_mass(iat), &
               tmp_chrg(iat), irpt
          if(ierr .ne. 0) call error(un_log,5,1)

          !Check to see if atoms in the molecule are for CV analysis
          do j=1,numconn
             if(upcase(tmp_nam(iat)) .eq. cnatnam(j)) then
                molcntyp(i)=1
             end if
          end do
          
          if( irpt .gt. 1 ) then
             tmp_nam(iat+1:iat+irpt-1)=tmp_nam(iat)
             tmp_mass(iat+1:iat+irpt-1)=tmp_mass(iat)
             tmp_chrg(iat+1:iat+irpt-1)=tmp_chrg(iat)
             iat=iat+irpt-1
          end if
          
       end do rpts
       

       !Construct alst entries for this molecule type
       
       jprev=sum(mdef%nmol(1:i-1))

       do j=1,mdef%nmol(i)
          do k=1,mdef%nmat(i)

             iatom=iatom+1
             atlst(iatom)%mtyp=i
             atlst(iatom)%imol=jprev+j
             atlst(iatom)%mass=tmp_mass(k)
             atlst(iatom)%chrg=tmp_chrg(k)

          end do
       end do

       deallocate(tmp_mass) 
       deallocate(tmp_chrg)
       deallocate(tmp_nam)
       cont=.true.
       
       !Locate the bonding lists

       do
          read(un_fld,*,iostat=ierr) line
          if(ierr .ne. 0) call error(un_log,5,1)
          line = adjustl(line)
          line = upcase(line)

          if(line(1:4) .eq. "BOND") exit
          if(line(1:4) .eq. "CONS") exit
          if(line(1:4) .eq. "FINI") then
             cont=.false.
             exit
          end if
       end do

       if(cont) then
          
          !Get the bond information
          backspace(un_fld,iostat=ierr)
          if(ierr .ne. 0) call error(un_log,5,1)
          read(un_fld,*,iostat=ierr) line, mdef%nbnd(i)
          if(ierr .ne. 0) call error(un_log,5,1)

          line = adjustl(line)
          line = upcase(line)
             
          !Assume the largest bondlist is in first molecule
          
          if(mxbnd .eq. 0) then
             mxbnd=mdef%nbnd(i)
             
             if(mxbnd .gt. 0) then
                allocate(mdef%lbnd(2,mxbnd,mdef%nmtp))
                mdef%lbnd=0
             endif
             
          else
             if(mdef%nbnd(i) .gt. mxbnd) then
                
                !Reallocate storage when the assumption is wrong
                mxbnd_old=mxbnd
                mxbnd=mdef%nbnd(i)
                
                allocate(tmp_lbnd(2,mxbnd,mdef%nmtp))
                tmp_lbnd = 0
                tmp_lbnd(:,1:mxbnd_old,1:i-1)=mdef%lbnd(:,:,1:i-1)
                deallocate(mdef%lbnd)
                
                allocate(mdef%lbnd(2,mxbnd,mdef%nmtp))
                mdef%lbnd=tmp_lbnd
                deallocate(tmp_lbnd)
                
             endif
          endif

          select case(line(1:4))
          case( "BOND" )
             formspec = "(4x,2i5)"
          case( "CONS" )
             formspec = "(2i5)"
          end select

          do j=1,mdef%nbnd(i)
             read(un_fld,fmt=formspec,iostat=ierr) mdef%lbnd(:,j,i)
             if(ierr .ne. 0) call error(un_log,5,1)
          end do

         !Read to finish
          do
             read(un_fld,*,iostat=ierr) line
             line =adjustl(line)
             line = upcase(line)
             if(line .eq. "FINISH") exit
          end do
          
       end if !cont - i.e. read bond info
       
    end do


!---Close the FIELD file
    call closech(un_fld)


!---Evaluate the total number of molecules
    numol=sum(mdef%nmol)


!---Construct molecule list arrays
    allocate(mlst(numol))
    mlst%mtyp=0
    mlst%imtp=0
    mlst%n1=0

    im=0
    nfirst=1

    do i=1,mdef%nmtp
       do j = 1,mdef%nmol(i)

          im=im+1
          mlst(im)%mtyp=i
          mlst(im)%imtp=j
          mlst(im)%n1=nfirst
          nfirst=nfirst+mdef%nmat(i)

       end do

    end do

!---Allocate and construct mdef atom-lists

    nmax=maxval(mdef%nmat)
    allocate(mdef%imat(nmax,mdef%nmtp))
    mdef%imat=0
    do i=1,mdef%nmtp
       mdef%imat(1:mdef%nmat(i),i) = (/ (j, j=1,mdef%nmat(i)) /)
    end do

!---Check the number of atoms

    natm_f=dot_product(mdef%nmol,mdef%nmat)
    if(natm_f .ne. natms) call error(un_log,13,1)

!---Report to log
    nmcn=0
    mdef%nclmtp=sum(molcntyp(:))
    allocate(mdef%cntyp(mdef%nclmtp))
    allocate(mdef%cnnam(mdef%nclmtp))
    mdef%cntyp=0
    mdef%cnnam=""
   
    write(un_log,'(a,i7)') "Total number of molecule types:", mdef%nmtp
    write(un_log,'(a,i7)') "Molecule types containing CV atoms:", mdef%nclmtp
    j=0
    do i=1, mdef%nmtp
       if(molcntyp(i) .eq. 1) then
          j=j+1
          mdef%cntyp(j) = i
          mdef%cnnam(j) = mdef%mname(i)
          nmcn=nmcn+mdef%nmol(i)
       end if
    end do
    
    write(un_log,'(a,i7)') "Total number of molecules with CV atoms:", nmcn
    write(un_log,'("  ",4a15)') mdef%cnnam(:)

!---Close the FIELD

    deallocate(molcntyp)
    call closech(un_fld)

    
  end subroutine fieldscan
  
  
  
  
end module readfld
