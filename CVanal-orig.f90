program main

  use glob_var
  use variables
  use channels
  use error_close
  use util
  use readinp
  use readcon
  use readfld
  use readhist
  use cvatdefine
  use coordinate


  implicit none


!===============================================================================
!
! Program to analyse MD simulations and calculate order parameters
!
! A. R. Finney (July 2017)
!
!===============================================================================


  integer :: i,ierr,rec,ird,nsteps
  logical :: skip
  real(dp) :: tau,tzero,tfinal,avncl
  real(dp), save :: otime
  integer :: maxclmol

  real(dp),dimension(:),allocatable :: maxmol

  logical, parameter :: DEBUG=.false.


!******************************************************************************
!
!   PUT YOUR VARIABLES HERE
!
  integer :: j,k,l,aj,ak,al,ncon,nox,nca,ncc
  integer :: kimol,limol
  real(dp),dimension(3) :: rd
  real(dp) :: nc,now,dist,dista,distc
  integer :: bin,molcount,imolcount,distbin
  integer,dimension(6) :: fscc
  integer,dimension(20) :: tmpimol
  real(dp),dimension(1:10) :: dstrbn,dentdstrbn
  real(dp),dimension(50) :: histodist,fsccdist,distara,distcal,rdara,rdcal
!
!
!
!
!******************************************************************************

  ird=0
  avncl=0
  dstrbn = 0.0_dp
  dentdstrbn = 0.0_dp

  if(DEBUG) write(*,*) "STARTING THE APP"

  !---Open output log file
  call opench(un_log,'CVanal.log','replace','write','formatted')

  write(un_log,'(/,a)',iostat=ierr) "******************************&
       &**************************************************"
  if(ierr .ne. 0) then
     write(*,*) &
          'Unable to write to the output log file'
  end if
  write(un_log,'(a)') "**                            &
       &                                                **"
  write(un_log,'(a)') "**               CV ANALYSIS O&
       &F MD TRAJECTORIES -- LOG FILE                   **"
  write(un_log,'(a)') "**                            &
       &                                                **"
  write(un_log,'(a)') "**                            &
       &       CVanal                                   **"
  write(un_log,'(a)') "**                            &
       &                                                **"
  write(un_log,'(a)') "**              Written by A. &
       &R. Finney (finneyar@gmail.com)                  **"
  write(un_log,'(a)') "**                            &
       &Version 1.0 07/2017                             **"
  write(un_log,'(a)') "**                            &
       &                                                **"
  write(un_log,'(a,/)') "******************************&
       &**************************************************"

  write(un_log,'(a)') "Units are those of input data"

  if(DEBUG) WRITE(*,*) "OUTPUT CHANNEL:", un_log

!---Read Canano-anal.inp
  call inputscan(nframes,nomit,nskip,numconn,cnatnam,&
       dcoord,xyzcl,trjfm)


  allocate(maxmol(nframes))
  maxmol=0.0_dp
  maxclmol=0

  if(DEBUG) write (*,*) "READ THE INPUT"


!---Read CONFIG file
  call atcnt(natms)

  allocate(atlst(natms))
  allocate(rxyz(3,natms))
  allocate(clxyz(3,natms))
  allocate(conn(natms))
  rxyz=0.0_dp
  clxyz=0.0_dp


  if(DEBUG) write (*,*) "COUNTED ATOMS IN CONFIG:", natms

  call configscan(natms,nuniq,auniq,atlst%ain,atlst%atnm,atlst%iatp,&
       atlst%iclu,atlst%aicl,nacn,cnatnam)


  if(DEBUG) write (*,*) "READ CONFIG"


!---Read FIELD file
  call fieldscan(mdef,atlst,mlst,numol,natms,numconn,cnatnam,maxclmol)

  if(DEBUG) write(*,*) "READ FIELD"

!---Define atoms for CV analysis
  call cvatdef(natms,atlst,nacn,cvlst)
  if(DEBUG) write(*,*) "DEFINED CV ATOMS"


!******************************************************************************
!---Initialise some variables
  nca = 0
  histodist = 0

    do j=1, nacn !Loop from 1 to the total number of CV atoms you have identified
        aj = cvlst(j)%ain ! get the atomic index of this CV atom
        if(cvlst(j)%atnm .eq. "CA      ") then ! if this is a calcium
          nca = nca + 1 ! count the calcium
        end if
    end do

    allocate(cafs(nca))

!---Read histogram data of Calcite and Aragonite
  open(unit=1,file="distara.dat")
  read(1,*)distara
  open(unit=2,file="distcal.dat")
  read(2,*)distcal
  open(unit=3,file="CV-distances.log")
  write(3,*) "# Central atom, central atom index, Aragonite Manhatten distance, &
       &Calcite Manhatten distance, Indeces for coordinated atoms"
  !******************************************************************************



!===Loop through the frames and analyse....

!---Read in the trajectory data

    traj: do rec=1, nframes
     if(DEBUG) WRITE(*,*) "REC", rec

     if(mod(rec,10) .eq. 0) write (*,'(a,i7)') &
          "Frames analysed: ", rec

     skip=(rec .lt. nframes)

     ! Call the appropriate routine depending on DL version
     if(trjfm%frm) then ! formatted
        if(trjfm%vrsn .eq. 2) then
           call readtraj_f(ird,natms,atlst%atnm,rxyz,&
                nstp,tstp,cell,imcon,nomit,nskip,skip)
        elseif(trjfm%vrsn .eq. 4) then
           call readtraj_dl4_f(ird,natms,atlst%atnm,rxyz,&
                nstp,tstp,cell,imcon,nomit,nskip,skip)
        end if
     else
        if(trjfm%vrsn .eq. 2) then ! unformatted
           call readtraj_u(ird,natms,atlst%atnm,rxyz,&
                nstp,tstp,cell,imcon,nomit,nskip,skip)
        elseif(trjfm%vrsn .eq. 4) then
           call readtraj_dlp4_u(ird,natms,atlst%atnm,rxyz,&
                nstp,tstp,cell,imcon,nomit,nskip,skip)
        end if
     end if

     nsteps=rec

     if(ird .lt. 0) then
        nsteps=rec-1
        exit traj
     end if

     if(DEBUG) write(*,*) " Read frame", nstp*tstp

     !Record time zero
     if(rec .eq. 1) then
        tzero=nstp*tstp
        otime=tzero
        tau=tzero
     end if
     if(rec .eq. 2) then
        tau = nstp*tstp - otime
     end if

     !Get av cell
     do i=1,9
        avcell(i)=avcell(i)+cell(i)
     end do


     !Create the connection lists for CV atoms
     call connect(nacn,atlst,rxyz,cvlst,conn,dcoord,imcon,cell)
     if(DEBUG) write(*,*) " Connected atoms"




!******************************************************************************
!
! This is the place to do your CV analysis on each frame of the trajectory file.
! Look in the file variables.f90 to check how to reference data with
! derived types.
! For example, if you want to calculate the average number of carbons around
! a calcium in each frame, you could do the following...


     !Initialise some variables
     nca = 0
     now = 0
     cafs(:)%nmol=0
     nc = 0.0_dp
     distbin = 0

     do j=1, nacn !Loop from 1 to the total number of CV atoms you have identified

        aj = cvlst(j)%ain ! get the atomic index of this CV atom

        if(cvlst(j)%atnm .eq. "CA      ") then !if this is a calcium then...

           nca = nca + 1 ! count this as a Ca using the nca counter
           cafs(nca)%imol(:) = 0 ! initialise the array of calcium firstshell
           ncon = conn(aj)%ncn  ! get the number of atom connections to this atom

           ncc = 0
           nox = 0
           tmpimol = 0

           do k=1, ncon !Loop through the connections

              ak = conn(aj)%ain(k) ! get the atomic index of this connected atom

              if(atlst(ak)%atnm .eq. "OC      ") then ! if this is an oxygen in carbonate then ...

                 nox = nox + 1 ! count the number of oxygens
                 tmpimol(nox) = atlst(ak)%imol ! get the molecular number of the oxygens counted

                 !***** print ak and check the coordination in VMD
              else if(atlst(ak)%atnm .eq. "OW      ") then ! if this is an oxygen in water then ...

                 if(conn(aj)%dcn(k) .lt. 3.0_dp) then ! if the oxygen is within 3A of the Calcium

                    now = now + 1.0_dp ! count the number of oxygens

                 end if
              end if
           end do

           ! Check for monodentate vs. bidentate vs. tridentate
           molcount = count(mask = tmpimol(:) .gt. 0) ! count the occurance of carbonate molecules
           !***** print *, "molcount", molcount

           do k=1,molcount
              !***** print *, tmpimol(1:molcount)

              if(tmpimol(k) .gt. 0) then

                 imolcount = count(mask = tmpimol(:) .eq. tmpimol(k)) ! check the denticity
                 !***** print *, k, imolcount

                 dentdstrbn(imolcount) = dentdstrbn(imolcount) + 1.0_dp ! count all denticity

                 nc = nc + 1.0_dp
                 ncc = ncc + 1

                 !***** print *, ncc

                 ! Count the number of carbonates in the first shell of this Ca
                 cafs(nca)%nmol = ncc
                 cafs(nca)%imol(ncc) = tmpimol(k)

              end if

              !***** print *, cafs(nca)%nmol
              !***** print *, cafs(nca)%imol(:)

              where(tmpimol .eq. tmpimol(k)) tmpimol=0
           end do

           ! count the CA-OC coordination distribution
           bin = ncc + 1
           dstrbn(bin) = dstrbn(bin) + 1.0_dp

        end if
     end do

     !print *, ncc, nca
     ! Write the data to the log file
     if(rec .eq. 1) then
        write(un_log,'(/,a)') "***** Average number of carbons per calcium"
     end if
     write(un_log,'(f15.3,f6.3,f6.3)') nstp*tstp, real(nc)/real(nca) !now/real(nca)

     ! GO through the Ca first shell carbponates
     do j=1, nca

        !print *, cafs(j)%nmol
        fsccdist = 0

        if(cafs(j)%nmol .eq. 6)then

           !print *, "***** Index of Calcium", j-1
           do k=1, cafs(j)%nmol-1

           !Get the index for this carbonate in the shell
           kimol = cafs(j)%imol(k)
           ak = mlst(kimol)%n1 ! the first one in 4 atoms is a carbon
           fscc(k) = ak

              do l=k+1, cafs(j)%nmol ! loop through other carbons

              limol = cafs(j)%imol(l)
              al = mlst(limol)%n1
              fscc(l) = al

              ! Get the distance between the carbnate centres
              rd = rxyz(:,al) - rxyz(:,ak)
              rd = minimg(imcon,cell,rd(1),rd(2),rd(3))

              dist = sqrt(dot_product(rd,rd)) ! calculate the distance by array dot product

              distbin = int(dist*5) ! set bin for distance distribution to 0.1A
              histodist(distbin) = histodist(distbin) + 1
              fsccdist(distbin) = fsccdist(distbin) + 1

              !print *,ak,al,dist
              end do
           end do
           !print '(f11.6)', fsccdist/sum(fsccdist)

           rdara = fsccdist/sum(fsccdist) - distara
           rdcal = fsccdist/sum(fsccdist) - distcal

           dista = sum(abs(rdara))
           distc = sum(abs(rdcal))

           write(3,*) "Calcium",j-1,dista,distc,fscc-1

        end if
     end do




!******************************************************************************




  end do traj ! main trajectory analysis loop



  write(un_log,'(/,a)') "***** Distribution of CA-OC Coordination"
  do i=1,10
     write(un_log,'(i5,f11.6)') i,dstrbn(i)/(nca*nsteps)
  end do
  
  write(un_log,'(/,a)') "***** Distribution of Carbonate denticity"
  do i=1,10
     write(un_log,'(i5,f11.6)') i,dentdstrbn(i)/sum(dentdstrbn)
  end do
  
  write(un_log,'(/,a)') "***** Distribution of firstshell CC-CC Distance"
  do i=1,50
     write(un_log,'(f8.3,f11.6)') 0.1*i,histodist(i)/sum(histodist)
  end do

  if(DEBUG) WRITE(*,*) "FINISHED ITERATING THROUGH FRAMES"
  write(un_log,'(/,a)') " "
  write(un_log,'(/,a)') "***** ANALYSIS COMPLETE *****"


!---Report time
  write(un_log,'(/,a)') "Trajectory times (ps)"
  write(un_log,'(a,f15.5)') "t0", tzero
  tfinal=nstp*tstp
  write(un_log,'(a,f15.5)') "tf", tfinal
  write(un_log,'(a,f15.5)') "tau (between frames read)", tau

!---Calculate avcell
  write(un_log,'(/,a9,3f12.3)') "Av. cell:", avcell(1:3)/nsteps
  write(un_log,'(a9,3f12.3)') "         ",avcell(4:6)/nsteps
  write(un_log,'(a9,3f12.3)') "         ",avcell(7:9)/nsteps


end program main
