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
  use dstrb_smooth

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
  integer :: j,k,l,aj,ak,al,ncon,nox,nca,ncc,nc,now
  integer :: kimol,limol,cain,pres
  integer,parameter :: distbins=50,angbins=32,coordbin=15
  real(dp),parameter :: distbindelta=0.2,angbindelta=0.1,smfact=2.23548
  real(dp) :: mandistara,mandistcal,manangara,manangcal,mandistcaara,mandistcacal
  real(dp),dimension(3) :: rd,vk,vl
  real(dp) :: dist,dvk,dvl,theta,frac
  real(dp) :: var_distara,var_distcal,var_angara,var_angcal,var_distcaara,var_distcacal
  integer :: molcount,imolcount,distbin
  integer,dimension(10) :: fscc
  integer,dimension(20) :: tmpimol
  real(dp),dimension(10) :: dentdstrbn
  real(dp),dimension(coordbin) :: coordox,coordcc,coordow
  real(dp),dimension(distbins) :: distfscc,distararef,distcalref,distara,distcal
  real(dp),dimension(distbins) :: distfscac,distcaararef,distcacalref,distcaara,distcacal
  real(dp),dimension(angbins) :: angfscc,angara,angcal,angararef,angcalref
  real(dp),dimension(:,:),allocatable :: avdistfscc,avangfscc,avdistfscac
  real(dp),dimension(:,:,:),allocatable :: avMD
!
!
!
!
!******************************************************************************

  ird=0
  avncl=0
  coordox = 0.0_dp
  coordcc = 0.0_dp
  coordow = 0.0_dp
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

    do j=1, nacn !Loop from 1 to the total number of CV atoms you have identified
        aj = cvlst(j)%ain ! get the atomic index of this CV atom
        if(cvlst(j)%atnm .eq. "CA      ") then ! if this is a calcium
          nca = nca + 1 ! count the calcium
        end if
    end do

    allocate(cafs(nca))
    allocate(avdistfscc(distbins,nca))
    allocate(avdistfscac(distbins,nca))
    allocate(avangfscc(angbins,nca))
    allocate(avMD(nca,7,2))

!---Read histogram data of Calcite and Aragonite
    open(unit=1,file="dist-cc-araref.dat")
    read(1,*)distararef
    open(unit=2,file="dist-cc-calref.dat")
    read(2,*)distcalref
    open(unit=10,file="angararef.dat")
    read(10,*)angararef
    open(unit=11,file="angcalref.dat")
    read(11,*)angcalref
    open(unit=12,file="dist-cac-araref.dat")
    read(12,*)distcaararef
    open(unit=13,file="dist-cac-calref.dat")
    read(13,*)distcacalref

    open(unit=3,file="CV-distances.log")
    write(3,*) "# Central atom index, 6-coord flag, x, y, z, Aragonite C-C dist. MD, &
         &Calcite C-C dist. MD, Aragonite C-Ca-C ang. MD, &
         & Calcite C-Ca-C ang. MD, Aragonite Ca-C dist. MD, Calcite Ca-C dist. MD"

  !******************************************************************************
  
  avdistfscc = 0.0_dp
  avangfscc = 0.0_dp
  avdistfscac = 0.0_dp
  avMD = 0.0_dp

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
     cafs(:)%ain = 0
     cafs(:)%nmol = 0
     nc = 0
     distbin = 0

     do j=1, nacn !Loop from 1 to the total number of CV atoms you have identified
        
        aj = cvlst(j)%ain ! get the atomic index of this CV atom
        
        if(cvlst(j)%atnm .eq. "CA      ") then !if this is a calcium then...
           
           nca = nca + 1 ! count this as a Ca using the nca counter
           cafs(nca)%ain = aj ! record this atom index
           cafs(nca)%imol(:) = 0
           ncon = conn(aj)%ncn  ! get the number of atom connections to this atom
           
           ncc = 0
           nox = 0
           tmpimol = 0
           
           do k=1, ncon !Loop through the connections

              ak = conn(aj)%ain(k) ! get the atomic index of this connected atom
              
              if(atlst(ak)%atnm .eq. "OC      ") then ! if this is an oxygen in carbonate then ...
                 
                 nox = nox + 1 ! count the number of oxygens
                 tmpimol(nox) = atlst(ak)%imol ! get the molecular number of the oxygens counted

             else if(atlst(ak)%atnm .eq. "OW      ") then ! if this is an oxygen in water then ...
                 
                 if(conn(aj)%dcn(k) .lt. 3.0_dp) then ! if the oxygen is within 3A of the Calcium
                    
                    now = now + 1 ! count the number of oxygens
                    
                 end if
              end if
           end do

           ! Check for monodentate vs. bidentate vs. tridentate
           molcount = count(mask = tmpimol(:) .gt. 0) ! count the occurance of carbonate molecules

           do k=1,molcount

              if(tmpimol(k) .gt. 0) then

                 imolcount = count(mask = tmpimol(:) .eq. tmpimol(k)) ! check the denticity

                 dentdstrbn(imolcount) = dentdstrbn(imolcount) + 1.0_dp ! count all denticity

                 nc = nc + 1
                 ncc = ncc + 1

                 ! Count the number of carbonates in the first shell of this Ca
                 cafs(nca)%nmol = ncc
                 cafs(nca)%imol(ncc) = tmpimol(k)

              end if

              where(tmpimol .eq. tmpimol(k)) tmpimol=0
           end do

           ! Record the coordination numbers
           coordox(nox+1) = coordox(nox+1) + 1.0_dp
           coordow(now+1) = coordow(now+1) + 1.0_dp
           coordcc(ncc+1) = coordcc(ncc+1) + 1.0_dp


        end if
     end do

     ! Write the data to the log file
     if(rec .eq. 1) then
        write(un_log,'(/,a)') "***** Average number of carbons per calcium"
     end if
     write(un_log,'(f15.3,f6.3,f6.3)') nstp*tstp, real(nc)/real(nca) !now/real(nca)


     write(3,*) "#  Frame", rec
     ! Go through the Ca first shell carbonates
     do j=1, nca

        distfscc = 0.0_dp
        angfscc = 0.0_dp
        distfscac = 0.0_dp
        cain = cafs(j)%ain !Ca atom index
        
        !If this Ca has six carbons....
        if(cafs(j)%nmol .eq. 6)then
           pres = 1
        !Record the atom as 6-coordinate
        avMD(j,1,1) = avMD(j,1,1) + 1.0_dp
        else
           pres = 0
        end if
        !Count the number of times this Ca appears
        avMD(j,1,2) = avMD(j,1,2) + 1.0_dp

        !Loop through the Carbons....
        do k=1, cafs(j)%nmol-1
           
           !Get the index for this carbonate in the shell
           kimol = cafs(j)%imol(k)
           ak = mlst(kimol)%n1 ! the first one in 4 atoms is a carbon
           fscc(k) = ak
           
           !Ca-C(k) distance
           vk=rxyz(:,ak) - rxyz(:,cain)
           vk = minimg(imcon,cell,vk(1),vk(2),vk(3))
           dvk=sqrt(dot_product(vk,vk))
           ! Get the Ca-C distance distribution
           if(k .eq. 1) then
              distbin=int(dvk/distbindelta) + 1 !Ca-C(1)
              if(distbin .gt. distbins) stop "Increase the number of distbins"
              
              distfscac(distbin) = distfscac(distbin) + 1.0_dp
              avdistfscac(distbin,j) = avdistfscac(distbin,j) + 1.0_dp
           end if
           
           
           do l=k+1, cafs(j)%nmol ! loop through its carbon pair
              
              limol = cafs(j)%imol(l)
              al = mlst(limol)%n1
              fscc(l) = al
              
              ! Get the distance between the carbonate centres
              rd = rxyz(:,al) - rxyz(:,ak)
              rd = minimg(imcon,cell,rd(1),rd(2),rd(3))
              dist = sqrt(dot_product(rd,rd)) ! calculate the distance by array dot product
              
              distbin = int(dist/distbindelta) + 1
              if(distbin .gt. distbins) stop "Increase the number of distbins"
              
              distfscc(distbin) = distfscc(distbin) + 1.0_dp
              avdistfscc(distbin,j) = avdistfscc(distbin,j) + 1.0_dp
              
              ! Get the C-Ca-C angles
              vl=rxyz(:,al) - rxyz(:,cain)
              vl = minimg(imcon,cell,vl(1),vl(2),vl(3))
              dvl=sqrt(dot_product(vl,vl))
              
              theta=acos(dot_product(vl,vk)/(dvl*dvk))
              
              distbin=int(theta/angbindelta)+1
              angfscc(distbin)=angfscc(distbin)+1.0_dp
              avangfscc(distbin,j)=avangfscc(distbin,j)+1.0_dp
              
              ! Get the Ca-C distance distribution
              distbin=int(dvl/distbindelta) + 1    !Ca-C(5-6)
              if(distbin .gt. distbins) stop "Increase the number of distbins"
              
              distfscac(distbin) = distfscac(distbin) + 1.0_dp
              avdistfscac(distbin,j) = avdistfscac(distbin,j) + 1.0_dp
              
           end do
        end do
        
        
        call smooth(distfscc,smfact)
        call smooth(angfscc,smfact)
        call smooth(distfscac,smfact)
        
        !pres = 1
        
        !====Manhatten distance
        distara(:) = distfscc(:)/sum(distfscc) - distararef(:)
        distcal(:) = distfscc(:)/sum(distfscc) - distcalref(:)
        angara(:) = angfscc(:)/sum(angfscc) - angararef(:)
        angcal(:) = angfscc(:)/sum(angfscc) - angcalref(:)
        distcaara(:) = distfscac(:)/sum(distfscac) - distcaararef(:)
        distcacal(:) = distfscac(:)/sum(distfscac) - distcacalref(:)
        
        !           if(j .eq. 1) then
        !              do k=1,distbins
        !                 print *, distfscac(k)/sum(distfscac),distcaararef(k),distcaara(k)
        !              end do
        !           end if
        
        mandistara = sum(abs(distara))
        mandistcal = sum(abs(distcal))
        manangara = sum(abs(angara))
        manangcal = sum(abs(angcal))
        mandistcaara = sum(abs(distcaara))
        mandistcacal = sum(abs(distcacal))
        
        !====Euclidean distance
        !distara(:) = (distfscc(:)/sum(distfscc) - distararef(:))**2
        !distcal(:) = (distfscc(:)/sum(distfscc) - distcalref(:))**2
        !angara(:) = (angfscc(:)/sum(angfscc) - angararef(:))**2
        !angcal(:) = (angfscc(:)/sum(angfscc) - angcalref(:))**2
        
        !mandistara = sum(abs(distara))
        !mandistcal = sum(abs(distcal))
        !manangara = sum(abs(angara))
        !manangcal = sum(abs(angcal))
        
        !====Canberra distance
        !mandistara = sum(distfscc(:))
        !mandistcal = sum(angfscc(:))
        !distfscc(:)=distfscc(:)/mandistara
        !angfscc(:)=angfscc(:)/mandistcal
        
        !distfscc(:)=distfscc(:)+0.00000001_dp
        !angfscc(:)=angfscc(:)+0.00000001_dp
        !if(rec .eq. 1 .and. j .eq. 1) then
        !   distararef(:)=distararef(:)+0.00000001_dp
        !   distcalref(:)=distcalref(:)+0.00000001_dp
        !   angararef(:)=angararef(:)+0.00000001_dp
        !   angcalref(:)=angcalref(:)+0.00000001_dp
        !end if
        
        !distara(:)=abs(distfscc(:)-distararef(:))/(abs(distfscc(:))+abs(distararef(:)))
        !distcal(:)=abs(distfscc(:)-distcalref(:))/(abs(distfscc(:))+abs(distcalref(:)))
        !angara(:)=abs(angfscc(:)-angararef(:))/(abs(angfscc(:))+abs(angararef(:)))
        !angcal(:)=abs(angfscc(:)-angcalref(:))/(abs(angfscc(:))+abs(angcalref(:)))
        
        !mandistara = sum(abs(distara))
        !mandistcal = sum(abs(distcal))
        !manangara = sum(abs(angara))
        !manangcal = sum(abs(angcal))
        
        write(3,'(i6,i3,3f10.4,6f10.5)') cain,pres,rxyz(:,cain),mandistara,mandistcal,&
             manangara,manangcal,mandistcaara,mandistcacal
        
        !Record the sum and sum squared variables
        avMD(j,2,1) = avMD(j,2,1) + mandistara
        avMD(j,3,1) = avMD(j,3,1) + mandistcal
        avMD(j,4,1) = avMD(j,4,1) + manangara
        avMD(j,5,1) = avMD(j,5,1) + manangcal
        avMD(j,6,1) = avMD(j,6,1) + mandistcaara
        avMD(j,7,1) = avMD(j,7,1) + mandistcacal
        
        avMD(j,2,2) = avMD(j,2,2) + mandistara*mandistara
        avMD(j,3,2) = avMD(j,3,2) + mandistcal*mandistcal
        avMD(j,4,2) = avMD(j,4,2) + manangara*manangara
        avMD(j,5,2) = avMD(j,5,2) + manangcal*manangcal
        avMD(j,6,2) = avMD(j,6,2) + mandistcaara*mandistcaara
        avMD(j,7,2) = avMD(j,7,2) + mandistcacal*mandistcacal
        
        
        !else
        
        !  pres=0
        !  write(3,'(i6,i3,3f10.4)') cain,pres,rxyz(:,cain)
        
        !end if
        
        !====Manhatten output
        
        
        !====Euclidean output
        !write(3,*) j-1,sqrt(mandistara),sqrt(mandistcal),sqrt(manangara),sqrt(manangcal),fscc-1
        
        !====Canberra output
        !write(3,*) j-1,mandistara,mandistcal,manangara,manangcal,fscc-1
        
     end do
     
!******************************************************************************
     
     
     
     
  end do traj ! main trajectory analysis loop


!******************************************************************************
!
  write(un_log,'(/,a)') "***** Distribution of CA-OC CA-OW and CA-CC Coordination"
  do i=1,coordbin
     write(un_log,'(i5,3f11.6)') i-1,coordox(i)/sum(coordox),coordow(i)/sum(coordow),coordcc(i)/sum(coordcc)
  end do
  
  write(un_log,'(/,a)') "***** Distribution of Carbonate denticity"
  do i=1,3
     write(un_log,'(i5,f11.6)') i,dentdstrbn(i)/sum(dentdstrbn)
  end do



  !====Average distance and angle distributions for all Ca
  distfscc(:) = 0.0_dp
  angfscc(:) = 0.0_dp
  distfscac(:) = 0.0_dp

  do i=1,nca
     distfscc(:) = distfscc(:) + avdistfscc(:,i)
     angfscc(:) = angfscc(:) + avangfscc(:,i)
     distfscac(:) = distfscac(:) + avdistfscac(:,i)
  end do

  distfscc(:) = distfscc(:)/sum(distfscc(:))
  angfscc(:) = angfscc(:)/sum(angfscc(:))
  distfscac(:) = distfscac(:)/sum(distfscac(:))

  call smooth(distfscc(:),smfact)
  call smooth(angfscc(:),smfact)
  call smooth(distfscac(:),smfact)

  write(un_log,'(/,a)') "***** Distribution of firstshell CC-CC Distances"
  do i=1,distbins
     write(un_log,'(f8.3,f11.6)') (i-0.5)*distbindelta,distfscc(i)
  end do

  write(un_log,'(/,a)') "***** Distribution of firstshell CC--CA--CC Angles"
  do i=1,angbins
     write(un_log,'(f8.3,f11.6)') (i-0.5)*angbindelta,angfscc(i)
  end do  

  write(un_log,'(/,a)') "***** Distribution of firstshell CA-CC Distances"
  do i=1,distbins
     write(un_log,'(f8.3,f11.6)') (i-0.5)*distbindelta,distfscac(i)
  end do

!====Average Manhatten distance based on average distributions
  write(un_log,'(/,a)') "***** Average Manhatten distances based on average distributions"
  write(un_log,'(a)') "Atom index, trajectory fraction, Aragonite C-C dist. MD, &
       &Calcite C-C dist. MD, Aragonite C-Ca-C ang. MD, &
       & Calcite C-Ca-C ang. MD, Aragonite Ca-C dist. MD, Calcite Ca-C dist MD"
  
  do j=1, nca
     if(avMD(j,1,2) .gt. 0.0_dp) then !Ca has to be at least coordinated to C
        !Smooth the average distributions for each Ca
        call smooth(avdistfscc(:,j),smfact)
        call smooth(avangfscc(:,j),smfact)
        call smooth(avdistfscac(:,j),smfact)
        
        distara(:) = avdistfscc(:,j)/sum(avdistfscc(:,j)) - distararef(:)
        distcal(:) = avdistfscc(:,j)/sum(avdistfscc(:,j)) - distcalref(:)
        angara(:) = avangfscc(:,j)/sum(avangfscc(:,j)) - angararef(:)
        angcal(:) = avangfscc(:,j)/sum(avangfscc(:,j)) - angcalref(:)
        distcaara(:) = avdistfscac(:,j)/sum(avdistfscac(:,j)) - distcaararef(:)
        distcacal(:) = avdistfscac(:,j)/sum(avdistfscac(:,j)) - distcacalref(:)

        mandistara = sum(abs(distara))
        mandistcal = sum(abs(distcal))
        manangara = sum(abs(angara))
        manangcal = sum(abs(angcal))
        mandistcaara = sum(abs(distcaara))
        mandistcacal = sum(abs(distcacal))

        frac = avMD(j,1,1)/real(nsteps)
        
        !====Manhatten output
        write(un_log,'(i6,f8.3,6f16.10)') cafs(j)%ain,frac,mandistara,mandistcal,manangara,manangcal,&
             mandistcaara,mandistcacal
     else
        !write(un_log,'(i6,f8.3,4f16.10)') cafs(j)%ain,frac
     end if
  end do

!  do j=1, nca
!     avangfscc(:,j) = avangfscc(:,j)/sum(avangfscc(:,j))
!  end do
!  do j=1, angbins
!     print *, avangfscc(j,:)
!  end do



!====Average Manhatten distance based on averages of stepwise MD values
  write(un_log,'(/,a)') "***** Average Manhatten distances based on average MD values per timestep"
  write(un_log,'(a)') "Atom index, trajectory fraction, Aragonite C-C dist. MD, Standard Dev., &
       &Calcite C-C dist. MD, Standard Dev. "
  
  !C-C Distances
  do j=1, nca
     if(avMD(j,1,2) .gt. 1.0_dp) then !Ca has to be at least coordinated to C

        !Mean values
        !mean = sumx / N
        mandistara = avMD(j,2,1) / avMD(j,1,2)
        mandistcal = avMD(j,3,1) / avMD(j,1,2)

        !Standard deviations
        !sd = (sumx^2 - (sumx*sumx)/N) / (N-1)
        var_distara = (avMD(j,2,2) - (avMD(j,2,1)*avMD(j,2,1))/avMD(j,1,2))/(avMD(j,1,2)-1)
        var_distcal = (avMD(j,3,2) - (avMD(j,3,1)*avMD(j,3,1))/avMD(j,1,2))/(avMD(j,1,2)-1)
        
        frac = avMD(j,1,1)/real(nsteps)

        !====Manhatten output
        write(un_log,'(i6,f8.3,4f16.10)') cafs(j)%ain,frac,mandistara,sqrt(var_distara),&
             mandistcal,sqrt(var_distcal)

     else

        !write(un_log,'(i6,f8.3)') cafs(j)%ain,frac

     end if
  end do

  !Angles
  write(un_log,'(/,a)') " "
  write(un_log,'(a)') "Atom index, trajectory fraction, Aragonite C-Ca-C ang. MD, &
       & Standard Dev., Calcite C-Ca-C ang. MD, Standard Dev."

  do j=1, nca
     if(avMD(j,1,2) .gt. 1.0_dp) then !Ca has to be at least coordinated to C

        !Mean values
        !mean = sumx / N
        manangara = avMD(j,4,1) / avMD(j,1,2)
        manangcal = avMD(j,5,1) / avMD(j,1,2)

        !Standard deviations
        !sd = (sumx^2 - (sumx*sumx)/N) / (N-1)
        var_angara = (avMD(j,4,2) - (avMD(j,4,1)*avMD(j,4,1))/avMD(j,1,2))/(avMD(j,1,2)-1)
        var_angcal = (avMD(j,5,2) - (avMD(j,5,1)*avMD(j,5,1))/avMD(j,1,2))/(avMD(j,1,2)-1)
        
        
        frac = avMD(j,1,1)/real(nsteps)

        !====Manhatten output
        write(un_log,'(i6,f8.3,4f16.10)') cafs(j)%ain,frac,manangara,sqrt(var_angara),&
             manangcal,sqrt(var_angcal)

     else

        !write(un_log,'(i6,f8.3)') cafs(j)%ain,frac

     end if
  end do

  write(un_log,'(/,a)') " "
  write(un_log,'(a)') "Atom index, trajectory fraction, Aragonite Ca-C dist. MD, Standard Dev., &
       &Calcite Ca-C dist. MD, Standard Dev. "
  
  !Ca-C Distances
  do j=1, nca
     if(avMD(j,1,1) .gt. 1.0_dp) then !Ca has to be at least coordinated to C

        !Mean values
        !mean = sumx / N
        mandistcaara = avMD(j,6,1) / avMD(j,1,2)
        mandistcacal = avMD(j,7,1) / avMD(j,1,2)

        !Standard deviations
        !sd = (sumx^2 - (sumx*sumx)/N) / (N-1)
        var_distcaara = (avMD(j,6,2) - (avMD(j,6,1)*avMD(j,6,1))/avMD(j,1,1))/(avMD(j,1,2)-1)
        var_distcacal = (avMD(j,7,2) - (avMD(j,7,1)*avMD(j,7,1))/avMD(j,1,1))/(avMD(j,1,2)-1)
        
        frac = avMD(j,1,1)/real(nsteps)

        !====Manhatten output
        write(un_log,'(i6,f8.3,4f16.10)') cafs(j)%ain,frac,mandistcaara,sqrt(var_distcaara),&
             mandistcacal,sqrt(var_distcacal)

     else

        !write(un_log,'(i6,f8.3)') cafs(j)%ain,frac

     end if
  end do



!
!******************************************************************************


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
