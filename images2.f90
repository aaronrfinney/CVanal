module images
  use util
  
  implicit none
  
contains

  subroutine minimg(imcon,natms,cell,xxx,yyy,zzz)
      
!=======================================================================
! calculates the minimum image for atom pairs
!
! adapted from the dl_poly routine images
!
! imcon=0 no boundary conditions apply
! imcon=1 standard cubic boundaries apply
! imcon=2 orthorhombic boundaries apply
! imcon=3 parallelepiped boundaries apply
! imcon=4 truncated octahedron boundaries apply
! imcon=5 rhombic dodecahedron boundaries apply
! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
!
! note: in all cases the centre of the cell is at (0,0,0)
!=======================================================================

    integer :: imcon,natms
    real(dp),dimension(*) :: xxx,yyy,zzz
    real(dp),dimension(9) :: cell,rcell

    intent(in) :: imcon,natms,cell
    intent(inout) :: xxx,yyy,zzz

    integer :: i
    real(dp) :: aaa,bbb,ccc,rt2,det
    real(dp) :: ssx,ssy,ssz,xss,yss,zss
    
    ! Calculate according to image convention
      
    if(imcon.eq.1) then 

       ! IMCON=1 Standard cubic boundary conditions
       
       aaa=1.0_dp/cell(1)
       
       do i=1,natms
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
       end do
       
    else if(imcon.eq.2) then
       
       ! IMCON=2 Rectangular (slab) boundary conditions
       
       aaa=1.0_dp/cell(1)
       bbb=1.0_dp/cell(5)
       ccc=1.0_dp/cell(9)
       
       do i=1,natms
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(5)*nint(bbb*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(ccc*zzz(i))
          
       end do
        
    else if(imcon.eq.3) then

       ! IMCON=3 Parallelepiped boundary conditions
        
       call invert(cell,rcell,det)
        
       do i=1,natms
          
          ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
          ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
          ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))
          
          xss=ssx-nint(ssx)
          yss=ssy-nint(ssy)
          zss=ssz-nint(ssz)
          
          xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
       end do
        
      else if(imcon.eq.4) then

         ! IMCON=4 Truncated octahedral boundary conditions
        
         aaa=1.0_dp/cell(1)
        
         do i=1,natms
          
            xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
            yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
            zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
          
            if((abs(xxx(i))+abs(yyy(i))+abs(zzz(i))).ge. &
                 (0.75_dp*cell(1)))then
            
               xxx(i)=xxx(i)-0.5_dp*sign(cell(1),xxx(i))
               yyy(i)=yyy(i)-0.5_dp*sign(cell(1),yyy(i))
               zzz(i)=zzz(i)-0.5_dp*sign(cell(1),zzz(i))
               
            endif
          
         end do
        
      else if(imcon.eq.5) then

         ! IMCON=5 Rhombic dodecahedral boundary conditions
        
         rt2=sqrt(2.0_dp)
         
         aaa=1.0_dp/cell(1)
         bbb=1.0_dp/cell(9)
        
         do i=1,natms
          
            xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
            yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
            zzz(i)=zzz(i)-cell(9)*nint(bbb*zzz(i))
          
            if((abs(xxx(i))+abs(yyy(i))+abs(rt2*zzz(i))).ge. &
                 cell(1)) then
            
               xxx(i)=xxx(i)-0.5_dp*sign(cell(1),xxx(i))
               yyy(i)=yyy(i)-0.5_dp*sign(cell(1),yyy(i))
               zzz(i)=zzz(i)-0.5_dp*sign(cell(9),zzz(i))
            
            endif
          
         end do
        
      else if(imcon.eq.6) then

         ! IMCON=6 x-y boundary conditions 

         det = cell(1)*cell(5) - cell(2)*cell(4)

         det = 1.0_dp/det

         rcell(1) =  det*cell(5)
         rcell(2) = -det*cell(2)
         rcell(4) = -det*cell(4)
         rcell(5) =  det*cell(1)
        
         do i=1,natms

            ssx = rcell(1)*xxx(i) + rcell(4)*yyy(i)
            ssy = rcell(2)*xxx(i) + rcell(5)*yyy(i)

            xss = ssx - nint(ssx)
            yss = ssy - nint(ssy)

            xxx(i)=cell(1)*xss + cell(4)*yss
            yyy(i)=cell(2)*xss + cell(5)*yss

         end do

      end if

    end subroutine minimg

  end module images

