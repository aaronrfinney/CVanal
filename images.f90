module images
  use util
  
  implicit none
  
contains

  function minimg(imcon,cell,xxx,yyy,zzz)
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

    real(dp), dimension(3) :: minimg

    integer :: imcon
    real(dp) :: xxx,yyy,zzz
    real(dp),dimension(9) :: cell,rcell

    intent(in) :: imcon,cell
    intent(inout) :: xxx,yyy,zzz

    real(dp) :: aaa,bbb,ccc,rt2,det
    real(dp) :: ssx,ssy,ssz,xss,yss,zss

    minimg = 0.0_dp
    
    ! Calculate according to image convention
      
    if(imcon.eq.1) then 
       
       ! IMCON=1 Standard cubic boundary conditions
       
       aaa=1.0_dp/cell(1)
       
       xxx=xxx-cell(1)*nint(aaa*xxx)
       yyy=yyy-cell(1)*nint(aaa*yyy)
       zzz=zzz-cell(1)*nint(aaa*zzz)
       
       
    else if(imcon.eq.2) then
       
       ! IMCON=2 Rectangular (slab) boundary conditions
       
       aaa=1.0_dp/cell(1)
       bbb=1.0_dp/cell(5)
       ccc=1.0_dp/cell(9)
       
       
       xxx=xxx-cell(1)*nint(aaa*xxx)
       yyy=yyy-cell(5)*nint(bbb*yyy)
       zzz=zzz-cell(9)*nint(ccc*zzz)
       
       
    else if(imcon.eq.3) then
       
       ! IMCON=3 Parallelepiped boundary conditions
       
       call invert(cell,rcell,det)
       
       
       ssx=(rcell(1)*xxx+rcell(4)*yyy+rcell(7)*zzz)
       ssy=(rcell(2)*xxx+rcell(5)*yyy+rcell(8)*zzz)
       ssz=(rcell(3)*xxx+rcell(6)*yyy+rcell(9)*zzz)
       
       xss=ssx-nint(ssx)
       yss=ssy-nint(ssy)
       zss=ssz-nint(ssz)
       
       xxx=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
       yyy=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
       zzz=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
       
       
    else if(imcon.eq.4) then
       
       ! IMCON=4 Truncated octahedral boundary conditions
       
       aaa=1.0_dp/cell(1)

       
       xxx=xxx-cell(1)*nint(aaa*xxx)
       yyy=yyy-cell(1)*nint(aaa*yyy)
       zzz=zzz-cell(1)*nint(aaa*zzz)
       
       if((abs(xxx)+abs(yyy)+abs(zzz)).ge. &
            (0.75_dp*cell(1)))then
          
          xxx=xxx-0.5_dp*sign(cell(1),xxx)
          yyy=yyy-0.5_dp*sign(cell(1),yyy)
          zzz=zzz-0.5_dp*sign(cell(1),zzz)
          
       endif
       
       
       
    else if(imcon.eq.5) then
       
       ! IMCON=5 Rhombic dodecahedral boundary conditions
       
       rt2=sqrt(2.0_dp)
       
       aaa=1.0_dp/cell(1)
       bbb=1.0_dp/cell(9)
       
       
       xxx=xxx-cell(1)*nint(aaa*xxx)
       yyy=yyy-cell(1)*nint(aaa*yyy)
       zzz=zzz-cell(9)*nint(bbb*zzz)
       
       if((abs(xxx)+abs(yyy)+abs(rt2*zzz)).ge. &
            cell(1)) then
          
          xxx=xxx-0.5_dp*sign(cell(1),xxx)
          yyy=yyy-0.5_dp*sign(cell(1),yyy)
          zzz=zzz-0.5_dp*sign(cell(9),zzz)
          
       endif
       
       
    else if(imcon.eq.6) then
       
       ! IMCON=6 x-y boundary conditions 
       
       det = cell(1)*cell(5) - cell(2)*cell(4)
       
       det = 1.0_dp/det
       
       rcell(1) =  det*cell(5)
       rcell(2) = -det*cell(2)
       rcell(4) = -det*cell(4)
       rcell(5) =  det*cell(1)
       
       
       ssx = rcell(1)*xxx + rcell(4)*yyy
       ssy = rcell(2)*xxx + rcell(5)*yyy
       
       xss = ssx - nint(ssx)
       yss = ssy - nint(ssy)
       
       xxx=cell(1)*xss + cell(4)*yss
       yyy=cell(2)*xss + cell(5)*yss
       

    end if

    minimg(1) = xxx
    minimg(2) = yyy
    minimg(3) = zzz
    
  end function minimg
  
end module images

