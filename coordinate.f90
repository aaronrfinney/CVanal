module coordinate

    use glob_var
    use images
    use error_close
    use toxyz
    
    implicit none


  contains
    subroutine connect(nacn,atlst,rxyz,cvlst,conn,dcoord,imcon,cell)
      
!===============================================================================
!
! A routine to identify coordinated atoms defined by a cut-off distance
!
!===============================================================================

      integer :: nacn,imcon
      type(atom_list),dimension(:) :: atlst
      type(atom_list),dimension(:) :: cvlst 
      type(connection),dimension(:) :: conn
      real(dp) :: dcoord
      real(dp), dimension(9) :: cell
      real(dp), dimension(:,:) :: rxyz
 
      intent(in) :: nacn,dcoord
      intent(inout) :: atlst,rxyz,cvlst,conn,imcon,cell
      
      integer :: i, j, ai, aj
      real(dp), dimension(3) :: rd
      real(dp) :: r2coord,rd2
      
      logical, parameter :: DEBUG=.false.
      
      

!---Initialise the connection arrays
      conn(:)%ncn=0
      do i=1,maxcoord
         
         conn(:)%ain(i)=0
         conn(:)%iatp(i)=0
         conn(:)%dcn(i)=0.0_dp
         
      end do
 

!---Get the coordination distance cut-off squared
      r2coord=dcoord*dcoord
      
      !---Loop over CV atom coordinates with all other
      !   atom coordinates to identify connectivity


      do i=1, nacn-1
         
         ai = cvlst(i)%ain
         

         do j=i+1, nacn

            aj = cvlst(j)%ain
            
            
            !Calculate the distance between atoms
            rd = rxyz(:,aj) - rxyz(:,ai)
            rd = minimg(imcon,cell,rd(1),rd(2),rd(3))
            
            rd2 = dot_product(rd,rd)
            

            !If connected, log the coordination
            if(rd2 .lt. r2coord) then
               
               conn(ai)%ncn = conn(ai)%ncn + 1
               conn(aj)%ncn = conn(aj)%ncn + 1
               
               if(conn(ai)%ncn .gt. maxcoord) call error(un_log,16,1)
               if(conn(aj)%ncn .gt. maxcoord) call error(un_log,16,1)
               
               conn(ai)%ain(conn(ai)%ncn) = atlst(aj)%ain
               conn(ai)%iatp(conn(ai)%ncn) = atlst(aj)%iatp
               conn(ai)%dcn(conn(ai)%ncn) = sqrt(rd2)
               
               conn(aj)%ain(conn(aj)%ncn) = atlst(ai)%ain
               conn(aj)%iatp(conn(aj)%ncn) = atlst(ai)%iatp
               conn(aj)%dcn(conn(aj)%ncn) = sqrt(rd2)
                           
            end if
         end do !j loop
      end do ! i loop
      
      return      
      
    
    end subroutine connect
    
    
  end module coordinate
  
  
  
  

