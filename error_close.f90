module error_close

  implicit none


contains 
  subroutine error(unitno,eid,kid)
    
!===============================================================================
! Subroutine to indicate the cause of an error and to close the 
! analysis package if necessary 
!===============================================================================

    integer, intent(in) :: unitno, eid, kid

    if(eid .gt. 0) then
       write(unitno,'(//,a)') &
            'CVanal has terminated due to error:'

       if(eid .eq. 1) then
          write(*,*) &
               'Unable to write to the output log file'
          
       elseif(eid .eq. 3) then
          write(unitno,'(a)') &
               'CVanal.inp read error'
          
       elseif(eid .eq. 5) then
          write(unitno,'(a)') &
               'Unable to read from FIELD'
          
       elseif(eid .eq. 6) then
          write(unitno,'(a)') &
               'The number of unique atom types is greater&
               & than the maximum allowed'
          write(unitno,'(/,a)') &
               'Action: Increase the size of the variable limit'

       elseif(eid .eq. 7) then
          write(unitno,'(a)') &
               'HISTORY file is inconsistent with CONFIG'

       elseif(eid .eq. 8) then
          write(unitno,'(a)') &
               'Error reading HISTORY'

       elseif(eid .eq. 9) then
          write(unitno,'(a)') &
               'Number of atoms in CONFIG and HISTORY&
               & do not macth'          

       elseif(eid .eq. 10) then
          write(unitno,'(a)') &
               'EOF found during frame read in HISTORY'    

       elseif(eid .eq. 12) then
          write(unitno,'(a)') &
               'Error writing to CVanal.xyz' 

       elseif(eid .eq. 13) then
          write(unitno,'(a)') &
               'Number of atoms in CONFIG and FIELD&
               & do not macth'

       elseif(eid .eq. 15) then
          write(unitno,'(a)') &
               'Number of iterations to reconstruct &
               &bonds across the boundary exceeds &
               &the limit' 

       elseif(eid .eq. 16) then
          write(unitno,'(a)') &
               'Coordination of atoms exceed that of &
               & the maximum.' 
          write(unitno,'(/,a)') &
               'Action: Increase the size of the variable maxcoord'
       end if

       write(*,'(a,i4)') 'CVanal has terminated due to error:',eid

    end if



    !Close the package
    if(kid .gt. 0) then
       
       stop

    end if


  end subroutine error


 

end module error_close
