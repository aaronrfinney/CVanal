module dstrb_smooth
  use glob_var
  
  implicit none


  contains

    subroutine gaussian(fwhm,nkb,kernel)

!===============================================================================
!
! Calculation of the normalised discrete Gaussian for smoothing. Provides a
! real array, kernel, of dimension nkb, which is half of a symmetrical,
! normal distribution.
!
!===============================================================================


      integer :: nkb
      real(dp) :: fwhm
      real(dp), dimension(nkb) :: kernel

      intent(in) :: nkb, fwhm
      intent(out) :: kernel

      real(dp) :: sigma, mu
      real(dp) :: pre_exp, k_sum
   
      integer :: kb
 
      !Standard distribution
      mu = 1.0_dp

      kernel = 0.0_dp



!---P(x;mu,pre_exp)= 1 / (sigma*sqrt(2*pi)) * exp(-(xb-mu) / (2*sigma**2))
!---Where sigma=fwhm/sqrt(8*log(2))

      sigma = fwhm / (sqrt(8.0_dp * log(2.0_dp))) 
      
      pre_exp = 1.0_dp / (sigma * sqrt(2.0_dp*pi))

      k_sum = 0.0_dp

!---Calculate half a Gaussian

      do kb=1, nkb 

       kernel(kb) = pre_exp * exp((-1.0_dp*(real(kb) - mu)**2.0_dp) / &
            (2*sigma*sigma))

       if(kb .eq. 1) then
          k_sum = k_sum + kernel(kb)
       else
          k_sum = k_sum + (2 * kernel(kb))
       end if

    end do   
    
!---Normalise the Gaussian

    kernel = kernel / k_sum

  end subroutine gaussian




  subroutine smooth(fnct,fwhm)

!===============================================================================
!
! Routine which smooths a function using Gaussian kernel regression
!
!===============================================================================

    real(dp), dimension(:) :: fnct
    real(dp) :: fwhm
    intent(inout) :: fnct
    

    !Gaussian parameters
    integer, parameter :: nkb = 10                   !# kernel bins
    real(dp), dimension(nkb) :: kernel       

    
    real(dp), dimension(:), allocatable :: h_calc
    real(dp), dimension(:), allocatable :: s_fnct
    real(dp) :: wght
    integer :: tl
    integer :: lngth
    integer :: i, kb

!---Calculate the Gaussian distribution

    call gaussian(fwhm,nkb,kernel)


!---Set up the arrays

    lngth = (size(fnct))

    allocate(s_fnct(lngth))
    allocate(h_calc(nkb*2-1))

    s_fnct = 0.0_dp


!---Smooth the function sequentially for each bin
    do i=1, lngth

       if(fnct(i) .gt. 0.0) then

          !Zero the temporary array
          h_calc = 0.0_dp          
          wght = 0.0_dp

          !Ensure conservation in the initial and final regions of fnct
          
          !Re-weighting of the Gaussian bins in these regions:
          ! g(bin1) = g(bin1) / (w(1) + 2w(2) + 2w(3) ... +2w(nbins))

          !peak about mu
          h_calc(nkb) = fnct(i) * kernel(1)
          wght = wght + h_calc(nkb)

          !tails outside of mu
          do kb=2, nkb

             tl = kb - 1

             if(i-tl .gt. 0) then

                h_calc(nkb - tl) = fnct(i) * kernel(kb)
                wght = wght + h_calc(nkb - tl)

             end if


             if(i+tl .lt. lngth) then
                
                h_calc(nkb + tl) = fnct(i) * kernel(kb)
                wght = wght + h_calc(nkb + tl)

             end if

          end do

          !Normalise
          h_calc(:) = h_calc(:) * (fnct(i) / wght)

          !Now add the smoothed delta peak to the smoothed function

          s_fnct(i) = s_fnct(i) + h_calc(nkb)

          do kb = 2, nkb

             tl =kb - 1

             if(i-tl .gt. 0) then

                s_fnct(i - tl) = s_fnct(i - tl) + h_calc(nkb - tl)

             end if

             if(i+tl .lt. lngth) then

                s_fnct(i + tl) = s_fnct(i + tl) + h_calc(nkb + tl)

             end if

          end do

       end if

    end do


    !Make the smoothed function equal to the function
    fnct = s_fnct


    deallocate(s_fnct)
    deallocate(h_calc)


  end subroutine smooth


end module dstrb_smooth
