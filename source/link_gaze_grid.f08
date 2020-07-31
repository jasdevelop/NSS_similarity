!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   gaze_grid_link.f08
! Author: jorge
!
! Created on August 11, 2019, 6:50 PM
!

MODULE link_gaze_grid
    use input_module    
    use fft_convolve_3d
    implicit none
    private
    
    public setup_mask_values, setup_observation_points
    
    contains
    subroutine setup_mask_values(user_list, nu, ti, tf )        
        integer ::  tp ! timepoint to add
        integer :: ti, tf , tc !
        integer :: p(3)  !store the point to add 
        integer :: u , nu  ! user index and number of users in user_list 
        logical :: user_list(nu)   ! indicates which users are toe be included (True) or excluded (False)
        double precision :: temp_val  !store the value here momentarily to asses if its NaN
        double precision :: amplitude  !amplitude of the Impulse added to the mask
        integer :: status
        integer :: k
        integer :: count_valid
        count_valid=0
        tc=(ti+tf)/2   !center of the time window
        do u=1,nu
            if ( user_list(u) .eqv. .false. ) cycle  
            do tp=ti,tf           
                p=0  !reset
                do k=1,3
                    !convert coordinates to gridpoints
                    temp_val= (gaze_data(tp,u,k)-box_bounds(k,1) )/box_spacing(k) 
                    !check for NaN
                    if ( isnan(temp_val) .eqv. .TRUE. ) then
                        !give a voided value, it will be rejected later upon bounds checking
                        p(k)=-1
                    else
                        !store the value of the coordinate in points into p
                        p(k)=NINT( temp_val  )+1
                    end if
                    
                end do    
                
                !We now check if planar mode is active
                if (planar_mode .eqv. .TRUE.) then 
                    !Override the location in time and force it onto the plane
                    !amplitude=exp( -0.5*( tp-tc   )**2  /sigma_in_points(3)**2   )
                    p(3)=1  ! ensures the Gaussian is placed on the grid
                    
                end if
                
                call mask_add_point(p, 1.d0, status)
                
                if (status .eq. 0 ) then
                    count_valid=count_valid+1
                end if
            end do
            
        end do
        !write(*,*) "Added points: ", count_valid
    end subroutine setup_mask_values
    
    subroutine setup_observation_points( obs_array , nobs, nobs_valid, u, ti, tf )        
        integer, intent(in) :: nobs
        integer, intent(inout) :: obs_array(nobs,3)
        integer , intent(out)  :: nobs_valid
        
        integer, intent(in) :: u       !  user_index 
        integer, intent(in) :: ti, tf   !  timepoint at start and end of box
        integer :: count_valid
        integer :: p(3)
        integer :: k , tp   ! indexes for iterating over coordinate and timepoint
        double precision :: temp_val
        integer :: tc  !center of the window, in points
        
        !Returns the valid values, converted into gridpoints, in the obs_array 
        
        tc=(ti+tf)/2
        count_valid=0
        do tp=ti,tf           
            p=0  !reset
            do k=1,3
                !convert coordinates to gridpoints
                temp_val= (gaze_data(tp,u,k)-box_bounds(k,1) )/box_spacing(k) 
                !check for NaN                 
                if ( isnan(temp_val) .eqv. .TRUE. ) then
                    ! if NAN set the result to -1 which is void, it will be rejected later 
                    p(k)=-1
                else 
                    !store the value of the coordinate in points into p
                    p(k)=NINT( temp_val  )+1
                end if
                !write(*,*) "pk" , p(k), k, box_shape(k)
                !Override if it is out of bounds
                if ((p(k) .lt. 1) .or. ( p(k) .gt. box_shape(k) ) ) then
                    !write(*,*) "Voided ", p(k) , k
                    p(k)=-1
                    
                end if
            end do         
            
            if (planar_mode .eqv. .TRUE. ) then
                !Override the third dimension if doing planar, perhaps overturning the bound checkers veredict for 3rd dim
                p(3)=1
            end if
            
            !Test the array for any bad values, marked as -1
            !Store the value if clean for NANs 
            if ( ANY(p==-1 ) .eqv. .FALSE. ) then
                count_valid=count_valid+1
                obs_array( count_valid,: ) =p(:)
                                
            else
                !write(*,*) "Discarded ", p(1), p(2), p(3)
                continue
            end if

        end do

       !Export the number of valid observations found    
       nobs_valid=count_valid 
       !write(*,*)  "Valid observation points ", nobs_valid , nobs 
    end subroutine setup_observation_points
    

END MODULE link_gaze_grid
