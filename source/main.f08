!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   main.f08
! Author: jorge
!
! Created on July 28, 2019, 2:14 PM
!

!include 'mkl_vsl.f90'

program nss_calculator
  use fft_convolve_3d
  use input_module
  use link_gaze_grid
  use lm_exact
  use d_mcarlo
  implicit none
  call initialize_random_mkl(777)
  call mc_randomizer_g()  
  call read_parameters()
  call setup_axis_all()
  call read_data()
  call read_groups()
  
  write(*,*) "Done reading data"
  !call grid_calculation()
  
  select case(mode_code)
      case(1) 
          write(*,*) "Executing numerical grid calculation"
          call grid_calculation( box_shape, kernel_shape, sigma_in_points  )
          
      case(2)
          write(*,*) "Executing analytical solution calculator"      
          call exact_calculation(.FALSE.)
      
      case(3)    
          write(*,*) "Executing analytical solution with erf approximation"
          call exact_calculation(.TRUE.)
      case default 
          write(*,*) "Bad mode code"
  end select  
  write(*,*) samp_freq
  write(*,*) axis_x(1), axis_x(2) , axis_x( box_shape(1)  )
  write(*,*) axis_all(1,1) , axis_all(2,1), axis_all(3,1)
  write(*,*) axis_all(1,2) , axis_all(2,2), axis_all(3,2)
  write(*,*) "The end"
  
contains 

subroutine grid_calculation(mask_box , kernel_box, sigma_vals  )
  
  double precision :: time_g1=0.d0, time_g2=0.d0, time_g_tot=0.d0
  integer            :: counter =0
  double precision   :: NSS_score=0.d0  , NSS_single=0.d0
  integer            :: nobs=1
  integer            :: nobs_max 
  integer, allocatable :: obs(:,:)
  integer  :: mask_box(3)  ,   kernel_box(3)  !shapes of the mask and kernel
  double precision :: sigma_vals(3)  !values of sigma UNITS of GRIDPOINTS
  integer :: u  !user being analyzed
  integer :: ti , tf  !initial and final timepoint of window  
  logical, allocatable :: user_selector(:)  !indicates if a user is to be used for the current window
  integer :: g !group number
  integer :: k ! Monte Carlo iteration
  integer :: up  !user prime
  integer :: i
  
  integer, parameter :: fid_NSS_out=16
  character*255       :: filename_NSS_out
  
  !open file where score will be exported
  filename_NSS_out= 'NSS_score_export.dat'
  open(fid_NSS_out, file=trim(filename_NSS_out), action='write' )

  
  !initialize observation-related variables
  nobs_max=points_t  
  allocate( obs(nobs_max,3)   )
  obs=0
  nobs=0
  
  
  !initialize user selection variables
  allocate( user_selector(num_users) )
  user_selector=.TRUE.
  
  
  write(*,*) "Initializing the grid"
 
  
  call initialize_3D( mask_nx=mask_box(1)  , mask_ny=mask_box(2), mask_nz=mask_box(3), & 
        kernel_nx=kernel_box(1) , kernel_ny=kernel_box(2) , kernel_nz=kernel_box(3)  , &
        sigma_x= sigma_vals(1)  ,   sigma_y= sigma_vals(2)  , sigma_z= sigma_vals(3) , planar=planar  )
  
                
  NSS_single=calc_NSS_single()
 
  
  write(*,*) "Sigma vals ", sigma_vals
  
  write(*,*)  "Calculating up to ", num_timepoints
  
  call CPU_time(time_g1)
  
 do ti=1,num_timepoints, window_stride      
     tf=ti+points_t-1   !end of the time interval
     
     !If the window falls out of bounds then don't keep executing
     if (.not. (tf .le. num_timepoints)) cycle
     write(*,*) " Iteration " , ti, tf , num_timepoints
     
     !update time axis bounds
     box_bounds(3,1)= (ti-1)*time_step
     box_bounds(3,2)= (tf-1)*time_step
     
     !write(*,*)  " Iteration  " , ti, tf, box_bounds(3,1) , box_bounds(3,2)
    do k=1, mc_iters_max 
        !prepare shuffled groups
        if (k .gt. 1) then
            call prepare_groups_random(.TRUE.)
        else
            call prepare_groups_random(.FALSE.)
        end if
        
        do u =1,num_users
            do g=1,num_groups
                !call reset_random_number( 314 )

                call mask_clear()
                !reset user mask
                user_selector=.FALSE.
                do i=1,group_sizes(g)
                    !select a user from the group
                    up=groups_shuffled(g,i)
                    !skip the user itself
                    if (up .eq. u) cycle
                    !activate the mask for that user
                    user_selector(up)=.TRUE.
                end do
                
                call setup_mask_values(user_selector, num_users, ti, tf)

                call convolve_3d()
                call setup_observation_points(  obs, nobs_max, nobs, u ,ti, tf  )

                !obs=1    !  gaze_in_points(ti:tf ,u,:  )
                NSS_score=calc_NSS_grid(obs,nobs  )
                if (export_result_file .eqv. .TRUE.) then
                    write(fid_NSS_out,*) NSS_score , NSS_single, NSS_score/NSS_single , u, g, ti, k
                end if
                !call copy_mask_to_real()
                !write(*,*) "NSS score ", NSS_score , NSS_single, NSS_score/NSS_single
            end do
        end do
    end do
 end do
 
  call CPU_time(time_g2)
  write(*,*) "NSS score ", NSS_score , NSS_single, NSS_score/NSS_single  
    
  !record execution time
  
  time_g_tot=time_g2-time_g1
  write(*,*) "Total time ",  time_g_tot
    
  call  write_time_file(time_g_tot )
  
  call copy_mask_to_real()
  call export_array_double2D(kernel_space_r, "kernel_2D.dat")  
  call export_array_double2D( conv_space_r, "convolved_2D.dat")  
  call export_array_3D(conv_space_r, "fixation_map", output_single=.TRUE.)
  call export_array_double2D(kernel_original, "kernel_original_2D.dat")  
    
end subroutine grid_calculation






subroutine exact_calculation( approx  )  
  double precision :: time_g1=0.d0, time_g2=0.d0, time_g_tot=0.d0 , time_rep1=0.d0, time_rep2=0.d0
  double precision :: time_setup1=0.d0, time_setup2=0.d0 , time_setup=0.d0
  double precision :: time_mc_prep1=0.d0, time_mc_prep2=0.d0
  
  integer            :: counter =0
  double precision   :: NSS_score=0.d0  , NSS_single=0.d0
  
  integer :: u  !user being analyzed
  integer :: ti , tf  !initial and final timepoint of window  
  logical, allocatable :: user_selector(:)  !indicates if a user is to be used for the current window
  logical, allocatable :: user_selector2d(:,:)
  double precision :: ave_field !average of the field
  double precision :: std_field !standard deviation of the field
  double precision :: ovlp  !overlap matrix element contributions
  double precision :: sim
  double precision :: volume 
  double precision :: v(3)
  double precision :: margin_time=0.d0
  integer :: i, k
  integer :: eff_dim     !effective dimensionality 
  logical :: is_planar = .FALSE.
  logical, intent(in) :: approx 
  integer :: g  !group
  integer :: gsize  !group size
  integer :: g_uind1, g_uind2  !index us user within group
  integer :: up, upp  !user prime and double prime
  !allocate( user_selector(num_users) , user_selector2d(num_users, num_users) )
    
  integer, parameter :: fid_NSS_out=16
  character*255       :: filename_NSS_out
  
    write(*,*) "Exact solution chosen"
    if (planar .eq. 1) then
        !we use only two dimensions, for we are on a plane
        eff_dim=2
        is_planar=.TRUE.
        write(*,*) "2D calculation"
    else
        eff_dim=3
        is_planar=.FALSE.
        write(*,*) "3D calculation"
    end if
    
    !open file where score will be exported
    filename_NSS_out= 'NSS_score_export.dat'
    open(fid_NSS_out, file=trim(filename_NSS_out), action='write' )

    call initialize_arrays_exact()

    !Apply automatic padding in time only if using the approximate form
    if (approx .eqv. .TRUE.) then
        margin_time=5.d0 !number of standard deviations
        write(*,*) "Applying padding to time axis of sigma=", margin_time
        !This will prevent the approximation from failing 
    else
        margin_time=0.d0
    end if
    
    box_bounds(3,1)=   -margin_time*sigma_array(3)
    box_bounds(3,2)= (points_t-1 )*time_step  +margin_time*sigma_array(3)    
    
    v=0.5d0/sigma_array**2
    volume=1.d0
    
    
    !we only use the effective dimensions for volume 
    do i=1,eff_dim
        volume=volume*( box_bounds(i,2) -box_bounds(i,1)  )
        
    end do
        
    NSS_single=NSS_calc_single(is_planar, approx)
    
    call CPU_time(time_g1)
        
    
 do ti=1,num_timepoints, window_stride   
     call CPU_time(time_setup1)
     tf=ti+points_t-1   !end of the time interval
     
     !If the window falls out of bounds then don't keep executing
     if (.not. (tf .le. num_timepoints)) cycle
     write(*,*) " Iteration " , ti, tf , num_timepoints
     
     !update time axis bounds
     box_bounds(3,1)= (ti-1)*time_step  -margin_time*sigma_array(3)
     box_bounds(3,2)= (tf-1)*time_step  +margin_time*sigma_array(3)
     
     
    call setup_arrays_exact(ti, tf, is_planar, approx)
    !write(44,*) similarity
    !write(45,*) overlap*norm**2
    !write(46, *) average_field*norm/volume
    call CPU_time( time_setup2 )
    call CPU_time( time_rep1 ) 
    do k=1, mc_iters_max
        !prepare new groups
        if (k .gt. 1) then
            call CPU_time( time_mc_prep1 ) 
            call prepare_groups_random(.TRUE.)
            call CPU_time( time_mc_prep2 ) 
            !write(177, *)k, time_mc_prep2-time_mc_prep1
        else
            call prepare_groups_random(.FALSE.)
        end if
        
        do u =1,num_users
            do g=1,num_groups
                !user_selector= group_masks_1d(g,:)       
                !user_selector2d=group_masks_2d(g,:,:)  

                !remove user from field
                
                !user_selector(u)=.FALSE.
                !user_selector2d(u,:)=.FALSE.
                !user_selector2d(:,u)=.FALSE.
                
                gsize=group_sizes(g)
                ave_field=0.d0
                ovlp=0.d0
                sim=0.d0
                do g_uind1=1,gsize
                    !define reference user u prime from group
                    up=groups(g,g_uind1)
                    if (u .eq. up) cycle
                    !contribution of average field due to all the other users in group
                    ave_field=ave_field+average_field(up)
                                        
                    !similarity between uprime and u
                    sim=sim+similarity( up,u )
                    do g_uind2=1,gsize
                        upp= groups(g,g_uind2)
                        if (u .eq. upp) cycle
                        ovlp=ovlp+overlap(up, upp)
                    end do
                    
                end do
                ave_field=ave_field/volume
                std_field=( (ovlp/volume) -ave_field**2 )**0.5
                sim=sim/valid_points(u)
                !write(134,*) ave_field, std_field, sim
                
                !ave_field=SUM( average_field, MASK=user_selector )/volume
                !std_field=( SUM( overlap , MASK=user_selector2d   )/volume -ave_field**2 )**0.5
                !sim=(SUM( similarity( :,u ) , MASK=user_selector )   )/valid_points(u)
                
                
                NSS_score=(sim-ave_field)/std_field
                if (export_result_file .eqv. .TRUE.) then
                    write(fid_NSS_out,*) NSS_score , NSS_single, NSS_score/NSS_single , u, g, ti,k
                end if
                !write(*,*)  "sim ave std  ", sim, ave_field, std_field
                !write(*,*) "NSS score ", NSS_score , NSS_single, NSS_score/NSS_single
                i=i+k+u
            end do
        end do
    end do
    call CPU_time(time_rep2)
 end do
    call CPU_time(time_g2)
    time_g_tot=time_g2-time_g1
    time_setup=time_setup2-time_setup1
    

    write(*,*) "Total time ", time_g_tot , " repetition time ", (time_rep2-time_rep1)/num_users,"setup ", time_setup
    write(77,*) "Total time ", time_g_tot, " repetition time ", (time_rep2-time_rep1)/num_users, "setup ", time_setup 
    
    call write_time_file(time_g_tot )
    
    write(*,*) "NSS score ", NSS_score , NSS_single, NSS_score/NSS_single
    write(*,*) " i ", i
end subroutine exact_calculation


subroutine write_time_file(time_measured)
    integer , parameter :: fid_time_meas=778
    character*255       :: filename_time_meas
    double precision , intent(in) :: time_measured
    
    filename_time_meas= 'measured_time.dat'
    open(fid_time_meas, file=trim(filename_time_meas), action='write' )
    write(fid_time_meas,*) "Total time ", time_measured
    close(fid_time_meas)
    
end subroutine write_time_file

subroutine prepare_groups_random(shuffle_active)
    implicit none
    integer :: i,j
    integer :: gindex
    logical, intent(in) :: shuffle_active
    gindex=1
    !unload groups onto flattened array
    do i=1,num_groups
        do j=1,group_sizes(i)
            groups_flat(gindex)=groups(i,j)
            gindex=gindex+1
        end do
    end do
    
    !only if shuffle is active 
    if (shuffle_active .eqv. .TRUE.) then
        call Shuffle_buf(groups_flat)
    end if
    
    !load the flattened array into the shuffled group array
    gindex=1
    do i=1,num_groups
        do j=1,group_sizes(i)
            groups_shuffled(i,j)=groups_flat(gindex)
            gindex=gindex+1
        end do
    end do    
    
end subroutine prepare_groups_random



end program nss_calculator

