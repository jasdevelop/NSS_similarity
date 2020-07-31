!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   data_input.f08
! Author: jorge
!
! Created on August 8, 2019, 12:56 AM
!

MODULE input_module

implicit none
private
save
integer, public ::  num_users=1
integer, public ::  num_groups=1
integer :: num_iterations=1

!Monte Carlo parameters
integer, public :: mc_iters_max =1  !number of Monte Carlo iterations
logical, public :: monte_carlo_enabled=.FALSE.

double precision, allocatable , public   :: gaze_data(:,:,:)
integer , allocatable,  public        ::  gaze_in_points(:,:,:)
!integer, allocatable    :: group_mat_ids(:,:) !user ids to consider
!integer, allocatable    ::  group_mat_n(:) !number of users in group
!integer, allocatable    ::  iter_info(:,:) !information on iteration with frame, user and group 

logical, public :: export_result_file =.FALSE.


!Store user group data and size of each group
integer , allocatable, public :: groups(:,:) , group_sizes(:)


!Stores all users in all groups  and backup for shuffling 
integer, allocatable, public :: groups_flat(:), groups_shuffled(:,:)


!input parameters
double precision, public    ::  samp_freq
double precision, public    ::  time_step
integer, public             ::  points_t                
integer, private            ::  points_t_in
double precision, public    ::  fsamp_factor=1.d0               

integer              ::  points_x
integer              ::  points_y
                
double precision            ::  start_x                
double precision            ::  end_x
                
double precision            ::  start_y                  
double precision            ::  end_y  

double precision            ::  space_x                
double precision            ::  space_y


integer, public             ::  window_stride

integer, public             ::  planar=0
integer, public             ::  mode_code
integer , public            ::  num_timepoints


double precision, public , allocatable   :: axis_x(:) , axis_y(:), axis_t(:) , axis_all(:,:)
integer, public              ::  box_shape(3)   
double precision, public     ::  box_bounds(3,2) , box_spacing(3)
double precision, public     ::  sigma_array(3) , sigma_in_points(3)
double precision, public     ::  kernel_size_std(3)
integer , public             ::  kernel_shape(3)

double precision , parameter :: min_space_box =1.0/1d6

public :: read_parameters , setup_axis_all , read_data , conv_coord_to_points, read_groups

CONTAINS 

subroutine read_data
    implicit none
    character*255,parameter ::  filename_data="gaze_data.dat"
    integer, parameter      ::  fid_data=100
    integer                 ::  num_lines , line_index=0 
    double precision        ::  x, y, t
    integer                 ::  user
    integer                 ::  it    
    integer                 ::  timepoint
        
    
    open(fid_data,file=filename_data ,action='read')
    read(fid_data,*) num_lines , num_timepoints , num_users      
    write(*,*) 'Lines: ', num_lines, ' ', ' time_points: ', num_timepoints , 'num_users: ', num_users
    
    allocate( gaze_data(num_timepoints,num_users, 3) )
    gaze_data=0.d0
    read(fid_data, *) 
    do it=1,num_lines 
        read(fid_data, *) line_index, timepoint,  user, x,y,t
        !write(*,*) 'line ', line_index, time, frame, user, x,y,z , '\n'
        
        gaze_data(timepoint,user,1)=x
        gaze_data(timepoint,user,2)=y
        gaze_data(timepoint,user,3)=t
                
    end do
        
    close(fid_data)
    
end subroutine read_data



subroutine read_parameters()
    implicit none
    integer, parameter :: num_lines_max=20
    integer, parameter :: fid_param=250
    character*255,parameter :: filename_data="parameters.txt"
    character*255 :: param_name, param_val
    integer :: i
    integer :: export_result_file_int=1
    integer ::  num_lines 
    
    num_lines=count_lines(filename_data, num_lines_max)
    
    open(fid_param,file=filename_data ,action='read')
    
    !the order of the lines does not matter 
    do i=1,num_lines
        read(fid_param, * ) param_name, param_val
                
        !parse the read in value
        select case( trim(param_name)  )  
            !1
            case('samp_freq')
                read( param_val,* ) samp_freq
            !2
            case('points_t')
                read( param_val,* ) points_t_in                
            !3
            case('sigma_t')
                read( param_val,* ) sigma_array(3)    
            !4
            case('window_stride')
                read( param_val,* ) window_stride    
            !5    
            case('sigma_xy')
                read( param_val,* ) sigma_array(1)         
                sigma_array(2)=sigma_array(1)

            !6    
            case('kernel_size_xy')
                read( param_val,* ) kernel_size_std(1)
                kernel_size_std(2)=kernel_size_std(1)
            !7    
            case('kernel_size_t')
                read( param_val,* ) kernel_size_std(3)                 

            !8    
            case('space_xy')
                read( param_val,* ) space_x
                space_y=space_x
                            
            !9    
            case('start_x')  
                read( param_val,* ) start_x                
            !10
            case('end_x')  
                read( param_val,* ) end_x
                
            !11    
            case('start_y')  
                read( param_val,* ) start_y                  
            !12
            case('end_y')  
                read( param_val,* ) end_y                 
            !13    
            case('mode')
                read(param_val, *) mode_code 
               
            !14
            case('planar')
                read( param_val,* ) planar         
                
            !15
            case('mc_iters')
                read( param_val,* ) mc_iters_max   
                if (mc_iters_max > 1) then
                    monte_carlo_enabled = .TRUE.
                else 
                    monte_carlo_enabled = .FALSE.
                end if
            
            !16    
            case('export')
                read(param_val,*) export_result_file_int
                export_result_file=export_result_file_int >0  
            
            !17
            case('fsamp_factor' ) 
                read( param_val,* ) fsamp_factor
                if (fsamp_factor .ne. 1) then
                    call failure_sub(1, trim("Only the value 1.0 is currently permitted for the fsamp_factor " )  )
                end if
                
            case default                
                !write(*,*) 'unknown parameter ', param_name
                call failure_sub(1, trim("Unknown parameter in input file "//param_name )  )
                       
        end select
               
    end do
    

end subroutine read_parameters  

function parse_mode(mode_text)
    implicit none
    character*255 :: mode_text
    character*255 :: text_err
    integer :: parse_mode 
    character*4 , parameter :: calc_type_grid='grid'
    character*5 , parameter :: calc_type_exact='exact'
    
    select case( trim(mode_text) )
        case(calc_type_grid)
            parse_mode=1
        case(calc_type_exact )
            parse_mode=2
        case default
            write(text_err,fmt='(A)') trim('Unknown value for mode: '//mode_text ) 
            call failure_sub(1, trim(text_err) )
    end select
    
end function parse_mode

function count_lines(filename_data, maxlines)
    implicit none
    character(len=*) :: filename_data   
    integer, intent(in) :: maxlines
    integer:: count_lines 
    integer :: i , counter
    integer, parameter :: fid_param=251    
    character*255 :: param_name, param_val
    integer :: io
    
    !This function counts the number of contiguous lines that are valid inputs
    counter=0
    open(fid_param,file=trim(filename_data) ,action='read')
    do i = 1, maxlines
        read(fid_param, *, IOSTAT=io)  param_name, param_val
        
        if (io .eq. 0 ) then
            counter=counter+1
        else
            EXIT
        end if
    end do
    close(fid_param)
    write(*,*) "Detected ", counter, " valid lines."
    count_lines=counter
    
end function count_lines

subroutine failure_sub( err_code, err_message )
    ! gives an error message an halts execution
    ! it is up to the user to decide what yo give as en error code
    implicit none
    integer , intent(in):: err_code   !error code from the routine that originated the issue
    character(len=*) :: err_message    !error message
    print '("  Error, status = ",I0)', err_code
    print *, trim(err_message)
    print *, "Terminating"
    call exit(1)
    
end subroutine failure_sub

subroutine linspace( array, p_start, p_stop  , npoints )
    !Similar to the python and Matlab function
    !This one allocates an array of size npoints, starting at p_start 
    ! and increasing in regular intervals up to p_end.  
    double precision, allocatable :: array(:) 
    integer          :: i
    integer          :: npoints
    double precision :: p_start , p_stop
    double precision :: delta
    integer :: status
    status=0
    allocate( array(npoints) , STAT=status )
    if (status.ne.0) call failure_sub(1, "Could not allocate memory for axis ")
    
    array(1)=p_start
    
    delta= (p_stop-p_start)/(npoints -1)
    do i=2, npoints
        array(i) = p_start+  delta*(i-1)
    end do
    
end subroutine linspace


subroutine setup_axis_all()
    implicit none
    integer :: i , j , max_n , n
    integer :: status
    status=0
    write(*,*) "Setting up x axis "
    call generate_range( axis_x, start_x, end_x, space_x )
    points_x=size(axis_x,1)
    write(*,*) "Setting up y axis "
    call generate_range( axis_y, start_y, end_y, space_y )
    points_y=size(axis_y,1)
    
    write(*,*) "Setting up time axis"
    points_t=NINT(points_t_in*fsamp_factor)
    time_step=1/(samp_freq*fsamp_factor)
    
    allocate( axis_t(points_t) )
    
    do i=1, points_t
        axis_t( i )= time_step*(i-1)
    end do
    write(*,*) "Time points ", points_t,  " Time step ", time_step
    box_shape(1)=points_x
    box_shape(2)=points_y
    box_shape(3)=points_t
    
    !Store the start and end points for each axis
    box_bounds(1,1)=start_x
    box_bounds(1,2)=end_x
    box_bounds(2,1)=start_y
    box_bounds(2,2)=end_y
    
    box_bounds(3,1)=axis_t(1)
    box_bounds(3,2)=axis_t(points_t)
    
    
    ! Store the shapes of all the axis into one array, zero padded for the shorter dimensions
    max_n=MAXVAL( box_shape  )
    write(*,*) "Maximum size is " , max_n    
    allocate( axis_all(3, max_n   ) )    
    axis_all=0.d0    
    axis_all( 1 ,1:box_shape(1) ) = axis_x(:)        
    axis_all( 2 ,1:box_shape(2) ) = axis_y(:)   
    axis_all( 3 ,1:box_shape(3) ) = axis_t(:)   
      
    
    !Store the grid spacing values (pixel or time units)
    do i=1,3
        box_spacing(i)=axis_all(i, 2 )-axis_all( i, 1)
    end do
    
    do i=1,3
        if ( box_spacing(i) > min_space_box ) then
            ! case of spacing that is not too small
            n=CEILING(kernel_size_std(i)*sigma_array(i)/box_spacing(i) )
        else
            !For cases with singularities 
            n=1
        end if
        if (modulo(n,2) .eq.0) n=n+1
        kernel_shape(i)=n
        write(*,*) "Number of gridpoints in kernel ", n
    end do
    
    !Calculate the sigmas in units of gridpoints
    do i=1,3
        if (box_spacing(i) > min_space_box ) then
            sigma_in_points(i)=sigma_array(i)/box_spacing(i)
        else
            write(*,*) "Singularity case in box spacing ", i
            sigma_in_points(i)=1.d0/min_space_box
        end if
        
    end do
    
end subroutine setup_axis_all


function round_n_decimals(value_in,n)
    double precision, intent(in) :: value_in 
    integer, intent(in) :: n
    double precision  :: round_n_decimals
        
    round_n_decimals=ANINT(value_in*10**n, 8)/10**n
    
end function round_n_decimals

subroutine generate_range(array, p_start, p_end, p_step  )
    ! generate a range, with given spacing 
    ! the endpoint is included
    double precision, allocatable :: array(:)
    integer :: n 
    integer :: i
    double precision :: p_start, p_end, p_step
    
    n= CEILING( (p_end-p_start)/p_step  ) +1
    allocate(array(n) )
    
    do i=1,n
        array(i) = p_start + p_step*(i-1)
    end do
    
end subroutine generate_range


subroutine conv_coord_to_points(data_in, data_out )
    implicit none
    double precision :: data_in(:,:,:)
    integer, allocatable, intent(out) :: data_out(:,:,:)    
    integer :: sx, sy, sz , i, j , k

    sx=size(data_in, 1)
    sy=size(data_in, 2)
    sz=size(data_in, 3)
    if ( allocated( data_out  ) ) then
        call failure_sub(1,"The output data array is already allocated ")
    end if        
    allocate(  data_out(  sx, sy, sz )   )
    
    do k=1,3
        data_out(:,:,k)=NINT( (data_in(:,:,k)-box_bounds(k,1) )/box_spacing(k)  )+1
    end do
    
end subroutine conv_coord_to_points

subroutine read_groups()
    implicit none    
    integer, parameter :: fid_group=222
    character*255,parameter :: filename_group="groups.txt"   
    integer :: index_group, size_group
    integer :: i, j

    logical, allocatable :: mask_temp(:,:)
    integer  :: max_group_size
    integer :: groups_enabled
    integer :: gindex
    integer :: tot_group_users
    
    !no group is allowed to repeat users
    max_group_size=num_users
    
    open(fid_group,file=filename_group ,action='read')
    read(fid_group,*) num_groups  
       
    
    write(*,*) "We have ", num_groups , " groups " 
    
    allocate( groups(num_groups, max_group_size ) )
    allocate( groups_shuffled(num_groups, max_group_size ) )
    allocate( group_sizes(num_groups)  )

    
    
    
    write(*,*) "Allocated memory for groups ", num_users
    groups=0
    groups_shuffled=0
    
    
    do i=1,num_groups
        read(fid_group,*) size_group
        group_sizes(i)=size_group
        !Store the user id of the members of the group
        read(fid_group,*) (groups(i,j),j=1,size_group)


    end do
    
    close(fid_group)
    !write(*,*) groups
    !write(*,*) group_masks_1d
    !write(*,*) "2D masks "
    !write(*,*) group_masks_2d
    groups_shuffled=groups
    
    tot_group_users=SUM( group_sizes )
    allocate(groups_flat(tot_group_users))
    
    !mc_mask_temp=.FALSE.
    !do i=1, num_groups
    !    groups_all(gindex:gindex+group_sizes(i)-1 )=groups(i,1:group_sizes(i) )
    !    !write(*,*) "From ", gindex , gindex+group_sizes(i) , i
    !    gindex=gindex+group_sizes(i)
    !end do    
    !groups_shuffled=groups_all
    !write(*,*) "All groups"
    !write(*,*) groups_all
end subroutine read_groups







END MODULE input_module
