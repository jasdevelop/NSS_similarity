!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   convolutionfft3d.f08
! Author: jorge
!
! Created on August 7, 2019, 11:35 PM
!



MODULE fft_convolve_3d
  private
  save
  include 'fftw3.f'
  
    
  !Array sizes
  
  integer  :: L1=0   , L2=0  ,  L3=0    !size of the unpadded mask for each dimension
  integer  :: N1=0   , N2=0  ,  N3=0   !size of first and second dimension of padded mask
  integer ::  M1=0   , M2=0  ,  M3=0   ! size of kernel 
  integer :: pad1=0, pad2=0, pad3=0
  
  integer, public  :: calc_dim=0             !effective dimensionality for grid-purposes, takes values 1,2 or 3 
  logical, public ::  planar_mode=.FALSE.   ! if active the grid is planar , but the time windows are still kept    
  
  double precision , parameter :: min_sigma =0.1 !minimum value for sigma to be valid
  
  ! Data arrays 
  ! the suffix r indicates real values, while c is for complex
  
  double precision, allocatable, public :: kernel_original(:,:,:)
    
  !arrays for the mask, padded with zeros
  complex*16, allocatable :: mask_array_c(:,:,:)
  double precision, allocatable :: mask_array_r(:,:,:)

  complex*16, allocatable :: kernel_c(:,:,:)
  double precision, allocatable, public :: kernel_space_r(:,:,:)
  
  double precision, allocatable, public :: conv_space_r(:,:,:)
  
  ! FFTW plans
  integer*8 :: plan_fwd_mask = 0
  integer*8 :: plan_fwd_kernel = 0
  
  integer*8 :: plan_back_mask = 0   
  integer*8 :: plan_back_kernel =0  !Not really needed

  integer :: status =0
  
  
  public :: initialize_3D , mask_set_random_values_3D, convolve_3D, random_observations_3D, calc_NSS_grid, &
            copy_mask_to_real, export_array_double2D, export_array_3D , mask_clear, mask_add, mask_add_multiple, &
            get_max_from_array , mask_add_point , calc_NSS_single
  
  
contains

function get_max_from_array()
    double precision :: get_max_from_array
    
     get_max_from_array=MAXVAL( conv_space_r  )
    
end function get_max_from_array

subroutine set_grid_sizes( mask_nx , mask_ny, mask_nz, kernel_nx,  kernel_ny, kernel_nz, planar )
    implicit none
    integer :: kernel_nx,  kernel_ny, kernel_nz
    integer :: mask_nx , mask_ny, mask_nz
    integer :: planar
    !initialize padding variables
    !These paddings are for masks that are too small for the kernel
    pad1=0
    pad2=0
    pad3=0
    
  !Set array sizes 
    
    L1=  mask_nx
    L2=  mask_ny
    L3=  mask_nz
    
    M1=kernel_nx
    M2=kernel_ny
    M3=kernel_nz

!Check if the planar mode is activated and adjust L3    
    if (planar .eq. 1 ) then
        planar_mode =.TRUE.
        L3=1
    end if
    
    ! Set and adjust the size of the kernel grid 
    
    ! We have 3 cases supported: 1D, 2D and 3D
    
    if ( (  L2 .eq. 1  ) .and. (L3 .eq. 1 ) ) then
        ! For 1D  case L2 and L3 are 1
        !Flatten y and z and make the kernel a line
        M2=1
        M3=1
        calc_dim=1
        write(*,*) "Performing a 1D calculation"
    else if ( ( L1 .gt. 1  ) .and. ( L2 .gt. 1 ) .and. ( L3 .eq. 1) ) then
        !For 2D we have that L3 is 1, and L1 and L2 greater than 1        
        !We use flattened z 
        M3=1
        calc_dim=2
        write(*,*) "Performing a 2D calculation"
    else if (  ( L1 .gt. 1  ) .and. ( L2 .gt. 1 ) .and. ( L3 .gt. 1 )  ) then    
        !3D case  
        calc_dim=3
        write(*,*) "Performing a 3D calculation"
        !Provide possible extra padding to any dimension that is too short
        if (L1 .lt. (M1/2)+1  ) then
            pad1= M1/2+1
            write(*,*) "Padding 1 with ", pad1
        end if
        
         if (L2 .lt. (M2/2)+1  ) then
            pad2= M2/2+1
            write(*,*) "Padding 2 with ", pad2
        end if       
        
        if (L3 .lt. (M3/2)+1  ) then
            pad3= M3/2+1
            write(*,*) "Padding 3 with " , pad3
        end if        
        
    else 
        call failure_sub(1, "Dimensions incorrectly set. The grid must be along x (1D) , xy (2D) , or xyz (3D) ")
        
    end if
    
    !M2 and M1 will always be identical 
    
    !Verification
    !Make sure the kernel has an odd number of points
    if ( (modulo(M1,2) .ne. 1).or. (modulo(M2,2) .ne. 1) .or. ( modulo(M3,2) .ne.1 ) )  then
        call failure_sub(status, "Kernel must have an odd number of points")
    end if
    !Stop initialization if the grid sizes are not valid
    if ( (M1 < 1 ) .or. (M2 < 1) .or. (M3 <1) .or. (L1 < 1) .or. (L2 <1) .or. ( L3<1 )  ) then
        call failure_sub( 1, "Array sizes need to be set to valid sizes before initializing" )
    end if    
    
end subroutine set_grid_sizes


subroutine initialize_3d( mask_nx , mask_ny, mask_nz, kernel_nx, kernel_ny ,kernel_nz,  sigma_x, sigma_y , sigma_z , planar )    
    implicit none
    integer :: kernel_nx, kernel_ny , kernel_nz
    integer :: mask_nx , mask_ny, mask_nz
    double precision :: sigma_x, sigma_y ,sigma_z
    integer :: planar
    

    
    !set the grid variables, which are global
    call set_grid_sizes(  mask_nx , mask_ny, mask_nz, kernel_nx,  kernel_ny, kernel_nz, planar )
    
    !Validating sigmas
    if ( ( sigma_x.lt. min_sigma ) .or. (  sigma_y .lt. min_sigma ) .or. ( sigma_z .lt. min_sigma ) ) then
        call failure_sub( 1, "Sigma is too small." )
    end if

    !Allocate memory for kernel array
    allocate( kernel_original(M1, M2, M3 ) ,    STAT=status )
    if (status .ne. 0) call failure_sub(status, "Bad allocation of core kernel array")
  
  !Initialize array
  kernel_original=0.d0
   
  !Define enlarged sizes, padding with extra zeros 
  !M1 is always odd and symmetric, so we pad with M/2 (integer division) zeros, plus some extra (sometimes, for small L)
  !N1, N2 and N3 are global variables in the module
  N1=L1+M1/2 +pad1
  N2=L2+M2/2 +pad2
  N3=L3+M3/2 +pad3  
  
  write(*,*) "Physical size of box : " , N1, N2, N3
  write(*,*) "Logical size of box : " , L1, L2, L3
  write(*,*) "Kernel box size : " , M1, M2, M3
  
 ! initialize arrays and prepare plans
  call setup_plans_allocate(N1, N2, N3, mask_array_c , plan_fwd_mask, plan_back_mask )
  call setup_plans_allocate(N1, N2, N3,  kernel_c , plan_fwd_kernel, plan_back_kernel )

  !Allocate array for the real values answer output
  allocate( kernel_space_r(N1, N2, N3 )  )  

  kernel_space_r=0.d0
  !Prepare kernel arrays
  write(*,*) "Kernel setup initiated ..."
  call init_gaussian3D(kernel_original, sigma_x, sigma_y,sigma_z) 
  write(*,*) "Rolling ..."
  call kernel_pad_roll3D(kernel_original, kernel_space_r)
  
  kernel_c=kernel_space_r
  write(*,*) "Preparing kernel in Fourier domain ..."
  call dfftw_execute_dft(plan_fwd_kernel, kernel_c, kernel_c )
  write(*,*) "Kernel setup complete"  
end subroutine initialize_3d

subroutine copy_mask_to_real()
  implicit none
  integer :: i,j,k
  double precision :: thres=1.d0
  !Allocate array for the real values answer output, if not done so already
  if (allocated(conv_space_r) .eqv. .false.  ) then
      allocate(conv_space_r(L1, L2, L3)  )
  end if
  
  !Copy the array , taking the real part since the source is complex
  conv_space_r( 1:L1, 1:L2, 1:L3 )=REALPART( mask_array_c( 1:L1, 1:L2, 1:L3 ) )

  ! Get the value of the machine epsilon
  thres=5*EPSILON(thres) 
  !Purge the array from small value errors, especially negative numbers
  forall (i=1:L1, j=1:L2, k=1:L3 ,  conv_space_r(i,j,k) .lt. thres  )
      !Set the value to zero 
      conv_space_r(i,j,k)=0.d0
  end forall
  
  
end subroutine

subroutine set_mask_value(x,y,z  )
    !deprecated
    implicit none
    !sets the specified position of the mask to one
    integer :: x,y,z
    !Check the bounds 
    if ( (x >0).and.(y>0).and.(z>0).and.(x .le. L1) .and. (y .le. L2) .and. ( z .le. L3 ) ) then
        !Set the value
        mask_array_c(x,y , z )=1.d0
    else
        call failure_sub(1,"Index out of bounds for the mask" )
    end if
    
end subroutine set_mask_value

subroutine mask_add(x,y,z ,amplitude )
    !deprecated
    implicit none
    double precision , intent(in):: amplitude
    !adds a specified value to the  position of the mask 
    integer :: x,y,z
    !Check the bounds 
    if ( (x >0).and.(y>0).and.(z>0).and.(x .le. L1) .and. (y .le. L2) .and. ( z .le. L3 ) ) then
        !Set the value
        mask_array_c(x,y , z )= REALPART( mask_array_c(x,y , z )) + amplitude
    else
        !call failure_sub(1,"Index out of bounds for the mask" )
        continue
    end if
    
end subroutine mask_add

subroutine mask_add_point(p ,amplitude , status)
    implicit none
    double precision , intent(in):: amplitude
    !adds a specified value to the  position of the mask 
    integer :: x,y,z
    integer :: p(3)
    integer , intent(out) :: status
    
    
    x=p(1)
    y=p(2)
    z=p(3)
    
    !Check the bounds 
    if ( (x .ge. 1 ).and.(y .ge.1 ).and.(z .ge. 1).and.(x .le. L1) .and. (y .le. L2) .and. ( z .le. L3 ) ) then
        !Set the value
        mask_array_c(x,y , z )= REALPART( mask_array_c(x,y , z )) + amplitude
        !No problem
        status=0
    else
        !call failure_sub(1,"Index out of bounds for the mask" )
        !out ob bounds
        status=1
        
    end if
    
end subroutine mask_add_point

subroutine mask_add_multiple(points , user_exclude )
    !the input unit is gridpoints 
    implicit none
    integer :: points(:,:,:)  !  (point index, user, coordinate)
    integer :: user_exclude   !selected users to add
    integer :: n , i
    integer :: u , num_u
    
    n=size(points,1)
    num_u=size(points, 2)
    do u=1, num_u
        if ( u .eq. user_exclude) cycle   
        do i=1, n                      
            call mask_add( points(i, u ,1 ) , points(i, u ,2) , points(i, u ,3) , 1.d0   )        
        end do
    end do
end subroutine mask_add_multiple


subroutine mask_clear()
    implicit none
    mask_array_c=0.d0
end subroutine


subroutine mask_set_random_values_3D()
    !Creates random Gaussians to put into the box 
    implicit none
    integer , parameter :: nr=20
    real    :: array1(nr) , array2(nr),  array3(nr)
    real  :: numx, numy, numz
    integer :: i =1 , pos1, pos2, pos3
    
    !Reset the mask
    mask_array_c=0.d0
    
    call RANDOM_NUMBER(array1)
    call RANDOM_NUMBER(array2)
    call RANDOM_NUMBER(array3)
    
    do i=1, nr
        numx=array1(i)
        numy=array2(i)
        numz=array3(i)
        pos1=INT( FLOOR(  numx*(L1)   )   ) +1
        pos2=INT( FLOOR(  numy*(L2)  )   )  +1
        pos3=INT( FLOOR(  numz*(L3)  )   )  +1
        
        mask_array_c( pos1, pos2, pos3 )=1.d0
        !write( 99,* ) pos1, pos2, pos3        
        
    end do
    
    
end subroutine mask_set_random_values_3D    

subroutine random_observations_3D(obs ,  nr)
    integer :: nr
    integer :: obs(:,:)
    
    real    :: array1(nr) , array2(nr),  array3(nr)
    real  :: numx, numy, numz
    integer :: i =1 , pos1, pos2, pos3
    
    !Reset the mask
    
    
    call RANDOM_NUMBER(array1)
    call RANDOM_NUMBER(array2)
    call RANDOM_NUMBER(array3)
    
    do i=1, nr
        numx=array1(i)
        numy=array2(i)
        numz=array3(i)
        pos1=INT( FLOOR(  numx*(L1)   )   ) +1
        pos2=INT( FLOOR(  numy*(L2)  )   )  +1
        pos3=INT( FLOOR(  numz*(L3)  )   )  +1
        
        obs(i,1 )=pos1
        obs(i,2 )=pos2
        obs(i,3) =pos3
        
    end do   
    
end subroutine random_observations_3D

function calc_NSS_grid(obs_points, npoints, bounds_in )
    implicit none
    double precision :: calc_NSS_grid   
    double precision :: mean 
    double precision :: std
    double precision :: F, NSS, NSS_ave
      
    integer :: obs_points(:,:)  !has dimensions npoints x 3  (ie x, y and z)  
    integer , intent(in) :: npoints   !number of observation points
    
    integer :: i    
    double precision :: FS  !sum of F
    double precision :: volume
    integer , optional :: bounds_in(3,2)   
    integer :: b(3,2)   !stores the boundaries for calculation
    
    if (present(bounds_in).eqv. .TRUE.) then
        !Override default bounds
        b(:,:)=bounds_in(:,:)
        
    else
        !Use the standard box boundaries
        b(:,1)=1  !starting boundary
        b(1,2)=L1
        b(2,2)=L2
        if (planar_mode .eqv. .TRUE.) then
            b(3,2)=1   !Is the calculation is a planar we only have 1 point
        else
            b(3,2)=L3  !This is the case for full 3d calculations
        end if
    end if
    volume=1.d0*PRODUCT( (b(:,2)-b(:,1) )+1 )
    
    mean=SUM(REALPART( mask_array_c( b(1,1):b(1,2) , b(2,1):b(2,2)  , b(3,1):b(3,2)  ) ) )/volume
    std=(SUM( ( REALPART( mask_array_c(  b(1,1):b(1,2) , b(2,1):b(2,2) , b(3,1):b(3,2)  ) ) - mean )**2 )/volume )**0.5
    
    FS=0.d0
    NSS_ave=0.d0
    do i=1, npoints
        !value of the field F at the observation locations 
        F= REALPART( mask_array_c( obs_points(i,1) ,   obs_points(i,2),  obs_points(i,3)    ) )
        FS=FS+F
        !convert to Z score
        NSS= (F-mean )/std
        NSS_ave=NSS_Ave+NSS        
    end do
    NSS_ave=NSS_ave/npoints
    FS=FS/npoints
    calc_NSS_grid=NSS_ave
    
end function calc_NSS_grid

function calc_NSS_single()    
    integer :: status 
    integer :: p(3) , p_aug(1,3)    
    double precision :: calc_NSS_single
    status=0
    p(1)=L1/2+1
    p(2)=L2/2+1    
    p(3)=L3/2+1
    
    
    p_aug(1,:)=p(:)
    
    call mask_clear()    
    call mask_add_point(p, 1.d0 , status)
    call convolve_3d()    
    calc_NSS_single=calc_NSS_grid(p_aug, 1  )
    
end function calc_NSS_single

subroutine convolve_3d()    
  ! Execute FFT space -> frequency domain
  call dfftw_execute_dft(plan_fwd_mask  ,  mask_array_c, mask_array_c )
  
  !Multiply and normalize
  mask_array_c=mask_array_c*kernel_c/(N1*N2*N3)
  
  ! Execute inverse FFT frequency -> space
  call dfftw_execute_dft(plan_back_mask, mask_array_c, mask_array_c)
  
end subroutine convolve_3d

subroutine export_array_double(array, filename_data)
    integer , parameter :: fid_param=333
    double precision :: array(:,:,:)
    character(len=*) :: filename_data !"data.txt"        
    open(fid_param,file=trim(filename_data) ,action='write')
    write(fid_param,*), array
    close(fid_param)        
end subroutine export_array_double

subroutine export_array_complex(array, filename_data)
    implicit none
    integer , parameter :: fid_param=333
    complex*16 :: array(:,:,:)
    character(len=*) :: filename_data !"data.txt"        
    open(fid_param,file=trim(filename_data) ,action='write')
    write(fid_param,*), REALPART(array), IMAGPART(array)
    close(fid_param)        
end subroutine export_array_complex

subroutine export_array_complex2D(array, filename_root, extension)
    !export z=1 plane for complex array
    implicit none
    integer :: n1,n2
    complex*16 :: array(:,:,:)
    !file names will be of the form "root_suffix.ext" , eg: data_imag.dat
    character(len=*) :: filename_root   ! root for the file name
    character(len=*) :: extension       ! extension for the file eg ".dat"
    character*255, parameter :: suffix_r="real" , suffix_i="imag"
    character*255 :: name_real_file, name_imag_file  ! we store the new names here
    integer :: fid_export=33
    integer :: i,j
    
    n1=size(array,1)
    n2=size(array,2)
    
    !prepare the name for the real part file    
    write(name_real_file,fmt='(A)') trim(filename_root)//'_'//trim(suffix_r)//'.'//trim(extension)
    !export it as if it as a double precision array
    call export_array_double2D( REALPART(array) , trim(name_real_file)  )
    
    !we repeat for the imaginary part
    write(name_imag_file, fmt='(A)') trim(filename_root)//'_'//trim(suffix_i)//'.'//trim(extension)
    call export_array_double2D( IMAGPART(array) , trim(name_imag_file)  )
           
    
end subroutine export_array_complex2D

! 

subroutine export_array_double2D(array, filename_data )
    !export the z=1 plane in double precision
    implicit none
    integer :: n1,n2
    double precision :: array(:,:,:)
    character(len=*) :: filename_data !"data.txt"   
    integer :: fid_export=33
    integer :: i,j
    n1=size(array,1)
    n2=size(array,2)
    
    !The file is written with the first dimension as rows and the second as columns in the file
    !The values are separated by spaces
    open(fid_export, file=trim(filename_data) , action="write")
    do i=1,n1
        do j=1,n2
            write(fid_export,fmt='(es15.6e3)', advance="no") array(i,j,1) 
        end do
        
        write(fid_export,* )
    end do
    
    close(fid_export)
    
end subroutine export_array_double2D

subroutine export_array_3d(array, filename_root, output_single )
    double precision :: array(:,:,:)
    logical :: output_single 
    integer, parameter :: gen_purpose_bovfile= 44  
    integer, parameter :: fid_file_bin=55
    character(len=*) :: filename_root   ! root for the file name
    character*255 :: fname_bov
    character*255 :: fname_dat
    character*255 :: variable_name
    character*255 :: output_precision
    double precision, parameter :: time_to_write =0.d0
    integer :: npoints(3) , ncenter(3)
    integer :: i
    
    
    variable_name="fixation_map"
    do i=1,3
        npoints(i)=size(array,i)
    end do
    
    ncenter=1  
    write(fname_dat,fmt='(A)') trim(filename_root)//'.bin'
    open(fid_file_bin , file=trim(adjustl(fname_dat)),form="unformatted",status='replace')
    
    if (output_single .eqv. .TRUE. ) then
        write(fid_file_bin) REAL(array, 4)
        write( output_precision, fmt='(A)' )  trim( 'FLOAT')        
    else 
        write(fid_file_bin) array
        write( output_precision, fmt='(A)' )  trim( 'DOUBLE')
    end if
    
    close(fid_file_bin)
    
    
    write(fname_bov,fmt='(A)') trim(filename_root)//'.bov'
    
    open(gen_purpose_bovfile,file=trim(adjustl(fname_bov)),status='replace')
    write(gen_purpose_bovfile,'(1x,a,1x,es23.15)') 'TIME:',time_to_write
    write(gen_purpose_bovfile,'(1x,a,1x,a)') 'DATA_FILE:',trim(adjustl(fname_dat))
    write(gen_purpose_bovfile,'(1x,a,3(1x,i6))') 'DATA_SIZE:',npoints(1), npoints(2), npoints(3)
    write(gen_purpose_bovfile,'(1x,a)') 'DATA_FORMAT:  '//trim(output_precision)
    write(gen_purpose_bovfile,'(1x,a,1x,a)') 'VARIABLE:',trim(adjustl(variable_name))
    !write(gen_purpose_bovfile,'(1x,a)') 'DATA_COMPONENTS: 2'
    write(gen_purpose_bovfile,'(1x,a)') 'DATA_ENDIAN:  LITTLE'
    write(gen_purpose_bovfile,'(1x,a)') 'CENTERING:  zonal'
    write(gen_purpose_bovfile,'(1x,a,3(1x,i6))') 'BRICK_ORIGIN:', ncenter(1) , ncenter(2) , ncenter(3)
    write(gen_purpose_bovfile,'(1x,a,3(1x,i6))') 'BRICK_SIZE:', npoints(1), npoints(2), npoints(3)
    write(gen_purpose_bovfile,'(1x,a)') 'BYTE_OFFSET: 4'
    close(gen_purpose_bovfile)

    
end subroutine export_array_3d



subroutine export_array_boolean2D(array, filename_data )
    !export the z=1 plane in double precision
    implicit none
    integer :: n1,n2
    logical :: array(:,:)
    character(len=*) :: filename_data !"data.txt"   
    integer :: fid_export=33
    integer :: i,j
    n1=size(array,1)
    n2=size(array,2)
    
    !The file is written with the first dimension as rows and the second as columns in the file
    !The values are separated by spaces
    open(fid_export, file=trim(filename_data) , action="write")
    do i=1,n1
        do j=1,n2
            write(fid_export,fmt='(es15.6e3)', advance="no") array(i,j) 
        end do
        
        write(fid_export,* )
    end do
    
    close(fid_export)
    
end subroutine export_array_boolean2D


subroutine export_all()
  call copy_mask_to_real()
  call export_array_double( kernel_original, "kernel.dat" )
  
  call export_array_double( conv_space_r, "convolved_full.dat" )
  call export_array_complex( mask_array_c, "mask_freq.dat" )
  call export_array_double2D( conv_space_r(1:L1, 1:L2, 1:L3 )   , "convolved_2D.dat" )
  call export_array_double2D( kernel_space_r   , "kernel_rotated_2D.dat" ) 
  call export_array_double2D( kernel_original   , "kernel_2D.dat" )
  call export_array_complex2D( mask_array_c, "data", "dat" )
      
end subroutine export_all


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

subroutine roll_array3D(array_in, shifts , array_out )
    implicit none
    !roll the array by a number of points determined by shift, in the specified axis
    !the roll is done towards the origin for positive shift values
    ! this subroutine does the operation out-of-place , so array_in is unaltered
    integer , intent(in) :: shifts(3) !amount of points to roll for each axis
    
    double precision, intent(in) :: array_in(:,:,:)
    !the result is stored in array_out, which is also used as temporary storage
    double precision             :: array_out(:,:,:)
    double precision, allocatable :: temp(:,:,:)
    integer :: nx, ny ,nz !size of input array
    integer :: i   !indices for iteration
    integer :: ir   !rolled indices
        
    nx=size(array_in, 1)
    ny=size(array_in, 2)
    nz=size(array_in, 3)
    allocate(temp( nx,ny,nz  )  )
    
    temp=0.d0    
    !store the default output for no transformation
    array_out(:,:,:)=array_in(:,:,:)
        
    do i=0,nx-1
        !determine the target grid-point, "ir", to copy
        !we select the point that is located a certain distance away, given by "shifts" , with wrap-around 
        !we use the modulo operator instead of mod to allow the correct output with negative shifts
        !the index is chosen to be 0-based for simplicity of calculation with modulo, 
        !   and corrected later by adding 1 
        ir=modulo( i+shifts(1) , nx )
        !we copy the element into out temporary array, from out current working state
        temp( i+1, : ,: )=array_out(ir+1 , :,:)
    end do
    !update the result
    array_out(:,:,:)=temp(:,:,:)
    
    !repeat for y axis
    do i=0,ny-1
        ir=modulo( i+shifts(2), ny )
        temp( :, i+1 ,: )=array_out( :, ir+1 , : )
    end do
    
    array_out(:,:,:)=temp(:,:,:)    
    
    !repeat for z
    do i=0,nz-1
        ir=modulo( i+shifts(3), nz )
        temp( :,:, i+1  )=array_out( :, : , ir+1  )
    end do    
    
    array_out(:,:,:)=temp(:,:,:)    
    
    ! we discard the temporary array
    deallocate(temp)
    
end subroutine roll_array3D


subroutine kernel_pad_roll3D(kernel, padded)
    !
    ! Copy the kernel to the padded array and then roll it
    ! This puts the kernel segments at the edges, in the "wrapped-up" form needed for FFT
    implicit none
    double precision :: kernel(:,:,:)
    double precision :: padded(:,:,:)
   
    integer          :: mx, my , mz      ! sizes of the kernel array
    integer          :: shifts(3)

    mx=size(kernel, 1)
    my=size(kernel, 2)
    mz=size(kernel, 3)
    
    write(*,*)  "Kernel sizes ", mx, my, mz
    !if ( ( mx .gt. N1 ) .or. (my .gt. N2) ) then
    !    call failure_sub( 1, "The kernel size cannot be greater than the padded array." )
    !end if
    
    padded=0.d0
    write(*,*)  "Initiating kernel copy to padded array.. "
    write(*,*) "Padded ", size(padded,3)  , " kernel ", size(kernel,3)
    write(*,*) " values ", kernel(4,4,1)
    write(*,*) " value check done .."
    !copy the kernel to the origin of the padded array
    padded(1:mx, 1:my, 1:mz)=kernel(1:mx, 1:my, 1:mz)
    
    write(*,*) "Copy complete. "
    shifts(1)=  mx/2
    shifts(2)=  my/2
    shifts(3)=  mz/2
    write(*,*) "Shifts ",  shifts
    !roll the array to make the wrap-around form
    !we store the output on the same input array
    call roll_array3D(padded, shifts, padded)
    write(*,*) "Rolling complete"
end subroutine kernel_pad_roll3D




subroutine setup_plans_allocate(k1, k2, k3, x_cmplx, plan_fwd, plan_back)
  implicit none
    !set up execution plans for inverse FFT and allocate array
  ! the transformation is done in-place for array x_cmplx  
  integer  :: k1, k2 , k3
  ! Working precision is double precision
  integer, parameter :: WP = selected_real_kind(15,307)
  ! Data arrays
  complex(WP), allocatable :: x_cmplx(:,:,:)
  

  ! FFTW plan
  integer*8 :: plan_fwd, plan_back
  integer :: status
  status=0 
  
  print *,"Allocate data array "
  allocate ( x_cmplx (k1, k2, k3), STAT = status)
  if (0 /= status) call failure_sub(status, "Bad real array allocation")
  
  !Initialize
  x_cmplx(:,:,:)=0.d0
  print *,"Create FFTW complex-to-complex plans"
  
  ! prepare plan for forward transform
  call dfftw_plan_dft_3d(plan_fwd, k1, k2, k3, x_cmplx, x_cmplx, FFTW_FORWARD, FFTW_ESTIMATE)  
  !confirm 
  if (0 == plan_fwd) call failure_sub(status, "Bad forward DFT plan")
  
  !Prepare plan for inverse transform
  call dfftw_plan_dft_3d(plan_back , k1, k2, k3, x_cmplx, x_cmplx, FFTW_BACKWARD, FFTW_ESTIMATE)    
  ! confirm
  if (0 == plan_back) call failure_sub(status, "Bad backward DFT plan")
  write(*,*) "Plans all set up. "
end subroutine setup_plans_allocate


subroutine init_gaussian3D(array, sigma_x,  sigma_y, sigma_z )
  implicit none
  integer :: nx, ny, nz  !sizes of the array 
  integer :: kx, ky, kz,  k0x, k0y, k0z
  double precision , intent(in) :: sigma_x,sigma_y ,sigma_z  !the standard deviation of the Gaussian
  !The standard deviations must be given in units of grid-points within this entire module
  ! Working precision is double precision
  integer, parameter :: WP = selected_real_kind(15,307)
  real(WP)   :: array(:,:,:)
  array=0.d0
  nx=size(array,1)
  ny=size(array,2)
  nz=size(array,3)
  !locate the center of the array
  k0x=nx/2+1 
  k0y=ny/2+1
  k0z=nz/2+1
!evaluate the Gaussian function at all the points on the whole array     
   
  forall( kx=1:nx, ky=1:ny  , kz=1:nz)
      array(kx,ky, kz)=exp( -( (kx-k0x)**2)/(2.0*sigma_x**2) -((ky-k0y)**2 )/(2.0*sigma_y**2) -( (kz-k0z)**2 )/( 2.0*sigma_z**2 ) )
  end forall
  !write(*,*) "Maximum value ", MAXVAL(array)
end subroutine init_gaussian3D



subroutine init_sphere3D(array, radius )
   !Generates a grid with a sphere of the given radius  
  implicit none
  integer :: nx, ny, nz  !sizes of the array 
  integer :: kx, ky, kz,  k0x, k0y, k0z
  double precision , intent(in) :: radius  !the radius of the sphere
  ! Working precision is double precision
  integer, parameter :: WP = selected_real_kind(15,307)
  real(WP)   :: array(:,:,:)
  logical :: inside_radius 
  array=0.d0
  nx=size(array,1)
  ny=size(array,2)
  nz=size(array,3)
  !locate the center of the array
  k0x=nx/2+1 
  k0y=ny/2+1
  k0z=nz/2+1
!evaluate the Gaussian function at all the points on the whole array  
  forall( kx=1:nx, ky=1:ny  , kz=1:nz,  (  (kx-k0x)**2  + (ky-k0y)**2 +(kz-k0z)**2   ) .le. radius**2      )     
     array(kx,ky, kz)=1.d0  
  end forall
  !write(*,*) x_real
end subroutine init_sphere3D

subroutine color_by_user(positions )
    integer :: i,j,k
    double precision :: thres !threshold to apply    
    integer :: positions(:,:)  !stores the user that "owns" each gaze point     
    !the first dimension is the point number, the second one has the values of x,y and z 
    integer, allocatable :: distances2(:,:)
    integer :: numpoints
    
    if (allocated( conv_space_r ).eqv. .false.) then
        call copy_mask_to_real()    
    end if
    
    numpoints=size( positions,1 )    
    allocate(distances2(numpoints,3)  )
    
    !we set the theshold to the value of the Gaussian at 1 standard deviation
    thres=exp( -0.5 )
    !forall (i=1:L1, j=1:L2, k=1:L3 , conv_space_r >thres )
        !calc_dist2( positions, i,j,k )
    !end forall
    
    deallocate(distances2)
    
    write(*,*) "Color"
end subroutine color_by_user


subroutine mask_set_random_walk(steps)
    implicit none
    integer :: steps
    
end subroutine mask_set_random_walk





END MODULE fft_convolve_3d
