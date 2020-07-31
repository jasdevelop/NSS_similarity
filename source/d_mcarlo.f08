!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   d_mcarlo.f08
! Author: jorge
!
! Created on January 9, 2020, 7:28 AM
!
include 'mkl_vsl.f90'
MODULE d_mcarlo
    !implicit none
    !private
    !save
    
    
    USE MKL_VSL_TYPE
    USE MKL_VSL
    implicit none
    private
    save
    TYPE (VSL_STREAM_STATE) :: stream

    integer :: brng,method,seed
    integer, parameter :: r_buffer_size=10000
    double precision :: r_buffer(r_buffer_size)
    
    public :: mc_randomizer, mc_randomizer_g, initialize_random_mkl, reset_random_mkl, Shuffle_buf
    
    contains
    
    subroutine mc_randomizer()
        integer :: i
        
        do i=1,10
            write(*,*) i
        end do
        
    end subroutine
    
    subroutine initialize_random_mkl(seed_new)
        implicit none
        integer, intent(in)  :: seed_new
        integer(kind=4) errcode    
    
        brng=VSL_BRNG_MT19937
        method=VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
        seed=seed_new

        !     ***** Initializing *****
        errcode=vslnewstream( stream, brng,  seed )          
        if (errcode .ne. VSL_STATUS_OK ) then            
            call failure_sub_mc( errcode, "Error with random generator initialization" )
        end if

        
    end subroutine initialize_random_mkl
    
    
    subroutine close_random_mkl()
        implicit none
        integer(kind=4) errcode
        !     ***** Deinitialize *****
        errcode=vsldeletestream( stream )
    end subroutine close_random_mkl
    
    subroutine reset_random_mkl( new_seed )
        implicit none
        integer, intent(in) :: new_seed
        
        call close_random_mkl()
        call initialize_random_mkl(new_seed)
        
    end subroutine reset_random_mkl
    
    subroutine random_number_mkl(r,n)
        implicit none
        integer(kind=4) :: errcode
        integer(kind=8) :: n  ! number of random numbers to generate
        real(kind=8) :: a, b ! parameters of distribution
        real(kind=8) , intent(inout):: r(:)  ! buffer for random numbers
        
        a=0.d0
        b=1.d0        
        
        errcode=vdrnguniform( method, stream, n, r, a, b )
        if (errcode .ne. VSL_STATUS_OK ) then            
            call failure_sub_mc( errcode, "Error with random generator" )
        end if
        
        
    end subroutine random_number_mkl
    
    subroutine mc_randomizer_g()
        implicit none
        integer, parameter :: q=6
        double precision :: array(q)
        integer :: array_int(q)
        integer :: n
        integer :: i,j
        n=size(array,1)
        array=0.d0
        !call initialize_random_mkl
        call random_number_mkl(array,n)
        
        !write(*,*) array
        
        do j=1,5
            do i=1,q
                array_int(i)=i
            end do
            call Shuffle(array_int)
            write(*,*) array_int
        end do
        !call close_random_mkl
 
    end subroutine mc_randomizer_g    
    
      
    subroutine Shuffle(a)
        integer, intent(inout) :: a(:)
        integer :: i, randpos, temp        
        double precision :: r(1)
        integer :: n
        n=size(a)
        do i = n, 2, -1          
          call random_number_mkl(r,1)
          randpos = int(r(1) * i) + 1
          temp = a(randpos)
          a(randpos) = a(i)
          a(i) = temp
        end do

    end subroutine Shuffle

    ! Adapted from  https://www.rosettacode.org/wiki/Knuth_shuffle#Fortran

    
    subroutine Shuffle_buf(a)
        integer, intent(inout) :: a(:)
        integer :: i, randpos, temp             
        integer :: n
        n=size(a)
        !prepares a buffer with random numbers for the entire shuffle 
        call random_number_mkl(r_buffer,n)
        
        !Shuffle using the buffered random values
        do i = n, 2, -1           
          randpos = int(r_buffer(i) * i) + 1
          temp = a(randpos)
          a(randpos) = a(i)
          a(i) = temp
        end do

    end subroutine Shuffle_buf
    
      
    subroutine failure_sub_mc( err_code, err_message )
        ! gives an error message an halts execution
        ! it is up to the user to decide what yo give as en error code
        implicit none
        integer(kind=4) , intent(in):: err_code   !error code from the routine that originated the issue
        character(len=*) :: err_message    !error message
        print '("  Error, status = ",I0)', err_code
        print *, trim(err_message)
        print *, "Terminating"
        call exit(1)
    
    end subroutine failure_sub_mc
    
    
    
END MODULE d_mcarlo
