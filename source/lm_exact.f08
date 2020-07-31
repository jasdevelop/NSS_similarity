!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   exact_module.f08
! Author: jorge
!
! Created on August 12, 2019, 7:04 PM
!

MODULE lm_exact
use input_module , ONLY : num_users, box_bounds,  gaze_data, sigma_array
implicit none
private
save
double precision, parameter  , public      :: pi=3.1415926535897932384d0
double precision, public, allocatable :: overlap(:,:)  , similarity(:,:), average_field(:) 

integer , allocatable , public :: valid_points(:)  !number of valid gaze points by each user 

public :: initialize_arrays_exact , setup_arrays_exact, reset_arrays_exact, NSS_calc_single, calc_norm_infinity

contains

subroutine initialize_arrays_exact(  )
    
allocate(  overlap(num_users, num_users  )  )
allocate(  similarity(num_users, num_users ) )
allocate(  average_field(num_users)  ) 
allocate(  valid_points(num_users) )

call reset_arrays_exact()

end subroutine initialize_arrays_exact

subroutine reset_arrays_exact()
    
    overlap=0.d0
    similarity=0.d0
    average_field=0.d0
    valid_points=0
    
end subroutine reset_arrays_exact

subroutine setup_arrays_exact( ti, tf, planar_mode, approx)
    implicit none
    integer :: t1, t2 
    integer, intent(in) :: ti, tf
    integer :: u1, u2
    integer :: k 
    double precision :: p, q
    double precision :: x1(3)  , x2(3) , xp(3), xm(3)
    double precision :: vdist
    double precision ::  simi_prod , sfac(3)      
    double precision ::  afac(3), afac_prod  !average_field partial sum and factors
    double precision ::   ofac(3), ofac_prod  !overlap factors and product 
    double precision :: v(3)  !array with the Gaussian width parameter nu (greek letter) 
    integer :: eff_dim 
    logical, intent(in) :: planar_mode
    logical, intent(in) :: approx
    
    if (planar_mode .eqv. .TRUE.) then
        !Keep the calculation to a plane
        eff_dim=2
    else
        !Use all three dimensions
        eff_dim=3
    end if
    
    x1=0.d0
    x2=0.d0
    call reset_arrays_exact()
    
    v=0.5d0/sigma_array**2
    
    valid_points=0
    
    do u1=1,num_users
        do u2=u1,num_users
                        
            do t1=ti, tf
                x1(:)=gaze_data(t1, u1, :   )
       
                             
                do t2=ti, tf
                    x2(:)=gaze_data(t2, u2, :)
                    xm=(x1-x2)
                    xp=( x1+x2 )/2.d0
                    
                    !similarity  calculation
                    
                    
                    vdist=0.d0
                    do k=1,eff_dim
                        vdist=vdist+ v(k)*(  xm(k)   )**2
                    end do
                    simi_prod=exp( -vdist   )
                    !sum(     )
                    
                    !write(99,*) simi_prod, u1, u2, t1, t2
                    if (.not. isnan(simi_prod) ) then
                        similarity(u1,u2)=similarity(u1,u2) +simi_prod
                        similarity(u2,u1)=similarity(u1,u2)
                    end if
                    
                    !overlap calculation
                    
                    ofac_prod=1.d0
                    do k=1,eff_dim
                        p=box_bounds(k,1)
                        q=box_bounds(k,2)
                        ofac_prod =ofac_prod* exp(-0.5*v(k)*xm(k)**2 )*integral_single( p,q, 2*v(k) , xp(k) , approx  )
                        !write(98,*) ofac(k), p,q, xp(k), k
                    end do
                    
                    if ( .not.  isnan( ofac_prod )    ) then
                        overlap(u1,u2)=overlap(u1,u2) +ofac_prod
                        overlap(u2,u1)=overlap(u1,u2)
                    end if
                    
                    
                end do  !end of loop for time of user 2
            end do    ! end of loop for time of user 1
            
        end do   !end of loop for user 2     
    end do  ! end of loop for user 1

    !Calculation for average field
    do u1=1,num_users
        do t1=ti, tf
            x1(:)=gaze_data(t1, u1, :   )

            ! average_field calculation
            afac_prod=1.d0
            do k=1,eff_dim
                p=box_bounds(k,1)
                q=box_bounds(k,2)
                afac_prod=afac_prod* integral_single(p,q,v(k), x1(k) , approx )
            end do
            
            if ( .not.  isnan( afac_prod )   ) then
                !Test if it has any NANs
                ! if valid then add in the contribution of the Gaussian to the average field
                average_field(u1)=average_field(u1)+ afac_prod
                !we also add 1 to the valid gaze point counter of this user
                valid_points(u1)=valid_points(u1)+1
            end if
            
        end do
        
    end do
    
end subroutine setup_arrays_exact


function integral_single( p,q , v , x0, approx )
    double precision :: p,q
    double precision :: x0
    double precision :: v
    logical , intent(in) :: approx 
    double precision :: r1, r2
    double precision, parameter        :: pi=3.1415926535897932384d0
    double precision :: integral_single
    r1=0.5d0*( pi/v )**0.5d0
    if (approx .eqv. .TRUE. ) then
        r2=2.d0
    else 
        r2=erf(  ( q-x0 )*(v )**0.5d0  ) - erf(  ( p-x0)*( v)**0.5d0 )    
    end if
    integral_single=r1*r2
end function integral_single




function NSS_calc_single(planar_mode, approx)
    double precision :: ofac(3) , afac(3)
    double precision :: p,q , v(3)
    double precision :: std , volume, sim, avg, variance
    integer :: k
    logical, intent(in) :: planar_mode
    integer :: eff_dim 
    logical, intent(in) :: approx 
    double precision NSS_calc_single
    
    v=0.5/sigma_array**2
    
    ofac=1.d0
    afac=1.d0
    
    if (planar_mode .eqv. .TRUE. ) then
        eff_dim=2        
    else
        eff_dim=3
    end if
    
    volume=1.d0    
    do k=1,eff_dim
        volume=volume*( box_bounds(k,2) -box_bounds(k,1)  )        
    end do
    
    !overlap calculation
    do k=1,eff_dim
        p=box_bounds(k,1)
        q=box_bounds(k,2)
        !we put the Gaussian in the center of the box
        ofac(k) =  integral_single( p,q, 2*v(k) , (p+q)/2 , approx )
        afac(k)  =  integral_single( p,q,   v(k) , (p+q)/2 , approx )
    end do    
    sim=1.0  ! it is similar to itself
    avg=PRODUCT( afac  )/volume   !The integral to infinity of all components
    variance=PRODUCT( ofac  )/volume -avg**2 
    std=variance**0.5
    write(*,* ) "Single Gaussian avg, std, var ", avg, std, variance    
    NSS_calc_single=(sim-avg)/std
end function NSS_calc_single

function calc_norm_infinity( w )
    double precision :: w
    double precision :: calc_norm_infinity
    calc_norm_infinity=( (2**0.5)*( w /pi)**0.5 )**0.5
end function calc_norm_infinity


END MODULE lm_exact
