program estimate
use aux
use estimate_module
implicit none
real :: vector(l)
 call initialize_variables
 call calculate_contrast_and_phase
 call calculate_phase_moments_driver
 !call calculate_factorial_moments_driver
 call write_estimated_error
end program estimate
