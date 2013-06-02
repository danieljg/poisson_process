program estimate
use aux
use estimate_module
implicit none
real :: vector(l)
 call initialize_variables
 call random_seed
 call calculate_histogram_and_phases
 call phase_moments_driver
 call factorial_moments_driver
end program estimate
