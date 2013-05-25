program estimate
use aux
use estimate_module
implicit none
 call initialize_variables
 call calculate_cosine(cosine)
 call random_seed
 call calculate_histogram_and_phases
 call phase_moments_driver
 call factorial_moments_driver
end program estimate
