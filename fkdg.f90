! Fake Data Generator
! 5 de Agosto de 2012

! #########################
! ### Programa principal ##
! #########################

program fkdg
use aux
implicit none

 call random_seed

 call make_data

 call write_data

end program fkdg


