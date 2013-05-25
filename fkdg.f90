! Fake Data Generator
! 5 de Agosto de 2012
module dg_module
use aux
implicit none
contains

 subroutine initialize_variables
 implicit none
  phase  = pi
  Nbar   = 4000
  vis    = 0.9
  dt     = 1.3107e-3   ! time interval for a measurement
 end subroutine initialize_variables

 subroutine write_data
 implicit none
 integer :: i
  open(15, file="fakedata.dat")
  do i =1, l
   write(15,*) i,dataset(i)
  enddo
  close(15)
  write(*,*)'mean N:    ',Nbar
  write(*,*)'observed N:',Ncum
 end subroutine write_data
end module dg_module

! #########################
! ### Programa principal ##
! #########################

program fkdg
use dg_module
implicit none
 call initialize_variables
 call random_seed
 call make_fringe_dataset
 call write_data
end program fkdg
