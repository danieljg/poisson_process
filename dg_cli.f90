! Fake Data Generator
! 5 de Agosto de 2012
module dg_module
use aux
implicit none
contains

 subroutine read_command_arguments
 implicit none
 integer :: k
  if(command_argument_count().ne.2)then
   write(*,*) 'The program must be run with two arguments "./dg_cli -nbar -vis"'
   write(*,*) 'or none at all (defaultsn: bar=4000,vis=0.9)'
   stop
  endif
  k=1
  k=read_real(k,Nbar)
  k=read_real(k,vis)
 end subroutine read_command_arguments

 subroutine initialize_variables
 implicit none
  phase  = pi
  Nbar   = 4000
  vis    = 0.9
  dt     = 1.3107e-3   ! time interval for a measurement
  if (command_argument_count().ne.0)then
   call read_command_arguments
  endif
 end subroutine initialize_variables

 integer function read_real(k,variable)
 implicit none
 integer :: k
 real    :: variable
 character(12) :: value
  call get_command_argument(k,value)
  read(value,*) variable
  read_real=k+1
 end function read_real
 
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
