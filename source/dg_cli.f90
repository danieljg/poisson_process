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
  Nbar=10*(2**Nbar)
  k=read_real(k,vis)
  vis=vis/100
 end subroutine read_command_arguments

 subroutine initialize_variables
 implicit none
 integer :: date_time(8), time_seed(2), M=2
 character (len = 12) :: real_clock(3)
 real :: ran
  phase  = pi
  Nbar   = 4000
  vis    = 0.9
  dt     = 1.3107e-3   ! time interval for a measurement
  if (command_argument_count().ne.0)then
   call read_command_arguments
  endif
  CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
                      REAL_CLOCK (3), DATE_TIME)

 call date_and_time(real_clock(1), real_clock(2), &
                    real_clock(3), date_time)
 time_seed(1) = date_time(6)+date_time(7)+date_time(8)
 time_seed(2) = date_time(1)+date_time(2)+date_time(8)
 call random_seed(size=M)
 call random_seed(put=time_seed(1:M))
 call random_number(ran)
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
 integer :: i,j
  open(15, file="fakedata.dat")
  do i =1, l
   j=nint(dataset(i))
   write(15,*) i,j
  enddo
  close(15)
 end subroutine write_data
end module dg_module

! #########################
! ### Programa principal ##
! #########################

program fkdg
use dg_module
implicit none
 call initialize_variables
 call make_fringe_dataset
 call write_data
end program fkdg
