! #######################################################
! ### Test suite for periodic poisson process simulator #
! #######################################################

module estimate_module
use aux
implicit none
integer :: nn, nbin     ! # of realizations and bins
integer :: fact_moments ! # number of factorial moments
                        ! to be calculated
integer :: phase_moments! # of moments about the mean to
                      ! be calculated for the phase data
parameter(nn=4e3,nbin=l, fact_moments=20, phase_moments=8)
real :: phase_n(nn), phasevec(phase_moments), phase_hist(nbin)
real :: counts_hist(nbin), moments(fact_moments)
real :: phase_stddev
character(6) :: file_nbar, file_vis
character(40) ::filepath
contains

subroutine initialize_variables
implicit none
 phase  = pi
 dt     = 1.3107e-3   ! time interval for a measurement
 if(command_argument_count().eq.0)then
  write(*,*)'this should be run from the script, fool'
  write(*,*)"unless you now what i'm doing"
  write(*,*)'my notebook has the instructions'
  stop
 else
  call get_command_argument(1,file_nbar)
  call get_command_argument(2,file_vis)
  filepath='data/nbar_'//trim(file_nbar)//'/vis_'//trim(file_vis)//'/'
  read(file_nbar,*)Nbar
  read(file_vis,*)vis
  Nbar   = 10*(2**nbar)
  vis    = vis/100
 endif
end subroutine initialize_variables

subroutine int_to_char(tmp,n)
implicit none
 integer n
 character(len=6) :: tmp
 write(tmp,'(I6)')n
 tmp = adjustl(tmp)
end subroutine int_to_char

subroutine read_histogram(ggzz)
implicit none
integer :: dummy,k,ggzz
character(6) :: ccgz
 call int_to_char(ccgz,ggzz)
 call system('gunzip '//trim(filepath)//trim(ccgz)//'.gz -N')
 open(unit=15,file=trim(filepath)//'fakedata.dat')
 do k = 1, nbin
  read(15,*) dummy,dataset(k)
 end do
 Ncum=sum(dataset)
 close(15)
 call system('gzip '//trim(filepath)//'fakedata.dat')
 call system('mv '//trim(filepath)//'fakedata.dat.gz '//trim(filepath)//trim(ccgz)//'.gz')
end subroutine read_histogram

subroutine calculate_contrast_and_phase
implicit none
integer k, diff         ! difference of Nbar and Ncum
 ! initialize histogram recording array
 counts_hist(1:nbin) = 0
 phase_hist(1:nbin) = 0
 ! run data-reading routine nn times
 do k=1,nn
  call read_histogram(k) ! reads dataset and defines Ncum
  ! update histogram 
   diff=Ncum-Nbar
   if( abs(diff) .gt. (nbin-1)/2 ) then 
    write(*,*)'this should not happen',Ncum,Nbar
   else
    counts_hist((nbin-1)/2+diff) = counts_hist( (nbin-1)/2 + diff ) + 1
   end if
  call estimate_phase(dataset,l,phase_n(k))
  ! update phase histogram
  diff=phase_n(k)-(phase*nbin/(2*pi))
  phase_hist( nbin/2+diff ) = phase_hist(nbin/2+diff) + 1
 end do
 ! writeout
 call write_counts_bin_data
 call write_phase_bin_data
end subroutine calculate_contrast_and_phase

subroutine write_phase_bin_data
implicit none
integer k
 ! write to file
 open(15, file="phase_histogram.dat")
 write(15,*)"# error, relative frequency"
 do k=1,nbin
  write(15,*) k-nbin/2,phase_hist(k)/nn
 end do
 close(15) 
end subroutine write_phase_bin_data

subroutine write_counts_bin_data
implicit none
integer k
 ! normalize histogram
 counts_hist(1:nbin) = counts_hist(1:nbin)/nn
 ! write to file
 open(15, file="counts_histogram.dat")
 write(15,*)"#"
 do k=1,nbin
  write(15,*) k-1 -(nbin-1)/2 + nbar , counts_hist(k)
 end do
 close(15)
end subroutine write_counts_bin_data

subroutine calculate_phase_moments_driver
implicit none
 ! calculate moments
 call calculate_phase_moments(phase_moments,nn,&
                        phase_n,phasevec,phase_stddev)
 ! save phase data to file
 call write_phase_moments()
end subroutine calculate_phase_moments_driver

subroutine write_phase_moments
implicit none
integer k
 open(15,file="phase.dat")
 write(15,*) '# ','"k" ','" phase(k)" '
 do k=1,nn
  write(15,*)k,phase_n(k)
 enddo
 close(15)
 open(16,file="phase_moments.dat")
 write(16,*) '# ','"k" ','"moments about the mean" ',&
             '"gaussian moments" '
 do k=1,phase_moments
  write(16,*) k,phasevec(k),(sqrt(phasevec(2)))**k *dfact(k-1)
 end do
 close(16)
end subroutine write_phase_moments

subroutine calculate_factorial_moments_driver
implicit none
integer  :: k
 call factorial_moments(fact_moments,nbin,counts_hist,moments)
 open(15, file="fact_moments.dat")
  write(15,*) '# "r" ','"<n>^r" ',' "fact. moments"',' "% error (theo)"',' "% error (emp)"'
  do k=1,fact_moments
   write(15,'(I2,4E16.6E4)') k, nbar**(k), moments(k), 1-moments(k)/nbar**k, 1-moments(k)/moments(1)**k
  end do
 close(15)
end subroutine calculate_factorial_moments_driver

subroutine write_estimated_error
implicit none
 write(*,*)vis,Nbar,phase_stddev
end subroutine write_estimated_error

end module estimate_module
