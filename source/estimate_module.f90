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
parameter(nn=1e5,nbin=l, fact_moments=20, phase_moments=8)
real :: phase_n(nn), phasevec(phase_moments)
real :: hist(nbin), moments(fact_moments)
contains

subroutine initialize_variables
implicit none
 phase  = pi
 Nbar   = 400
 vis    = 0.7
 dt     = 1.3107e-3   ! time interval for a measurement
end subroutine initialize_variables

subroutine calculate_histogram_and_phases
implicit none
integer k, diff         ! difference of Nbar and Ncum
 ! initialize histogram recording array
 hist(1:nbin) = 0
 ! run data-generating routine nn times
 do k=1,nn
  ! simulate random process
  call make_fringe_dataset ! returns dataset o defines Nbar (expected) and Ncum (observed)
  ! update histogram 
   diff=Ncum-Nbar
   if( abs(diff) .gt. (nbin-1)/2 ) then 
    write(*,*)'this should not happen'
   else
    hist((nbin-1)/2+diff) = hist( (nbin-1)/2 + diff ) + 1
   end if
  ! estimate phase
  call estimate_phase(dataset,l,phase_n(k))
 end do
 ! writeout
 call write_bin_data
end subroutine calculate_histogram_and_phases

subroutine write_bin_data
implicit none
integer k
 open(15, file="histogram.dat")
 write(15,*)"#"
 do k=1,nbin
  write(15,*) k-1 -(nbin-1)/2 + nbar , hist(k)
 end do
 close(15)
 ! normalize histogram
 hist(1:nbin) = hist(1:nbin)/nn
end subroutine write_bin_data

subroutine phase_moments_driver
implicit none
 ! calculate moments
 call calculate_moments(phase_moments,nn,&
                        phase_n,phasevec)
 ! save phase data to file
 call write_phase_moments()
end subroutine phase_moments_driver

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

subroutine factorial_moments_driver
implicit none
integer  :: k
 call factorial_moments(fact_moments,nbin,hist,moments)
 open(15, file="fact_moments.dat")
  write(15,*) '# "r" ','"<n>^r" ',' "fact. moments"',' "% error (theo)"',' "% error (emp)"'
  do k=1,fact_moments
   write(15,'(I2,4E16.6E4)') k, nbar**(k), moments(k), 1-moments(k)/nbar**k, 1-moments(k)/moments(1)**k
  end do
 close(15)
end subroutine

end module estimate_module
