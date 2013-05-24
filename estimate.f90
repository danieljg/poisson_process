! #########################################################
! ### Test suite for periodic poisson process simulator ###
! #########################################################

program estimate
use aux
implicit none
integer :: k            ! do loop counter
integer :: diff         ! difference of Nbar and Ncum
integer :: nn, nbin     ! # of realizations and bins
integer :: fact_moments ! # number of factorial moments
                        ! to be calculated
integer :: phase_moments! # of moments about the mean to
                        ! be calculated for the phase data
parameter(nn=1e5,nbin=1001, fact_moments=24, phase_moments=8)
real,  &
dimension(nbin) :: hist ! histogram of differences
real, &
dimension(fact_moments) :: moments ! factorial moments vector
real, &
dimension(nn)  :: phase_n
real, &
dimension(phase_moments) :: phasevec
real, &
dimension(l)   :: cosine

 ! initialize variables
 call init

 ! calculate cosine for correlation
 call calculate_cosine(cosine)

 ! initialize random number generator
 call random_seed

 ! initialize histogram recording array
 hist(1:nbin) = 0

 ! run data-generating routine nn times
 do k=1,nn
  ! simulate random process
  call make_data ! defines Nbar (expected) and Ncum (observed)
  ! update histogram 
   diff=Ncum-Nbar
   if( abs(diff) .gt. (nbin-1)/2 ) then 
    write(*,*)'this should not happen'
   else
    hist( (nbin-1)/2 + diff) = hist( (nbin-1)/2 + diff ) + 1
   end if
  ! estimate phase
  call estimate_phase(dataset,cosine,l,phase_n(k))
 end do

 ! save bin data to file
 ! call bin_data( )
 ! write data to file
 open(15, file="histogram.dat")
 write(15,*)"#"
 do k=1,nbin
  write(15,*) k-1 -(nbin-1)/2 + nbar , hist(k)
 end do
 close(15)
 ! normalize histogram
 hist(1:nbin) = hist(1:nbin)/nn

 ! process phase data
 ! call phase_data( )
 ! calculate moments
 call calculate_moments(phase_moments,nn,&
                              phase_n,phasevec)
 ! save phase data to file
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
  write(16,*) k,phasevec(k),(sqrt(phasevec(2)))**k * dfact(k-1)
 end do
 close(16)

 ! calculate the factorial moments
 ! and write data to file
 ! call fact_driver()
 call factorial_moments(fact_moments,nbin,hist,moments)
 open(15, file="fact_moments.dat")
  write(15,*) '# "r" ','"<n>^r" ',' "fact. moments"',' "% error (theo)"',' "% error (emp)"'
  do k=1,fact_moments
   write(15,'(I2,4E16.6E4)') k, nbar**(k), moments(k), 1-moments(k)/nbar**k, 1-moments(k)/moments(1)**k
  end do
 close(15)

contains

end program estimate
