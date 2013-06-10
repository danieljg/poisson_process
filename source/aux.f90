! #########################
! ### Modulo auxiliar #####
! #########################

module aux
use MKL_DFTI
implicit none
 integer, parameter :: l=1024
 real,    parameter :: pi=4.0*atan(1.0)
 real, dimension(l) :: dataset ! Realizations of periodic Poisson process
 integer :: Ncum=0     ! Cumulative number of points
 real :: y             ! random number dummies
 real :: phase, vis, Nbar, Period, dt! Physical constants 

 contains

 ! Generar datos de acuerdo a un proceso periodico de Poisson
 subroutine make_fringe_dataset
 implicit none
 integer :: k          ! contador del ciclo
 real    :: ycum       ! tiempo acumulado
 real    :: lambda
  Period = l*dt        ! total experiment length
  ! reset contador de N
  Ncum   = 0
  ! Generar los l "bins" del proceso periodico de Poisson
  ! cada uno es estacionario
  do k=1,l
   ! calcular la intensidad por este "bin"
   lambda = &
   ( 1.0 + (vis*Period/(2.0*pi*dt))*&
    (sin(2.0*pi*(k+0.5)*dt/Period + phase)&
     - sin(2.0*pi*(k-0.5)*dt/Period + phase) )& 
   ) * (Nbar/Period)
   ! llenar el "bin" de puntos
   ycum = 0.0
   dataset(k) = 0
   do
    call random_number(y)
    y = -log(1-y)/lambda
    ycum = ycum+y
    if(ycum.gt.dt)exit
    dataset(k) = dataset(k)+1
    Ncum = Ncum+1
   end do
  end do
 end subroutine make_fringe_dataset

 ! We want to calculate the first rmax factorial moments
 ! of the empirical distribution given by hist(nbin).
 ! The output will be given by the vector f_moments(rmax)
 subroutine factorial_moments&
            (rmax,nbin,histogram,f_moments)
 implicit none
 integer, intent(in)  :: rmax,nbin
 real, dimension(nbin), intent(in):: histogram
 real, dimension(rmax), intent(out) :: f_moments
 integer :: k,l
 real , dimension(nbin,rmax) :: Ss
 ! The first row of Ss is a sequence from nmin to nmax
 Ss(1,1) = nbar - (nbin-1)/2
 do k=2,nbin
  Ss(k,1) = Ss(k-1,1) + 1
 end do
 ! The following rows are given by multiplying the
 ! previous row element by (n+1)
 do k=1,nbin
  do l=2,rmax
   Ss(k,l) = Ss(k,l-1) * ( nbar - (nbin-1)/2 + k - l)
  end do
 end do
 ! We now multiply the newly created matrix by
 ! the empirical histogram to get the factorial moments
 f_moments = Matmul(histogram,Ss)
 end subroutine factorial_moments

 ! Estimar fase y visibilidad
 subroutine estimate_phase(dataset,l,phaseval)
 implicit none
 integer, intent(in) :: l
 real, dimension(l)  &
     , intent(in)    :: dataset
 real, intent(out)   :: phaseval
 integer :: i,tt=1
  ! Calcular fase por standard quantum limit
  call std_quantum(dataset,phaseval)
  phaseval=nint(phaseval*l/(2*pi))
 end subroutine estimate_phase

 subroutine std_quantum(dataset,phase)
 implicit none
 real, dimension(l), intent(in) :: dataset
 real, intent(out)  :: phase
 real    :: cosine(l), sine(l), re, im
 integer :: j
 re = 0.0
 im = 0.0
 call calculate_cosine(cosine)
 call calculate_sine(sine)
 re = dot_product(dataset,cosine)
 im = dot_product(dataset,sine)
 if(atan2(im,re).ge.0)then
  phase = atan2(im,re)
 else
  phase = atan2(im,re)+2*pi
 endif
 end subroutine

 ! Calcular el vector coseno para la correlacion cruzada
 subroutine calculate_cosine(cosine)
 implicit none
 integer :: k
 real, intent(out), &
       dimension(l) :: cosine
  do k=1, l
   cosine(k) = cos(2.0*pi*(k-1)/real(l))
  end do
 end subroutine calculate_cosine

 subroutine calculate_sine(sine)
 implicit none
 integer :: k
 real, intent(out), &
       dimension(l) :: sine
  do k=1, l
   sine(k) = sin(2.0*pi*(k-1)/real(l))
  end do
 end subroutine calculate_sine

 ! Calcular momentos normalizados alrededor de la media
 subroutine calculate_phase_moments(n,k,vec,phase_moments,stddev)
 implicit none
 integer, intent(in) :: n,k
 real, intent(in)    :: vec(k)
 real, intent(out)   :: phase_moments(n),stddev
 integer :: kk,jj
 real, dimension(k)  :: temp
 real :: mean
 ! First moment about the origin
 mean = sum(vec)/k
 ! Calculate std deviation
 stddev = sqrt( sum((vec-mean)**2)/k )
 ! Moments about the origin
 do kk=1,n
  do jj=1,k
   temp(jj) = (vec(jj)-mean)**kk
  end do
  phase_moments(kk) = sum( temp ) / (k*stddev**kk)
 end do
 end subroutine calculate_phase_moments

 ! Calculates double factorial
 function dfact(n)
 integer ::n
 real    ::dfact
 integer :: kk
 dfact=n
  do
   if(n.lt.2) exit
   dfact = dfact*(n-2)
   n=n-2
  end do
 end function dfact

end module aux
