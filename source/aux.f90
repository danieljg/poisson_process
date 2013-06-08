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
 subroutine estimate_phase(dataset,distribution,l,phaseval)
 implicit none
 integer, intent(in) :: l
 real, dimension(l)  &
     , intent(in)    :: dataset,distribution
 real, dimension(l)  :: dist_copy
 real, intent(out)   :: phaseval
 real, dimension(l)  :: crosscorr
 integer :: i,tt=1
 real    :: mc,md
  ! Crear una copia de la distribucion
  dist_copy = distribution
  ! Calcular funcion de correlacion cruzada
  crosscorr = ccf_dft(dataset, dist_copy)
  ! Encontrar la posicion del maximo
  phaseval = maxloc(crosscorr,1)
  !write(*,*)"phase",phase
  !max_loc  = 1
  if(tt.eq.1) then
   mc = maxval(crosscorr,1)
   md = maxval(dataset,1)
   open(15,file="ccf.dat")
   write(*,*)"expected_phase",phase*l/(2*pi)
   write(*,*)"estimated_phase",phaseval
   do i=1,l
    write(15,'(I4,3F6.2)')i,crosscorr(i)/mc,&
                dataset(i)/md,distribution(i)
   end do
   tt=0
  end if
 end subroutine estimate_phase

 function ccf_dft(array1,array2)
 implicit none
 real, dimension(l) :: array1, array2
 real, dimension(l) :: ccf_dft
 complex, dimension(l) :: carray1, carray2,&
                          farray1, farray2
 complex, dimension(l) :: cout1, cout2
 type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
 integer :: stat
  carray1=array1
  carray2=array2
  stat = DftiCreateDescriptor( My_Desc1_Handle,&
                               DFTI_SINGLE,&
                               DFTI_COMPLEX,&
                               1, l )
  stat = DftiCommitDescriptor( My_Desc1_Handle )
!  stat = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT,&
!                       DFTI_NOT_INPLACE)
  stat = DftiComputeForward(My_Desc1_Handle,&
                            carray1)!,farray1)
  stat = DftiComputeForward(My_Desc1_Handle,&
                            carray2)!,farray2)
!  farray1 = farray1*conjg(farray2)
  carray1=carray1*conjg(carray2)
  stat = DftiComputeBackward(My_Desc1_Handle,&
                             carray1)!farray1,carray1)
  ccf_dft = real(carray1)
  stat = DftiFreeDescriptor(My_Desc1_Handle)
 end function ccf_dft

 ! Calcular la correlacion cruzada
 function ccf(array1,array2)
 implicit none
 real, dimension(l) :: array1, array2
 real, dimension(l) :: ccf
 integer :: i
  ! Calcular la ccf
  do i=1,l
   ccf(i) = dot_product(array1,array2)
   ! shift the array
   array2 = shift(array2)
  end do
 end function ccf

 ! Deslizar el pulso periodico en un dt
 function shift(array)
 implicit none
 real :: keeper
 real, dimension(l), intent(in) :: array
 real, dimension(l) :: shift
  ! save first value
  keeper       = array(1)
  shift(1:l-1) = array(2:l)
  shift(l)     = keeper
 end function shift

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

 ! Calcular el vector gaussiano para la correlacion cruzada
 subroutine calculate_gaussian(gaussian,rad_width)
 implicit none
 integer :: k
 real, intent(in) :: rad_width
 real, intent(out) :: gaussian(l)
 real :: width
 width = rad_width*l/(2*pi)
 do k=1, l/2
  gaussian(k) = (exp( (k-1)/width))**2
 end do
 do k=l/2+1,l
  gaussian(k) = (exp( (k-l-1)/width))**2
 end do
 end subroutine calculate_gaussian

 ! Calcular momentos normalizados alrededor de la media
 subroutine calculate_phase_moments(n,k,vec,phase_moments)
 implicit none
 integer, intent(in) :: n,k
 real, intent(in)    :: vec(k)
 real, intent(out)   :: phase_moments(n)
 integer :: kk,jj
 real, dimension(k)  :: temp
 real :: mean, stddev
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
