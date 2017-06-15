!
! Finite Well calculation calculation using matrix approach
! with harmonic oscillator basis set
! developed by Dr. Stoyan Pisov 2012
!
program finitewell
implicit none
include "fftw3.f"

double precision, parameter :: pi = 3.14159265358979323846264338327d0
! N - Size of basis set, should be odd number
integer, parameter :: n = 2048
! lwork - work variable required by ZHEEV subroutine
integer, parameter :: lwork = 3*n - 1
! Length - half of the distance of finite well
double precision, parameter :: length = 1.d0
! Define imaginary number
complex*16, parameter :: img = ( 0.d0, 1.d0)
! Define complex number one
complex*16, parameter :: one = ( 1.d0, 1.d0)


! X - Vector with discrete points
double precision, dimension(n) :: x

! HF - array containing first 0 - N harmonic oscillator basis
! function discretized over M points
! BWF - calculated basis function of finite well
complex*16, dimension(n,n) :: bwf, pbwf, sbwf
double precision, dimension(n,n) :: prbf

!Work array for calculating the potential energy matrix with FFT
complex*16, dimension(n) :: invec, outvec

! work variable for creating FFT plan
! plan - forward plan
! iplan - backward plan
integer*8 :: plan, iplan

! Epot - N x N potential energy matrix calculated in oscillator basis
! Ekin - N x N kinetic energy matrix calculated in oscillator basis
! Etot - N x N total energy matrix calculated in oscillator basis
complex*16, dimension(1:n,1:n) :: epot,ekin,etot
! Energies - Eigenvalue vector of Etot
double precision, dimension(1:n) :: energies, epotf
! work - array required by ZHEEV
double precision :: work(lwork), rwork(lwork)

!Potential energy level V0
double precision, parameter :: V0 = 10.d0

!work variables
double precision :: xmin, xmax, Delta, step
integer :: i, j, info, np, nbo
complex*16 :: integ

character*1 :: jobz, uplo

! We request eigenvector as well
jobz = 'V'
uplo = 'U'

! Define function range

Delta = 20.d0
xmin = 0.d0
xmax = Delta

! Calculate discretization step
step =  Delta /n

! Initialize coordinate vector X
do i = 1, n
  x(i) = (i-1) * step + xmin
end do

! Calculate potential energy profile finite well with limit equal V0
invec(:) = cmplx(V0, 0.d0)
do i = 1, n
  if(abs(x(i)-0.5*Delta).le.(length)) then
    invec(i) = (0.d0, 0.d0)
  end if
end do

epotf(:) = real(invec(:))

! Create FFT forward plan
call dfftw_plan_dft_1d( plan, n,invec, outvec,fftw_forward, FFTW_ESTIMATE)

! Create FFT backward plan
call dfftw_plan_dft_1d(iplan, n,outvec, invec,fftw_backward, FFTW_ESTIMATE)

! Execute the plan itself
call dfftw_execute(plan)

! Norm the array
outvec(:) = outvec(:) / n

do i = 1, n
  do j = i, n
    epot(i,j) = outvec(j - i + 1)
    if (i.ne.j) then
      epot(j, i) = dconjg(outvec(j - i + 1))
    end if
  end do
end do


ekin(:, :) = cmplx(0.d0, 0.d0)
do i = 1, n 
  ekin(i, i) = (one - 0.5d0*exp(2*img*pi*i/n) - 0.5d0*exp(-2*img*pi*i/n))*n**2/Delta**2
end do

! Calculate Epot, Ekin, Etot matrices
etot(:, :) = ekin(:, :) + epot(:, :)

! External LAPACK subroutine call return eigenvalues and eigenvectors of Etot
call zheev (jobz, uplo, n, etot, n, energies, work, lwork, rwork, info)

! Calculate finite well eigenfunctions
do i = 1, n
  outvec(:) = etot(:, i)
  call dfftw_execute(iplan)
  sbwf(:, i) = invec(:) / sqrt(Delta)
  prbf(:, i) = cdabs(sbwf(:, i))**2 + energies(i)
end do

! Print out the data
! Print first five eigenvalues to standard output
do i = 1, 5
  write(0,'(A1,I1,A3,F15.7)')'E',i,' = ',energies(i)
end do

! Print first four eigenvectors and potential energy profile
! to standard error output
do i = 1, n
  write(6,'(6F15.7)')x(i),prbf(i,1),prbf(i,2),prbf(i,3),prbf(i,4),epotf(i)
end do

end program finitewell
