program fd1d_heat_explicit_prb
 use :: types_mod, only: dp

 implicit none

 integer,parameter :: T_NUM = 201
 integer,parameter :: X_NUM = 21

 real (kind=dp) :: cfl
 real (kind=dp) :: dt
! real (kind=dp) :: h(x_num)
! real (kind=dp) :: h_new(x_num)
 real (kind=dp), dimension(:), allocatable :: h, h_new
! the "matrix" stores all x-values for all t-values
! remember Fortran is column major, meaning that rows are contiguous
! real (kind=dp) :: hmat(x_num, t_num)
 real (kind=dp), dimension(:,:), allocatable :: hmat
 integer :: i
 integer :: j
 real (kind=dp) :: k

 !real (kind=dp) :: t(t_num)
 real (kind=dp) :: t_max
 real (kind=dp) :: t_min
 !real (kind=dp) :: x(x_num)
 real (kind=dp) :: x_max
 real (kind=dp) :: x_min
 real (kind=dp), dimension(:), allocatable :: t, x

 write (*, '(a)') ' '
 write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
 write (*, '(a)') '  FORTRAN77 version.'
 write (*, '(a)') '  Test the FD1D_HEAT_EXPLICIT library.'

 write (*, '(a)') ' '
 write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
 write (*, '(a)') '  Normal end of execution.'
 write (*, '(a)') ' '

 write (*, '(a)') ' '
 write (*, '(a)') 'FD1D_HEAT_EXPLICIT_TEST01:'
 write (*, '(a)') '  Compute an approximate solution to the time-dependent'
 write (*, '(a)') '  one dimensional heat equation:'
 write (*, '(a)') ' '
 write (*, '(a)') '    dH/dt - K * d2H/dx2 = f(x,t)'
 write (*, '(a)') ' '
 write (*, '(a)') '  Run a simple test case.'

! allocate variables
 allocate( h(1:X_NUM) )
 allocate( h_new(1:X_NUM) )
 allocate( hmat(1:X_NUM,1:T_NUM) )
 allocate( t(1:T_NUM) )
 allocate( x(1:X_NUM) )

! heat coefficient
 k = 0.002e+00_dp

! the x-range values
 x_min = 0.0e+00_dp
 x_max = 1.0e+00_dp
! x_num is the number of intervals in the x-direction
 call r8vec_linspace(x_min, x_max, x)

! the t-range values. integrate from t_min to t_max
 t_min = 0.0e+00_dp
 t_max = 80.0e+00_dp

! t_num is the number of intervals in the t-direction
 dt = (t_max-t_min)/real(t_num-1, kind=dp)
 call r8vec_linspace(t_min, t_max, t)

! get the CFL coefficient
 call fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, x_max, cfl)

 if (0.5e+00_dp<=cfl) then
  write (*, '(a)') ' '
  write (*, '(a)') 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
  write (*, '(a)') '  CFL condition failed.'
  write (*, '(a)') '  0.5 <= K * dT / dX / dX = CFL.'
  stop
 end if

! set the initial condition
 do j = 1, x_num
  h(j) = 50.0e+00_dp
 end do

! set the bounday condition
 h(1) = 90.0e+00_dp
 h(x_num) = 70.0e+00_dp

! initialise the matrix to the initial condition
 do i = 1, x_num
  hmat(i, 1) = h(i)
 end do

! the main time integration loop 
 do j = 2, t_num
  call fd1d_heat_explicit(x, t(j-1), dt, cfl, h, h_new)

  do i = 1, x_num
   hmat(i, j) = h_new(i)
   h(i) = h_new(i)
  end do
 end do

! write data to files
 call r8mat_write('h_test01.txt', hmat)
 call r8vec_write('t_test01.txt', t)
 call r8vec_write('x_test01.txt', x)

 deallocate( h, h_new, hmat, t, x)

contains

 function func(j, x) result (d)
  implicit none

  integer,intent(in) :: j
  real (kind=dp) :: d
  real (kind=dp), dimension(:), intent(in) :: x

  d = 0.0e+00_dp
 end function

 subroutine fd1d_heat_explicit(x, t, dt, cfl, h, h_new)
  implicit none

!  integer,intent(in) :: x_num

  real (kind=dp), intent(in) :: cfl
  real (kind=dp), intent(in) :: dt
  real (kind=dp), dimension(:), intent(in) :: h
  real (kind=dp), dimension(:), intent(out) :: h_new
  integer :: j
  real (kind=dp), intent(in) :: t
  real (kind=dp), dimension(:), intent(in) :: x
  real (kind=dp) :: f(size(x))

  do j = 1, x_num
   f(j) = func(j, x)
  end do

  h_new(1) = 0.0e+00_dp

  do j = 2, x_num - 1
   h_new(j) = h(j) + dt*f(j) + cfl*(h(j-1)-2.0e+00_dp*h(j)+h(j+1))
  end do

! set the boundary conditions again
  h_new(1) = 90.0e+00_dp
  h_new(x_num) = 70.0e+00_dp
 end subroutine

 subroutine fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, &
   x_max, cfl)

  implicit none

  real (kind=dp),intent(out) :: cfl
  real (kind=dp) :: dx
  real (kind=dp) :: dt
  real (kind=dp),intent(in) :: k
  real (kind=dp),intent(in) :: t_max
  real (kind=dp),intent(in) :: t_min
  integer,intent(in) :: t_num
  real(kind=dp), intent(in) :: x_max
  real (kind=dp),intent(in) :: x_min
  integer,intent(in) :: x_num

  dx = (x_max-x_min)/real(x_num-1, kind=dp)
  dt = (t_max-t_min)/real(t_num-1, kind=dp)

  cfl = k*dt/dx/dx

  write (*, '(a)') ' '
  write (*, '(a,g14.6)') '  CFL stability criterion value = ', cfl

 end subroutine

 subroutine r8mat_write(output_filename, table)
  implicit none

  integer :: m
  integer :: n

  integer :: j
  character (len=*), intent(in) :: output_filename
  integer :: output_unit_id
  character (len=30) :: string
  real (kind=dp), dimension(:,:), intent(inout) :: table

  m = size(table(:,:), 1)
  n = size(table(:,:), 2)

  output_unit_id = 10
  open (unit=output_unit_id, file=output_filename, status='replace')

  write (string, '(a1,i8,a1,i8,a1,i8,a1)') '(', m, 'g', 24, '.', 16, ')'

  do j = 1, n
   write (output_unit_id, string) table(1:m, j)
  end do

  close (unit=output_unit_id)
 end subroutine

 subroutine r8vec_linspace(a_first, a_last, a)

  implicit none

  integer :: n
  real (kind=dp), dimension(:), intent(inout) :: a
  real (kind=dp), intent(in) :: a_first
  real (kind=dp), intent(in) :: a_last
  integer :: i

  n=size(a)

  do i = 1, n
   a(i) = (real(n-i,kind=dp)*a_first+real(i-1,kind=dp)*a_last)/ &
     real(n-1, kind=dp)
  end do

 end subroutine

 subroutine r8vec_write(output_filename, x)

  implicit none

  integer :: m
  integer :: n

  integer :: j
  character (len=*),intent(in) :: output_filename
  integer :: output_unit_id
  real (kind=dp), dimension(:), intent(inout) :: x

  n = size(x)

  output_unit_id = 11
  open (unit=output_unit_id, file=output_filename, status='replace')

  do j = 1, n
   write (output_unit_id, '(2x,g24.16)') x(j)
  end do

  close (unit=output_unit_id)
 end subroutine

end program
