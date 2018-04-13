  module IO_mod
    use  Types_mod
    use  netcdf
    implicit none
    !everything is private unless otherwise stated
    private
    public :: r8mat_write,r8vec_linspace,r8vec_write 

   contains
    subroutine r8mat_write(x, t, output_filename, table)
     implicit none

     integer :: m
     integer :: n

     integer :: j
     real (kind=dp), dimension(:), intent(in) :: x
     real (kind=dp), dimension(:), intent(in) :: t
     character (len=*), intent(in) :: output_filename
     integer :: output_unit_id
     character (len=30) :: string
     real (kind=dp), dimension(:,:), intent(inout) :: table

     integer :: ierr, ncid, x_dimid, t_dimid, var_x_id, var_t_id, var_sol_id

     m = size(table(:,:), 1)
     n = size(table(:,:), 2)

     ! output_unit_id = 10
     ! open (unit=output_unit_id, file=output_filename, status='replace')

     ! write (string, '(a1,i8,a1,i8,a1,i8,a1)') '(', m, 'g', 24, '.', 16, ')'

     ! do j = 1, n
      ! write (output_unit_id, string) table(1:m, j)
     ! end do

     ! close (unit=output_unit_id)

     ! write NetCDF files
     ierr = NF90_CREATE(output_filename, NF90_CLOBBER, ncid)
     ierr = NF90_DEF_DIM(ncid, "x", size(x), x_dimid)
     ierr = NF90_DEF_DIM(ncid, "t", size(t), t_dimid) 
     ierr = NF90_DEF_VAR(ncid, "x-range", NF90_REAL8, [x_dimid], var_x_id)
     ierr = NF90_DEF_VAR(ncid, "t-range", NF90_REAL8, [t_dimid], var_t_id)
     ierr = NF90_DEF_VAR(ncid, "solution", NF90_REAL8,[x_dimid, t_dimid], var_sol_id)
     ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, "purpose", "Fortran workshop")
     ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, "name", "Zhong-Nan Wang")
     ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, "insitution", "University of Cambridge")
     ierr = NF90_PUT_ATT(ncid, var_x_id, "units", "metres")
     ierr = NF90_PUT_ATT(ncid, var_t_id, "units", "seconds")
     ierr = NF90_PUT_ATT(ncid, var_sol_id, "units", "Celsius")
     ierr = NF90_ENDDEF(ncid)
     ierr = NF90_PUT_VAR(ncid, var_x_id, x)
     ierr = NF90_PUT_VAR(ncid, var_t_id, t)
     ierr = NF90_PUT_VAR(ncid, var_sol_id, table)
     ierr = NF90_CLOSE(ncid)

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

   end module IO_mod
