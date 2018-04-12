  module CFL_mod
    use  Types_mod
    
    implicit none
    !everything is private unless otherwise stated
    private
    public :: fd1d_heat_explicit_cfl 

   contains
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
   end module CFL_mod
