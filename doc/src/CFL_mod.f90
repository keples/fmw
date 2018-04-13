!> this module calculates the CFL number
  module CFL_mod
    use  Types_mod
    
    implicit none
    !everything is private unless otherwise stated
    private
    public :: fd1d_heat_explicit_cfl 

   contains
    !> calculates the CFL number
    !> \( \text{CFL} = \kappa\frac{\Delta t}{\Delta x^2} \)    
    subroutine fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, &
     x_max, cfl)

     implicit none
     
     !> calculate CFL number
     real (kind=dp),intent(out) :: cfl
     real (kind=dp) :: dx
     real (kind=dp) :: dt
     !> the heat constant \( \kappa \)
     real (kind=dp),intent(in) :: k
     !> uppper bound of t-axis
     real (kind=dp),intent(in) :: t_max
     !> lower bound of t-axis
     real (kind=dp),intent(in) :: t_min
     !> number of intervals in t-axis
     integer,intent(in) :: t_num
     !> upper bound of x-axis
     real(kind=dp), intent(in) :: x_max
     !> lower bound of x-axis
     real (kind=dp),intent(in) :: x_min
     !> number of intervals in x-axis
     integer,intent(in) :: x_num

     dx = (x_max-x_min)/real(x_num-1, kind=dp)
     dt = (t_max-t_min)/real(t_num-1, kind=dp)

     cfl = k*dt/dx/dx

     write (*, '(a)') ' '
     write (*, '(a,g14.6)') '  CFL stability criterion value = ', cfl

    end subroutine
  end module CFL_mod
