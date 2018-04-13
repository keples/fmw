  module Solver_mod
    use  Types_mod
    use  RHS_mod
    implicit none
    !everything is private unless otherwise stated
    private
    public :: fd1d_heat_explicit

   contains
    subroutine fd1d_heat_explicit(x, t, dt, cfl, h, h_new)
     implicit none

     integer :: x_num
     real (kind=dp), intent(in) :: cfl
     real (kind=dp), intent(in) :: dt
     real (kind=dp), dimension(:), intent(in) :: h
     real (kind=dp), dimension(:), intent(out) :: h_new
     integer :: j
     real (kind=dp), intent(in) :: t
     real (kind=dp), dimension(:), intent(in) :: x
     real (kind=dp) :: f(size(x))

     x_num = size(x)

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
   end module Solver_mod
