  module RHS_mod
    use  Types_mod
    implicit none
    !everything is private unless otherwise stated
    private
    public :: func

   contains
    function func(j, x) result (d)
     implicit none

     integer,intent(in) :: j
     real (kind=dp) :: d
     real (kind=dp), dimension(:), intent(in) :: x

     d = 0.0e+00_dp
    end function
   end module RHS_mod
