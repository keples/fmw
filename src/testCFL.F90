!@test
subroutine testCFL()
  use pFUnit_mod
  use CFL_mod
  use Types_mod

  integer, parameter :: t_num = 201
  integer, parameter :: x_num = 21
  real(kind=dp) :: k, x_min, x_max, t_min, t_max
  real(kind=dp) :: cfl, cfl_exact, tol

  tol = 0.0000001_DP
  cfl_exact = 0.32_DP
  k = 0.002_DP

  x_min = 0.0_DP
  x_max = 1.0_DP
  t_min = 0.0_DP
  t_max = 80.0_DP

  call fd1d_heat_explicit_cfl( k, t_num, t_min, t_max, &
                               x_num, x_min, x_max, cfl )
#line 23 "testCFL.pf"
  call assertEqual( cfl, cfl_exact, tol, &
 & location=SourceLocation( &
 & 'testCFL.pf', &
 & 23) )
  if (anyExceptions()) return
#line 24 "testCFL.pf"
end subroutine testCFL

module WraptestCFL
   use pFUnit_mod
   implicit none
   private

contains


end module WraptestCFL

function testCFL_suite() result(suite)
   use pFUnit_mod
   use WraptestCFL
   type (TestSuite) :: suite

   external testCFL


   suite = newTestSuite('testCFL_suite')

   call suite%addTest(newTestMethod('testCFL', testCFL))


end function testCFL_suite

