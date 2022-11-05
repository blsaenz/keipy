module kei_ocncommon

	use kei_parameters

	implicit none

  real, save :: hmixd(0:1)         ! storage arrays for extrapolations
  real, save :: Us(NZP1,NVEL ,0:1) !  ..      ..     ..  ..
  real, save :: Xs(NZP1,NSCLR,0:1) ! ..      ..     ..  ..
	integer, save :: old,new         ! extrapolation index for Us,Xs,hmixd

  real, save :: qsw_absorbed(NZ)

	! Common tridiagonal matrix factors (set in "init ocn")
  real, save :: tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix	

end module kei_ocncommon