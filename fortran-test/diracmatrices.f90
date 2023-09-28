subroutine diracmatrices(xmn_in)
	use onebody_currents_sf
	IMPLICIT NONE
	real*8 :: xmn_in
	
	call dirac_matrices_in(xmn_in)

	return

end subroutine diracmatrices
