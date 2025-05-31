!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine read_harmonic_ser
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, j
	integer :: harmonic_ser_id_dummy
	character(len=20) :: ctemp1
	! jw
	
	open(pw_harmonic_ser,file=id_harmonic_ser,form='formatted',status='old')

	! jw
	call skip_header_lines(pw_harmonic_ser, id_harmonic_ser)
	
	! jw
	do i=1,num_harmonic_ser
		read(pw_harmonic_ser,*) harmonic_ser_id_dummy, no_tidal_constituent(i), tidal_phase_shift(i)
		do j=1,no_tidal_constituent(i)
			read(pw_harmonic_ser,*) tidal_period(j,i), tidal_amplitude(j,i), tidal_phase(j,i), &
			&	tidal_nodal_factor(j,i), equilibrium_argument(j,i), ctemp1
		end do
	end do
end subroutine read_harmonic_ser