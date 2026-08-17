!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine scan_harmonic_ser
	use mod_global_variables
	use mod_file_definition, only : pw_harmonic_ser, id_harmonic_ser
	implicit none
	
	character(len = 200) :: line
	integer :: j, k, i1, i2
	! jw

	! jw
	open(pw_harmonic_ser, file=id_harmonic_ser, form='formatted', status = 'old')

	! jw
	call skip_header_lines(pw_harmonic_ser,id_harmonic_ser)
	
	max_no_tidal_constituent = 1 ! jw
	do k = 1, num_harmonic_ser
		read(pw_harmonic_ser,*) i1, i2
		max_no_tidal_constituent = MAX(max_no_tidal_constituent,i2) ! jw
		do j=1,i2
			read(pw_harmonic_ser,*) line
		end do
	end do	

	close(pw_harmonic_ser)
end subroutine scan_harmonic_ser