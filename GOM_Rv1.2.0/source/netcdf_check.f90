!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine netcdf_check(netcdf_status)
	use netcdf
	implicit none

	integer, intent(in) :: netcdf_status
	! jw

	if (netcdf_status /= nf90_noerr) then
		print *, "NetCDF ERROR: ", nf90_strerror(netcdf_status)
		stop "NetCDF Error"
	end if
end subroutine netcdf_check