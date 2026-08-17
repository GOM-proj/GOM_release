!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Tecplot & VTK files require different file formats at the header lines.
!! Thus, this subroutine will write different header lines for each file format.
!! 
subroutine open_tser_outfiles_tec(pw_tser,id_tser,tser_list)
	use mod_global_variables
	use mod_file_definition
	implicit none
	
   integer, intent(in) :: pw_tser, tser_list
   character(len=200), intent(in) :: id_tser
   
   character(len=200) :: tec_title
   ! jw

	! jw
	open(pw_tser, file=trim(id_tser), form='formatted', status='replace')
	
	! jw
	! jw
	if(tser_list == 1) then ! jw
		tec_title = 'Title = "Water elevation time series from msl'
	else if(tser_list == 2) then ! jw
		tec_title = 'Title = "Total water depth (H) time series'
	else if(tser_list == 3) then ! jw
		tec_title = 'Title = "Velocity 2d u time series'
	else if(tser_list == 4) then ! jw
		tec_title = 'Title = "Velocity 2d v time series'
	else if(tser_list == 5) then ! jw
		tec_title = 'Title = "Salinity time series'
	else if(tser_list == 6) then ! jw
		tec_title = 'Title = "Temperature time series'
	else if(tser_list == 7) then ! jw
	   tec_title = 'Title = "Air pressure time series'
	end if

	! jw
	if(tser_hloc == 1) then
		tec_title = trim(tec_title)//' at cell"'
	else if(tser_hloc == 2) then
		tec_title = trim(tec_title)//' at node"'
	end if
			
	! jw
	write(pw_tser,'(A)') trim(tec_title)
	! jw
	
	! jw
   if(tser_time == 1) then ! jw
   	write(pw_tser,'(A)', advance = 'no') 'Variables = "Elapsed Time [sec]"'
   else if(tser_time == 2) then ! jw
   	write(pw_tser,'(A)', advance = 'no') 'Variables = "Elapsed Time [min]"'
   else if(tser_time == 3) then ! jw
   	write(pw_tser,'(A)', advance = 'no') 'Variables = "Elapsed Time [hr]"'
   else if(tser_time == 4) then ! jw
   	write(pw_tser,'(A)', advance = 'no') 'Variables = "Elapsed Time [day]"'
   end if


	! jw
	call open_tser_outfiles_sub(pw_tser,id_tser,tser_list)
end subroutine open_tser_outfiles_tec