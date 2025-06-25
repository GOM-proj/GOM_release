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
   ! End of local variables ==================================================!

	! Tecplot outfiles uses ***.dat, but VKT files use ***.txt
	open(pw_tser, file=trim(id_tser), form='formatted', status='replace')
	
	! Only tecplot files need this title line =================================!
	! write tecplot title =====================================================!
	if(tser_list == 1) then ! tser_eta
		tec_title = 'Title = "Water elevation time series from msl'
	else if(tser_list == 2) then ! tser_H
		tec_title = 'Title = "Total water depth (H) time series'
	else if(tser_list == 3) then ! tser_u
		tec_title = 'Title = "Velocity 2d u time series'
	else if(tser_list == 4) then ! tser_v
		tec_title = 'Title = "Velocity 2d v time series'
	else if(tser_list == 5) then ! tser_salt
		tec_title = 'Title = "Salinity time series'
	else if(tser_list == 6) then ! tser_temp
		tec_title = 'Title = "Temperature time series'
	else if(tser_list == 7) then ! tser_airp
	   tec_title = 'Title = "Air pressure time series'
	end if

	! add variable horizontal location at the end of the tecplot title
	if(tser_hloc == 1) then
		tec_title = trim(tec_title)//' at cell"'
	else if(tser_hloc == 2) then
		tec_title = trim(tec_title)//' at node"'
	end if
			
	! write title for tecplot file
	write(pw_tser,'(A)') trim(tec_title)
	! =========================================================================!
	
	! Tecplot requires "Variables =" but not in VTK file
   if(tser_time == 1) then ! [sec]
   	write(pw_tser,'(A)', advance = 'no') 'Variables = "Elapsed Time [sec]"'
   else if(tser_time == 2) then ! [min]
   	write(pw_tser,'(A)', advance = 'no') 'Variables = "Elapsed Time [min]"'
   else if(tser_time == 3) then ! [hr]
   	write(pw_tser,'(A)', advance = 'no') 'Variables = "Elapsed Time [hr]"'
   else if(tser_time == 4) then ! [day]
   	write(pw_tser,'(A)', advance = 'no') 'Variables = "Elapsed Time [day]"'
   end if


	! now call identical part =================================================!	
	call open_tser_outfiles_sub(pw_tser,id_tser,tser_list)
end subroutine open_tser_outfiles_tec