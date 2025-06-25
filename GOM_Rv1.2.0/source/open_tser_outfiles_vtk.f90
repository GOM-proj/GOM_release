!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Tecplot & VTK files require different file formats at the header lines.
!! Thus, this subroutine will write different header lines for each file format.
!! 
subroutine open_tser_outfiles_vtk(pw_tser,id_tser,tser_list)
	use mod_global_variables
	use mod_file_definition
	implicit none
	
   integer, intent(in) :: pw_tser, tser_list
   character(len=200), intent(in) :: id_tser	
	! End of local variables ==================================================!
	
	! Tecplot outfiles uses ***.dat, but VKT files use ***.txt
	open(pw_tser, file=trim(id_tser), form='formatted', status='replace')

	! Tecplot requires "Variables =" but not in VTK file
   if(tser_time == 1) then ! [sec]
   	write(pw_tser,'(A)', advance = 'no') '"Elapsed Time [sec]"'
   else if(tser_time == 2) then ! [min]
   	write(pw_tser,'(A)', advance = 'no') '"Elapsed Time [min]"'
   else if(tser_time == 3) then ! [hr]
   	write(pw_tser,'(A)', advance = 'no') '"Elapsed Time [hr]"'
   else if(tser_time == 4) then ! [day]
   	write(pw_tser,'(A)', advance = 'no') '"Elapsed Time [day]"'
   end if
   
	! now call identical part =================================================!	
	call open_tser_outfiles_sub(pw_tser,id_tser,tser_list)
end subroutine open_tser_outfiles_vtk