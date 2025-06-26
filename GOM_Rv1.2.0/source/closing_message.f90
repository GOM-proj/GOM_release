!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine closing_message
	use mod_global_variables
	use mod_file_definition	
	implicit none
	
	integer :: i
	logical :: fopened
	
	! jw
	character(len=8)	:: date
	character(len=10)	:: time
	character(len=5)	:: zone
	integer,dimension(8) :: values
	real(dp):: total_elapsed_time
	integer :: hh,mm,ss
	! jw
	

	! jw
	! jw
	call date_and_time(date,time,zone,values)	! jw
	call cpu_time(finish)							! jw
	
	gom_finish_year 	= values(1)
	gom_finish_month 	= values(2)
	gom_finish_day 	= values(3)
	gom_finish_hour 	= values(5)
	gom_finish_minute = values(6)
	gom_finish_second = values(7)
	
	total_elapsed_time = (gom_finish_day - gom_start_day) * 86400.0		&
	&						 + (gom_finish_hour - gom_start_hour) * 3600.0		&
	&						 + (gom_finish_minute - gom_start_minute) * 60.0	&
	&						 + (gom_finish_second - gom_start_second)
	
	hh = INT(total_elapsed_time/3600.0_dp)
	mm = INT((total_elapsed_time - hh*3600.0_dp)/60.0_dp)
	ss = INT(total_elapsed_time - hh*3600.0_dp - mm*60.0_dp)

	!=== Write an end message on 'run.log' file ===============================!
	write(pw_run_log,*) "!==============================================================================!"
	write(pw_run_log,*) 'Good Job ~!'
	write(pw_run_log,*) 'GOM reaches the normal termination.'
	write(pw_run_log,*) " "
	write(pw_run_log,'(A35,I4,A1,I2,A1,I2, A3, I2,A1,I2,A1,I2)') &
	&                "Today is (yyyy:mm:dd - hh:mm:ss): ", &
	&                 values(1),":",values(2),":",values(3)," - ",values(5),":",values(6),":",values(7) 
	write(pw_run_log,'(A22,F10.0,A10)') &
	&						"Total cpu time     = ", finish - start,"seconds" 		! jw
	write(pw_run_log,'(A22,F10.0,A10, A1,I2,A1,I2,A1,I2,A1)') &
	&						"Total elapsed time = ", total_elapsed_time,"seconds", '[',hh,':',mm,':',ss,']' ! jw
	write(pw_run_log,*) " "
	write(pw_run_log,*) 'The End.'
	write(pw_run_log,*) "!==============================================================================!"
	
	!=== Write a welcome message on the screen ================================!
	write(*,*) "!==============================================================================!"
	write(*,*)'Good Job ~!'
	write(*,*)'GOM reaches the normal termination.'
	write(*,*) " "
	write(*,'(A35,I4,A1,I2,A1,I2, A3, I2,A1,I2,A1,I2)') &
	&                "Today is (yyyy:mm:dd - hh:mm:ss): ", &
	&                 values(1),":",values(2),":",values(3)," - ",values(5),":",values(6),":",values(7) 
	
	write(*,'(A22,F10.0,A10)') &
	&						"Total cpu time     = ", finish - start,"seconds" 		! jw
	write(*,'(A22,F10.0,A10, A1,I2,A1,I2,A1,I2,A1)') &
	&						"Total elapsed time = ", total_elapsed_time,"seconds", '[',hh,':',mm,':',ss,']' ! jw
	write(*,*) " "
	write(*,*) 'The End.'
	write(*,*) "!==============================================================================!"

	!=== close all opened files ===============================================!
	! jw
	do i=11,999
		inquire(unit=i,opened=fopened)
		if(fopened) then
			close(unit=i)
		end if
	end do	
	!=== end of closing files =================================================!	
end subroutine closing_message