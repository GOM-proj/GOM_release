!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine welcome_message
	use mod_global_variables
	use mod_file_definition	
	implicit none

	! jw
	character(len=8)	:: date
	character(len=10)	:: time
	character(len=5)	:: zone
	integer,dimension(8) :: values
	! jw
	
	! jw
	call date_and_time(date,time,zone,values)	! jw
	call cpu_time(start)								! jw
	gom_start_year 	= values(1)
	gom_start_month 	= values(2)
	gom_start_day 		= values(3)
	gom_start_hour 	= values(5)
	gom_start_minute 	= values(6)
	gom_start_second 	= values(7)
	
	!==== Create 'run.log' file for writing log information. ==================!
	open(pw_run_log, file = id_run_log, status='replace')
	write(pw_run_log,*) "!==============================================================================!"
	write(pw_run_log,*) "Step 0: Welcome Message!"
	write(pw_run_log,*) "Welcome to GOM."
	write(pw_run_log,*) "GOM is developed by Jungwoo Lee & Jun Lee."
	write(pw_run_log,*) "Then, good luck!"
	write(pw_run_log,*) " "
	write(pw_run_log,'(A35,I4,A1,I4,A1,I2, A3, I2,A1,I2,A1,I2)') 					&
	&                "Today is (yyyy:mm:dd - hh:mm:ss): ", 							&
	&                 values(1),":",values(2),":",values(3)," - ",values(5),":",values(6),":",values(7) 
	
	write(pw_run_log,*) "!==============================================================================!"
	write(pw_run_log,*) "Now, GOM is scanning main input file - main.inp"
	write(pw_run_log,*) "."
	write(pw_run_log,*) "."
	write(pw_run_log,*) "."	
	!=== End of writing welcome message on the log file =======================!
	
	!=== Write a welcome message on the screen ================================!
	write(*,*) "!==============================================================================!"
	write(*,*) "Step 0: Welcome Message!"
	write(*,*) "Welcome to GOM."
	write(*,*) "GOM is developed by Jungwoo Lee & Jun Lee."
	write(*,*) "Then, good luck!"
	write(*,*) " "
	write(*,'(A35,I4,A1,I2,A1,I2, A3, I2,A1,I2,A1,I2)') &
	&                "Today is (yyyy:mm:dd - hh:mm:ss): ", &
	&                 values(1),":",values(2),":",values(3)," - ",values(5),":",values(6),":",values(7) 
	write(*,'(A22,F10.0,A10)') "Total elapsed time = ", start,"seconds"
	write(*,*) "!==============================================================================!"
	write(*,*) "Now, GOM is scanning main input file - main.inp"
	write(*,*) "."
	write(*,*) "."
	write(*,*) "."
	!=== End of writing welcome message on the screen =========================!	

end subroutine welcome_message