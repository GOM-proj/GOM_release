!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine find_omp_range(g_start,g_end,nthreads,irank,i_start,i_end)
	implicit none
	integer,intent(in) :: g_start, g_end, nthreads,irank
	integer,intent(inout):: i_start,i_end
	integer :: share, rest
	! jw
	
	share = (g_end - g_start + 1)/nthreads
	rest = mod(g_end - g_start + 1, nthreads)
	i_start = irank * share + g_start + min(irank,rest)
	i_end = i_start + share - 1
	if(rest > irank) then
		i_end = i_end + 1
	end if
end subroutine find_omp_range