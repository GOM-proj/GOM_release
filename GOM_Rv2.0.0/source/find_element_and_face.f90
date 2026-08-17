!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! GIVEN TWO ADJACENT NODES (NUMBERS), THIS ROUTINE DETERMINES
!! THE CORRESPONDING FACE (OR SIDE) NUMBER AND ELEMENT NUMBER
!!
!! Note: 
!! 		n1 and n2 should be given in counterclockwise order
 subroutine find_element_and_face(n1,n2)
	use mod_global_variables
	implicit none
	
	integer, intent(in) :: n1, n2
! jw
	integer :: i, j, j1, j2, nod1, nod2
	
	! jw
	face_id = 0
	element_id = 0
	
	! jw
outer:	do i = 1, MAXELE
	inner:	do j1 = 1, 4
					nod1 = nodenum_at_cell(j1,i)	! jw
					j2 = j1 + 1
	
					if(j2 > 4) then
						j2 = 1
					end if
			
					nod2 = nodenum_at_cell(j2,i)	! jw

					if(nod2 /= 0) then	! jw
						if(nod1 == n1 .AND. nod2 == n2)then
							element_id = i
							! jw
							! jw
							
							! jw
							! jw
							! jw
							exit outer
						end if
					else					! jw
						nod2 = nodenum_at_cell(1,i)
						if(nod1 == n1 .AND. nod2 == n2)then
							element_id = i
							! jw
							! jw
							
							! jw
							! jw
							exit outer
						end if
					end if
				end do inner ! jw
			end do outer ! jw
	
	! jw
	! jw
	do j=1,maxface
		if(n1 == nodenum_at_face(1,j) .and. n2 == nodenum_at_face(2,j)) then
			face_id = j
			exit
		end if
	end do
end subroutine find_element_and_face
