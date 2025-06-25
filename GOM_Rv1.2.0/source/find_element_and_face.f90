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
! 	integer, intent(inout) :: face_id, element_id
	integer :: i, j, j1, j2, nod1, nod2
	
	! initialize face_id and element_id before beginning
	face_id = 0
	element_id = 0
	
	! Find element_id =========================================================!
outer:	do i = 1, MAXELE
	inner:	do j1 = 1, 4
					nod1 = nodenum_at_cell(j1,i)	! first node in an element
					j2 = j1 + 1
	
					if(j2 > 4) then
						j2 = 1
					end if
			
					nod2 = nodenum_at_cell(j2,i)	! second node in an element

					if(nod2 /= 0) then	! quadrilateral grid
						if(nod1 == n1 .AND. nod2 == n2)then
							element_id = i
							! face_id = j1
							! here, j1 is the face count number (order) at this element: 1, 2, 3, 4th face in this element
							
							! this approach is not clear, so I will use next approach
							! jw, I don't know why face number and face_id is not aligned...
							! face_id = facenum_at_cell(j1-1,element_id)
							exit outer
						end if
					else					! triangular grid
						nod2 = nodenum_at_cell(1,i)
						if(nod1 == n1 .AND. nod2 == n2)then
							element_id = i
							! face_id = j1
							! here, j1 is the face count number (order) at this element: 1, 2, 3, 4th face in this element
							
							! this approach is not clear, so I will use next approach
							! face_id = facenum_at_cell(j1-1,element_id)
							exit outer
						end if
					end if
				end do inner ! j1 = 1,4
			end do outer ! i = 1,MAXELE	
	
	! Find face_id ============================================================!
	! This is the new robust approach
	do j=1,maxface
		if(n1 == nodenum_at_face(1,j) .and. n2 == nodenum_at_face(2,j)) then
			face_id = j
			exit
		end if
	end do
end subroutine find_element_and_face
