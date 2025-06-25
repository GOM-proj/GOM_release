!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================!
!!
!! This subroutine creates grid information in netcdf filr format
!!
subroutine write_netcdf_grid
	use mod_global_variables
	use mod_file_definition
	use mod_netcdf_utils
	use netcdf
	implicit none

	! NetCDF IDs
	integer :: ncid0 ! file ID; this is for distinguishing from ncid
! 	integer :: node_dimid, elem_dimid, layer_dimid ! dimension IDs ! move to mod_netcdf_utils.f90
! 	integer :: four_dimid ! move to mod_netcdf_utils.f90
	integer :: nodenum_at_cell_dimid(2)
	integer :: varid_xn, varid_yn, varid_hn, varid_tri_or_quad, varid_nodenum_at_cell ! variable IDs
	integer :: netcdf_status
	! End of local variables ==================================================!
	
	! Create NetCDF file (overwrite if exists) ================================!
	netcdf_status = nf90_create(id_netcdf_grid, nf90_clobber, ncid0) ! this will create: ./output/netcdf_grid.nc
	call netcdf_check(netcdf_status)

	! Define dimensions =======================================================!
	call netcdf_check(nf90_def_dim(ncid0, "maxnod", maxnod, node_dimid))
	call netcdf_check(nf90_def_dim(ncid0, "maxele", maxele, cell_dimid))
! 	call netcdf_check(nf90_def_dim(ncid0, "maxface", maxface, face_dimid))
	call netcdf_check(nf90_def_dim(ncid0, "maxlayer", maxlayer, layer_dimid))
	call netcdf_check(nf90_def_dim(ncid0, "four", 4, four_dimid))

	! Define variabls =========================================================!
	! Define node_info
	call netcdf_check(nf90_def_var(ncid0, "x_node", nf90_double, (/node_dimid/), varid_xn))
	call netcdf_check(nf90_put_att(ncid0, varid_xn, "units", "[m]"))
	call netcdf_check(nf90_put_att(ncid0, varid_xn, "long_name", "node x_coordinate"))
	
	call netcdf_check(nf90_def_var(ncid0, "y_node", nf90_double, (/node_dimid/), varid_yn))
	call netcdf_check(nf90_put_att(ncid0, varid_yn, "units", "[m]"))
	call netcdf_check(nf90_put_att(ncid0, varid_yn, "long_name", "node y_coordinate"))

	call netcdf_check(nf90_def_var(ncid0, "h_node", nf90_double, (/node_dimid/), varid_hn))
	call netcdf_check(nf90_put_att(ncid0, varid_hn, "units", "[m]"))
	call netcdf_check(nf90_put_att(ncid0, varid_hn, "positive", "down from MSL"))
	call netcdf_check(nf90_put_att(ncid0, varid_hn, "long_name", "node water_depth"))
	
	! Define cell_info
	call netcdf_check(nf90_def_var(ncid0, "tri_or_quad", nf90_int,(/cell_dimid/), varid_tri_or_quad))
	call netcdf_check(nf90_put_att(ncid0, varid_tri_or_quad, "3", "triangle"))
	call netcdf_check(nf90_put_att(ncid0, varid_tri_or_quad, "4", "rectangle"))
	call netcdf_check(nf90_put_att(ncid0, varid_tri_or_quad, "long_name", "element shape: 3 = trianlge, 4 = rectangle"))
	
	! here, both methods should work:
 	! call netcdf_check(nf90_def_var(ncid0, "nodenum_at_cell", nf90_int, (/four_dimid,cell_dimid/), varid_nodenum_at_cell))
	nodenum_at_cell_dimid = (/four_dimid, cell_dimid/)
 	call netcdf_check(nf90_def_var(ncid0, "nodenum_at_cell", nf90_int, nodenum_at_cell_dimid, varid_nodenum_at_cell))
 	call netcdf_check(nf90_put_att(ncid0, varid_nodenum_at_cell, "Fill_Value", "0"))
 	call netcdf_check(nf90_put_att(ncid0, varid_nodenum_at_cell, "long_name", "cell connectivity: e.g., (nd_1, nd_2, nd_3, nd_4)"))

	! End define mode
	call netcdf_check(nf90_enddef(ncid0))

	! Write grid information: =================================================!
	call netcdf_check(nf90_put_var(ncid0, varid_xn, x_node))
	call netcdf_check(nf90_put_var(ncid0, varid_yn, y_node))
	call netcdf_check(nf90_put_var(ncid0, varid_hn, h_node))
	
	! write cell information:
	call netcdf_check(nf90_put_var(ncid0, varid_tri_or_quad, tri_or_quad))
 	call netcdf_check(nf90_put_var(ncid0, varid_nodenum_at_cell, nodenum_at_cell))

	! Close the NetCDF file ===================================================!
	call netcdf_check(nf90_close(ncid0))

end subroutine write_netcdf_grid