	INTEGER*8 :: H5_STRING_T
	INTEGER*8 :: H5_INT16_T
	INTEGER*8 :: H5_INT32_T
	INTEGER*8 :: H5_INT64_t
	INTEGER*8 :: H5_FLOAT32_T
	INTEGER*8 :: H5_FLOAT64_T
	INTEGER*8 :: H5_MAX_NAME_LEN
	INTEGER*8 :: H5_SUCCESS
	INTEGER*8 :: H5_OK
	INTEGER*8 :: H5_NOK
	INTEGER*8 :: H5_FAILURE
	INTEGER*8 :: H5_ERR_BADF
	INTEGER*8 :: H5_ERR_NOMEM
	INTEGER*8 :: H5_ERR_INVAL
	INTEGER*8 :: H5_ERR_BADFD
	INTEGER*8 :: H5_ERR_LAYOUT
	INTEGER*8 :: H5_ERR_NOENTRY
	INTEGER*8 :: H5_ERR_MPI
	INTEGER*8 :: H5_ERR_HDF5
	INTEGER*8 :: H5_ERR_H5
	INTEGER*8 :: H5_ERR_H5PART
	INTEGER*8 :: H5_ERR_H5BLOCK
	INTEGER*8 :: H5_ERR_H5FED
	INTEGER*8 :: H5_ERR_INTERNAL
	INTEGER*8 :: H5_ERR_NOT_IMPLEMENTED
	PARAMETER (H5_STRING_T =  1)
	PARAMETER (H5_INT16_T =   2)
	PARAMETER (H5_INT32_T =   3)
	PARAMETER (H5_INT64_T =   4)
	PARAMETER (H5_FLOAT32_T = 5)
	PARAMETER (H5_FLOAT64_T = 6)
	PARAMETER (H5_MAX_NAME_LEN = 64)
	PARAMETER (H5_SUCCESS = 		0)
	PARAMETER (H5_OK =			H5_SUCCESS)
	PARAMETER (H5_NOK =			-1)
	PARAMETER (H5_FAILURE = 		-2)
	PARAMETER (H5_ERR_BADF =		-9)
	PARAMETER (H5_ERR_NOMEM	=       	-12)
	PARAMETER (H5_ERR_INVAL	=       	-22)
	PARAMETER (H5_ERR_BADFD	=       	-77)
	PARAMETER (H5_ERR_LAYOUT =		-100)
	PARAMETER (H5_ERR_NOENTRY =		-101)
	PARAMETER (H5_ERR_MPI = 		-201)
	PARAMETER (H5_ERR_HDF5 =		-202)
	PARAMETER (H5_ERR_H5 =  		-203)
	PARAMETER (H5_ERR_H5PART =		-204)
	PARAMETER (H5_ERR_H5BLOCK =		-205)
	PARAMETER (H5_ERR_H5FED	=       	-206)
	PARAMETER (H5_ERR_INTERNAL =		-253)
	PARAMETER (H5_ERR_NOT_IMPLEMENTED =	-254)
	INTEGER*8 h5_openr
	INTEGER*8 h5_openw
	INTEGER*8 h5_opena
	INTEGER*8 h5_openr_par
	INTEGER*8 h5_openw_par
	INTEGER*8 h5_opena_par
	INTEGER*8 h5_openr_align
	INTEGER*8 h5_openw_align
	INTEGER*8 h5_opena_align
	INTEGER*8 h5_openr_par_align
	INTEGER*8 h5_openw_par_align
	INTEGER*8 h5_opena_par_align
	INTEGER*8 h5_close
	INTEGER*8 h5_finalize
	INTEGER*8 h5_check
	INTEGER*8 h5_setstep
	INTEGER*8 h5_getstep
	INTEGER*8 h5_getnsteps
	INTEGER*8 h5_set_verbosity_level
	INTEGER*8 h5_getnfileattribs
	INTEGER*8 h5_getfileattribinfo
	INTEGER*8 h5_writefileattrib_string
	INTEGER*8 h5_readfileattrib_string
	INTEGER*8 h5_writefileattrib_r8
	INTEGER*8 h5_readfileattrib_r8
	INTEGER*8 h5_writefileattrib_r4
	INTEGER*8 h5_readfileattrib_r4
	INTEGER*8 h5_writefileattrib_i8
	INTEGER*8 h5_readfileattrib_i8
	INTEGER*8 h5_writefileattrib_i4
	INTEGER*8 h5_readfileattrib_i4
	INTEGER*8 h5_getnstepattribs
	INTEGER*8 h5_getstepattribinfo
	INTEGER*8 h5_writestepattrib_string
	INTEGER*8 h5_readstepattrib_string
	INTEGER*8 h5_writestepattrib_r8
	INTEGER*8 h5_readstepattrib_r8
	INTEGER*8 h5_writestepattrib_r4
	INTEGER*8 h5_readstepattrib_r4
	INTEGER*8 h5_writestepattrib_i8
	INTEGER*8 h5_readstepattrib_i8
	INTEGER*8 h5_writestepattrib_i4
	INTEGER*8 h5_readstepattrib_i4
	INTEGER*8 h5pt_setnpoints
	INTEGER*8 h5pt_setnpoints_strided
	INTEGER*8 h5pt_getndatasets
	INTEGER*8 h5pt_getnpoints
	INTEGER*8 h5pt_getdatasetname
	INTEGER*8 h5pt_getdatasetinfo
	INTEGER*8 h5pt_setview
	INTEGER*8 h5pt_setview_indices
	INTEGER*8 h5pt_resetview
	INTEGER*8 h5pt_hasview
	INTEGER*8 h5pt_getview
	INTEGER*8 h5pt_writedata_r8
	INTEGER*8 h5pt_writedata_r4
	INTEGER*8 h5pt_writedata_i8
	INTEGER*8 h5pt_writedata_i4
	INTEGER*8 h5pt_readdata_r8
	INTEGER*8 h5pt_readdata_r4
	INTEGER*8 h5pt_readdata_i8
	INTEGER*8 h5pt_readdata_i4
	INTEGER*8 h5bl_3d_setview
	INTEGER*8 h5bl_3d_getview
	INTEGER*8 h5bl_3d_getreducedview
	INTEGER*8 h5bl_3d_hasview
	INTEGER*8 h5bl_3d_setchunk
	INTEGER*8 h5bl_getnumfields
	INTEGER*8 h5bl_getfieldinfo
	INTEGER*8 h5bl_getnfieldattribs
	INTEGER*8 h5bl_getfieldattribinfo
	INTEGER*8 h5bl_writefieldattrib_string
	INTEGER*8 h5bl_readfieldattrib_string
	INTEGER*8 h5bl_writefieldattrib_r8
	INTEGER*8 h5bl_readfieldattrib_r8
	INTEGER*8 h5bl_writefieldattrib_r4
	INTEGER*8 h5bl_readfieldattrib_r4
	INTEGER*8 h5bl_writefieldattrib_i8
	INTEGER*8 h5bl_readfieldattrib_i8
	INTEGER*8 h5bl_writefieldattrib_i4
	INTEGER*8 h5bl_readfieldattrib_i4
	INTEGER*8 h5bl_3d_get_field_spacing
	INTEGER*8 h5bl_3d_set_field_spacing
	INTEGER*8 h5bl_3d_get_field_origin
	INTEGER*8 h5bl_3d_set_field_origin
	INTEGER*8 h5bl_3d_write_scalar_field_r8
	INTEGER*8 h5bl_3d_read_scalar_field_r8
	INTEGER*8 h5bl_3d_write_vector3d_field_r8
	INTEGER*8 h5bl_3d_read_vector3d_field_r8
	INTEGER*8 h5bl_3d_write_scalar_field_r4
	INTEGER*8 h5bl_3d_read_scalar_field_r4
	INTEGER*8 h5bl_3d_write_vector3d_field_r4
	INTEGER*8 h5bl_3d_read_vector3d_field_r4
	INTEGER*8 h5bl_3d_write_scalar_field_i8
	INTEGER*8 h5bl_3d_read_scalar_field_i8
	INTEGER*8 h5bl_3d_write_vector3d_field_i8
	INTEGER*8 h5bl_3d_read_vector3d_field_i8
	INTEGER*8 h5bl_3d_write_scalar_field_i4
	INTEGER*8 h5bl_3d_read_scalar_field_i4
	INTEGER*8 h5bl_3d_write_vector3d_field_i4
	INTEGER*8 h5bl_3d_read_vector3d_field_i4
