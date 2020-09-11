# python gradunwarp-1.1.0/bin/gradient_unwarp.py \
# 		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM\16_t1_mprage_sag_CORR_DISTOR_3D_pos1_ND.nii') \
# 		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM_DC\16_t1_mprage_sag_CORR_DISTOR_3D_pos1_ND_DC.nii') \
# 		siemens \
# 		-g $(wpc 'C:\Users\Lucas\Desktop\test_nii\coeff_grad\Verio_coeff_AS097.grad') \
# 		-n --numpoints 160 --interp_order 1 --fovmin -.2 --fovmax .13 \
# 		--verbose 


python bin/gradient_unwarp.py \
		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM\19_t1_mprage_sag_CORR_DISTOR_2D_pos1_ND.nii') \
		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM_DC\19_t1_mprage_sag_CORR_DISTOR_2D_pos1_ND_DC_N160_Jacob_FOV_-0122_0198.nii') \
		siemens \
		-g $(wpc 'C:\Users\Lucas\Desktop\Matlab\Projects\2019-08_ihMT_PostProcess_bash\Tools\gradunwarp\coeff_grad\Verio_coeff_AS097.grad') \
		--numpoints 160 --interp_order 1 --fovmin -.300 --fovmax .300 \
		--verbose 
rm -f fullWarp_abs.nii.gz



# python gradunwarp-1.1.0/bin/gradient_unwarp.py \
# 		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM\19_t1_mprage_sag_CORR_DISTOR_2D_pos1_ND.nii') \
# 		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM_DC\19_t1_mprage_sag_CORR_DISTOR_2D_pos1_ND_DC_N160_noJacob.nii') \
# 		siemens \
# 		-g $(wpc 'C:\Users\Lucas\Desktop\test_nii\coeff_grad\Verio_coeff_AS097.grad') \
# 		-n --numpoints 160 --interp_order 1 --fovmin -.2 --fovmax .13 \
# 		--verbose 


# python gradunwarp-1.1.0/bin/gradient_unwarp.py \
# 		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM\19_t1_mprage_sag_CORR_DISTOR_2D_pos1_ND.nii') \
# 		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM_DC\19_t1_mprage_sag_CORR_DISTOR_2D_pos1_ND_DC_N160_Jacob_FOV_-032_032.nii') \
# 		siemens \
# 		-g $(wpc 'C:\Users\Lucas\Desktop\test_nii\coeff_grad\Verio_coeff_AS097.grad') \
# 		--numpoints 160 --interp_order 1 --fovmin -.32 --fovmax .32 \
# 		--verbose 

# python gradunwarp-1.1.0/bin/gradient_unwarp.py \
# 		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM\19_t1_mprage_sag_CORR_DISTOR_2D_pos1_ND.nii') \
# 		$(wpc 'C:\Users\Lucas\Desktop\test_nii\SAM_DC\19_t1_mprage_sag_CORR_DISTOR_2D_pos1_ND_DC_N160_Jacob_FOV_-0198_0122.nii') \
# 		siemens \
# 		-g $(wpc 'C:\Users\Lucas\Desktop\test_nii\coeff_grad\Verio_coeff_AS097.grad') \
# 		--numpoints 160 --interp_order 1 --fovmin -.198 --fovmax .122 \
# 		--verbose 