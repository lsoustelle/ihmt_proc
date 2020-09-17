# ///////////////////////////////////////////////////////////////////////////////////////////
# //                 Université d’Aix Marseille (AMU) -
# //                 Centre National de la Recherche Scientifique (CNRS) –
# //                 Copyright © 2018 AMU, CNRS All Rights Reserved.
# //
# //     These computer program listings and specifications, herein, are
# //     the property of Université d’Aix Marseille and CNRS
# //     shall not be reproduced or copied or used in whole or in part as
# //     the basis for manufacture or sale of items without written permission.
# //     For a license agreement, please contact: <mailto: licensing@sattse.com>
# //
# //	 Developed by Lucas Soustelle, PhD, Aix Marseille Univ, CNRS, CRMBM
# // 	 Mail: lucas.soustelle@univ-amu.fr
# //////////////////////////////////////////////////////////////////////////////////////////

import sys
import numpy as np
import numpy.matlib
import nibabel as nib

print('---------------------')
print('Proceeding to apodization')
print('---------------------')

nii_path  = sys.argv[1]
interpfac = 2

## load nifti
img = nib.load(nii_path)

## build cosine kernel                          
wx  = np.array( np.cos( ( np.arange(img.shape[0]) - np.floor(img.shape[0]/2) ) * np.pi / img.shape[0] ) )
wy  = np.array( np.cos( ( np.arange(img.shape[1]) - np.floor(img.shape[1]/2) ) * np.pi / img.shape[1] ) )
wz  = np.array( np.cos( ( np.arange(img.shape[2]) - np.floor(img.shape[2]/2) ) * np.pi / img.shape[2] ) )
# wx  = np.array( np.cos( ( np.arange(img.shape[0]//2+1) - np.floor(img.shape[0]//4) ) * np.pi / (img.shape[0]/2) ) )
# wx 	= np.concatenate( (np.zeros(img.shape[0]//4), wx, np.zeros(img.shape[0]//4)) )
# wy  = np.array( np.cos( ( np.arange(img.shape[1]//2) - np.floor(img.shape[1]//4) ) * np.pi / (img.shape[1]/2) ) )
# wy 	= np.concatenate( (np.zeros(img.shape[1]//4), wy, np.zeros(img.shape[1]//4)) )
# wz  = np.array( np.cos( ( np.arange(img.shape[2]//2) - np.floor(img.shape[2]//4) ) * np.pi / (img.shape[2]/2) ) )
# wz 	= np.concatenate( (np.zeros(img.shape[2]//4), wz, np.zeros(img.shape[2]//4)) )

wx  = wx[np.newaxis,:]
wy  = wy[np.newaxis,:]
wz  = wz[np.newaxis,np.newaxis,:]
w2D = np.dot(wx.T,wy)
w2D = w2D[:,:,np.newaxis]  
w3D = np.squeeze(np.dot(w2D, wz))

## process fft & apodize
img_apo  = np.zeros(( img.shape[0]*interpfac,img.shape[1]*interpfac,img.shape[2]*interpfac,img.shape[3] ))
img_data = img.dataobj
for ii in range(img.shape[3]):
    img_data_tmp    = img_data[:,:,:,ii]
    img_data_fft    = np.fft.fftshift( np.fft.fftn(img_data_tmp) ) / ( img.shape[0] * img.shape[1] * img.shape[2] )
    img_data_ifft   = np.fft.ifftn( img_data_fft * w3D, \
                                   [ img.shape[0]*interpfac,img.shape[1]*interpfac,img.shape[2]*interpfac ] ) * \
                                     img.shape[0]*interpfac*img.shape[1]*interpfac*img.shape[2]*interpfac
    img_apo[:,:,:,ii] = abs(img_data_ifft)

## modify pixdim & save nifti
img.header['pixdim'][1:4] = img.header['pixdim'][1:4] / interpfac
new_img = nib.Nifti1Image(img_apo, img.affine, img.header)
nib.save(new_img, nii_path)

print('Proceeding to apodization: done')
print('---------------------')
