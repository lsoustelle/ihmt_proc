% Denoise test
clear; close all;
addpath ../NIfTI_20140122_Toolbox

% vol_nii = load_nii('C:\Users\Luka\Desktop\DiffUTE_Appli\CTL06\CTL06_DiffUTE1_reg.nii');
% vol_nii = load_nii('\\ipb\dskcommun\soustelle\HeadCTL06_10814394_1_Default_DiffUTE_20170721_11142075_6.0.1.PvDatasets_FILES\Unregistered\UTE\1\1_clipped_0001.nii');
vol_nii = load_nii('\\ipb\dskcommun\soustelle\HeadCTL09_11404685_1_Default_DiffUTE_20171019_11699194_6.0.1\Unregistered\UTE\1\1_clipped_0001.nii');
std = 4000;
% std= 4650/sqrt(2-pi/2);
% std = 0;
tic
[y_est, sigma_est] = bm4d(vol_nii.img, 'Rice', std, 'mp');
toc

figure
imshow(y_est(:,:,9),[1000 5e4]);
set(gcf,'position',[1 1 1000 1000])

% javaaddpath 'C:\Program Files\MATLAB\R2016a\java\mij.jar'
% javaaddpath 'C:\Program Files\MATLAB\R2016a\java\ij.jar'
% MIJ.start
% MIJ.createImage(Img)
% MIJ.createImage(y_est)