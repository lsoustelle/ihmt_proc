# Minimal requirement: 
* ANTs (https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS)
* Bash > 4.0 (tested on Bash 4.4 & 5.0)

# Additional features & dependencies:
* ihMT-MoCo cannot be included in this repository since under license. The script is freely available, provided that the End User License Agreement is signed, at: https://crmbm.univ-amu.fr/ihmt-moco/
	* ihMT-MoCo requires the BET tool from FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)

* DICOMs handling: dcm2niix package (https://github.com/rordenlab/dcm2niix) & set in the PATH variable

* Degibbsing:
	1.  unring (https://bitbucket.org/reisert/unring): provided binary compatible on Debian 9 / Ubuntu 18.04 --- to be set in your path in any case
	2.  cosine apodization (python3 with nibabel & numpy): provided .py script

* Denoising:
	1.  MP-PCA (http://userdocs.mrtrix.org/en/latest/installation/linux_install.html): dwidenoise by MRtrix3 
	2.  BM4D (Matlab with Parallel Computing Toolbox)


