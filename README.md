# Minimal requirement: 
* ANTs (https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS)
* Bash > 4.0 (tested on Bash 4.4 & 5.0)

# Additional features & dependencies:
* DICOMs handling: dcm2niix package & set in the PATH variable

* Degibbsing:
1.  unring (https://bitbucket.org/reisert/unring): provided binary compatible on Debian 9 / Ubuntu 18.04 --- to be set in your path in any case
2.  cosine apodization (python3 with nibabel & numpy): provided py script

* Denoising:
1.  MP-PCA (http://userdocs.mrtrix.org/en/latest/installation/linux_install.html): dwidenoise through MRtrix3 
2.  BM4D (Matlab with Parallel Computing Toolbox)

