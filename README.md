# Minimal requirement: 
* ANTs (https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS)
* Bash > 4.0 (tested on Bash 4.4 & 5.0)

# Additional features & dependencies:
* ihMT-MoCo [1] cannot be included in this repository since under license. The script is freely available, provided that the End User License Agreement is signed, at: https://crmbm.univ-amu.fr/ihmt-moco/
	* ihMT-MoCo requires the BET tool from [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)

* DICOMs handling: dcm2niix package & set in the PATH variable

* Degibbsing:
1.  unring: mrdegibbs through [MRtrix3](http://userdocs.mrtrix.org/en/latest/installation/linux_install.html) 
2.  cosine apodization (python3 with nibabel & numpy): provided *apodize_cos.py* script

* Denoising: MP-PCA: dwidenoise through [MRtrix3](http://userdocs.mrtrix.org/en/latest/installation/linux_install.html)  

* ihMTsat computation:
	- ihMTsat computation [2] is made possible (*fit_ihMTsat_CLI.py*), providing pre-processed ihMT and co-registered T1 map as minimal image inputs. 
	- See -h option for further description about inputs/outputs.
	- The script requires up-to-date Python (developed on 3.8.2) with numpy (≥ 1.19.4), scipy (≥ 1.5.4) and nibabel (≥ 3.2.0) packages.

* References:
[1] L. Soustelle et al., A Motion Correction Strategy for Multi-Contrast based 3D parametric imaging : Application to Inhomogeneous Magnetization Transfer (ihMT), bioRxiv 2020

[2] F. Munsch et al., Characterization of the cortical myeloarchitecture with inhomogeneous magnetization transfer imaging (ihMT), NeuroImage 2020;225:117442
 
[3] L. Soustelle et al., A strategy to reduce the sensitivity of inhomogeneous magnetization transfer (ihMT) imaging to radiofrequency transmit field variations at 3 T, Magnetic Resonance in Medicine, 2022 (DOI: 10.1002/mrm.29055)