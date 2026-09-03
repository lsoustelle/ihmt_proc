# ihMT processing
A self-contained package for ihMT image pre-/post-processing.

## Installation
The package requires a dedicated [conda/miniconda](https://anaconda.org/) environment.
```bash
git clone https://github.com/lsoustelle/ihmt_proc.git && cd ihmt_proc
conda env create -f environment.yml
conda activate ihmtproc
pip install .
```

## ihMT pre-processing
The (mostly-optional) pre-processing steps are as follows:
1. **MP-PCA denoising**: uses up-to-date [(t)MP-PCA denoiser](https://github.com/lsoustelle/tMPPCA) [1].
2. **Gibbs-ringing removal**: either (i) [`svsdegibbs`](https://github.com/lsoustelle/svsdegibbs) (port of `mrdegibbs` from MRtrix3 [dev](https://github.com/MRtrix3/mrtrix3/tree/dev) branch; for 3D volumes [2,3]), or (ii) cosine apodization+2-fold upscaling (legacy).
3. **Gradient non-linearity distortion correction**: gradient unwarping ([`gradunwarp`](https://github.com/Washington-University/gradunwarp/tree/master)); requires to provide the associated \*.grad coil file to your MR system (Siemens only).
4. **Motion correction**: ihMT-MoCo [4-6], or [ANTs'](https://github.com/antsx/antspy) antsMotionCorr onto MT<sub>0</sub> or average target [7]. ihMT-MoCo uses niprep' implementation of [`mri_synthstrip`](https://github.com/nipreps/synthstrip) and [FreeSurfer's weights](https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/) (automatically pulled) [8]. When using ihMT-MoCo, you accept the terms stated in the following [EULA](https://crmbm.univ-amu.fr/resources/ihmt-moco/).
5. **ihMT and derived maps**: outputs pre-processed ihMT, ihMTR, MTR-single, MTR-dual, normalized MT-single and/or normalized MT-dual.

### Usage
```bash
proc-ihMT	${FLD_DATA}/ihMT_2p5.nii \
			${FLD_DATA}/res_ \
			--maps ihMTp,ihMTR \
			--mppca \
			--unring 1 \
			--gnldc --gnldc_grad coeff.grad \
			--moco 1 \
			--out_int --out_gz \
			--nthreads 16
```

## ihMT post-processing
The package supports the calculation of ihMT-sat maps as performed in Ref. [9], providing a T<sub>1</sub> map (mandatory) and a B<sub>1</sub> map (highly recommended).

### Usage
```bash
fit-ihMTsat 	${FLD_DATA}/res_ihMT.nii \
				${FLD_DATA}/T1map.nii \
				${FLD_DATA}/ihMTsat.nii \
				5.0,1,10,100.0,10.0 \
 				8.26,7,96,0,2000.0 \
 				--B1 ${FLD_DATA}/B1map.nii.gz \
				--mask ${FLD_DATA}/mask.nii.gz \
				--ihMTsatB1sq ${FLD_DATA}/ihMTsatB1sq.nii.gz \
				--nthreads 16 \
				-R 1 -S 2,4 -D 3,5
```


## References:
1. Olesen et al., Tensor denoising of multidimensional MRI data, Magn Reson Med. 2022 [(doi: 10.1002/mrm.29478)](https://doi.org/10.1002/mrm.29478)
2. Kellner et al., Gibbs-ringing artifact removal based on local subvoxel-shifts. Magn Reson Med. 2015 [(doi: 10.1002/mrm.26054)](https://doi.org/10.1002/mrm.26054)
3. Bautista et al., Removal of Gibbs ringing artefacts for 3D acquisitions using subvoxel shifts. Proc. ISMRM, 2021 [(p. 3535)](https://cds.ismrm.org/protected/21MProceedings/PDFfiles/3535.html)
4. Soustelle et al., A Motion Correction Strategy for Multi-Contrast based 3D parametric imaging: Application to Inhomogeneous Magnetization Transfer (ihMT), bioRxiv, 2020 [(doi: 10.1101/2020.09.11.292649)](https://doi.org/10.1101/2020.09.11.292649)
5. Soustelle et al., A strategy to reduce the sensitivity of inhomogeneous magnetization transfer (ihMT) imaging to radiofrequency transmit field variations at 3 T, Magn Reson Med., 2022 [(doi: 10.1002/mrm.29055)](https://doi.org/10.1002/mrm.29478)
6. Anderson et al., SNR‐Efficient Inhomogeneous Magnetization Transfer (ihMT) for Clinical Applications at 7 T, Magn Reson Med. 2026 [(doi: 10.1002/mrm.70419)](https://onlinelibrary.wiley.com/doi/10.1002/mrm.70419)
7. Avants et al., A reproducible evaluation of ANTs similarity metric performance in brain image registration, NeuroImage 2011 [(doi: 10.1016/j.neuroimage.2010.09.025)](http://dx.doi.org/10.1016/j.neuroimage.2010.09.025)
8. Hoopes et al., SynthStrip: skull-stripping for any brain image, NeuroImage 2022 [(doi: 10.1016/j.neuroimage.2022.119474)](https://doi.org/10.1016/j.neuroimage.2022.119474)
9. Munsch et al., Characterization of the cortical myeloarchitecture with inhomogeneous magnetization transfer imaging (ihMT), NeuroImage 2020 [(doi: 10.1016/j.neuroimage.2020.117442)](https://doi.org/10.1016/j.neuroimage.2020.117442)
