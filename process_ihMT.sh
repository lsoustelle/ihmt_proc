#!/usr/bin/env bash
export LC_CTYPE=C
trap control_c SIGINT # delete temp folder upon CTRL+C hit
BASEDIR=$(dirname "$0")
TOOLSDIR=${BASEDIR}/Tools

function cleanup
{ 
	rm -rf ${TMP_FLD} 
}
function control_c
{
  if [[ ${FLAG_KEEP_TMP} -eq 0 ]]; then cleanup; fi
  exit
}

function Usage {
    cat <<USAGE
Usage:
`basename $0` -i input4D -o output_prefix -c ihMTR

Example Case:
`basename $0` -i path/to/input/in_ihMT.nii -o path/to/output/out_ \\
		-c ihMTR,MTRd \\
		-n 8 -d 1 \\
		-R 1 -S 2,4 -D 3,5

Compulsory arguments:
	-i: input 4D ihMT NIfTI image OR DICOM folder
	-o: output path-prefix to 3D ihMT-derived NIfTI images
	-c: comma-separated ihMT-derived maps outputs (ihMT,ihMTR,MTRs,MTRd,ihMTRinv,MTRsinv,MTRdinv)
		ihMT:		raw data in NIfTI format
		MTRs: 		1 - MTs/MT0
		MTRd: 		1 - MTd/MT0
		ihMTR: 		2 * (MTRd - MTRs)
		MTRsinv: 	MT0 / MTs - 1
		MTRdinv: 	MT0 / MTd - 1
		ihMTRinv: 	2 * (MTRdinv - MTRsinv)

Optional arguments:
	-n: number of threads to be used for subsequent processes (default = 1)
	-d: MP-PCA denoising of raw ihMT images ((0)/1)
	-e: MP-PCA kernel extent (default = 3,3,3) 
	-w: gradient non-linearity distortion correction ((0)/1)
	-g: *.grad file name (Avanto_coeff_AS05.grad/Verio_coeff_AS097.grad so far) for distortion correction
	-N: number of points constituting the grid for distortion dorrection (default = 60)
	-t: BM4D denoising of raw MT0/MTs/MTd images ((0)/1)
	-b: BM4D denoising of ihMT-derived maps ((0)/1)
	-p: full path to Matlab binary for BM4D denoising
	-u: unring native ihMT images (default: none; 1: unring routine; 2: cos-kernel apodization & zero-filling x2)
	-m: ihMT Motion Correction (default: none; 1: ihMT-MoCo; 2: antsMotionCorr (MT0); 3: antsMotionCorr (average))
	-R: comma-separated indices of MT reference volume(s) in the 4D input stack (default is 1)
	-S: comma-separated indices of MT single volumes in the 4D input stack (default is 2,4,6,...,N-1)
	-D: comma-separated indices of MT dual volumes in the 4D input stack (default is 3,5,7,...,N)
	-k: keep temporary files ((0)/1)

Notes: 
	- ANTSPATH variable has to be properly set (minimal settings)
	- If a DICOM folder is provided, dcm2niix package has to be installed & in your PATH variable 
	- Performs: 
		1) MP-PCA denoising (if set; EXTENT should be >= amount of raw 3D volumes)
		2) BM4D denoising of raw MT0/MTw (if set)
		3) Gradient non-linearity distortion correction (if set)
		4) Unringing (if set)
		5) Motion Correction (if set)
		6) Maps computation
		7) BM4D "denoising" of computed maps (if set)
	- MP-PCA denoising (Veerart 2016) necessitates an up-to-date MRtrix3 package & in your PATH variable
	- gradient non-linearity distortion correction necessitates the right *.grad file and a large enough grid
	- unring (Kellner 2016) necessitates the unring binary in your PATH variable
	- cos-kernel apodization necessitates python3 with up-to-date numpy & nibabel libs
	- BM4D denoising (Maggioni 2013) necessitates a provided path to Matlab binary (WITH set up Parallel Computing Toolbox)
	- BM4D on computed maps is discouraged (the functionnality necessitating raw images), and should be employed as a last resort on raw maps
	- please cite the subsequent sources depending on the used functionnalities

USAGE
    exit 1
}

function func_denoiseBM4D {
	local -n _TMPNIIPATH=$1
	local -n _FLAG_DENOISE=$2

	if [[ ${_FLAG_DENOISE} -eq 1 ]]; then
		echo -e "\nBM4D denoising"
		${MATLAB_BIN} -nodisplay <<EOF
		addpath(genpath('Tools'));
		cell_in_niipaths 	= 	strsplit('${_TMPNIIPATH[@]}');
		if numel(cell_in_niipaths) > 1
			parpool(${N_THREADS});
			parfor ii = 1 : numel(cell_in_niipaths)
				tmp_nii    		= load_untouch_nii(cell_in_niipaths{ii});
			    tmp_nii.img 	= bm4d(tmp_nii.img,'Rice');
				save_untouch_nii(tmp_nii,cell_in_niipaths{ii})
			end
		else
			tmp_nii    		= load_untouch_nii(cell_in_niipaths{1});
			tmp_nii.img 	= bm4d(tmp_nii.img,'Rice');
			save_untouch_nii(tmp_nii,cell_in_niipaths{1})
		end
		exit
EOF
		echo "BM4D denoising: done"
	fi 
}

function func_denoiseMPPCA {
	local -n _TMPNIIPATH=$1
	local -n _FLAG_DENOISE=$2

	if [[ ${_FLAG_DENOISE} -eq 1 ]]; then
		echo "MP-PCA denoising"
		dwidenoise ${_TMPNIIPATH} ${_TMPNIIPATH/.nii/_MPPCA.nii} -nthreads ${N_THREADS} -force -extent ${EXTENT}
		echo "MP-PCA denoising: done"
		_TMPNIIPATH=${_TMPNIIPATH/.nii/_MPPCA.nii}
	fi
}



## Reading command line arguments
FLAG_DCM=0
FLAG_MAP_ihMT=0
FLAG_MAP_ihMTR=0
FLAG_MAP_MTRs=0
FLAG_MAP_MTRd=0
FLAG_MAP_ihMTRinv=0
FLAG_MAP_MTRsinv=0
FLAG_MAP_MTRdinv=0
FLAG_IDX=0
FLAG_DENOISE_RAW=0
FLAG_DENOISE_MT=0
FLAG_GNLDC=0
FLAG_UNRING=0
FLAG_MoCo=0
FLAG_DENOISE_MAPS=0
N_THREADS=1
if [[ $# -eq 0 ]] ; then Usage && echo $USAGE ; exit 0 ; fi
while getopts "i:o:c:n:d:e:w:g:N:t:b:p:u:m:R:S:D:k:" OPT; do
	case $OPT in
		i)	
	INPUT_PATH=$OPTARG
	;;
		o)	
	OUTPUT_NIIPATH=$OPTARG
	;;
		c)	
	OUT_MAPS=( ${OPTARG//,/' '} )
	;;
		n)
	N_THREADS=$OPTARG
	;;
		d)
	FLAG_DENOISE_RAW=$OPTARG
	;;
		e)
	EXTENT=$OPTARG
	;;
		w)
	FLAG_GNLDC=$OPTARG
	;;
		g)
	GRAD_PATH=$OPTARG
	GRAD_PATH=${TOOLSDIR}/gradunwarp/coeff_grad/${GRAD_PATH}
	;;
		N)
	N_GRID=$OPTARG
	;;
		t)
	FLAG_DENOISE_MT=$OPTARG
	;;
		b)
	FLAG_DENOISE_MAPS=$OPTARG
	;;
		p)
	MATLAB_BIN=$OPTARG
	;;
		u)
	FLAG_UNRING=$OPTARG
	;;
		m)
	FLAG_MoCo=$OPTARG
	;;
		R)
	IDX_MT0=( ${OPTARG//,/' '} )
	IDX_MT0_MoCo=$OPTARG
	;;
		S)
	IDX_MTs=( ${OPTARG//,/' '} )
	IDX_MTs_MoCo=$OPTARG
	;;
		D)
	IDX_MTd=( ${OPTARG//,/' '} )
	IDX_MTd_MoCo=$OPTARG
	;;
		k)
	FLAG_KEEP_TMP=$OPTARG
	;;
		\?)
   	Usage && echo $USAGE
  	exit 1
	;;
	esac
done


## Check inconsistencies & set default stuff
if [[ -z ${INPUT_PATH} ]]; then
	echo "----------------"
    echo "Input not defined."
    exit
fi
CHECK_INPUT=$(file -b ${INPUT_PATH})
if [[ ${CHECK_INPUT} == *"cannot open"* ]]; then
	echo "----------------"
    echo "Input not found."
    exit
elif [[ ${CHECK_INPUT} == *"directory"* ]]; then
	FLAG_DCM=1
	echo "----------------"
    echo "Input is a DICOM folder."
elif [[ ${CHECK_INPUT} == *"data"* ]]; then
	FLAG_DCM=0
	echo "----------------"
    echo "Input is a NIfTI."
fi

if [[ -z ${OUTPUT_NIIPATH} ]]; then
	echo "----------------"
    echo "Output prefix not defined."
    exit
fi

if [[ -z ${OUT_MAPS} ]]; then
	echo "----------------"
    echo "No output map defined."
    exit
fi

if [[ ! -z ${OUT_MAPS} ]]; then
	for STR_TMP in ${OUT_MAPS[@]}; do
		case $STR_TMP in
			'ihMT')
		FLAG_MAP_ihMT=1
		;;
			'ihMTR')
		FLAG_MAP_ihMTR=1
		;;
			'MTRs')
		FLAG_MAP_MTRs=1
		;;
			'MTRd')
		FLAG_MAP_MTRd=1
		;;
			'ihMTRinv')
		FLAG_MAP_ihMTRinv=1
		;;
			'MTRsinv')
		FLAG_MAP_MTRsinv=1
		;;
			'MTRdinv')
		FLAG_MAP_MTRdinv=1
		;;
		esac
	done

	CHECKSUM=$(( $FLAG_MAP_ihMT+$FLAG_MAP_ihMTR+$FLAG_MAP_MTRs+$FLAG_MAP_MTRd+$FLAG_MAP_ihMTRinv+$FLAG_MAP_MTRsinv+$FLAG_MAP_MTRdinv ))
	if [[ ${CHECKSUM} -lt ${#OUT_MAPS[@]} ]]; then
		echo "----------------"
	    echo "Inconsistencies detected in the -c option; please redefined."
	    exit
	fi
fi

if [[ -z ${ANTSPATH} ]]; then
	echo "----------------"
    echo "ANTSPATH variable not defined."
    exit
fi


if [[ ${N_THREADS} =~ ^-?[0-9]+$  ]]; then 
	echo Processing with ${N_THREADS} threads
else
	N_THREADS=1
	echo "----------------"
	echo "-n input should be an integer; setting number of threads to 1"
fi

if [[ -z `which dwidenoise` && ${FLAG_DENOISE_RAW} -eq 1 ]]; then
	echo "----------------"
	echo "dwidenoise program from MRtrix3 can't be found. MP-PCA denoising impossible."
	exit
fi

if [[ -z ${EXTENT} ]]; then
	EXTENT=3,3,3
fi

if [[ -z `which python` && ${FLAG_GNLDC} -eq 1 ]]; then
	echo "----------------"
    echo "python program can't be found. Please install the package & dependencies (numpy, scipy & nibabel)."
    exit
fi

if [[ ${FLAG_GNLDC} -eq 1 ]] && [[ ! -e ${GRAD_PATH} || -z ${GRAD_PATH} ]]; then
	echo "----------------"
    echo "Gradient file for distortion correction can't be found. Please check."
    exit
fi

if [[ ${FLAG_GNLDC} -eq 1 && -z ${N_GRID} ]]; then
	echo "----------------"
    echo "Grid size for distortion correction not provided: setting it to 60^3."
    N_GRID=60
fi

if [[ ${FLAG_MoCo} -lt 0 || ${FLAG_MoCo} -gt 3 ]]; then
	echo "----------------"
    echo "MoCo option not recognized: MoCo impossible."
    exit
fi



if [[ ${FLAG_DENOISE_MAPS} -eq 1 || ${FLAG_DENOISE_MT} -eq 1 ]] && [[ ! -e ${MATLAB_BIN} || -z ${MATLAB_BIN} ]]; then
	echo "----------------"
	echo "Matlab binary path not provided or nonexistent: BM4D denoising impossible."
	exit
fi

if [[ ${FLAG_UNRING} -eq 1 && -z `which unring` ]]; then
	echo "----------------"
    echo "unring program can't be found. Please install package & set your environment path."
    exit
fi
if [[ -z `which python3` && ${FLAG_UNRING} -eq 2 ]]; then
	echo "----------------"
    echo "python3 program can't be found. Please install the package & dependencies (numpy & nibabel)."
    exit
elif [[ ! -z `which python3` && ${FLAG_UNRING} -eq 2 ]]; then
	CHECK_LIB=$(pip3 freeze | grep "numpy\|nibabel")
	if [[ ! ${CHECK_LIB} == *"numpy"* || ! ${CHECK_LIB} == *"nibabel"* ]]; then
		echo "----------------"
	    echo "numpy or nibabel libraries not found. Make sure that both libs are installed (run 'pip3 freeze')."
	    exit
	fi
fi

[ ! -z "$IDX_MT0" ] && BOOL_MT0=1 || BOOL_MT0=0
[ ! -z "$IDX_MTs" ] && BOOL_MTs=1 || BOOL_MTs=0
[ ! -z "$IDX_MTd" ] && BOOL_MTd=1 || BOOL_MTd=0
CHECKSUM=$(( $BOOL_MT0+$BOOL_MTs+$BOOL_MTd ))
if [[ ${CHECKSUM} -eq 3 ]]; then
	FLAG_IDX=1
	IDX_CAT+=( "${IDX_MT0[@]}" "${IDX_MTs[@]}" "${IDX_MTd[@]}" )
	DUPLICATE=$(printf '%s\n' "${IDX_CAT[@]}"|awk '!($0 in seen){seen[$0];next} 1')
	if [[ ! -z $DUPLICATE ]]; then
		echo "----------------"
		echo "Inconsistencies in provided -R/-S/-D options: at least one index is duplicated."
		exit
	fi
elif [[ ${CHECKSUM} -eq 1 || ${CHECKSUM} -eq 2 ]]; then
	echo "----------------"
	echo "Inconsistencies in provided -R/-S/-D options: at least one option does not exist."
	exit
fi

if [ -z ${FLAG_KEEP_TMP} ] 	; then FLAG_KEEP_TMP=0 ; fi


## Declare tmp stuff & paths
RDM_STR='tmp'
if [[ ${OUTPUT_NIIPATH} == *'/'* ]]; then
	TMP_FLD=$(dirname ${OUTPUT_NIIPATH}) 
else
	TMP_FLD=$(readlink -f ${OUTPUT_NIIPATH}) # root is basically current directory
	TMP_FLD=$(dirname ${TMP_FLD})
fi
TMP_FLD=${TMP_FLD}/process_ihMT_${RDM_STR}/
mkdir -p ${TMP_FLD}
if [[ ! -d ${TMP_FLD} ]]; then
	echo "----------------"
    echo "Can't write in provided output folder. Please provide a proper one."
    exit
fi

# takes care of dcm2niix if necessary
if [[ -z `which dcm2niix` && FLAG_DCM -eq 1 ]]; then
	echo "----------------"
    echo "dcm2niix program not found. Please install package & set your environment path."
    exit
elif [[ ! -z `which dcm2niix` && FLAG_DCM -eq 1 ]]; then
	echo "----------------"
    echo "Proceeding to dcm2niix conversion"
    dcm2niix -f %s_%d -o ${TMP_FLD} ${INPUT_PATH}
    readarray -t ARRAY_CONV_NII < <(find ${TMP_FLD} -name \*.nii)
    if [[ ${#ARRAY_CONV_NII[@]} -ge 2 ]]; then
    	echo "----------------"
	    echo "Multiple DICOMs were contained in the provided input path. Please provide only one DICOM dataset."
	    cleanup
	    exit
	else
		INPUT_PATH=${ARRAY_CONV_NII[@]}
    fi
fi

ihMT_TMPNIIPATH=${TMP_FLD}/$(basename ${INPUT_PATH})
cp -f ${INPUT_PATH} ${ihMT_TMPNIIPATH} 2>/dev/null
if [[ ${ihMT_TMPNIIPATH} == *'.gz'* ]]; then  # keep it nii
	gzip -d ${ihMT_TMPNIIPATH}
	ihMT_TMPNIIPATH=${ihMT_TMPNIIPATH/.gz/}
fi

# get nii header info
SIZESTRING=$( ${ANTSPATH}/PrintHeader ${ihMT_TMPNIIPATH} 2 ) # header information to catch 4th dimension amount
SIZESTRING="${SIZESTRING%\\n}"
SIZE=( `echo ${SIZESTRING##*x}` )
IDX_END=$((1000+$(($SIZE-1))))
IDX_SUFF=( `seq 1000 ${IDX_END}` ) # as is with ANTs split

if [[ FLAG_IDX -eq 0 ]]; then # simplest case [MT0, MTs, MTd, ..., MTd]
	IDX_MT0=( 1 ) # first volume indexed as "1"
	IDX_MTs=( $(seq 2 2 $(($SIZE))) )
	IDX_MTd=( $(seq 3 2 $(($SIZE))) )
	IDX_MT0_MoCo=${IDX_MT0}
	IDX_MTs_MoCo=${IDX_MTs}
	IDX_MTd_MoCo=${IDX_MTd}
fi

INPUT_PREF=${ihMT_TMPNIIPATH/%.*} # delete extension
INPUT_PREF=${INPUT_PREF##*/}
SPLIT_PREF=${TMP_FLD}${INPUT_PREF}_splitted_.nii
SPLIT_NIIS=(		${IDX_SUFF[@]/%/.nii} 									)
SPLIT_NIIPATHS=( 	${SPLIT_NIIS[@]/#/${TMP_FLD}${INPUT_PREF}_splitted_} 	)	

MT0_AVER_TMPNIIPATH=${TMP_FLD}/MT0_aver.nii
MTs_AVER_TMPNIIPATH=${TMP_FLD}/MTs_aver.nii
MTd_AVER_TMPNIIPATH=${TMP_FLD}/MTd_aver.nii
for ii in ${IDX_MT0[@]}; do
	MT0_INDIV_TMPNIIPATH[$ii]=${SPLIT_NIIPATHS[$ii-1]}
done
for ii in ${IDX_MTs[@]}; do
	MTs_INDIV_TMPNIIPATH[$ii]=${SPLIT_NIIPATHS[$ii-1]}
done
for ii in ${IDX_MTd[@]}; do
	MTd_INDIV_TMPNIIPATH[$ii]=${SPLIT_NIIPATHS[$ii-1]}
done
MASK_TMPNIIPATH=${TMP_FLD}/mask_0.nii
ONES_TMPNIIPATH=${TMP_FLD}/ones.nii
MTRs_TMPNIIPATH=${TMP_FLD}/MTRs.nii
MTRd_TMPNIIPATH=${TMP_FLD}/MTRd.nii
ihMTR_TMPNIIPATH=${TMP_FLD}/ihMTR.nii
ihMTRinv_TMPNIIPATH=${TMP_FLD}/ihMTRinv.nii
MTRsinv_TMPNIIPATH=${TMP_FLD}/MTRsinv.nii
MTRdinv_TMPNIIPATH=${TMP_FLD}/MTRdinv.nii

## Denoise native ihMT images (MP-PCA)
func_denoiseMPPCA ihMT_TMPNIIPATH FLAG_DENOISE_RAW

# Denoise MT0/MTw (BM4D)
func_denoiseBM4D ihMT_TMPNIIPATH FLAG_DENOISE_MT

## Unring
if [[ ${FLAG_UNRING} -eq 1 && ! -z `which unring` ]]; then 
	echo 'Unringing (unring)'
	ihMT_TMPNIIGZPATH=${ihMT_TMPNIIPATH/.nii/.nii.gz}
	${BASEDIR}/unring ${ihMT_TMPNIIPATH} 						${ihMT_TMPNIIGZPATH/.nii.gz/_d1.nii.gz} -d 1
	${BASEDIR}/unring ${ihMT_TMPNIIGZPATH/.nii.gz/_d1.nii.gz} 	${ihMT_TMPNIIGZPATH/.nii.gz/_d1d3.nii.gz} -d 3
	mv ${ihMT_TMPNIIGZPATH/.nii.gz/_d1d3.nii.gz} ${ihMT_TMPNIIGZPATH}
	gzip -d -f ${ihMT_TMPNIIGZPATH}
	echo 'Unringing (unring): done'
elif [[ ${FLAG_UNRING} -eq 2 ]]; then
	echo 'Unringing (cos-kernel apodization)'
	python3 ${BASEDIR}/apodize_cos.py ${ihMT_TMPNIIPATH}
	${ANTSPATH}/ImageMath 4 ${ihMT_TMPNIIPATH} + ${ihMT_TMPNIIPATH} 0  # refresh header
	echo 'Unringing (cos-kernel apodization): done'
fi

## Gradient non-linearity distortion correction
if [[ ${FLAG_GNLDC} -eq 1 ]]; then
	echo 'Gradient non-linearity distortion correction'
	${ANTSPATH}/ImageMath 4 ${SPLIT_PREF} TimeSeriesDisassemble ${ihMT_TMPNIIPATH}
	ihMT_INDIV_TMPNIIPATH+=( ${MT0_INDIV_TMPNIIPATH[@]} ${MTs_INDIV_TMPNIIPATH[@]} ${MTd_INDIV_TMPNIIPATH[@]} ) # correct only volumes of interest

	for TMP_NIIPATH in ${ihMT_INDIV_TMPNIIPATH[@]}; do
		echo ${TMP_NIIPATH}
		python ${TOOLSDIR}/gradunwarp/bin/gradient_unwarp.py \
				${TMP_NIIPATH} \
				${TMP_NIIPATH} \
				siemens -g ${GRAD_PATH} \
				--numpoints ${N_GRID} --interp_order 1 \
				--fovmin -.300 --fovmax .300 \
				--verbose
	done

	${ANTSPATH}/ImageMath 4 ${ihMT_TMPNIIPATH} TimeSeriesAssemble 1 0 ${SPLIT_NIIPATHS[@]}
	rm fullWarp_abs.nii.gz
	rm ${ihMT_INDIV_TMPNIIPATH[@]}

	echo 'Gradient non-linearity distortion correction: done'
fi

## MoCo
if [[ ${FLAG_MoCo} -eq 1 ]]; then
	echo 'Motion Correction (MC-MoCo)'
	${BASEDIR}/ihMT_MoCo.sh 	-i ${ihMT_TMPNIIPATH} 	-o ${ihMT_TMPNIIPATH} 	-n ${N_THREADS} -f 0.6 \
								-R ${IDX_MT0_MoCo} 		-S ${IDX_MTs_MoCo} 		-D ${IDX_MTd_MoCo}
	echo 'Motion Correction (MC-MoCo): done'
elif [[ ${FLAG_MoCo} -eq 2 ]]; then
	echo 'Motion Correction (antsMotionCorr onto MT0)'
	TARGET_TMPNIIPATH=${TMP_FLD}/MT0_extracted100.nii
	${ANTSPATH}/ImageMath 4 ${TARGET_TMPNIIPATH/100.nii/.nii} TimeSeriesSubset ${ihMT_TMPNIIPATH} ${IDX_MT0[0]}
	${ANTSPATH}/antsMotionCorr 	-d 3 \
								-o [${ihMT_TMPNIIPATH/.nii/},${ihMT_TMPNIIPATH}] \
								-m MI[${TARGET_TMPNIIPATH},${ihMT_TMPNIIPATH},1,32,Regular,0.25] \
								-t Rigid[0.005] \
								-i 50x50 -u 1 -s 1x0 -f 2x1 -e 1 -v 1
	echo 'Motion Correction (antsMotionCorr onto MT0): done'
elif [[ ${FLAG_MoCo} -eq 3 ]]; then
	echo 'Motion Correction (antsMotionCorr onto stack averaged)'
	TARGET_TMPNIIPATH=${TMP_FLD}/ihMT_avg.nii
	${ANTSPATH}/antsMotionCorr  -d 3 \
								-a ${ihMT_TMPNIIPATH} \
								-o ${TARGET_TMPNIIPATH}

	${ANTSPATH}/antsMotionCorr 	-d 3 \
								-o [${ihMT_TMPNIIPATH/.nii/},${ihMT_TMPNIIPATH}] \
								-m MI[${TARGET_TMPNIIPATH},${ihMT_TMPNIIPATH},1,32,Regular,0.25] \
								-t Rigid[0.005] \
								-i 50x50 -u 1 -s 1x0 -f 2x1 -e 1 -v 1
	echo 'Motion Correction (antsMotionCorr onto stack averaged): done'
fi

## Compute maps
${ANTSPATH}/ImageMath 		4 ${SPLIT_PREF} TimeSeriesDisassemble ${ihMT_TMPNIIPATH}
${ANTSPATH}/AverageImages 	3 ${MT0_AVER_TMPNIIPATH} 0 ${MT0_INDIV_TMPNIIPATH[@]}
${ANTSPATH}/AverageImages 	3 ${MTs_AVER_TMPNIIPATH} 0 ${MTs_INDIV_TMPNIIPATH[@]}
${ANTSPATH}/AverageImages 	3 ${MTd_AVER_TMPNIIPATH} 0 ${MTd_INDIV_TMPNIIPATH[@]}
${ANTSPATH}/CreateImage 	3 ${MT0_AVER_TMPNIIPATH} ${ONES_TMPNIIPATH} 1


# MTRs
${ANTSPATH}/ImageMath 3 ${MTRs_TMPNIIPATH} 		/ ${MTs_AVER_TMPNIIPATH} 	${MT0_AVER_TMPNIIPATH}
${ANTSPATH}/ImageMath 3 ${MTRs_TMPNIIPATH} 		- ${ONES_TMPNIIPATH} 		${MTRs_TMPNIIPATH}

# MTRd
${ANTSPATH}/ImageMath 3 ${MTRd_TMPNIIPATH} 		/ ${MTd_AVER_TMPNIIPATH} 	${MT0_AVER_TMPNIIPATH}
${ANTSPATH}/ImageMath 3 ${MTRd_TMPNIIPATH} 		- ${ONES_TMPNIIPATH} 		${MTRd_TMPNIIPATH}

# ihMTR
${ANTSPATH}/ImageMath 3 ${ihMTR_TMPNIIPATH} 	- ${MTRd_TMPNIIPATH} 		${MTRs_TMPNIIPATH}
${ANTSPATH}/ImageMath 3 ${ihMTR_TMPNIIPATH} 	+ ${ihMTR_TMPNIIPATH} 		${ihMTR_TMPNIIPATH}

# MTRsinv
${ANTSPATH}/ImageMath 3 ${MTRsinv_TMPNIIPATH} 	/ ${MT0_AVER_TMPNIIPATH} 	${MTs_AVER_TMPNIIPATH}
${ANTSPATH}/ImageMath 3 ${MTRsinv_TMPNIIPATH} 	- ${MTRsinv_TMPNIIPATH} 	${ONES_TMPNIIPATH}

# MTRdinv
${ANTSPATH}/ImageMath 3 ${MTRdinv_TMPNIIPATH} 	/ ${MT0_AVER_TMPNIIPATH} 	${MTd_AVER_TMPNIIPATH}
${ANTSPATH}/ImageMath 3 ${MTRdinv_TMPNIIPATH} 	- ${MTRdinv_TMPNIIPATH} 	${ONES_TMPNIIPATH}
 
# ihMTRinv
${ANTSPATH}/ImageMath 3 ${ihMTRinv_TMPNIIPATH} 	- ${MTRdinv_TMPNIIPATH} 	${MTRsinv_TMPNIIPATH}
${ANTSPATH}/ImageMath 3 ${ihMTRinv_TMPNIIPATH} 	+ ${ihMTRinv_TMPNIIPATH} 	${ihMTRinv_TMPNIIPATH}


## Outputs
if [[ ${FLAG_MAP_ihMT} -eq 1 ]]; then
	OUTPUT_ihMTNIIPATH=${OUTPUT_NIIPATH}ihMT.nii
	cp ${ihMT_TMPNIIPATH} ${OUTPUT_ihMTNIIPATH}
fi

declare -a NII_TO_DENOISE
if [[ ${FLAG_MAP_MTRs} -eq 1 ]]; then
	OUTPUT_MTRsNIIPATH=${OUTPUT_NIIPATH}MTRs.nii
	${ANTSPATH}/ThresholdImage  3 ${MTRs_TMPNIIPATH}  				${MASK_TMPNIIPATH} 		0 1
	${ANTSPATH}/ImageMath 		3 ${OUTPUT_MTRsNIIPATH} 		m 	${MTRs_TMPNIIPATH} 		${MASK_TMPNIIPATH}
	NII_TO_DENOISE=( ${NII_TO_DENOISE[@]} ${OUTPUT_MTRsNIIPATH} )
fi

if [[ ${FLAG_MAP_MTRd} -eq 1 ]]; then
	OUTPUT_MTRdNIIPATH=${OUTPUT_NIIPATH}MTRd.nii
	${ANTSPATH}/ThresholdImage  3 ${MTRd_TMPNIIPATH}  				${MASK_TMPNIIPATH} 		0 1
	${ANTSPATH}/ImageMath 		3 ${OUTPUT_MTRdNIIPATH} 		m 	${MTRd_TMPNIIPATH} 		${MASK_TMPNIIPATH}
	NII_TO_DENOISE=( ${NII_TO_DENOISE[@]} ${OUTPUT_MTRdNIIPATH} )
fi

if [[ ${FLAG_MAP_ihMTR} -eq 1 ]]; then
	OUTPUT_ihMTR_NIIPATH=${OUTPUT_NIIPATH}ihMTR.nii
	${ANTSPATH}/ThresholdImage  3 ${ihMTR_TMPNIIPATH}  				${MASK_TMPNIIPATH} 		0 1
	${ANTSPATH}/ImageMath 		3 ${OUTPUT_ihMTR_NIIPATH} 		m 	${ihMTR_TMPNIIPATH} 	${MASK_TMPNIIPATH}
	NII_TO_DENOISE=( ${NII_TO_DENOISE[@]} ${OUTPUT_ihMTR_NIIPATH} )
fi

if [[ ${FLAG_MAP_MTRsinv} -eq 1 ]]; then
	OUTPUT_MTRsinvNIIPATH=${OUTPUT_NIIPATH}MTRsinv.nii
	${ANTSPATH}/ThresholdImage  3 ${MTRsinv_TMPNIIPATH}  			${MASK_TMPNIIPATH} 		0 1
	${ANTSPATH}/ImageMath 		3 ${OUTPUT_MTRsinvNIIPATH} 		m 	${MTRsinv_TMPNIIPATH} 	${MASK_TMPNIIPATH}
	NII_TO_DENOISE=( ${NII_TO_DENOISE[@]} ${OUTPUT_MTRsinvNIIPATH} )
fi

if [[ ${FLAG_MAP_MTRdinv} -eq 1 ]]; then
	OUTPUT_MTRdinvNIIPATH=${OUTPUT_NIIPATH}MTRdinv.nii
	${ANTSPATH}/ThresholdImage  3 ${MTRdinv_TMPNIIPATH}  			${MASK_TMPNIIPATH} 		0 1
	${ANTSPATH}/ImageMath 		3 ${OUTPUT_MTRdinvNIIPATH} 		m 	${MTRdinv_TMPNIIPATH} 	${MASK_TMPNIIPATH}
	NII_TO_DENOISE=( ${NII_TO_DENOISE[@]} ${OUTPUT_MTRdinvNIIPATH} )
fi

if [[ ${FLAG_MAP_ihMTRinv} -eq 1 ]]; then
	OUTPUT_ihMTRinvNIIPATH=${OUTPUT_NIIPATH}ihMTRinv.nii
	${ANTSPATH}/ThresholdImage  3 ${ihMTRinv_TMPNIIPATH}  			${MASK_TMPNIIPATH} 		0 1
	${ANTSPATH}/ImageMath 		3 ${OUTPUT_ihMTRinvNIIPATH} 	m 	${ihMTRinv_TMPNIIPATH} 	${MASK_TMPNIIPATH}
	NII_TO_DENOISE=( ${NII_TO_DENOISE[@]} ${OUTPUT_ihMTRinvNIIPATH} )
fi

# Denoise ihMT-derived maps (BM4D)
func_denoiseBM4D NII_TO_DENOISE FLAG_DENOISE_MAPS


## Clean (or not)
if [ ${FLAG_KEEP_TMP} -eq 0 ]; then
	cleanup
fi