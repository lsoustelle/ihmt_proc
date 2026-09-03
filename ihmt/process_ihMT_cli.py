#!/usr/bin/env python3
"""
process_ihMT_cli.py - Processing pipeline for inhomogeneous Magnetization Transfer (ihMT) MRI data.

Performs (in order, when enabled):
    1) MP-PCA denoising of raw ihMT images
    2) Gibbs-ringing removal
    3) Gradient non-linearity distortion correction
    4) Motion correction
    5) ihMT-derived map computation

References:
  - Soustelle et al., bioRxiv 2020. DOI: 10.1101/2020.09.11.292649
  - Soustelle et al., Magn Reson Med 2022. DOI: 10.1002/mrm.29055
"""

import argparse
import os
import shutil
import signal
import subprocess
import sys
import tempfile
from argparse import RawTextHelpFormatter
from pathlib import Path
from datetime import datetime

import nibabel
import numpy as np

# Globals (populated after argument parsing)
tmp_fld: Path | None = None
flag_keep_tmp: bool = False

# Cleanup helpers
def cleanup() -> None:
    global tmp_fld
    if tmp_fld is not None and tmp_fld.is_dir():
        shutil.rmtree(tmp_fld, ignore_errors=True)

def _signal_handler(sig, frame) -> None: 
    if not flag_keep_tmp:
        cleanup()
    sys.exit(0)

signal.signal(signal.SIGINT, _signal_handler)
signal.signal(signal.SIGTERM, _signal_handler)

###################################################################
############## Argument parsing
################################################################### 
VALID_MAPS = {"ihMTp", "ihMTR", "MTRs", "MTRd", "MTNs", "MTNd"}

MAP_DESCRIPTIONS = """\
  ihMTp     : pre-processed ihMT stack
  ihMTR     : 2 * (MTRd - MTRs)
  MTRs      : 1 - MTs/MT0
  MTRd      : 1 - MTd/MT0
  MTNs      : MTs/MT0
  MTNd      : MTd/MT0"""

text_description = f"""\
Processing pipeline for inhomogeneous Magnetization Transfer (ihMT) MRI data.

Example:
  process_ihMT.py path/to/in_ihMT.nii \\
                  path/to/out_ \\
                  --maps ihMTR,MTRd \\
                  --mppca \\
                  --unring 1 \\
                  --moco 1 \\
                  --idx_mt0 1 --idx_mts 2,4 --idx_mtd 3,5 \\
                  --nthreads 8

Available output maps:
{MAP_DESCRIPTIONS}
"""

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=text_description, formatter_class=RawTextHelpFormatter)

    # Positional
    parser.add_argument("input", help="Input 4D ihMT NIfTI image.")
    parser.add_argument("output_prefix", help="Output path-prefix for 3D ihMT-derived NIfTI images\n(e.g. /path/to/out_ -> /path/to/out_ihMTR.nii).")
    parser.add_argument("--maps", "-c", required=True, metavar="MAP1,MAP2,...", help="Comma-separated list of maps to compute.\n" f"Valid values: {', '.join(VALID_MAPS)}")

    # Optional
    parser.add_argument("--mppca", "-d", action="store_true", help="Perform MP-PCA denoising of raw ihMT images.")
    parser.add_argument("--mppca_window", "-e", help="MP-PCA kernel extent as comma-separated integers (default: 5,5,5).")
    parser.add_argument("--unring", "-u", type=int, choices=[0, 1, 2], default=0, help= "Gibbs-ringing removal method:\n"  
                                                                                        "  0: none (default)\n" 
                                                                                        "  1: 3D Sub-voxel shift unringing\n" 
                                                                                        "  2: cos-kernel apodization & zero-filling ×2")
    parser.add_argument("--gnldc", "-w", action="store_true", help="Perform gradient non-linearity distortion correction.")
    parser.add_argument("--gnldc_grad", "-g", default=None, help="Path to *.grad file for distortion correction.")
    parser.add_argument("--gnldc_ngrid", "-N", type=int, default=60, help="Number of grid points for distortion correction (default: 60).")
    parser.add_argument("--moco", "-m", type=int, choices=[0, 1, 2, 3], default=0, help="Motion correction method:\n"
                                                                                        "  0: none (default)\n"
                                                                                        "  1: ihMT-MoCo\n"
                                                                                        "  2: antsMotionCorr (register to first MT0)\n"
                                                                                        "  3: antsMotionCorr (register to time-series average)")
    parser.add_argument("--idx_mt0", "-R", default=None, help="Comma-separated 1-based indices of MT Reference (MT0) volumes\n(default: 1).")
    parser.add_argument("--idx_mts", "-S", default=None, help="Comma-separated 1-based indices of MT Single (MTs) volumes\n(default: 2,4,6,...,N-1).")
    parser.add_argument("--idx_mtd", "-D", default=None, help="Comma-separated 1-based indices of MT Dual (MTd) volumes\n(default: 3,5,7,...,N).")
    parser.add_argument("--nthreads", "-n", type=int, default=1, help="Number of threads for ANTs operations (default: 1).")
    parser.add_argument("--out_int", "-I", action="store_true", help="x1000 scaling of ihMT-derived maps and convert to int16 (data compression benefit).")
    parser.add_argument("--out_gz",  "-G", action="store_true", help="Write outputs as .nii.gz file(s).")
    parser.add_argument("--keep_tmp", "-k", action="store_true", help="Keep temporary files (default: delete on exit).")
    parser.add_argument("--verbose", "-v", action="store_true",help="High verbosity mode.")

    return parser.parse_args()

###################################################################
############## Validation helpers
################################################################### 
def _parse_int_list(s: str, name: str) -> list[int]:
    # Parse a comma-separated string of integers, raising argparse errors on failure.
    try:
        return [int(v) for v in s.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError(f"--{name} must be a comma-separated list of integers (e.g. 1,2,3).")

def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> dict:
    # Validate and resolve all arguments.  Returns a dict of parsed/resolved values
    # so the rest of the script works with clean Python objects.
    v: dict = {}

    # Input
    v["input_path"] = Path(args.input)

    # Output prefix
    output_prefix = Path(args.output_prefix)
    out_dir = output_prefix.parent
    if not out_dir.is_dir():
        parser.error(f"Output directory does not exist: {out_dir}")
    v["output_prefix"] = output_prefix

    # Maps
    requested_maps = [m.strip() for m in args.maps.split(",")]
    invalid = set(requested_maps) - VALID_MAPS
    if invalid:
        parser.error(
            f"Unrecognised map(s): {sorted(invalid)}. "
            f"Valid values: {sorted(VALID_MAPS)}"
        )
    v["maps"] = set(requested_maps)

    # Threads
    v["nthreads"] = min(get_physCPU_number(), args.nthreads)

    # MP-PCA denoising
    v["flag_mppca"] = True if args.mppca else False
    if v["flag_mppca"] is False and args.mppca_window is not None: 
        parser.error("--mppca_window set but not --mppca, ignoring.")
    if args.mppca_window is None:
        args.mppca_window = "5,5,5" # default
    try:
        v["mppca_window"] = [int(x) for x in args.mppca_window.split(",")]
    except ValueError:
        parser.error("--mppca_window must be comma-separated integers (e.g. 5,5,5).")

    # Unring
    if args.unring not in (0, 1, 2):
        parser.error("--unring must be 0, 1 or 2.")
    v["flag_unring"] = args.unring

    # Gradient non-linearity distortion correction
    v["flag_gnldc"] = True if args.gnldc else False
    if v["flag_gnldc"]:
        if args.gnldc_grad is None:
            parser.error("--gnldc_grad is required when --gnldc is set.")
        gnldc_grad_path = Path(args.gnldc_grad)
        if not gnldc_grad_path.is_file():
            parser.error(f"Gradient file for distortion correction not found: {args.grad}")
        v["gnldc_grad_path"] = gnldc_grad_path
        v["gnldc_n_grid"] = args.gnldc_ngrid

    # MoCo
    if args.moco not in (0, 1, 2, 3):
        parser.error("--moco must be 0, 1, 2, or 3.")
    v["flag_moco"] = args.moco

    # Volume indices
    n_provided = sum([args.idx_mt0 is not None,
                      args.idx_mts is not None,
                      args.idx_mtd is not None])
    if n_provided not in (0, 3):
        parser.error(
            "Either all three of --idx_mt0, --idx_mts, --idx_mtd must be provided, "
            "or none of them (auto-derive from volume count)."
        )
    v["custom_idx"] = n_provided == 3
    if v["custom_idx"]:
        idx_mt0 = _parse_int_list(args.idx_mt0, "idx_mt0")
        idx_mts = _parse_int_list(args.idx_mts, "idx_mts")
        idx_mtd = _parse_int_list(args.idx_mtd, "idx_mtd")
        all_idx = idx_mt0 + idx_mts + idx_mtd
        if len(all_idx) != len(set(all_idx)):
            parser.error("Duplicate indices detected across --idx_mt0/--idx_mts/--idx_mtd.")
        v["idx_mt0"] = idx_mt0
        v["idx_mts"] = idx_mts
        v["idx_mtd"] = idx_mtd
    else:
        # deferred: set after we know the 4th-dim size
        v["idx_mt0"] = None
        v["idx_mts"] = None
        v["idx_mtd"] = None

    # Output as int16
    v["out_int"] = True if args.out_int else False

    # Output as .nii.gz
    v["out_gz"] = True if args.out_gz else False

    # Keep temporary files
    global flag_keep_tmp 
    flag_keep_tmp = True if args.keep_tmp else False

    # Verbosity
    v["verbose"] = True if args.verbose else False

    return v

# nibabel helpers
def _nib_load(path: Path) -> tuple[np.ndarray, nibabel.Nifti1Image]:
    nii = nibabel.load(str(path))
    return nii.get_fdata(dtype=np.float32), nii
def _nib_save(data: np.ndarray, ref_nii: nibabel.Nifti1Image, path: Path) -> None:
    data    = np.asarray(data).copy() # prevent issues with memory-mapping (Bus error (core-dumped))
    img     = nibabel.Nifti1Image(data, ref_nii.affine, ref_nii.header)
    img.set_data_dtype(data.dtype)
    nibabel.save(img, str(path))

# Step helpers
def step_denoise_mppca(nii_path: Path, extent: list[int], nthreads: int) -> None:
    # MP-PCA denoising via denoise-tmppca & returns path to the denoised NIfTI.
    window_str = ",".join(str(e) for e in extent)
    print(f"  MP-PCA denoising (window={window_str})...")
    if shutil.which("denoise-tmppca") is None:
        sys.exit(f"denoise-tmppca wrapper not found; please check installation.")
    subprocess.run(["denoise-tmppca",str(nii_path), str(nii_path),"--window", window_str,"--num_threads", str(nthreads)],check=True)
    print("  MP-PCA denoising: done\n")

def step_unring_svsdegibbs(nii_path: Path, nthreads: int) -> None:
    # Gibbs-ringing removal via svsdegibbs (port from MRtrix3 dev branch for 3D mrdegibbs).
    print("  Unringing (svsdegibbs)...")
    if shutil.which("svsdegibbs") is None:
        sys.exit("svsdegibbs wrapper not found; please check installation.")
    subprocess.run(["svsdegibbs", "-dimensionality", "3", "-nthreads", str(nthreads), str(nii_path), str(nii_path)],check=True)
    print("  Unringing (svsdegibbs): done\n")

def step_unring_apodize(nii_path: Path) -> None:
    # Cos-kernel apodization
    print('  Unringing (cos-kernel apodization)')
    interpfac= 2

    ## load nifti
    img, ref_nii = _nib_load(nii_path)

    ## build cosine kernel                          
    wx  = np.array( np.cos( ( np.arange(img.shape[0]) - np.floor(img.shape[0]/2) ) * np.pi / img.shape[0] ) )
    wy  = np.array( np.cos( ( np.arange(img.shape[1]) - np.floor(img.shape[1]/2) ) * np.pi / img.shape[1] ) )
    wz  = np.array( np.cos( ( np.arange(img.shape[2]) - np.floor(img.shape[2]/2) ) * np.pi / img.shape[2] ) )
    wx  = wx[np.newaxis,:]
    wy  = wy[np.newaxis,:]
    wz  = wz[np.newaxis,np.newaxis,:]
    w2D = np.dot(wx.T,wy)
    w2D = w2D[:,:,np.newaxis]  
    w3D = np.squeeze(np.dot(w2D, wz))

    ## process fft & apodize
    img_apo  = np.zeros(( img.shape[0]*interpfac,img.shape[1]*interpfac,img.shape[2]*interpfac,img.shape[3] ))
    for ii in range(img.shape[3]):
        img_data_tmp    = img[:,:,:,ii]
        img_data_fft    = np.fft.fftshift( np.fft.fftn(img_data_tmp) ) / ( img.shape[0] * img.shape[1] * img.shape[2] )
        img_data_ifft   = np.fft.ifftn( img_data_fft * w3D, \
                                       [ img.shape[0]*2,img.shape[1]*2,img.shape[2]*2 ] ) * \
                                         img.shape[0]*2*img.shape[1]*2*img.shape[2]*2
        img_apo[:,:,:,ii] = abs(img_data_ifft)

    ## modify pixdim & sform in ref & save nifti
    ref_nii.header['pixdim'][1:4] = ref_nii.header['pixdim'][1:4] / interpfac
    _nib_save(img_apo, ref_nii, nii_path)
    subprocess.run( ["ImageMath", "4", nii_path, "+", nii_path, "0"], check=True) # refresh header
    print('  Unringing (cos-kernel apodization): done\n')

def step_gnldc(nii_path: Path, grad_path: Path, n_grid: int, tmp_fld: Path, n_vols: int, verbose: bool = False) -> None:
    # Gradient non-linearity distortion correction (gradient_unwarp.py, in-place).
    print("  Gradient non-linearity distortion correction...")
    if shutil.which("gradient_unwarp.py") is None:
        sys.exit("gradient_unwarp.py wrapper not found; please check installation.")
    data, ref_nii = _nib_load(nii_path)
    for vi in range(n_vols):
        vol_path = tmp_fld / f"vol_{vi:04d}.nii"
        _nib_save(data[..., vi], ref_nii, vol_path)
        print(f"    Correcting volume {vi + 1}/{n_vols}...")
        subprocess.run( [   "gradient_unwarp.py",
                            str(vol_path), str(vol_path),
                            "siemens", "-g", str(grad_path),
                            "--numpoints", str(n_grid),
                            "--interp_order", "1",
                            "--fovmin", "-.300", "--fovmax", ".300",
                            "--verbose"],
                        check=True,
                        stdout=None if verbose else subprocess.DEVNULL, 
                        stderr=None if verbose else subprocess.DEVNULL)
        data[..., vi] = nibabel.load(str(vol_path)).get_fdata()
        vol_path.unlink(missing_ok=True)
    _nib_save(data, ref_nii, nii_path)

    # gradient_unwarp creates this file whatsoever - remove
    warp_path = Path.cwd() / Path("fullWarp_abs.nii.gz")
    warp_path.unlink(missing_ok=True)

    print("  Gradient non-linearity distortion correction: done\n")

def step_moco_ihmt(nii_path: Path, idx_mt0: list[int], idx_mts: list[int], idx_mtd: list[int], nthreads: int, verbose: bool = False) -> None:
    # Motion correction via ihMT-MoCo
    print("  Motion Correction (ihMT-MoCo)...")
    if shutil.which("moco-ihMT") is None:
        sys.exit("moco-ihMT wrapper not found; please check installation.")
    print("  When using ihMT-MoCo, you acknowledge the following EULA: https://crmbm.univ-amu.fr/resources/ihmt-moco/")
    cmd =   [   "moco-ihMT", str(nii_path), str(nii_path),
                "--nthreads", str(nthreads),
                "-R", ",".join(str(i) for i in idx_mt0),
                "-S", ",".join(str(i) for i in idx_mts),
                "-D", ",".join(str(i) for i in idx_mtd),
            ]
    if flag_keep_tmp:
        cmd.append("--keep_tmp")
    if verbose:
        cmd.append("--verbose")
    subprocess.run(cmd, check=True)
    print("  Motion Correction (ihMT-MoCo): done\n")

def step_moco_ants_mt0(nii_path: Path, idx_mt0: list[int], tmp_fld: Path, verbose: bool = False) -> None:
    # antsMotionCorr registration of all volumes onto the first MT0 volume. 
    print("  Motion Correction (antsMotionCorr -> MT0)...")
    if shutil.which("antsMotionCorr") is None:
        sys.exit("antsMotionCorr not found; please check installation.")
    target_path = tmp_fld / "target_mt0.nii"
    prefix_path = tmp_fld / "antsmotioncorr_"

    # save first MT0
    data, ref_nii = _nib_load(nii_path)
    _nib_save(data[..., idx_mt0[0]-1], ref_nii, target_path)

    # run antsMotionCorr
    cmd =   [   "antsMotionCorr", 
                "-d", "3", 
                "-o", f"[{prefix_path},{nii_path}]",
                "-m", f"MI[{target_path},{nii_path},1,32,Regular,0.25]",
                "-t", "Rigid[0.005]",
                "-i", "50x50", "-u", "1", "-s", "1x0", "-f", "2x1", "-e", "1", 
                "--verbose", "1" if verbose else "0",
            ]
    # cmd.extend(["--random-seed", "12345"])
    subprocess.run(cmd, check=True)
    print("  Motion Correction (antsMotionCorr -> MT0): done\n")


def step_moco_ants_avg(nii_path: Path, tmp_fld: Path, verbose: bool = False) -> None:
    # antsMotionCorr registration of all volumes onto the time-series average.
    # Step 1: compute average; Step 2: register each volume to average.
    print("  Motion Correction (antsMotionCorr -> average)...")
    if shutil.which("antsMotionCorr") is None:
        sys.exit("antsMotionCorr not found; please check installation.")
    target_path = tmp_fld / "target_aver.nii"
    prefix_path = tmp_fld / "antsmotioncorr_"

    # save stack's average
    data, ref_nii = _nib_load(nii_path)
    _nib_save(np.mean(data,axis=3), ref_nii, target_path)

    # run antsMotionCorr
    cmd =   [   "antsMotionCorr", 
                "-d", "3", 
                "-o", f"[{prefix_path},{nii_path}]",
                "-m", f"MI[{target_path},{nii_path},1,32,Regular,0.25]",
                "-t", "Rigid[0.005]",
                "-i", "50x50", "-u", "1", "-s", "1x0", "-f", "2x1", "-e", "1", "-v", "1", 
                "--verbose", "1" if verbose else "0",
            ]
    print(cmd)
    # cmd.extend(["--random-seed", "12345"])
    subprocess.run(cmd, check=True)
    print("  Motion Correction (antsMotionCorr -> average): done\n")

def _average_volumes(data: np.ndarray, indices_1based: list[int]) -> np.ndarray:
    # Return the voxel-wise mean of selected volumes (1-based indices).
    vols = np.stack([data[..., i - 1] for i in indices_1based], axis=-1)
    return vols.mean(axis=-1)

def step_compute_maps(nii_path: Path, idx_mt0: list[int], idx_mts: list[int], idx_mtd: list[int],
                      requested_maps: set[str], output_prefix: Path, out_gz: bool, out_int: bool) -> None:
    # Compute all requested ihMT-derived maps and write them to output_prefix<map>.nii.{gz}.
    # ANTs ImageMath / AverageImages / ThresholdImage operations -> numpy.
    print("  Writing outputs...")
    ihMTp, ref_nii = _nib_load(nii_path)

    # Average per contrast
    MT0_avg = _average_volumes(ihMTp, idx_mt0)
    MTs_avg = _average_volumes(ihMTp, idx_mts)
    MTd_avg = _average_volumes(ihMTp, idx_mtd)

    # Maps
    MTRs  = 1.0 - np.divide(MTs_avg, MT0_avg, out=np.zeros_like(MTs_avg), where=MT0_avg != 0)  # 1 - MTs/MT0
    MTRd  = 1.0 - np.divide(MTd_avg, MT0_avg, out=np.zeros_like(MTs_avg), where=MT0_avg != 0)  # 1 - MTd/MT0
    ihMTR = 2.0 * (MTRd - MTRs)                                                                # 2*(MTRd - MTRs)
    MTNs  = np.divide(MTs_avg, MT0_avg,out=np.zeros_like(MTs_avg),where=MT0_avg != 0)          # MTs/MT0
    MTNd  = np.divide(MTd_avg, MT0_avg,out=np.zeros_like(MTs_avg),where=MT0_avg != 0)          # MTd/MT0

    # Save
    def _write_outputs(arr: np.ndarray, name: str) -> None:
        if name is not "ihMTp":
            arr  = arr * (arr > 0).astype(np.float32) # keep >0 values
            arr  = arr * (arr < 1).astype(np.float32) # keep <1 values
            arr  = (arr * 1000).astype(np.int16) if out_int else arr
        ext      = ".nii.gz" if out_gz else ".nii" 
        out_path = Path(str(output_prefix) + name + ext)
        _nib_save(arr, ref_nii, out_path)
        print(f"    Saved: {out_path}")

    if "ihMTp"  in requested_maps: _write_outputs(ihMTp, "ihMTp")
    if "MTRs"   in requested_maps: _write_outputs(MTRs,  "MTRs")
    if "MTRd"   in requested_maps: _write_outputs(MTRd,  "MTRd")
    if "ihMTR"  in requested_maps: _write_outputs(ihMTR, "ihMTR")
    if "MTNs"   in requested_maps: _write_outputs(MTNs,  "MTNs")
    if "MTNd"   in requested_maps: _write_outputs(MTNd,  "MTNd")

    print("  Writing outputs: done")

###################################################################
############## Get CPU info
################################################################### 
def get_physCPU_number():
    # from joblib source code (commit d5c8274)
    # https://github.com/joblib/joblib/blob/master/joblib/externals/loky/backend/context.py#L220-L246
    if sys.platform == "linux":
        cpu_info = subprocess.run(
            "lscpu --parse=core".split(" "), capture_output=True)
        cpu_info = cpu_info.stdout.decode("utf-8").splitlines()
        cpu_info = {line for line in cpu_info if not line.startswith("#")}
        cpu_count_physical = len(cpu_info)
    elif sys.platform == "win32":
        cpu_info = subprocess.run(
            "wmic CPU Get NumberOfCores /Format:csv".split(" "),
            capture_output=True)
        cpu_info = cpu_info.stdout.decode('utf-8').splitlines()
        cpu_info = [l.split(",")[1] for l in cpu_info
                    if (l and l != "Node,NumberOfCores")]
        cpu_count_physical = sum(map(int, cpu_info))
    elif sys.platform == "darwin":
        cpu_info = subprocess.run(
            "sysctl -n hw.physicalcpu".split(" "), capture_output=True)
        cpu_info = cpu_info.stdout.decode('utf-8')
        cpu_count_physical = int(cpu_info)
    else:
        raise NotImplementedError(
            "unsupported platform: {}".format(sys.platform))
    if cpu_count_physical < 1:
            raise ValueError(
                "found {} physical cores < 1".format(cpu_count_physical))
    return cpu_count_physical

###################################################################
############## main
################################################################### 
def main() -> None:
    global tmp_fld, flag_keep_tmp

    # Re-add all arguments & validate
    args    = parse_args()
    v       = validate_args(args, argparse.ArgumentParser(description=text_description, formatter_class=RawTextHelpFormatter))

    print("ihMT preprocessing...")
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(v["nthreads"])
    print(f"Processing with {v['nthreads']} thread(s)")

    # Temporary directory
    tmp_fld = Path(f"/tmp/tmp_ihMTproc_{datetime.now().strftime('%y%m%d-%H%M%S-%f')[:-3]}")
    tmp_fld.mkdir(parents=True, exist_ok=True)
    if not tmp_fld.is_dir():
        sys.exit(f"Cannot create temporary folder: {tmp_fld}")
    print(f"Temporary folder: {tmp_fld}")

    # Input path
    input_path: Path = v["input_path"]

    # Copy to tmp & ensure uncompressed NIfTI
    tmp_nii = tmp_fld / input_path.name
    shutil.copy2(input_path, tmp_nii)
    if tmp_nii.suffix == ".gz":
        import gzip
        uncompressed = tmp_nii.with_suffix("")  # removes .gz, keeps .nii
        with gzip.open(tmp_nii, "rb") as f_in, open(uncompressed, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        tmp_nii.unlink()
        tmp_nii = uncompressed

    # Determine 4th-dim size and default indices
    nii_obj = nibabel.load(str(tmp_nii))
    n_vols = nii_obj.shape[3] if nii_obj.ndim >= 4 else 1
    print(f"Input shape: {nii_obj.shape}\n")

    if not v["custom_idx"]:
        # Default: [MT0, MTs, MTd, ..., MTd]  (1-indexed)
        v["idx_mt0"] = [1]
        v["idx_mts"] = list(range(2, n_vols + 1, 2))
        v["idx_mtd"] = list(range(3, n_vols + 1, 2))
        print(f"Auto-derived indices:")
        print(f"  MT0 : {v['idx_mt0']}")
        print(f"  MTs : {v['idx_mts']}")
        print(f"  MTd : {v['idx_mtd']}")

    # Step 1 - MP-PCA denoising
    if v["flag_mppca"]:
        print("--- Step 1: MP-PCA denoising")
        step_denoise_mppca(tmp_nii, v["mppca_window"], v["nthreads"])
    else:
        print("--- Step 1: MP-PCA denoising (skipped)")

    # Step 2 - Unringing
    if v["flag_unring"] == 1:
        print("--- Step 2: Unringing")
        step_unring_svsdegibbs(tmp_nii, v["nthreads"])
    elif v["flag_unring"] == 2:
        print("--- Step 2: Unringing")
        step_unring_apodize(tmp_nii)
    else:
        print("--- Step 2: Unringing (skipped)")

    # Step 3 - Gradient non-linearity distortion correction
    if v["flag_gnldc"]:
        print("--- Step 3: Gradient non-linearity distortion correction")
        step_gnldc(tmp_nii, v["gnldc_grad_path"], v["gnldc_n_grid"], tmp_fld, n_vols, verbose=v["verbose"])
    else:
        print("--- Step 3: Gradient non-linearity distortion correction (skipped)")

    # Step 4 - Motion correction
    if v["flag_moco"] == 1:
        print("--- Step 4: Motion Correction (ihMT-MoCo)")
        step_moco_ihmt(tmp_nii,v["idx_mt0"], v["idx_mts"], v["idx_mtd"], v["nthreads"], verbose=v["verbose"])
    elif v["flag_moco"] == 2:
        print("--- Step 4: Motion Correction (antsMotionCorr -> MT0)")
        step_moco_ants_mt0(tmp_nii, v["idx_mt0"], tmp_fld, verbose=v["verbose"])
    elif v["flag_moco"] == 3:
        print("--- Step 4: Motion Correction (antsMotionCorr -> average)")
        step_moco_ants_avg(tmp_nii, tmp_fld, verbose=v["verbose"])
    else:
        print("--- Step 4: Motion Correction (skipped)")

    # Step 5 - Map computation
    print("--- Step 5: Map computation")
    step_compute_maps(tmp_nii, v["idx_mt0"], v["idx_mts"], v["idx_mtd"], v["maps"], v["output_prefix"], v["out_gz"], v["out_int"])

    # Cleanup
    if not flag_keep_tmp:
        cleanup()
        print("\nTemporary files removed.")
    else:
        tmp_dest = v["output_prefix"].parent / tmp_fld.name
        shutil.move(tmp_fld, tmp_dest)
        print(f"\nTemporary files kept at: {tmp_dest}")

    print("ihMT preprocessing: Done.\n")

    print("References:\n"
    "- Soustelle et al., A Motion Correction Strategy for Multi-Contrast based 3D parametric imaging: Application to Inhomogeneous Magnetization Transfer (ihMT) bioRxiv, 2020.\n"
    "DOI: 10.1101/2020.09.11.292649 \n"
    "- Soustelle et al., A strategy to reduce the sensitivity of inhomogeneous magnetization transfer (ihMT) imaging to radiofrequency transmit field variations at 3 T Magnetic Resonance in Medicine, 2022. \n"
    "DOI: 10.1002/mrm.29055\n")

if __name__ == "__main__":
    main()