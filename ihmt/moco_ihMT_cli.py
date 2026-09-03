#!/usr/bin/env python3
"""
ihMT_MoCo.py - Motion Correction for inhomogeneous Magnetization Transfer (ihMT) MRI data.

Performs:
  1) N4 bias field correction on each 3D volume
  2) Skull-stripping (nipreps-synthstrip from niprep)
  3) Iterative rigid registration of MTw volumes to build a target
  4) Registration of all volumes onto the final target
  5) Apply transforms to original (non-N4, non-stripped) volumes
  6) Reassemble corrected 4D stack

ANTs operations are implemented via ANTsPy.
nipreps-synthstrip is called as a subprocess.

References:
  - Soustelle et al., A Motion Correction Strategy for Multi-Contrast based 3D parametric imaging: Application to Inhomogeneous Magnetization Transfer (ihMT), bioRxiv, 2020. 
    DOI: 10.1101/2020.09.11.292649
"""

import argparse
import os
import shutil
import signal
import subprocess
import sys
from argparse import RawTextHelpFormatter
from pathlib import Path
from datetime import datetime
import urllib.request

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

def _signal_handler(sig, frame) -> None:  # noqa: ANN001
    if not flag_keep_tmp:
        cleanup()
    sys.exit(0)


signal.signal(signal.SIGINT, _signal_handler)
signal.signal(signal.SIGTERM, _signal_handler)


###################################################################
############## Argument parsing
################################################################### 
DESCRIPTION = """\
Motion Correction pipeline for inhomogeneous Magnetization Transfer (ihMT) MRI data.

Example:
  ihMT_MoCo.py path/to/ihMT_in.nii path/to/ihMT_out.nii \\
      --idx_mt0 1 --idx_mts 2,4 --idx_mtd 3,5 \\
      --nthreads 20
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=RawTextHelpFormatter)

    # Positional
    parser.add_argument("input",help="Input 4D ihMT NIfTI image.")
    parser.add_argument("output",help="Output 4D motion-corrected ihMT NIfTI image.")

    # Optional
    parser.add_argument("--idx_mt0", "-R",default=None,help="Comma-separated 1-based indices of MT0 (reference) volumes\n (default: 1).")
    parser.add_argument("--idx_mts", "-S",default=None,help="Comma-separated 1-based indices of MTs (single) volumes\n (default: 2,4,6,...,N-1).")
    parser.add_argument("--idx_mtd", "-D",default=None,help="Comma-separated 1-based indices of MTd (dual) volumes\n (default: 3,5,7,...,N).")
    parser.add_argument("--nthreads", "-n", type=int, default=1, help="Number of threads for ANTs operations (default: 1).")
    parser.add_argument("--keep_tmp", "-k",action="store_true",help="Keep temporary files.")
    parser.add_argument("--verbose", "-v",action="store_true",help="High verbosity mode.")

    return parser.parse_args()


###################################################################
############## Validation helpers
################################################################### 
def _parse_int_list(s: str, name: str, parser: argparse.ArgumentParser) -> list[int]:
    try:
        return [int(v) for v in s.split(",")]
    except ValueError:
        parser.error(f"--{name} must be a comma-separated list of integers (e.g. 1,2,3).")


def _parse_path_list(s: str, name: str, parser: argparse.ArgumentParser) -> list[Path]:
    paths = [Path(p.strip()) for p in s.split(",")]
    if len(paths) == 2 or len(paths) > 3:
        parser.error(
            f"--{name}: provide either 1 (common) or 3 (one per contrast) paths "
            f"(found {len(paths)})."
        )
    for p in paths:
        if not p.is_file():
            parser.error(f"--{name}: file not found: {p}")
    return paths


def validate_args(args: argparse.Namespace,
                  parser: argparse.ArgumentParser) -> dict:
    v: dict = {}

    # Input
    input_path = Path(args.input)
    if not input_path.is_file():
        parser.error(f"Input not found: {input_path}")
    v["input_path"] = input_path

    # Output
    output_path = Path(args.output)
    if not output_path.parent.is_dir():
        parser.error(f"Output directory does not exist: {output_path.parent}")
    v["output_path"] = output_path

    # Threads
    v["nthreads"] = min(get_physCPU_number(), args.nthreads)

    # Skull-stripping method flag
    v["synthstrip_weights"] = get_synthstrip_weights()

    # Volume indices
    n_provided = sum([args.idx_mt0 is not None,
                      args.idx_mts is not None,
                      args.idx_mtd is not None])
    if n_provided not in (0, 3):
        parser.error("Either all three of --idx_mt0, --idx_mts, --idx_mtd must be provided, or none (auto-derived from volume count).")
    v["custom_idx"] = n_provided == 3
    if v["custom_idx"]:
        idx_mt0 = _parse_int_list(args.idx_mt0, "idx_mt0", parser)
        idx_mts = _parse_int_list(args.idx_mts, "idx_mts", parser)
        idx_mtd = _parse_int_list(args.idx_mtd, "idx_mtd", parser)
        all_idx = idx_mt0 + idx_mts + idx_mtd
        if len(all_idx) != len(set(all_idx)):
            parser.error("Duplicate indices detected across --idx_mt0/--idx_mts/--idx_mtd.")
        v["idx_mt0"] = idx_mt0
        v["idx_mts"] = idx_mts
        v["idx_mtd"] = idx_mtd
    else:
        v["idx_mt0"] = None
        v["idx_mts"] = None
        v["idx_mtd"] = None

    # Keep tmp
    global flag_keep_tmp 
    flag_keep_tmp = True if args.keep_tmp else False

    # Verbosity
    v["verbose"] = True if args.verbose else False

    return v

###################################################################
############## NIfTI helpers
################################################################### 
def _nib_load(path: Path) -> tuple[np.ndarray, nibabel.Nifti1Image]:
    nii = nibabel.load(str(path))
    return nii.get_fdata(dtype=np.float32), nii
def _nib_save(data: np.ndarray, ref_nii: nibabel.Nifti1Image, path: Path) -> None:
    data    = np.asarray(data, dtype=np.float32).copy()
    img     = nibabel.Nifti1Image(data, ref_nii.affine, ref_nii.header)
    img.set_data_dtype(data.dtype)
    nibabel.save(img, str(path))

# Step: N4 bias field correction (native antsRegistration suite)
def step_n4(vol: np.ndarray, ref_nii: nibabel.Nifti1Image, vol_path: Path, out_path: Path, verbose: bool = False) -> np.ndarray:
    # N4 bias field correction on a single 3D volume.
    # Equivalent to: N4BiasFieldCorrection -d 3 -s 2 -v 1 -c [50x50x50x50,1e-10]
    _nib_save(vol.astype(np.float32), ref_nii, vol_path)
    subprocess.run(
        ["N4BiasFieldCorrection", "-d", "3",
                                  "-i", str(vol_path), "-o", str(out_path),
                                  "-s", "2", "-c", "[50x50x50x50,0.0000000001]",
                                  "-v", "1" if verbose else "0",
        ],
        check=True,
    )
    n4_vol, _ = _nib_load(out_path)
    return n4_vol

# Step: skull-stripping
def step_synthstrip(vol: np.ndarray, ref_nii: nibabel.Nifti1Image, vol_path: Path, out_path: Path, synthstrip_weights: Path, nthreads: int, verbose: bool = False) -> np.ndarray:
    # Brain extraction via nipreps-synthstrip (subprocess, reads/writes NIfTI).
    _nib_save(vol.astype(np.float32), ref_nii, vol_path)
    subprocess.run(
        ["nipreps-synthstrip", "-i", str(vol_path), "-o", str(out_path), 
                                "--model", synthstrip_weights, 
                                "--num-threads", str(nthreads)],
        check=True,
        stdout=None if verbose else subprocess.DEVNULL, 
        stderr=None if verbose else subprocess.DEVNULL
    )
    extracted, _ = _nib_load(out_path)
    return extracted

# Step: rigid registration (native antsRegistration)
def _rigid_register(fixed_path: Path, moving_path: Path, out_prefix: str, verbose: bool = False) -> dict:
    # Rigid registration of moving onto fixed.
    # Returns a dict with 'warpedmovout' (np.ndarray) and 'fwdtransforms' (list[str]).
    warped_path = f"{out_prefix}Warped.nii.gz"
    cmd =   [
                "antsRegistration", "--dimensionality", "3",
                "--output", f"[{out_prefix},{warped_path}]",
                "--interpolation", "Linear",
                "--winsorize-image-intensities", "[0.005,0.995]",
                "--use-histogram-matching", "0",
                "--initial-moving-transform", f"[{fixed_path},{moving_path},1]",
                "--transform", "Rigid[0.01]",
                "--metric", f"MI[{fixed_path},{moving_path},1,32,Regular,0.25]",
                "--convergence", "[500x250x100,1e-6,10]",
                "--shrink-factors", "4x2x1",
                "--smoothing-sigmas", "2x1x0vox",
                "--write-composite-transform", "0",
                "--verbose", "1" if verbose else "0", 
            ]
    # cmd.extend(["--random-seed", "12345"])
    subprocess.run(cmd, check=True)
    warped, _ = _nib_load(Path(warped_path))
    return {"warpedmovout": warped, "fwdtransforms": [f"{out_prefix}0GenericAffine.mat"]}

def _apply_transform(fixed_path: Path,moving_path: Path,transformlist: list[str], out_path: Path, verbose: bool = False) -> np.ndarray:
    cmd = ["antsApplyTransforms", "-d", "3", "-i", str(moving_path), "-r", str(fixed_path),
                                  "-o", str(out_path), "-n", "Linear", "--verbose", "1" if verbose else "0"]
    for t in transformlist:
        cmd += ["-t", t]
    subprocess.run(cmd, check=True)
    warped, _ = _nib_load(out_path)
    return warped

###################################################################
############## Get mri_synthstrip weights path
################################################################### 
def get_synthstrip_weights() -> Path:
    # Make sure synthstrip.1.pt is present next to the nipreps-synthstrip binary.
    # Try downloading it on first use if it isn't there yet.
    # Returns the path to the weights file.
    SYNTHSTRIP_WEIGHTS_URL = ("https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/requirements/synthstrip.1.pt")
    SYNTHSTRIP_WEIGHTS_NAME = "synthstrip.1.pt"
    bin_path = shutil.which("nipreps-synthstrip")
    if bin_path is None:
        raise RuntimeError("nipreps-synthstrip not found in PATH; check installation.")

    dest = Path(bin_path).parent / SYNTHSTRIP_WEIGHTS_NAME
    if dest.is_file():
        return dest

    print(
        "\n[ihMT-MoCo] -------------------------------------------------------\n"
        "[ihMT-MoCo] synthstrip.1.pt (niprep-synthstrip requirement) weight not found...\n"
        "[ihMT-MoCo] synthstrip.1.pt is part of the FreeSurfer software.\n"
        "[ihMT-MoCo] By downloading this file, you acknowledge the FreeSurfer\n"
        "[ihMT-MoCo] licence (MIT): https://choosealicense.com/licenses/mit/\n"
        "[ihMT-MoCo] More information at:\n"
        "[ihMT-MoCo]   https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/\n"
        "[ihMT-MoCo] -------------------------------------------------------\n"
    )
    print(f"[ihMT-MoCo] Downloading synthstrip weights -> {dest}...")
    try:
        urllib.request.urlretrieve(SYNTHSTRIP_WEIGHTS_URL, dest)
        print(f"[ihMT-MoCo] synthstrip weights downloaded: {dest}")
    except Exception as exc:
        raise RuntimeError(
            f"Could not download synthstrip weights: {exc}\n"
            f"  Download manually with:\n"
            f"    wget {SYNTHSTRIP_WEIGHTS_URL} -O {dest}"
        ) from exc

    return dest

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

    parser  = argparse.ArgumentParser(description=DESCRIPTION,formatter_class=RawTextHelpFormatter)
    args    = parse_args()
    v       = validate_args(args, parser)

    print("ihMT-MoCo...")
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(v["nthreads"])
    print(f"Processing with {v['nthreads']} thread(s)")

    # Temporary directory
    input_path: Path  = v["input_path"]
    output_path: Path = v["output_path"]
    input_stem = input_path.name.replace(".nii.gz", "").replace(".nii", "")
    tmp_fld = input_path.parent / f"tmp_MoCo_{datetime.now().strftime('%y%m%d-%H%M%S-%f')[:-3]}"
    tmp_fld.mkdir(parents=True, exist_ok=True)
    print(f"Temporary folder: {tmp_fld}")

    # Load 4D input
    data_4d, ref_nii = _nib_load(input_path)
    n_vols = data_4d.shape[3]
    print(f"Input shape: {data_4d.shape}")

    # Derive indices if not provided
    if not v["custom_idx"]:
        v["idx_mt0"] = [1]
        v["idx_mts"] = list(range(2, n_vols + 1, 2))
        v["idx_mtd"] = list(range(3, n_vols + 1, 2))
        print(f"Auto-derived indices:")
        print(f"  MT0 : {v['idx_mt0']}")
        print(f"  MTs : {v['idx_mts']}")
        print(f"  MTd : {v['idx_mtd']}")

    idx_mt0: list[int] = v["idx_mt0"]
    idx_mts: list[int] = v["idx_mts"]
    idx_mtd: list[int] = v["idx_mtd"]

    # Target volumes are MTs + MTd (1-based)
    idx_target = idx_mts + idx_mtd

    # Step 1 - N4 + Brain extraction for each volume
    print("\n--- ihMT-MoCo - Step 1: N4 bias field correction + skull-stripping")
    n4_extracted: list[np.ndarray] = []
    n4_extracted_paths: list[Path] = []

    for ii in range(n_vols):
        print(f"  Volume {ii + 1}/{n_vols}...")

        # N4
        vol_raw_path = tmp_fld / f"vol_{ii:04d}_raw.nii"
        vol_nii_path = tmp_fld / f"vol_{ii:04d}_N4.nii"
        n4_vol = step_n4(data_4d[...,ii], ref_nii, vol_raw_path, vol_nii_path, verbose=v["verbose"])

        # Skull-stripping
        vol_ext_path = tmp_fld / f"vol_{ii:04d}_N4E.nii"
        extracted = step_synthstrip(n4_vol, ref_nii, vol_nii_path, vol_ext_path, v["synthstrip_weights"], v["nthreads"], verbose=v["verbose"])

        n4_extracted.append(extracted)
        n4_extracted_paths.append(vol_ext_path)

    print("--- ihMT-MoCo - Step 1: done\n")

    # Step 2 - Iterative rigid registration of target volumes -> build target
    print("--- ihMT-MoCo - Step 2: Iterative registration to build target")
    target_path = tmp_fld / f"{input_stem}_N4E_target.nii"
    warped_to_avg: list[np.ndarray] = []

    for cnt, ii_1based in enumerate(idx_target):
        ii = ii_1based - 1  # 0-based
        moving_path = n4_extracted_paths[ii]

        if cnt == 0:
            fixed_path = moving_path
        else:
            # Average of previously warped volumes
            avg_data = np.stack(warped_to_avg, axis=-1).mean(axis=-1)
            _nib_save(avg_data.astype(np.float32), ref_nii, target_path)
            fixed_path = target_path

        out_prefix  = str(tmp_fld / f"vol_{ii:04d}_N4E_")

        print(f"  Registering volume {ii_1based} (target iteration {cnt + 1}) ...")
        reg = _rigid_register(fixed_path, moving_path, out_prefix, verbose=v["verbose"])

        warped_to_avg.append(reg["warpedmovout"])

    # Final target = average of all warped target volumes
    avg_data = np.stack(warped_to_avg, axis=-1).mean(axis=-1)
    _nib_save(avg_data.astype(np.float32), ref_nii, target_path)
    print("--- ihMT-MoCo - Step 2: done\n")

    # Step 3 - Register ALL volumes onto final target; apply to original vols
    print("--- ihMT-MoCo - Step 3: Register all volumes onto final target")
    corrected_vols: list[np.ndarray] = []

    for ii in range(n_vols):
        print(f"  Volume {ii + 1}/{n_vols} ...")
        moving_n4e_path = n4_extracted_paths[ii]
        out_prefix      = str(tmp_fld / f"vol_{ii:04d}_N4E_final_")

        # Register N4+extracted volume -> target (to get the transform)
        reg = _rigid_register(target_path, moving_n4e_path, out_prefix, verbose=v["verbose"])

        # Apply the same transform to the ORIGINAL (non-N4, non-stripped) volume
        vol_orig_path    = tmp_fld / f"vol_{ii:04d}_orig.nii"
        warped_orig_path = tmp_fld / f"vol_{ii:04d}_orig_warped.nii"
        _nib_save(data_4d[...,ii].astype(np.float32), ref_nii, vol_orig_path)
        warped_orig = _apply_transform(target_path, vol_orig_path, reg["fwdtransforms"], warped_orig_path, verbose=v["verbose"])
        corrected_vols.append(warped_orig)

    print("--- ihMT-MoCo - Step 3: done\n")

    # Step 4 - Reassemble corrected 4D stack
    print("--- ihMT-MoCo - Step 4: Reassembling 4D output")
    corrected_4d = np.stack(corrected_vols, axis=-1).astype(np.float32)
    _nib_save(corrected_4d, ref_nii, output_path)
    print(f"  Output written: {output_path}")
    print("--- ihMT-MoCo - Step 4: done\n")

    # Cleanup
    if not flag_keep_tmp:
        cleanup()
        print("Temporary files removed.")
    else:
        print(f"Temporary files kept at: {tmp_fld}")

    print("\nihMT-MoCo: done")


if __name__ == "__main__":
    main()