"""
calc_ihMTmaps_cli.py - Computes ihMT-derived maps: ihMTR, ihMT, MTRs, MTRd, MTNs, MTNd

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

###################################################################
############## Argument parsing
################################################################### 
VALID_MAPS = {"ihMT", "ihMTR", "MTRs", "MTRd", "MTNs", "MTNd"}

text_description = f"""\
  calc-ihMTmaps path/to/ihMT_preproc.nii \\
                --ihMTR path/to/ihMTR.nii.gz \\
                --ihMT path/to/ihMT.nii.gz \\
                --MTRs path/to/MTRs.nii.gz \\
                --MTRd path/to/MTRd.nii.gz \\
                --MTNs path/to/MTNs.nii.gz \\
                --MTNd path/to/MTNd.nii.gz \\
                --idx_mt0 1 --idx_mts 2,4 --idx_mtd 3,5 --scale_int 1000

Available output maps:
    ihMTR  : 2 * (MTRd - MTRs)
    ihMT   : MTs - MTd
    MTRs   : 1 - MTs/MT0
    MTRd   : 1 - MTd/MT0
    MTNs   : MTs/MT0
    MTNd   : MTd/MT0
"""

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=text_description, formatter_class=RawTextHelpFormatter)

    # Positional
    parser.add_argument("input",    help="Input 4D ihMT NIfTI image.")
    parser.add_argument("--ihMTR",  default=None, help="Output path to ihMTR map NIfTI image.")
    parser.add_argument("--MTRd",   default=None, help="Output path to Dual MTR map NIfTI image.")
    parser.add_argument("--MTRs",   default=None, help="Output path to Single MTR map NIfTI image.")
    parser.add_argument("--MTNs",   default=None, help="Output path to Normalized MT-Single NIfTI image.")
    parser.add_argument("--MTNd",   default=None, help="Output path to Normalized MT-Dual NIfTI image.")
    parser.add_argument("--ihMT",   default=None, help="Output path to ihMT NIfTI image.")

    # Optional
    parser.add_argument("--idx_mt0", "-R",  default=None, help="Comma-separated 1-based indices of MT reference (MT0) volumes.")
    parser.add_argument("--idx_mts", "-S",  default=None, help="Comma-separated 1-based indices of MT single volumes.")
    parser.add_argument("--idx_mtd", "-D",  default=None, help="Comma-separated 1-based indices of MT dual volumes.")
    parser.add_argument("--scale_int", "-I", default=1.0, type=np.float32, help="Scaling factor of ihMT-derived maps -- convert to uint16 if >1 for data compression benefit (default: 1).")

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

    # Outputs
    requested_maps = []
    for name in VALID_MAPS:
        value = getattr(args, name)
        if value is not None:
            requested_maps.append(name)
            v[name] = value
        else:
            v[name] = None
    if not requested_maps:
        parser.error("At least one map (ihMTR, etc.) must be requested.")
    v["requested_maps"] = set(requested_maps)

    # Volume indices
    n_provided = sum([args.idx_mt0 is not None,
                      args.idx_mts is not None,
                      args.idx_mtd is not None])
    if n_provided not in (0, 3):
        parser.error("Either all three of --idx_mt0, --idx_mts, --idx_mtd must be provided, or none of them.")
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

    # Output as uint16
    if args.scale_int < 1.0:
        parser.error("--scale_int should have value >= 1.0.")
    v["scale_int"] = args.scale_int

    return v

# Helpers
def _nib_save(data: np.ndarray, ref_nii: nibabel.Nifti1Image, path: Path) -> None:
    data    = np.asarray(data).copy() # prevent issues with memory-mapping (Bus error (core-dumped))
    img     = nibabel.Nifti1Image(data, ref_nii.affine, ref_nii.header)
    img.set_data_dtype(data.dtype)
    nibabel.save(img, str(path))
def _average_volumes(data: np.ndarray, indices_1based: list[int]) -> np.ndarray:
    # Return the voxel-wise mean of selected volumes (1-based indices).
    vols = np.stack([data[..., i - 1] for i in indices_1based], axis=-1)
    return vols.mean(axis=-1)

###################################################################
############## main
################################################################### 
def main() -> None:
    # Re-add all arguments & validate
    args    = parse_args()
    v       = validate_args(args, argparse.ArgumentParser(description=text_description, formatter_class=RawTextHelpFormatter))

    print("ihMT map calculation...")

    # Determine 4th-dim size and default indices
    nii_obj = nibabel.load(v["input_path"])
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

    # Map computation
    ihMTp = nii_obj.get_fdata(dtype=np.float32)

    ## Average per contrast
    MT0_avg = _average_volumes(ihMTp, v["idx_mt0"])
    MTs_avg = _average_volumes(ihMTp, v["idx_mts"])
    MTd_avg = _average_volumes(ihMTp, v["idx_mtd"])

    ## Maps
    MTNs  = np.divide(MTs_avg, MT0_avg, out=np.zeros_like(MTs_avg), where=MT0_avg != 0)  # MTs/MT0
    MTNd  = np.divide(MTd_avg, MT0_avg, out=np.zeros_like(MTd_avg), where=MT0_avg != 0)  # MTd/MT0
    MTRs  = 1.0 - MTNs                                                                   # 1 - MTs/MT0
    MTRd  = 1.0 - MTNd                                                                   # 1 - MTd/MT0
    ihMTR = 2.0 * (MTRd - MTRs)                                                          # 2*(MTRd - MTRs)
    ihMT  = 2.0 * (MTs_avg - MTd_avg)                                                    # 2*(MTs - MTd)

    ## Save
    def _write_outputs(arr: np.ndarray, path: Path) -> None:
        arr  = arr * (arr > 0).astype(np.float32) # keep >0 values
        arr  = arr * (arr < 1).astype(np.float32) # keep <1 values
        arr  = (arr * v["scale_int"]).astype(np.uint16) if v["scale_int"] > 1.0 else arr
        _nib_save(arr, nii_obj, path)
        print(f"    Saved: {path}")

    if "MTRs"   in v["requested_maps"]: _write_outputs(MTRs,  v["MTRs"])
    if "MTRd"   in v["requested_maps"]: _write_outputs(MTRd,  v["MTRd"])
    if "ihMTR"  in v["requested_maps"]: _write_outputs(ihMTR, v["ihMTR"])
    if "MTNs"   in v["requested_maps"]: _write_outputs(MTNs,  v["MTNs"])
    if "MTNd"   in v["requested_maps"]: _write_outputs(MTNd,  v["MTNd"])
    if "ihMT"   in v["requested_maps"]: _nib_save(ihMT, nii_obj, v["ihMT"])


    print("ihMT map calculation: Done.\n")

    print("References:\n"
    "- Soustelle et al., A Motion Correction Strategy for Multi-Contrast based 3D parametric imaging: Application to Inhomogeneous Magnetization Transfer (ihMT) bioRxiv, 2020.\n"
    "DOI: 10.1101/2020.09.11.292649 \n"
    "- Soustelle et al., A strategy to reduce the sensitivity of inhomogeneous magnetization transfer (ihMT) imaging to radiofrequency transmit field variations at 3 T Magnetic Resonance in Medicine, 2022. \n"
    "DOI: 10.1002/mrm.29055\n")

if __name__ == "__main__":
    main()