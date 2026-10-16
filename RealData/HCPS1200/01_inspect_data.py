#!/usr/bin/env python3
"""Inspect HCP S1200 ZIP layout and CIFTI metadata for this project."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import struct
import sys
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


MOTOR_RE = re.compile(r"^(\d+)_3T_tfMRI_MOTOR_preproc\.zip$")
REST_RE = re.compile(r"^(\d+)_3T_rfMRI_REST_fix\.zip$")

REST_RUNS = (
    "rfMRI_REST1_LR",
    "rfMRI_REST1_RL",
    "rfMRI_REST2_LR",
    "rfMRI_REST2_RL",
)
MOTOR_RUNS = ("tfMRI_MOTOR_LR", "tfMRI_MOTOR_RL")


@dataclass
class CiftiSummary:
    subject_id: str
    condition: str
    run: str
    zip_path: str
    member: str
    datatype: int
    bitpix: int
    nifti_dim: Tuple[int, ...]
    vox_offset: int
    series_points: int
    series_step: float
    brain_dim_size: int
    brain_model_count: int
    first_structures: List[Dict[str, str]]


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("inspect_hcp")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)
    return logger


def collect_subjects(folder: Path, pattern: re.Pattern) -> Dict[str, Path]:
    out: Dict[str, Path] = {}
    for path in sorted(folder.glob("*.zip")):
        m = pattern.match(path.name)
        if m:
            out[m.group(1)] = path
    return out


def detect_endian(prefix: bytes) -> str:
    if len(prefix) < 4:
        raise ValueError("Header prefix too short.")
    le = struct.unpack("<i", prefix[0:4])[0]
    be = struct.unpack(">i", prefix[0:4])[0]
    if le == 540:
        return "<"
    if be == 540:
        return ">"
    raise ValueError("Not a NIfTI-2 header (sizeof_hdr != 540).")


def parse_nifti2_header(prefix: bytes) -> Dict[str, object]:
    if len(prefix) < 552:
        raise ValueError("Need at least 552 bytes to parse NIfTI-2 + extension header.")
    endian = detect_endian(prefix)
    datatype = struct.unpack(f"{endian}h", prefix[12:14])[0]
    bitpix = struct.unpack(f"{endian}h", prefix[14:16])[0]
    dim = struct.unpack(f"{endian}8q", prefix[16:80])
    vox_offset = struct.unpack(f"{endian}q", prefix[168:176])[0]
    ext_flag = prefix[540:544]
    if ext_flag[0] != 1:
        raise ValueError("No NIfTI extension flag set; expected CIFTI XML extension.")
    esize, ecode = struct.unpack(f"{endian}ii", prefix[544:552])
    if ecode != 32:
        raise ValueError(f"Unexpected extension code {ecode}; expected CIFTI ecode=32.")
    xml_start = 552
    xml_end = xml_start + esize - 8
    if xml_end > len(prefix):
        raise ValueError("Prefix too short to include complete CIFTI XML extension.")
    xml_bytes = prefix[xml_start:xml_end].rstrip(b"\x00")
    return {
        "endian": endian,
        "datatype": int(datatype),
        "bitpix": int(bitpix),
        "dim": tuple(int(x) for x in dim),
        "vox_offset": int(vox_offset),
        "xml_bytes": xml_bytes,
    }


def find_member(zf: zipfile.ZipFile, run: str, condition: str) -> str:
    if condition == "rest":
        suffix = f"{run}/{run}_Atlas_MSMAll_hp2000_clean.dtseries.nii"
    else:
        suffix = f"{run}/{run}_Atlas_MSMAll.dtseries.nii"
    for name in zf.namelist():
        if name.endswith(suffix):
            return name
    raise FileNotFoundError(f"No CIFTI member ending with '{suffix}' found in ZIP.")


def summarize_member(subject: str, zip_path: Path, condition: str, run: str) -> CiftiSummary:
    with zipfile.ZipFile(zip_path, "r") as zf:
        member = find_member(zf, run=run, condition=condition)
        with zf.open(member, "r") as handle:
            prefix = handle.read(900_000)
    hdr = parse_nifti2_header(prefix)
    root = ET.fromstring(hdr["xml_bytes"])

    series_map = None
    brain_map = None
    for mim in root.findall(".//MatrixIndicesMap"):
        map_type = mim.attrib.get("IndicesMapToDataType", "")
        if map_type == "CIFTI_INDEX_TYPE_SERIES":
            series_map = mim
        elif map_type == "CIFTI_INDEX_TYPE_BRAIN_MODELS":
            brain_map = mim

    if series_map is None or brain_map is None:
        raise ValueError("Could not find required SERIES and BRAIN_MODELS maps in CIFTI XML.")

    series_points = int(series_map.attrib.get("NumberOfSeriesPoints", "0"))
    series_step = float(series_map.attrib.get("SeriesStep", "0"))
    brain_models = brain_map.findall("BrainModel")
    first_structures = []
    for bm in brain_models[:10]:
        first_structures.append(
            {
                "BrainStructure": bm.attrib.get("BrainStructure", ""),
                "IndexOffset": bm.attrib.get("IndexOffset", ""),
                "IndexCount": bm.attrib.get("IndexCount", ""),
                "ModelType": bm.attrib.get("ModelType", ""),
            }
        )

    dim = hdr["dim"]
    brain_dim_size = int(dim[6]) if len(dim) > 6 else 0
    return CiftiSummary(
        subject_id=subject,
        condition=condition,
        run=run,
        zip_path=str(zip_path),
        member=member,
        datatype=int(hdr["datatype"]),
        bitpix=int(hdr["bitpix"]),
        nifti_dim=tuple(int(x) for x in dim),
        vox_offset=int(hdr["vox_offset"]),
        series_points=series_points,
        series_step=series_step,
        brain_dim_size=brain_dim_size,
        brain_model_count=len(brain_models),
        first_structures=first_structures,
    )


def write_inventory_csv(
    path: Path, subjects: List[str], motor: Dict[str, Path], rest: Dict[str, Path]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["subject_id", "has_motor_zip", "has_rest_zip", "motor_zip", "rest_zip"])
        for sid in subjects:
            writer.writerow(
                [
                    sid,
                    int(sid in motor),
                    int(sid in rest),
                    str(motor[sid]) if sid in motor else "",
                    str(rest[sid]) if sid in rest else "",
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Inspect HCP ZIP/CIFTI data layout.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent, help="Application data/output root (default: script directory).")
    parser.add_argument(
        "--inventory-csv",
        type=Path,
        default=Path("data_processed/hcp_zip_inventory.csv"),
        help="Output CSV for subject ZIP inventory.",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=Path("data_processed/hcp_data_inventory.json"),
        help="Output JSON summary file.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/01_inspect_data.log"),
        help="Log file path.",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    logger = setup_logger(root / args.log_file)

    motor_dir = root / "HCP_Data" / "MotorTask"
    rest_dir = root / "HCP_Data" / "RestingState"
    if not motor_dir.exists() or not rest_dir.exists():
        raise FileNotFoundError("Expected HCP_Data/MotorTask and HCP_Data/RestingState directories.")

    motor = collect_subjects(motor_dir, MOTOR_RE)
    rest = collect_subjects(rest_dir, REST_RE)

    all_subjects = sorted(set(motor) | set(rest))
    overlap = sorted(set(motor) & set(rest))
    logger.info("Motor ZIPs: %d", len(motor))
    logger.info("Rest ZIPs: %d", len(rest))
    logger.info("Overlap subjects: %d", len(overlap))

    write_inventory_csv(root / args.inventory_csv, all_subjects, motor, rest)
    logger.info("Wrote ZIP inventory: %s", root / args.inventory_csv)

    sample_entries: List[CiftiSummary] = []
    if overlap:
        sid = overlap[0]
        try:
            sample_entries.append(summarize_member(sid, motor[sid], "motor", MOTOR_RUNS[0]))
            sample_entries.append(summarize_member(sid, rest[sid], "rest", REST_RUNS[0]))
        except Exception as exc:
            logger.exception("Failed to parse sample CIFTI files: %s", exc)
            raise

    summary = {
        "motor_zip_count": len(motor),
        "rest_zip_count": len(rest),
        "overlap_count": len(overlap),
        "first_overlap_subjects": overlap[:20],
        "sample_cifti_summaries": [
            {
                "subject_id": s.subject_id,
                "condition": s.condition,
                "run": s.run,
                "zip_path": s.zip_path,
                "member": s.member,
                "datatype": s.datatype,
                "bitpix": s.bitpix,
                "nifti_dim": s.nifti_dim,
                "vox_offset": s.vox_offset,
                "series_points": s.series_points,
                "series_step": s.series_step,
                "brain_dim_size": s.brain_dim_size,
                "brain_model_count": s.brain_model_count,
                "first_structures": s.first_structures,
            }
            for s in sample_entries
        ],
    }

    (root / args.summary_json).parent.mkdir(parents=True, exist_ok=True)
    with (root / args.summary_json).open("w") as f:
        json.dump(summary, f, indent=2)
    logger.info("Wrote summary JSON: %s", root / args.summary_json)
    logger.info("Inspection completed.")


if __name__ == "__main__":
    main()
