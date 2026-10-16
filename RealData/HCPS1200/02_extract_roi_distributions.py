#!/usr/bin/env python3
"""Extract ROI-wise distribution summaries from HCP CIFTI dtseries ZIP archives."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import random
import re
import struct
import sys
import xml.etree.ElementTree as ET
import zipfile
from array import array
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


MOTOR_RE = re.compile(r"^(\d+)_3T_tfMRI_MOTOR_preproc\.zip$")
REST_RE = re.compile(r"^(\d+)_3T_rfMRI_REST_fix\.zip$")

RUN_SPECS: Sequence[Tuple[str, str, float]] = (
    ("rest", "rfMRI_REST1_LR", 0.10),
    ("rest", "rfMRI_REST1_RL", 0.30),
    ("rest", "rfMRI_REST2_LR", 0.50),
    ("rest", "rfMRI_REST2_RL", 0.70),
    ("motor", "tfMRI_MOTOR_LR", 0.85),
    ("motor", "tfMRI_MOTOR_RL", 1.00),
)

ROI_NAMES = ("cortex_left", "cortex_right", "cerebellum", "subcortical")


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("extract_hcp")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    fh = logging.FileHandler(log_file, mode="w")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    return logger


def collect_subjects(folder: Path, pattern: re.Pattern) -> Dict[str, Path]:
    out: Dict[str, Path] = {}
    for path in sorted(folder.glob("*.zip")):
        m = pattern.match(path.name)
        if m:
            out[m.group(1)] = path
    return out


def detect_endian(prefix: bytes) -> str:
    le = struct.unpack("<i", prefix[0:4])[0]
    be = struct.unpack(">i", prefix[0:4])[0]
    if le == 540:
        return "<"
    if be == 540:
        return ">"
    raise ValueError("Not a valid NIfTI-2 header (sizeof_hdr != 540).")


def parse_cifti_header(prefix: bytes) -> Dict[str, object]:
    if len(prefix) < 552:
        raise ValueError("Prefix too short for NIfTI-2 + extension parsing.")
    endian = detect_endian(prefix)
    datatype = struct.unpack(f"{endian}h", prefix[12:14])[0]
    bitpix = struct.unpack(f"{endian}h", prefix[14:16])[0]
    dim = struct.unpack(f"{endian}8q", prefix[16:80])
    vox_offset = struct.unpack(f"{endian}q", prefix[168:176])[0]
    ext_flag = prefix[540:544]
    if ext_flag[0] != 1:
        raise ValueError("CIFTI XML extension flag not found in NIfTI header.")
    esize, ecode = struct.unpack(f"{endian}ii", prefix[544:552])
    if ecode != 32:
        raise ValueError(f"Unexpected extension code {ecode}; expected 32 for CIFTI.")
    xml_start = 552
    xml_end = xml_start + esize - 8
    if xml_end > len(prefix):
        raise ValueError("Prefix too short to capture full CIFTI XML extension.")
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
    raise FileNotFoundError(f"CIFTI dtseries member not found with suffix: {suffix}")


def read_member_prefix(zf: zipfile.ZipFile, member: str, nbytes: int = 900_000) -> bytes:
    with zf.open(member, "r") as handle:
        return handle.read(nbytes)


def parse_brain_models(xml_bytes: bytes) -> Tuple[int, float, List[Dict[str, object]]]:
    root = ET.fromstring(xml_bytes)
    series_map = None
    brain_map = None
    for mim in root.findall(".//MatrixIndicesMap"):
        t = mim.attrib.get("IndicesMapToDataType", "")
        if t == "CIFTI_INDEX_TYPE_SERIES":
            series_map = mim
        elif t == "CIFTI_INDEX_TYPE_BRAIN_MODELS":
            brain_map = mim
    if series_map is None or brain_map is None:
        raise ValueError("Missing SERIES or BRAIN_MODELS map in CIFTI XML.")

    n_series = int(series_map.attrib.get("NumberOfSeriesPoints", "0"))
    series_step = float(series_map.attrib.get("SeriesStep", "0.0"))
    models: List[Dict[str, object]] = []
    for bm in brain_map.findall("BrainModel"):
        offset = int(bm.attrib.get("IndexOffset", "0"))
        count = int(bm.attrib.get("IndexCount", "0"))
        structure = bm.attrib.get("BrainStructure", "")
        models.append(
            {
                "offset": offset,
                "count": count,
                "structure": structure,
            }
        )
    return n_series, series_step, models


def ranges_from_models(models: Sequence[Dict[str, object]], structures: Iterable[str]) -> List[int]:
    target = set(structures)
    out: List[int] = []
    for bm in models:
        structure = str(bm["structure"])
        if structure in target:
            start = int(bm["offset"])
            count = int(bm["count"])
            out.extend(range(start, start + count))
    return out


def build_roi_indices(models: Sequence[Dict[str, object]]) -> Dict[str, List[int]]:
    cortex_left = ranges_from_models(models, ["CIFTI_STRUCTURE_CORTEX_LEFT"])
    cortex_right = ranges_from_models(models, ["CIFTI_STRUCTURE_CORTEX_RIGHT"])
    cerebellum = ranges_from_models(
        models, ["CIFTI_STRUCTURE_CEREBELLUM_LEFT", "CIFTI_STRUCTURE_CEREBELLUM_RIGHT"]
    )
    excluded = {
        "CIFTI_STRUCTURE_CORTEX_LEFT",
        "CIFTI_STRUCTURE_CORTEX_RIGHT",
        "CIFTI_STRUCTURE_CEREBELLUM_LEFT",
        "CIFTI_STRUCTURE_CEREBELLUM_RIGHT",
    }
    subcortical = []
    for bm in models:
        structure = str(bm["structure"])
        if structure in excluded:
            continue
        start = int(bm["offset"])
        count = int(bm["count"])
        subcortical.extend(range(start, start + count))

    roi = {
        "cortex_left": sorted(cortex_left),
        "cortex_right": sorted(cortex_right),
        "cerebellum": sorted(cerebellum),
        "subcortical": sorted(subcortical),
    }
    for k in ROI_NAMES:
        if len(roi[k]) == 0:
            raise ValueError(f"ROI {k} has zero indices; cannot proceed.")
    return roi


def sample_roi_indices(
    roi_indices: Dict[str, List[int]], sample_per_roi: int, seed: int
) -> Dict[str, List[int]]:
    rng = random.Random(seed)
    out: Dict[str, List[int]] = {}
    for roi, idx in roi_indices.items():
        n_take = min(sample_per_roi, len(idx))
        out[roi] = sorted(rng.sample(idx, n_take))
    return out


def read_exact(handle, n: int) -> bytes:
    chunks = []
    remaining = n
    while remaining > 0:
        chunk = handle.read(remaining)
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    data = b"".join(chunks)
    if len(data) != n:
        raise EOFError(f"Expected {n} bytes, got {len(data)} bytes.")
    return data


def skip_exact(handle, n: int) -> None:
    remaining = n
    while remaining > 0:
        chunk = handle.read(min(4_194_304, remaining))
        if not chunk:
            raise EOFError(f"Unexpected EOF while skipping {n} bytes.")
        remaining -= len(chunk)


def timeseries_sd(ts_bytes: bytes, n_time: int) -> float:
    vals = array("f")
    vals.frombytes(ts_bytes)
    if sys.byteorder != "little":
        vals.byteswap()
    if len(vals) != n_time:
        raise ValueError(f"Expected {n_time} floats, got {len(vals)}.")
    mean = math.fsum(vals) / n_time
    var = math.fsum((v - mean) * (v - mean) for v in vals) / n_time
    return math.sqrt(var)


def quantile_linear(sorted_vals: Sequence[float], p: float) -> float:
    n = len(sorted_vals)
    if n == 0:
        raise ValueError("Cannot compute quantile of empty sequence.")
    if n == 1:
        return float(sorted_vals[0])
    pos = p * (n - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(sorted_vals[lo])
    w = pos - lo
    return float((1.0 - w) * sorted_vals[lo] + w * sorted_vals[hi])


def quantiles(values: Sequence[float], probs: Sequence[float]) -> List[float]:
    sv = sorted(values)
    return [quantile_linear(sv, p) for p in probs]


def process_member(
    zf: zipfile.ZipFile,
    member: str,
    selected_indices: Dict[str, List[int]],
    block_grayords: int,
) -> Tuple[Dict[str, List[float]], int, int]:
    prefix = read_member_prefix(zf, member, nbytes=900_000)
    hdr = parse_cifti_header(prefix)
    n_series, _series_step, models = parse_brain_models(hdr["xml_bytes"])
    dim = hdr["dim"]
    n_gray = int(dim[6]) if len(dim) > 6 else 0
    if n_gray <= 0 or n_series <= 0:
        raise ValueError("Invalid n_series/n_gray parsed from CIFTI header.")

    index_to_roi: Dict[int, List[str]] = {}
    for roi, idx_list in selected_indices.items():
        for idx in idx_list:
            index_to_roi.setdefault(idx, []).append(roi)
    selected_union = sorted(index_to_roi.keys())
    if not selected_union:
        raise ValueError("No selected grayordinates available.")

    ts_bytes_len = n_series * 4
    roi_values: Dict[str, List[float]] = {k: [] for k in ROI_NAMES}

    with zf.open(member, "r") as handle:
        skip_exact(handle, int(hdr["vox_offset"]))
        sel_ptr = 0
        n_sel = len(selected_union)

        for block_start in range(0, n_gray, block_grayords):
            gcount = min(block_grayords, n_gray - block_start)
            block = read_exact(handle, gcount * ts_bytes_len)
            block_end = block_start + gcount

            while sel_ptr < n_sel and selected_union[sel_ptr] < block_end:
                idx = selected_union[sel_ptr]
                if idx < block_start:
                    sel_ptr += 1
                    continue
                local = idx - block_start
                byte_lo = local * ts_bytes_len
                byte_hi = byte_lo + ts_bytes_len
                ts = block[byte_lo:byte_hi]
                sd_val = timeseries_sd(ts, n_series)
                for roi in index_to_roi[idx]:
                    roi_values[roi].append(sd_val)
                sel_ptr += 1

    # Basic consistency check with brain models
    if sum(int(m["count"]) for m in models) != n_gray:
        raise ValueError("Brain model index counts do not match n_gray dimension.")
    return roi_values, n_series, n_gray


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract ROI quantiles from HCP CIFTI dtseries.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent, help="Application data/output root (default: script directory).")
    parser.add_argument(
        "--max-subjects",
        type=int,
        default=20,
        help="Pilot subset size. Use <=0 to process all intersected subjects.",
    )
    parser.add_argument(
        "--sample-per-roi",
        type=int,
        default=1500,
        help="Number of grayordinates sampled per ROI in each run.",
    )
    parser.add_argument(
        "--quantile-grid-size",
        type=int,
        default=101,
        help="Number of quantile probabilities in [0,1].",
    )
    parser.add_argument(
        "--block-grayords",
        type=int,
        default=256,
        help="Grayordinate block size for streaming reads.",
    )
    parser.add_argument("--seed", type=int, default=20260304, help="Random seed.")
    parser.add_argument(
        "--output-quantiles",
        type=Path,
        default=Path("data_processed/hcp_roi_quantiles_raw.csv"),
        help="Output CSV for ROI quantiles.",
    )
    parser.add_argument(
        "--output-runs",
        type=Path,
        default=Path("data_processed/hcp_run_inventory.csv"),
        help="Output CSV for run-level inventory/status.",
    )
    parser.add_argument(
        "--output-sampling",
        type=Path,
        default=Path("data_processed/hcp_roi_sampling.json"),
        help="Output JSON for ROI sampling details.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/02_extract_roi_distributions.log"),
        help="Log file path.",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    logger = setup_logger(root / args.log_file)

    motor_dir = root / "HCP_Data" / "MotorTask"
    rest_dir = root / "HCP_Data" / "RestingState"
    if not motor_dir.is_dir() or not rest_dir.is_dir():
        raise FileNotFoundError("Supply HCP_Data/MotorTask and HCP_Data/RestingState under --root; see the application README.")
    motor = collect_subjects(motor_dir, MOTOR_RE)
    rest = collect_subjects(rest_dir, REST_RE)
    overlap = sorted(set(motor) & set(rest))

    if args.max_subjects is not None and args.max_subjects > 0:
        subjects = overlap[: args.max_subjects]
    else:
        subjects = overlap

    logger.info("Motor ZIPs: %d | Rest ZIPs: %d | Overlap: %d", len(motor), len(rest), len(overlap))
    logger.info("Subjects selected for extraction: %d", len(subjects))
    if not subjects:
        raise ValueError("No overlapping subjects available.")

    # Build ROI indices from a reference CIFTI file.
    ref_sid = subjects[0]
    with zipfile.ZipFile(rest[ref_sid], "r") as zf:
        ref_member = find_member(zf, run="rfMRI_REST1_LR", condition="rest")
        ref_prefix = read_member_prefix(zf, ref_member, nbytes=900_000)
    ref_hdr = parse_cifti_header(ref_prefix)
    _n_series_ref, _series_step_ref, ref_models = parse_brain_models(ref_hdr["xml_bytes"])
    roi_indices = build_roi_indices(ref_models)
    sampled = sample_roi_indices(roi_indices, sample_per_roi=args.sample_per_roi, seed=args.seed)

    sampling_meta = {
        "reference_subject": ref_sid,
        "reference_zip": str(rest[ref_sid]),
        "reference_member": ref_member,
        "sample_per_roi_requested": int(args.sample_per_roi),
        "seed": int(args.seed),
        "roi_total_counts": {k: len(v) for k, v in roi_indices.items()},
        "roi_sampled_counts": {k: len(v) for k, v in sampled.items()},
        "roi_sampled_indices_head": {k: v[:30] for k, v in sampled.items()},
    }
    (root / args.output_sampling).parent.mkdir(parents=True, exist_ok=True)
    with (root / args.output_sampling).open("w") as f:
        json.dump(sampling_meta, f, indent=2)
    logger.info("Saved ROI sampling metadata: %s", root / args.output_sampling)

    probs = [i / (args.quantile_grid_size - 1) for i in range(args.quantile_grid_size)]

    quant_rows: List[List[object]] = []
    run_rows: List[List[object]] = []

    for sid in subjects:
        logger.info("Processing subject %s", sid)
        for condition, run_name, t_val in RUN_SPECS:
            zip_path = rest[sid] if condition == "rest" else motor[sid]
            status = "ok"
            message = ""
            member = ""
            n_series = ""
            n_gray = ""
            try:
                with zipfile.ZipFile(zip_path, "r") as zf:
                    member = find_member(zf, run=run_name, condition=condition)
                    roi_vals, n_series_i, n_gray_i = process_member(
                        zf=zf,
                        member=member,
                        selected_indices=sampled,
                        block_grayords=args.block_grayords,
                    )
                n_series = n_series_i
                n_gray = n_gray_i

                for roi in ROI_NAMES:
                    vals = roi_vals[roi]
                    qvals = quantiles(vals, probs)
                    for p, q in zip(probs, qvals):
                        quant_rows.append(
                            [
                                sid,
                                condition,
                                run_name,
                                t_val,
                                roi,
                                p,
                                q,
                                len(vals),
                                n_series_i,
                                n_gray_i,
                                str(zip_path),
                                member,
                            ]
                        )
            except Exception as exc:  # pylint: disable=broad-except
                status = "error"
                message = str(exc)
                logger.exception(
                    "Failed subject=%s condition=%s run=%s: %s", sid, condition, run_name, exc
                )

            run_rows.append(
                [
                    sid,
                    condition,
                    run_name,
                    t_val,
                    status,
                    message,
                    str(zip_path),
                    member,
                    n_series,
                    n_gray,
                ]
            )

    out_quant = root / args.output_quantiles
    out_runs = root / args.output_runs
    out_quant.parent.mkdir(parents=True, exist_ok=True)
    out_runs.parent.mkdir(parents=True, exist_ok=True)

    with out_quant.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "subject_id",
                "condition",
                "run",
                "t",
                "roi",
                "p",
                "q_raw",
                "n_sampled_grayords",
                "n_timepoints",
                "n_total_grayords",
                "source_zip",
                "source_member",
            ]
        )
        w.writerows(quant_rows)

    with out_runs.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "subject_id",
                "condition",
                "run",
                "t",
                "status",
                "message",
                "source_zip",
                "source_member",
                "n_timepoints",
                "n_total_grayords",
            ]
        )
        w.writerows(run_rows)

    n_ok = sum(1 for r in run_rows if r[4] == "ok")
    n_total = len(run_rows)
    logger.info("Run extraction status: %d / %d succeeded", n_ok, n_total)
    logger.info("Wrote quantile file: %s", out_quant)
    logger.info("Wrote run inventory: %s", out_runs)
    logger.info("Extraction step completed.")


if __name__ == "__main__":
    main()
