#!/usr/bin/env python3
"""Build a source-only public ZIP without modifying the working repository."""

from __future__ import annotations

import argparse
import hashlib
import io
import re
import sys
import zipfile
from pathlib import Path


EXCLUDED_DIRECTORIES = {
    "results", "result", "output", "outputs", "plots", "logs", "archive",
    "archives", "audit", "audits", "research_archive", "__pycache__",
    ".git", ".svn", ".cache", ".pytest_cache", ".ipynb_checkpoints",
}
SOURCE_SUFFIXES = {".r", ".py", ".sbatch"}
ASSET_SUFFIXES = {".pdf", ".png", ".tex"}
DATA_SUFFIXES = {
    ".rds", ".rdata", ".rda", ".csv", ".tsv", ".mat", ".h5", ".hdf5",
    ".nii", ".npy", ".npz", ".pkl", ".pickle", ".parquet", ".feather",
    ".xlsx", ".xls", ".sqlite", ".db", ".zip", ".gz", ".bz2", ".xz",
}
ROOT_DOCUMENTS = {
    "readme.md", "dependencies.tsv", "check_environment.r",
    "reproduce.md", "verification.md",
}
REQUIRED_FILES = {
    "dependencies.tsv", "check_environment.R", "REPRODUCE.md",
    "VERIFICATION.md", "tools/build_release.py",
}

# Generic patterns avoid retaining any author's machine or cluster names.
# Scheduler-specific allocations belong in user-supplied submission options.
PRIVATE_PATTERNS = (
    re.compile(r"/(?:Users|home|nas|Volumes)/[^\s\"'<>]+"),
    re.compile(r"^\s*#SBATCH\s+(?:--(?:account|partition|nodelist)(?:[=\s]|$)|-[Awp](?:[=\s]|[A-Za-z0-9]))"),
)


def is_readme(path: Path) -> bool:
    return path.name.lower() in {"readme", "readme.md"}


def collect_files(source: Path) -> tuple[dict[str, bytes], list[str]]:
    """Apply an explicit allowlist, then reject unsafe included content."""
    files: dict[str, bytes] = {}
    failures: list[str] = []
    for path in sorted(source.rglob("*")):
        relative = path.relative_to(source)
        parts_lower = [part.lower() for part in relative.parts]
        if any(part in EXCLUDED_DIRECTORIES for part in parts_lower[:-1]):
            continue
        if any(part.startswith(".") for part in relative.parts):
            continue
        if path.is_dir():
            continue
        if path.is_symlink():
            failures.append(f"Symlinks are not permitted in the release: {relative.as_posix()}")
            continue
        if not path.is_file():
            continue

        name = relative.as_posix()
        data_directory = any(part in {"data", "hcp_data", "data_processed"}
                             for part in parts_lower[:-1])
        if data_directory and not is_readme(path):
            failures.append(f"Data directory contains a non-README file: {name}")
            continue
        if path.suffix.lower() in DATA_SUFFIXES and name != "dependencies.tsv":
            failures.append(f"Data or archive file is not permitted: {name}")
            continue

        root_document = len(relative.parts) == 1 and path.name.lower() in ROOT_DOCUMENTS
        source_file = path.suffix.lower() in SOURCE_SUFFIXES
        asset = (len(relative.parts) == 2 and relative.parts[0] == "manuscript_assets"
                 and path.suffix.lower() in ASSET_SUFFIXES)
        if not (root_document or source_file or is_readme(path) or asset):
            continue
        content = path.read_bytes()
        if not (asset and path.suffix.lower() in {".pdf", ".png"}):
            try:
                text = content.decode("utf-8")
            except UnicodeDecodeError:
                failures.append(f"Distributed text must be UTF-8: {name}")
                continue
            for line_number, line in enumerate(text.splitlines(), 1):
                if any(pattern.search(line) for pattern in PRIVATE_PATTERNS):
                    failures.append(f"Personal path or machine-specific scheduler directive at {name}:{line_number}")
        files[name] = content

    for required in sorted(REQUIRED_FILES - files.keys()):
        failures.append(f"Missing required release file: {required}")
    if not any(Path(name).name.lower() == "readme.md" and "/" not in name for name in files):
        failures.append("Missing repository-root README.md or README.MD")
    for folder in ("Simulation", "RealData"):
        if not any(name.startswith(folder + "/") for name in files):
            failures.append(f"Missing distributed source directory: {folder}")
    return files, failures


def make_archive(files: dict[str, bytes]) -> bytes:
    checksum_lines = [f"{hashlib.sha256(content).hexdigest()}  {name}"
                      for name, content in sorted(files.items())]
    manifest = ("\n".join(checksum_lines) + "\n").encode("utf-8")
    payload = dict(files)
    # SHA256SUMS covers every other archive member; it cannot hash itself.
    payload["SHA256SUMS"] = manifest
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=9) as archive:
        for name, content in sorted(payload.items()):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, content, compresslevel=9)
    return buffer.getvalue()


def main() -> int:
    default_source = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=default_source,
                        help="Distributed code root containing Simulation/ and RealData/ (default: parent of tools/).")
    parser.add_argument("--output", type=Path,
                        help="Destination ZIP outside the code root (default: ../output/reproducibility/FPCADVP.zip).")
    parser.add_argument("--replace", action="store_true",
                        help="Explicitly allow replacing an existing destination ZIP.")
    args = parser.parse_args()
    source = args.source.expanduser().resolve()
    output = (args.output.expanduser().resolve() if args.output else
              source.parent / "output" / "reproducibility" / "FPCADVP.zip")
    if not source.is_dir():
        parser.error("--source must name an existing code directory")
    if output == source or source in output.parents:
        parser.error("--output must be outside the source directory")
    if output.suffix.lower() != ".zip":
        parser.error("--output must end in .zip")
    if output.exists() and not args.replace:
        parser.error("Destination exists; choose another output or pass --replace explicitly")
    files, failures = collect_files(source)
    if failures:
        print("Release validation failed:", file=sys.stderr)
        for failure in failures:
            print("- " + failure, file=sys.stderr)
        return 1
    archive_bytes = make_archive(files)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation prevents a race from overwriting an existing package.
    with output.open("wb" if args.replace else "xb") as handle:
        handle.write(archive_bytes)
    print(f"Created {output}")
    print(f"Included {len(files)} source/document/asset files plus SHA256SUMS; no data files.")
    print(f"ZIP SHA256: {hashlib.sha256(archive_bytes).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
