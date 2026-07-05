from pathlib import Path
from subprocess import run
from typing import Any

from rna_masshunter.models import RunConfig
from rna_masshunter.warnings_manager import add_warning


def prepare_input_file(config: RunConfig, logger, warnings: list[dict[str, Any]]) -> Path | None:
    mzml_path = config.input.get("mzml_path")
    raw_path = config.input.get("raw_path")
    if mzml_path:
        path = Path(mzml_path)
        if path.exists():
            logger.info("Using mzML input: %s", path)
            return path
        raise FileNotFoundError(f"Configured mzML file does not exist: {path}")
    if not raw_path:
        add_warning(warnings, "WARNING", "conversion", "No input file configured.")
        return None
    path = Path(raw_path)
    if not path.exists():
        raise FileNotFoundError(f"Configured raw file does not exist: {path}")
    if path.suffix.lower() == ".mzml":
        return path
    if path.suffix.lower() in {".wiff", ".wiff2"}:
        return run_msconvert(path, config.input.get("msconvert_path"), logger)
    raise ValueError(f"Unsupported input file extension for MVP-1: {path.suffix}")


def run_msconvert(raw_path: str | Path, msconvert_path: str | None, logger) -> Path:
    if not msconvert_path:
        raise ValueError("msconvert_path is required to convert WIFF/WIFF2 files in MVP-1.")
    raw = Path(raw_path)
    out_dir = raw.parent
    command = [msconvert_path, str(raw), "--mzML", "-o", str(out_dir)]
    logger.info("Running msconvert for %s", raw)
    result = run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"msconvert failed: {result.stderr or result.stdout}")
    mzml_path = out_dir / f"{raw.stem}.mzML"
    if not mzml_path.exists():
        raise FileNotFoundError(f"msconvert completed but mzML was not found: {mzml_path}")
    return mzml_path
