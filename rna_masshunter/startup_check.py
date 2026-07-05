import importlib.util
import sys
from pathlib import Path
from typing import Any

from rna_masshunter.models import RunConfig
from rna_masshunter.warnings_manager import add_warning


REQUIRED_PACKAGES = ["yaml", "pandas", "openpyxl", "numpy", "pyteomics", "lxml", "tqdm"]


def run_startup_check(project_root: str | Path, config: RunConfig, logger, warnings: list[dict[str, Any]]) -> None:
    root = Path(project_root)
    logger.info("Startup check: Python %s", sys.version.split()[0])
    if sys.version_info < (3, 10):
        add_warning(warnings, "WARNING", "startup_check", "Python 3.10 or newer is recommended.")

    for package in REQUIRED_PACKAGES:
        if importlib.util.find_spec(package) is None:
            add_warning(warnings, "ERROR", "startup_check", f"Required package is not installed: {package}")
            logger.error("Required package is not installed: %s", package)
        else:
            logger.info("Package OK: %s", package)

    required_paths = [
        root / "config.yaml",
        root / "data" / "modifications.yaml",
        root / "data" / "rule_sets",
        root / "data" / "pathways",
    ]
    for path in required_paths:
        if path.exists():
            logger.info("Path OK: %s", path)
        else:
            add_warning(warnings, "ERROR", "startup_check", f"Required path is missing: {path}")
            logger.error("Required path is missing: %s", path)

    for key in ("output_dir", "log_dir", "cache_dir"):
        path = Path(config.project[key])
        try:
            path.mkdir(parents=True, exist_ok=True)
            logger.info("Writable directory OK: %s", path)
        except OSError as exc:
            add_warning(warnings, "ERROR", "startup_check", f"Cannot create directory: {path}", str(exc))

    mzml_path = config.input.get("mzml_path")
    raw_path = config.input.get("raw_path")
    if not mzml_path and not raw_path:
        add_warning(warnings, "WARNING", "startup_check", "No mzML or raw input configured. Edit config.yaml before real analysis.")
    else:
        for label, value in (("mzML", mzml_path), ("raw", raw_path)):
            if value and not Path(value).exists():
                add_warning(warnings, "ERROR", "startup_check", f"{label} input path does not exist.", value)

    raw_suffix = Path(raw_path).suffix.lower() if raw_path else ""
    if raw_suffix in {".wiff", ".wiff2"} and not config.input.get("msconvert_path"):
        add_warning(warnings, "ERROR", "startup_check", "msconvert_path is required for WIFF/WIFF2 conversion.")
