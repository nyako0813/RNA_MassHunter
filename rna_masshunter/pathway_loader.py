from pathlib import Path
from typing import Any

import yaml

from rna_masshunter.warnings_manager import add_warning


def load_pathways(pathway_dir: str | Path, warnings: list[dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    directory = Path(pathway_dir)
    pathways = []
    if not directory.exists():
        if warnings is not None:
            add_warning(warnings, "ERROR", "pathway_loader", "Pathway directory is missing.", str(directory))
        return pathways
    for path in sorted(directory.glob("*.yaml")):
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
        data.setdefault("_file", path.name)
        pathways.append(data)
    return pathways


def validate_pathways(pathways: list[dict[str, Any]], warnings: list[dict[str, Any]] | None = None) -> None:
    if not pathways and warnings is not None:
        add_warning(warnings, "WARNING", "pathway_loader", "No pathway YAML files were loaded.")
    for pathway in pathways:
        if not (pathway.get("id") or pathway.get("name")) and warnings is not None:
            add_warning(warnings, "WARNING", "pathway_loader", "Pathway lacks id/name.", pathway.get("_file"))
