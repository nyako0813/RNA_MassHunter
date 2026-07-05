import json
from pathlib import Path
from typing import Any


def _checkpoint_path(checkpoint_dir: str | Path, name: str) -> Path:
    return Path(checkpoint_dir) / f"{name}.json"


def save_checkpoint(checkpoint_dir: str | Path, name: str, data: dict[str, Any]) -> Path:
    directory = Path(checkpoint_dir)
    directory.mkdir(parents=True, exist_ok=True)
    path = _checkpoint_path(directory, name)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)
    return path


def load_checkpoint(checkpoint_dir: str | Path, name: str) -> dict[str, Any] | None:
    path = _checkpoint_path(checkpoint_dir, name)
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def checkpoint_exists(checkpoint_dir: str | Path, name: str) -> bool:
    return _checkpoint_path(checkpoint_dir, name).exists()
