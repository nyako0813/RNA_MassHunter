import hashlib
import pickle
from pathlib import Path
from typing import Any


def get_cache_key(*parts: Any) -> str:
    digest = hashlib.sha256()
    for part in parts:
        digest.update(repr(part).encode("utf-8"))
    return digest.hexdigest()


def load_cache(cache_dir: str | Path, key: str) -> Any:
    path = Path(cache_dir) / f"{key}.pkl"
    if not path.exists():
        return None
    with path.open("rb") as handle:
        return pickle.load(handle)


def save_cache(cache_dir: str | Path, key: str, value: Any) -> Path:
    directory = Path(cache_dir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{key}.pkl"
    with path.open("wb") as handle:
        pickle.dump(value, handle)
    return path
