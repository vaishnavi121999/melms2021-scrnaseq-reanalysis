"""Shared utilities for the Melms 2021 reanalysis pipeline."""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict

import pandas as pd
import yaml

LOGGER_NAME = "melms2021"


def setup_logger() -> logging.Logger:
    """Configure a module-level logger only once."""
    logger = logging.getLogger(LOGGER_NAME)
    if logger.handlers:
        return logger
    handler = logging.StreamHandler()
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


def load_config(path: str | Path = "workflow/config.yaml") -> Dict[str, Any]:
    """Load YAML config and keep absolute root for downstream modules."""
    path = Path(path)
    with path.open() as handle:
        config = yaml.safe_load(handle)
    config["project_root"] = str(path.parent.parent.resolve())
    return config


def ensure_parent(path: str | Path) -> None:
    """Create parent directories for the provided path if they do not exist."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_sample_sheet(path: str | Path) -> pd.DataFrame:
    """Read the metadata/sample sheet and enforce required columns."""
    df = pd.read_csv(path)
    required = {"sample_id", "file_name", "condition"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Sample sheet missing columns: {', '.join(sorted(missing))}")
    return df


def write_json(data: Dict[str, Any], path: str | Path) -> None:
    ensure_parent(path)
    with Path(path).open("w") as handle:
        json.dump(data, handle, indent=2)


def resolve_path(path: str | Path, root: str | Path | None = None) -> Path:
    """Resolve relative paths against the provided project root."""
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    base = Path(root) if root else Path.cwd()
    return (base / candidate).resolve()
