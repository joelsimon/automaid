"""Thin package wrapper for the legacy automaid processing scripts."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _legacy_version() -> str:
    setup_path = Path(__file__).resolve().parents[1] / "scripts" / "setup.py"
    if not setup_path.exists():
        return "v4.5.3"

    spec = importlib.util.spec_from_file_location("_automaid_legacy_setup", setup_path)
    if spec is None or spec.loader is None:
        return "v4.5.3"

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.get_version()


__version__ = _legacy_version()
