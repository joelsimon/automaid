"""Preflight validation for the legacy automaid processing workflow."""

from __future__ import annotations

import os
from pathlib import Path


REQUIRED_ICDF_BINARIES = (
    "icdf24_v103_test",
    "icdf24_v103ec_test",
)


class PreflightError(RuntimeError):
    """Raised when automaid cannot safely start processing."""

    def __init__(self, messages: list[str]):
        self.messages = messages
        super().__init__("\n".join(messages))


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def scripts_dir() -> Path:
    return repo_root() / "scripts"


def default_bin_dir() -> Path:
    return scripts_dir() / "bin"


def validate_process_inputs(
    *,
    server: str | os.PathLike[str],
    processed: str | os.PathLike[str],
    database: str | os.PathLike[str],
    bin_dir: str | os.PathLike[str] | None = None,
) -> None:
    """Validate paths and native binaries needed by ``automaid process``."""

    errors: list[str] = []
    server_path = Path(server).expanduser()
    processed_path = Path(processed).expanduser()
    database_path = Path(database).expanduser()
    icdf_bin_dir = Path(bin_dir).expanduser() if bin_dir else default_bin_dir()

    _check_existing_directory(errors, "server", server_path)
    _check_creatable_directory(errors, "processed", processed_path)
    _check_creatable_directory(errors, "database", database_path)
    _check_icdf_binaries(errors, icdf_bin_dir)

    if errors:
        raise PreflightError(errors)


def _check_existing_directory(errors: list[str], label: str, path: Path) -> None:
    if not path.exists():
        errors.append(
            f"{label} path does not exist: {path}\n"
            f"  Fix: create it or pass --{label} PATH pointing to an existing directory."
        )
        return

    if not path.is_dir():
        errors.append(
            f"{label} path is not a directory: {path}\n"
            f"  Fix: pass --{label} PATH pointing to a directory."
        )


def _check_creatable_directory(errors: list[str], label: str, path: Path) -> None:
    if path.exists():
        if not path.is_dir():
            errors.append(
                f"{label} path exists but is not a directory: {path}\n"
                f"  Fix: remove the file or pass --{label} PATH pointing to a directory."
            )
        elif not os.access(path, os.W_OK):
            errors.append(
                f"{label} directory is not writable: {path}\n"
                f"  Fix: choose a writable directory or update permissions."
            )
        return

    parent = path.parent if path.parent != Path("") else Path(".")
    if not parent.exists():
        errors.append(
            f"{label} directory cannot be created because its parent is missing: {path}\n"
            f"  Fix: create {parent} first or pass --{label} PATH under an existing directory."
        )
        return

    if not parent.is_dir():
        errors.append(
            f"{label} directory cannot be created because its parent is not a directory: {path}\n"
            f"  Fix: pass --{label} PATH under an existing directory."
        )
        return

    if not os.access(parent, os.W_OK):
        errors.append(
            f"{label} directory cannot be created because parent is not writable: {path}\n"
            f"  Fix: choose a writable parent directory or update permissions."
        )


def _check_icdf_binaries(errors: list[str], bin_dir: Path) -> None:
    if not bin_dir.exists():
        errors.append(
            f"icdf binary directory does not exist: {bin_dir}\n"
            "  Fix: compile the native wavelet inversion programs with:\n"
            f"    cd {scripts_dir() / 'src' / 'V103'} && make && cp icdf24_v103_test ../../bin/\n"
            f"    cd {scripts_dir() / 'src' / 'V103EC'} && make && cp icdf24_v103ec_test ../../bin/"
        )
        return

    if not bin_dir.is_dir():
        errors.append(
            f"icdf binary path is not a directory: {bin_dir}\n"
            "  Fix: pass a directory containing icdf24_v103_test and icdf24_v103ec_test."
        )
        return

    for binary in REQUIRED_ICDF_BINARIES:
        binary_path = bin_dir / binary
        if not binary_path.exists():
            errors.append(
                f"required icdf binary is missing: {binary_path}\n"
                "  Fix: compile and copy the native binary:\n"
                f"    cd {scripts_dir() / 'src' / _source_dir_for(binary)} && make && cp {binary} ../../bin/"
            )
        elif not os.access(binary_path, os.X_OK):
            errors.append(
                f"required icdf binary is not executable: {binary_path}\n"
                f"  Fix: run chmod +x {binary_path} or rebuild it with make."
            )


def _source_dir_for(binary: str) -> str:
    if binary == "icdf24_v103ec_test":
        return "V103EC"
    return "V103"
