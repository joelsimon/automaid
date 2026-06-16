"""Command line interface for automaid."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .preflight import PreflightError, scripts_dir, validate_process_inputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="automaid",
        description="Process MERMAID files with the legacy automaid workflow.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"automaid {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command")

    process = subparsers.add_parser(
        "process",
        help="process MERMAID server files",
        description="Process MERMAID files from a server directory into processed outputs.",
    )
    process.add_argument(
        "-s",
        "--server",
        required=True,
        help="directory containing raw MERMAID server files",
    )
    process.add_argument(
        "-p",
        "--processed",
        required=True,
        help="directory for processed automaid outputs",
    )
    process.add_argument(
        "-d",
        "--database",
        required=True,
        help="directory for downloaded or cached MERMAID database files",
    )
    process.set_defaults(func=_process)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not hasattr(args, "func"):
        parser.print_help()
        return 0

    try:
        return args.func(args)
    except PreflightError as exc:
        print("automaid preflight failed:", file=sys.stderr)
        for message in exc.messages:
            print(f"- {message}", file=sys.stderr)
        return 2


def _process(args: argparse.Namespace) -> int:
    validate_process_inputs(
        server=args.server,
        processed=args.processed,
        database=args.database,
    )

    legacy_scripts = str(scripts_dir())
    if legacy_scripts not in sys.path:
        sys.path.insert(0, legacy_scripts)

    import main as legacy_main

    legacy_main.main(
        server=str(Path(args.server).expanduser()),
        processed=str(Path(args.processed).expanduser()),
        database=str(Path(args.database).expanduser()),
    )
    return 0
