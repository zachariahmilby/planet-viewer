from pathlib import Path

_project_directory = Path(__file__).resolve().parent
"""Absolute path to this package's installation directory."""

kernel_path = Path(_project_directory, 'anc', 'kernels')
"""Package SPICE kernel path."""
