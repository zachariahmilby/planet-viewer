from pathlib import Path

"""Absolute path to this package's installation directory."""
_project_directory = Path(__file__).resolve().parent

"""Package SPICE kernel path."""
kernel_path = Path(_project_directory, 'anc', 'kernels')
