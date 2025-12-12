# scripts/silence_warnings.py
"""
Silence known harmless warnings from Python and third-party libraries.
Import this at the top of every Snakemake script to keep console output clean.
"""

import warnings
import os

# ---  General Python warnings ---
warnings.filterwarnings("ignore", category=SyntaxWarning)

# ---  pkg_resources deprecation (from louvain / setuptools) ---
warnings.filterwarnings("ignore", message=".*pkg_resources.*deprecated.*")

# ---  Scanpy sparse densification ---
warnings.filterwarnings("ignore", message=".*zero-centering a sparse array.*")

# ---  Prevent propagation to stderr in spawned processes (optional safety) ---
os.environ["PYTHONWARNINGS"] = "ignore"
