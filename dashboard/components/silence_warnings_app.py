# components/silence_warnings_app.py
"""
Silence known harmless warnings from Python and third-party libraries.
Import this at the top of every Snakemake script to keep console output clean.
"""

import warnings
import os

# --- General Python warnings ---
warnings.filterwarnings("ignore", category=SyntaxWarning)

# --- pkg_resources deprecation (from louvain / setuptools) ---
warnings.filterwarnings("ignore", message=".*pkg_resources.*deprecated.*")

# --- Scanpy sparse densification ---
warnings.filterwarnings("ignore", message=".*zero-centering a sparse array.*")

# --- pandas FutureWarning about observed=False in groupby ---
warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    message=".*The default of observed=False is deprecated.*"
)

# --- MuData FutureWarnings about .update() ---
warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    message=".*From 0\\.4 \.update\\(\\) will not pull obs/var columns.*"
)

# --- anndata / scanpy __version__ FutureWarning ---
warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    message=".*`__version__` is deprecated.*"
)

# --- Prevent propagation to stderr in spawned processes ---
os.environ["PYTHONWARNINGS"] = "ignore"
