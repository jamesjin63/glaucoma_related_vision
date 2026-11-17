"""Shared paths for the command-line analysis modules."""
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
INPUT = Path(os.environ.get('GLAUCOMA_INPUT', ROOT / 'data/local')).expanduser().resolve()
OUTPUT = Path(os.environ.get('GLAUCOMA_OUTPUT', ROOT / 'results')).expanduser().resolve()
if INPUT == OUTPUT or INPUT in OUTPUT.parents or OUTPUT in INPUT.parents:
    raise ValueError('Input and output directories must be separate and must not contain each other.')
