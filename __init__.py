"""
Open Ballistic Library - Package Initialization

This package provides tools for ballistic calculations, trajectory predictions,
and orbital mechanics for projectiles and space objects.
"""

from .main import BallisticCalculator
from .orbital_prediction import OrbitalPredictor

__version__ = "0.1.0"
__author__ = "Open Ballistic Library Contributors"
__license__ = "MIT"

def get_library_info():
    """Return basic information about the library."""
    return {
        "name": "Open Ballistic Library",
        "version": __version__,
        "author": __author__,
        "description": "An open-source library for ballistic calculations and orbital mechanics"
    }