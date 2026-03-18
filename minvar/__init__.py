"""
MinVar: Minority Variant Detection Pipeline

Author: Germán Vallejo Palma
Institution: Instituto de Salud Carlos III (ISCIII)

This package provides tools for detecting and filtering minority variants
from IRMA-generated viral sequencing data.
"""

from .main import run_minvar

__all__ = ["run_minvar"]
__version__ = "0.1.0"
