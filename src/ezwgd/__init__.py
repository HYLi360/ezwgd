#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

""" """

from ._base import version, nickname
from . import evo
from . import tidy
from . import utils


__version__ = f"{version} {nickname}"
__author__ = "HYLi360"

__all__ = [
    "__version__",
    "__author__",
    "version",
    "evo",
    "tidy",
    "utils"
]
