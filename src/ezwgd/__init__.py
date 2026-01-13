#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

""" """
from rich.console import Console
from rich.traceback import install

from ezwgd._base import version, nickname
from ezwgd import evo
from ezwgd import tidy
from ezwgd import utils

# Start Rich Engine.
console = Console()
install(show_locals=True)

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
