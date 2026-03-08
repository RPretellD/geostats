from importlib.metadata import version as _version
__version__ = _version("geostats")

from .Kriging import Model
from .Kriging import Site
from .Kriging import Kriging

import geostats.build_Mrho as build_Mrho
import geostats.geostats_tools as geostats_tools
