# Get the version number from the setup file
import os

DIRECTORY = os.path.dirname(os.path.abspath(__file__))
ABOUT_DIR = os.path.join(DIRECTORY, "about")
__version__ = open(os.path.join(ABOUT_DIR,"version.txt")).read().strip()

from .balltree import *
