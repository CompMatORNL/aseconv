"""  Default plugin modules

    Loaded by defaults file names with 'plug*.py' or 'io*.py'
    The custom plugin file path need to be set in 'ASEC_PLUGIN_PATH' enviornment variable.
"""
import os
import glob
from pathlib import Path
import traceback
from importlib import util
import importlib
# https://packaging.python.org/en/latest/guides/creating-and-discovering-plugins/

import sys, imp

def _load_module(path, name):
    module_name = os.path.splitext(os.path.basename(path))[0]
    spec = util.spec_from_file_location(name, path)
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _search_modules(dirpath):
    files = sorted(glob.glob(os.path.join(dirpath, "*.py")))
    for full in files:
        pfull = Path(full)
        fname = pfull.name
        fpath = pfull.parent
        if fname.startswith("plug") or fname.startswith("io"):
            try:
                _load_module(pfull, fname)
            except Exception:
                traceback.print_exc()


# Get current path
path = os.path.abspath(__file__)
dirpath = os.path.dirname(path)
_search_modules(dirpath)

custom_path = os.environ.get("ASEC_PLUGIN_PATH", "")
if custom_path != "":
    for x in custom_path.split(":"):
        _search_modules(x)
# additional plugin paths
