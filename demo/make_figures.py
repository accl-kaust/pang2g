#!/usr/bin/env python3
"""
demo/make_figures.py — compatibility shim.

The paper's figure generator now lives at the repository root
(../make_figures.py) and reads the committed measurement tables by default.
This shim forwards all arguments to it so the `python3 demo/make_figures.py
--<x>-csv ...` hints printed by the demo harnesses keep working.
"""
import os
import runpy
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.argv[0] = os.path.join(ROOT, "make_figures.py")
runpy.run_path(sys.argv[0], run_name="__main__")
