#!/usr/bin/env python3
"""Render the migrated, source-specific mois four-panel climate figure.

Requires numpy/matplotlib; preserves the original FCP rendering implementation.
Original figures/global_moisture_colour.svg remains the immutable historical copy.
"""
from __future__ import annotations
import argparse
import importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
LEGACY=ROOT/"scripts/analysis/render_fcp_original_macro_figure_20261014.py"
spec=importlib.util.spec_from_file_location("mois_original_figure",LEGACY)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",default=str(ROOT/"results/figure_source_data.json"))
    parser.add_argument("--outdir",default=str(ROOT/"figures/recreated"))
    args=parser.parse_args()
    created=mod.make_figure(args.source,args.outdir)
    for path in created:
        print(path)

if __name__=="__main__":
    main()
