#!/usr/bin/env python3
"""Run one Namibia model scenario with portable data and output paths."""

from __future__ import annotations

import argparse
import importlib
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODELS = {
    "full": "namibia_pde_model_total_death",
    "nofood": "namibia_pde_model_total_death_nofood",
    "notoxin": "namibia_pde_model_total_death_notoxin",
    "uniform_low_food": "namibia_pde_model_uniform_low_food",
}


def required_file(parser: argparse.ArgumentParser, data_dir: Path, name: str, value: str) -> None:
    """Check the input before importing heavy simulation dependencies."""
    file_path = Path(value).expanduser()
    if not file_path.is_absolute():
        file_path = data_dir / file_path
    if not file_path.is_file():
        parser.error(f"{name} file not found: {file_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", choices=MODELS.keys())
    parser.add_argument(
        "--data-dir", type=Path,
        default=Path(os.environ.get("NAMIBIA_DATA_DIR", ROOT / "data" / "namibia")),
        help="directory with the RData grid, mortality CSV, and mask raster",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        help="output folder for this scenario; defaults to RESULTS_DIR/scenario",
    )
    parser.add_argument("--grid", default="ud_pois_plot.rdata", help="RData file name or path")
    parser.add_argument("--mortality", default="Poisoning_GIS.csv", help="CSV file name or path")
    parser.add_argument("--mask", default="combinedgeo.tif", help="Namibia raster file name or path")
    parser.add_argument(
        "--parks", default="Namibia_National_Parks_Boundaries_2001.shp",
        help="optional park boundaries shapefile for figures",
    )
    parser.add_argument("--n-steps", type=int, help="override the default 3600 integration steps")
    parser.add_argument("--no-plots", action="store_true", help="save CSV/JSON outputs without figures")
    parser.add_argument("--quiet", action="store_true", help="hide the full JSON summary")
    parser.add_argument(
        "--uniform-food-level", type=float,
        help="override the uniform_low_food scenario's calibrated default of 0.025737",
    )
    args = parser.parse_args()

    data_dir = args.data_dir.expanduser().resolve()
    results_dir = Path(os.environ.get("NAMIBIA_RESULTS_DIR", ROOT / "results" / "namibia"))
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir else (results_dir / args.scenario).expanduser().resolve()
    )
    if args.n_steps is not None and args.n_steps < 1:
        parser.error("--n-steps must be positive")
    if args.uniform_food_level is not None:
        if args.scenario != "uniform_low_food":
            parser.error("--uniform-food-level applies only to uniform_low_food")
        if not 0 <= args.uniform_food_level <= 1:
            parser.error("--uniform-food-level must be between 0 and 1")
    for name, value in (("grid", args.grid), ("mortality", args.mortality), ("mask", args.mask)):
        required_file(parser, data_dir, name, value)

    model = importlib.import_module(MODELS[args.scenario])
    params = model.Params()
    params.base_dir = str(data_dir)
    params.grid_rdata_name = args.grid
    params.mortality_csv_name = args.mortality
    params.namibia_mask_tif_name = args.mask
    params.parks_shp_name = args.parks
    # os.path.join(data_dir, absolute_output_dir) preserves the absolute output.
    params.output_dir_name = str(output_dir)
    params.make_plots = not args.no_plots
    params.verbose = not args.quiet
    if args.n_steps is not None:
        params.n_steps = args.n_steps
    if args.uniform_food_level is not None:
        params.uniform_food_level = args.uniform_food_level
    model.main(params)
    print(f"Saved {args.scenario} results to {output_dir}")


if __name__ == "__main__":
    main()
