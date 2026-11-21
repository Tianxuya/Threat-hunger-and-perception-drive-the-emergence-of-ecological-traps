#!/usr/bin/env python3
"""Optuna calibration of the cognition model with fixed mixed initial condition w=0.8.

Initial condition:
    u0 = (1 - w) * UD + w * Uniform
with:
    w = 0.8

Therefore:
    u0 = 0.2 * UD + 0.8 * Uniform

Objective:
    maximize Pearson + 0.2*Top20 hotspot overlap - 0.5*MSE

Example:
    python scripts/python/calibrate_params.py --n-trials 10
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]


def load_model(path: str):
    path = os.path.abspath(path)
    spec = importlib.util.spec_from_file_location("namibia_total_death_model", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import model file: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_params(model, args, trial=None):
    p = model.Params()
    p.base_dir = args.base_dir
    p.grid_rdata_name = args.grid
    p.mortality_csv_name = args.mortality
    p.namibia_mask_tif_name = args.mask
    p.parks_shp_name = args.parks
    p.save_outputs = False
    p.make_plots = False
    p.verbose = False

    # Fixed mixed initial condition:
    # u0 = (1-w)*UD + w*Uniform, with w=0.8
    p.initial_condition = "mixed"
    p.uniform_weight = 0.8

    if hasattr(p, "condition"):
        p.condition = "Cognition"

    if trial is not None:
        p.search_radius_km_cognition = trial.suggest_float(
            "search_radius_km", 20.0, 80.0
        )
        p.sigma_factor = trial.suggest_float(
            "sigma_factor", 0.25, 1.0
        )
        p.diffusion_rate = trial.suggest_float(
            "diffusion_rate", 1e-4, 3e-2, log=True
        )
        p.directed_rate = trial.suggest_float(
            "directed_rate", 1e-2, 1.5, log=True
        )
        p.d1 = trial.suggest_float(
            "d1",
            2.7392206258276093e-6,
            4.3851814415556833e-4,
            log=True,
        )
        p.decision_epsilon = trial.suggest_float(
            "decision_epsilon", 0.0, 0.10
        )
        p.hunger_threshold = trial.suggest_float(
            "hunger_threshold", 0.5, 2.0
        )
        p.toxin_threshold = trial.suggest_float(
            "toxin_threshold", 0.02, 0.40
        )
        p.eta = trial.suggest_float(
            "eta", 5.0, 45.0, log=True
        )

    return p


def get_metrics(summary: dict):
    m = summary["normalized_spatial_fit"]
    return (
        float(m["pearson_minmax01_spatial_fit"]),
        float(m["mse_minmax01_spatial_fit"]),
        float(m["spearman_minmax01_spatial_fit"]),
        float(m["top20_hotspot_overlap_fraction"]),
        float(m["rmse_minmax01_spatial_fit"]),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-file", default=str(Path(__file__).with_name("namibia_pde_model_total_death.py")))
    ap.add_argument(
        "--data-dir", "--base-dir", dest="base_dir",
        default=os.environ.get("NAMIBIA_DATA_DIR", str(ROOT / "data" / "namibia")),
        help="directory containing the grid, mortality and mask inputs",
    )
    ap.add_argument("--grid", default="ud_pois_plot.rdata")
    ap.add_argument("--mortality", default="Poisoning_GIS.csv")
    ap.add_argument("--mask", default="combinedgeo.tif")
    ap.add_argument("--parks", default="Namibia_National_Parks_Boundaries_2001.shp")
    ap.add_argument("--n-trials", type=int, default=60)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--study-name", default="w08_initial_total_death")
    ap.add_argument("--storage", help="Optuna storage URL; defaults to a SQLite DB in output directory")
    ap.add_argument(
        "--output-dir",
        default=str(Path(os.environ.get("NAMIBIA_RESULTS_DIR", ROOT / "results" / "namibia")) / "calibration"),
    )
    args = ap.parse_args()

    if args.n_trials < 1:
        ap.error("--n-trials must be positive")
    args.base_dir = str(Path(args.base_dir).expanduser().resolve())
    args.output_dir = str(Path(args.output_dir).expanduser().resolve())
    if not Path(args.model_file).expanduser().is_file():
        ap.error(f"model file not found: {args.model_file}")
    for label, name in (("grid", args.grid), ("mortality", args.mortality), ("mask", args.mask)):
        candidate = Path(name).expanduser()
        if not candidate.is_absolute():
            candidate = Path(args.base_dir) / candidate
        if not candidate.is_file():
            ap.error(f"{label} file not found: {candidate}")

    import optuna

    model = load_model(args.model_file)
    os.makedirs(args.output_dir, exist_ok=True)
    if args.storage is None:
        args.storage = "sqlite:///" + str(Path(args.output_dir) / "w08_initial_total_death.db")

    def objective(trial: optuna.Trial) -> float:
        p = make_params(model, args, trial)

        try:
            summary = model.main(p)
            pearson, mse, spearman, top20, rmse = get_metrics(summary)

            if not np.isfinite(pearson) or not np.isfinite(mse):
                print(f"Trial {trial.number:3d} failed: non-finite metric")
                return -1e6

            score = pearson + 0.2 * top20 - 0.5 * mse

            trial.set_user_attr("pearson", pearson)
            trial.set_user_attr("mse", mse)
            trial.set_user_attr("rmse", rmse)
            trial.set_user_attr("spearman", spearman)
            trial.set_user_attr("top20", top20)
            trial.set_user_attr(
                "final_population_mass",
                float(summary["final_population_mass"]),
            )
            trial.set_user_attr(
                "final_cumulative_total_mortality_mass",
                float(summary["final_cumulative_total_mortality_mass"]),
            )
            trial.set_user_attr(
                "final_cumulative_poison_mortality_mass",
                float(summary["final_cumulative_poison_mortality_mass"]),
            )

            print(
                f"Trial {trial.number:3d} | "
                f"Pearson={pearson:.4f} | "
                f"Spearman={spearman:.4f} | "
                f"Top20={top20:.4f} | "
                f"MSE={mse:.5f} | "
                f"RMSE={rmse:.5f} | "
                f"Score={score:.5f}"
            )

            return float(score)

        except Exception as exc:
            trial.set_user_attr("error", repr(exc))
            print(f"Trial {trial.number:3d} failed: {exc}")
            return -1e6

    sampler = optuna.samplers.TPESampler(
        seed=args.seed,
        multivariate=True,
    )

    study = optuna.create_study(
        study_name=args.study_name,
        storage=args.storage,
        load_if_exists=True,
        direction="maximize",
        sampler=sampler,
    )

    study.optimize(
        objective,
        n_trials=args.n_trials,
        n_jobs=1,
        gc_after_trial=True,
    )

    best = study.best_trial

    if best.value is None or best.value <= -1e6:
        raise RuntimeError("All calibration trials failed; inspect the trial errors in the Optuna database.")

    print("\nBEST TRIAL")
    print(f"trial={best.number}")
    print(f"score={best.value:.6f}")
    print(f"Pearson={best.user_attrs.get('pearson')}")
    print(f"Spearman={best.user_attrs.get('spearman')}")
    print(f"Top20={best.user_attrs.get('top20')}")
    print(f"MSE={best.user_attrs.get('mse')}")
    print(f"RMSE={best.user_attrs.get('rmse')}")
    print(f"Final population mass={best.user_attrs.get('final_population_mass')}")
    print(
        "Final cumulative total mortality mass="
        f"{best.user_attrs.get('final_cumulative_total_mortality_mass')}"
    )
    print(
        "Final cumulative poison mortality mass="
        f"{best.user_attrs.get('final_cumulative_poison_mortality_mass')}"
    )
    print("\nBest parameters:")
    print(json.dumps(best.params, indent=2))

    study.trials_dataframe().to_csv(
        os.path.join(args.output_dir, "w08_initial_calibration_trials.csv"),
        index=False,
    )

    with open(
        os.path.join(args.output_dir, "best_w08_initial_calibration.json"),
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(
            {
                "trial": best.number,
                "objective": best.value,
                "metrics": best.user_attrs,
                "best_parameters": best.params,
                "fixed_initial_condition": "mixed",
                "uniform_weight": 0.8,
                "initial_condition_definition": "u0 = 0.2*UD + 0.8*Uniform",
                "objective_definition": "Pearson + 0.2*Top20 hotspot overlap - 0.5*MSE",
            },
            f,
            indent=2,
        )

    # Re-run best parameterization with normal outputs and plots.
    p = make_params(model, args, trial=None)
    p.search_radius_km_cognition = best.params["search_radius_km"]
    p.sigma_factor = best.params["sigma_factor"]
    p.diffusion_rate = best.params["diffusion_rate"]
    p.directed_rate = best.params["directed_rate"]
    p.d1 = best.params["d1"]
    p.decision_epsilon = best.params["decision_epsilon"]
    p.hunger_threshold = best.params["hunger_threshold"]
    p.toxin_threshold = best.params["toxin_threshold"]
    p.eta = best.params["eta"]
    p.output_dir_name = args.output_dir
    p.save_outputs = True
    p.make_plots = True
    p.verbose = True

    final_summary = model.main(p)

    with open(
        os.path.join(args.output_dir, "best_w08_initial_full_summary.json"),
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(final_summary, f, indent=2)


if __name__ == "__main__":
    main()
