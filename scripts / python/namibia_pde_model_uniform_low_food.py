"""Namibia hunger-structured PDE / graph model: uniform low-food scenario.

Uses uniform_food_level=0.025737 on the normalized resource scale and keeps
the heterogeneous toxin field. Use run_model.py for portable paths.
"""

from __future__ import annotations

import json

import os

import warnings

from dataclasses import asdict, dataclass

from typing import Literal

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

import pyreadr

import rasterio

from rasterio.crs import CRS

from rasterio.warp import transform as rio_transform

from scipy.sparse import csr_matrix, eye

from scipy.spatial import cKDTree

from scipy.stats import spearmanr

try:

    import geopandas as gpd

except ImportError:

    gpd = None


DomainMode = Literal["ud_positive", "all_valid"]

Condition = Literal["Random", "NoCognition", "Cognition"]


@dataclass

class Params:

    # Final parameterization selected from Pearson+MSE calibration.

    # Trial 97: Pearson=0.7157, MSE=0.0355, Top20=0.706,

    # Spearman=0.6035 (diagnostic).

    # Inputs / outputs

    base_dir: str = "."

    grid_rdata_name: str = "ud_pois_plot.rdata"

    mortality_csv_name: str = "Poisoning_GIS.csv"

    namibia_mask_tif_name: str = "combinedgeo.tif"

    parks_shp_name: str = "Namibia_National_Parks_Boundaries_2001.shp"

    output_dir_name: str = "pde_mortality_namibia_output_uniform_low_food"

    # Set these False in calibration to speed up Optuna.

    save_outputs: bool = True

    make_plots: bool = True

    verbose: bool = True

    # Domain

    # "all_valid" = all valid Namibia cells; UD determines initial condition only.

    # "ud_positive" = restrict PDE domain to cells with UD > 0.

    domain_mode: DomainMode = "all_valid"

    # Spatial fields
    toxin_col: str = "pois_hwc_all"
    food_col: str = "food_resc"

    # Scenario: uniform low food on the normalized [0,1] scale,
    # and original spatially heterogeneous toxin field.
    uniform_food_level: float = 0.025737

    # Movement condition used for the Namibia application.

    # "Cognition" means directed nonlocal perception is active.

    condition: Condition = "Cognition"

    k_neighbors: int = 8

    # Local crowding / competition radius R_u (km).

    # This is biologically distinct from the nonlocal perception radius R.

    # The 30-km value is a provisional baseline and can later be calibrated

    # or subjected to sensitivity analysis.

    crowding_radius_km: float = 30.0

    # Movement / perception.

    # Radius is specified directly in physical units (km), not "number of cells".

    # Final cognition radius selected from calibration.

    search_radius_km_cognition: float = 53.595283559290

    search_radius_km_nocognition: float = 15.0

    # sigma_R = sigma_factor * R.

    sigma_factor: float = 0.8851218263931865

    # Behavioural thresholds.

    # v_h = 1 corresponds to roughly 3 days of fasting when a0 = 1/3 day^-1.

    hunger_threshold: float = 0.808976066888251

    toxin_threshold: float = 0.22105985007781886

    decision_epsilon: float = 0.032274814350080394

    # Graph transition rates (day^-1).

    diffusion_rate: float = 0.00014073696431122817

    directed_rate: float = 0.023689607034513824 #0.02669786875925223

    # Hunger-structured dynamics:

    # G(x,t,v,U) = a0 - eta*s(x)/(1 + beta*C(x,t)) * v

    #

    # eta and beta follow the current continuum formulation.

    # NOTE: because the current continuum equation contains the factor v,

    # this is a model update rather than a purely algebraic rewrite of the

    # old a0 - s/(b + u_local) equation.

    # Hunger state is dimensionless; time is measured in days.

    # v = 0 is interpreted as fully satiated, and v_h = 1 as the food-seeking

    # threshold reached after about 3 food-free days.

    v_start: float = 0.0

    v_min: float = 0.0

    v_max: float = 12.0

    n_hunger_bins: int = 61

    newborn_hunger: float = 0.0

    # Hunger accumulation rate: 1/3 day^-1 -> v rises from 0 to 1 in 3 days

    # when the food-mediated term is absent.

    a0: float = 1.0 / 3.0

    # Current hunger-dynamics parameterization.

    eta: float = 40.927512727237236

    beta: float = 1.0

    # Population is kept entirely on a normalized mass scale.

    # Initial total population mass = 1.

    # No arbitrary conversion to animal counts is used.

    # Birth / death reaction. Rates are per day.

    # Default recruitment is centered within the literature-informed search range.

    r: float = 0.0004944139192066642

    # Effective background mortality hazard.

    # The final calibration range is derived from annual survival

    # S = 0.852--0.999 using d1 = -ln(S)/365.25.

    d1: float = 7.646949137034071e-06

    # Shared behavioural/mortality sensitivities.

    # Hunger remains active with the calibrated full-model value.

    # Toxin sensitivity remains active in behaviour and mortality.

    d_v: float = 0.6360304243122281

    d_tau: float = 2.133000054600129

    # Food-dependent birth saturation

    birth_K_min: float = 0.0010977030337361

    birth_K_max: float = 0.009152941020177215

    # Initial condition

    initial_condition: str = "mixed"       # "ud", "uniform", or "mixed"

    uniform_weight: float = 0.8

    # Time integration.

    # One model time unit = one day.

    # dt = 0.05 day = 1.2 hours; 3600 steps = 180 days.

    time_unit: str = "day"

    dt: float = 0.05

    n_steps: int = 3600

    # Evaluation radius used to compare local modeled and observed mortality

    # burdens. It is fixed during biological calibration.

    event_compare_radius_km: float = 50.0

    # Plot limits copied from the original Namibia R plotting window

    xlim_min: float = 14

    xlim_max: float = 20.5

    ylim_min: float = -29.0

    ylim_max: float = -19

# ---------------------------------------------------------------------

# General helpers

# ---------------------------------------------------------------------

def normalize_nonnegative(x: np.ndarray) -> np.ndarray:

    x = np.asarray(x, dtype=float)

    y = np.maximum(x, 0.0)

    total = float(y.sum())

    return y / total if total > 0 else y


def rescale01(x: np.ndarray) -> np.ndarray:

    x = np.asarray(x, dtype=float)

    lo = np.nanmin(x)

    hi = np.nanmax(x)

    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:

        return np.zeros_like(x, dtype=float)

    return (x - lo) / (hi - lo)


def normalize_minmax01(x: np.ndarray) -> np.ndarray:

    """

    Rescale values to [0, 1].

    0 = lowest value

    1 = highest value

    """

    x = np.asarray(x, dtype=float)

    lo = float(np.nanmin(x))

    hi = float(np.nanmax(x))

    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:

        return np.zeros_like(x, dtype=float)

    return (x - lo) / (hi - lo)


def lonlat_to_km(

    lon: np.ndarray,

    lat: np.ndarray,

    lon0: float | None = None,

    lat0: float | None = None,

) -> tuple[np.ndarray, np.ndarray, float, float]:

    lon = np.asarray(lon, dtype=float)

    lat = np.asarray(lat, dtype=float)

    if lon0 is None:

        lon0 = float(np.nanmean(lon))

    if lat0 is None:

        lat0 = float(np.nanmean(lat))

    x = (lon - lon0) * 111.32 * np.cos(np.deg2rad(lat0))

    y = (lat - lat0) * 110.57

    return x, y, lon0, lat0


def load_rdata_frame(path: str) -> pd.DataFrame:

    objects = pyreadr.read_r(path)

    if "df" in objects:

        return objects["df"].copy()

    if "data" in objects:

        return objects["data"].copy()

    raise ValueError(

        f"Expected an object named 'df' or 'data' in {path}; found {list(objects.keys())}"

    )


# ---------------------------------------------------------------------

# Input preparation

# ---------------------------------------------------------------------

def sample_namibia_mask(grid: pd.DataFrame, mask_path: str) -> np.ndarray:

    """

    Reproduce the R logic:

        ex <- raster::extract(combinedgeo.tif, points)

        keep points where ex is not NA

    """

    if not os.path.exists(mask_path):

        raise FileNotFoundError(f"Could not find Namibia mask raster: {mask_path}")

    xs = grid["x"].to_numpy(dtype=float)

    ys = grid["y"].to_numpy(dtype=float)

    with rasterio.open(mask_path) as src:

        # Grid coordinates are lon/lat. Transform only if the mask raster is not lon/lat.

        if src.crs is not None and src.crs != CRS.from_epsg(4326):

            tx, ty = rio_transform("EPSG:4326", src.crs, xs.tolist(), ys.tolist())

            sample_xy = list(zip(tx, ty))

        else:

            sample_xy = list(zip(xs, ys))

        vals = np.array([sample[0] for sample in src.sample(sample_xy)], dtype=float)

        nodata = src.nodata

    valid = np.isfinite(vals)

    if nodata is not None and np.isfinite(nodata):

        valid &= vals != nodata

    return valid


def prepare_grid(raw_df: pd.DataFrame, params: Params) -> pd.DataFrame:

    required = {"ud", "x", "y", params.toxin_col, params.food_col}

    missing = required.difference(raw_df.columns)

    if missing:

        raise ValueError(f"Grid file is missing required columns: {sorted(missing)}")

    keep = raw_df[["ud", "x", "y", params.toxin_col, params.food_col]].copy()

    keep = keep.dropna(subset=["x", "y", params.toxin_col, params.food_col])

    # Sum UD by cell, keep raster covariates as cell-level means.

    grid = (

        keep.groupby(["x", "y"], as_index=False)

        .agg(

            ud=("ud", "sum"),

            toxin_raw=(params.toxin_col, "mean"),

            food_raw=(params.food_col, "mean"),

        )

        .sort_values(["x", "y"])

        .reset_index(drop=True)

    )

    # Keep Namibia cells only.

    mask_path = os.path.join(params.base_dir, params.namibia_mask_tif_name)

    grid = grid[sample_namibia_mask(grid, mask_path)].copy().reset_index(drop=True)

    grid["ud"] = grid["ud"].fillna(0.0)

#    grid["food_model"] = rescale01(grid["food_raw"].to_numpy())
    # Constant resource level at every valid Namibia spatial cell.
    grid["food_model"] = np.full(
        len(grid), params.uniform_food_level, dtype=float
    )

    # Preserve the calibrated, spatially heterogeneous toxin field.
    grid["toxin_model"] = rescale01(
        grid["toxin_raw"].to_numpy(dtype=float)
    )

    if params.domain_mode == "ud_positive":

        grid = grid[grid["ud"] > 0].copy().reset_index(drop=True)

    elif params.domain_mode == "all_valid":

        grid = grid.copy().reset_index(drop=True)

    else:

        raise ValueError(f"Unknown domain_mode: {params.domain_mode}")

    if grid.empty:

        raise ValueError("No Namibia cells remain after applying the domain filter.")

    if float(grid["ud"].sum()) <= 0:

        raise ValueError("UD sums to zero after filtering; cannot build initial condition.")

    return grid


def prepare_mortality(path: str) -> pd.DataFrame:

    if not os.path.exists(path):

        raise FileNotFoundError(f"Could not find mortality CSV: {path}")

    mort = pd.read_csv(path)

    required = {"X", "Y", "N_Dead"}

    missing = required.difference(mort.columns)

    if missing:

        raise ValueError(f"Mortality CSV is missing required columns: {sorted(missing)}")

    mort = mort.dropna(subset=["X", "Y", "N_Dead"])

    mort = mort[mort["N_Dead"] > 0].reset_index(drop=True)

    if mort.empty:

        raise ValueError(

            "No mortality records remain after filtering valid coordinates and N_Dead > 0."

        )


    return mort


# ---------------------------------------------------------------------

# Movement operators

# ---------------------------------------------------------------------

def build_diffusion_matrix(

    x_km: np.ndarray,

    y_km: np.ndarray,

    params: Params,

) -> tuple[csr_matrix, np.ndarray, np.ndarray, dict[str, float]]:

    """

    Build the local graph random-walk operator used for spatial diffusion.

    Two distinct spatial scales are reported:

    * median_nearest_neighbor_km:

        median distance from each valid spatial node to its nearest other node.

        This is the most useful diagnostic for the empirical grid resolution.

    * graph_kernel_sigma_km:

        median distance across all graph edges used in the k-nearest-neighbor

        diffusion operator. This controls the Gaussian edge weights and should

        not be interpreted as the physical "cell size".

    """

    coords = np.column_stack([x_km, y_km])

    tree = cKDTree(coords)

    k = min(params.k_neighbors + 1, len(coords))

    dists, inds = tree.query(coords, k=k)

    if dists.ndim == 1:

        dists = dists[:, None]

        inds = inds[:, None]

    # First column is the node itself; second column is the nearest other node.

    nearest = dists[:, 1] if dists.shape[1] > 1 else np.asarray([], dtype=float)

    nearest = nearest[np.isfinite(nearest) & (nearest > 0)]

    median_nearest_neighbor_km = (

        float(np.median(nearest)) if nearest.size else 1.0

    )

    # Graph edge-length scale used only for the random-walk weights.

    positive_edges = dists[:, 1:][

        np.isfinite(dists[:, 1:]) & (dists[:, 1:] > 0)

    ]

    graph_kernel_sigma_km = (

        float(np.median(positive_edges)) if positive_edges.size else 1.0

    )

    graph_kernel_sigma_km = max(graph_kernel_sigma_km, 1e-6)

    rows: list[int] = []

    cols: list[int] = []

    vals: list[float] = []

    for i in range(len(coords)):

        js = np.asarray(inds[i, 1:], dtype=int)

        ds = np.asarray(dists[i, 1:], dtype=float)

        valid = np.isfinite(ds) & (ds > 0) & (js >= 0) & (js < len(coords))

        js = js[valid]

        ds = ds[valid]

        if js.size == 0:

            rows.append(i)

            cols.append(i)

            vals.append(1.0)

            continue

        w = np.exp(-0.5 * (ds / graph_kernel_sigma_km) ** 2)

        sw = float(np.sum(w))

        if sw <= 0 or not np.isfinite(sw):

            w = np.ones_like(w, dtype=float) / len(w)

        else:

            w = w / sw

        rows.extend([i] * len(js))

        cols.extend(js.tolist())

        vals.extend(w.tolist())

    P_diff = csr_matrix(

        (vals, (rows, cols)),

        shape=(len(coords), len(coords)),

    )

    spacing_meta = {

        "median_nearest_neighbor_km": median_nearest_neighbor_km,

        "graph_kernel_sigma_km": graph_kernel_sigma_km,

    }

    return P_diff, dists, inds, spacing_meta

def precompute_weighted_targets(

    x_km: np.ndarray,

    y_km: np.ndarray,

    field: np.ndarray,

    search_radius_km: float,

    sigma_km: float,

    optimize_type: Literal["max", "min"],

    score_mode: Literal["weighted", "distance_penalized"] = "weighted",

    tol: float = 1e-10,

) -> np.ndarray:

    """

    Precompute the preferred nonlocal target at every spatial cell.

    Food seeking uses the distance-weighted food score

        Phi_s(z;x) = w_R(z;x) * s(z),

    and therefore maximizes ``field * weight``.

    Toxin avoidance uses the distance-penalized toxin score

        Phi_tau(z;x) = tau(z) / w_R(z;x),

    and therefore minimizes ``field / weight``.

    The latter prevents a distant candidate from appearing artificially safe

    simply because its Gaussian perception weight is small.

    """

    coords = np.column_stack([x_km, y_km])

    tree = cKDTree(coords)

    targets = np.arange(len(coords), dtype=int)

    sigma_km = max(float(sigma_km), 1e-12)

    for i, xy in enumerate(coords):

        cand = np.asarray(tree.query_ball_point(xy, r=search_radius_km), dtype=int)

        if cand.size == 0:

            continue

        # A spatially uniform environmental cue has no preferred
        # directional target, even though Gaussian distance weights vary.
        if np.ptp(field[cand]) < tol:
            targets[i] = i
            continue

        d2 = np.sum((coords[cand] - xy) ** 2, axis=1)

        weights = np.exp(-d2 / (2.0 * sigma_km**2))

        if score_mode == "weighted":

            score = field[cand] * weights

        elif score_mode == "distance_penalized":

            score = field[cand] / np.maximum(weights, 1e-12)

        else:

            raise ValueError(

                "score_mode must be 'weighted' or 'distance_penalized'"

            )

        if score.size == 0:

            continue

        if abs(float(np.max(score)) - float(np.min(score))) < tol:

            targets[i] = i

            continue

        if optimize_type == "max":

            optimum = float(np.max(score))

            tied = np.flatnonzero(np.isclose(score, optimum, rtol=0.0, atol=tol))

        elif optimize_type == "min":

            optimum = float(np.min(score))

            tied = np.flatnonzero(np.isclose(score, optimum, rtol=0.0, atol=tol))

        else:

            raise ValueError("optimize_type must be 'max' or 'min'")

        # If several candidates have the same optimal score, choose the nearest

        # one to avoid an arbitrary dependence on cKDTree return order.

        if tied.size > 1:

            best_local = int(tied[int(np.argmin(d2[tied]))])

        else:

            best_local = int(tied[0])

        targets[i] = int(cand[best_local])

    return targets


def build_directional_transition_matrix(

    x_km: np.ndarray,

    y_km: np.ndarray,

    neighbor_inds: np.ndarray,

    targets: np.ndarray,

) -> csr_matrix:

    """

    Convert nonlocal ABM targets into one-neighbour graph transport steps.

    """

    coords = np.column_stack([x_km, y_km])

    rows: list[int] = []

    cols: list[int] = []

    vals: list[float] = []

    for i in range(len(coords)):

        target = int(targets[i])

        if target == i:

            rows.append(i)

            cols.append(i)

            vals.append(1.0)

            continue

        target_vec = coords[target] - coords[i]

        target_norm = float(np.linalg.norm(target_vec))

        if target_norm <= 1e-12:

            rows.append(i)

            cols.append(i)

            vals.append(1.0)

            continue

        target_unit = target_vec / target_norm

        js = np.asarray(neighbor_inds[i, 1:], dtype=int)

        if js.size == 0:

            rows.append(i)

            cols.append(i)

            vals.append(1.0)

            continue

        edge_vecs = coords[js] - coords[i]

        edge_norms = np.linalg.norm(edge_vecs, axis=1)

        valid = edge_norms > 1e-12

        if not np.any(valid):

            rows.append(i)

            cols.append(i)

            vals.append(1.0)

            continue

        edge_units = edge_vecs[valid] / edge_norms[valid, None]

        align = edge_units @ target_unit

        valid_js = js[valid]

        best_idx = int(np.argmax(align))

        if float(align[best_idx]) <= 0:

            rows.append(i)

            cols.append(i)

            vals.append(1.0)

        else:

            rows.append(i)

            cols.append(int(valid_js[best_idx]))

            vals.append(1.0)

    return csr_matrix((vals, (rows, cols)), shape=(len(coords), len(coords)))


def prepare_movement_operators(

    x_km: np.ndarray,

    y_km: np.ndarray,

    food: np.ndarray,

    tau: np.ndarray,

    params: Params,

) -> tuple[csr_matrix, csr_matrix, csr_matrix, dict[str, float | str]]:

    """

    Prepare spatial diffusion and directed-movement graph operators.

    The biological perception radius R is specified directly in kilometres.

    This avoids tying the interpretation of R to the numerical grid spacing.

    """

    P_diff, _dists, inds, spacing_meta = build_diffusion_matrix(

        x_km, y_km, params

    )

    if params.condition == "Cognition":

        search_radius_km = float(params.search_radius_km_cognition)

    elif params.condition == "NoCognition":

        search_radius_km = float(params.search_radius_km_nocognition)

    elif params.condition == "Random":

        search_radius_km = 0.0

    else:

        raise ValueError(f"Unknown condition: {params.condition}")

    sigma_km = float(

        max(search_radius_km * params.sigma_factor, 1e-12)

    )

    if params.condition == "Random" or search_radius_km <= 0:

        identity = eye(len(x_km), format="csr")

        Q_food = identity

        Q_toxin = identity

    else:

        food_targets = precompute_weighted_targets(

            x_km,

            y_km,

            food,

            search_radius_km=search_radius_km,

            sigma_km=sigma_km,

            optimize_type="max",

            score_mode="weighted",

        )

        toxin_targets = precompute_weighted_targets(

            x_km,

            y_km,

            tau,

            search_radius_km=search_radius_km,

            sigma_km=sigma_km,

            optimize_type="min",

            score_mode="distance_penalized",

        )

        Q_food = build_directional_transition_matrix(

            x_km, y_km, inds, food_targets

        )

        Q_toxin = build_directional_transition_matrix(

            x_km, y_km, inds, toxin_targets

        )

    median_nn = float(spacing_meta["median_nearest_neighbor_km"])

    radius_in_median_steps = (

        search_radius_km / median_nn if median_nn > 0 else float("nan")

    )

    meta = {

        **spacing_meta,

        "search_radius_km": search_radius_km,

        "search_radius_in_median_nn_steps": radius_in_median_steps,

        "sigma_km": sigma_km,

        "sigma_factor": float(params.sigma_factor),

        "sigma_definition": "sigma_R = sigma_factor * R",

        "food_perception_score": "w_R_times_s",

        "toxin_perception_score": "tau_over_w_R",

        "distance_unit": "km",

    }

    return P_diff, Q_food, Q_toxin, meta

def build_local_population_matrix(

    x_km: np.ndarray,

    y_km: np.ndarray,

    params: Params,

) -> tuple[csr_matrix, dict[str, float | str]]:

    r"""

    Build the binary local-neighborhood matrix used for crowding.

    The crowding variable is

        C_i(t) = sum_{j: d_ij <= R_u} U_j(t),

    where U_j(t) is the total population mass in spatial cell j after

    summing over all hunger states.

    Important:

    - The matrix is NOT row-normalized.

    - Every included neighbor contributes weight 1.

    - The focal cell i is included because d_ii = 0 <= R_u.

    - Since the initial total population mass is normalized to 1,

      C_i(t) is interpretable as the fraction of the initial total

      population currently contained in the R_u-neighborhood of cell i.

    """

    coords = np.column_stack([x_km, y_km])

    tree = cKDTree(coords)

    radius_km = float(params.crowding_radius_km)

    if radius_km <= 0:

        raise ValueError("crowding_radius_km must be positive.")

    neighborhoods = tree.query_ball_point(coords, r=radius_km)

    rows: list[int] = []

    cols: list[int] = []

    vals: list[float] = []

    neighbor_counts = np.zeros(len(coords), dtype=int)

    for i, js in enumerate(neighborhoods):

        js = np.asarray(js, dtype=int)

        if js.size == 0:

            js = np.asarray([i], dtype=int)

        neighbor_counts[i] = int(js.size)

        rows.extend([i] * len(js))

        cols.extend(js.tolist())

        vals.extend([1.0] * len(js))

    A_local = csr_matrix(

        (vals, (rows, cols)),

        shape=(len(coords), len(coords)),

    )

    meta = {

        "crowding_definition": "C_i(t) = sum_{j: d_ij <= R_u} U_j(t)",

        "crowding_radius_km": radius_km,

        "crowding_matrix_row_normalized": False,

        "crowding_focal_cell_included": True,

        "crowding_min_cells_per_neighborhood": int(np.min(neighbor_counts)),

        "crowding_median_cells_per_neighborhood": float(np.median(neighbor_counts)),

        "crowding_max_cells_per_neighborhood": int(np.max(neighbor_counts)),

    }

    return A_local, meta


# ---------------------------------------------------------------------

# Hunger-structured helpers

# ---------------------------------------------------------------------

def build_hunger_grid(params: Params) -> np.ndarray:

    """Uniform hunger grid used to approximate the structured variable v."""

    if params.n_hunger_bins < 2:

        raise ValueError("n_hunger_bins must be at least 2.")

    if params.v_max <= params.v_min:

        raise ValueError("v_max must be greater than v_min.")

    return np.linspace(

        float(params.v_min),

        float(params.v_max),

        int(params.n_hunger_bins),

    )


def point_mass_on_hunger_grid(v_grid: np.ndarray, value: float) -> np.ndarray:

    """

    Discrete approximation of a delta distribution at a prescribed hunger value.

    """

    p = np.zeros_like(v_grid, dtype=float)

    idx = int(np.argmin(np.abs(v_grid - float(value))))

    p[idx] = 1.0

    return p


def transport_hunger_by_characteristics(

    u_xv: np.ndarray,

    v_grid: np.ndarray,

    food: np.ndarray,

    crowd_count: np.ndarray,

    params: Params,

) -> np.ndarray:

    r"""

    Conservative characteristic transport for

        \partial_t u + \partial_v(G u) = 0,

    with

        G(x,t,v,U)

        = a0 - eta*s(x)/(1 + beta*Crowd[U](x,t)) * v.

    During one operator-split time step, the crowding-dependent coefficient is

    held fixed at each spatial cell. The characteristic ODE

        dv/dt = a0 - k(x,t) v

    is then solved analytically and the population mass in each hunger bin is

    remapped to the destination hunger bins. This avoids a restrictive explicit

    CFL condition in hunger space and preserves mass (up to floating-point

    roundoff and clipping at the finite v-domain boundaries).

    """

    u_xv = np.asarray(u_xv, dtype=float)

    v_grid = np.asarray(v_grid, dtype=float)

    food = np.asarray(food, dtype=float)

    crowd_count = np.asarray(crowd_count, dtype=float)

    n_x, n_v = u_xv.shape

    if n_v != len(v_grid):

        raise ValueError("u_xv second dimension must match len(v_grid).")

    dv = float(v_grid[1] - v_grid[0])

    v_min = float(v_grid[0])

    v_max = float(v_grid[-1])

    # k >= 0 because food, eta, beta, and crowding are nonnegative.

    k = (

        params.eta

        * np.maximum(food, 0.0)

        / (1.0 + params.beta * np.maximum(crowd_count, 0.0))

    )

    decay = np.exp(-k * params.dt)

    v0 = v_grid[None, :]

    small_k = k <= 1e-14

    equilibrium = np.zeros_like(k)

    equilibrium[~small_k] = params.a0 / k[~small_k]

    v_dest = equilibrium[:, None] + (

        v0 - equilibrium[:, None]

    ) * decay[:, None]

    if np.any(small_k):

        v_dest[small_k, :] = v0 + params.a0 * params.dt

    # Finite hunger domain: accumulate any characteristic reaching outside the

    # numerical interval in the closest boundary bin.

    v_dest = np.clip(v_dest, v_min, v_max)

    j_float = (v_dest - v_min) / dv

    j_lo = np.floor(j_float).astype(np.int64)

    j_lo = np.clip(j_lo, 0, n_v - 1)

    j_hi = np.minimum(j_lo + 1, n_v - 1)

    frac_hi = np.clip(j_float - j_lo, 0.0, 1.0)

    frac_hi[j_lo == j_hi] = 0.0

    frac_lo = 1.0 - frac_hi

    row_offset = (np.arange(n_x, dtype=np.int64) * n_v)[:, None]

    flat_lo = (row_offset + j_lo).ravel()

    flat_hi = (row_offset + j_hi).ravel()

    out_flat = np.bincount(

        flat_lo,

        weights=(u_xv * frac_lo).ravel(),

        minlength=n_x * n_v,

    )

    out_flat += np.bincount(

        flat_hi,

        weights=(u_xv * frac_hi).ravel(),

        minlength=n_x * n_v,

    )

    return out_flat.reshape(n_x, n_v)


# ---------------------------------------------------------------------

# PDE simulation

# ---------------------------------------------------------------------

def simulate_pde(

    u0: np.ndarray,

    food: np.ndarray,

    tau: np.ndarray,

    P_diff: csr_matrix,

    Q_food: csr_matrix,

    Q_toxin: csr_matrix,

    A_local: csr_matrix,

    params: Params,

) -> dict[str, np.ndarray | float]:

    r"""

    Run the hunger-structured spatial model.

    State

    -----

    ``u_xv[i, j]`` is the normalized population mass in spatial cell i and

    hunger bin j, approximating u(x,t,v).

    The total spatial population density is

        U(x,t) = integral u(x,t,v) dv,

    represented numerically by summing the hunger-bin masses.

    Hunger drift

    ------------

        G = a0 - eta*s/(1 + beta*C) * v, with C_i(t) = sum_{j:d_ij<=R_u} U_j(t).

    Directed movement

    -----------------

    Food seeking is active when

        v > v_h

        and d_v*v - d_tau*tau >= epsilon,

    while toxin avoidance is active when

        tau > tau_h

        and d_tau*tau - d_v*v >= epsilon.

    Hence the epsilon interval is a diffusion-only decision deadband.

    """

    spatial_u0 = np.asarray(u0, dtype=float).copy()

    food = np.asarray(food, dtype=float)

    tau = np.asarray(tau, dtype=float)

    v_grid = build_hunger_grid(params)

    initial_hunger_distribution = point_mass_on_hunger_grid(

        v_grid, params.v_start

    )

    newborn_hunger_distribution = point_mass_on_hunger_grid(

        v_grid, params.newborn_hunger

    )

    # Structured initial condition u(x,0,v).

    u_xv = spatial_u0[:, None] * initial_hunger_distribution[None, :]

    # u0 is normalized in main(), so the initial total population mass

    # should be 1. This is the reference population against which U and C

    # are interpreted; the state is NOT renormalized during the simulation.

    initial_total_mass = float(u_xv.sum())

    if not np.isclose(initial_total_mass, 1.0, rtol=1e-8, atol=1e-10):

        warnings.warn(

            f"Initial total population mass is {initial_total_mass:.12g}, "

            "not 1. The intended interpretation of C_i(t) as a fraction "

            "of the initial total population assumes normalized u0.",

            RuntimeWarning,

            stacklevel=2,

        )

    # Cumulative mortality is stored spatially after integrating over hunger.

    M_total = np.zeros_like(spatial_u0)

    M_baseline_hunger = np.zeros_like(spatial_u0)

    M_poison_excess = np.zeros_like(spatial_u0)

    time_history = [0.0]

    total_population_history = [

        float(u_xv.sum())

    ]

    deaths_this_step_history = [0.0]

    cumulative_total_deaths_history = [0.0]

    cumulative_poison_deaths_history = [0.0]

    warned_about_reaction_step = False

    warned_about_total_outflow = False

    v_row = v_grid[None, :]

    tau_col = tau[:, None]

    # Mortality hazards depend on hunger state and spatial toxin exposure.

    # The same d_v and d_tau used in the behavioural decision rule are also

    # used here as mortality multipliers (parsimonious shared-sensitivity

    # assumption).

    baseline_hunger_rate = params.d1 * (1.0 + params.d_v * v_row)

    poison_multiplier = params.d_tau * tau_col

    crowding_initial = np.asarray(

        A_local @ u_xv.sum(axis=1),

        dtype=float,

    ).ravel()

    for step in range(params.n_steps):

        # U(x,t): total mass at each spatial cell after integrating over hunger.

        U = u_xv.sum(axis=1)

        # Local population mass within radius R_u:

        #

        #     C_i(t) = sum_{j: d_ij <= R_u} U_j(t).

        #

        # A_local is binary and is NOT row-normalized. Because the initial

        # total population mass is normalized to 1, C_i(t) is the fraction

        # of the initial total modeled population currently present in the

        # local neighborhood around cell i.

        crowd_count = np.asarray(A_local @ U, dtype=float).ravel()

        # 1) Hunger-space transport:

        #    -partial_v [G(x,t,v,U) u].

        u_h = transport_hunger_by_characteristics(

            u_xv=u_xv,

            v_grid=v_grid,

            food=food,

            crowd_count=crowd_count,

            params=params,

        )

        # Hunger transport is conservative within each spatial cell, so U is

        # unchanged apart from floating-point roundoff.

        U_h = u_h.sum(axis=1)

        # 2) Recruitment B(x,U) pi_b(v).

        K_local = params.birth_K_min + (

            params.birth_K_max - params.birth_K_min

        ) * food

        birth_density_factor = K_local / (

            K_local + crowd_count + 1e-12

        )

        B_total = params.r * U_h * birth_density_factor

        birth = B_total[:, None] * newborn_hunger_distribution[None, :]

        # 3) Mortality M(x,v,tau)u.

        poison_excess_rate = baseline_hunger_rate * poison_multiplier

        death_rate = baseline_hunger_rate + poison_excess_rate

        max_fraction_removed = float(params.dt * np.max(death_rate))

        if max_fraction_removed > 1.0 and not warned_about_reaction_step:

            warnings.warn(

                "Explicit reaction step may be unstable: "

                "dt * max(death_rate) > 1. Reduce Params.dt.",

                RuntimeWarning,

                stacklevel=2,

            )

            warned_about_reaction_step = True

        # A sufficient (conservative) positivity check for the explicit

        # spatial/reaction update. Directed movement is active for at most one

        # behavioural mode at a time because of the epsilon deadband.

        directed_out_rate = (

            0.0 if params.condition == "Random" else params.directed_rate

        )

        max_outflow_fraction = params.dt * (

            params.diffusion_rate

            + directed_out_rate

            + float(np.max(death_rate))

        )

        if max_outflow_fraction > 1.0 and not warned_about_total_outflow:

            warnings.warn(

                "Explicit spatial/reaction step may lose positivity: "

                "dt * (diffusion_rate + directed_rate + max death rate) > 1. "

                "Reduce Params.dt or the transition rates.",

                RuntimeWarning,

                stacklevel=2,

            )

            warned_about_total_outflow = True

        death_total_xv = death_rate * u_h

        death_baseline_hunger_xv = baseline_hunger_rate * u_h

        death_poison_excess_xv = poison_excess_rate * u_h

        # 4) Random spatial movement acts on every hunger class.

        diffusion_move = params.diffusion_rate * (P_diff.T @ u_h - u_h)

        # 5) Hunger- and toxin-dependent directed movement.

        if params.condition == "Random":

            directed_move = np.zeros_like(u_h)

        else:

            hunger_pressure = params.d_v * v_row

            toxin_pressure = params.d_tau * tau_col

            food_mask = (

                (v_row > params.hunger_threshold)

                & (

                    hunger_pressure - toxin_pressure

                    >= params.decision_epsilon

                )

            )

            toxin_mask = (

                (tau_col > params.toxin_threshold)

                & (

                    toxin_pressure - hunger_pressure

                    >= params.decision_epsilon

                )

            )

            food_mass = u_h * food_mask

            toxin_mass = u_h * toxin_mask

            active_mass = food_mass + toxin_mass

            incoming_directed = (

                Q_food.T @ food_mass

                + Q_toxin.T @ toxin_mass

            )

            directed_move = params.directed_rate * (

                incoming_directed - active_mass

            )

        # Operator-split explicit spatial/reaction update.

        u_xv = u_h + params.dt * (

            diffusion_move

            + directed_move

            + birth

            - death_total_xv

        )

        u_xv[u_xv < 0] = 0.0

        # Integrate mortality over hunger state to obtain spatial mortality.

        death_total_x = death_total_xv.sum(axis=1)

        death_baseline_hunger_x = death_baseline_hunger_xv.sum(axis=1)

        death_poison_excess_x = death_poison_excess_xv.sum(axis=1)

        M_total += params.dt * death_total_x

        M_baseline_hunger += params.dt * death_baseline_hunger_x

        M_poison_excess += params.dt * death_poison_excess_x

        time_history.append((step + 1) * params.dt)

        total_population_history.append(

            float(u_xv.sum())

        )

        deaths_this_step_history.append(

            float(

                params.dt

                * death_total_x.sum()


            )

        )

        cumulative_total_deaths_history.append(

            float(M_total.sum() )

        )

        cumulative_poison_deaths_history.append(

            float(M_poison_excess.sum() )

        )

        if float(u_xv.sum()) <= 1e-15:

            break

    U_final = u_xv.sum(axis=1)

    crowding_final = np.asarray(

        A_local @ U_final,

        dtype=float,

    ).ravel()

    total_final_mass = float(u_xv.sum())

    upper_hunger_boundary_mass_fraction = (

        float(u_xv[:, -1].sum()) / total_final_mass

        if total_final_mass > 0

        else 0.0

    )

    if upper_hunger_boundary_mass_fraction > 0.01:

        warnings.warn(

            "More than 1% of the final population mass is in the highest "

            "hunger bin. The finite v_max may be clipping the hunger "

            "distribution; consider increasing v_max or revisiting the "

            "hunger dynamics.",

            RuntimeWarning,

            stacklevel=2,

        )

    mean_hunger = np.divide(

        u_xv @ v_grid,

        U_final,

        out=np.zeros_like(U_final),

        where=U_final > 1e-15,

    )

    return {

        # Spatially integrated population quantities.

        "final_alive_raw": U_final,

        "final_alive_prob": normalize_nonnegative(U_final),

        "final_hunger": mean_hunger,

        # Full hunger-structured final state.

        "hunger_grid": v_grid,

        "final_structured_density_raw": u_xv,

        # Spatial mortality integrated over hunger state.

        "total_mortality_raw": M_total,

        "total_mortality_prob": normalize_nonnegative(M_total),

        "baseline_hunger_mortality_raw": M_baseline_hunger,

        "baseline_hunger_mortality_prob": normalize_nonnegative(

            M_baseline_hunger

        ),

        "poison_excess_mortality_raw": M_poison_excess,

        "poison_excess_mortality_prob": normalize_nonnegative(

            M_poison_excess

        ),

        "final_total_mass": float(U_final.sum()),

        "upper_hunger_boundary_mass_fraction": float(

            upper_hunger_boundary_mass_fraction

        ),

        "initial_total_population_mass": float(initial_total_mass),

        "crowding_definition": "C_i(t) = sum_{j: d_ij <= R_u} U_j(t)",

        "crowding_radius_km": float(params.crowding_radius_km),

        "crowding_interpretation": (

            "local population mass as a fraction of the initial total "

            "population mass (initial total mass is normalized to 1)"

        ),

        "crowding_initial_quantiles": {

            "q00": float(np.quantile(crowding_initial, 0.00)),

            "q25": float(np.quantile(crowding_initial, 0.25)),

            "q50": float(np.quantile(crowding_initial, 0.50)),

            "q75": float(np.quantile(crowding_initial, 0.75)),

            "q90": float(np.quantile(crowding_initial, 0.90)),

            "q99": float(np.quantile(crowding_initial, 0.99)),

            "q100": float(np.quantile(crowding_initial, 1.00)),

        },

        "crowding_final_quantiles": {

            "q00": float(np.quantile(crowding_final, 0.00)),

            "q25": float(np.quantile(crowding_final, 0.25)),

            "q50": float(np.quantile(crowding_final, 0.50)),

            "q75": float(np.quantile(crowding_final, 0.75)),

            "q90": float(np.quantile(crowding_final, 0.90)),

            "q99": float(np.quantile(crowding_final, 0.99)),

            "q100": float(np.quantile(crowding_final, 1.00)),

        },

        "population_scale": "normalized mass; initial total population = 1",

        "time": np.asarray(time_history),

        "total_population_mass": np.asarray(total_population_history),

        "mortality_this_step_mass": np.asarray(deaths_this_step_history),

        "cumulative_total_mortality_mass": np.asarray(

            cumulative_total_deaths_history

        ),

        "cumulative_poison_mortality_mass": np.asarray(

            cumulative_poison_deaths_history

        ),

    }


# ---------------------------------------------------------------------

# Evaluation

# ---------------------------------------------------------------------

def safe_pearson(x: np.ndarray, y: np.ndarray) -> float:

    x = np.asarray(x, dtype=float)

    y = np.asarray(y, dtype=float)

    if len(x) < 2 or np.std(x) <= 0 or np.std(y) <= 0:

        return float("nan")

    return float(np.corrcoef(x, y)[0, 1])


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:

    x = np.asarray(x, dtype=float)

    y = np.asarray(y, dtype=float)

    if len(x) < 2 or np.std(x) <= 0 or np.std(y) <= 0:

        return float("nan")

    return float(spearmanr(x, y).correlation)


def top_fraction_hotspot_overlap(

    modeled_values: np.ndarray,

    observed_values: np.ndarray,

    fraction: float = 0.20,

) -> dict[str, float | int]:

    """

    Compare the highest-risk fraction of locations in modeled and observed

    mortality patterns.

    For fraction=0.20 and n=83 locations:

        k = ceil(0.20 * 83) = 17.

    The overlap fraction is:

        |H_model intersect H_obs| / k.

    Stable descending sorting makes ties reproducible.

    """

    modeled_values = np.asarray(modeled_values, dtype=float)

    observed_values = np.asarray(observed_values, dtype=float)

    if modeled_values.shape != observed_values.shape:

        raise ValueError(

            "modeled_values and observed_values must have the same shape."

        )

    n = int(modeled_values.size)

    if n == 0:

        return {

            "fraction": float(fraction),

            "n_locations": 0,

            "n_hotspots": 0,

            "overlap_count": 0,

            "overlap_fraction": float("nan"),

            "jaccard": float("nan"),

            "random_expected_overlap_fraction": float("nan"),

        }

    if not (0.0 < fraction <= 1.0):

        raise ValueError("fraction must lie in (0, 1].")

    k = max(1, int(np.ceil(fraction * n)))

    model_order = np.argsort(-modeled_values, kind="mergesort")

    obs_order = np.argsort(-observed_values, kind="mergesort")

    model_top = model_order[:k]

    obs_top = obs_order[:k]

    overlap_count = int(

        np.intersect1d(

            model_top,

            obs_top,

            assume_unique=True,

        ).size

    )

    overlap_fraction = float(overlap_count / k)

    union_count = 2 * k - overlap_count

    jaccard = (

        float(overlap_count / union_count)

        if union_count > 0

        else float("nan")

    )

    return {

        "fraction": float(fraction),

        "n_locations": n,

        "n_hotspots": k,

        "overlap_count": overlap_count,

        "overlap_fraction": overlap_fraction,

        "jaccard": jaccard,

        "random_expected_overlap_fraction": float(k / n),

    }


def compare_normalized_spatial_fit(

    modeled_local_deaths: np.ndarray,

    observed_local_deaths: np.ndarray,

) -> tuple[dict[str, float], np.ndarray, np.ndarray]:

    """

    Compare spatial death patterns after scaling both modeled and observed

    local mortality burdens to [0, 1].

    This ignores total death magnitude and asks whether the model predicts

    relatively high-death vs low-death locations correctly.

    """

    modeled_norm = normalize_minmax01(modeled_local_deaths)

    observed_norm = normalize_minmax01(observed_local_deaths)

    residual = modeled_norm - observed_norm

    mse = float(np.mean(residual**2))

    top20 = top_fraction_hotspot_overlap(

        modeled_norm,

        observed_norm,

        fraction=0.20,

    )

    metrics = {

        "pearson_minmax01_spatial_fit": safe_pearson(

            modeled_norm,

            observed_norm,

        ),

        "spearman_minmax01_spatial_fit": safe_spearman(

            modeled_norm,

            observed_norm,

        ),

        "top20_hotspot_overlap_fraction": float(

            top20["overlap_fraction"]

        ),

        "top20_hotspot_overlap_count": int(

            top20["overlap_count"]

        ),

        "top20_hotspot_n_locations": int(

            top20["n_hotspots"]

        ),

        "top20_hotspot_jaccard": float(

            top20["jaccard"]

        ),

        "top20_random_expected_overlap_fraction": float(

            top20["random_expected_overlap_fraction"]

        ),

        "mae_minmax01_spatial_fit": float(

            np.mean(np.abs(residual))

        ),

        "mse_minmax01_spatial_fit": mse,

        "rmse_minmax01_spatial_fit": float(np.sqrt(mse)),

    }

    return metrics, modeled_norm, observed_norm


def observed_sum_near_events(

    event_xy_km: np.ndarray,

    observed_deaths: np.ndarray,

    radius_km: float,

) -> tuple[np.ndarray, np.ndarray]:

    """

    Sum observed deaths inside a radius around each observed event.

    Returns:

      observed_deaths_within_radius

      n_observed_events_within_radius

    """

    tree = cKDTree(event_xy_km)

    neighborhoods = tree.query_ball_point(event_xy_km, r=radius_km)

    observed_sums = np.asarray(

        [float(observed_deaths[idx].sum()) if len(idx) else 0.0 for idx in neighborhoods],

        dtype=float,

    )

    n_events = np.asarray([len(idx) for idx in neighborhoods], dtype=int)

    return observed_sums, n_events


def local_sum_near_events(

    x_km: np.ndarray,

    y_km: np.ndarray,

    grid_values: np.ndarray,

    event_xy_km: np.ndarray,

    radius_km: float,

) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    """

    Sum grid_values inside a radius around each observed event.

    Returns:

      local_sums, n_cells_in_radius, nearest_grid_row, nearest_grid_distance_km

    """

    coords = np.column_stack([x_km, y_km])

    tree = cKDTree(coords)

    neighborhoods = tree.query_ball_point(event_xy_km, r=radius_km)

    local_sums = np.asarray(

        [float(grid_values[idx].sum()) if len(idx) else 0.0 for idx in neighborhoods],

        dtype=float,

    )

    n_cells = np.asarray([len(idx) for idx in neighborhoods], dtype=int)

    nearest_dist, nearest_idx = tree.query(event_xy_km, k=1)

    return local_sums, n_cells, nearest_idx.astype(int), nearest_dist.astype(float)


def compare_modeled_deaths_to_observed(

    modeled_local_deaths: np.ndarray,

    observed_deaths: np.ndarray,

    modeled_total_deaths: float,

) -> dict[str, float]:

    obs_total = float(np.sum(observed_deaths))

    residual = modeled_local_deaths - observed_deaths

    return {

        "observed_total_deaths": obs_total,

        "modeled_total_deaths": float(modeled_total_deaths),

        "modeled_minus_observed_total_deaths": float(modeled_total_deaths - obs_total),

        "modeled_to_observed_total_death_ratio": (

            float(modeled_total_deaths / obs_total) if obs_total > 0 else float("nan")

        ),

        "pearson_local_modeled_vs_observed_deaths": safe_pearson(

            modeled_local_deaths, observed_deaths

        ),

        "spearman_local_modeled_vs_observed_deaths": safe_spearman(

            modeled_local_deaths, observed_deaths

        ),

        "mae_local_modeled_vs_observed_deaths": float(np.mean(np.abs(residual))),

        "rmse_local_modeled_vs_observed_deaths": float(np.sqrt(np.mean(residual**2))),

    }


# ---------------------------------------------------------------------

# Plotting

# ---------------------------------------------------------------------

def maybe_plot_parks(ax, params: Params) -> None:

    if gpd is None:

        return

    shp_path = os.path.join(params.base_dir, params.parks_shp_name)

    if not os.path.exists(shp_path):

        return

    parks = gpd.read_file(shp_path)

    if parks.crs is not None and parks.crs.to_epsg() != 4326:

        parks = parks.to_crs("EPSG:4326")

    parks.boundary.plot(ax=ax, linewidth=0.7)


def format_namibia_axes(ax, params: Params) -> None:

    ax.set_xlim(params.xlim_min, params.xlim_max)

    ax.set_ylim(params.ylim_min, params.ylim_max)

    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("Longitude")

    ax.set_ylabel("Latitude")


def add_observed_points(ax, mort: pd.DataFrame) -> None:

    ax.scatter(

        mort["X"],

        mort["Y"],

        s=np.sqrt(mort["N_Dead"].to_numpy(dtype=float)) * 18,

        facecolors="none",

        edgecolors="black",

        linewidths=0.8,

        zorder=5,

    )

def save_field_map(

    outdir: str,

    grid: pd.DataFrame,

    value_col: str,

    title: str,

    colorbar_label: str,

    filename: str,

    params: Params,

) -> None:

    fig, ax = plt.subplots(figsize=(8.2, 7.2))

    sc = ax.scatter(

        grid["x"],

        grid["y"],

        c=grid[value_col],

        s=16,

        marker="s",

        linewidths=0,

    )

    maybe_plot_parks(ax, params)

    format_namibia_axes(ax, params)

    ax.set_title(title)

    cbar = fig.colorbar(sc, ax=ax, fraction=0.065, pad=0.035)

    cbar.set_label(colorbar_label)

    fig.tight_layout()

    fig.savefig(os.path.join(outdir, filename), dpi=300)

    plt.close(fig)


def save_plots(

    outdir: str,

    grid: pd.DataFrame,

    mort: pd.DataFrame,

    population_death_ts: pd.DataFrame,

    metrics_normalized_spatial: dict[str, float],

    params: Params,

) -> None:

    """Save figures using normalized model mass only."""

    os.makedirs(outdir, exist_ok=True)

    save_field_map(

        outdir=outdir,

        grid=grid,

        value_col="food_model",

        title="Namibia: food resource distribution",

        colorbar_label="Food resource, rescaled to [0, 1]",

        filename="namibia_food_resource_map.png",

        params=params,

    )

    save_field_map(

        outdir=outdir,

        grid=grid,

        value_col="toxin_model",

        title="Namibia: poisoning-risk distribution",

        colorbar_label="Poisoning-risk index, rescaled to [0, 1]",

        filename="namibia_toxin_distribution_map.png",

        params=params,

    )

    fig, ax = plt.subplots(figsize=(7.2, 7.2))

    sc = ax.scatter(

        grid["x"],

        grid["y"],

        c=grid["final_population_mass"],

        s=16,

        marker="s",

        linewidths=0,

    )

    maybe_plot_parks(ax, params)

    add_observed_points(ax, mort)

    format_namibia_axes(ax, params)

    ax.set_title("Namibia: final normalized population mass")

    fig.colorbar(sc, ax=ax, label="Normalized population mass per cell")

    fig.tight_layout()

    fig.savefig(

        os.path.join(outdir, "namibia_final_population_map.png"),

        dpi=300,

    )

    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 7.2))

    sc = ax.scatter(

        grid["x"],

        grid["y"],

        c=grid["cumulative_total_mortality_mass"],

        s=16,

        marker="s",

        linewidths=0,

    )

    maybe_plot_parks(ax, params)

    add_observed_points(ax, mort)

    format_namibia_axes(ax, params)

    ax.set_title("Namibia: cumulative normalized mortality mass")

    fig.colorbar(sc, ax=ax, label="Cumulative mortality mass per cell")

    fig.tight_layout()

    fig.savefig(

        os.path.join(outdir, "namibia_cumulative_total_mortality_map.png"),

        dpi=300,

    )

    plt.close(fig)

    fig, axes = plt.subplots(

        1, 2, figsize=(12.0, 4.8), constrained_layout=True

    )

    axes[0].plot(

        population_death_ts["time"],

        population_death_ts["total_population_mass"],

    )

    axes[0].set_xlabel("Time (days)")

    axes[0].set_ylabel("Normalized population mass")

    axes[0].set_title("Total living population mass")

    axes[1].plot(

        population_death_ts["time"],

        population_death_ts["cumulative_total_mortality_mass"],

        label="Total mortality",

    )

    axes[1].plot(

        population_death_ts["time"],

        population_death_ts["cumulative_poison_mortality_mass"],

        label="Poison-related mortality",

    )

    axes[1].set_xlabel("Time (days)")

    axes[1].set_ylabel("Normalized mortality mass")

    axes[1].set_title("Cumulative mortality mass")

    axes[1].legend()

    fig.savefig(

        os.path.join(outdir, "population_and_mortality_timeseries.png"),

        dpi=300,

    )

    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 5.2))

    ax.scatter(

        mort["modeled_total_mortality_mass_within_radius_minmax01"],

        mort["observed_deaths_within_radius_minmax01"],

        s=28,

        alpha=0.75,

    )

    ax.set_xlabel(

        f"Modeled mortality mass within "

        f"{params.event_compare_radius_km:g} km, scaled 0-1"

    )

    ax.set_ylabel(

        f"Observed mortality within "

        f"{params.event_compare_radius_km:g} km, scaled 0-1"

    )

    ax.set_title("Min-max normalized spatial fit of mortality")

    txt = "\n".join(

        [

            f"Pearson = "

            f"{metrics_normalized_spatial['pearson_minmax01_spatial_fit']:.3f}",

            f"Spearman = "

            f"{metrics_normalized_spatial['spearman_minmax01_spatial_fit']:.3f}",

            f"Top-20% overlap = "

            f"{metrics_normalized_spatial['top20_hotspot_overlap_fraction']:.3f}",

            f"MAE = "

            f"{metrics_normalized_spatial['mae_minmax01_spatial_fit']:.4f}",

            f"MSE = "

            f"{metrics_normalized_spatial['mse_minmax01_spatial_fit']:.4f}",

            f"RMSE = "

            f"{metrics_normalized_spatial['rmse_minmax01_spatial_fit']:.4f}",

        ]

    )

    ax.text(0.05, 0.95, txt, transform=ax.transAxes, va="top")

    fig.tight_layout()

    fig.savefig(

        os.path.join(outdir, "minmax01_spatial_fit_event_radius.png"),

        dpi=300,

    )

    plt.close(fig)


# ---------------------------------------------------------------------

# Main

# ---------------------------------------------------------------------

def main(params: Params | None = None) -> dict:

    params = Params() if params is None else params

    grid_path = os.path.join(params.base_dir, params.grid_rdata_name)

    mortality_path = os.path.join(params.base_dir, params.mortality_csv_name)

    outdir = os.path.join(params.base_dir, params.output_dir_name)

    if params.save_outputs:

        os.makedirs(outdir, exist_ok=True)

    raw_df = load_rdata_frame(grid_path)

    grid = prepare_grid(raw_df, params)

    mort = prepare_mortality(mortality_path)

    # Coordinates in approximate kilometers for graph geometry.

    x_km, y_km, lon0, lat0 = lonlat_to_km(grid["x"].to_numpy(), grid["y"].to_numpy())

    mx_km, my_km, _, _ = lonlat_to_km(

        mort["X"].to_numpy(),

        mort["Y"].to_numpy(),

        lon0=lon0,

        lat0=lat0,

    )

    event_xy_km = np.column_stack([mx_km, my_km])

    observed_deaths = mort["N_Dead"].astype(float).to_numpy()

    observed_deaths_within_radius, n_observed_events_within_radius = observed_sum_near_events(

        event_xy_km=event_xy_km,

        observed_deaths=observed_deaths,

        radius_km=params.event_compare_radius_km,

    )

    # Model fields: uniform resource availability and original spatially
    # heterogeneous toxin exposure. Both influence the model dynamically.

    food = grid["food_model"].to_numpy(dtype=float)

    tau = grid["toxin_model"].to_numpy(dtype=float)

# Initial condition: UD, uniform, or mixture.

    ud0 = grid["ud"].to_numpy(dtype=float)

    ud0 = ud0 / ud0.sum()

    uniform0 = np.ones(len(grid), dtype=float)

    uniform0 = uniform0 / uniform0.sum()

    if params.initial_condition == "ud":

        u0 = ud0

    elif params.initial_condition == "uniform":

        u0 = uniform0;

    elif params.initial_condition == "mixed":

        w = params.uniform_weight

        u0 = (1.0 - w) * ud0 + w * uniform0

        u0 = u0 / u0.sum()

    else:

        raise ValueError(

            f"Unknown initial_condition: {params.initial_condition}. "

            "Use 'ud', 'uniform', or 'mixed'."

        )


    # Movement operators.

    P_diff, Q_food, Q_toxin, movement_meta = prepare_movement_operators(

        x_km,

        y_km,

        food,

        tau,

        params,

    )

    # Local crowding operator:

    # binary radius-neighborhood SUM, not a row-normalized smoother.

    A_local, crowding_meta = build_local_population_matrix(

        x_km,

        y_km,

        params,

    )

    # Run the full hunger-structured spatial model.

    full = simulate_pde(

        u0,

        food,

        tau,

        P_diff,

        Q_food,

        Q_toxin,

        A_local,

        params,

    )

    # Keep all spatial outputs on the normalized population-mass scale.

    initial_population_mass = u0

    final_population_mass = full["final_alive_raw"]

    cumulative_total_mortality_mass = full["total_mortality_raw"]

    cumulative_baseline_hunger_mortality_mass = (

        full["baseline_hunger_mortality_raw"]

    )

    cumulative_poison_mortality_mass = (

        full["poison_excess_mortality_raw"]

    )

    # Compare modeled cumulative total mortality with observed mortality.

    (

        modeled_total_mortality_mass_within_radius,

        n_grid_cells_within_radius,

        nearest_grid_row,

        nearest_grid_distance_km,

    ) = local_sum_near_events(

        x_km=x_km,

        y_km=y_km,

        grid_values=cumulative_total_mortality_mass,

        event_xy_km=event_xy_km,

        radius_km=params.event_compare_radius_km,

    )

    (

        metrics_normalized_spatial_fit,

        modeled_total_mortality_mass_within_radius_minmax01,

        observed_deaths_within_radius_minmax01,

    ) = compare_normalized_spatial_fit(

        modeled_local_deaths=modeled_total_mortality_mass_within_radius,

        observed_local_deaths=observed_deaths_within_radius,

    )

    # Population and mortality time series remain on the normalized mass scale.


    population_death_ts = pd.DataFrame(

        {

            "time": full["time"],

            "total_population_mass": full["total_population_mass"],

            "mortality_this_step_mass": full["mortality_this_step_mass"],

            "cumulative_total_mortality_mass": full["cumulative_total_mortality_mass"],

            "cumulative_poison_mortality_mass": full["cumulative_poison_mortality_mass"],

        }

    )

    # Grid output.

    grid_out = grid.copy()

    grid_out["initial_population_mass"] = initial_population_mass

    grid_out["final_population_mass"] = final_population_mass

    grid_out["final_population_prob"] = full["final_alive_prob"]

    grid_out["final_hunger"] = full["final_hunger"]

    grid_out["cumulative_total_mortality_mass"] = cumulative_total_mortality_mass

    grid_out["cumulative_baseline_hunger_mortality_mass"] = cumulative_baseline_hunger_mortality_mass

    grid_out["cumulative_poison_mortality_mass"] = cumulative_poison_mortality_mass

    # Event-level local comparison.

    mort_out = mort.copy().reset_index(drop=True)

    mort_out["nearest_grid_row"] = nearest_grid_row

    mort_out["nearest_grid_distance_km"] = nearest_grid_distance_km

    mort_out["n_grid_cells_within_radius"] = n_grid_cells_within_radius

    mort_out["comparison_radius_km"] = params.event_compare_radius_km

    # Normalized modeled mortality mass and empirical observed mortality.

    mort_out["modeled_total_mortality_mass_within_radius"] = modeled_total_mortality_mass_within_radius

    mort_out["observed_deaths_at_event"] = observed_deaths

    mort_out["observed_deaths_within_radius"] = observed_deaths_within_radius

    mort_out["n_observed_events_within_radius"] = n_observed_events_within_radius

    # Min-max normalized 0-1 spatial burden.

    mort_out["modeled_total_mortality_mass_within_radius_minmax01"] = (

        modeled_total_mortality_mass_within_radius_minmax01

    )

    mort_out["observed_deaths_within_radius_minmax01"] = (

        observed_deaths_within_radius_minmax01

    )

    # Local covariates at nearest model grid cell.

    mort_out["food_at_nearest_grid_cell"] = food[nearest_grid_row]

    mort_out["toxin_at_nearest_grid_cell"] = tau[nearest_grid_row]

    summary = {

        "model": "Namibia hunger-structured PDE/graph model: uniform low food and heterogeneous toxin",

        "scenario": "uniform_low_food_with_toxin",

        "ablation_definition": f"food_model = {params.uniform_food_level:g} across all cells; heterogeneous toxin and calibrated d_tau retained",

        "toxin_effect_enabled": True,

        "parameters": asdict(params),

        "population_scale": "normalized mass; initial total population = 1",

        "absolute_population_count_used": False,

        "movement_metadata": movement_meta,

        "crowding_metadata": crowding_meta,

        "shared_behavior_mortality_sensitivity": True,

        "mortality_model": "multiplicative",

        "mortality_equation": "M = d1 * (1 + d_v*v) * (1 + d_tau*tau(x))",

        "simulation_horizon_days": float(params.dt * params.n_steps),

        "evaluation_radius_km": float(params.event_compare_radius_km),

        "n_domain_cells": int(len(grid_out)),

        "n_records_used": int(len(mort_out)),

        "initial_population_mass": float(

            population_death_ts["total_population_mass"].iloc[0]

        ),

        "final_population_mass": float(

            population_death_ts["total_population_mass"].iloc[-1]

        ),

        "final_cumulative_total_mortality_mass": float(

            population_death_ts["cumulative_total_mortality_mass"].iloc[-1]

        ),

        "final_cumulative_poison_mortality_mass": float(

            population_death_ts["cumulative_poison_mortality_mass"].iloc[-1]

        ),

        "normalized_spatial_fit": metrics_normalized_spatial_fit,

        "hotspot_evaluation": {

            "fraction": 0.20,

            "overlap_fraction": metrics_normalized_spatial_fit[

                "top20_hotspot_overlap_fraction"

            ],

            "overlap_count": metrics_normalized_spatial_fit[

                "top20_hotspot_overlap_count"

            ],

            "n_hotspots": metrics_normalized_spatial_fit[

                "top20_hotspot_n_locations"

            ],

            "random_expected_overlap_fraction": metrics_normalized_spatial_fit[

                "top20_random_expected_overlap_fraction"

            ],

        },

    }

    if params.save_outputs:

        population_death_ts.to_csv(

            os.path.join(outdir, "population_death_timeseries.csv"),

            index=False,

        )

        grid_out.to_csv(os.path.join(outdir, "pde_mortality_grid.csv"), index=False)

        mort_out.to_csv(os.path.join(outdir, "mortality_event_comparison.csv"), index=False)

        np.savez_compressed(

            os.path.join(outdir, "final_hunger_structured_state.npz"),

            v_grid=full["hunger_grid"],

            u_xv=full["final_structured_density_raw"],

        )

        with open(os.path.join(outdir, "metrics.json"), "w", encoding="utf-8") as f:

            json.dump(summary, f, indent=2)

        if params.make_plots:

            save_plots(

                outdir,

                grid_out,

                mort_out,

                population_death_ts,

                metrics_normalized_spatial_fit,

                params,

            )

    if params.verbose:

        print(json.dumps(summary, indent=2))

    return summary


if __name__ == "__main__":

    main()
