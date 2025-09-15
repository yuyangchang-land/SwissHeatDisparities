# =========================================
# Two-Pass GWR (mgwr): re-run + full diagnostics
#   PASS 1: without socioassis
#   PASS 2: with socioassis (non-NA only)
# =========================================
import geopandas as gpd
import pandas as pd
import numpy as np
from mgwr.gwr import GWR
from mgwr.sel_bw import Sel_BW
from libpysal.weights import DistanceBand
from esda.moran import Moran
import statsmodels.api as sm
import matplotlib.pyplot as plt
import os

# -----------------------------
# 0) Paths
# -----------------------------
shp_path = "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/muni_with_zonal_means.shp"
base_out_dir = "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/results_gwr_multi"
os.makedirs(base_out_dir, exist_ok=True)

# -----------------------------
# 1) Load shapefile & project to meters (LV95 / EPSG:2056)
# -----------------------------
gdf_poly = gpd.read_file(shp_path)
if gdf_poly.crs is None:
    gdf_poly.set_crs(epsg=4326, inplace=True)  # assume WGS84 if missing
gdf_poly = gdf_poly.to_crs(epsg=2056)

# Representative points for GWR
gdf_pts = gdf_poly.copy()
gdf_pts["geometry"] = gdf_pts.representative_point()

# -----------------------------
# 2) Helpers
# -----------------------------
def z(x: pd.Series) -> pd.Series:
    m = np.nanmean(x)
    s = np.nanstd(x)
    return (x - m) / s if s not in (0, np.nan) else x * 0

def mm01(x_arr: np.ndarray) -> np.ndarray:
    rmin, rmax = np.nanmin(x_arr), np.nanmax(x_arr)
    return (x_arr - rmin) / (rmax - rmin) if rmax != rmin else np.zeros_like(x_arr)

def zscore_cols(df, cols):
    df = df.copy()
    for c in cols:
        if c in df.columns:
            df[c] = z(pd.to_numeric(df[c], errors="coerce"))
    return df

def ensure_numeric(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

# -----------------------------
# 3) Column name harmonization (use your printed names)
# -----------------------------
alias_map = {
    # CHEI components
    "lst_urban": "lst_rbn",   # LST
    "heat_d":    "avgt25_",   # heat days
    "heat_p":    "htwv_p_",   # heatwave probability

    # physical vars
    "dem_mean":  "dem_men",
    "imp_mean":  "imprv_m",
    "evi_mean":  "evi_men",
    "bldg_mean": "bhght_m",
    "albd_mean": "albd_mn",

    # socio vars
    "foreign19": "forgn19",
    "socioassis":"socisss",
    "ave_living":"av_lvng",
    "num_reside":"num_rsd",

    # others used directly
    "income": "income",
    "age65_": "age65_",
    "age80_": "age80_",
    "for_eu": "for_eu"
}
for std, real in alias_map.items():
    if std not in gdf_pts.columns:
        if real in gdf_pts.columns:
            gdf_pts[std] = pd.to_numeric(gdf_pts[real], errors="coerce")
        else:
            raise KeyError(f"Expected '{real}' for {std}, but not found.")

# -----------------------------
# 4) Derived fields (using standard names)
# -----------------------------
to_numeric = [
    "lst_urban","heat_d","heat_p",
    "dem_mean","imp_mean","evi_mean","bldg_mean","albd_mean",
    "income","age65_","age80_","foreign19","for_eu",
    "socioassis","ave_living","num_reside"
]
gdf_pts = ensure_numeric(gdf_pts, to_numeric)

gdf_pts["for_eu_prop"]  = np.where(gdf_pts["for_eu"] > 1, gdf_pts["for_eu"] / 100, gdf_pts["for_eu"])
gdf_pts["age65_79"]     = gdf_pts["age65_"] - gdf_pts["age80_"]
gdf_pts["percen_eu"]    = gdf_pts["foreign19"] * gdf_pts["for_eu_prop"]
gdf_pts["percen_noneu"] = gdf_pts["foreign19"] * (1 - gdf_pts["for_eu_prop"])
gdf_pts["w"]            = np.log(np.maximum(gdf_pts["num_reside"], 0) + 1)
gdf_pts["income_ln"]    = np.log(np.maximum(gdf_pts["income"], 0) + 1)

# CHEI components (z-score of raw indicators)
gdf_pts["LST_z"] = z(gdf_pts["lst_urban"])
gdf_pts["HWD_z"] = z(gdf_pts["heat_d"])
gdf_pts["HWP_z"] = z(gdf_pts["heat_p"])
# Equal-weight composite as in your main spec
gdf_pts["CHEI_pca"] = mm01(
    np.dot(gdf_pts[["LST_z","HWD_z","HWP_z"]], np.repeat(1/3, 3))
)

# -----------------------------
# 5) Define predictors for each pass
# -----------------------------
socio_vars_base = ["income_ln","age65_79","age80_","percen_eu","percen_noneu","ave_living"]
physical_vars   = ["dem_mean","evi_mean","bldg_mean","albd_mean","imp_mean"]

# PASS 1: exclude socioassis
rhs_pass1 = socio_vars_base + physical_vars
# PASS 2: include socioassis (drop NA socioassis rows)
rhs_pass2 = socio_vars_base + ["socioassis"] + physical_vars

# -----------------------------
# 6) Dependent variables
# -----------------------------
dep_vars = ["CHEI_pca", "LST_z", "HWD_z", "HWP_z"]

# -----------------------------
# 7) GWR runner with full diagnostics
# -----------------------------
def run_gwr_block(df_in, dep_vars, rhs_vars, tag):
    """Fit GWR for each DV; save results under base_out_dir/tag_<DV>/..., plus rich diagnostics."""
    diag_rows = []

    for dv in dep_vars:
        print(f"\n=== [{tag}] Running GWR for {dv} ===")
        out_dir = os.path.join(base_out_dir, f"{tag}_{dv}")
        os.makedirs(out_dir, exist_ok=True)

        # Drop NA rows for this DV + predictors
        cols_needed = [dv] + rhs_vars
        df_model = df_in.dropna(subset=cols_needed).copy()
        if df_model.empty:
            print(f"[{tag}:{dv}] No data after dropping NAs. Skipping.")
            continue

        # Scale predictors (already z-scored outside if desired; safe to re-z)
        df_model = zscore_cols(df_model, rhs_vars)

        coords = np.array([(geom.x, geom.y) for geom in df_model.geometry])
        y = df_model[[dv]].to_numpy()
        X = df_model[rhs_vars].to_numpy()

        # Bandwidth
        selector = Sel_BW(coords, y, X, kernel="bisquare", fixed=False, spherical=False)
        bw = selector.search()
        print(f"[{tag}:{dv}] Bandwidth (neighbors): {bw}")

        # Fit
        gwr = GWR(coords, y, X, bw, kernel="bisquare", fixed=False, spherical=False)
        res = gwr.fit()

        # Predictions & residuals
        yhat = np.asarray(res.predy).ravel()
        resid = (y.ravel() - yhat)
        rmse = float(np.sqrt(np.nanmean(resid**2)))
        mae  = float(np.nanmean(np.abs(resid)))

        # Global OLS baseline (for AICc comparison, same sample)
        ols = sm.OLS(y, sm.add_constant(X)).fit()
        n, k = X.shape
        AICc_global = float(ols.aic + (2 * k * (k + 1)) / max(n - k - 1, 1))

        # Spatial residual autocorrelation
        w = DistanceBand(coords, threshold=50000, binary=True, silence_warnings=True)
        mi = Moran(resid, w)

        # Local results
        param_names = ["Intercept"] + rhs_vars
        params = pd.DataFrame(res.params, columns=param_names)
        tvals  = pd.DataFrame(res.tvalues, columns=[f"{v}_t" for v in param_names])
        localR2 = pd.Series(np.asarray(res.localR2).ravel(), name="localR2")

        out = pd.concat(
            [df_model.reset_index(drop=True)[["geometry"]], params, tvals, localR2,
             pd.Series(resid, name="gw_residual")],
            axis=1
        )
        out_gdf = gpd.GeoDataFrame(out, geometry="geometry", crs="EPSG:2056")

        # Diagnostics row
        diag = {
            "tag": tag,
            "dv": dv,
            "n_obs": int(len(df_model)),
            "k_predictors": int(k),
            "bandwidth_neighbors": int(bw),
            # Fit quality
            "pseudoR2_GWR": float(res.R2),
            "AIC_GWR": float(res.aic),
            "AICc_GWR": float(res.aicc),
            "AICc_global_OLS": AICc_global,
            "AICc_diff_OLS_minus_GWR": float(AICc_global - res.aicc),
            "RMSE_GWR": rmse,
            "MAE_GWR": mae,
            # Spatial autocorrelation of residuals
            "MoranI_resid": float(mi.I),
            "MoranI_p_norm": float(mi.p_norm),
            # Local R2 distribution
            "localR2_q25": float(np.nanquantile(localR2, 0.25)),
            "localR2_median": float(np.nanmedian(localR2)),
            "localR2_q75": float(np.nanquantile(localR2, 0.75)),
        }

        # Save diagnostics (per DV)
        pd.DataFrame([diag]).to_csv(os.path.join(out_dir, f"{dv}_diagnostics.csv"), index=False)
        print(diag)

        # Save local outputs
        out_gdf.to_file(os.path.join(out_dir, f"gwr_results_{dv}.gpkg"),
                        layer=f"gwr_local_results_{dv}", driver="GPKG")
        out_gdf.drop(columns="geometry").to_csv(os.path.join(out_dir, f"{dv}_gwr_results.csv"), index=False)

        # DV summary & plots (optional)
        desc = df_model[[dv]].describe().T
        desc.to_csv(os.path.join(out_dir, f"{dv}_descriptives.csv"))

        plt.figure(figsize=(6,4))
        pd.Series(df_model[dv]).plot(kind="hist", bins=40)
        plt.title(f"{dv} distribution [{tag}]")
        plt.xlabel(dv); plt.ylabel("Count"); plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{dv}_histogram.png"), dpi=300)
        plt.close()

        df_map = gdf_poly.reset_index(drop=True).join(df_model[[dv]].reset_index(drop=True))
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        df_map.plot(column=dv, cmap="viridis", legend=True,
                    legend_kwds={"label": f"{dv} [{tag}]"}, ax=ax,
                    edgecolor="black", linewidth=0.2)
        ax.set_title(f"{dv} map [{tag}]", fontsize=14)
        ax.axis("off")
        fig.savefig(os.path.join(out_dir, f"{dv}_map.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

        diag_rows.append(diag)

    return pd.DataFrame(diag_rows)

# -----------------------------
# 8) PASS 1 — exclude socioassis
# -----------------------------
df_pass1 = gdf_pts.copy()
df_pass1 = zscore_cols(df_pass1, rhs_pass1)  # z-score predictors
diag1 = run_gwr_block(df_pass1, dep_vars, rhs_pass1, tag="no_socio")

# -----------------------------
# 9) PASS 2 — include socioassis (drop NA socioassis)
# -----------------------------
df_pass2 = gdf_pts.dropna(subset=["socioassis"]).copy()
df_pass2 = zscore_cols(df_pass2, rhs_pass2)  # z-score predictors (incl socioassis)
diag2 = run_gwr_block(df_pass2, dep_vars, rhs_pass2, tag="with_socio")

# -----------------------------
# 10) Master diagnostics table
# -----------------------------
diag_all = pd.concat([d for d in [diag1, diag2] if d is not None], ignore_index=True)
diag_all = diag_all.sort_values(["tag","dv"]).reset_index(drop=True)
master_csv = os.path.join(base_out_dir, "gwr_diagnostics_master.csv")
diag_all.to_csv(master_csv, index=False)
print(f"\nAll GWR runs complete.\n- Per-DV outputs in: {base_out_dir}/<tag>_<DV>/\n- Master diagnostics: {master_csv}\n")