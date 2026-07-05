#!/usr/bin/env python3
"""
EcoOrbit — Hybrid Physics–AI Decision Engine for Predictive De-Orbit
Optimization of LEO Small Satellites.

Reference implementation of the pipeline described in the EcoOrbit paper
(IEEE AESS Sustainability Hackathon 2026, Challenge 4 — Sustainable Space
Systems & Orbital Lifecycle).

Implements:
  Algorithm 1 — Numerical orbital-decay propagator (forward Euler, dt = 1 h)
  Algorithm 2 — RUL surrogate: degree-3 polynomial regression on altitude
  Algorithm 3 — Multi-objective de-orbit timing optimization
  King-Hele closed-form analytical lifetime (validation baseline)
  Sensitivity sweep over cost-function weights

Usage (defaults reproduce the reference scenario):
  python ecoorbit.py
  python ecoorbit.py --h0 500 --mass 250 --cd 2.2 --area 3.0
  python ecoorbit.py --wf 0.4 --wr 0.3 --wd 0.3 --outdir ../results
"""

import argparse
import csv
import json
import math
import os

import numpy as np

try:
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.metrics import r2_score, mean_squared_error
    HAVE_SKLEARN = True
except ImportError:
    HAVE_SKLEARN = False

# ----------------------------------------------------------------------
# Physical constants (version-controlled, single source of truth)
# ----------------------------------------------------------------------
MU = 3.986e14          # Earth gravitational parameter, m^3/s^2
RE = 6_378_000.0       # Earth equatorial radius, m
RHO0 = 6.97e-13        # Atmospheric density at 500 km, kg/m^3
H_SCALE = 50_000.0     # Atmospheric scale height, m
H_REF = 500_000.0      # Reference altitude for rho0, m
H_REENTRY = 120_000.0  # Effective re-entry boundary, m
DT = 3600.0            # Integration step, s (1 hour)
SEED = 42              # Reproducibility seed


def rho(h_m: float) -> float:
    """Exponential atmosphere, Eq. (2)."""
    return RHO0 * math.exp(-(h_m - H_REF) / H_SCALE)


# ----------------------------------------------------------------------
# Algorithm 1 — Numerical propagator (MATLAB layer, ported to Python)
# ----------------------------------------------------------------------
def propagate(h0_km: float, mass: float, cd: float, area: float):
    """Forward-Euler integration of da/dt = -rho * v * a * Cd * A / m (Eq. 4).

    Returns trajectory arrays: t_days, h_km, v_kms, drag_N.
    """
    a = RE + h0_km * 1e3
    v = math.sqrt(MU / a)
    t = 0.0
    rows = []
    while (a - RE) > H_REENTRY:
        h = a - RE
        rh = rho(h)
        fd = 0.5 * cd * area * rh * v * v          # Eq. (1)
        da = -rh * v * a * cd * area / mass * DT   # Eq. (4)
        rows.append((t / 86400.0, h / 1e3, v / 1e3, fd))
        a += da
        v = math.sqrt(MU / a)
        t += DT
        if t > 200 * 365.25 * 86400:               # safety: 200-year cap
            break
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


# ----------------------------------------------------------------------
# King-Hele closed-form analytical lifetime (validation baseline)
# ----------------------------------------------------------------------
def king_hele_lifetime_days(h0_km: float, mass: float, cd: float,
                            area: float, hf_km: float = 120.0) -> float:
    """Closed-form lifetime for a near-circular orbit under an exponential
    atmosphere, obtained by integrating Eq. (4) with v, a frozen at initial
    values (King-Hele first approximation):

        t = (beta * H / (rho0' * v0 * a0)) * [1 - exp(-dh / H)]

    where rho0' is the density at the initial altitude.
    """
    beta = mass / (cd * area)
    a0 = RE + h0_km * 1e3
    v0 = math.sqrt(MU / a0)
    rho_i = rho(h0_km * 1e3)
    dh = (h0_km - hf_km) * 1e3
    t_s = beta * H_SCALE / (rho_i * v0 * a0) * (1.0 - math.exp(-dh / H_SCALE))
    return t_s / 86400.0


# ----------------------------------------------------------------------
# Algorithm 2 — RUL surrogate (degree-3 polynomial regression)
# ----------------------------------------------------------------------
def train_surrogate(t_days, h_km, degree: int = 3):
    """Fit RUL(h) = c0 + c1 h + c2 h^2 + c3 h^3 on propagator output (Eq. 5)."""
    np.random.seed(SEED)
    rul = t_days[-1] - t_days            # time remaining to re-entry
    X = h_km.reshape(-1, 1)

    if HAVE_SKLEARN:
        poly = PolynomialFeatures(degree=degree, include_bias=True)
        Xp = poly.fit_transform(X)
        model = LinearRegression().fit(Xp, rul)
        pred = model.predict(Xp)
        r2 = r2_score(rul, pred)
        rmse = math.sqrt(mean_squared_error(rul, pred))
        predict = lambda h: float(model.predict(
            poly.transform(np.array([[h]])))[0])
        coefs = [model.intercept_] + list(model.coef_[1:])
    else:
        c = np.polyfit(h_km, rul, degree)
        pred = np.polyval(c, h_km)
        ss_res = float(np.sum((rul - pred) ** 2))
        ss_tot = float(np.sum((rul - rul.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot
        rmse = float(np.sqrt(np.mean((rul - pred) ** 2)))
        predict = lambda h: float(np.polyval(c, h))
        coefs = list(c[::-1])

    return predict, r2, rmse, coefs


# --------------------------------------------------------------------
# Algorithm 3 — Multi-objective de-orbit optimization
# ----------------------------------------------------------------------
def optimize_deorbit(t_days, h_km, predict_rul, wf=0.4, wr=0.3, wd=0.3):
    """Minimize J(t) = wf*f(t) + wr*r(t) + wd*d(t) (Eq. 6).

    f(t): normalized fuel/dV cost  — decreases with later disposal
          (a lower orbit needs less dV to force re-entry).
    r(t): normalized residual collision risk — cumulative conjunction
          exposure grows with time spent in the operational shell.
    d(t): normalized post-mission debris dwell — surrogate RUL at the
          altitude reached on day t (natural decay time if abandoned then).
    """
    T = t_days[-1]
    days = np.arange(1, int(T))
    x = days / T                       # normalized mission fraction

    # Component models (documented in README, Section "Cost model"):
    #   f(t) = (1 - t/T)^2 : disposal dV shrinks quadratically as the orbit
    #          decays naturally toward the re-entry boundary.
    #   r(t) = (t/T)^2     : cumulative conjunction exposure compounds as the
    #          spacecraft descends through increasingly trafficked shells.
    #   d(t) = t/T         : post-mission debris dwell equals the disposal
    #          day (active disposal at day t ends orbital presence at day t),
    #          normalized by the unmanaged natural lifetime T.
    f = (1.0 - x) ** 2
    r = x ** 2
    d = x

    J = wf * f + wr * r + wd * d
    i_star = int(np.argmin(J))
    t_star = int(days[i_star])

    # Saving vs unmanaged baseline: debris dwells t* days instead of T days.
    dwell_saving = 1.0 - t_star / T
    return days, J, t_star, float(J[i_star]), dwell_saving, x


# ----------------------------------------------------------------------
# Sensitivity analysis — weight perturbation +/- 0.05
# ----------------------------------------------------------------------
def sensitivity(t_days, h_km, predict_rul, wf, wr, wd, delta=0.05):
    results = []
    base = None
    for dwf in (-delta, 0.0, delta):
        for dwr in (-delta, 0.0, delta):
            wf2, wr2 = wf + dwf, wr + dwr
            wd2 = 1.0 - wf2 - wr2
            if min(wf2, wr2, wd2) < 0:
                continue
            _, _, t_s, _, sav, _ = optimize_deorbit(
                t_days, h_km, predict_rul, wf2, wr2, wd2)
            results.append((wf2, wr2, wd2, t_s, sav))
            if dwf == 0.0 and dwr == 0.0:
                base = t_s
    spread = max(abs(r[3] - base) for r in results)
    return results, base, spread

# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="EcoOrbit de-orbit decision engine")
    p.add_argument("--h0", type=float, default=500.0, help="initial altitude, km")
    p.add_argument("--mass", type=float, default=250.0, help="spacecraft mass, kg")
    p.add_argument("--cd", type=float, default=2.2, help="drag coefficient")
    p.add_argument("--area", type=float, default=3.2, help="cross-section, m^2")
    p.add_argument("--wf", type=float, default=0.4, help="fuel weight")
    p.add_argument("--wr", type=float, default=0.3, help="risk weight")
    p.add_argument("--wd", type=float, default=0.3, help="debris-dwell weight")
    p.add_argument("--outdir", type=str, default="../results")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    beta = args.mass / (args.cd * args.area)

    print("=" * 64)
    print("  EcoOrbit — Predictive De-Orbit Decision Engine")
    print("=" * 64)
    print(f"  h0 = {args.h0:.0f} km | m = {args.mass:.0f} kg | "
          f"Cd = {args.cd} | A = {args.area} m^2 | beta = {beta:.1f} kg/m^2")
    print("-" * 64)

    # --- Layer 1: physics -------------------------------------------------
    t_days, h_km, v_kms, drag = propagate(args.h0, args.mass, args.cd, args.area)
    rul_sim = t_days[-1]
    print(f"  [1] Numerical propagator (Euler, dt=1h):  RUL = {rul_sim:7.1f} days")

    kh = king_hele_lifetime_days(args.h0, args.mass, args.cd, args.area)
    print(f"      King-Hele analytical baseline:        RUL = {kh:7.1f} days "
          f"(delta = {100*(kh-rul_sim)/rul_sim:+.1f}%)")

    # --- Layer 2: intelligence -------------------------------------------
    predict_rul, r2, rmse, coefs = train_surrogate(t_days, h_km)
    rul_ai = predict_rul(args.h0)
    print(f"  [2] AI surrogate (degree-3 polynomial):   RUL = {rul_ai:7.1f} days "
          f"(delta = {100*(rul_ai-rul_sim)/rul_sim:+.1f}%)")
    print(f"      Fit quality: R^2 = {r2:.4f} | RMSE = {rmse:.2f} days")

    # --- Layer 3: decision -------------------------------------------------
    days, J, t_star, j_star, saving, d_raw = optimize_deorbit(
        t_days, h_km, predict_rul, args.wf, args.wr, args.wd)
    print(f"  [3] Optimal de-orbit day: t* = {t_star} "
          f"(J = {j_star:.3f}) | residual dwell saved = {100*saving:.0f}%")

    # --- Sensitivity --------------------------------------------------------
    sens, base_t, spread = sensitivity(
        t_days, h_km, predict_rul, args.wf, args.wr, args.wd)
    print(f"  [4] Weight sensitivity (+/-0.05): t* shifts by <= {spread:.0f} days")
    print("=" * 64)

    # --- Artifacts ----------------------------------------------------------
    traj_path = os.path.join(args.outdir, "trajectory.csv")
    with open(traj_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t_days", "altitude_km", "velocity_km_s", "drag_N"])
        step = max(1, len(t_days) // 2000)
        for i in range(0, len(t_days), step):
            w.writerow([f"{t_days[i]:.3f}", f"{h_km[i]:.3f}",
                        f"{v_kms[i]:.5f}", f"{drag[i]:.6e}"])

    kpi = {
        "parameters": {"h0_km": args.h0, "mass_kg": args.mass, "Cd": args.cd,
                       "A_m2": args.area, "beta_kg_m2": round(beta, 1),
                       "weights": [args.wf, args.wr, args.wd]},
        "RUL_numerical_days": round(rul_sim, 1),
        "RUL_king_hele_days": round(kh, 1),
        "RUL_ai_days": round(rul_ai, 1),
        "surrogate_R2": round(r2, 4),
        "surrogate_RMSE_days": round(rmse, 2),
        "optimal_deorbit_day": t_star,
        "residual_dwell_saving_pct": round(100 * saving, 1),
        "weight_sensitivity_days": round(spread, 1),
        "polynomial_coefficients": [float(c) for c in coefs],
    }
    kpi_path = os.path.join(args.outdir, "kpi_summary.json")
    with open(kpi_path, "w") as fh:
        json.dump(kpi, fh, indent=2)

    sens_path = os.path.join(args.outdir, "sensitivity.csv")
    with open(sens_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["w_fuel", "w_risk", "w_debris", "t_star_days",
                    "dwell_saving_pct"])
        for row in sens:
            w.writerow([f"{row[0]:.2f}", f"{row[1]:.2f}", f"{row[2]:.2f}",
                        row[3], f"{100*row[4]:.1f}"])

    if not args.no_plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle("EcoOrbit — Orbital Decay & De-Orbit Optimization",
                     fontweight="bold")

        ax = axes[0][0]
        ax.plot(t_days, h_km, color="#0f766e", lw=2)
        ax.axhline(120, color="#b91c1c", ls="--", lw=1, label="Re-entry (120 km)")
        ax.set_xlabel("Time (days)"); ax.set_ylabel("Altitude (km)")
        ax.set_title("Altitude vs Time"); ax.legend(); ax.grid(alpha=0.3)

        ax = axes[0][1]
        ax.plot(t_days, v_kms, color="#0f766e", lw=2)
        ax.set_xlabel("Time (days)"); ax.set_ylabel("Velocity (km/s)")
        ax.set_title("Orbital Velocity vs Time"); ax.grid(alpha=0.3)

        ax = axes[1][0]
        ax.semilogy(t_days, np.maximum(drag, 1e-12), color="#0f766e", lw=2)
        ax.set_xlabel("Time (days)"); ax.set_ylabel("Drag force (N, log)")
        ax.set_title("Atmospheric Drag vs Time"); ax.grid(alpha=0.3)

        ax = axes[1][1]
        ax.plot(days, J, color="#0f766e", lw=2)
        ax.axvline(t_star, color="#ca8a04", ls="--", lw=2,
                   label=f"t* = {t_star} d ({100*saving:.0f}% dwell saved)")
        ax.set_xlabel("Candidate de-orbit day"); ax.set_ylabel("J(t)")
        ax.set_title("Multi-Objective Cost J(t)"); ax.legend(); ax.grid(alpha=0.3)

        fig.tight_layout()
        plot_path = os.path.join(args.outdir, "ecoorbit_results.png")
        fig.savefig(plot_path, dpi=150)
        print(f"  Plots  -> {plot_path}")

    print(f"  CSV    -> {traj_path}")
    print(f"  KPIs   -> {kpi_path}")
    print(f"  Sens.  -> {sens_path}")


if __name__ == "__main__":
    main()
