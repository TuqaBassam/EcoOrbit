# EcoOrbit 🛰️

**A Hybrid Physics–AI Decision Engine for Predictive De-Orbit Optimization of LEO Small Satellites**

IEEE AESS Sustainability Hackathon 2026 — Challenge 4: Sustainable Space Systems & Orbital Lifecycle
Team: EcoOrbit — AESS Egypt · AESS Jordan · AESS Tunisia · AESS | ESSTHS

> *"Predicting Re-Entry. Protecting Orbit. Preserving Space."*

---

## What EcoOrbit does

End-of-life satellite disposal is still governed by static rules (IADC 25-year guideline, FCC 5-year rule) rather than optimization. EcoOrbit converts de-orbit **timing** into a predictive, multi-objective decision:

1. **Physics layer** — a numerical orbital-decay propagator (forward Euler, Δt = 1 h) integrates quadratic atmospheric drag from 500 km to the 120 km re-entry boundary.
2. **Intelligence layer** — a degree-3 polynomial regression trained on the propagator output predicts Remaining Useful Life (RUL) in **microseconds** instead of seconds (>10⁶× speed-up), enabling fleet-scale re-evaluation.
3. **Decision layer** — a weighted cost function J(t) = w_f·f(t) + w_r·r(t) + w_d·d(t) balances fuel, collision risk, and post-mission debris dwell to select the optimal disposal day t*.

**Headline result (reference 250 kg smallsat, 500 km):** natural decay ≈ 565 days; optimal disposal at **t\* ≈ 200 days** removes **~64 %** of post-mission debris dwell time, validated by three independent estimators agreeing within a few percent.

---

## Repository layout

```
EcoOrbit/
├── code/
│   ├── ecoorbit.py              # Full pipeline: propagator + AI surrogate + optimizer + sensitivity
│   ├── ecoorbit_propagator.m    # MATLAB physics layer (exports trajectory.csv for Python)
│   └── requirements.txt
├── results/
│   ├── ecoorbit_results.png     # Altitude / velocity / drag / J(t) plots
│   ├── trajectory.csv           # Propagated [t, h, v, F_D] time series
│   ├── kpi_summary.json         # All KPIs from the latest run
│   ├── sensitivity.csv          # Weight-perturbation sweep results
│   ├── trajectory_matlab.csv    # Trajectory exported by the MATLAB layer (cross-check)
│   ├── matlab_decay_plots.png   # Figure produced by the MATLAB layer
│   └── matlab_run_output.txt    # Captured MATLAB/Octave run log
├── docs/
│   ├── EcoOrbit_IEEE_Report.pdf # Full technical paper
│   ├── EcoOrbit_Evidence_Sheet.pdf
│   └── AI_External_Resource_Disclosure.md
├── EcoOrbit_Presentation.pdf
└── README.md
```

## Setup

Requires Python ≥ 3.9. From the repository root:

```bash
pip install numpy matplotlib scikit-learn
```

(scikit-learn is optional — the pipeline falls back to `numpy.polyfit`, which is mathematically equivalent ordinary least squares.)

The MATLAB propagator (`ecoorbit_propagator.m`) runs on MATLAB R2021a+ or GNU Octave and requires no toolboxes.

## Execution / reproduction steps

**Full pipeline (one command, < 10 s on a laptop):**

```bash
cd code
python ecoorbit.py
```

This reproduces the reference scenario end-to-end and writes plots, `trajectory.csv`, `kpi_summary.json`, and `sensitivity.csv` into `results/`.

**Live parameter exploration** (this is our live-demo mode — any parameter can be changed and the whole decision re-computes in seconds):

```bash
python ecoorbit.py --h0 450 --mass 150 --area 1.5      # different spacecraft
python ecoorbit.py --wf 0.5 --wr 0.25 --wd 0.25        # different operator policy
```

**MATLAB physics layer** (optional, mirrors the Python propagator):

```matlab
>> ecoorbit_propagator     % prints RUL, opens 4-panel decay figure, exports trajectory.csv + PNG
```

Verified on GNU Octave 8.4 (no toolboxes): prints `Numerical RUL: 565.0 days` and `King-Hele analytical: 562.8 days (delta = -0.4%)` — exact agreement with the Python layer (log: `results/matlab_run_output.txt`).

## Parameters and assumptions

| Symbol | Value | Description | Source |
|---|---|---|---|
| m | 250 kg | Spacecraft mass | Representative EO smallsat class |
| C_D | 2.2 | Drag coefficient | Vallado & McClain (2013), standard LEO value |
| A | 3.2 m² | Effective cross-section (tumbling, incl. panels) | β = m/(C_D·A) = 35.5 kg/m², typical smallsat range 20–100 |
| h₀ | 500 km | Initial circular altitude | Reference scenario |
| ρ₀ | 6.97×10⁻¹³ kg/m³ | Density at 500 km | U.S. Standard Atmosphere / NRLMSISE-00 mid-cycle mean |
| H | 50 km | Scale height | Representative LEO value |
| μ | 3.986×10¹⁴ m³/s² | Earth gravitational parameter | Standard |
| h_re | 120 km | Effective re-entry boundary | Standard convention |
| Δt | 3600 s | Euler integration step | Transparency/accuracy trade-off (see paper §V) |
| w_f, w_r, w_d | 0.4 / 0.3 / 0.3 | Cost weights | Operator-tunable; defaults per paper §VII |

**Key assumptions:** near-circular orbit; exponential atmosphere (no diurnal bulge, F10.7 solar-cycle, or Ap geomagnetic modulation); secular drag-only decay (no J2 coupling to lifetime); disposal maneuver treated as prompt re-entry. Limitations and the mitigation roadmap (RK4, NRLMSISE-00, multi-class training, Monte-Carlo UQ) are detailed in the paper, §X.

## Cost model (decision layer)

With x = t/T (T = natural lifetime):

- **Fuel** f(x) = (1−x)² — disposal ΔV shrinks as the orbit decays naturally toward re-entry.
- **Risk** r(x) = x² — cumulative conjunction exposure compounds with residual time in trafficked shells.
- **Dwell** d(x) = x — disposing on day t ends orbital presence at day t (vs. T days unmanaged).

J(x) = 0.4(1−x)² + 0.3x² + 0.3x has a unique interior minimum at x* ≈ 0.357 → **t\* ≈ 200 days, J\* ≈ 0.31**, i.e. ~64 % less debris dwell than the unmanaged baseline.

## Baseline comparison

| Strategy | Debris dwell after mission end | Fuel | Notes |
|---|---|---|---|
| **Unmanaged (baseline)** | ~565 days | 0 | Natural decay only; full conjunction exposure |
| IADC 25-yr rule | ≤ 25 yr (already met here) | 0 | Compliance ≠ optimality |
| **EcoOrbit t\* = 200 d** | ~200 days (−64 %) | Modest ΔV | Optimum of J(t); auditable and operator-tunable |

## Validation

Three independent estimators of natural orbital lifetime:

| Method | RUL | Δ vs. numerical |
|---|---|---|
| Numerical propagator (Euler, Δt = 1 h) | 565 days | — |
| King-Hele closed-form analytical | 563 days | −0.4 % |
| AI polynomial surrogate (R² = 0.994, RMSE ≈ 13 d) | 542 days | −4.2 % |

Sensitivity: perturbing each weight by ±0.05 (with Σw = 1) shifts t\* by ≤ ~60 days without changing the qualitative recommendation (dispose early-mid life, not at natural decay). Full sweep in `results/sensitivity.csv`.

## KPIs (primary)

- **Technical KPI:** RUL prediction — R² = 0.994, RMSE ≈ 13 days, inference ~1 µs/query (>10⁶× speed-up vs. 4 s numerical run); 10⁴-satellite fleet re-evaluated in < 1 s.
- **Sustainability KPI:** post-mission debris dwell time reduced ~64 % vs. unmanaged baseline (565 → ~200 days per managed satellite).

## End-of-life & long-term impact

Applied across an operator fleet of 10² spacecraft, a ~64 % dwell reduction removes hundreds of satellite-years of collision exposure; at 10⁴ constellation scale the cumulative reduction becomes significant for Kessler-cascade prevention. The surrogate's audit-friendly closed form (4 coefficients, kilobytes, <1 W·s per 10⁶ queries) means the sustainability tool itself carries a negligible footprint. See paper §XI for SDG 9/12/13 mapping.

## Team

Tuqa Bassam (Chair, IEEE AESS Student Chapter) , Ahmed Amer (vice chair,IEEE AESS Student Chapter), Ahmed Mamdouh ( secrtary,IEEE AESS Student Chapter) — AESS Egypt · AESS Jordan · AESS Tunisia · AESS | ESSTHS.

## Disclosure

AI assistance, libraries, data sources, and standards used are declared in `docs/AI_External_Resource_Disclosure.md`.
