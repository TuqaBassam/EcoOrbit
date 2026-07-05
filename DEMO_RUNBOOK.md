# AI / External Resource Disclosure — Team EcoOrbit

IEEE AESS Sustainability Hackathon 2026 — Challenge 4 Final Phase

Per Section 2 of the Finalist Guide, we disclose all meaningful AI assistance, libraries, datasets, prior work, and external technical sources. Every item below can be explained line-by-line by registered team members.

## 1. AI assistance

| Tool | Used for | Not used for |
|---|---|---|
| Claude (Anthropic) | Repository packaging: refactoring our pipeline into a documented, parameterized Python reference implementation (`code/ecoorbit.py`) reproducing the algorithms specified in our paper (§V–§VII, Algorithms 1–3); drafting README, evidence sheet, and presentation-support documents; cross-checking that published headline figures reproduce from the stated equations and constants (this identified and corrected a cross-section units typo, A = 3.2 m², β = 35.5 kg/m²) | The scientific concept, system architecture, equations, algorithm design, cost-function formulation, and the technical paper — these are the team's original work |

The physics (Eqs. 1–4), the hybrid MATLAB→Python architecture, the RUL surrogate approach, and the multi-objective J(t) formulation predate and are independent of any AI assistance.

## 2. Software libraries

| Library | Version | Role | License |
|---|---|---|---|
| NumPy | ≥1.24 | Numerics, OLS fallback | BSD-3 |
| scikit-learn | ≥1.3 | PolynomialFeatures + LinearRegression (Eq. 5) | BSD-3 |
| matplotlib | ≥3.7 | Result plots | PSF-based |
| MATLAB R2024a / GNU Octave | — | Physics-layer propagator | Commercial / GPL |

## 3. Datasets

No external datasets. Training data for the AI surrogate is generated internally by our own deterministic propagator (paper §VI-A). Atmospheric reference density is a published constant (see §4).

## 4. Technical sources, standards, and datasheets

- King-Hele, *Theory of Satellite Orbits in an Atmosphere* (1964) — closed-form lifetime baseline.
- Vallado & McClain, *Fundamentals of Astrodynamics and Applications*, 4th ed. — drag modeling, C_D = 2.2 convention.
- U.S. Standard Atmosphere (1976) / NRLMSISE-00 (Picone et al., 2002) — ρ₀ = 6.97×10⁻¹³ kg/m³ at 500 km, H ≈ 50 km.
- IADC Space Debris Mitigation Guidelines (IADC-02-01 Rev. 3, 2021); FCC Report & Order 22-74 (5-year rule); UNOOSA LTS Guidelines (2019) — regulatory framing.
- ESA Annual Space Environment Report 2024 — debris population statistics.
- Full citation list: references [1]–[23] of our IEEE-format paper (`docs/EcoOrbit_IEEE_Report.pdf`).

## 5. Prior work by the team

The Phase 1 concept submission (same challenge, same core project) — refined, implemented, and validated for the Final Phase per Rule 8.1 of the Finalist Guide. No external or third-party prior projects were reused.

---
*Team EcoOrbit — AESS Egypt · AESS Jordan · AESS Tunisia · AESS | ESSTHS. All registered members can explain, run, and modify every file in this repository.*
