# ftsim — N point charges under Finite Theory's time-dilation factor

![ftsim in action](ftsim.gif)

## Build

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/ftsim
```

Requires Qt6 Widgets, OpenGL and OpenGLWidgets (Debian/Ubuntu:
`qt6-base-dev libqt6opengl6-dev`; a Qt5 fallback is commented in
`CMakeLists.txt`). The view is a `QOpenGLWidget` on a fixed-function OpenGL
2.1 compatibility context, which every desktop driver provides.

## 3D view

All positions, velocities, accelerations and fields are 3-vectors (`Vec3`);
the integrators, the FT kernels, Biot–Savart, the spin dipole terms, the
gravitomagnetic field and the photon/pair model all run in three dimensions.
The planets are placed on their real orbital inclinations and ascending nodes
(J2000) with their spin axes at the real obliquities; atomic/quark scenes fill
a cube or a ball; Bohr stationary states stay in the xy plane by construction.

The camera orbits a target point: **left-drag** orbits (yaw about z, pitch
above the xy plane), **right-drag / Shift-drag** pans in the camera plane,
**wheel** zooms about the cursor, **double-click** resets (pitch 62°; pitch 90°
is the old top-down view). Grid rings lie in the xy plane; the axes are drawn
x amber, y cyan, z green; the crossover sphere as three great circles. The
scene is handed to the GPU in "target pixels" — `(r − target)·scale` — so a
10³⁰ a0 galaxy and a 10⁻⁵ a0 hadron are rendered with the same float
precision. Trails, rings and discs are OpenGL; symbols, HUD, legend and the
field compass are a 2D overlay blitted on top.

`FTSIM_SHOT=file.png[:tab]` runs for 1.5 s, saves a screenshot and quits
(for headless checks, e.g. under `xvfb-run` with `LIBGL_ALWAYS_SOFTWARE=1`).

Field dials: E and g are set in the xy plane (angle from +x); B and B_g in the
xz plane (tilt from +z), and their in-plane component now acts.
`core_test.cpp` is a standalone check of the acceleration kernel and needs no Qt:

```sh
g++ -O2 -std=c++17 core_test.cpp -o core_test && ./core_test
```

## Physics

`f(r) = eta*r / (eta*r + S)`, with `S = |q_source|`.

| route | acceleration | behaviour |
|---|---|---|
| Potential | `-C·eta²(S − eta·r)/(eta·r + S)³` | sign change at `r = S/eta` |
| Centripetal | `C·eta²/(eta·r + S)²` | sign-definite, no crossover |
| Coulomb | `C/r²` | reference |

`C = K·q_i·q_j/m_i`, and `K = 1` in the unit system below.

Both FT routes are **finite at r = 0**, so the integrator needs no softening
length — the `f` factor regularises the Coulomb singularity by itself.

## Units

| quantity | unit | value |
|---|---|---|
| length | a0 | 5.29177210903e-11 m |
| speed | v_Bohr | 2.18769126e6 m/s |
| time | a0/v_Bohr | 2.4189e-17 s |
| charge | e | 1.602176634e-19 C |
| mass | m_e | 9.1093837015e-31 kg |

`ETA = eta·a0/e`. **ETA = 8.7705e5** corresponds to `eta = c·eps0`; the
crossover is then at 1.14e-6 a0 and the run is indistinguishable from Coulomb
(ratio to Coulomb at 1 a0 = 0.999995439). Set ETA to roughly 0.1–10 to make
the FT structure visible on screen.

## Pair Production tab: the photon as a bound e⁻e⁺ pair (internal doc, Section 3)

With *photon = bound e⁻e⁺ pair* ticked, each photon is two real charged point
particles (an electron and a positron) co-moving at `c` along the beam at a
fixed separation `d0` across it (default `2 r_e`, Sec. 3.2) — no wave and no
internal oscillation; the frequency only sets the energy `ε = hf/m_e c²`.
Their motion is prescribed; they are real sources (other charges feel the
dipole field) and the ambient E/B sliders act on them.

The binding is the FT kernel of the current route with the current `eta`.
On the Potential route `|F(d)| = K e² eta² (eta d − 1)/(eta d + 1)³` for
`d > r_c = 1/eta`, and its maximum sits at exactly `d = 2 r_c`; with the
Section 3 preset `eta = e/r_e` that is `d0 = 2 r_e`, so

| | FT (Potential, eta = e/r_e) | Coulomb |
|---|---|---|
| F_bind at 2 r_e | `K e²/(27 r_e²)` = 1.08 N | `K e²/(4 r_e²)` = 7.26 N |
| U_bind at 2 r_e | `2/9 m_e c²` = 113.6 keV | `½ m_e c²` = 255.5 keV |
| E_split = F_bind / e | 6.7e18 V/m = 5.1 E_S | 4.5e19 V/m = 35 E_S |

Each lepton feels `e E_eff` outward (opposite charges, opposite ways) against
`F_bind` inward, so the pair splits when `e E_eff > F_bind`; since `d0` is the
force maximum, no barrier remains beyond it. `E_eff = |E + c k × B|`, any
direction (two point charges are sheared apart along the beam just as across it). The
capacitor of Sec. 3.2 (2.555 MV/m) gives `e E / F_bind = 3.8e-13`.

## Solar tab: body spin

Each body carries its rotational kinetic **potential** (per unit mass, no mass
in it): `Φ_spin = K_spin/m = ½·k·v_eq²`, `k = I/(mR²)`, from the sidereal
rotation period, equatorial radius and moment-of-inertia factor (NASA fact
sheets). Two options:

- **spin speed on the axis** (on by default) — `v_s = √(2Φ_spin)` is added on the rotation
  axis (z, out of the orbital plane): the kinematic factor sees `v² + v_s²`
  and the energy counts `m·Φ_spin`. The in-plane momentum inversion keeps its
  closed form with `1 → 1 − v_s²/c²` (SR, FT) or `1 − v_s²/2c²` (FT ½).
- **lose Φ_spin from the Sun's potential** — the Sun's potential at the body
  is the offset `G·M/d − Φ_spin` (per unit mass). A constant offset has no
  gradient, so the Newtonian pull is unchanged; it enters only the FT factor
  `f = 1/(1 + Φ/c²)`, as `η_g → η_g·(1 − Φ_spin/c²)` for that body (2.1e-10
  at Jupiter).

| body | v_s km/s | Φ_spin/c² |
|---|---|---|
| Sun | 0.534 | 1.6e-12 |
| Earth | 0.267 | 4.0e-13 |
| Jupiter | 6.20 | 2.1e-10 |
| Saturn | 4.37 | 1.1e-10 |

Both effects are ≤ 2e-10 at gain 1 — raise the gravity gain to see them.

- **axial force from the spin threshold law** (on by default) — the law
  fitted in FX-CDLIFT-2026-10 Rev D, per unit mass an acceleration along the
  spin vector, `a = a_A·[(ω/ω0)^(5/2) − 1]` for ω > ω0 and zero below, where
  a_A = 6.35 m/s² and ω0 = 243.8 rad/s (2330 rpm). It replaces the earlier
  spin-down force `√k·R·dω/dt`. Every solar body spins far below ω0 (Jupiter
  1.8e-4 rad/s, 7e-7 of ω0), so the force is zero for all of them; the spin
  table's `a_ax` column shows it per body. Scaling by equatorial speed instead
  (v0 = ω0 × 60 mm = 14.6 m/s) would give Earth 3.6e4 m/s² (6e6 times the
  Sun's pull) and push Jupiter past c within 20 minutes, so it is not offered.
  Counted in E as −m·a·(axis·r).

### Shadow planets

A darker copy of the Sun and planets is integrated in parallel with mainstream
orbital mechanics — `a = −G Σ m_j r_ij / r³` (× gravity gain), `p = m·v`, no FT
factor, no spin terms — from the same starting positions on the Newtonian
circular speed, and drawn under the live bodies with dashed trails. The gap
between a body and its shadow is what the selected route, η_g, kinematics and
spin options change: at η_g = c²/G it is 6e-7 AU for Mercury after ten orbits
(the FT perihelion-type drift); the spin options add ~2e-10 relative to that. The body matrix in the group box lists, per planet, the
live-minus-shadow spherical offset about the Sun — Δr in AU, Δθ (azimuth) and
Δφ (polar angle) in degrees — refreshed every frame.

## Section 2.8 (perihelion shift) in the simulator

Manuscript Eq. 66 is the GR effective potential with its metric factor
replaced by f²: `V_eff = −GMm/r + L²f²/(2mr²)`, and Sec. 2.8 quotes the GR
result `δφ = 6πGM/(c²a(1−e²))` for it. Three readings were integrated (an
eccentric Mercury, e = 0.2056, G amplified so the angle resolves):

| reading | φ̇ | δφ per orbit | Mercury |
|---|---|---|---|
| A — Eq. 66 as written, Newtonian `L = m r² φ̇` | `L/(m r²)` | **6π g** | 42.98″/cy |
| B — the Lagrangian whose effective potential is Eq. 66: `½mṙ² + ½mr²φ̇²/f² + GMm/r` | `L f²/(m r²)` | **2π g** | 14.3″/cy |
| C — `1 − 2GM/c²r → f²` everywhere in the Schwarzschild radial equation | `L/(m r²)` | **3π g** | 21.5″/cy |

`g = GM/(c² a(1−e²))`. Reading A reproduces GR, but it pairs the f²-modified
potential with the unmodified angular momentum, and no Lagrangian does that:
once `L²f²/(2mr²)` is the centrifugal term, the conserved momentum is
`m r² φ̇/f²` and the azimuth advances as `L f²/(m r²)` — that costs −4π g and
leaves 2π g. Reading A is GR's radial equation transcribed, not an equation of
motion of FT. A pure time reparametrisation `dτ = f·dt` of Kepler motion gives
no precession at all (the orbit shape does not depend on the clock).

### The velocity-dependent reading

Sec. 2.8 closes by saying the Eq. 66 term "can appear as the result of the
velocity-dependent correction". Taking that literally — the body's clock runs
at f, so its own velocity is u = v/f and its mechanics is written in its own
time, S = ∫ L(u) dτ with dτ = f·dt — gives Lagrangians of the form
`f·L_kin(v/f) + (gravity)`. Their precession (same integrator, per orbit, in
units of g):

| Lagrangian (per unit mass) | δφ |
|---|---|
| metric `−c²√(f² − v²/c²)` (gravity from f alone, no U) | π |
| `v²/(2f) + U` | 2π |
| `v²/(2f²) + U` | 4π |
| `f·L_41(v/f) + U`, Eq. 41 kinetics `L_41 = −(c²/2) ln(1 − u²/c²)` | 4π |
| `f·L_½(v/f) + U`, FT ½ kinetics | 3π |
| `f·L_41(v/f) + U·f²` (the simulator's Potential route) | 0 |
| **`f·L_41(v/f) + U/f`** | **6π** |

So the velocity-dependent correction with FT's own Eq. 41 momentum does
produce two thirds of GR by itself, and the full 6π follows if the potential
term is `U/f` — i.e. `U/f²` in the body's proper time, `S = ∫[L_41(u) + U/f²] dτ`.
That combination is implemented as the momentum factor **"FT proper time"**
(use with the Coulomb route): `p = m(v/f)/(1 − v²/f²c²)`, Hamiltonian
`H = f·K_41(u) − U/f`, forces Newton's plus `−K_41(u)∇f + ∇(U/f − U)`.
Measured +18.99 g per orbit at g = 2.7·10⁻³ (6π = 18.85, the rest is second
order), energy conserved to 5·10⁻⁷ and scaling with dt.

The simulator also implements reading B as the momentum factor **"Sec. 2.8
Lagrangian"** (use with the Coulomb route): `p_r = m v_r`, `p_t = m v_t/f²`
w.r.t. the heaviest body, with the two extra Hamiltonian forces
`−(p_t²/2m)·d(f²)/dr` (radial) and `(f²−1) p_r p_t/(m r)` (tangential, keeps
`L = r p_t`). Measured: +6.28 g per orbit with E conserved to leapfrog order.

## Freefall tab: steel vs NdFeB

Two 66 g balls, steel and NdFeB, released together from height h in Earth's g
onto a floor at z = 0 (Newtonian kinematics, Coulomb route; side view, x right,
z up; dt = 1e13 = 0.24 ms). The magnet carries two options:

- **FT magnetic acceleration** `a_B = 2·B·c/√(K/G)` along its magnetisation
  axis (S→N), `√(K/G) = 1.1606·10¹⁰ kg/C`: 0.0672 m/s² at B = 1.3 T, 0.69 % of
  g. Default N pole up: a = g − a_B, the magnet falls slower; N pole down adds
  it; sideways it drifts sideways. Dial sets the axis tilt from straight down.
- **field energy as kinetic potential** — the spinning-disk analogue:
  `Φ_B = W/m`, `W = B²V/(6μ₀)` for a uniformly magnetised sphere (66 g NdFeB:
  2.0 J, 29.9 J/kg), `v_B = √(2Φ_B) = 7.7 m/s` on the axis in the kinematic
  factor, and the potential offset Φ_g − Φ_B. Both are ~10⁻¹⁶ of g; B is
  constant during the fall, so the axial force m·dv_B/dt is zero.

The HUD gives predicted and run impact times for comparison with a drop-test
regression. From 1.5 m with the default N up: steel 0.55310 s, magnet
0.55500 s, the magnet lands 1.90 ms later and is 10.3 mm higher when the steel
ball lands (Δa/g = −6.85·10⁻³). N down: 1.88 ms earlier. Run impact times are interpolated
inside the last step and match the analytic values to 10 µs.

**The magnet swings.** A 66 g, 1.3 T sphere is a rigid dipole of m = 9.1 A·m²
with I = 4.3·10⁻⁶ kg·m²; in Earth's field (5.4·10⁻⁵ T, inclination 71°) its
libration period is 2π√(I/mB) = 0.59 s — the fall time itself. N down sits near
the stable direction (21° off vertical) and oscillates about it; N up is the
unstable equilibrium with an e-folding time of 0.09 s, so a 1° release error
grows to a flip before the floor. Since a_B rides on the axis, the deficit it
produces is not the fixed-axis value: with the torque on (default, 1°
perturbation at release) the N-up lead of +1.90 ms becomes +0.87 ms, the N-down
lag −1.88 → −1.74 ms, and a horizontal release gives −1.10 ms instead of zero.
A fixed-axis prediction therefore overstates what a_B can produce, and the
pole orientation at impact is not the one at release.

Two extra knobs for fitting a measured deficit: the coefficient β in
`a_B = β·B·c/√(K/G)` (2 as postulated, 2/3 in the Oct 2025 note;
β = Δa_meas·√(K/G)/(B·c) = Δa_meas × 29.78 at 1.3 T), and a velocity drag
`a = −κ·v` on the magnet (eddy currents in nearby conductors act this way).
The two forms separate with drops from several heights: a constant a_B gives a
fixed Δt/t at every h, a drag gives a Δt/t that grows with h (tuned to the same
deficit at 1.5 m: 0.6 % at 0.5 m vs 1.4 % at 3 m). A fitted β is a measurement
only if it is the same for two magnets of different B.

Unit fix made here: the atomic unit of acceleration is a0/t² = 9.0442·10²² m/s²
(the earlier rounded 9.0·10²² mis-scaled g by 0.5 %).

## Energy audit

`E = kinetic() + potential()` is the quantity the integrator conserves, with
`kinetic()` the exact `∫ v·dp` of the selected momentum law. With the body spin
on, `p = g(v² + v_s²)·m·v` at fixed `v_s`, and that integral is **not**
`K_kind(√(v² + v_s²))`: `dp/dv` carries a factor `(1 − v_s²/c²)`, so

| kind | K (per body) |
|---|---|
| SR | `m c² (1 − β_s²)[γ(u) − γ(v_s)] + m c²[γ(v_s) − 1]` |
| FT | `m c² [(1 − β_s²)(G(u) − G_s) − ½ ln(G(u)/G_s)] + m c²[G_s − 1 − ½ ln G_s]` |
| FT ½ | `m c² [(2 − β_s²)(1/D(u) − 1/D_s) + ln(D(u)/D_s)] + m c²[2/D_s − 2 + ln D_s]` |

with `u² = v² + v_s²`, `G = 1/(1 − u²/c²)`, `D = 1 − u²/2c²`, the second term
being the spin's own (constant) energy. All differences are evaluated in
cancellation-free form (`log1p`, `b²/(√(1−b²)(1+√(1−b²)))`): the Sun's
`(γ − 1) M c²` computed naively is round-off of 2e64 that swamps the planets'
1e51 orbital energies by 1e-4.

Max |ΔE| / (|K| + |U|) at the tab's dt and at dt/2 (leapfrog: ×¼ expected):

| scene | dt | dt/2 | note |
|---|---|---|---|
| Solar, spin on, circular | 3.8e-11 | 8.4e-12 | |
| Solar, Mercury at 0.7 v_circ, v_s = 0.6c on Mercury | 2.5e-3 | 6.3e-4 | the naive `K_kind(u)` drifts by 0.6 of Mercury's KE |
| Atomic, 60 bodies, electric + gravity | 1.9e-4 | 4.6e-5 | |
| Atomic, + Biot–Savart, spin moments, gravitomagnetic | 2.0e-1 | 1.1e-2 | velocity-dependent forces have no potential: not conserved by construction |
| Quark, uud Cornell, SR | 1.2e-13 | 3.1e-14 | was a dt-independent 6.4e-9: the force and potential used unit constants differing in the 7th digit, now one constant |
| Pair, freed e⁻e⁺ under the FT potential | 2.3e-3 | 5.7e-4 | photons are not in `E`; it jumps at a split or recombination |
| Galactic, bulge + 200 tracers | 1.1e-5 | 2.7e-5 | does not scale: unresolved tracer–tracer encounters at 5–9 pc (one step ≈ 24 Myr), not a conservation bug |

Every momentum factor was rechecked with a large test spin (v_s = 0.3c on an
eccentric Mercury): the drift equals the spin-free drift and scales as dt²
for all nine. SR, FT Eq. 41, FT ½, FT f²×(m0 + K/c²) and FT proper time carry
the spin in the momentum law (`u² = v² + v_s²`; for FT proper time
`u² = v²/f² + v_s²`, the spin speed in the body's own frame, so the spin
energy also weighs through f). Newtonian, the two unlimited f² kinds and the
Sec. 2.8 Lagrangian keep it as the separate constant `½ m v_s²` — conserved,
but with no effect on the orbit.

With the spin offset on, one η_g per pair (the lighter body's) is used for both
members, so the pair force stays equal and opposite and derives from
`potential()`.

## License

Copyright (C) 2026 Phil Bouchard <phil@fornux.com>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See `COPYING` for the full text.
