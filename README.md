# ftsim — N point charges under Finite Theory's time-dilation factor

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
  is `G·M/d − Φ_spin` and its acceleration `(G·M/d − Φ_spin)/d`, i.e. `Φ_spin/d`
  less inward, so the net potential is kept. For Jupiter `Φ_spin/(GM/d) =
  Φ_spin/v_orb² = 20/171 (km/s)² = 12 %` less pull — a circular start drifts
  outward (Restart re-circularises on the new law). The same Φ enters
  `f = 1/(1 + Φ/c²)` as `η_g → η_g·(1 − Φ_spin/c²)` for that body (~2e-10).

| body | v_s km/s | Φ_spin/c² |
|---|---|---|
| Sun | 0.534 | 1.6e-12 |
| Earth | 0.267 | 4.0e-13 |
| Jupiter | 6.20 | 2.1e-10 |
| Saturn | 4.37 | 1.1e-10 |

The axis-speed effect alone is ≤ 2e-10 at gain 1; the potential subtraction is
visible at once on Jupiter and Saturn.

### Shadow planets

A darker copy of the Sun and planets is integrated in parallel with mainstream
orbital mechanics — `a = −G Σ m_j r_ij / r³` (× gravity gain), `p = m·v`, no FT
factor, no spin terms — from the same starting positions on the Newtonian
circular speed, and drawn under the live bodies with dashed trails. The gap
between a body and its shadow is what the selected route, η_g, kinematics and
spin options change: at η_g = c²/G it is 6e-7 AU for Mercury after ten orbits
(the FT perihelion-type drift), and 0.5 AU for Jupiter after 3000 steps with the
spin subtraction on. The body matrix in the group box lists, per planet, the
live-minus-shadow spherical offset about the Sun — Δr in AU, Δθ (azimuth) and
Δφ (polar angle) in degrees — refreshed every frame.

## License

Copyright (C) 2026 Phil Bouchard <phil@fornux.com>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See `COPYING` for the full text.
