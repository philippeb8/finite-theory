# ftsim — N point charges under Finite Theory's time-dilation factor

## Build

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/ftsim
```

Requires Qt6 Widgets (a Qt5 fallback is commented in `CMakeLists.txt`).
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

## License

Copyright (C) 2026 Phil Bouchard <phil@fornux.com>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See `COPYING` for the full text.
